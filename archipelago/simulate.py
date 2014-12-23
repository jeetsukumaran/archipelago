#! /usr/bin/env python

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import random
import collections
import argparse
import pprint
import copy
from distutils.util import strtobool

import dendropy
from dendropy.utility import textprocessing

import archipelago
from archipelago import utility

def weighted_choice(seq, weights, rng=None):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `weights` (which must be at least as long as the
    length of `seq` - 1).
    """
    if weights is None:
        weights = [1.0/len(seq) for count in range(len(seq))]
    else:
        weights = list(weights)
    if len(weights) < len(seq) - 1:
        raise Exception("Insufficient number of weights specified")
    if len(weights) == len(seq) - 1:
        weights.append(1 - sum(weights))
    return seq[weighted_index_choice(weights, rng)]

def weighted_index_choice(weights, rng=None):
    """
    (From: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    The following is a simple function to implement weighted random choice in
    Python. Given a list of weights, it returns an index randomly, according
    to these weights [1].
    For example, given [2, 3, 5] it returns 0 (the index of the first element)
    with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5.
    The weights need not sum up to anything in particular, and can actually be
    arbitrary Python floating point numbers.
    If we manage to sort the weights in descending order before passing them
    to weighted_choice_sub, it will run even faster, since the random call
    returns a uniformly distributed value and larger chunks of the total
    weight will be skipped in the beginning.
    """
    rnd = rng.uniform(0, 1) * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

class InsufficientFocalAreaLineagesSimulationException(Exception):
    pass

class TotalExtinctionException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class IndexGenerator(object):

    def __init__(self, start=0):
        self.start = start
        self.index = start

    def __next__(self):
        c = self.index
        self.index += 1
        return c
    next = __next__

    def reset(self, start=None):
        if start is None:
            start = self.start
        self.index = start

class StatesVector(object):
    """
    A vector in which each element is an integer represents the state of a
    trait.

    E.g.,

        [1,0,1,2]

    is a 4-trait vector, where trait 0 is in state 1, trait 1 is in
    state 0, and so on.
    """

    def __init__(self,
            nchar,
            nstates=None,
            values=None,
            ):
        """
        Parameters
        ----------
        nchar : integer
            The number of traits to be tracked.
        nstates : list of integers
            The number of states for each trait. If not specified, defaults
            to binary (i.e., 2 states, 0 and 1). If specifed, must be a list of
            length `nchar`, with each element in the list being integer > 0.
        values : iterable of ints
            Vector of initial values. If not specified, defaults to all 0's.
        """
        self._nchar = nchar
        if nstates is None:
            self._nstates = nstates
        else:
            self._nstates = [2] * nchar
        if not values:
            self._states = [0] * nchar
        else:
            assert len(values) == nchar
            self._states = list(values)

    def clone(self):
        s = self.__class__(
                nchar=self._nchar,
                nstates=self._nstates,
            )
        s._states = list(self._states)
        return s

    @property
    def nchar(self):
        return len(self)

    def __len__(self):
        return self._nchar

    def __getitem__(self, trait_index):
        return self._states[trait_index]

    def __setitem__(self, trait_index, v):
        self._states[trait_index] = v

class DistributionVector(StatesVector):

    def __init__(self, num_areas):
        StatesVector.__init__(self,
                nchar=num_areas,
                nstates=[2] * num_areas,
                )

    def presences(self):
        """
        Returns list of indexes in which lineage is present.
        """
        return [idx for idx, s in enumerate(self._states) if s == 1]

    def clone(self):
        s = self.__class__(num_areas=self._nchar)
        s._states = list(self._states)
        return s

class TraitType(object):

    def __init__(self,
            index=None,
            label=None,
            nstates=None,
            transition_rate=None,
            transition_weights=None,
            ):
        self.index = index
        self.label = label
        self.nstates = nstates
        self.transition_rate = transition_rate
        self.transition_weights = transition_weights
        self.transition_rate_matrix = []

    def compile_matrix(self):
        self.transition_rate_matrix = []
        for i, tw in enumerate(self.transition_weights):
            self.transition_rate_matrix.append([])
            for j, w in enumerate(tw):
                self.transition_rate_matrix[i].append( w * self.transition_rate )

class TraitTypes(object):

    def __init__(self):
        self.normalize_transition_weights = True

    def __iter__(self):
        return iter(self.trait_types)

    def __getitem__(self, idx):
        return self.trait_types[idx]

    def parse_definition(self,
            trait_types,
            run_logger,
            verbose=True):
        self.trait_types = []
        self.trait_label_index_map = collections.OrderedDict()
        for trait_idx, trait_d in enumerate(trait_types):
            trait = TraitType(
                index=trait_idx,
                label=str(trait_d.pop("label", trait_idx)),
                nstates=trait_d.pop("nstates", 2),
                transition_rate=trait_d.pop("transition_rate", 0.01),
            )
            self.trait_types.append(trait)
            transition_weights = trait_d.pop("transition_weights", None) # delay processing until all traits have been defined
            if not transition_weights:
                trait.transition_weights = [[1.0 for j in range(trait.nstates)] for i in range(trait.nstates)]
            total_transition_weight = 0.0
            for a1_idx in range(trait.nstates):
                for a2_idx in range(trait.nstates):
                    if a1_idx == a2_idx:
                        trait.transition_weights[a1_idx][a2_idx] = 0.0
                    else:
                        trait.transition_weights[a1_idx][a2_idx] = float(trait.transition_weights[a1_idx][a2_idx])
                        total_transition_weight += trait.transition_weights[a1_idx][a2_idx]
            if self.normalize_transition_weights and total_transition_weight:
                for a1_idx in range(trait.nstates):
                    for a2_idx in range(trait.nstates):
                        if a1_idx == a2_idx:
                            continue
                        trait.transition_weights[a1_idx][a2_idx] /= total_transition_weight
            self.trait_label_index_map[trait.label] = trait.index
            if verbose:
                run_logger.info("[ECOLOGY] Configuring trait {idx}: '{label}': {nstates} states, transition rate of {trate} with {normalized}transition weights of {tweights}".format(
                    idx=trait.index,
                    label=trait.label,
                    normalized="normalized " if self.normalize_transition_weights else "",
                    nstates=trait.nstates,
                    trate=trait.transition_rate,
                    tweights=trait.transition_weights,
                    ))
            trait.compile_matrix()
            if trait_d:
                raise TypeError("Unsupported trait keywords: {}".format(trait_d))
        # if len(self.trait_types) < 1:
        #     raise ValueError("No traits defined")
        if verbose:
            run_logger.info("[ECOLOGY] {} traits defined: {}".format(
                len(self.trait_types),
                ", ".join("'{}'".format(a.label) for a in self.trait_types),
                ))
        self.trait_nstates = [trait.nstates for trait in self.trait_types]

    def new_traits_vector(self, values=None):
        s = StatesVector(
                nchar=len(self.trait_types),
                nstates=self.trait_nstates)
        return s

class Area(object):

    def __init__(self,
            index=None,
            label=None,
            is_supplemental=False,
            relative_diversity=None,
           transition_weights=None,
            ):
        self.index = index
        self.label = label
        self.is_supplemental = is_supplemental
        self.relative_diversity = relative_diversity

class Geography(object):

    def __init__(self):
        self.normalize_dispersal_weights = True

    def __iter__(self):
        return iter(self.areas)

    def __getitem__(self, idx):
        return self.areas[idx]

    def parse_definition(self,
            areas,
            run_logger,
            verbose=True):
        self.areas = []
        self.area_label_index_map = collections.OrderedDict()
        self.area_indexes = []
        self.focal_area_indexes = []
        self.supplemental_area_indexes = []
        for area_idx, area_d in enumerate(areas):
            area = Area(
                index=area_idx,
                label=str(area_d.pop("label", area_idx)),
                relative_diversity=area_d.pop("relative_diversity", 1.0),
                is_supplemental=area_d.pop("is_supplemental", False)
            )
            area._dispersal_weights_d = area_d.pop("dispersal_weights", {}) # delay processing until all areas have been defined
            self.areas.append(area)
            self.area_label_index_map[area.label] = area.index
            self.area_indexes.append(area.index)
            if area.is_supplemental:
                self.supplemental_area_indexes.append(area.index)
            else:
                self.focal_area_indexes.append(area.index)
            if verbose:
                run_logger.info("[GEOGRAPHY] Configuring area {idx}: '{label}' ({is_supplemental}; relative diversity: {relative_diversity})".format(
                    idx=area.index,
                    label=area.label,
                    is_supplemental="supplemental" if area.is_supplemental else "primary",
                    relative_diversity=area.relative_diversity,
                    ))
            if area_d:
                raise TypeError("Unsupported area keywords: {}".format(area_d))
        if len(self.areas) < 1:
            raise ValueError("No areas defined")
        if verbose:
            run_logger.info("[GEOGRAPHY] {} areas defined: {}".format(
                len(self.areas),
                ", ".join("'{}'".format(a.label) for a in self.areas),
                ))
        self.dispersal_weights = []
        total_dispersal_weight = 0.0
        for a1_idx, area1 in enumerate(self.areas):
            self.dispersal_weights.append([])
            for a2_idx, area2 in enumerate(self.areas):
                if a1_idx == a2_idx:
                    self.dispersal_weights[a1_idx].append(0.0)
                else:
                    d = float(area1._dispersal_weights_d.pop(area2.label, 1.0))
                    self.dispersal_weights[a1_idx].append(d)
                    total_dispersal_weight += d
            if area1._dispersal_weights_d:
                raise ValueError("Undefined dispersal targets in '{}': '{}'".format(area1.label, area1._dispersal_weights_d))
            del area1._dispersal_weights_d
        if self.normalize_dispersal_weights and total_dispersal_weight:
            for a1_idx, area1 in enumerate(self.areas):
                for a2_idx, area2 in enumerate(self.areas):
                    self.dispersal_weights[a1_idx][a2_idx] /= total_dispersal_weight
        if verbose:
            if self.normalize_dispersal_weights:
                weight_type = "Normalized dispersal"
            else:
                weight_type = "Dispersal"
            for a1, area1 in enumerate(self.areas):
                run_logger.info("[GEOGRAPHY] {} weights from area '{}': {}".format(weight_type, area1.label, self.dispersal_weights[a1]))

        # instead of recalculating every time
        self.area_nstates = [2 for i in self.areas]

    def new_distribution_vector(self):
        s = DistributionVector(num_areas=len(self.areas))
        return s

class Lineage(dendropy.Node):

    def __init__(self,
            index,
            distribution_vector=None,
            traits_vector=None,
            ):
        dendropy.Node.__init__(self)
        self.index = index
        self.distribution_vector = distribution_vector
        self.traits_vector = traits_vector
        self.is_extant = True
        self.edge.length = 0

    def __repr__(self):
        return utility.encode_lineage(self)

class Phylogeny(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return Lineage(**kwargs)
    node_factory = classmethod(node_factory)

    def __init__(self, *args, **kwargs):
        if "system" in kwargs:
            self.system = kwargs.pop("system")
            self.lineage_indexer = IndexGenerator(0)
            if "seed_node" not in kwargs:
                seed_node = self.node_factory(
                        index=next(self.lineage_indexer),
                        distribution_vector=self.system.geography.new_distribution_vector(),
                        traits_vector=self.system.trait_types.new_traits_vector(),
                        )
                seed_node.distribution_vector[0] = 1
                kwargs["seed_node"] = seed_node
            dendropy.Tree.__init__(self, *args, **kwargs)
            self.current_lineages = set([self.seed_node])
        else:
            dendropy.Tree.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        memo[id(self.system)] = self.system
        memo[id(self.taxon_namespace)] = self.taxon_namespace
        return dendropy.Tree.__deepcopy__(self, memo)

    def iterate_current_lineages(self):
        for lineage in self.current_lineages:
            yield lineage

    def split_lineage(self, lineage):

        # speciation modes
        # 0:  single-area sympatric speciation
        #     -   ancestral range copied to both daughter species
        # 1:  sympatric subset: multi-area sympatric speciation
        #     -   d1: inherits complete range
        #     -   d2: inherits single area in ancestral range
        # 2:  vicariance
        #     -   ancestral range divided up between two daughter species
        # 3:  jump dispersal
        #     -   single
        presences = lineage.distribution_vector.presences()
        num_presences = len(presences)
        num_areas = len(self.system.geography.area_indexes)
        if num_presences <= 1:
            speciation_mode = 0
        elif num_presences == num_areas:
            speciation_mode = self.system.rng.randint(0, 2)
        else:
            speciation_mode = self.system.rng.randint(0, 3)
        if speciation_mode == 0:
            dist1 = lineage.distribution_vector.clone()
            dist2 = lineage.distribution_vector.clone()
        elif speciation_mode == 1:
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.system.geography.new_distribution_vector()
            # TODO: area diversity base speciation
            dist2[ self.system.rng.choice(presences) ] = 1
        elif speciation_mode == 2:
            dist1 = self.system.geography.new_distribution_vector()
            dist2 = self.system.geography.new_distribution_vector()
            if num_presences == 2:
                dist1[presences[0]] = 1
                dist2[presences[1]] = 1
            else:
                n1 = self.system.rng.randint(1, num_presences-1)
                n2 = num_presences - n1
                if n2 == n1:
                    n1 += 1
                    n2 -= 1
                sample1 = set(self.system.rng.sample(presences, n1))
                for idx in self.system.geography.area_indexes:
                    if idx in sample1:
                        dist1[idx] = 1
                    else:
                        dist2[idx] = 1
        elif speciation_mode == 3:
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.system.geography.new_distribution_vector()
            absences = [idx for idx in self.system.geography.area_indexes if idx not in presences]
            dist2[ self.system.rng.choice(absences) ] = 1
        else:
            raise ValueError(speciation_mode)

        c1 = self.node_factory(
                index=next(self.lineage_indexer),
                distribution_vector=dist1,
                traits_vector=lineage.traits_vector.clone(),
                )
        c2 = self.node_factory(
                index=next(self.lineage_indexer),
                distribution_vector=dist2,
                traits_vector=lineage.traits_vector.clone(),
                )
        if self.system.debug_mode:
            self.system.run_logger.debug("Splitting {} with distribution {} under speciation mode {} to: {} (distribution: {}) and {} (distribution: {})".format(
                lineage,
                lineage.distribution_vector.presences(),
                speciation_mode,
                c1,
                dist1.presences(),
                c2,
                dist2.presences(),
                ))
            assert len(dist1.presences()) > 0
            assert len(dist2.presences()) > 0

        lineage.is_extant = False
        self.current_lineages.remove(lineage)
        lineage.add_child(c1)
        lineage.add_child(c2)
        self.current_lineages.add(c1)
        self.current_lineages.add(c2)

    def extirpate_lineage(self, lineage):
        if self.system.is_lineage_death_global:
            self._make_lineage_extinct_on_phylogeny(lineage)
        else:
            presences = lineage.distribution_vector.presences()
            assert len(presences) > 0
            if len(presences) == 1:
                self._make_lineage_extinct_on_phylogeny(lineage)
            else:
                lineage.distribution_vector[ self.system.rng.choice(presences) ] = 0

    def _make_lineage_extinct_on_phylogeny(self, lineage):
        if len(self.current_lineages) == 1:
            self.total_extinction_exception("No extant lineages remaining")
        lineage.is_extant = False
        self.current_lineages.remove(lineage)
        self.prune_subtree(lineage)

    def total_extinction_exception(self, msg):
        # self.run_logger.info("Total extinction: {}".format(msg))
        raise TotalExtinctionException(msg)

    def evolve_trait(self, lineage, trait_idx, state_idx):
        lineage.traits_vector[trait_idx] = state_idx

    def disperse_lineage(self, lineage, dest_area_idx):
        lineage.distribution_vector[dest_area_idx] = 1

    def focal_area_lineages(self):
        focal_area_lineages = set()
        for lineage in self.iterate_current_lineages():
            for area_idx in self.system.geography.focal_area_indexes:
                if lineage.distribution_vector[area_idx] == 1:
                    focal_area_lineages.add(lineage)
                    break
        return focal_area_lineages

    def num_focal_area_lineages(self):
        count = 0
        for lineage in self.iterate_current_lineages():
            for area_idx in self.system.geography.focal_area_indexes:
                if lineage.distribution_vector[area_idx] == 1:
                    count += 1
                    break
        return count

    def extract_focal_area_tree(self):
        # tcopy = Phylogeny(self)
        tcopy = copy.deepcopy(self)
        focal_area_lineages = tcopy.focal_area_lineages()
        if len(focal_area_lineages) < 2:
            raise InsufficientFocalAreaLineagesSimulationException("Insufficient lineages in focal area at termination: {}".format(len(focal_area_lineages)))
        try:
            tcopy.filter_leaf_nodes(filter_fn=lambda x: x in focal_area_lineages)
        except dendropy.SeedNodeDeletionException:
            raise InsufficientFocalAreaLineagesSimulationException("Insufficient lineages in focal area at termination: {}".format(len(focal_area_lineages)))
        return tcopy

class ArchipelagoSimulator(object):

    @staticmethod
    def get_fixed_value_function(v, description):
        f = lambda x: v
        f.__doc__ = description
        return f

    def __init__(self,
            model_d=None,
            config_d=None,
            verbose_setup=True
            ):

        # configure
        self.elapsed_time = 0.0
        if config_d is None:
            config_d = {}
        else:
            config_d = dict(config_d)
        self.configure_simulator(config_d, verbose=verbose_setup)

        # set up model
        if model_d is None:
            model_d = {}
        else:
            model_d = dict(model_d)
        self.set_model(model_d, verbose=verbose_setup)

        # start
        self.phylogeny = Phylogeny(system=self)

        # begin logging generations
        self.run_logger.system = self

    def configure_simulator(self, config_d, verbose=True):

        self.name = config_d.pop("name", None)
        if self.name is None:
            self.name = str(id(self))
        self.output_prefix = config_d.pop("output_prefix", "archipelago-{}".format(self.name))

        self.run_logger = config_d.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(
                    name="archipelago",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_path=self.output_prefix + ".log",
                    file_logging_level=config_d.pop("file_logging_level", "info"),
                    )
        self.run_logger.system = self

        if verbose:
            self.run_logger.info("Configuring simulation '{}'".format(self.name))

        if config_d.pop("store_focal_area_trees", True):
            self.focal_areas_tree_log = open(self.output_prefix + ".focal-areas.trees", "w")
            if verbose:
                self.run_logger.info("Focal area trees filepath: {}".format(self.focal_areas_tree_log.name))
        else:
            self.focal_areas_tree_log = None
            if verbose:
                self.run_logger.info("Focal area trees will not be stored")

        if config_d.pop("store_full_area_trees", True):
            self.all_areas_tree_log = open(self.output_prefix + ".all-areas.trees", "w")
            if verbose:
                self.run_logger.info("Full area trees filepath: {}".format(self.all_areas_tree_log.name))
        else:
            self.all_areas_tree_log = None
            if verbose:
                self.run_logger.info("Full area trees will not be stored")

        if not self.focal_areas_tree_log and not self.all_areas_tree_log:
            self.run_logger.warning("No trees will be stored!")

        self.is_suppress_internal_node_labels = config_d.pop("suppress_internal_node_labels", False)
        if verbose:
            self.run_logger.info("Internal node labels will{} be written on trees".format(" not" if self.is_suppress_internal_node_labels else ""))

        self.rng = config_d.pop("rng", None)
        if self.rng is None:
            self.random_seed = config_d.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            if verbose:
                self.run_logger.info("Initializing with random seed {}".format(self.random_seed))
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in config_d:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            if verbose:
                self.run_logger.info("Using existing random number generator")

        termination_conditions_d = config_d.pop("termination_conditions", {})
        self.target_num_tips = termination_conditions_d.pop("target_num_tips", 50)
        # self.gsa_termination_num_tips = termination_conditions_d.pop("gsa_termination_num_tips", 500)
        self.gsa_termination_num_tips = termination_conditions_d.pop("gsa_termination_num_tips", 0)
        self.max_time = termination_conditions_d.pop("max_time", 0)
        if termination_conditions_d:
            raise TypeError("Unsupported configuration keywords: {}".format(termination_conditions_d))
        if self.gsa_termination_num_tips and not self.target_num_tips:
            raise ValueError("Cannot specify 'gsa_termination_num_tips' without specifying 'target_num_tips'")
        if not self.target_num_tips and self.max_time:
            desc = "Simulation will terminate after {} time units".format(self.max_time)
        elif self.target_num_tips and not self.gsa_termination_num_tips and not self.max_time:
            desc = "Simulation will terminate when {} tips are generated (no time limit)".format(self.target_num_tips)
        elif self.target_num_tips and not self.gsa_termination_num_tips and self.max_time:
            desc = "Simulation will terminate when {} tips are generated or after {} time units".format(self.target_num_tips, self.max_time)
        elif self.target_num_tips and self.gsa_termination_num_tips and not self.max_time:
            desc = "Simulation will terminate when {} tips are generated, with a random snapshot of the phylogeny when there were {} extant tips will be sampled and returned".format(self.target_num_tips, self.gsa_termination_num_tips)
        elif self.target_num_tips and self.gsa_termination_num_tips and self.max_time:
            desc = "Simulation will terminate when {} tips are generated or after {} time units, with a random snapshot of the phylogeny when there were {} extant tips will be sampled and returned".format(self.target_num_tips, self.max_time, self.gsa_termination_num_tips)
        elif not self.target_num_tips and not self.max_time:
            raise ValueError("Unspecified termination condition")
        else:
            raise ValueError("Unsupported termination condition(s)")
        if verbose:
            self.run_logger.info(desc)

        if self.target_num_tips:
            default_log_frequency = 1
        else:
            default_log_frequency = self.max_time/100
        self.log_frequency = config_d.pop("log_frequency", default_log_frequency)

        self.debug_mode = config_d.pop("debug_mode", False)
        if verbose and self.debug_mode:
            self.run_logger.info("Running in DEBUG mode")

        if config_d:
            raise TypeError("Unsupported configuration keywords: {}".format(config_d))

    def set_model(self, model_d, verbose=True):

        # Geography
        if "areas" not in model_d:
            raise ValueError("No areas defined")
        self.geography = Geography()
        self.geography.parse_definition(
                model_d.pop("areas"),
                run_logger=self.run_logger,
                verbose=verbose)

        # Ecology
        self.trait_types = TraitTypes()
        self.trait_types.parse_definition(
                model_d.pop("traits", {}),
                run_logger=self.run_logger,
                verbose=verbose)

        # Diversification: speciation
        diversification_d = model_d.pop("diversification", {})
        if "lineage_birth_probability_function" in diversification_d:
            self.lineage_birth_probability_function = diversification_d.pop("lineage_birth_probability_function")
        else:
            self.lineage_birth_probability_function = ArchipelagoSimulator.get_fixed_value_function(
                    0.01,
                    "Fixed speciation probability: {}".format(0.01)
            )
        if verbose:
            desc = getattr(self.lineage_birth_probability_function, "__doc__", None)
            if desc is None:
                desc = "(no description available)"
            self.run_logger.info("[DIVERSIFICATION] Setting lineage speciation probability function: {}".format(desc,))

        # Diversification: death/extinction/extirpation
        if "lineage_death_probability_function" in diversification_d:
            self.lineage_death_probability_function = diversification_d.pop("lineage_death_probability_function")
        else:
            self.lineage_death_probability_function = ArchipelagoSimulator.get_fixed_value_function(
                    0.01,
                    "Fixed death probability: {}".format(0.01)
            )
        self.is_lineage_death_global = strtobool(str(model_d.pop("is_lineage_death_global", 0)))
        if verbose:
            desc = getattr(self.lineage_death_probability_function, "__doc__", None)
            if desc is None:
                desc = "(no description available)"
            self.run_logger.info("[DIVERSIFICATION] Setting lineage death (= {is_global}) probability function: {desc}".format(
                is_global="global extinction" if self.is_lineage_death_global else "local extirpation",
                desc=desc,
                ))

        # Dispersal submodel
        if "lineage_dispersal_probability_function" in model_d:
            self.lineage_dispersal_probability_function = model_d.pop("lineage_dispersal_probability_function")
        else:
            self.lineage_dispersal_probability_function = ArchipelagoSimulator.get_fixed_value_function(
                    0.01,
                    "Fixed dispersal probability: {}".format(0.01)
            )
        if verbose:
            desc = getattr(self.lineage_dispersal_probability_function, "__doc__", None)
            if desc is None:
                desc = "(no description available)"
            self.run_logger.info("[DISPERSAL] Setting lineage dispersal probability function: {}".format(desc,))

        if model_d:
            raise TypeError("Unsupported model keywords: {}".format(model_d))

    def run(self):

        ### Initialize time
        self.elapsed_time = 0.0
        if self.log_frequency:
            if self.target_num_tips:
                last_logged_num_tips = 0
            else:
                last_logged_time = 0.0

        ### Initialize debugging
        if self.debug_mode:
            num_events = 0

        ### Initialize termination conditiong checking
        ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
        ntips = len(self.phylogeny.current_lineages)

        while True:

            ### DEBUG
            if self.debug_mode:
                num_events += 1
                self.run_logger.debug("Pre-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    assert len(lineage.distribution_vector.presences()) > 0

            ### LOGGING
            if self.log_frequency:
                if self.target_num_tips:
                    if ntips_in_focal_areas - last_logged_num_tips >= self.log_frequency:
                        last_logged_num_tips = ntips_in_focal_areas
                        self.run_logger.info("{} lineages occurring in focal areas, {} lineages across all areas".format(ntips_in_focal_areas, ntips))
                else:
                    if self.elapsed_time - last_logged_time >= self.log_frequency:
                        last_logged_time = self.elapsed_time
                        self.run_logger.info("{} lineages occurring in focal areas, {} lineages across all areas".format(ntips_in_focal_areas, ntips))

            ### EVENT SCHEDULING
            event_calls, event_rates, sum_of_event_rates = self.schedule_events()

            if self.debug_mode:
                if sum_of_event_rates == 0:
                    self.run_logger.debug("Sum of event rates is 0: {}".format(event_rates))

            time_till_event = self.rng.expovariate(sum_of_event_rates)
            self.elapsed_time += time_till_event
            if self.max_time and self.elapsed_time > max_time:
                raise NotImplementedError
            for lineage in self.phylogeny.iterate_current_lineages():
                lineage.edge.length += time_till_event

            ### EVENT SELECTION AND EXECUTION
            event_idx = weighted_index_choice(event_rates, self.rng)
            if self.debug_mode:
                self.run_logger.debug("Event {}: {}".format(num_events, event_calls[event_idx]))
            event_calls[event_idx][0](*event_calls[event_idx][1:])

            ### DEBUG
            if self.debug_mode:
                self.run_logger.debug("Post-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    assert len(lineage.distribution_vector.presences()) > 0


            ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
            ntips = len(self.phylogeny.current_lineages)
            if self.gsa_termination_num_tips and ntips_in_focal_areas >= self.gsa_termination_num_tips:
                raise NotImplementedError
            elif self.gsa_termination_num_tips and ntips_in_focal_areas == self.target_num_tips:
                raise NotImplementedError
            elif self.target_num_tips and ntips_in_focal_areas >= self.target_num_tips:
                if self.focal_areas_tree_log is not None:
                    focal_area_tree = self.phylogeny.extract_focal_area_tree()
                    n = len(focal_area_tree.seed_node._child_nodes)
                    if n < 2:
                        raise FailedSimulationException("Insufficient lineages in focal area: {}".format(n))
                    self.write_tree(
                            out=self.focal_areas_tree_log,
                            tree=focal_area_tree,
                            focal_areas_only_labeling=True,
                            )
                if self.all_areas_tree_log is not None:
                    self.write_tree(
                            out=self.all_areas_tree_log,
                            tree=self.phylogeny,
                            focal_areas_only_labeling=False,
                            )
                break

    def schedule_events(self):
        event_calls = []
        event_rates = []

        # if self.debug_mode:
        #     num_current_lineages = len(self.phylogeny.current_lineages)
        #     self.run_logger.debug("Scheduling events for {} current lineages".format(
        #         num_current_lineages))

        for lineage in self.phylogeny.iterate_current_lineages():

            # if self.debug_mode:
            #     self.run_logger.debug("Scheduling events for lineage {}".format(lineage))

            # speciation
            event_calls.append( (self.phylogeny.split_lineage, lineage) )
            event_rates.append(self.lineage_birth_probability_function(lineage))
            # extinction
            extinction_prob = self.lineage_death_probability_function(lineage)
            if extinction_prob:
                event_calls.append( (self.phylogeny.extirpate_lineage, lineage) )
                event_rates.append(extinction_prob)
            # trait evolution
            for trait_idx, trait_state in enumerate(lineage.traits_vector):
                for state_idx in range(self.trait_types[trait_idx].nstates):
                    if state_idx == trait_idx:
                        continue
                    trait_transition_rate = self.trait_types[trait_idx].transition_rate_matrix[trait_state][state_idx]
                    if trait_transition_rate:
                        event_calls.append( (self.phylogeny.evolve_trait, lineage, trait_idx, state_idx) )
                        event_rates.append(trait_transition_rate)
            # dispersal
            for area_idx, occurs in enumerate(lineage.distribution_vector):
                if not occurs:
                    continue
                for dest_idx in self.geography.area_indexes:
                    if dest_idx == area_idx:
                        continue
                    dispersal_weight = self.geography.dispersal_weights[area_idx][dest_idx]
                    dispersal_prob = self.lineage_dispersal_probability_function(lineage)
                    dispersal_rate = dispersal_weight * dispersal_prob
                    if dispersal_rate:
                        event_calls.append( (self.phylogeny.disperse_lineage, lineage, dest_idx) )
                        event_rates.append(dispersal_rate)
        sum_of_event_rates = sum(event_rates)
        return event_calls, event_rates, sum_of_event_rates

    def write_tree(self,
            out,
            tree,
            focal_areas_only_labeling,
            ):
        if focal_areas_only_labeling:
            labelf = lambda x: utility.encode_lineage(x, exclude_areas=self.geography.supplemental_area_indexes)
        else:
            labelf = lambda x: utility.encode_lineage(x, exclude_areas=None)
        tree.write_to_stream(
                out,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_func=labelf,
                suppress_internal_node_labels=self.is_suppress_internal_node_labels,
                )

    def debug_compose_tree(self, tree):
        s = tree.as_string(
                "newick",
                node_label_compose_func=utility.encode_lineage,
                suppress_edge_lengths=True)
        return s.replace("\n", "")

def repeat_run(
        nreps,
        model_d,
        config_d,
        random_seed=None,
        output_prefix=None,
        stderr_logging_level="info",
        file_logging_level="debug",
        maximum_num_reruns_per_replicates=1000):
    """
    Executes multiple runs of the Archipelago simulator under identical
    parameters to produce the specified number of replicates, discarding failed
    runs.

    Parameters
    ----------
    nreps : integer
        Number of replicates to produce.
    config_d : dict
        Simulator configuration parameters as keyword-value pairs. To be
        re-used for each replicate.
    model_d : dict
        Simulator model parameters as keyword-value pairs. To be re-used for
        each replicate.
    random_seed : integer
        Random seed to be used (for single random number generator across all
        replicates).
    stderr_logging_level : string or None
        Message level threshold for screen logs; if 'none' or `None`, screen
        logs will be supprsed.
    file_logging_level : string or None
        Message level threshold for file logs; if 'none' or `None`, file
        logs will be supprsed.
    maximum_num_reruns_per_replicates : int
        A failed replicate (due to e.g., total extinction of all taxa) will be
        re-run. This limits the number of re-runs.
    """
    if output_prefix is None:
        output_prefix = config_d.pop("output_prefix", "archipelago")
    config_d["output_prefix"] = output_prefix
    if stderr_logging_level is None or stderr_logging_level.lower() == "none":
        log_to_stderr = False
    else:
        log_to_stderr = True
    if file_logging_level is None or file_logging_level.lower() == "none":
        log_to_file = False
    else:
        log_to_file = True
    if "run_logger" not in config_d:
        config_d["run_logger"] = utility.RunLogger(
                name="archipelago",
                log_to_stderr=log_to_stderr,
                stderr_logging_level=stderr_logging_level,
                log_to_file=log_to_file,
                log_path=output_prefix + ".log",
                file_logging_level=file_logging_level,
                )
    run_logger = config_d["run_logger"]
    run_logger.info("||ARCHIPELAGO-META|| Starting: {}".format(archipelago.description()))
    if "rng" not in config_d:
        if random_seed is None:
            random_seed = config_d.pop("random_seed", None)
            if random_seed is None:
                random_seed = random.randint(0, sys.maxsize)
        run_logger.info("||ARCHIPELAGO-META|| Initializing with random seed: {}".format(random_seed))
        config_d["rng"] = random.Random(random_seed)
    else:
        run_logger.info("||ARCHIPELAGO-META|| Using existing RNG: {}".format(config_d["rng"]))
    header_written = False
    current_rep = 0
    while current_rep < nreps:
        simulation_name="Run{}".format((current_rep+1))
        run_output_prefix = "{}.R{:04d}".format(output_prefix, current_rep+1)
        run_logger.info("||ARCHIPELAGO-META|| Run {} of {}: starting".format(current_rep+1, nreps))
        num_reruns = 0
        while True:
            archipelago_simulator = ArchipelagoSimulator(
                model_d=model_d,
                config_d=config_d,
                verbose_setup=num_reruns == 0)
            try:
                archipelago_simulator.run()
                run_logger.system = None
            except TotalExtinctionException as e:
                run_logger.system = None
                run_logger.info("||ARCHIPELAGO-META|| Replicate {} of {}: total extinction of all lineages before termination condition".format(current_rep+1, nreps, num_reruns))
                num_reruns += 1
                if num_reruns > maximum_num_reruns_per_replicates:
                    run_logger.info("||ARCHIPELAGO-META|| Replicate {} of {}: maximum number of re-runs exceeded: aborting".format(current_rep+1, nreps, num_reruns))
                else:
                    run_logger.info("||ARCHIPELAGO-META|| Replicate {} of {}: re-running replicate (number of re-runs: {})".format(current_rep+1, nreps, num_reruns))
            else:
                run_logger.system = None
                run_logger.info("||ARCHIPELAGO-META|| Replicate {} of {}: completed to termination condition".format(current_rep+1, nreps, num_reruns))
                num_reruns = 0
                break
        current_rep += 1
