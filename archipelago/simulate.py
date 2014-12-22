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
from distutils.util import strtobool

import dendropy
from dendropy.utility import textprocessing

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
        """
        self._states = [2] * nchar

    @property
    def nchar(self):
        return len(self)

    def __len__(self):
        return self._nchar

    def __getitem__(self, trait_index):
        return self._states[trait_index]

    def __setitem__(self, trait_index, v):
        self._states[trait_index] = v

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
        if values is None:
            values = [0] * len(self.trait_types)
        for idx, v in enumerate(values):
            if v is None:
                v = 0
            assert v >= 0 and v < self.trait_types[idx]
            s[idx] = v
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
        return iter(self.geography)

    def __getitem__(self, idx):
        return self.geography[idx]

    def parse_definition(self,
            areas,
            run_logger,
            verbose=True):
        self.geography = []
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
            self.geography.append(area)
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
        if len(self.geography) < 1:
            raise ValueError("No areas defined")
        if verbose:
            run_logger.info("[GEOGRAPHY] {} areas defined: {}".format(
                len(self.geography),
                ", ".join("'{}'".format(a.label) for a in self.geography),
                ))
        self.dispersal_weights = []
        total_dispersal_weight = 0.0
        for a1_idx, area1 in enumerate(self.geography):
            self.dispersal_weights.append([])
            for a2_idx, area2 in enumerate(self.geography):
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
            for a1_idx, area1 in enumerate(self.geography):
                for a2_idx, area2 in enumerate(self.geography):
                    self.dispersal_weights[a1_idx][a2_idx] /= total_dispersal_weight
        if verbose:
            if self.normalize_dispersal_weights:
                weight_type = "Normalized dispersal"
            else:
                weight_type = "Dispersal"
            for a1, area1 in enumerate(self.geography):
                run_logger.info("[GEOGRAPHY] {} weights from area '{}': {}".format(weight_type, area1.label, self.dispersal_weights[a1]))

        # instead of recalculating every time
        self.area_nstates = [2 for i in self.geography]

    def new_occurrence_vector(self, incidences=None):
        s = StatesVector(
                nchar=len(self.geography),
                nstates=self.area_nstates,
                )
        if incidences is None:
            incidences = [0] * len(self.geography)
        for idx, v in enumerate(incidences):
            if v is None:
                v = 0
            assert v >= 0 and v < 2
            s[idx] = v
        return s

class Lineage(dendropy.Node):

    def __init__(self,
            index,
            area_occurrences=None,
            trait_suite=None,
            ):
        dendropy.Node.__init__(self)
        self.index = index
        self.area_occurrences = area_occurrences
        self.trait_suite = trait_suite
        self.extant = True
        self.edge.length = 0

class Phylogeny(dendropy.Tree):

    class TotalExtinctionException(Exception):
        def __init__(self, *args, **kwargs):
            Exception.__init__(self, *args, **kwargs)

    class TargetNumberOfTipsException(Exception):
        def __init__(self, num_extant_tips_exception_trigger, num_extant_tips, *args, **kwargs):
            self.num_extant_tips = num_extant_tips
            self.num_extant_tips_exception_trigger = num_extant_tips_exception_trigger
            Exception.__init__(self, *args, **kwargs)

    def node_factory(cls, **kwargs):
        return Lineage(**kwargs)
    node_factory = classmethod(node_factory)

    def __init__(self, system):
        self.system = system
        self.lineage_indexer = IndexGenerator(0)
        seed_node = self.node_factory(
                index=next(self.lineage_indexer),
                area_occurrences=self.system.geography.new_occurrence_vector(),
                trait_suite=self.system.trait_types.new_traits_vector(),
                )
        seed_node.area_occurrences[0] = 1
        dendropy.Tree.__init__(self, seed_node=seed_node)
        self.tips = set([self.seed_node])

    def event(self):

        # set up event phenomenology
        event_calls = []
        event_rates = []
        for lineage in self.tips:
            # speciation
            event_calls.append( (self.split_lineage, lineage) )
            event_rates.append(self.system.lineage_speciation_probability_function(lineage))
            # extinction
            extinction_prob = self.system.lineage_death_probability_function(lineage)
            if extinction_prob:
                event_calls.append( (self.extinguish_lineage, lineage) )
                event_rates.append(extinction_prob)
            # trait evolution
            for trait_idx, trait_state in enumerate(lineage.trait_suite):
                for state_idx in range(self.system.trait_types[trait_idx].nstates):
                    if state_idx == trait_idx:
                        continue
                    trait_transition_rate = self.system.trait_types[trait_idx].transition_rate_matrix[trait_state][state_idx]
                    if trait_transition_rate:
                        event_calls.append( (self.evolve_trait, lineage, trait_idx, state_idx) )
                        event_rates.append(trait_transition_rate)
            # dispersal
            for area_idx, occurrence in enumerate(lineage.area_occurrences):
                for dest_idx in self.system.geography.area_indexes:
                    if dest_idx == area_idx:
                        continue
                    dispersal_weight = self.system.geography.dispersal_weights[area_idx][dest_idx]
                    dispersal_prob = self.system.lineage_dispersal_probability_function(lineage)
                    dispersal_rate = dispersal_weight * dispersal_prob
                    if dispersal_rate:
                        event_calls.append( (self.disperse_lineage, lineage, dest_idx) )
                        event_rates.append(dispersal_rate)
        sum_of_event_rates = sum(event_rates)

        # execute event
        time_till_event = self.system.rng.expovariate(sum_of_event_rates)
        self.system.current_time += time_till_event
        for lineage in self.tips:
            lineage.edge.length += time_till_event
        event_idx = weighted_index_choice(event_rates, self.system.rng)
        event_calls[event_idx][0](*event_calls[event_idx][1:])

    def split_lineage(self, lineage):
        lineage.extant = False
        self.tips.remove(lineage)
        c1 = self.node_factory(
                index=next(self.lineage_indexer),
                area_occurrences=self.system.geography.new_occurrence_vector(),
                trait_suite=self.system.trait_types.new_traits_vector(),
                )
        c2 = self.node_factory(
                index=next(self.lineage_indexer),
                area_occurrences=self.system.geography.new_occurrence_vector(),
                trait_suite=self.system.trait_types.new_traits_vector(),
                )
        # to do:
        # - handle sympatrix, allopatric, para-allopatric speciation
        lineage.add_child(c1)
        lineage.add_child(c2)
        self.tips.add(c1)
        self.tips.add(c2)

    def extinguish_lineage(self, lineage):
        pass

    def evolve_trait(self, lineage, trait_idx, state_idx):
        lineage.trait_suite[trait_idx] = state_idx

    def disperse_lineage(self, lineage, dest_area_idx):
        lineage.area_occurrences[dest_area_idx] = 1

class ArchipelagoSimulator(object):

    @staticmethod
    def get_fixed_value_function(v, description):
        f = lambda x: v
        f.__doc__ = description
        return f

    def __init__(self,
            config_d=None,
            model_d=None,
            verbose_setup=True
            ):

        # system globals
        self.current_gen = 0
        self.phylogeny = None

        # run configuration
        self.output_prefix = None
        self.run_logger = None
        self.name = None
        self.tree_log = None
        self.general_stats_log = None
        self.rng = None
        self.track_extinct_lineages = None
        self.debug_mode = None
        self.log_frequency = None
        self.report_frequency = None

        # configure
        if config_d is None:
            config_d = {}
        self.configure_simulator(config_d, verbose=verbose_setup)
        if model_d is None:
            model_d = {}
        self.set_model(model_d, verbose=verbose_setup)

        # start
        self.current_time = 0.0
        self.phylogeny = Phylogeny(self)

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
                    log_path=self.output_prefix + ".log")
        self.run_logger.system = self

        if verbose:
            self.run_logger.info("Configuring simulation '{}'".format(self.name))

        self.tree_log = config_d.pop("tree_log", None)
        if self.tree_log is None:
            self.tree_log = open(self.output_prefix + ".trees", "w")
        if verbose:
            self.run_logger.info("Tree log filepath: {}".format(self.tree_log.name))

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

        self.track_extinct_lineages = config_d.pop("track_extinct_lineages", False)
        if verbose:
            if self.track_extinct_lineages:
                self.run_logger.info("Extinct lineages will be tracked: lineages will be retained in the tree even if they are extirpated from all habitats in all islands")
            else:
                self.run_logger.info("Extinct lineages will not be tracked: lineages will be pruned from the tree if they are extirpated from all habitats in all islands")
        self.debug_mode = config_d.pop("debug_mode", False)
        if verbose and self.debug_mode:
            self.run_logger.info("Running in DEBUG mode")

        self.log_frequency = config_d.pop("log_frequency", 1000)
        self.report_frequency = config_d.pop("report_frequency", None)
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
        if "lineage_speciation_probability_function" in diversification_d:
            self.lineage_speciation_probability_function = diversification_d.pop("lineage_speciation_probability_function")
        else:
            self.lineage_speciation_probability_function = ArchipelagoSimulator.get_fixed_value_function(
                    0.01,
                    "Fixed speciation probability: {}".format(0.01)
            )
        if verbose:
            desc = getattr(self.lineage_speciation_probability_function, "__doc__", None)
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
