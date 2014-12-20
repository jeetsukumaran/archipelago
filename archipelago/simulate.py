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

    def reset(self, start=None):
        if start is None:
            start = self.start
        self.index = start

class StatesVector(object):
    """
    A vector in which each element is an integer represents the state of a
    character.

    E.g.,

        [1,0,1,2]

    is a 4-character vector, where character 0 is in state 1, character 1 is in
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
            The number of characters to be tracked.
        nstates : list of integers
            The number of states for each character. If not specified, defaults
            to binary (i.e., 2 states, 0 and 1). If specifed, must be a list of
            length `nchar`, with each element in the list being integer > 0.
        """
        self._states = [0] * nchar
        self._nchar = len(self.states)

    @property
    def nchar(self):
        return len(self)

    def __len__(self):
        return self._nchar

class Area(object):

    def __init__(self,
            index=None,
            label=None,
            is_supplemental=False,
            relative_diversity=None,
            ):
        self.index = index
        self.label = label
        self.is_supplemental = is_supplemental
        self.relative_diversity = relative_diversity

class Geography(object):

    def __init__(self):
        self.normalize_dispersal_weights = True

    def parse_definition(self,
            area_definitions,
            run_logger,
            verbose=True):
        self.areas = []
        self.area_label_index_map = collections.OrderedDict()
        self.focal_area_indexes = []
        self.supplemental_area_indexes = []
        for area_idx, area_d in enumerate(area_definitions):
            area = Area(
                index=area_idx,
                label=str(area_d.pop("label", area_idx)),
                relative_diversity=area_d.pop("relative_diversity", 1.0),
                is_supplemental=area_d.pop("is_supplemental", False)
            )
            area._dispersal_weights_d = area_d.pop("dispersal_weights", {}) # delay processing until all areas have been defined
            self.areas.append(area)
            self.area_label_index_map[area.label] = area.index
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

class Lineage(dendropy.Node):

    def __init__(self,
            index,
            geography=None,
            traits=None,
            ):
        dendropy.Node.__init__(self)
        self.index = index
        self.geography = geography
        self.traits = traits

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

    def __init__(self):
        self.lineage_indexer = IndexGenerator(0)

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
        self.traits = collections.OrderedDict()
        if "traits" in model_d:
            trait_d_list = model_d.pop("traits")
            for trait_idx, trait_d in enumerate(trait_d_list):
                label = trait_d.pop("label", "trait{}".format(trait_idx))
                self.traits[label] = {}
                self.traits[label]["index"] = trait_idx
                self.traits[label]["nstates"] = trait_d.pop("nstates", 2)
                self.traits[label]["transition_probability"] = trait_d.pop("nstates", None)
                if self.traits[label]["transition_probability"] is None:
                    self.traits[label]["transition_probability"] = ArchipelagoSimulator.get_fixed_value_function(
                            0.0,
                            "Fixed uniform transition probability: {}".format(0.01)
                    )
                if self.traits[label]["transition_probability"].__doc__ is None:
                    self.traits[label]["transition_probability"].__doc__ = "undescribed transition probabilities"
                if verbose:
                    self.run_logger.info("[ECOLOGY] Configuring trait {idx}: '{label}' ({nstates} states, {trans_prob})".format(
                        idx=trait_idx,
                        label=label,
                        nstates=self.traits[label]["nstates"],
                        trans_prob=self.traits[label]["transition_probability"].__doc__,
                        ))
            if trait_d:
                raise TypeError("Unsupported trait keywords: {}".format(trait_d))
            if verbose:
                self.run_logger.info("[ECOLOGY] Total number of traits defined: {}".format(len(self.traits)))
        else:
            if verbose:
                self.run_logger.info("No traits defined")

        # Diversification submodel
        diversification_d = model_d.pop("diversification", {})
        if "lineage_birth_probability" in diversification_d:
            self.lineage_birth_probability = diversification_d.pop("lineage_birth_probability")
        else:
            self.lineage_birth_probability = ArchipelagoSimulator.get_fixed_value_function(
                    0.01,
                    "Fixed birth probability: {}".format(0.01)
            )
        if verbose:
            desc = getattr(self.lineage_birth_probability, "__doc__", None)
            if desc is None:
                desc = "(no description available)"
            self.run_logger.info("[DIVERSIFICATION] Setting lineage birth probability function: {}".format(desc,))

        if "lineage_death_probability" in diversification_d:
            self.lineage_death_probability = diversification_d.pop("lineage_death_probability")
        else:
            self.lineage_death_probability = ArchipelagoSimulator.get_fixed_value_function(
                    0.01,
                    "Fixed death probability: {}".format(0.01)
            )
        if verbose:
            desc = getattr(self.lineage_death_probability, "__doc__", None)
            if desc is None:
                desc = "(no description available)"
            self.run_logger.info("[DIVERSIFICATION] Setting lineage death probability function: {}".format( desc,))

        # Dispersal submodel

        # Trait evolution submodel

        if model_d:
            raise TypeError("Unsupported model keywords: {}".format(model_d))
