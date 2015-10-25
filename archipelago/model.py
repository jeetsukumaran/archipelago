#! /usr/bin/env python

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import os
import random
import collections
import argparse
import pprint
import copy
import json
from distutils.util import strtobool
import dendropy

from archipelago import utility
from archipelago import error

def weighted_choice(seq, weights, rng):
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
    sow = sum(weights)
    if len(weights) == len(seq) - 1:
        weights.append(1 - sow)
    return seq[weighted_index_choice(weights, sow, rng)]

def weighted_index_choice(weights, sum_of_weights, rng):
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
    rnd = rng.uniform(0, 1) * sum_of_weights
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

class ArchipelagoModel(object):

    _TRAITS_SEPARATOR = "."
    _LABEL_COMPONENTS_SEPARATOR = "^"
    _NULL_TRAITS = "NA"

    @classmethod
    def create(
            cls,
            model_definition_source,
            model_definition_type,
            interpolate_missing_model_values=False,
            run_logger=None,
            ):
        """
        Create and return a model under which to run a simulation.

        Parameters
        ----------
        model_definition_source : object
            See 'model_definition_type' argument for values this can take.
        model_definition_type : str
            Whether 'model_definition_source' is:

                - 'python-dict' : a Python dictionary defining the model.
                - 'python-dict-str' : a string providing a Python dictionary
                  defining the model.
                - 'python-dict-filepath' : a path to a Python file to be evaluated;
                  the file should be a valid Python script containing nothing but a
                  dictionary defining the model.
                - 'json-filepath': a path to a JSON file containing a dictionary
                  defining the model.

        Returns
        -------
        m : ArchipelagoModel
            A fully-specified Archipelago model.

        """
        if model_definition_type == "python-dict-filepath":
            src = open(model_definition_source, "r")
            model_definition = eval(src.read())
        elif model_definition_type == "python-dict-str":
            model_definition = eval(model_definition_source)
        elif model_definition_type == "python-dict":
            model_definition = model_definition_source
        elif model_definition_type == "json-filepath":
            src = open(model_definition_source, "r")
            model_definition = json.load(src)
        else:
            raise ValueError("Unrecognized model definition type: '{}'".format(model_definition_type))
        return cls.from_definition_dict(
                model_definition=model_definition,
                run_logger=run_logger,
                interpolate_missing_model_values=interpolate_missing_model_values)

    @classmethod
    def from_definition_dict(cls,
            model_definition,
            interpolate_missing_model_values=False,
            run_logger=None):
        archipelago_model = cls()
        archipelago_model.parse_definition(
                model_definition=model_definition,
                interpolate_missing_model_values=interpolate_missing_model_values,
                run_logger=run_logger,
        )
        return archipelago_model

    @staticmethod
    def compose_encoded_label(
            lineage,
            excluded_area_indexes=None):
        if lineage.traits_vector:
            traits_v = ArchipelagoModel._TRAITS_SEPARATOR.join(str(i) for i in lineage.traits_vector)
        else:
            traits_v = ArchipelagoModel._NULL_TRAITS
        if excluded_area_indexes is None:
            areas_v = "".join(str(i) for i in lineage.distribution_vector)
        else:
            areas_v = "".join(str(i) for idx, i in enumerate(lineage.distribution_vector) if idx not in excluded_area_indexes)
        encoding = "s{lineage_index}{sep}{traits_v}{sep}{areas_v}".format(
                lineage_index=lineage.index,
                traits_v=traits_v,
                areas_v=areas_v,
                sep=ArchipelagoModel._LABEL_COMPONENTS_SEPARATOR)
        return encoding

    @staticmethod
    def decode_label(label):
        parts = label.split(ArchipelagoModel._LABEL_COMPONENTS_SEPARATOR)
        traits_string = parts[1]
        if not traits_string or traits_string == ArchipelagoModel._NULL_TRAITS:
            traits_vector = StatesVector(nchar=0)
        else:
            traits_string_parts = traits_string.split(ArchipelagoModel._TRAITS_SEPARATOR)
            traits_vector = StatesVector(
                    nchar=len(traits_string_parts),
                    # The trait states need to be an integer if
                    # archipelago-summarize.py coerces the user input to
                    # integers
                    # values=[int(i) for i in traits_string_parts],
                    # The reason we do NOT want it parsed to an integer value
                    # is to allow null traits 'NA', 'null', etc.
                    values=[i for i in traits_string_parts],
                    )
        distribution_string = parts[2]
        distribution_vector = DistributionVector(
                num_areas=len(distribution_string),
                values=[int(i) for i in distribution_string],)
        return traits_vector, distribution_vector

    @staticmethod
    def set_lineage_data(
            tree,
            leaf_nodes_only=False,
            lineage_data_source="node",
            traits_filepath=None,
            areas_filepath=None,
            ):
        if lineage_data_source == "node":
            _decode = lambda x: ArchipelagoModel.decode_label(x.label)
        elif lineage_data_source == "taxon":
            _decode = lambda x: ArchipelagoModel.decode_label(x.taxon.label)
        else:
            raise ValueError("'lineage_data_source' must be 'node' or 'taxon'")
        for nd in tree:
            if (not leaf_nodes_only or not nd._child_nodes) and (lineage_data_source == "node" or nd.taxon is not None):
                traits_vector, distribution_vector = _decode(nd)
                nd.traits_vector = traits_vector
                nd.distribution_vector = distribution_vector
            else:
                nd.traits_vector = None
                nd.distribution_vector = None

    def __init__(self):
        pass

    def parse_definition(self,
            model_definition,
            run_logger=None,
            interpolate_missing_model_values=True):

        # initialize
        if model_definition is None:
            model_definition = {}
        else:
            model_definition = dict(model_definition)

        # Geography
        if "areas" not in model_definition:
            if interpolate_missing_model_values:
                model_definition["areas"] = [
                        {'is_supplemental': False, 'label': 'a1'},
                        {'is_supplemental': False, 'label': 'a2'},
                        {'is_supplemental': False, 'label': 'a3'},
                        {'is_supplemental': False, 'label': 'a4'},
                        {'is_supplemental': True, 'label': 's1'}
                ]
            else:
                raise ValueError("No areas defined")
        self.geography = Geography()
        self.geography.parse_definition(
                copy.deepcopy(model_definition.pop("areas", [])),
                run_logger=run_logger)

        # Ecology
        self.trait_types = TraitTypes()
        self.trait_types.parse_definition(
                copy.deepcopy(model_definition.pop("traits", [])),
                run_logger=run_logger)

        # Diversification
        ## speciation
        diversification_d = dict(model_definition.pop("diversification", {}))
        if "lineage_birth_rate" in diversification_d:
            self.lineage_birth_rate_function = RateFunction.from_definition_dict(diversification_d.pop("lineage_birth_rate"), self.trait_types)
        else:
            self.lineage_birth_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.01",
                    description="fixed: 0.01",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Setting lineage birth rate function: {desc}".format(
                desc=self.lineage_birth_rate_function.description,))
        if diversification_d:
            raise TypeError("Unsupported diversification model keywords: {}".format(diversification_d))

        # Dispersal submodel
        anagenetic_range_evolution_d = dict(model_definition.pop("anagenetic_range_evolution", {}))
        if "global_area_gain_rate" not in anagenetic_range_evolution_d and "mean_area_gain_rate" not in anagenetic_range_evolution_d:
            if interpolate_missing_model_values:
                anagenetic_range_evolution_d["global_area_gain_rate"] = 0.01
            else:
                raise TypeError("Exactly one of 'global_area_gain_rate' or 'mean_area_gain_rate' must be specified")
        if "global_area_gain_rate" in anagenetic_range_evolution_d and "mean_area_gain_rate" in anagenetic_range_evolution_d:
            raise TypeError("No more than one of 'global_area_gain_rate' or 'mean_area_gain_rate' can be specified")
        elif "global_area_gain_rate" in anagenetic_range_evolution_d:
            self.global_area_gain_rate = float(anagenetic_range_evolution_d.pop("global_area_gain_rate"))
            self.mean_area_gain_rate = None
            if run_logger is not None:
                run_logger.info("(ANAGENETIC RANGE EVOLUTION) Global dispersal rate is: {}".format(self.global_area_gain_rate))
            self.geography.set_global_area_gain_rate(self.global_area_gain_rate)
        else:
            self.mean_area_gain_rate = float(anagenetic_range_evolution_d.pop("mean_area_gain_rate"))
            self.global_area_gain_rate = None
            run_logger.info("(ANAGENETIC RANGE EVOLUTION) Mean dispersal rate is: {}".format(self.mean_area_gain_rate))
            self.geography.set_mean_area_gain_rate(self.mean_area_gain_rate)
        if run_logger is not None:
            for a1, area1 in enumerate(self.geography.areas):
                run_logger.info("(ANAGENETIC RANGE EVOLUTION) Effective rate of area gain from area '{}': {}".format(area1.label, self.geography.effective_area_gain_rates[a1]))

        if "lineage_area_gain_weight" in anagenetic_range_evolution_d:
            self.lineage_area_gain_weight_function = RateFunction.from_definition_dict(anagenetic_range_evolution_d.pop("lineage_area_gain_weight"), self.trait_types)
        else:
            self.lineage_area_gain_weight_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.01",
                    description="fixed: 0.01",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC RANGE EVOLUTION) Setting lineage area gain rate function: {desc}".format(
                desc=self.lineage_area_gain_weight_function.description,))

        ## extinction
        self.treat_area_loss_rate_as_lineage_death_rate = strtobool(str(anagenetic_range_evolution_d.pop("treat_area_loss_rate_as_lineage_death_rate", 0)))
        if "lineage_area_loss_rate" in anagenetic_range_evolution_d:
            self.lineage_area_loss_rate_function = RateFunction.from_definition_dict(anagenetic_range_evolution_d.pop("lineage_area_loss_rate"), self.trait_types)
        else:
            self.lineage_area_loss_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.0",
                    description="fixed: 0.0",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            # run_logger.info("(ANAGENETIC RANGE EVOLUTION) Setting lineage area loss (= {is_global}) rate function: {desc}".format(
            #     is_global="global extinction" if self.treat_area_loss_rate_as_lineage_death_rate else "local extirpation",
            #     desc=self.lineage_area_loss_rate_function.description,
            #     ))
            if self.treat_area_loss_rate_as_lineage_death_rate:
                run_logger.info("(ANAGENETIC RANGE EVOLUTION) Setting lineage global extinction rate function: {desc}".format(
                    desc=self.lineage_area_loss_rate_function.description,
                    ))
            else:
                run_logger.info("(ANAGENETIC RANGE EVOLUTION) Setting lineage area loss rate function: {desc}".format(
                    desc=self.lineage_area_loss_rate_function.description,
                    ))

        if anagenetic_range_evolution_d:
            raise TypeError("Unsupported dispersal model keywords in anagenetic dispersal submodel: {}".format(anagenetic_range_evolution_d))

        # Cladogenetic range inheritance submodel
        cladogenesis_d = dict(model_definition.pop("cladogenetic_range_evolution", {}))
        self.cladogenesis_sympatric_subset_speciation_weight = float(cladogenesis_d.pop("sympatric_subset_speciation_weight", 1.0))
        self.cladogenesis_single_area_vicariance_speciation_weight = float(cladogenesis_d.pop("single_area_vicariance_speciation_weight", 1.0))
        self.cladogenesis_widespread_vicariance_speciation_weight = float(cladogenesis_d.pop("widespread_vicariance_speciation_weight", 1.0))
        self.cladogenesis_founder_event_speciation_weight = float(cladogenesis_d.pop("founder_event_speciation_weight", 0.0))
        if cladogenesis_d:
            raise TypeError("Unsupported keywords in cladogenesis submodel: {}".format(cladogenesis_d))
        if run_logger is not None:
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of sympatric subset speciation mode: {}".format(self.cladogenesis_sympatric_subset_speciation_weight))
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of single area vicariance speciation mode: {}".format(self.cladogenesis_single_area_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of widespread vicariance speciation mode: {}".format(self.cladogenesis_widespread_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of founder event speciation ('jump dispersal') mode: {} (note that the effective weight of this event for each lineage is actually the product of this and the lineage-specific dispersal weight)".format(self.cladogenesis_founder_event_speciation_weight))

        termination_conditions_d = dict(model_definition.pop("termination_conditions", {}))
        self.target_focal_area_lineages = termination_conditions_d.pop("target_focal_area_lineages", None)
        self.gsa_termination_focal_area_lineages = termination_conditions_d.pop("gsa_termination_focal_area_lineages", None)
        self.max_time = termination_conditions_d.pop("max_time", None)
        if termination_conditions_d:
            raise TypeError("Unsupported termination condition model keywords: {}".format(termination_conditions_d))
        if self.gsa_termination_focal_area_lineages and not self.target_focal_area_lineages:
            raise ValueError("Cannot specify 'gsa_termination_focal_area_lineages' without specifying 'target_focal_area_lineages'")
        if self.target_focal_area_lineages is None and self.max_time is None:
            if run_logger is not None:
                run_logger.info("Termination conditions not specified: default termination conditions applied")
            self.target_focal_area_lineages = 50
        if not self.target_focal_area_lineages and self.max_time:
            desc = "Simulation will terminate at time t = {}".format(self.max_time)
        elif self.target_focal_area_lineages and not self.gsa_termination_focal_area_lineages and not self.max_time:
            desc = "Simulation will terminate when there are {} lineages in focal areas (no time limit)".format(self.target_focal_area_lineages)
        elif self.target_focal_area_lineages and not self.gsa_termination_focal_area_lineages and self.max_time:
            desc = "Simulation will terminate at time t = {} or when there are {} lineages in focal areas".format(self.max_time, self.target_focal_area_lineages)
        elif self.target_focal_area_lineages and self.gsa_termination_focal_area_lineages and not self.max_time:
            desc = "Simulation will terminate when there are {} lineages in focal areas (with the phylogeny sampled at a random slice of time when there were {} extant lineages in the focal areas)".format(self.gsa_termination_focal_area_lineages, self.target_focal_area_lineages)
        elif self.target_focal_area_lineages and self.gsa_termination_focal_area_lineages and self.max_time:
            desc = "Simulation will terminate at time t = {} or when there are {} lineages in focal areas (with the phylogeny sampled at a random slice of time when there were {} extant lineages in the focal areas)".format(self.max_time, self.gsa_termination_focal_area_lineages, self.target_focal_area_lineages)
        elif not self.target_focal_area_lineages and not self.max_time:
            raise ValueError("Unspecified termination condition")
        else:
            raise ValueError("Unsupported termination condition(s)")
        if run_logger is not None:
            run_logger.info(desc)

        if model_definition:
            raise TypeError("Unsupported model keywords: {}".format(model_definition))

    def encode_lineage(self,
            lineage,
            set_label=False,
            add_annotation=False,
            exclude_supplemental_areas=False,
            ):
        encoded_label = ArchipelagoModel.compose_encoded_label(
                lineage=lineage,
                excluded_area_indexes=self.geography.supplemental_area_indexes if exclude_supplemental_areas else None,
                )
        if set_label:
            lineage.label =encoded_label
        if add_annotation:
            lineage.annotations.drop()
            lineage.annotations.add_new("traits_v", traits_v)
            lineage.annotations.add_new("distribution", areas_v)
            for trait_idx, trait in enumerate(self.trait_types):
                lineage.annotations.add_new(trait.label, lineage.traits_vector[trait_idx])
            area_list = []
            for area_idx, area in enumerate(self.geography.areas):
                if exclude_supplemental_areas and area.is_supplemental:
                    continue
                if lineage.distribution_vector[area_idx] == 1:
                    area_list.append(area.label)
            lineage.annotations.add_new("areas", area_list)
        return encoded_label

    def write_model(self, out):
        model_definition = collections.OrderedDict()
        model_definition["areas"] = self.geography.as_definition()
        model_definition["traits"] = self.trait_types.as_definition()
        model_definition["diversification"] = self.diversification_as_definition()
        model_definition["anagenetic_range_evolution"] = self.anagenetic_range_evolution_as_definition()
        model_definition["cladogenetic_range_evolution"] = self.cladogenetic_range_evolution_as_definition()
        model_definition["termination_conditions"] = self.termination_conditions_as_definition()
        json.dump(model_definition, out, indent=4, separators=(',', ': '))

    def diversification_as_definition(self):
        d = collections.OrderedDict()
        d["lineage_birth_rate"] = self.lineage_birth_rate_function.as_definition()
        return d

    def anagenetic_range_evolution_as_definition(self):
        d = collections.OrderedDict()
        if self.global_area_gain_rate is not None and self.mean_area_gain_rate is not None:
            raise TypeError("Both 'global_area_gain_rate' and 'mean_area_gain_rate' are populated")
        elif self.global_area_gain_rate is None and self.mean_area_gain_rate is None:
            raise TypeError("Neither 'global_area_gain_rate' and 'mean_area_gain_rate' are populated")
        elif self.global_area_gain_rate is not None:
            d["global_area_gain_rate"] = self.global_area_gain_rate
        else:
            d["mean_area_gain_rate"] = self.mean_area_gain_rate
        d["lineage_area_gain_weight"] = self.lineage_area_gain_weight_function.as_definition()
        d["lineage_area_loss_rate"] = self.lineage_area_loss_rate_function.as_definition()
        return d

    def cladogenetic_range_evolution_as_definition(self):
        d = collections.OrderedDict()
        d["sympatric_subset_speciation_weight"] = self.cladogenesis_sympatric_subset_speciation_weight
        d["single_area_vicariance_speciation_weight"] = self.cladogenesis_single_area_vicariance_speciation_weight
        d["widespread_vicariance_speciation_weight"] = self.cladogenesis_widespread_vicariance_speciation_weight
        d["founder_event_speciation_weight"] = self.cladogenesis_founder_event_speciation_weight
        return d

    def termination_conditions_as_definition(self):
        d = collections.OrderedDict()
        d["target_focal_area_lineages"] = self.target_focal_area_lineages
        d["gsa_termination_focal_area_lineages"] = self.gsa_termination_focal_area_lineages
        d["max_time"] = self.max_time
        return d

class RateFunction(object):

    @classmethod
    def from_definition_dict(cls, rate_function_d, trait_types):
        rf = cls()
        rf.parse_definition(rate_function_d, trait_types)
        return rf

    def __init__(self,
            definition_type=None,
            definition_content=None,
            description=None,
            trait_types=None,
            ):
        self.definition_type = definition_type # value, lambda, function, map
        self.definition_content = definition_content
        self.description = description
        self._compute_rate = None
        if trait_types is not None:
            self.compile_function(trait_types)

    def __call__(self, lineage):
        return self._compute_rate(lineage)

    def parse_definition(self, rate_function_d, trait_types):
        rate_function_d = dict(rate_function_d)
        self.definition_type = rate_function_d.pop("definition_type").replace("-", "_")
        self.definition_content = rate_function_d.pop("definition")
        self.description = rate_function_d.pop("description", "")
        if rate_function_d:
            raise TypeError("Unsupported function definition keywords: {}".format(rate_function_d))
        self.compile_function(trait_types)

    def compile_function(self, trait_types):
        self.definition_type = self.definition_type.replace("-", "_")
        if self.definition_type == "fixed_value":
            self.definition_content = float(self.definition_content)
            self._compute_rate = lambda lineage: self.definition_content
        elif self.definition_type == "lambda_definition":
            self._compute_rate = eval(self.definition_content)
        elif self.definition_type.startswith("trait_state_index_map"):
            parts = self.definition_type.split(":")
            if len(parts) != 2:
                raise ValueError("Expecting definition type in form of 'trait_state_index_map:<TRAIT-LABEL>' but found: '{}'".format(self.definition_type))
            trait_label = parts[1]
            if trait_label not in trait_types.trait_label_index_map:
                raise ValueError("Trait '{}' not defined: {}".format(trait_label, trait_types.trait_label_index_map.keys()))
            trait = trait_types.get_by_label(trait_label)
            rate_list = list(self.definition_content)
            if len(rate_list) != trait.nstates:
                raise ValueError("Trait '{}' has {} states, but rate mapping only provides {} values".format(
                    trait_label, trait.nstates, len(rate_list)))
            rates = [float(v) for v in rate_list]
            self._compute_rate = lambda lineage: rates[lineage.traits_vector[trait.index]]
        elif self.definition_type == "function_object":
            self._compute_rate = self.definition_content
        else:
            raise ValueError("Unrecognized function definition type: '{}'".format(self.definition_type))

    def as_definition(self):
        d = collections.OrderedDict()
        d["definition_type"] = self.definition_type
        if d["definition_type"] == "function_object":
            d["definition"] = str(self.definition_content)
        else:
            d["definition"] = self.definition_content
        d["description"] = self.description
        return d

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
        if nstates is not None:
            self._nstates = list(nstates)
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

    def __getitem__(self, idx):
        return self._states[idx]

    def __setitem__(self, idx, v):
        self._states[idx] = v

    def __repr__(self):
        return str(self._states)

class DistributionVector(StatesVector):

    def __init__(self, num_areas, values=None):
        StatesVector.__init__(self,
                nchar=num_areas,
                nstates=[2] * num_areas,
                values=values,
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

    def as_definition(self):
        d = collections.OrderedDict()
        d["label"] = self.label
        d["nstates"] = self.nstates
        d["transition_rate"] = self.transition_rate
        d["transition_weights"] = self.transition_weights
        return d

class TraitTypes(object):

    def __init__(self):
        self.normalize_transition_weights = True

    def __len__(self):
        return len(self.trait_types)

    def __iter__(self):
        return iter(self.trait_types)

    def __getitem__(self, idx):
        return self.trait_types[idx]

    def get_by_label(self, label):
        idx = self.trait_label_index_map[label]
        return self.trait_types[idx]

    def parse_definition(self,
            trait_types,
            run_logger):
        self.trait_types = []
        self.trait_label_index_map = collections.OrderedDict()
        for trait_idx, trait_d in enumerate(trait_types):
            if "label" not in trait_d:
                raise ValueError("Trait requires 'label' to be defined")
            trait_label = str(trait_d.pop("label"))
            if trait_label in self.trait_label_index_map:
                raise ValueError("Trait with label '{}' already defined".format(trait_label))
            trait = TraitType(
                index=trait_idx,
                label=trait_label,
                nstates=trait_d.pop("nstates", 2),
                transition_rate=trait_d.pop("transition_rate", 0.01),
            )
            self.trait_types.append(trait)

            trait.transition_weights = list(trait_d.pop("transition_weights", []))
            if not trait.transition_weights:
                trait.transition_weights = [[1.0 for j in range(trait.nstates)] for i in range(trait.nstates)]
            else:
                trait.transition_weights = [list(i) for i in trait.transition_weights]
            assert len(trait.transition_weights) == trait.nstates
            for a1_idx in range(trait.nstates):
                assert len(trait.transition_weights[a1_idx]) == trait.nstates
                for a2_idx in range(trait.nstates):
                    if a1_idx == a2_idx:
                        trait.transition_weights[a1_idx][a2_idx] = 0.0
                    else:
                        trait.transition_weights[a1_idx][a2_idx] = float(trait.transition_weights[a1_idx][a2_idx])
            if self.normalize_transition_weights:
                for a1_idx in range(trait.nstates):
                    normalization_factor = sum(trait.transition_weights[a1_idx])
                    if normalization_factor:
                        for a2_idx in range(trait.nstates):
                            if a1_idx == a2_idx:
                                continue
                            trait.transition_weights[a1_idx][a2_idx] /= normalization_factor
            self.trait_label_index_map[trait.label] = trait.index
            if run_logger is not None:
                run_logger.info("(TRAIT EVOLUTION) configuring trait {idx}: '{label}': {nstates} states, transition rate of {trate} with {normalized}transition weights of {tweights}".format(
                    idx=trait.index,
                    label=trait.label,
                    normalized="normalized " if self.normalize_transition_weights else "",
                    nstates=trait.nstates,
                    trate=trait.transition_rate,
                    tweights=trait.transition_weights,
                    ))
            trait.compile_matrix()
            if trait_d:
                raise TypeError("Unsupported trait model keywords: {}".format(trait_d))
        # if len(self.trait_types) < 1:
        #     raise ValueError("No traits defined")
        if run_logger is not None:
            run_logger.info("(TRAIT EVOLUTION) Total of {} traits defined{}{}".format(
                len(self.trait_types),
                ": " if self.trait_types else "",
                ", ".join("'{}'".format(a.label) for a in self.trait_types),
                ))
        self.trait_nstates = [trait.nstates for trait in self.trait_types]

    def new_traits_vector(self, values=None):
        s = StatesVector(
                nchar=len(self.trait_types),
                nstates=self.trait_nstates)
        return s

    def as_definition(self):
        traits = [t.as_definition() for t in self.trait_types]
        return traits

class Area(object):

    def __init__(self,
            index=None,
            label=None,
            is_supplemental=False,
            relative_diversity=None,
            area_gain_weights=None,
            ):
        self.index = index
        self.label = label
        self.is_supplemental = is_supplemental
        self.relative_diversity = relative_diversity
        # this is here mainly for description purposes in
        # `Area.as_definition()`; actual usage is through
        # `Geography.area_gain_weights`
        self.area_gain_weights = area_gain_weights

    def as_definition(self):
        d = collections.OrderedDict()
        # d["index"] = self.index
        d["label"] = self.label
        d["is_supplemental"] = self.is_supplemental
        d["relative_diversity"] = self.relative_diversity
        d["area_gain_weights"] = self.area_gain_weights
        return d

class Geography(object):

    def __init__(self):
        self.normalize_area_gain_weights = False

    def __iter__(self):
        return iter(self.areas)

    def __getitem__(self, idx):
        return self.areas[idx]

    def parse_definition(self,
            areas,
            run_logger):
        self.areas = []
        self.area_label_index_map = collections.OrderedDict()
        self.area_indexes = []
        self.focal_area_indexes = []
        self.supplemental_area_indexes = []
        for area_idx, area_d in enumerate(areas):
            if "label" not in area_d:
                raise ValueError("Area requires 'label' to be defined")
            area_label = str(area_d.pop("label", None))
            if area_label in self.area_label_index_map:
                raise ValueError("Area with label '{}' already defined".format(area_label))
            area = Area(
                index=area_idx,
                label=area_label,
                relative_diversity=area_d.pop("relative_diversity", 1.0),
                is_supplemental=area_d.pop("is_supplemental", False)
            )
            area._area_gain_weights_d = list(area_d.pop("area_gain_weights", [])) # delay processing until all areas have been defined
            self.areas.append(area)
            self.area_label_index_map[area.label] = area.index
            self.area_indexes.append(area.index)
            if area.is_supplemental:
                self.supplemental_area_indexes.append(area.index)
            else:
                self.focal_area_indexes.append(area.index)
            if run_logger is not None:
                run_logger.info("(GEOGRAPHY) Configuring area {idx}: '{label}' ({is_supplemental}; relative diversity: {relative_diversity})".format(
                    idx=area.index,
                    label=area.label,
                    is_supplemental="supplemental" if area.is_supplemental else "primary",
                    relative_diversity=area.relative_diversity,
                    ))
            if area_d:
                raise TypeError("Unsupported area model keywords: {}".format(area_d))
        if len(self.areas) < 1:
            raise ValueError("No areas defined")
        if len(self.focal_area_indexes) < 1:
            raise ValueError("No focal areas defined")
        if run_logger is not None:
            run_logger.info("(GEOGRAPHY) Total of {} areas defined: {}".format(
                len(self.areas),
                ", ".join("'{}'".format(a.label) for a in self.areas),
                ))
        self.area_gain_weights = []
        for a1_idx, area1 in enumerate(self.areas):
            if len(area1._area_gain_weights_d) == 0:
                area1._area_gain_weights_d = [1.0] * len(self.areas)
                area1._area_gain_weights_d[a1_idx] = 0.0
            if len(area1._area_gain_weights_d) != len(self.areas):
                raise ValueError("Expecting exactly {} elements in dispersal weight vector for area '{}', but instead found {}: {}".format(
                    len(self.areas), area1.label, len(area1._area_gain_weights_d), area1._area_gain_weights_d))
            self.area_gain_weights.append([])
            for a2_idx, area2 in enumerate(self.areas):
                if a1_idx == a2_idx:
                    d = float(area1._area_gain_weights_d[a2_idx])
                    if d != 0:
                        raise ValueError("Self-dispersal weight from area {label} to {label} must be 0.0, but instead found: {dw}".format(
                            label=area1.label, dw=d))
                    self.area_gain_weights[a1_idx].append(0.0)
                else:
                    d = float(area1._area_gain_weights_d[a2_idx])
                    self.area_gain_weights[a1_idx].append(d)
            # if area1._area_gain_weights_d:
            #     raise ValueError("Undefined dispersal targets in '{}': '{}'".format(area1.label, area1._area_gain_weights_d))
            area1.area_gain_weights = self.area_gain_weights[a1_idx]
            del area1._area_gain_weights_d
        if self.normalize_area_gain_weights:
            for a1_idx, area1 in enumerate(self.areas):
                normalization_factor = sum(self.area_gain_weights[a1_idx])
                if normalization_factor:
                    for a2_idx, area2 in enumerate(self.areas):
                        self.area_gain_weights[a1_idx][a2_idx] /= normalization_factor
                    area1.area_gain_weights = self.area_gain_weights[a1_idx]
        if run_logger is not None:
            if self.normalize_area_gain_weights:
                weight_type = "Normalized dispersal"
            else:
                weight_type = "Dispersal"
            for a1, area1 in enumerate(self.areas):
                run_logger.info("(GEOGRAPHY) {} weights from area '{}': {}".format(weight_type, area1.label, self.area_gain_weights[a1]))

        # instead of recalculating every time
        self.area_nstates = [2 for i in self.areas]

    def as_definition(self):
        areas = [a.as_definition() for a in self.areas]
        return areas

    def new_distribution_vector(self):
        s = DistributionVector(num_areas=len(self.areas))
        return s

    def set_global_area_gain_rate(self, global_area_gain_rate):
        self.effective_area_gain_rates = []
        for src_area_idx in self.area_indexes:
            self.effective_area_gain_rates.append([])
            for dest_area_idx in self.area_indexes:
                self.effective_area_gain_rates[src_area_idx].append(self.area_gain_weights[src_area_idx][dest_area_idx] * global_area_gain_rate)

    def set_mean_area_gain_rate(self, mean_area_gain_rate):
        self.effective_area_gain_rates = []
        for src_area_idx in self.area_indexes:
            self.effective_area_gain_rates.append([])
            for dest_area_idx in self.area_indexes:
                weight = self.area_gain_weights[src_area_idx][dest_area_idx]
                if weight:
                    self.effective_area_gain_rates[src_area_idx].append(mean_area_gain_rate / weight)
                else:
                    self.effective_area_gain_rates[src_area_idx].append(0.0)

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

class Phylogeny(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return Lineage(**kwargs)
    node_factory = classmethod(node_factory)

    def __init__(self, *args, **kwargs):
        if kwargs:
            self.model = kwargs.pop("model")
            self.rng = kwargs.pop("rng")
            self.debug_mode = kwargs.pop("debug_mode")
            self.run_logger = kwargs.pop("run_logger")
            self.lineage_indexer = utility.IndexGenerator(0)
            if "seed_node" not in kwargs:
                seed_node = self.node_factory(
                        index=next(self.lineage_indexer),
                        distribution_vector=self.model.geography.new_distribution_vector(),
                        traits_vector=self.model.trait_types.new_traits_vector(),
                        )
                for trait_idx in range(len(self.model.trait_types)):
                    trait_states = [i for i in range(self.model.trait_types[trait_idx].nstates)]
                    seed_node.traits_vector[trait_idx] = self.rng.choice(trait_states)
                initial_area = self.rng.randint(0, len(seed_node.distribution_vector)-1)
                seed_node.distribution_vector[initial_area] = 1
                kwargs["seed_node"] = seed_node
            dendropy.Tree.__init__(self, *args, **kwargs)
            self.current_lineages = set([self.seed_node])
        else:
            dendropy.Tree.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        memo[id(self.model)] = self.model
        memo[id(self.rng)] = None #self.rng
        memo[id(self.run_logger)] = self.run_logger
        memo[id(self.taxon_namespace)] = self.taxon_namespace
        return dendropy.Tree.__deepcopy__(self, memo)

    def iterate_current_lineages(self):
        for lineage in self.current_lineages:
            yield lineage

    def split_lineage(self, lineage):
        dist1, dist2 = self._get_daughter_distributions(lineage)
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
        if self.debug_mode:
            self.run_logger.debug("Splitting {} with distribution {} under speciation mode {} to: {} (distribution: {}) and {} (distribution: {})".format(
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

    def contract_lineage_range(self, lineage):
        if self.model.treat_area_loss_rate_as_lineage_death_rate:
            self._make_lineage_extinct_on_phylogeny(lineage)
        else:
            presences = lineage.distribution_vector.presences()
            assert len(presences) > 0
            if len(presences) == 1:
                self._make_lineage_extinct_on_phylogeny(lineage)
            else:
                lineage.distribution_vector[ self.rng.choice(presences) ] = 0

    def _get_daughter_distributions(self, lineage):
        # speciation modes
        # 0:  single-area sympatric speciation
        #     -   ancestral range copied to both daughter species
        # 1:  sympatric subset: multi-area sympatric speciation
        #     -   d1: inherits complete range
        #     -   d2: inherits single area in ancestral range
        # 2:  (single-area) vicariance
        #     -   d1: single area
        #     -   d2: all other areas
        # 3:  (multi-area) vicariance
        #     -   ancestral range divided up unequally between two daughter
        #         species
        # 4:  founder-event jump dispersal
        #     -   single new area colonized
        presences = lineage.distribution_vector.presences()
        num_presences = len(presences)
        num_areas = len(self.model.geography.area_indexes)
        if num_presences <= 1:
            speciation_mode = 0
        else:
            if num_presences < num_areas:
                lineage_area_gain_weight = self.model.lineage_area_gain_weight_function(lineage)
                fes_weight = self.model.cladogenesis_founder_event_speciation_weight * lineage_area_gain_weight
            else:
                fes_weight = 0.0
            speciation_mode_weights = [
                self.model.cladogenesis_sympatric_subset_speciation_weight,
                self.model.cladogenesis_single_area_vicariance_speciation_weight,
                self.model.cladogenesis_widespread_vicariance_speciation_weight,
                fes_weight,
            ]
            sum_of_weights = sum(speciation_mode_weights)
            speciation_mode = 1 + weighted_index_choice(
                    weights=speciation_mode_weights,
                    sum_of_weights=sum_of_weights,
                    rng=self.rng)
        if speciation_mode == 0:
            # single-area sympatric speciation
            #     -   ancestral range copied to both daughter species
            dist1 = lineage.distribution_vector.clone()
            dist2 = lineage.distribution_vector.clone()
        elif speciation_mode == 1:
            # sympatric subset: multi-area sympatric speciation
            #     -   d1: inherits complete range
            #     -   d2: inherits single area in ancestral range
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.model.geography.new_distribution_vector()
            # TODO: area diversity base speciation
            dist2[ self.rng.choice(presences) ] = 1
        elif speciation_mode == 2:
            # (single-area) allopatric vicariance
            #     -   d1: single area
            #     -   d2: all other areas
            dist1 = self.model.geography.new_distribution_vector()
            dist2 = self.model.geography.new_distribution_vector()
            self.rng.shuffle(presences)
            dist1[presences[0]] = 1
            for idx in presences[1:]:
                dist2[idx] = 1
        elif speciation_mode == 3:
            dist1 = self.model.geography.new_distribution_vector()
            dist2 = self.model.geography.new_distribution_vector()
            if num_presences == 2:
                dist1[presences[0]] = 1
                dist2[presences[1]] = 1
            else:
                n1 = self.rng.randint(1, num_presences-1)
                n2 = num_presences - n1
                if n2 == n1:
                    n1 += 1
                    n2 -= 1
                sample1 = set(self.rng.sample(presences, n1))
                for idx in self.model.geography.area_indexes:
                    if idx in sample1:
                        dist1[idx] = 1
                    else:
                        dist2[idx] = 1
        elif speciation_mode == 4:
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.model.geography.new_distribution_vector()
            absences = [idx for idx in self.model.geography.area_indexes if idx not in presences]
            dist2[ self.rng.choice(absences) ] = 1
        else:
            raise ValueError(speciation_mode)
        return dist1, dist2

    def _make_lineage_extinct_on_phylogeny(self, lineage):
        if len(self.current_lineages) == 1:
            self.total_extinction_exception("no extant lineages remaining")
        lineage.is_extant = False
        self.current_lineages.remove(lineage)
        self.prune_subtree(lineage)

    def total_extinction_exception(self, msg):
        # self.run_logger.info("Total extinction: {}".format(msg))
        raise error.TotalExtinctionException(msg)

    def evolve_trait(self, lineage, trait_idx, state_idx):
        lineage.traits_vector[trait_idx] = state_idx

    def disperse_lineage(self, lineage, dest_area_idx):
        lineage.distribution_vector[dest_area_idx] = 1

    def focal_area_lineages(self):
        focal_area_lineages = set()
        for lineage in self.iterate_current_lineages():
            for area_idx in self.model.geography.focal_area_indexes:
                if lineage.distribution_vector[area_idx] == 1:
                    focal_area_lineages.add(lineage)
                    break
        return focal_area_lineages

    def num_focal_area_lineages(self):
        count = 0
        for lineage in self.iterate_current_lineages():
            for area_idx in self.model.geography.focal_area_indexes:
                if lineage.distribution_vector[area_idx] == 1:
                    count += 1
                    break
        return count

    def extract_focal_areas_tree(self):
        # tcopy = Phylogeny(self)
        tcopy = copy.deepcopy(self)
        focal_area_lineages = tcopy.focal_area_lineages()
        if len(focal_area_lineages) < 2:
            raise error.InsufficientFocalAreaLineagesSimulationException("insufficient lineages in focal area at termination".format(len(focal_area_lineages)))
        try:
            tcopy.filter_leaf_nodes(filter_fn=lambda x: x in focal_area_lineages)
        except dendropy.SeedNodeDeletionException:
            raise error.InsufficientFocalAreaLineagesSimulationException("no extant lineages in focal area at termination".format(len(focal_area_lineages)))
        return tcopy


