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

class ArchipelagoModel(object):

    @classmethod
    def diagnose_model_file_schema(cls, filepath):
        ext = os.path.splitext(filepath)[1]
        if ext ==  ".json":
            schema = "json"
        elif ext == ".py":
            schema = "python"
        else:
            raise ValueError("Schema cannot be diagnosed from extension: '{}'".format(ext))
        return schema

    @classmethod
    def get_model_definition_from_path(cls, filepath, schema=None):
        if schema is None:
            schema = cls.diagnose_model_file_schema(filepath)
        src = open(filepath, "r")
        return cls.get_model_definition_from_file(src, schema)

    @classmethod
    def get_model_definition_from_file(cls, src, schema):
        if schema == "json":
            return json.load(src)
        elif schema == "python":
            return eval(src.read())
        else:
            raise ValueError("Unrecognized format: '{}'".format(schema))

    @classmethod
    def from_path(cls,
            filepath,
            schema=None,
            run_logger=None):
        if schema is None:
            schema = cls.diagnose_model_file_schema(filepath)
        src = open(filepath, "r")
        return cls.from_file(
                src=src,
                schema=schema,
                run_logger=run_logger)

    @classmethod
    def from_file(cls,
            src,
            schema,
            run_logger=None):
        if isinstance(src, str):
            src = open(src, "r")
        model_definition = cls.get_model_definition_from_file(
                src=src,
                schema=schema)
        return cls.from_definition(
                model_definition=model_definition,
                run_logger=run_logger)

    @classmethod
    def from_definition(cls, model_definition, run_logger):
        archipelago_model = cls()
        archipelago_model.parse_definition(
                model_definition=model_definition,
                run_logger=run_logger,
        )
        return archipelago_model

    @staticmethod
    def compose_encoded_label(
            lineage,
            excluded_area_indexes=None):
        if lineage.traits_vector:
            traits_v = "".join(str(i) for i in lineage.traits_vector)
        else:
            traits_v = "x"
        if excluded_area_indexes is None:
            areas_v = "".join(str(i) for i in lineage.distribution_vector)
        else:
            areas_v = "".join(str(i) for idx, i in enumerate(lineage.distribution_vector) if idx not in excluded_area_indexes)
        encoding = "s{}.{}.{}".format(lineage.index, traits_v, areas_v)
        return encoding

    @staticmethod
    def decode_label(label):
        parts = label.split(".")
        traits_string = parts[1]
        if not traits_string or traits_string == "x":
            traits_vector = StatesVector(nchar=0)
        else:
            traits_vector = StatesVector(
                    nchar=len(traits_string),
                    values=[int(i) for i in traits_string],
                    )
        distribution_string = parts[2]
        distribution_vector = DistributionVector(
                num_areas=len(distribution_string),
                values=[int(i) for i in distribution_string],)
        return traits_vector, distribution_vector

    @staticmethod
    def decode_tree_lineages_from_labels(
            tree,
            leaf_nodes_only=False,
            encoded_source="node",
            ):
        if encoded_source == "node":
            _decode = lambda x: ArchipelagoModel.decode_label(x.label)
        elif encoded_source == "taxon":
            _decode = lambda x: ArchipelagoModel.decode_label(x.taxon.label)
        else:
            raise ValueError("'encoded_source' must be 'node' or 'taxon'")
        for nd in tree:
            if (not leaf_nodes_only or not nd._child_nodes) and (encoded_source == "node" or nd.taxon is not None):
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
            run_logger=None):

        # initialize
        if model_definition is None:
            model_definition = {}
        else:
            model_definition = dict(model_definition)

        # Geography
        if "areas" not in model_definition:
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
            self.lineage_birth_rate_function = RateFunction.from_definition(diversification_d.pop("lineage_birth_rate"), self.trait_types)
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
        ## extinction
        self.is_lineage_death_global = strtobool(str(diversification_d.pop("is_lineage_death_global", 0)))
        if "lineage_death_rate" in diversification_d:
            self.lineage_death_rate_function = RateFunction.from_definition(diversification_d.pop("lineage_death_rate"), self.trait_types)
        else:
            self.lineage_death_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.0",
                    description="fixed: 0.0",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Setting lineage death (= {is_global}) probability function: {desc}".format(
                is_global="global extinction" if self.is_lineage_death_global else "local extirpation",
                desc=self.lineage_death_rate_function.description,
                ))
        if diversification_d:
            raise TypeError("Unsupported diversification model keywords: {}".format(diversification_d))

        # Dispersal submodel
        dispersal_d = dict(model_definition.pop("dispersal", {}))
        if "lineage_dispersal_rate" in dispersal_d:
            self.lineage_dispersal_rate_function = RateFunction.from_definition(dispersal_d.pop("lineage_dispersal_rate"), self.trait_types)
        else:
            self.lineage_dispersal_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.01",
                    description="fixed: 0.01",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(DISPERSAL) Setting lineage dispersal rate function: {desc}".format(
                desc=self.lineage_dispersal_rate_function.description,))
        if dispersal_d:
            raise TypeError("Unsupported dispersal model keywords: {}".format(dispersal_d))

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
        model_definition["dispersal"] = self.dispersal_as_definition()
        model_definition["termination_conditions"] = self.termination_conditions_as_definition()
        json.dump(model_definition, out, indent=4, separators=(',', ': '))

    def diversification_as_definition(self):
        d = collections.OrderedDict()
        d["lineage_birth_rate"] = self.lineage_birth_rate_function.as_definition()
        d["lineage_death_rate"] = self.lineage_death_rate_function.as_definition()
        return d

    def dispersal_as_definition(self):
        d = collections.OrderedDict()
        d["lineage_dispersal_rate"] = self.lineage_dispersal_rate_function.as_definition()
        return d

    def termination_conditions_as_definition(self):
        d = collections.OrderedDict()
        d["target_focal_area_lineages"] = self.target_focal_area_lineages
        d["gsa_termination_focal_area_lineages"] = self.gsa_termination_focal_area_lineages
        d["max_time"] = self.max_time
        return d

class RateFunction(object):

    @classmethod
    def from_definition(cls, rate_function_d, trait_types):
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
        if d["definition_type"] == "function":
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

    def __getitem__(self, trait_index):
        return self._states[trait_index]

    def __setitem__(self, trait_index, v):
        self._states[trait_index] = v

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
                run_logger.info("(ECOLOGY) configuring trait {idx}: '{label}': {nstates} states, transition rate of {trate} with {normalized}transition weights of {tweights}".format(
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
            run_logger.info("(ECOLOGY) Total of {} traits defined{}{}".format(
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
            dispersal_weights=None,
            ):
        self.index = index
        self.label = label
        self.is_supplemental = is_supplemental
        self.relative_diversity = relative_diversity
        self.dispersal_weights = dispersal_weights

    def as_definition(self):
        d = collections.OrderedDict()
        # d["index"] = self.index
        d["label"] = self.label
        d["is_supplemental"] = self.is_supplemental
        d["relative_diversity"] = self.relative_diversity
        d["dispersal_weights"] = self.dispersal_weights
        return d

class Geography(object):

    def __init__(self):
        self.normalize_dispersal_weights = False

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
            area._dispersal_weights_d = list(area_d.pop("dispersal_weights", [])) # delay processing until all areas have been defined
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
        self.dispersal_weights = []
        for a1_idx, area1 in enumerate(self.areas):
            if len(area1._dispersal_weights_d) == 0:
                area1._dispersal_weights_d = [1.0] * len(self.areas)
                area1._dispersal_weights_d[a1_idx] = 0.0
            if len(area1._dispersal_weights_d) != len(self.areas):
                raise ValueError("Expecting exactly {} elements in dispersal weight vector for area '{}', but instead found {}: {}".format(
                    len(self.areas), area1.label, len(area1._dispersal_weights_d), area1._dispersal_weights_d))
            self.dispersal_weights.append([])
            for a2_idx, area2 in enumerate(self.areas):
                if a1_idx == a2_idx:
                    d = float(area1._dispersal_weights_d[a2_idx])
                    if d != 0:
                        raise ValueError("Self-dispersal weight from area {label} to {label} must be 0.0, but instead found: {dw}".format(
                            label=area1.label, dw=d))
                    self.dispersal_weights[a1_idx].append(0.0)
                else:
                    d = float(area1._dispersal_weights_d[a2_idx])
                    self.dispersal_weights[a1_idx].append(d)
            # if area1._dispersal_weights_d:
            #     raise ValueError("Undefined dispersal targets in '{}': '{}'".format(area1.label, area1._dispersal_weights_d))
            area1.dispersal_weights = self.dispersal_weights[a1_idx]
            del area1._dispersal_weights_d
        if self.normalize_dispersal_weights:
            for a1_idx, area1 in enumerate(self.areas):
                normalization_factor = sum(self.dispersal_weights[a1_idx])
                if normalization_factor:
                    for a2_idx, area2 in enumerate(self.areas):
                        self.dispersal_weights[a1_idx][a2_idx] /= normalization_factor
                    area1.dispersal_weights = self.dispersal_weights[a1_idx]
        if run_logger is not None:
            if self.normalize_dispersal_weights:
                weight_type = "Normalized dispersal"
            else:
                weight_type = "Dispersal"
            for a1, area1 in enumerate(self.areas):
                run_logger.info("(GEOGRAPHY) {} weights from area '{}': {}".format(weight_type, area1.label, self.dispersal_weights[a1]))

        # instead of recalculating every time
        self.area_nstates = [2 for i in self.areas]

    def as_definition(self):
        areas = [a.as_definition() for a in self.areas]
        return areas

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
                seed_node.distribution_vector[0] = 1
                kwargs["seed_node"] = seed_node
            dendropy.Tree.__init__(self, *args, **kwargs)
            self.current_lineages = set([self.seed_node])
        else:
            dendropy.Tree.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        memo[id(self.model)] = self.model
        memo[id(self.rng)] = self.rng
        memo[id(self.run_logger)] = self.run_logger
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
        num_areas = len(self.model.geography.area_indexes)
        if num_presences <= 1:
            speciation_mode = 0
        elif num_presences == num_areas:
            speciation_mode = self.rng.randint(0, 2)
        else:
            speciation_mode = self.rng.randint(0, 3)
        if speciation_mode == 0:
            dist1 = lineage.distribution_vector.clone()
            dist2 = lineage.distribution_vector.clone()
        elif speciation_mode == 1:
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.model.geography.new_distribution_vector()
            # TODO: area diversity base speciation
            dist2[ self.rng.choice(presences) ] = 1
        elif speciation_mode == 2:
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
        elif speciation_mode == 3:
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.model.geography.new_distribution_vector()
            absences = [idx for idx in self.model.geography.area_indexes if idx not in presences]
            dist2[ self.rng.choice(absences) ] = 1
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

    def extirpate_lineage(self, lineage):
        if self.model.is_lineage_death_global:
            self._make_lineage_extinct_on_phylogeny(lineage)
        else:
            presences = lineage.distribution_vector.presences()
            assert len(presences) > 0
            if len(presences) == 1:
                self._make_lineage_extinct_on_phylogeny(lineage)
            else:
                lineage.distribution_vector[ self.rng.choice(presences) ] = 0

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


