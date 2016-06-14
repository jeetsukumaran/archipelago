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

def set_node_times_ages_extancy(tree,
        is_set_age=True,
        is_annotate=True,
        time_to_add_to_extant_tips=0,
        ):
    max_node_time = None
    for nd in tree:
        if nd.parent_node:
            nd.time = nd.parent_node.time + nd.edge.length
        else:
            nd.time = nd.edge.length
        if max_node_time is None or max_node_time < nd.time:
            max_node_time = nd.time
    for nd in tree:
        age = max(max_node_time - nd.time, 0)
        if hasattr(nd, "is_extant"):
            if nd.is_extant:
                a = abs(age)
                assert a < 1e-8, a
                nd.time = max_node_time + time_to_add_to_extant_tips
                nd.edge.length += time_to_add_to_extant_tips
                age = 0.0
            else:
                age += time_to_add_to_extant_tips
        else:
            if abs(age) < 1e-8:
                nd.is_extant = True
                nd.time = max_node_time + time_to_add_to_extant_tips
                nd.edge.length += time_to_add_to_extant_tips
                age = 0.0
            else:
                nd.is_extant = False
                age += time_to_add_to_extant_tips
        if is_set_age:
            nd.age = age
        if is_annotate:
            nd.annotations["is_extant"] = nd.is_extant
            nd.annotations["time"] = nd.time
            if is_set_age:
                nd.annotations["age"] = nd.age
    tree.max_node_time = max_node_time


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

    def trait_index(self, label):
        idx = self.trait_label_index_map[label]
        return idx

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
                    idx=trait.index+1,
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

    def __call__(self, **kwargs):
        return self._compute_rate(**kwargs)

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
            self._compute_rate = lambda **kwargs: self.definition_content
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
            self._compute_rate = lambda **kwargs: rates[kwargs["lineage"].traits_vector[trait.index]]
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

class Area(object):

    def __init__(self,
            index=None,
            focal_area_index=None,
            label=None,
            is_supplemental=False,
            relative_diversity=None,
            area_connection_weights=None,
            ):
        self.index = index
        self.focal_area_index = focal_area_index # index only considering focal areas
        self.label = label
        self.is_supplemental = is_supplemental
        self.relative_diversity = relative_diversity
        # this is here mainly for description purposes in
        # `Area.as_definition()`; actual usage is through
        # `Geography.area_connection_weights`
        self.area_connection_weights = area_connection_weights
        self.lineages = set()

    def __str__(self):
        return "Area_{}".format(self.index)

    def __repr__(self):
        return "<archipelago.model.Area object at {} with index {}>".format(id(self), self.index)

    def add_lineage(self, lineage):
        assert lineage not in self.lineages
        self.lineages.add(lineage)

    def remove_lineage(self, lineage):
        assert lineage in self.lineages
        self.lineages.remove(lineage)

    def debug_check(self):
        for lineage in self.lineages:
           assert lineage is extant
           assert self in lineage.areas

    def as_definition(self):
        d = collections.OrderedDict()
        # d["index"] = self.index
        d["label"] = self.label
        d["is_supplemental"] = self.is_supplemental
        d["relative_diversity"] = self.relative_diversity
        d["area_connection_weights"] = self.area_connection_weights
        return d

class Lineage(dendropy.Node):

    _TRAITS_SEPARATOR = "."
    _LABEL_COMPONENTS_SEPARATOR = "^"
    _NULL_TRAITS = "NA"

    class NullDistributionException(Exception):

        def __init__(self, lineage):
            self.lineage = lineage

    def __init__(self, index, model, geography, log_event):
        dendropy.Node.__init__(self)
        self.index = index
        self.model = model
        self.geography = geography
        self.log_event = log_event
        self.areas = set()
        self.distribution_vector = None # not used in simulation; assigned and used in profile/summary statistic calculation
        self.traits_vector = self.model.trait_types.new_traits_vector()
        self.is_extant = True
        self.edge.length = 0
        if self.log_event is not None:
            self.starting_distribution_bitstring = None # used for history logging
            self.ending_distribution_bitstring = None # used for history logging
            self.starting_focal_area_distribution_bitstring = None # used for history logging
            self.ending_focal_area_distribution_bitstring = None # used for history logging

    def __str__(self):
        return "<{} with index {} (root: {}, children: {})>".format(self.__class__.__name__, self.index, self.parent_node is None, len(self._child_nodes) == 0)

    def register_current_distribution_as_starting_distribution(self):
        self.starting_distribution_bitstring = self.distribution_bitstring(exclude_supplemental_areas=False)
        self.starting_focal_area_distribution_bitstring = self.distribution_bitstring(exclude_supplemental_areas=True)

    def register_current_distribution_as_ending_distribution(self):
        self.ending_distribution_bitstring = self.distribution_bitstring(exclude_supplemental_areas=False)
        self.ending_focal_area_distribution_bitstring = self.distribution_bitstring(exclude_supplemental_areas=True)

    def copy_areas(self, other):
        for area in self.geography.areas:
            if area in other.areas and area not in self.areas:
                self.add_area(area)
            elif area not in other.areas and area in self.areas:
                self.remove_area(area)

    def add_area(self, area, is_log_event=False):
        assert area not in self.areas
        self.areas.add(area)
        area.lineages.add(self)
        if is_log_event and self.log_event is not None:
            self.log_event(
                    lineage=self,
                    event_type="geography_anagenesis",
                    event_subtype="area_gain",
                    state_idx=area.index,
                    child0_lineage=None,
                    child1_lineage=None)

    def add_areas(self, areas, is_log_event=False):
        for area in areas:
            self.add_area(area, is_log_event=is_log_event)

    def remove_area(self, area, is_log_event=False):
        assert area in self.areas
        self.areas.remove(area)
        area.lineages.remove(self)
        if is_log_event and self.log_event is not None:
            self.log_event(
                    lineage=self,
                    event_type="geography_anagenesis",
                    event_subtype="area_loss",
                    state_idx=area.index,
                    child0_lineage=None,
                    child1_lineage=None)
        if len(self.areas) == 0:
            raise Lineage.NullDistributionException(self)

    def deactivate(self):
        self.is_extant = False
        self.clear_areas()

    def clear_areas(self):
        for area in self.areas:
            area.lineages.discard(self)
        self.areas = set()

    def copy_traits(self, other_lineage):
        self.traits_vector = other_lineage.traits_vector.clone()

    def distribution_bitstring(self, exclude_supplemental_areas=False):
        d = []
        if exclude_supplemental_areas:
            areas = self.geography.focal_areas
        else:
            areas = self.geography.areas
        for a in areas:
            if a in self.areas:
                d.append("1")
            else:
                d.append("0")
        return "".join(d)

    def compose_encoded_label(self, exclude_supplemental_areas=None):
        if self.traits_vector:
            traits_v = Lineage._TRAITS_SEPARATOR.join(str(i) for i in self.traits_vector)
        else:
            traits_v = Lineage._NULL_TRAITS
        areas_v = self.distribution_bitstring(exclude_supplemental_areas=exclude_supplemental_areas)
        encoding = "s{self_index}{sep}{traits_v}{sep}{areas_v}".format(
                self_index=self.index,
                traits_v=traits_v,
                areas_v=areas_v,
                sep=Lineage._LABEL_COMPONENTS_SEPARATOR)
        return encoding

    def encode_lineage(self,
            set_label=False,
            add_annotation=False,
            exclude_supplemental_areas=False,
            ):
        label = self.compose_encoded_label(exclude_supplemental_areas=exclude_supplemental_areas)
        if set_label:
            self.label = label
        if add_annotation:
            self.annotations.drop()
            self.annotations.add_new("traits", traits_v)
            self.annotations.add_new("distribution", areas_v)
            for trait_idx, trait in enumerate(self.model.trait_types):
                self.annotations.add_new(trait.label, self.traits_vector[trait_idx])
            area_list = []
            for area_idx, area in enumerate(self.geography.areas):
                if exclude_supplemental_areas and area.is_supplemental:
                    continue
                if self.distribution_vector[area_idx] == 1:
                    area_list.append(area.label)
            self.annotations.add_new("areas", area_list)
        return label

    def debug_check(self):
        for area in self.geography.areas:
            if area in self.areas:
                assert self in area.lineages
            else:
                assert self not in area.lineages

    def trait_state(self, trait_label):
        return self.traits_vector[self.model.trait_types.trait_index(trait_label)]

class Phylogeny(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return Lineage(**kwargs)
    node_factory = classmethod(node_factory)

    def __init__(self, *args, **kwargs):
        self.model = kwargs.pop("model", None)
        self.geography = kwargs.pop("geography", None)
        self.model_id = self.model.model_id
        self.annotations.add_bound_attribute("model_id")
        self.rng = kwargs.pop("rng", None)
        self.log_event = kwargs.pop("log_event", None)
        self.debug_mode = kwargs.pop("debug_mode", None)
        self.run_logger = kwargs.pop("run_logger", None)
        self.lineage_indexer = utility.IndexGenerator(0)
        assert "seed_node" not in kwargs
        seed_node = self.node_factory(
                index=next(self.lineage_indexer),
                model=self.model,
                geography=self.geography,
                log_event=self.log_event,
                )
        kwargs["seed_node"] = seed_node
        dendropy.Tree.__init__(self, *args, **kwargs)
        self.is_rooted = True

    def bootstrap(self):
        for trait_idx in range(len(self.model.trait_types)):
            trait_states = [i for i in range(self.model.trait_types[trait_idx].nstates)]
            self.seed_node.traits_vector[trait_idx] = self.rng.choice(trait_states)
        initial_area = self.rng.choice(self.geography.areas)
        self.seed_node.add_area(initial_area, is_log_event=True)
        self.current_lineages = set([self.seed_node])
        if self.log_event:
            self.seed_node.register_current_distribution_as_starting_distribution()

    def iterate_current_lineages(self):
        for lineage in self.current_lineages:
            yield lineage

    def split_lineage(self, lineage, area=None):
        c1 = self.node_factory(index=next(self.lineage_indexer), model=self.model, geography=self.geography, log_event=self.log_event)
        c2 = self.node_factory(index=next(self.lineage_indexer), model=self.model, geography=self.geography, log_event=self.log_event)
        c1.copy_traits(lineage)
        c2.copy_traits(lineage)
        speciation_mode_idx = self._set_daughter_distributions(
                parent=lineage,
                child1=c1,
                child2=c2,
                speciation_area=area)
        if self.log_event:
            if speciation_mode_idx == 0:
                speciation_mode_desc = "single-area sympatric speciation"
            elif speciation_mode_idx == 1:
                speciation_mode_desc = "subset sympatry"
            elif speciation_mode_idx == 2:
                speciation_mode_desc = "single area vicariance"
            elif speciation_mode_idx == 3:
                speciation_mode_desc = "multi-area vicariance"
            elif speciation_mode_idx == 4:
                speciation_mode_desc = "founder-event jump dispersal"
            self.log_event(
                    lineage=lineage,
                    event_type="cladogenesis",
                    event_subtype=speciation_mode_desc,
                    state_idx=None,
                    child0_lineage=c1,
                    child1_lineage=c2)
            c1.register_current_distribution_as_starting_distribution()
            c2.register_current_distribution_as_starting_distribution()
            lineage.register_current_distribution_as_ending_distribution()
        lineage.deactivate()
        self.current_lineages.remove(lineage)
        lineage.add_child(c1)
        lineage.add_child(c2)
        self.current_lineages.add(c1)
        self.current_lineages.add(c2)

    def extinguish_lineage(self, lineage, extinction_type):
        if len(self.current_lineages) == 1:
            self.total_extinction_exception("No extant lineages remaining")
        lineage.deactivate()
        self.current_lineages.remove(lineage)
        if self.log_event:
            self.log_event(
                    lineage=lineage,
                    event_type="extinction",
                    event_subtype=extinction_type,
                    state_idx=None,
                    child0_lineage=None,
                    child1_lineage=None,
                    )

    def extinguish_lineage_OLD(self, lineage, event_type):
        # assert not lineage._child_nodes
        if len(self.current_lineages) == 1:
            self.total_extinction_exception("No extant lineages remaining")
        if lineage is not self.seed_node:
            s = lineage.parent_node._child_nodes
            assert len(s) > 1
            if s[0] is lineage:
                lineage_to_inherit = s[1]
            else:
                lineage_to_inherit = s[0]
        else:
            lineage_to_inherit = None
        # if self.log_event is not None:
        #     if lineage is not self.seed_node:
        #         s = lineage.parent_node._child_nodes
        #         assert len(s) > 1
        #         if s[0] is lineage:
        #             lineage_inheriting_events = s[1]
        #         else:
        #             lineage_inheriting_events = s[0]
        #     self.log_event(
        #             lineage=lineage,
        #             event_type="extinction",
        #             lineage_inheriting_events=lineage_inheriting_events,
        #             )
        lineage.deactivate()
        self.current_lineages.remove(lineage)
        self.prune_subtree(lineage, suppress_unifurcations=False, update_bipartitions=False)
        remapped_nodes = self.suppress_unifurcations(update_bipartitions=False)
        if lineage_to_inherit is not None:
            remapped_nodes.insert( 0, (lineage, lineage_to_inherit))
        self.log_event(
                lineage=lineage,
                event_type="extinction",
                remapped_nodes=remapped_nodes)

    def _set_daughter_distributions(self, parent, child1, child2, speciation_area):
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
        num_presences = len(parent.areas)
        num_areas = len(self.geography.areas)
        if num_presences <= 1:
            speciation_mode = 0
        else:
            if num_presences < num_areas:
                area_gain_event_parameters, area_gain_event_rates, area_gain_rates_marginalized_by_destination_area = self.geography.calculate_raw_area_gain_events(
                        lineage=parent,
                        lineage_area_gain_weight_fn=self.model.lineage_area_gain_weight_function,
                        simulation_elapsed_time=None)
                if area_gain_event_parameters and area_gain_event_rates:
                    jump_dispersal_target_area_rates = [0.0 for a in self.geography.areas]
                    normalization_factor = float(sum(area_gain_event_rates))
                    if normalization_factor:
                        for ag_event_parameters, ag_event_rate in zip(area_gain_event_parameters, area_gain_event_rates):
                            jump_dispersal_target_area_rates[ag_event_parameters["to_area"].index] += (ag_event_rate/normalization_factor)
                        fes_weight = self.model.global_area_gain_rate * self.model.cladogenesis_founder_event_speciation_weight * sum(jump_dispersal_target_area_rates)
                    else:
                        fes_weight = 0.0
                else:
                    fes_weight = 0.0
                # if self.model.is_area_specific_gain_rate:
                #     for jd_area in self.geography.areas:
                #         if jd_area in parent.areas:
                #             continue
                #         jump_dispersal_target_areas.append(jd_area)
                #         jump_dispersal_target_rates.append(self.model.lineage_area_gain_weight_function(lineage=parent, area=jd_area))
                # else:
                #     parent_area_gain_rate = self.model.lineage_area_gain_weight_function(lineage=parent, area=None)
                #     for jd_area in self.geography.areas:
                #         if jd_area in parent.areas:
                #             continue
                #         jump_dispersal_target_areas.append(jd_area)
                #     per_area_gain_rate = parent_area_gain_rate / len(jump_dispersal_target_areas)
                #     jump_dispersal_target_rates = [per_area_gain_rate] * len(jump_dispersal_target_areas)
                # fes_weight = self.model.cladogenesis_founder_event_speciation_weight * sum(jump_dispersal_target_rates)
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
            child1.copy_areas(parent)
            child2.copy_areas(parent)
        elif speciation_mode == 1:
            # sympatric subset: multi-area sympatric speciation
            #     -   d1: inherits complete range
            #     -   d2: inherits single area in ancestral range
            child1.copy_areas(parent)
            if speciation_area is None:
                child2.add_area( self.rng.choice(tuple(parent.areas)) )
            else:
                assert speciation_area in parent.areas
                child2.add_area(speciation_area)
        elif speciation_mode == 2:
            # (single-area) allopatric vicariance
            #     -   d1: single area
            #     -   d2: all other areas
            if speciation_area is None:
                selected_idx = self.rng.randrange(len(parent.areas))
                for aidx, area in enumerate(parent.areas):
                    if aidx == selected_idx:
                        child2.add_area(area)
                    else:
                        child1.add_area(area)
            else:
                assert speciation_area in parent.areas
                for area in parent.areas:
                    if area is speciation_area:
                        child1.add_area(area)
                    else:
                        child2.add_area(area)
        elif speciation_mode == 3:
            if num_presences == 2:
                dist = list(parent.areas)
                self.rng.shuffle(dist)
                child1.add_area(dist[0])
                child2.add_area(dist[1])
            else:
                n1 = self.rng.randint(1, num_presences-1)
                n2 = num_presences - n1
                if n2 == n1:
                    n1 += 1
                    n2 -= 1
                assert n1 > 0
                assert n2 > 0
                sample1 = set(self.rng.sample(parent.areas, n1))
                # assert len(sample1) < len(self.geography.areas)
                for s_area in parent.areas:
                    if s_area in sample1:
                        child1.add_area(s_area)
                    else:
                        child2.add_area(s_area)
                assert len(child1.areas) + len(child2.areas) == len(parent.areas)
                # print("---")
                # print([a.index for a in parent.areas])
                # print([a.index for a in child1.areas])
                # print([a.index for a in child2.areas])
        elif speciation_mode == 4:
            child1.copy_areas(parent)
            jump_target_idx = weighted_index_choice(
                    weights=jump_dispersal_target_area_rates,
                    sum_of_weights=fes_weight,
                    rng=self.rng)
            jump_target_area = self.geography.areas[jump_target_idx]
            child2.add_area(jump_target_area)
        else:
            raise ValueError(speciation_mode)
        if self.debug_mode:
            self.run_logger.debug("Splitting {} with distribution {} under speciation mode {} to: {} (distribution: {}) and {} (distribution: {})".format(
                parent,
                parent.distribution_bitstring(),
                speciation_mode,
                child1,
                child1.distribution_bitstring(),
                child2,
                child2.distribution_bitstring(),
                ))
            assert len(child1.areas) > 0
            assert len(child2.areas) > 0
        return speciation_mode

    def total_extinction_exception(self, msg):
        # self.run_logger.info("Total extinction: {}".format(msg))
        raise error.TotalExtinctionException(msg)

    def evolve_trait(self, lineage, trait_idx, state_idx):
        lineage.traits_vector[trait_idx] = state_idx
        if self.log_event is not None:
            self.log_event(
                    lineage=lineage,
                    event_type="trait_evolution",
                    event_subtype="{}".format(trait_idx),
                    state_idx=trait_idx,
                    child0_lineage=None,
                    child1_lineage=None)

    def focal_area_lineages(self):
        focal_area_lineages = set()
        for area in self.geography.focal_areas:
            focal_area_lineages.update(area.lineages)
        return focal_area_lineages

    def num_focal_area_lineages(self):
        return len(self.focal_area_lineages())

    def _compose_tree_string(self,
            tree,
            node_label_compose_fn,
            is_suppress_internal_node_labels,
            ):
        s = StringIO()
        tree.write_to_stream(
                s,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_fn=node_label_compose_fn,
                suppress_internal_node_labels=is_suppress_internal_node_labels,
                )
        return s.getvalue()

    def generate_tree_strings_for_serialization(self,
            is_encode_nodes,
            is_annotate_nodes,
            is_suppress_internal_node_labels,
            ):

        # storage
        results = {}

        # prep
        all_areas_node_labels = {}
        all_areas_node_annotations = {}
        focal_areas_node_labels = {}
        focal_areas_node_annotations = {}
        if is_encode_nodes:
            all_areas_labelf = lambda x: x.encode_lineage(
                    set_label=False,
                    add_annotation=is_annotate_nodes,
                    exclude_supplemental_areas=False)
            focal_areas_labelf = lambda x: x.encode_lineage(
                    set_label=False,
                    add_annotation=is_annotate_nodes,
                    exclude_supplemental_areas=True)
        else:
            all_areas_labelf = ArchipelagoSimulator.simple_node_label_function
        focal_area_lineages = set()
        extinct_lineages = set()
        for nd in self:
            focal_areas_node_labels[nd] = focal_areas_labelf(nd)
            focal_areas_node_annotations[nd] = dendropy.AnnotationSet(nd.annotations)
            all_areas_node_labels[nd] = all_areas_labelf(nd)
            all_areas_node_annotations[nd] = dendropy.AnnotationSet(nd.annotations)
            if nd.is_leaf() and not nd.is_extant:
                extinct_lineages.add(nd)
            for area in self.geography.focal_areas:
                if area in nd.areas:
                    focal_area_lineages.add(nd)
                    break

        # all areas, complete
        results["all-areas.complete"] = self._compose_tree_string(
                tree=self,
                node_label_compose_fn=lambda nd: all_areas_node_labels[nd],
                is_suppress_internal_node_labels=is_suppress_internal_node_labels,
                )

        # all areas, extant only
        all_areas_extant_only_tree = self.extract_tree(
                node_filter_fn=lambda x: x not in extinct_lineages,
                tree_factory=dendropy.Tree,
                node_factory=dendropy.Node,
                )
        for nd in all_areas_extant_only_tree:
            nd.annotations = all_areas_node_annotations[nd.extraction_source]
        results["all-areas.extant"] = self._compose_tree_string(
                tree=all_areas_extant_only_tree,
                node_label_compose_fn=lambda nd: all_areas_node_labels[nd.extraction_source],
                is_suppress_internal_node_labels=is_suppress_internal_node_labels,
                )

        # focal areas
        focal_areas_extant_only_tree = self.extract_tree(
                node_filter_fn=lambda x: x not in extinct_lineages and x in focal_area_lineages,
                tree_factory=dendropy.Tree,
                node_factory=dendropy.Node,
                )
        for nd in focal_areas_extant_only_tree:
            nd.annotations = nd.extraction_source.annotations
        results["focal-areas"] = self._compose_tree_string(
                tree=focal_areas_extant_only_tree,
                node_label_compose_fn=lambda nd: focal_areas_node_labels[nd.extraction_source],
                is_suppress_internal_node_labels=is_suppress_internal_node_labels,
                )

        return results

    # def extract_focal_areas_tree(self):
        # # tcopy = Phylogeny(self)
        # # tcopy = copy.deepcopy(self)
        # tree_factory = lambda **kwargs: Phylogeny(
        #         model=self.model,
        #         geography=self.geography,
        #         rng=None,
        #         **kwargs
        #         )
        # node_factory = lambda: Lineage(index=None, model=self.model, geography=self.geography, log_event=None)
        # tcopy = self.extract_tree(
        #         tree_factory=tree_factory,
        #         node_factory=node_factory)
        # focal_area_lineages = set()
        # for nd in tcopy:
        #     nd.index = nd.extraction_source.index
        #     nd.areas = set(nd.extraction_source.areas)
        #     if nd.extraction_source.distribution_vector:
        #         nd.distribution_vector = nd.extraction_source.distribution_vector.clone()
        #     if nd.extraction_source.traits_vector:
        #         nd.traits_vector = nd.extraction_source.traits_vector.clone()
        #     if self.log_event:
        #         nd.starting_distribution_bitstring = nd.extraction_source.starting_distribution_bitstring
        #         nd.ending_distribution_bitstring = nd.extraction_source.ending_distribution_bitstring
        #         nd.starting_focal_area_distribution_bitstring = nd.extraction_source.starting_focal_area_distribution_bitstring
        #         nd.ending_focal_area_distribution_bitstring = nd.extraction_source.ending_focal_area_distribution_bitstring
        #     for area in self.geography.focal_areas:
        #         if area in nd.areas:
        #             focal_area_lineages.add(nd)
        #             break
        # if len(focal_area_lineages) < 2:
        #     raise error.InsufficientFocalAreaLineagesSimulationException("insufficient lineages in focal area at termination".format(len(focal_area_lineages)))
        # try:
        #     tcopy.filter_leaf_nodes(filter_fn=lambda x: x in focal_area_lineages)
        # except dendropy.SeedNodeDeletionException:
        #     raise error.InsufficientFocalAreaLineagesSimulationException("no extant lineages in focal area at termination".format(len(focal_area_lineages)))
        # return tcopy

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        memo[id(self.model)] = self.model
        memo[id(self.rng)] = None #self.rng
        memo[id(self.run_logger)] = self.run_logger
        memo[id(self.taxon_namespace)] = self.taxon_namespace
        return dendropy.Tree.__deepcopy__(self, memo)

class Geography(object):

    def __init__(self,
            areas_definition,
            run_logger=None):
        self.normalize_area_connection_weights = False
        self.parse_definition(areas_definition=areas_definition, run_logger=run_logger)

    def __iter__(self):
        return iter(self.areas)

    def __getitem__(self, idx):
        return self.areas[idx]

    def parse_definition(self, areas_definition, run_logger=None):
        self.areas = []
        self.area_label_index_map = collections.OrderedDict()
        self.area_indexes = []
        self.focal_areas = []
        self.supplemental_areas = []
        for area_idx, area_d in enumerate(areas_definition):
            if "label" not in area_d:
                raise ValueError("Area requires 'label' to be defined")
            area_label = str(area_d.pop("label", None))
            if area_label in self.area_label_index_map:
                raise ValueError("Area with label '{}' already defined".format(area_label))
            area = Area(
                index=area_idx,
                focal_area_index=area_d.pop("focal_area_index", -1),
                label=area_label,
                relative_diversity=area_d.pop("relative_diversity", 1.0),
                is_supplemental=area_d.pop("is_supplemental", False)
            )
            if area.is_supplemental:
                assert area.focal_area_index == -1
            else:
                assert area.focal_area_index == len(self.focal_areas)
            area._area_connection_weights_d = list(area_d.pop("area_connection_weights", [])) # delay processing until all areas have been defined
            self.areas.append(area)
            self.area_label_index_map[area.label] = area.index
            self.area_indexes.append(area.index)
            if area.is_supplemental:
                self.supplemental_areas.append(area)
            else:
                self.focal_areas.append(area)
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
        if len(self.focal_areas) < 1:
            raise ValueError("No focal areas defined")
        if run_logger is not None:
            run_logger.info("(GEOGRAPHY) Total of {} areas defined: {}".format(
                len(self.areas),
                ", ".join("'{}'".format(a.label) for a in self.areas),
                ))
        self.area_connection_weights = []
        for a1_idx, area1 in enumerate(self.areas):
            if len(area1._area_connection_weights_d) == 0:
                area1._area_connection_weights_d = [1.0] * len(self.areas)
                area1._area_connection_weights_d[a1_idx] = 0.0
            if len(area1._area_connection_weights_d) != len(self.areas):
                raise ValueError("Expecting exactly {} elements in area gain weight vector for area '{}', but instead found {}: {}".format(
                    len(self.areas), area1.label, len(area1._area_connection_weights_d), area1._area_connection_weights_d))
            self.area_connection_weights.append([])
            for a2_idx, area2 in enumerate(self.areas):
                if a1_idx == a2_idx:
                    d = float(area1._area_connection_weights_d[a2_idx])
                    if d != 0:
                        raise ValueError("Area gain weight from area {label} to {label} must be 0.0, but instead found: {dw}".format(
                            label=area1.label, dw=d))
                    self.area_connection_weights[a1_idx].append(0.0)
                else:
                    d = float(area1._area_connection_weights_d[a2_idx])
                    self.area_connection_weights[a1_idx].append(d)
            # if area1._area_connection_weights_d:
            #     raise ValueError("Undefined dispersal targets in '{}': '{}'".format(area1.label, area1._area_connection_weights_d))
            area1.area_connection_weights = self.area_connection_weights[a1_idx]
            del area1._area_connection_weights_d
        if self.normalize_area_connection_weights:
            for a1_idx, area1 in enumerate(self.areas):
                normalization_factor = sum(self.area_connection_weights[a1_idx])
                if normalization_factor:
                    for a2_idx, area2 in enumerate(self.areas):
                        self.area_connection_weights[a1_idx][a2_idx] /= normalization_factor
                    area1.area_connection_weights = self.area_connection_weights[a1_idx]
        if run_logger is not None:
            if self.normalize_area_connection_weights:
                weight_type = "Normalized area gain"
            else:
                weight_type = "Area gain"
            for a1, area1 in enumerate(self.areas):
                run_logger.info("(GEOGRAPHY) {} weights from area '{}': {}".format(weight_type, area1.label, self.area_connection_weights[a1]))

        # instead of recalculating every time
        self.area_nstates = [2 for i in self.areas]

    def as_definition(self):
        areas = [a.as_definition() for a in self.areas]
        return areas

    def calculate_raw_area_gain_events(self,
            lineage,
            lineage_area_gain_weight_fn,
            simulation_elapsed_time=None):
        # Submodel Design Objectives:
        # 1.    We want to be able to specify that the rate of gaining a
        #       particular area are functions of:
        #       -   The area (e.g., the number of area lineages in the area; the
        #           number of symbiont lineages in the area; or the particular
        #           area or symbiont lineages present/absent from an area).
        #       -   The source and destination areas (e.g., the phylogenetic
        #           distances between the source and destination areas; the
        #           number of resident lineages in the destination area; etc.)
        #       -   The particular characters/trait of the dispersing lineage.
        # 2.    We want to be able to specify a mean per-lineage rate of
        #       transmission. Allows for estimating this from empirical data
        #       using some reasonable if simplified and low-fidelity model:
        #       e.g., as a per-lineage trait evolution rate, where the area set
        #       is a multistate character.
        # Objective (1) means that the rates must be calculated on a
        # per-source area per-destination area per area basis.
        # Objective (2) means that the transmission (area gain) rate
        # weights across all area gain events needs to sum to 1 (with the
        # actual rate obtained by multiplying with the system-wide mean
        # (per-lineage) area gain rate.)
        src_areas = []
        dest_areas = []
        area_gain_rates_marginalized_by_destination_area = []
        for area in self.areas:
            area_gain_rates_marginalized_by_destination_area.append(0.0)
            if area in lineage.areas:
                src_areas.append(area)
            else:
                dest_areas.append(area)
        area_gain_event_parameters = []
        area_gain_event_rates = []
        if src_areas and dest_areas:
            for src_area in src_areas:
                for dest_area in dest_areas:
                    lineage_area_gain_weight = lineage_area_gain_weight_fn(
                            lineage=lineage,
                            from_area=src_area,
                            to_area=dest_area,
                            simulation_elapsed_time=simulation_elapsed_time)
                    rate = lineage_area_gain_weight * self.area_connection_weights[src_area.index][dest_area.index]
                    if rate:
                        area_gain_event_parameters.append({"from_area": src_area, "to_area": dest_area})
                        area_gain_event_rates.append(rate)
                        area_gain_rates_marginalized_by_destination_area[dest_area.index] += rate
        return area_gain_event_parameters, area_gain_event_rates, area_gain_rates_marginalized_by_destination_area

class ArchipelagoModel(object):

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
    def decode_traits_and_distribution_from_label(label_to_decode):
        parts = label_to_decode.split(Lineage._LABEL_COMPONENTS_SEPARATOR)
        traits_string = parts[1]
        if not traits_string or traits_string == Lineage._NULL_TRAITS:
            traits_vector = StatesVector(nchar=0)
        else:
            traits_string_parts = traits_string.split(Lineage._TRAITS_SEPARATOR)
            traits_vector = StatesVector(
                    nchar=len(traits_string_parts),
                    # The trait states need to be an integer if
                    # archipelago-summarize.py coerces the user input to
                    # integers:
                    # values=[int(i) for i in traits_string_parts],
                    # On the other hand, the reason we do NOT want it parsed to
                    # an integer value is to allow for null traits to be
                    # expressed via 'NA', 'null', etc.
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
        ## Typically used when reading output tree for summary statistics or
        ## profile calculation
        if lineage_data_source == "node":
            _decode = lambda x: ArchipelagoModel.decode_traits_and_distribution_from_label(x.label)
        elif lineage_data_source == "taxon":
            _decode = lambda x: ArchipelagoModel.decode_traits_and_distribution_from_label(x.taxon.label)
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

        # mode identification
        if "model_id" not in model_definition:
            model_definition["model_id"] = "Model1"
            if run_logger is not None:
                run_logger.warning("Model identifier not specified: defaulting to '{}'".format(model_definition["model_id"]))
        self.model_id = model_definition.pop("model_id", "Model1")
        if run_logger is not None:
            run_logger.info("Setting up model with identifier: '{}'".format(self.model_id))

        # geography
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
        self.geographical_definition = copy.deepcopy(model_definition.pop("areas"))
        focal_area_index = 0
        for area in self.geographical_definition:
            if area["is_supplemental"]:
                area["focal_area_index"] = -1
            else:
                area["focal_area_index"] = focal_area_index
                focal_area_index += 1
        self._example_geography = self.new_geography(run_logger=run_logger)

        # Ecology
        self.trait_types = TraitTypes()
        self.trait_types.parse_definition(
                copy.deepcopy(model_definition.pop("traits", [])),
                run_logger=run_logger)

        # Diversification
        diversification_d = dict(model_definition.pop("diversification", {}))
        ## speciation
        self.mean_birth_rate = diversification_d.pop("mean_birth_rate", 0.10)
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Mean diversification birth rate: {}".format(self.mean_birth_rate))
        if "lineage_birth_weight" in diversification_d:
            self.lineage_birth_weight_function = RateFunction.from_definition_dict(diversification_d.pop("lineage_birth_weight"), self.trait_types)
        else:
            self.lineage_birth_weight_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda **kwargs: 1.00",
                    description="fixed: 1.00",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Setting lineage-specific birth weight function: {desc}".format(
                desc=self.lineage_birth_weight_function.description,))
        ## (global) extinction
        self.mean_death_rate = diversification_d.pop("mean_death_rate", 0.00)
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Mean diversification death rate: {}".format(self.mean_death_rate))
        if "lineage_death_weight" in diversification_d:
            self.lineage_death_weight_function = RateFunction.from_definition_dict(diversification_d.pop("lineage_death_weight"), self.trait_types)
        else:
            self.lineage_death_weight_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda **kwargs: 1.0",
                    description="fixed: 1.0",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Setting lineage-specific death weight function: {desc}".format(
                desc=self.lineage_death_weight_function.description,))
        if diversification_d:
            raise TypeError("Unsupported diversification model keywords: {}".format(diversification_d))

        # Dispersal submodel
        anagenetic_range_evolution_d = dict(model_definition.pop("anagenetic_range_evolution", {}))
        # if "global_area_gain_rate" not in anagenetic_range_evolution_d and "mean_area_gain_rate" not in anagenetic_range_evolution_d:
        #     if interpolate_missing_model_values:
        #         anagenetic_range_evolution_d["global_area_gain_rate"] = 1.0
        #     else:
        #         raise TypeError("Exactly one of 'global_area_gain_rate' or 'mean_area_gain_rate' must be specified")
        # if "global_area_gain_rate" in anagenetic_range_evolution_d and "mean_area_gain_rate" in anagenetic_range_evolution_d:
        #     raise TypeError("No more than one of 'global_area_gain_rate' or 'mean_area_gain_rate' can be specified")
        # elif "global_area_gain_rate" in anagenetic_range_evolution_d:
        #     self.global_area_gain_rate = float(anagenetic_range_evolution_d.pop("global_area_gain_rate"))
        #     self.mean_area_gain_rate = None
        #     if run_logger is not None:
        #         run_logger.info("(ANAGENETIC RANGE EVOLUTION) Global area gain rate is: {}".format(self.global_area_gain_rate))
        #     self.geography.set_global_area_gain_rate(self.global_area_gain_rate)
        # else:
        #     self.mean_area_gain_rate = float(anagenetic_range_evolution_d.pop("mean_area_gain_rate"))
        #     self.global_area_gain_rate = None
        #     run_logger.info("(ANAGENETIC RANGE EVOLUTION) Mean area gain rate is: {}".format(self.mean_area_gain_rate))
        #     self.geography.set_mean_area_gain_rate(self.mean_area_gain_rate)
        # if run_logger is not None:
        #     for a1, area1 in enumerate(self.geography.areas):
        #         run_logger.info("(ANAGENETIC RANGE EVOLUTION) Effective rate of area gain from area '{}': {}".format(area1.label, self.geography.effective_area_gain_rates[a1]))

        self.global_area_gain_rate = anagenetic_range_evolution_d.pop("global_area_gain_rate", 0.20)
        if run_logger is not None:
            run_logger.info("(ANAGENETIC RANGE EVOLUTION) Global area gain rate: {}".format(self.global_area_gain_rate))
        if "lineage_area_gain_weight" in anagenetic_range_evolution_d:
            self.lineage_area_gain_weight_function = RateFunction.from_definition_dict(anagenetic_range_evolution_d.pop("lineage_area_gain_weight"), self.trait_types)
        else:
            self.lineage_area_gain_weight_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda **kwargs: 1.00",
                    description="fixed: 1.00",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC RANGE EVOLUTION) Setting lineage-specific area gain weight function: {desc}".format(
                desc=self.lineage_area_gain_weight_function.description,))

        ## extinction
        # self.treat_area_loss_rate_as_lineage_death_weight = strtobool(str(anagenetic_range_evolution_d.pop("treat_area_loss_rate_as_lineage_death_weight", 0)))
        # self.is_area_specific_loss_rate = anagenetic_range_evolution_d.pop("is_area_specific_loss_rate", False)
        # if run_logger is not None:
        #     if self.is_area_specific_loss_rate:
        #         run_logger.info("(ANAGENETIC RANGE EVOLUTION) Area loss will be modeled on a per-area basis: area loss rates will be taken to be per lineage per area rather than per lineage")
        #     else:
        #         run_logger.info("(ANAGENETIC RANGE EVOLUTION) Area loss will be modeled on a per-area basis: area loss rates will be taken to be per lineage rather than per lineage per area")
        self.mean_area_loss_rate = anagenetic_range_evolution_d.pop("mean_area_loss_rate", 0.10)
        if run_logger is not None:
            run_logger.info("(ANAGENETIC RANGE EVOLUTION) Global area loss rate: {}".format(self.mean_area_loss_rate))
        if "lineage_area_loss_weight" in anagenetic_range_evolution_d:
            self.lineage_area_loss_weight_function = RateFunction.from_definition_dict(anagenetic_range_evolution_d.pop("lineage_area_loss_weight"), self.trait_types)
        else:
            self.lineage_area_loss_weight_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda **kwargs: 1.00",
                    description="fixed: 1.00",
                    trait_types=self.trait_types,
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC RANGE EVOLUTION) Setting lineage-specific area loss weight function: {desc}".format(
                desc=self.lineage_area_loss_weight_function.description,
                ))

        if anagenetic_range_evolution_d:
            raise TypeError("Unsupported keywords in anagenetic range evolution submodel: {}".format(anagenetic_range_evolution_d))

        # Cladogenetic range inheritance submodel
        cladogenesis_d = dict(model_definition.pop("cladogenetic_range_evolution", {}))
        self.cladogenesis_sympatric_subset_speciation_weight = float(cladogenesis_d.pop("sympatric_subset_speciation_weight", 1.0))
        self.cladogenesis_single_area_vicariance_speciation_weight = float(cladogenesis_d.pop("single_area_vicariance_speciation_weight", 1.0))
        self.cladogenesis_widespread_vicariance_speciation_weight = float(cladogenesis_d.pop("widespread_vicariance_speciation_weight", 1.0))
        self.cladogenesis_founder_event_speciation_weight = float(cladogenesis_d.pop("founder_event_speciation_weight", 0.0))
        if cladogenesis_d:
            raise TypeError("Unsupported keywords in cladogenetic range evolution submodel: {}".format(cladogenesis_d))
        if run_logger is not None:
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of sympatric subset speciation mode: {}".format(self.cladogenesis_sympatric_subset_speciation_weight))
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of single area vicariance speciation mode: {}".format(self.cladogenesis_single_area_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of widespread vicariance speciation mode: {}".format(self.cladogenesis_widespread_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC RANGE EVOLUTION) Base weight of founder event speciation ('jump dispersal') mode: {} (note that the effective weight of this event for each lineage is actually the product of this and the lineage-specific area gain weight)".format(self.cladogenesis_founder_event_speciation_weight))

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

    def new_geography(self, run_logger=None):
        geography = Geography(
                areas_definition=copy.deepcopy(self.geographical_definition),
                run_logger=run_logger)
        return geography

    def write_model(self, out):
        model_definition = collections.OrderedDict()
        model_definition["model_id"] = self.model_id
        model_definition["areas"] = self._example_geography.as_definition()
        model_definition["traits"] = self.trait_types.as_definition()
        model_definition["diversification"] = self.diversification_as_definition()
        model_definition["anagenetic_range_evolution"] = self.anagenetic_range_evolution_as_definition()
        model_definition["cladogenetic_range_evolution"] = self.cladogenetic_range_evolution_as_definition()
        model_definition["termination_conditions"] = self.termination_conditions_as_definition()
        json.dump(model_definition, out, indent=4, separators=(',', ': '))

    def diversification_as_definition(self):
        d = collections.OrderedDict()
        d["mean_birth_rate"] = self.mean_birth_rate
        d["lineage_birth_weight"] = self.lineage_birth_weight_function.as_definition()
        d["mean_death_rate"] = self.mean_death_rate
        d["lineage_death_weight"] = self.lineage_death_weight_function.as_definition()
        return d

    def anagenetic_range_evolution_as_definition(self):
        d = collections.OrderedDict()
        d["global_area_gain_rate"] = self.global_area_gain_rate
        d["lineage_area_gain_weight"] = self.lineage_area_gain_weight_function.as_definition()
        d["mean_area_loss_rate"] = self.mean_area_loss_rate
        d["lineage_area_loss_weight"] = self.lineage_area_loss_weight_function.as_definition()
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
