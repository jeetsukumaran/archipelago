#! /usr/bin/env pythorvals

import math
from decimal import Decimal
import collections
import json
import dendropy
from archipelago import model

class EventLog(object):

    @staticmethod
    def prepare_tree_for_event_serialization(tree,
            time_to_add_to_extant_tips=0,
            is_set_taxa=True):
        old_taxon_namespace = tree.taxon_namespace
        tree.taxon_namespace = dendropy.TaxonNamespace()
        lineages_on_tree = set()
        model.set_node_times_ages_extancy(
                tree=tree,
                is_set_age=True,
                is_annotate=True,
                time_to_add_to_extant_tips=time_to_add_to_extant_tips)
        tree.is_rooted = True
        if is_set_taxa:
            node_label_fn = lambda x: x.encode_lineage(
                    set_label=False,
                    add_annotation=False,
                    exclude_supplemental_areas=False)
            for nd in tree:
                lineages_on_tree.add(nd)
                assert nd.taxon is None
                nd.taxon = tree.taxon_namespace.require_taxon(label=node_label_fn(nd))
            tree.encode_bipartitions()
            for nd in tree:
                nd.annotations["lineage_id"] = int(nd.bipartition)
        return old_taxon_namespace, lineages_on_tree

    @staticmethod
    def restore_tree_from_event_serialization(tree, old_taxon_namespace):
        for nd in tree:
            nd.taxon = None
            nd.annotations.drop()
        tree.taxon_namespace = old_taxon_namespace

    @staticmethod
    def compose_lineage_definitions(tree):
        lineage_defs = []
        for nd in tree.preorder_node_iter():
            is_leaf = len(nd._child_nodes) == 0
            lineage_definition = collections.OrderedDict([
                    ("lineage_id", int(nd.bipartition)),
                    ("lineage_label", nd.taxon.label),
                    ("lineage_parent_id", int(nd.parent_node.bipartition) if nd.parent_node is not None else None),
                    ("leafset_bitstring", nd.bipartition.leafset_as_bitstring()),
                    ("split_bitstring", nd.bipartition.split_as_bitstring()),
                    ("lineage_start_time", nd.parent_node.time if nd.parent_node else -1.0),
                    ("lineage_end_time", nd.time),
                    ("lineage_duration", nd.time - nd.parent_node.time if nd.parent_node else -1.0),
                    ("lineage_start_distribution_bitstring", nd.starting_distribution_bitstring),
                    ("lineage_end_distribution_bitstring", nd.ending_distribution_bitstring if nd._child_nodes else nd.distribution_bitstring(exclude_supplemental_areas=False)),
                    ("is_seed_node", nd.parent_node is None),
                    ("is_leaf", is_leaf),
                    ("is_extant_leaf", nd.is_extant and is_leaf),
                    ("age", nd.age),
            ])
            lineage_defs.append(lineage_definition)
        return lineage_defs

    @staticmethod
    def compose_tree_data(tree, max_event_time):
        tree_data = collections.OrderedDict()
        tree_data["newick"] = tree.as_string("newick",
                suppress_internal_node_labels=True,
                suppress_leaf_node_labels=False,
                suppress_annotations=False,
                )
        tree_data["seed_node_age"] = tree.seed_node.age
        tree_data["max_event_time"] = max_event_time
        tree_data["end_time"] = max(max_event_time, tree.seed_node.age)
        return tree_data

    def __init__(self):
        self.lineage_events = {}
        self.max_event_time = None
        self.lineages_on_tree = set()

    def register_event(self,
            lineage,
            event_time,
            event_type,
            **kwargs
            ):
        ev = {
            "lineage": lineage,
            "event_time": event_time,
            "event_type": event_type,
            "event_subtype": kwargs.pop("event_subtype"),
            "state_idx": kwargs.pop("state_idx"),
            "child0_lineage": kwargs.pop("child0_lineage"),
            "child1_lineage": kwargs.pop("child1_lineage"),
            }
        try:
            self.lineage_events[lineage].append(ev)
        except KeyError:
            self.lineage_events[lineage] = [ev]
        self.max_event_time = event_time if self.max_event_time is None else max(self.max_event_time, event_time)

    def map_to_new_lineage(self, old_lineage, new_lineage):
        if old_lineage in self.lineage_events:
            events = self.lineage_events[old_lineage]
            if new_lineage not in self.lineage_events:
                self.lineage_events[new_lineage] = []
            for event in events:
                event["lineage"] = new_lineage
                self.lineage_events[new_lineage].append(event)
            del self.lineage_events[old_lineage]
        for events in self.lineage_events.values():
            for event in events:
                assert event["lineage"] is not old_lineage
                for k, v in event.items():
                    if event[k] is old_lineage:
                        event[k] = new_lineage

    def write_histories(self,
            out,
            tree,
            ):
        old_taxon_namespace, self.lineages_on_tree = EventLog.prepare_tree_for_event_serialization(tree=tree)
        history_data = collections.OrderedDict()
        history_data["tree"] = EventLog.compose_tree_data(tree=tree, max_event_time=self.max_event_time)
        history_data["leaf_labels"] = self._compose_taxon_namespace(tree=tree)
        history_data["lineages"] = EventLog.compose_lineage_definitions(tree=tree)
        history_data["events"] = self._compose_event_entries()
        json.dump(history_data, out, indent=4, separators=(',', ': '))
        EventLog.restore_tree_from_event_serialization(tree=tree, old_taxon_namespace=old_taxon_namespace)

#     ### not used???
#     def compose_tree_string(self, tree):
#         return tree.as_string("newick")

    def _compose_taxon_namespace(self, tree):
        return [t.label for t in tree.taxon_namespace]

    def _compose_event_entries(self):
        events = []
        for lineage in self.lineage_events:
            assert lineage in self.lineages_on_tree, lineage.index
            for event in self.lineage_events[lineage]:
                # if event["event_type"].startswith("geography") and not event["event_subtype"].startswith("focal"):
                #     continue
                if lineage.parent_node:
                    assert event["event_time"] >= lineage.parent_node.time
                # assert event["event_time"] <= lineage.time, "{}, {}, {} ({})".format(event["event_time"], lineage.time, lineage, event["event_type"])
                ev_t = Decimal("{:0.8}".format(event["event_time"]))
                ln_t = Decimal("{:0.8}".format(lineage.time))
                assert ev_t <= ln_t, "{} <= {}: False, for event {} on lineage {}".format(
                        event["event_time"],
                        lineage.time,
                        event["event_type"],
                        lineage,
                        )
                if event["child0_lineage"] is not None:
                    assert event["child0_lineage"] in self.lineages_on_tree, event["child0_lineage"].index
                if event["child1_lineage"] is not None:
                    assert event["child1_lineage"] in self.lineages_on_tree, event["child1_lineage"].index
                d = collections.OrderedDict([
                    ("event_time", event["event_time"]),
                    ("lineage_id", int(event["lineage"].bipartition)),
                    ("event_type", event["event_type"]),
                    ("event_subtype", event["event_subtype"]),
                    ("state_idx", event["state_idx"]),
                    ("child0_lineage_id", int(event["child0_lineage"].bipartition) if event["child0_lineage"] is not None else None),
                    ("child1_lineage_id", int(event["child1_lineage"].bipartition) if event["child1_lineage"] is not None else None),
                    ])
                events.append(d)
        events.sort(key=lambda x:x["event_time"])
        return events
