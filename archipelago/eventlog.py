#! /usr/bin/env python

import collections
import json
import dendropy

class EventLog(object):

    def __init__(self):
        self.lineage_events = {}
        self.max_event_time = None

    def register_event(self,
            lineage,
            event_time,
            event_type,
            event_subtype,
            state_idx,
            child0_lineage,
            child1_lineage,
            ):
        if event_type == "extinction":
            self.log_lineage_extinction(lineage)
        else:
            ev = {
                "lineage": lineage,
                "event_time": event_time,
                "event_type": event_type,
                "event_subtype": event_subtype,
                "state_idx": state_idx,
                "child0_lineage": child0_lineage,
                "child1_lineage": child1_lineage,
                }
            try:
                self.lineage_events[lineage].append(ev)
            except KeyError:
                self.lineage_events[lineage] = [ev]
            self.max_event_time = event_time if self.max_event_time is None else max(self.max_event_time, event_time)

    def log_lineage_extinction(self, lineage):
        del self.lineage_events[lineage]

    def write_focal_areas_histories(self,
            out,
            tree,
            node_label_fn):
        old_taxon_namespace = self._prepare_tree_for_event_serialization(
                tree=tree,
                node_label_fn=node_label_fn)
        tree.calc_node_ages()
        history_data = collections.OrderedDict()
        history_data["tree"] = self._compose_tree_data(tree=tree)
        history_data["leaf_labels"] = self._compose_taxon_namespace(tree=tree)
        history_data["lineages"] = self._compose_lineage_definitions(tree=tree)
        history_data["events"] = self._compose_event_entries()
        json.dump(history_data, out, indent=4, separators=(',', ': '))
        self._restore_tree_from_event_serialization(tree=tree, old_taxon_namespace=old_taxon_namespace)

    def _prepare_tree_for_event_serialization(self, tree, node_label_fn):
        old_taxon_namespace = tree.taxon_namespace
        tree.taxon_namespace = dendropy.TaxonNamespace()
        for nd in tree:
            if nd.parent_node:
                nd.time = nd.parent_node.time + nd.edge.length
            else:
                nd.time = 0
            if nd.is_leaf():
                assert nd.taxon is None
                nd.taxon = tree.taxon_namespace.require_taxon(label=node_label_fn(nd))
        tree.is_rooted = True
        tree.encode_bipartitions()
        # for nd in tree:
        #     nd.annotations["lineage_id"] = str(int(nd.bipartition))
        return old_taxon_namespace

    def _restore_tree_from_event_serialization(self, tree, old_taxon_namespace):
        for nd in tree:
            nd.taxon = None
        tree.taxon_namespace = old_taxon_namespace

    def _compose_tree_data(self, tree):
        tree_data = collections.OrderedDict()
        # tree_data["newick"] = tree.as_string("newick", suppress_annotations=False)
        tree_data["newick"] = tree.as_string("newick")
        tree_data["seed_node_age"] = tree.seed_node.age
        tree_data["max_event_time"] = self.max_event_time
        tree_data["end_time"] = max(self.max_event_time, tree.seed_node.age)
        return tree_data

    def _compose_taxon_namespace(self, tree):
        return [t.label for t in tree.taxon_namespace]

    def _compose_lineage_definitions(self, tree):
        lineage_defs = []
        for nd in tree.preorder_node_iter():
            lineage_definition = collections.OrderedDict([
                    ("lineage_id", int(nd.bipartition)),
                    ("lineage_parent_id", int(nd.parent_node.bipartition) if nd.parent_node is not None else None),
                    ("leafset_bitstring", nd.bipartition.leafset_as_bitstring()),
                    ("split_bitstring", nd.bipartition.split_as_bitstring()),
                    ("lineage_start_time", nd.parent_node.time if nd.parent_node else -1.0),
                    ("lineage_end_time", nd.time),
                    ("lineage_start_distribution_bitstring", nd.starting_distribution_bitstring),
                    ("lineage_end_distribution_bitstring", nd.ending_distribution_bitstring if nd._child_nodes else nd.distribution_bitstring(exclude_supplemental_areas=True)),
                    ("is_seed_node", nd.parent_node is None),
                    ("is_leaf", len(nd._child_nodes) == 0),
            ])
            lineage_defs.append(lineage_definition)
        return lineage_defs

    def _compose_event_entries(self):
        events = []
        for lineage in self.lineage_events:
            for event in self.lineage_events[lineage]:
                if lineage.parent_node:
                    assert event["event_time"] >= lineage.parent_node.time
                assert event["event_time"] <= lineage.time, "{}, {}, {} ({})".format(event["event_time"], lineage.time, lineage, event["event_type"])
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
