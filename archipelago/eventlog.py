#! /usr/bin/env python

import collections
import json
import dendropy

class EventLog(object):

    def __init__(self):
        self.lineage_events = {}

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
            self.lineage_events[lineage] = ev

    def log_lineage_extinction(self, lineage):
        del self.lineage_events[lineage]

    def write_focal_areas_histories(self,
            out,
            tree,
            node_label_fn):
        history_data = collections.OrderedDict()
        history_data["lineages"] = self._compose_lineage_definitions(tree=tree, node_label_fn=node_label_fn)
        json.dump(history_data, out, indent=4, separators=(',', ': '))

    def _compose_lineage_definitions(self, tree, node_label_fn):
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
        tree.encode_bipartitions()
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
        for nd in tree:
            nd.taxon = None
        tree.taxon_namespace = old_taxon_namespace
        return lineage_defs


