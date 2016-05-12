#! /usr/bin/env python

import dendropy
import json

class EventLog(object):

    def __init__(self):
        self.lineage_events = {}

    def register_event(self,
            lineage,
            event_time,
            event_type,
            event_subtype,
            area_idx,
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
                "area_idx": area_idx,
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
        d = {}
        d["lineages"] = []
        for nd in tree.postorder_node_iter():
            lineage_definition = {
                    "lineage_id": int(nd.bipartition),
                    "lineage_parent_id": int(nd.parent_node.bipartition) if nd.parent_node is not None else None,
                    "leafset_bitstring": nd.bipartition.leafset_as_bitstring(),
                    "split_bitstring": nd.bipartition.split_as_bitstring(),
                    "lineage_start_time": nd.parent_node.time if nd.parent_node else -1.0,
                    "lineage_end_time": nd.time,
                    "lineage_start_distribution_bitstring": None,
                    "lineage_end_distribution_bitstring": None,
                    "is_seed_node": nd.parent_node is None,
                    "is_leaf": len(nd._child_nodes) == 0,
            }

