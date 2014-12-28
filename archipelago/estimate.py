#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import os
import re
import collections
import tempfile
import subprocess
import dendropy
from dendropy.utility import processio
from archipelago import model

class TraitEvolutionRateEstimator(object):

    # STATE_SYMBOLS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

    def __init__(self):
        self.tree_file = tempfile.NamedTemporaryFile()
        self.data_file = tempfile.NamedTemporaryFile()
        self.tree_file_name = self.tree_file.name
        self.data_file_name = self.data_file.name

    def estimate_trait_evolution_rate(self,
            trees,
            tree_results_map,
            trait_estimated_transition_rate_field_names,
            is_trees_decoded=False,
            is_suppressed_taxa=False,
            ):
        self.preprocess_tree_lineages(
            trees=trees,
            is_trees_decoded=is_trees_decoded,
            is_suppressed_taxa=is_suppressed_taxa)
        assert len(trait_estimated_transition_rate_field_names) == trees.num_trait_types
        for tree_idx, tree in enumerate(trees):
            with open(self.tree_file_name, "w") as tf:
                tree.write_to_stream(
                        tf,
                        "nexus",
                        translate_tree_taxa=True,
                        suppress_internal_node_labels=True,
                        suppress_internal_taxon_labels=True)
                tf.flush()
                tf.close()
            for trait_idx in range(trees.num_trait_types):
                rate = self._analyze(tree,
                        tree.lineage_trait_state_set_map[trait_idx])
                tree_results_map[tree][trait_estimated_transition_rate_field_names[trait_idx]] = rate
        self.restore_tree_lineages(trees)

    def preprocess_tree_lineages(self,
            trees,
            is_trees_decoded=False,
            is_suppressed_taxa=False):
        if not is_trees_decoded:
            model.ArchipelagoModel.decode_tree_lineages_from_labels(
                    trees=trees,
                    is_suppressed_taxa=is_suppressed_taxa,
                    leaf_nodes_only=True)
        sample_node = next(trees[0].leaf_node_iter())
        num_trait_types = len(sample_node.traits_vector)
        trees.num_trait_types = num_trait_types
        for tree_idx, tree in enumerate(trees):
            tree.lineage_trait_state_set_map = []
            for trait_idx in range(trees.num_trait_types):
                tree.lineage_trait_state_set_map.append({})
            tree._original_taxon_namespace = tree.taxon_namespace
            tree.taxon_namespace = dendropy.TaxonNamespace()
            for lineage_idx, lineage in enumerate(tree.leaf_node_iter()):
                lineage._original_taxon = lineage.taxon
                taxon_label = "s{}".format(lineage_idx+1)
                lineage._original_taxon = lineage.taxon
                lineage.taxon = tree.taxon_namespace.new_taxon(label=taxon_label)
                for trait_idx in range(trees.num_trait_types):
                    tree.lineage_trait_state_set_map[trait_idx][lineage.taxon.label] = str(lineage.traits_vector[trait_idx])
        return trees

    def _analyze(self,
            tree,
            taxon_label_state_map):
        symbols = set()
        with open(self.data_file_name, "w") as dataf:
            for taxon in tree.taxon_namespace:
                row = [taxon.label]
                states = taxon_label_state_map[taxon.label]
                symbols.update(states)
                row.append("".join(states))
                dataf.write("{}\n".format("\t".join(row)))
            dataf.flush()
            dataf.close()
        symbols = sorted(symbols)
        if len(symbols) < 2:
            return 0.0
        bt_commands = []
        bt_commands.append("1") # multstate
        bt_commands.append("1") # ml; 2 == mcmc
        if True: #len(name_to_symbol_map.SYMBOLS) > 7:
            bt_commands.append("restrictall q{}{}".format(
                symbols[0],
                symbols[1]))
        bt_commands.append("run")
        # bt_commands = "\n".join(bt_commands)
        p = subprocess.Popen(
                ["BayesTraits", self.tree_file_name, self.data_file_name],
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p, bt_commands)
        stdout = stdout.split("\n")
        result = dict(zip(stdout[-3].split("\t"), stdout[-2].split("\t")))
        rate = float(result['q{}{}'.format(symbols[0],symbols[1])])
        return rate

    def restore_tree_lineages(self, trees):
        for tree_idx, tree in enumerate(trees):
            tree.taxon_namespace = tree._original_taxon_namespace
            del tree._original_taxon_namespace
            lineage_trait_state_set_map = []
            for lineage_idx, lineage in enumerate(tree.leaf_node_iter()):
                lineage.taxon = lineage._original_taxon
                del lineage._original_taxon
        return trees

