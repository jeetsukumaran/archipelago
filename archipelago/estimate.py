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

class TraitTransitionRateEstimator(object):

    # STATE_SYMBOLS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    GEIGER_STATE_SYMBOLS = "123456789"

    def __init__(self):
        self.tree_file = tempfile.NamedTemporaryFile()
        self.data_file = tempfile.NamedTemporaryFile()
        self.work_file = tempfile.NamedTemporaryFile()
        self.output_file = tempfile.NamedTemporaryFile()
        self.tree_file_name = self.tree_file.name
        self.data_file_name = self.data_file.name
        self.work_file_name = self.work_file.name
        self.output_file_name = self.output_file.name

        # self.tree_file_name = "x1.nex"
        # self.data_file_name = "x1.txt"
        # self.work_file_name = "x1.in"
        # self.output_file_name = "x1.out"

    # def write_nexus(self, tree, taxon_label_state_map, filepath):
    #     from dendropy.interop import paup
    #     ds = dendropy.DataSet()
    #     ds.add_taxon_namespace(tree.taxon_namespace)
    #     trees = ds.new_tree_list(taxon_namespace=tree.taxon_namespace)
    #     trees.append(tree)
    #     cm = ds.new_char_matrix(
    #             char_matrix_type="standard",
    #             taxon_namespace=tree.taxon_namespace)
    #     for taxon_label in taxon_label_state_map:
    #         cm[taxon_label] = taxon_label_state_map[taxon_label]
    #     paup_block = []
    #     paup_block.append("BEGIN PAUP;")
    #     paup_block.append(paup.STANDARD_PREAMBLE + ";")
    #     paup_block.append("gett file='{}' storeBrlens=yes;".format(os.path.basename(filepath)))
    #     paup_block.append("END;")
    #     paup_block = "\n".join(paup_block)
    #     ds.write_to_path(
    #             filepath,
    #             "nexus",
    #             supplemental_blocks=[paup_block])

    def estimate_trait_transition_rate(self,
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

        self._estimate_trait_transition_rate_using_geiger(
            trees=trees,
            tree_results_map=tree_results_map,
            trait_estimated_transition_rate_field_names=trait_estimated_transition_rate_field_names)
        assert len(trait_estimated_transition_rate_field_names) == trees.num_trait_types

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
            tree.lineage_label_trait_state_set_map = []
            for trait_idx in range(trees.num_trait_types):
                tree.lineage_label_trait_state_set_map.append({})
            tree._original_taxon_namespace = tree.taxon_namespace
            tree.taxon_namespace = dendropy.TaxonNamespace()
            for node_idx, node in enumerate(tree.leaf_node_iter()):
                node._original_taxon = node.taxon
                taxon_label = "s{}".format(node_idx+1)
                node._original_taxon = node.taxon
                node.taxon = tree.taxon_namespace.new_taxon(label=taxon_label)
                for trait_idx in range(trees.num_trait_types):
                    tree.lineage_label_trait_state_set_map[trait_idx][node.taxon.label] = str(node.traits_vector[trait_idx])
        return trees

    def _estimate_trait_transition_rate_using_geiger(self,
            trees,
            tree_results_map,
            trait_estimated_transition_rate_field_names,
            ):
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
            lineage_row = {}
            for trait_idx in range(trees.num_trait_types):
                state_symbols = {}
                for taxon in tree.taxon_namespace:
                    state = tree.lineage_label_trait_state_set_map[trait_idx][taxon.label]
                    try:
                        symbol = state_symbols[state]
                    except KeyError:
                        symbol = TraitTransitionRateEstimator.GEIGER_STATE_SYMBOLS[len(state_symbols)]
                        state_symbols[state] = symbol
                    try:
                        lineage_row[taxon.label].append(symbol)
                    except KeyError:
                        lineage_row[taxon.label] = [ symbol ]
            dataf = open(self.data_file_name, "w")
            for taxon in tree.taxon_namespace:
                traits = ",".join(lineage_row[taxon.label])
                dataf.write("{},{}\n".format(taxon.label, traits))
            dataf.flush()
            dataf.close()
            rcmds = []
            rcmds.append("library(parallel, quietly=T)")
            rcmds.append("library(ape, quietly=T)")
            rcmds.append("library(geiger, quietly=T)")
            # rcmds.append("sink('{}')".format(self.output_file_name))
            rcmds.append("tree1 <- read.nexus('{}')".format(self.tree_file_name))
            rcmds.append("traits <- read.csv('{}', header=F, row.names=1)".format(self.data_file_name))
            for trait_idx in range(trees.num_trait_types):
                trait_var = "traitx"
                rcmds.append("{} <- traits[,{}]".format(trait_var, trait_idx+1))
                rcmds.append("names({}) <- row.names(traits)".format(trait_var))
                rcmds.append("m = fitDiscrete(tree1, {})".format(trait_var))
                rcmds.append(r"cat(c(m$opt$q12), sep='\n')")
            rcmds = "\n".join(rcmds)
            rfile = open(self.work_file_name, "w")
            rfile.write(rcmds + "\n")
            rfile.flush()
            rfile.close()
            shell_cmd = ["R",
                    "--vanilla",
                    "--no-save",
                    "--slave",
                    "--silent",
                    "-f",
                    self.work_file_name]
            p = subprocess.Popen(
                    shell_cmd,
                    stdout=subprocess.PIPE,
                    )
            stdout, stderr = processio.communicate(p)
            rows = [row.strip() for row in stdout.split("\n")]
            rows = [float(row) for row in rows if row]
            assert len(rows) == trees.num_trait_types
            assert len(rows) == len(trait_estimated_transition_rate_field_names)
            for rate, field_name in zip(rows, trait_estimated_transition_rate_field_names):
                tree_results_map[tree][field_name] = rate

    def _run_Geiger(self,
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

    def _estimate_trait_transition_rate_using_bayes_traits(self,
            trees,
            tree_results_map,
            trait_estimated_transition_rate_field_names,
            ):
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
                rate = self._run_BayesTraits(tree,
                        tree.lineage_label_trait_state_set_map[trait_idx])
                tree_results_map[tree][trait_estimated_transition_rate_field_names[trait_idx]] = rate

    def _run_BayesTraits(self,
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
            lineage_label_trait_state_set_map = []
            for lineage_idx, lineage in enumerate(tree.leaf_node_iter()):
                lineage.taxon = lineage._original_taxon
                del lineage._original_taxon
        return trees

