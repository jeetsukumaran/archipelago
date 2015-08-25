#! /usr/bin/env python

import sys
import os
import argparse
import collections
import re
import dendropy
from archipelago import model
from archipelago.profile import ArchipelagoProfiler

def create_traits_data_file(traits_data_file_path, tree, num_trait_types):
    lineage_trait_states = collections.OrderedDict()
    for trait_idx in range(num_trait_types):
        state_symbols = {}
        for taxon in tree.taxon_namespace:
            lineage_label = taxon.label
            state = taxon.traits_vector[trait_idx]
            try:
                symbol = state_symbols[state]
            except KeyError:
                symbol = ArchipelagoProfiler.GEIGER_STATE_SYMBOLS[len(state_symbols)]
                state_symbols[state] = symbol
            try:
                lineage_trait_states[lineage_label].append(symbol)
            except KeyError:
                lineage_trait_states[lineage_label] = [ symbol ]
    with open(traits_data_file_path, "w") as dataf:
        for lineage_label in lineage_trait_states:
            traits = ",".join(lineage_trait_states[lineage_label])
            dataf.write("{},{}\n".format(lineage_label, traits))
        dataf.flush()
        dataf.close()

def create_range_data_file(output_path, tree):
    sep = "\t"
    area_labels = ["a{}".format(idx+1) for idx, a in enumerate(tree.taxon_namespace[0].distribution_vector)]
    dataf = open(output_path, "w")
    # dataf.write("{num_lineages}\t{num_areas}\t({area_labels})\n".format(
    #     num_lineages=len(tree.taxon_namespace),
    #     num_areas=len(area_labels),
    #     area_labels=" ".join(area_labels),
    #     ))
    dataf.write("{num_lineages}\t{num_areas}\n".format(
        num_lineages=len(tree.taxon_namespace),
        num_areas=len(area_labels),
        ))
    for taxon in tree.taxon_namespace:
        incidences = [str(i) for i in taxon.distribution_vector]
        dataf.write("{}{}{}\n".format(taxon.label, sep, "".join(incidences)))
    dataf.flush()
    dataf.close()


def main():
    parser = argparse.ArgumentParser()
    source_options = parser.add_argument_group("Source Options")
    source_options.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to tree files.")
    source_options.add_argument("-f", "--format",
            dest="schema",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="Input data format (default: '%(default)s').")
    # output_options = parser.add_argument_group("Output Options")
    # output_options.add_argument("-o", "--output-filepath",
    #         default=None,
    #         help="Path to profile_results file (default: standard output)."
    #         )
    # output_options.add_argument("-l", "--labels",
    #         action="append",
    #         help="Labels to append to output (in format <FIELD-NAME>:value;)")
    # output_options.add_argument( "--no-header-row",
    #         action="store_true",
    #         default=False,
    #         help="Do not write a header row.")
    # output_options.add_argument( "--append",
    #         action="store_true",
    #         default=False,
    #         help="Append to output file if it already exists instead of overwriting.")

    args = parser.parse_args()
    source_filepaths = list(args.source_paths)
    tree_yielder = dendropy.Tree.yield_from_files(
            files=source_filepaths,
            schema=args.schema,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=True,
            )
    source_tree_idx = 0
    prev_source_idx = None
    for tree_idx, tree in enumerate(tree_yielder):
        if prev_source_idx != tree_yielder.current_file_index:
            prev_source_idx = tree_yielder.current_file_index
            source_tree_idx = 0
        else:
            source_tree_idx += 1
        sys.stderr.write("-archipelago- Source {source_idx} of {num_sources} ({source_filepath}), Tree #{tree_idx}\n".format(
                source_idx=tree_yielder.current_file_index+1,
                num_sources=len(source_filepaths),
                source_filepath=source_filepaths[tree_yielder.current_file_index],
                tree_idx=source_tree_idx+1,
                ))
        model.ArchipelagoModel.set_lineage_data(
                    tree=tree,
                    leaf_nodes_only=True,
                    lineage_data_source="node",
                    traits_filepath=None,
                    areas_filepath=None,
                    )
        tree.original_taxon_namespace = tree.taxon_namespace
        tree.taxon_namespace = dendropy.TaxonNamespace()
        for node_idx, node in enumerate(tree.leaf_node_iter()):
            node.original_taxon = node.taxon
            node.taxon = tree.taxon_namespace.new_taxon(label=node.label)
            node.taxon.traits_vector = node.traits_vector
            node.taxon.distribution_vector = node.distribution_vector
        output_file_stem = "{}.{:04d}".format(source_filepaths[tree_yielder.current_file_index], source_tree_idx+1)
        create_range_data_file(output_path=output_file_stem + ".ranges.txt", tree=tree)
        num_trait_types = len(tree.taxon_namespace[0].traits_vector)
        create_traits_data_file(
                traits_data_file_path=output_file_stem + ".traits.csv",
                tree=tree,
                num_trait_types=num_trait_types,
                )


if __name__ == "__main__":
    main()


