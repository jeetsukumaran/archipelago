#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from archipelago import profile
from archipelago import model

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to tree files.")
    # parser.add_argument(
    #         "-o", "--output-prefix",
    #         default='archipelago-validation',
    #         help="Output prefix (default: '%(default)s').")
    parser.add_argument("-f", "--format",
            dest="schema",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="Input data format (default: '%(default)s').")
    parser.add_argument("-m", "--model-file",
            default=None,
            help="Model file(s) for the input tree file(s)."
                 " Parameters of the model will be added to the"
                 " profile results to facilitate analysis."
            )
    parser.add_argument("-o", "--output-file",
            default=None,
            help="Path to results file (default: standard output)."
            )
    parser.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress progress messages.")
    parser.add_argument( "--no-source-columns",
            action="store_true",
            default=False,
            help="Do not include columns indicating source path(s) and tree indexes.")
    parser.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    parser.add_argument( "--header-row-only",
            action="store_true",
            default=False,
            help="Write header row only and exit.")
    args = parser.parse_args()
    source_filepaths = list(args.source_paths)

    if args.model_file:
        archipelago_model = model.ArchipelagoModel.from_path(args.model_file)
    else:
        archipelago_model = None
    tree_profiler = profile.TreeProfiler()
    results = collections.OrderedDict()
    if not args.no_source_columns:
        source_fieldnames = [
            "source.path",
            "tree.index",
        ]
    else:
        source_fieldnames = []
    if archipelago_model is not None:
        model_fieldnames = [
            "num.areas",
            "num.focal.areas",
            "num.supplemental.areas",
            "lineage.birth.rate.definition",
            "lineage.birth.rate.description",
            "lineage.death.rate.definition",
            "lineage.death.rate.description",
            "lineage.dispersal.rate.definition",
            "lineage.dispersal.rate.description",
        ]

    else:
        model_fieldnames = [ ]
    data_fieldnames = [
        "num.tips",
        "root.age",
        "pure.birth.rate",
    ]
    trait_labels = []
    trait_estimated_transition_rate_field_names = []
    if archipelago_model is not None:
        for trait_idx, trait in enumerate(archipelago_model.trait_types):
            trait_labels.append(trait.label)
            model_fieldnames.append("trait.{}.true.transition.rate".format(trait.label))
            data_fieldnames.append("trait.{}.est.transition.rate".format(trait.label))
            trait_estimated_transition_rate_field_names.append(data_fieldnames[-1])
    else:
        trees = dendropy.TreeList.get_from_path(source_filepaths[0], args.schema)
        tree_profiler.diagnose_num_trait_types(trees)
        for trait_idx in range(trees.num_trait_types):
            trait_labels.append(str(trait_idx))
            # model_fieldnames.append("trait.{}.true.transition.rate".format(trait_idx))
            data_fieldnames.append("trait.{}.est.transition.rate".format(trait_idx))
            trait_estimated_transition_rate_field_names.append(data_fieldnames[-1])
    fieldnames = source_fieldnames + model_fieldnames + data_fieldnames

    if args.output_file is None or args.output_file == "-":
        out = sys.stdout
    else:
        out = open(args.output_file, "w")
    writer = csv.DictWriter(out,
            fieldnames=fieldnames)
    if not args.no_header_row:
        writer.writeheader()
    if args.header_row_only:
        sys.exit(0)
    for source_idx, source_filepath in enumerate(source_filepaths):
        if not args.quiet:
            sys.stderr.write("-profiler- Source {source_idx} of {num_sources}: {source_filepath}\n".format(
                    source_idx=source_idx+1,
                    num_sources=len(source_filepaths),
                    source_filepath=source_filepath,
                    ))
        trees = dendropy.TreeList.get_from_path(source_filepath, args.schema)
        for tree_idx, tree in enumerate(trees):
            # if not args.quiet:
            #     sys.stderr.write("-profiler- Source {} of {}: Tree {} of {}\n".format(source_idx+1, len(source_filepaths), tree_idx+1, len(trees)))
            results[tree] = {}
            if not args.no_source_columns:
                results[tree]["source.path"] = source_filepath
                results[tree]["tree.index"] = tree_idx
            if archipelago_model:
                results[tree]["num.areas"] = len(archipelago_model.geography.areas)
                results[tree]["num.focal.areas"] = len(archipelago_model.geography.focal_area_indexes)
                results[tree]["num.supplemental.areas"] = len(archipelago_model.geography.supplemental_area_indexes)
                results[tree]["lineage.birth.rate.definition"] = archipelago_model.lineage_birth_rate_function.definition_content
                results[tree]["lineage.birth.rate.description"] = archipelago_model.lineage_birth_rate_function.description
                results[tree]["lineage.death.rate.definition"] = archipelago_model.lineage_death_rate_function.definition_content
                results[tree]["lineage.death.rate.description"] = archipelago_model.lineage_death_rate_function.description
                results[tree]["lineage.dispersal.rate.definition"] = archipelago_model.lineage_dispersal_rate_function.definition_content
                results[tree]["lineage.dispersal.rate.description"] = archipelago_model.lineage_dispersal_rate_function.description
                for trait_idx, trait in enumerate(archipelago_model.trait_types):
                    results[tree]["trait.{}.true.transition.rate".format(trait.label)] = trait.transition_rate
            tree.calc_node_ages()
            results[tree]["root.age"] = tree.seed_node.age
            results[tree]["num.tips"] = len(list(nd for nd in tree.leaf_node_iter()))
            # for nd in tree.postorder_node_iter():
            #     if nd.edge.length is None:
            #         nd.edge.length = 0.0
            #     elif nd.edge.length < 0:
            #         if not args.quiet:
            #             sys.stderr.write("-profiler- Source {source_idx} of {num_sources}: Tree {tree_idx} of {num_trees}: setting negative branch length {brlen} for node {node} to 0.0\n".format(
            #                     source_idx=source_idx+1,
            #                     num_sources=len(source_filepaths),
            #                     tree_idx=tree_idx+1,
            #                     num_trees=len(trees),
            #                     brlen=nd.edge.length,
            #                     node=str(nd),
            #                     ))
            #         nd.edge.length = 0.0
        tree_profiler.estimate_pure_birth(
                trees=trees,
                tree_results_map=results,
                )
        tree_profiler.estimate_trait_transition_rates(
                trees=trees,
                tree_results_map=results,
                trait_estimated_transition_rate_field_names=trait_estimated_transition_rate_field_names,
                is_trees_decoded=False,
                is_suppressed_taxa=False,
                )
    for row in results.values():
        writer.writerow(row)

if __name__ == "__main__":
    main()


