#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from archipelago import summarize

# library(adegenet)
# summary.df = read.table("processed/summary.txt", header=T)
# summary.df = na.omit(summary.df)
# groups = summary.df$dispersal.model
# cols.to.drop <- c(
#                   "dispersal.model",
#                   "birth.rate",
#                   "death.rate",
#                   "dispersal.rate",
#                   "niche.evolution.prob",
#                   "edges",
#                   "est.birth.rate",
#                   "length",
#                   "size"
#                   )
# predictors = summary.df[,!(names(summary.df) %in% cols.to.drop)]
# result = dapc(predictors, groups)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to simulated tree files.")
    parser.add_argument("-l", "--labels",
            action="append",
            help="Labels to append to output (in format <FIELD-NAME>:value;)")
    parser.add_argument(
            "-o", "--output-filepath",
            default='processed',
            help="Path to output file.")
    parser.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    parser.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")
    args = parser.parse_args()
    args.quiet = False
    args.group_processed_trees_by_model = False

    tree_summarizer = summarize.TreeSummarizer(
        drop_trees_not_occupying_all_islands=True,
        drop_trees_not_occupying_all_habitats=True,
        drop_stunted_trees=True,
    )
    param_keys = collections.OrderedDict()
    summaries = []
    output_root_dir = "."
    output_dir = output_root_dir
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    stats_fields = set()
    try:
        for source_idx, tree_filepath in enumerate(args.source_paths):
            if not args.quiet:
                sys.stderr.write("Processing job {} of {}: {}\n".format(source_idx+1, len(args.source_paths), tree_filepath))
            trees = dendropy.TreeList.get_from_path(
                    tree_filepath,
                    schema="newick",
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True,
                    )
            processed_trees, sub_stats_fields = tree_summarizer.summarize_trees(
                    trees,
                    summaries=summaries)
            stats_fields.update(sub_stats_fields)
    except KeyboardInterrupt:
        pass

    param_fields = list(param_keys.keys())
    stats_fields = sorted(list(stats_fields))
    all_fields = param_fields + stats_fields
    if args.output_filepath:
        out = open(args.output_filepath,
                "a" if args.append else "w")
    else:
        out = sys.stdout
    with out:
        writer = csv.DictWriter(out,
                fieldnames=all_fields,
                restval="NA",
                delimiter="\t")
        if not args.no_header_row:
            writer.writeheader()
        writer.writerows(summaries)

if __name__ == "__main__":
    main()
