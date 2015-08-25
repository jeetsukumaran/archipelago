#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy
import re

from archipelago import summarize
from archipelago import utility

def parse_fieldname_and_value(labels):
    if not labels:
        return collections.OrderedDict()
    fieldname_value_map = collections.OrderedDict()
    for label in labels:
        match = re.match(r"\s*(.*)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError("Cannot parse fieldname and label (format required: fieldname:value): {}".format(label))
        fieldname, value = match.groups(0)
        fieldname_value_map[fieldname] = value
    return fieldname_value_map

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
            default=None,
            help="Path to output file.")
    parser.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    parser.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")
    parser.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress progress messages.")
    args = parser.parse_args()
    args.group_processed_trees_by_model = False
    tree_summarizer = summarize.TreeSummarizer(
        drop_trees_not_occupying_all_areas=True,
        drop_trees_not_occupying_all_traits=True,
        drop_stunted_trees=True,
    )
    summary_results = []
    output_root_dir = "."
    output_dir = output_root_dir
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    extra_fields = parse_fieldname_and_value(args.labels)
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
            processed_trees, sub_stats_fields, sub_results = tree_summarizer.summarize_trees(
                    trees,
                    # lineage_data_source=lineage_data_source,
                    # traits_filepath=traits_filepath,
                    # distribution_filepath=distribution_filepath,
                    )
            stats_fields.update(sub_stats_fields)
            if extra_fields:
                for r in sub_results:
                    r.update(extra_fields)
            summary_results.extend(sub_results)
    except KeyboardInterrupt:
        pass

    stats_fields = sorted(list(stats_fields))
    all_fields = list(extra_fields.keys()) + stats_fields
    out = utility.open_output_file_for_csv_writer(
            filepath=args.output_filepath,
            append=args.append)
    with out:
        writer = csv.DictWriter(
                out,
                fieldnames=all_fields,
                restval="NA",
                delimiter=",",
                lineterminator=os.linesep,
                )
        if not args.no_header_row:
            writer.writeheader()
        writer.writerows(summary_results)

if __name__ == "__main__":
    main()
