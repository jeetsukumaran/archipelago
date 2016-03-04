#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy
import re
import datetime

from archipelago import summarize
from archipelago import utility
from archipelago.utility import USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE

def parse_trait_states(labels):
    if not labels:
        return collections.OrderedDict()
    trait_states = []
    for label in labels:
        match = re.match(r"\s*(.*?)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError("Cannot parse fieldname and label (format required: fieldname:value): {}".format(label))
        fieldname, value = match.groups(0)
        # The trait states need to be an integer if
        # ArchipelagoModel.decode_label coerces the labels to integers.
        # The reason we do NOT want it parsed to an integer value
        # is to allow null traits 'NA', 'null', etc.
        trait_states.append( (int(fieldname), value,) )
    return trait_states

def parse_fieldname_and_value(labels):
    if not labels:
        return collections.OrderedDict()
    fieldname_value_map = collections.OrderedDict()
    for label in labels:
        match = re.match(r"\s*(.*?)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError("Cannot parse fieldname and label (format required: fieldname:value): {}".format(label))
        fieldname, value = match.groups(0)
        fieldname_value_map[fieldname] = value
    return fieldname_value_map

def main():
    parser = argparse.ArgumentParser()
    source_options = parser.add_argument_group("Source Options")
    source_options.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to simulated tree files.")
    summarization_options = parser.add_argument_group("Summarization Options")
    summarization_options.add_argument("-x", "--exclude-trait",
            action="append",
            help="Index of trait to exclude, with first trait indexed with value {}; multiple traits can be specified by repeating the option (e.g., '--exclude-trait {} --ingore-trait {}').".format(
                USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE,
                USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE,
                USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE+1,
                ))
    summarization_options.add_argument("-X", "--exclude-trait-state",
            action="append",
            help="States of traits to exclude, (in format <TRAIT-INDEX:STATE-INDEX>. Not that traits are {}-based indexed, and states are 0-based indexed. E.g. '--exclude-trait-state 1:0 --exclude-trait-state 1:3').".format(
                USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE,
                ))
    summarization_options.add_argument("--no-drop-trees-not-spanning-all-areas",
            action="store_true",
            default=False,
            help="Do NOT skip trees that do not span all areas.")
    output_options = parser.add_argument_group("Source Options")
    output_options.add_argument("-l", "--labels",
            action="append",
            help="Labels to append to output (in format <FIELD-NAME>:value;)")
    output_options.add_argument(
            "-o", "--output-filepath",
            default=None,
            help="Path to output file.")
    output_options.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")
    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress progress messages.")
    args = parser.parse_args()
    args.group_processed_trees_by_model = False
    if args.quiet:
        _progress_update_fn = None
    else:
        log_frequency_percentage = 1
        def _progress_update_fn(current_idx, total):
            if not (int(float(current_idx)/total * 10) % log_frequency_percentage):
                sys.stderr.write("  [{}] Tree {} of {}\n".format(datetime.datetime.now(), current_idx+1, total))
    if args.exclude_trait:
        trait_indexes_to_exclude = [int(i) - USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE for i in args.exclude_trait]
        assert -1 not in trait_indexes_to_exclude
    else:
        trait_indexes_to_exclude = None
    trait_states_to_exclude = parse_trait_states(args.exclude_trait_state)
    tree_summarizer = summarize.TreeSummarizer(
        drop_trees_not_spanning_all_areas=not args.no_drop_trees_not_spanning_all_areas,
        trait_indexes_to_exclude=trait_indexes_to_exclude,
        trait_states_to_exclude=trait_states_to_exclude
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
                    progress_update_fn=_progress_update_fn,
                    # lineage_data_source=lineage_data_source,
                    # traits_filepath=traits_filepath,
                    # areas_filepath=areas_filepath,
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
