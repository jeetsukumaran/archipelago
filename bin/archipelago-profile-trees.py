#! /usr/bin/env python

import sys
import os
import argparse
import collections
import re
from archipelago import model
from archipelago import profile
from archipelago import utility

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
    parser = argparse.ArgumentParser(
            parents=[
                profile.ArchipelagoProfiler.get_profile_options_parser(),
                ],
            )
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
    source_options.add_argument("--no-preserve-underscores",
            action="store_true",
            default=False,
            help="Convert unquoted underscores to spaces, as dictated by the Newick/NEXUS standards.")
    source_options.add_argument("-r", "--range-data-file",
            default=None,
            help="Lineage range data file (in Phylip/BayArea/BioGeoBEARS format)."
                 " If not provided, lineage range data must be encoded in tip labels."
                "  If provided, lienage traits data file must be provided as well using"
                "  '-t'/'--traits-data-file' option."
            )
    source_options.add_argument("-t", "--traits-data-file",
            default=None,
            help="Lineage traits data file (in Phylip format)."
                 " If not provided, lineage traits data must be encoded in tip labels."
                "  If provided, lienage range data file must be provided as well using"
                "  '-t'/'--range-data-file' option."
            )
    source_options.add_argument("-m", "--model-file",
            default=None,
            help="Model file(s) for the input tree file(s)."
                 " Parameters of the model will be added to the"
                 " profile profile_results to facilitate analysis."
            )
    source_options.add_argument("-M", "--model-file-type",
            default="json",
            choices=("json", "python"),
            help="Model file format (default: %(default)s)."
            )
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-o", "--output-filepath",
            default=None,
            help="Path to profile_results file (default: standard output)."
            )
    output_options.add_argument("-l", "--labels",
            action="append",
            help="Labels to append to output (in format <FIELD-NAME>:value;)")
    output_options.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")

    args = parser.parse_args()
    extra_fields = parse_fieldname_and_value(args.labels)
    source_filepaths = list(args.source_paths)

    if args.range_data_file is not None and args.traits_data_file is None:
        sys.exit("If specifying '--range-data-file', must also specify '--traits-data-file'")
    elif args.range_data_file is None and args.traits_data_file is not None:
        sys.exit("If specifying '--traits-data-file', must also specify '--range-data-file'")
    elif args.range_data_file is None and args.traits_data_file is not None:
        lineage_data_source="filepath"
    else:
        lineage_data_source="node"

    if args.model_file:
        if args.model_file_type == "json":
            model_file_type = "json-filepath"
        elif args.model_file_type == "python":
            model_file_type = "python-filepath"
        else:
            raise ValueError(args.model_file_type)
        archipelago_model = model.ArchipelagoModel.create(
                model_definition_source=args.model_file,
                model_definition_type=model_file_type,
                )
    else:
        archipelago_model = None
    profiler = profile.ArchipelagoProfiler.from_option_args(args)
    profiles = []
    for source_idx, source_filepath in enumerate(source_filepaths):
        if args.verbosity >= 2:
            sys.stderr.write("-profiler- Source {source_idx} of {num_sources}: {source_filepath}\n".format(
                    source_idx=source_idx+1,
                    num_sources=len(source_filepaths),
                    source_filepath=source_filepath,
                    ))
        results = profiler.profile_trees_from_path(
                trees_filepath=source_filepath,
                schema=args.schema,
                preserve_underscores=not args.no_preserve_underscores,
                generating_model=archipelago_model,
                lineage_data_source=lineage_data_source,
                traits_filepath=args.traits_data_file,
                areas_filepath=args.range_data_file,
                )
        if extra_fields:
            for r in results:
                r.update(extra_fields)
        profiles.extend(results)
    out = utility.open_output_file_for_csv_writer(
            filepath=args.output_filepath,
            append=args.append)
    profiler.write_profiles(
            dest=out,
            profiles=profiles,
            suppress_headers=False)

if __name__ == "__main__":
    main()


