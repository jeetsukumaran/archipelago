#! /usr/bin/env python

import sys
import os
import argparse
from archipelago import model
from archipelago import profile

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
    source_options.add_argument("-m", "--model-file",
            default=None,
            help="Model file(s) for the input tree file(s)."
                 " Parameters of the model will be added to the"
                 " profile profile_results to facilitate analysis."
            )
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-o", "--output-file",
            default=None,
            help="Path to profile_results file (default: standard output)."
            )
    output_options.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")

    args = parser.parse_args()
    source_filepaths = list(args.source_paths)

    if args.model_file:
        archipelago_model = model.ArchipelagoModel.from_path(args.model_file)
    else:
        archipelago_model = None
    profiler = profile.ArchipelagoProfiler.from_option_args(args)
    profiles = []
    for source_idx, source_filepath in enumerate(source_filepaths):
        if not args.quiet:
            sys.stderr.write("-profiler- Source {source_idx} of {num_sources}: {source_filepath}\n".format(
                    source_idx=source_idx+1,
                    num_sources=len(source_filepaths),
                    source_filepath=source_filepath,
                    ))
        results = profiler.profile_trees_from_path(
                trees_filepath=source_filepath,
                schema=args.schema,
                generating_model=archipelago_model,
                )
        profiles.extend(results)
    profiler.write_profiles(
            dest=sys.stdout,
            profiles=profiles,
            suppress_headers=False)

if __name__ == "__main__":
    main()


