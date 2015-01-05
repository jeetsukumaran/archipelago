#! /usr/bin/env python

import sys
import os
import argparse
from archipelago import model
from archipelago import profile

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
                 " profile profile_results to facilitate analysis."
            )
    parser.add_argument("-o", "--output-file",
            default=None,
            help="Path to profile_results file (default: standard output)."
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
    profiler = profile.ArchipelagoProfiler()
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
                generating_model=archipelago_model)
        profiles.extend(results)
    profiler.write_profiles(
            dest=sys.stdout,
            profiles=profiles,
            suppress_headers=False)

if __name__ == "__main__":
    main()


