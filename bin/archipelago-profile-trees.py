#! /usr/bin/env python

import sys
import os
import argparse
from archipelago import model
from archipelago import profile

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
    source_options.add_argument("-m", "--model-file",
            default=None,
            help="Model file(s) for the input tree file(s)."
                 " Parameters of the model will be added to the"
                 " profile profile_results to facilitate analysis."
            )

    profile_options = parser.add_argument_group("Profile Options")

    profile_options.add_argument("--no-estimate-pure-birth",
            action="store_true",
            default=False,
            help="Do NOT estimate birth rate under a pure-birth model.")
    profile_options.add_argument("--no-estimate-trait-transition",
            action="store_true",
            default=False,
            help="Do NOT estimate trait transition rate.")
    profile_options.add_argument("--no-estimate-area-transition",
            action="store_true",
            default=False,
            help="Do NOT estimate area transition rate.")
    profile_options.add_argument("--estimate-dec",
            action="store_true",
            default=False,
            help="Estimate parameters under Lagrange's DEC model (using BioGeoBears).")

    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-o", "--output-file",
            default=None,
            help="Path to profile_results file (default: standard output)."
            )
    output_options.add_argument( "--no-header-row",
            action="store_true",
            default=False,
            help="Do not write a header row.")

    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress progress messages.")
    run_options.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debug mode (work files will not be deleted).")
    run_options.add_argument("--ignore-estimation-errors",
            action="store_true",
            default=False,
            help="Ignore errors raised by estimation internally or by external programs")
    args = parser.parse_args()
    source_filepaths = list(args.source_paths)

    if args.model_file:
        archipelago_model = model.ArchipelagoModel.from_path(args.model_file)
    else:
        archipelago_model = None
    profiler = profile.ArchipelagoProfiler(
            is_estimate_pure_birth_rate=not args.no_estimate_pure_birth,
            is_estimate_trait_transition_rates=not args.no_estimate_trait_transition,
            is_estimate_area_transition_rates=not args.no_estimate_area_transition,
            is_estimate_dec=args.estimate_dec,
            quiet=args.quiet,
            fail_on_estimation_error=not args.ignore_estimation_errors,
            debug_mode=args.debug_mode,
            )
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


