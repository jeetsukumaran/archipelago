#! /usr/bin/env python

import os
import sys
import argparse
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import archipelago
from archipelago import simulate
from archipelago import model

def main():
    parser = argparse.ArgumentParser(
            description="{} Biogeographical Simulator".format(archipelago.description())
            )
    model_options = parser.add_argument_group("Simulation Model")
    model_options.add_argument("model_file",
            nargs="?",
            help="Path to file defining model dictionary")
    model_options.add_argument("-f", "--model-format",
            dest="model_file_schema",
            choices=["json", "python"],
            default=None,
            help="Format of model file.")

    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type=str,
        default='archipelago',
        metavar='OUTPUT-FILE-PREFIX',
        help="Prefix for output files (default: '%(default)s').")

    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-n", "--nreps",
            type=int,
            default=10,
            help="Number of replicates (default: %(default)s).")
    run_options.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    run_options.add_argument("--log-frequency",
            default=None,
            type=float,
            help="Frequency that background progress messages get written to the log (0: do not log informational messages).")
    run_options.add_argument("--file-logging-level",
            default="info",
            help="Message level threshold for file logs.")
    run_options.add_argument("--stderr-logging-level",
            default="info",
            help="Message level threshold for screen logs.")
    run_options.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")

    args = parser.parse_args()

    config_d = {}
    if args.model_file is None:
        model_definition = {}
        interpolate_missing_model_values = True
    else:
        model_definition = model.ArchipelagoModel.get_model_definition_from_path(
                filepath=args.model_file,
                schema=args.model_file_schema)
        interpolate_missing_model_values = True

    simulate.repeat_run(
            output_prefix=args.output_prefix,
            nreps=args.nreps,
            model_definition=model_definition,
            interpolate_missing_model_values=interpolate_missing_model_values,
            config_d=config_d,
            random_seed=args.random_seed,
            stderr_logging_level=args.stderr_logging_level,
            file_logging_level=args.file_logging_level,
            debug_mode=args.debug_mode)

if __name__ == "__main__":
    main()

