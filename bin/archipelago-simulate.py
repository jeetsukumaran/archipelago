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

def read_model_definition(filename):
    src = open(filename, "rb").read()
    return eval(src)

def main():
    parser = argparse.ArgumentParser(
            description="{} Biogeographical Simulator".format(archipelago.description())
            )

    run_options = parser.add_argument("model_file",
            help="Path to (Python) file defining model dictionary")

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
            type=int,
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

    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type=str,
        default='archipelago_run',
        metavar='OUTPUT-FILE-PREFIX',
        help="Prefix for output files (default: '%(default)s').")

    args = parser.parse_args()

    config_d = {}
    model_d = read_model_definition(args.model_file)

    simulate.repeat_run(
            output_prefix=args.output_prefix,
            nreps=args.nreps,
            model_d=model_d,
            config_d=config_d,
            random_seed=args.random_seed,
            stderr_logging_level=args.stderr_logging_level,
            file_logging_level=args.file_logging_level)

if __name__ == "__main__":
    main()

