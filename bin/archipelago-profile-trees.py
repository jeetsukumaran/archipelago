#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from archipelago import profile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to tree files.")
    parser.add_argument(
            "-o", "--output-prefix",
            default='archipelago-validation',
            help="Output prefix (default: '%(default)s').")
    parser.add_argument(
            "-x", "--exclude-first-island-as-continental-source-outside-of-analysis",
            action="store_true",
            default=False,
            help="Treat Island 0 as a continental source, and exclude it from analysis.")
    parser.add_argument("-f", "--format",
            dest="schema",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="input data format (default='%(default)s')")
    args = parser.parse_args()
    args.quiet = False
    # est = estimate.TraitEvolutionRateEstimator(
    #         exclude_first_island_as_continental_source_outside_of_analysis=args.exclude_first_island_as_continental_source_outside_of_analysis,
    #         drop_stunted_trees=True)
    for source_path in args.source_paths:
        trees = dendropy.TreeList.get_from_path(source_path, args.schema)



if __name__ == "__main__":
    main()


