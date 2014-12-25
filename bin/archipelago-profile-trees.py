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
    # parser.add_argument(
    #         "-o", "--output-prefix",
    #         default='archipelago-validation',
    #         help="Output prefix (default: '%(default)s').")
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
    tree_profiler = profile.TreeProfiler()
    results = collections.OrderedDict()
    fieldnames = [
        "source.path",
        "tree.index",
        "pure.birth.rate",
    ]
    out = sys.stdout
    for source_idx, source_path in enumerate(args.source_paths):
        if not args.quiet:
            sys.stderr.write("-profiler- Source {} of {}: {}\n".format(source_idx+1, len(args.source_paths), source_path,))
        trees = dendropy.TreeList.get_from_path(source_path, args.schema)
        for tree_idx, tree in enumerate(trees):
            if not args.quiet:
                sys.stderr.write("-profiler- Source {} of {}: Tree {} of {}\n".format(source_idx+1, len(args.source_paths), tree_idx+1, len(trees)))
            results[tree] = {}
            results[tree]["source.path"] = source_path
            results[tree]["tree.index"] = tree_idx
        tree_profiler.estimate_pure_birth(trees, results)
    writer = csv.DictWriter(out,
            fieldnames=fieldnames)
    writer.writeheader()
    for row in results.values():
        writer.writerow(row)


if __name__ == "__main__":
    main()


