#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from archipelago import profile
from archipelago import utility

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
            action="append",
            help="Model file(s) for each of the input tree file(s). If specified, this"
                 " option must be specified exactly once *or* once and exactly once for each"
                 " tree file given as input. Parameters of the model will be added to the"
                 " profile results to facilitate analysis. If only one model file is specified"
                 " then it will be used for all the tree files. If multiple model files are specified"
                 " then the number should match the number of source files and the order that they"
                 " should be matched to the source files."
            )
    parser.add_argument("-q", "--quiet",
            action="store_true",
            dest="quiet",
            default=False,
            help="Suppress progress messages.")
    args = parser.parse_args()
    tree_profiler = profile.TreeProfiler()
    results = collections.OrderedDict()
    fieldnames = [
        "source.path",
        "tree.index",
        "pure.birth.rate",
    ]
    out = sys.stdout
    source_filepaths = list(args.source_paths)
    if not args.model_file:
        model_filepaths = [None] * len(source_filepaths)
    elif len(args.model_file) == 1:
        model_filepaths = [args.model_file[0]] * len(source_filepaths)
    elif len(args.model_file) == len(source_filepaths):
        model_filepaths = list(args.model_file)
    else:
        sys.exit("{} source file paths specified; exactly none, one or {} model filepaths required, but {} given".format(len(source_filepaths), len(source_filepaths), len(args.model_file)))
    for source_idx, (source_filepath, model_filepath) in enumerate(zip(source_filepaths, model_filepaths)):
        if model_filepath is None:
            model_file_desc = ""
        else:
            model_file_desc = " (model: {})".format(model_filepath)
        #     if model_filepath.endswith(".py"):
        #         model_d = utility.load_model_from_python_path(model_filepath)
        #     elif model_filepath.endswith(".json"):
        #         model_d = utility.load_model_from_json_path(model_filepath)
        #     else:
        #         raise ValueError("Unrecognized file type: {}".format(model_filepath))
        if not args.quiet:
            sys.stderr.write("-profiler- Source {source_idx} of {num_sources}{model_file_desc}: {source_filepath}\n".format(
                    source_idx=source_idx+1,
                    num_sources=len(source_filepaths),
                    model_file_desc=model_file_desc,
                    source_filepath=source_filepath,
                    ))
        trees = dendropy.TreeList.get_from_path(source_filepath, args.schema)
        for tree_idx, tree in enumerate(trees):
            if not args.quiet:
                sys.stderr.write("-profiler- Source {} of {}: Tree {} of {}\n".format(source_idx+1, len(source_filepaths), tree_idx+1, len(trees)))
            results[tree] = {}
            results[tree]["source.path"] = source_filepath
            results[tree]["tree.index"] = tree_idx
        tree_profiler.estimate_pure_birth(trees, results)
    writer = csv.DictWriter(out,
            fieldnames=fieldnames)
    writer.writeheader()
    for row in results.values():
        writer.writerow(row)


if __name__ == "__main__":
    main()


