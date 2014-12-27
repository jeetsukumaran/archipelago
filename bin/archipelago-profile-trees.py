#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from archipelago import profile
from archipelago import model

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
                 " profile results to facilitate analysis."
            )
    parser.add_argument("-q", "--quiet",
            action="store_true",
            dest="quiet",
            default=False,
            help="Suppress progress messages.")
    args = parser.parse_args()
    source_filepaths = list(args.source_paths)

    if args.model_file:
        archipelago_model = model.ArchipelagoModel.from_path(args.model_file)
    else:
        archipelago_model = None
    tree_profiler = profile.TreeProfiler()
    results = collections.OrderedDict()
    source_fieldnames = [
        "source.path",
        "tree.index",
    ]
    if archipelago_model is not None:
        model_fieldnames = [
            "num.areas",
            "num.focal.areas",
            "num.supplemental.areas",
            "lineage.birth.rate.definition",
            "lineage.birth.rate.description",
            "lineage.death.rate.definition",
            "lineage.death.rate.description",
            "lineage.dispersal.rate.definition",
            "lineage.dispersal.rate.description",
        ]
        for trait_idx, trait in enumerate(archipelago_model.trait_types):
            model_fieldnames.append("trait.{}.transition.rate".format(trait.label))

    else:
        model_fieldnames = [ ]
    data_fieldnames = [
        "pure.birth.rate",
    ]
    fieldnames = source_fieldnames + model_fieldnames + data_fieldnames

    out = sys.stdout
    for source_idx, source_filepath in enumerate(source_filepaths):
        if not args.quiet:
            sys.stderr.write("-profiler- Source {source_idx} of {num_sources}: {source_filepath}\n".format(
                    source_idx=source_idx+1,
                    num_sources=len(source_filepaths),
                    source_filepath=source_filepath,
                    ))
        trees = dendropy.TreeList.get_from_path(source_filepath, args.schema)
        for tree_idx, tree in enumerate(trees):
            if not args.quiet:
                sys.stderr.write("-profiler- Source {} of {}: Tree {} of {}\n".format(source_idx+1, len(source_filepaths), tree_idx+1, len(trees)))
            results[tree] = {}
            results[tree]["source.path"] = source_filepath
            results[tree]["tree.index"] = tree_idx
            if archipelago_model:
                results[tree]["num.areas"] = len(archipelago_model.geography.areas)
                results[tree]["num.focal.areas"] = len(archipelago_model.geography.focal_area_indexes)
                results[tree]["num.supplemental.areas"] = len(archipelago_model.geography.supplemental_area_indexes)
                results[tree]["lineage.birth.rate.definition"] = archipelago_model.lineage_birth_rate_function.definition_content
                results[tree]["lineage.birth.rate.description"] = archipelago_model.lineage_birth_rate_function.description
                results[tree]["lineage.death.rate.definition"] = archipelago_model.lineage_death_rate_function.definition_content
                results[tree]["lineage.death.rate.description"] = archipelago_model.lineage_death_rate_function.description
                results[tree]["lineage.dispersal.rate.definition"] = archipelago_model.lineage_dispersal_rate_function.definition_content
                results[tree]["lineage.dispersal.rate.description"] = archipelago_model.lineage_dispersal_rate_function.description
                for trait_idx, trait in enumerate(archipelago_model.trait_types):
                    results[tree]["trait.{}.transition.rate".format(trait.label)] = trait.transition_rate
        tree_profiler.estimate_pure_birth(trees, results)
    writer = csv.DictWriter(out,
            fieldnames=fieldnames)
    writer.writeheader()
    for row in results.values():
        writer.writerow(row)


if __name__ == "__main__":
    main()


