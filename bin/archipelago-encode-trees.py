#! /usr/bin/env python

import os
import sys
import argparse
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import dendropy
import archipelago
from archipelago import model

def symbols_desc(symbols):
    return ", ".join(["'{}'".format(s) for s in symbols])

def validate_and_set_taxon_data(
        data,
        data_description,
        allowed_symbols,
        attr_name,
        char_separator,
        _log):
    allowed_symbols = set(allowed_symbols)
    for taxon in data:
        seq = data[taxon]
        seq_symbols = seq.symbols_as_list()
        seq_symbol_set = set(seq_symbols)
        diff_symbols = seq_symbol_set - allowed_symbols
        if diff_symbols:
            sys.stderr.write("ERROR: Following symbols are not allowed in {} data:\n    {}\nAllowed symbols are:\n    {}\n".format(
                data_description,
                symbols_desc(diff_symbols),
                symbols_desc(allowed_symbols)))
        setattr(taxon, attr_name, char_separator.join(seq_symbols))

def read_taxon_data(
        trees,
        data_description,
        data_path,
        data_format,
        _log):
    pre_data_taxa = set([t.label for t in trees.taxon_namespace])
    pre_data_taxa_desc = ",".join(["'{}'".format(t) for t in pre_data_taxa])
    data_path = os.path.expanduser(os.path.expandvars(data_path))
    data = dendropy.StandardCharacterMatrix.get(
            path=data_path,
            schema=data_format,
            taxon_namespace=trees.taxon_namespace)
    post_data_taxa = set([t.label for t in trees.taxon_namespace])
    if pre_data_taxa != post_data_taxa:
        diff_taxa = post_data_taxa - pre_data_taxa
        diff_taxa_desc = ",".join(["'{}'".format(t) for t in diff_taxa])
        sys.stderr.write("ERROR: Following taxon names defined in {} data:\n    {}\nbut not in phylogeny:\n    {}\n".format(
            data_description,
            diff_taxa_desc,
            pre_data_taxa_desc))
        sys.exit(1)
    seq_length = None
    num_seqs = 0
    for taxon in data:
        num_seqs += 1
        pre_data_taxa.remove(taxon.label)
        seq = data[taxon]
        if seq_length is None:
            seq_length = len(seq)
        elif seq_length != len(seq):
            sys.stderr.write("Expecting sequence length of {} for {} data, but found {} for taxon '{}'\n".format(
                seq_length,
                data_description,
                len(seq),
                taxon.label))
    if pre_data_taxa:
        diff_taxa_desc = ",".join(["'{}'".format(t) for t in pre_data_taxa])
        sys.stderr.write("ERROR: Following taxon names defined in trees but not in {} data:\n    {}\n".format(
            data_description,
            diff_taxa_desc,
            ))
        sys.exit(1)
    _log("{} characters read for {} taxa for {} data from: {}".format(
        seq_length,
        num_seqs,
        data_description,
        data_path))
    return data

def main():
    parser = argparse.ArgumentParser(
            description="{} Data Encoder".format(archipelago.description())
            )
    parser.add_argument("-p", "--phylogeny",
            metavar="PHYLOGENY-FILE",
            help="Path to file defining the phylogeny.")
    parser.add_argument("-g", "--geography",
            metavar="GEOGRAPHY-FILE",
            help="Path to file defining the geography.")
    parser.add_argument("-t", "--traits",
            metavar="TRAITS-FILE",
            help="Path to file defining the traits.")
    parser.add_argument("-P", "--phylogeny-format",
            choices=("newick", "nexus"),
            default="newick",
            help="Format of the phylogeny file (default: '%(default)s').")
    parser.add_argument("--preserve-underscores",
            action="store_true",
            default=False,
            help="Do not convert unquoted underscores to blanks/spaces in labels of phylogeny.")
    parser.add_argument("-G", "--geography-format",
            choices=("phylip", "nexus",),
            default="phylip",
            help="Format of the geography file (default: '%(default)s').")
    parser.add_argument("-T", "--traits-format",
            choices=("phylip", "nexus",),
            default="phylip",
            help="Format of the traits file (default: '%(default)s').")
    parser.add_argument("-O", "--output-format",
            choices=("newick", "nexus"),
            default="newick",
            help="Format of the output data (default: '%(default)s').")
    parser.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Run in silent mode.")
    args = parser.parse_args()
    if args.phylogeny is None:
        sys.exit("Must specify path to phylogeny file using '-p'/'--phylogeny' option.")
    if args.geography is None:
        sys.exit("Must specify path to geography file using '-g'/'--geography' option.")
    # if args.traits is None:
    #     sys.exit("Must specify path to traits file using '-t'/'--traits' option.")
    if args.quiet:
        _log = lambda x: None
    else:
        _log = lambda x: sys.stderr.write("-archipelago- {}\n".format(x))

    taxon_namespace = dendropy.TaxonNamespace()

    trees_path = os.path.expanduser(os.path.expandvars(args.phylogeny))
    trees = dendropy.TreeList.get(
            path=trees_path,
            schema=args.phylogeny_format,
            taxon_namespace=taxon_namespace,
            preserve_underscores=args.preserve_underscores,
            )
    _log("{} tree(s), {} taxa read from: '{}'".format(len(trees), len(taxon_namespace), trees_path))

    geography_path = os.path.expanduser(os.path.expandvars(args.geography))
    geography_data = read_taxon_data(
            trees=trees,
            data_description="geography",
            data_path=geography_path,
            data_format=args.geography_format,
            _log=_log,
            )
    validate_and_set_taxon_data(
            data=geography_data,
            data_description="geography",
            allowed_symbols=("0", "1"),
            attr_name="geography_data_string",
            char_separator="",
            _log=_log)

    if args.traits is not None:
        traits_path = os.path.expanduser(os.path.expandvars(args.traits))
        traits_data = read_taxon_data(
                trees=trees,
                data_description="traits",
                data_path=traits_path,
                data_format=args.traits_format,
                _log=_log,
                )
        validate_and_set_taxon_data(
                data=traits_data,
                data_description="traits",
                allowed_symbols="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ",
                attr_name="traits_data_string",
                char_separator=model.Lineage._TRAITS_SEPARATOR,
                _log=_log)
    else:
        _log("No traits data specified")
        for taxon in trees.taxon_namespace:
            taxon.traits_data_string = model.Lineage._NULL_TRAITS

    for taxon in trees.taxon_namespace:
        taxon.old_label = taxon.label
        taxon.label = "{idx}{sep}{traits}{sep}{areas}".format(
                idx=taxon.label,
                sep=model.Lineage._LABEL_COMPONENTS_SEPARATOR,
                traits=taxon.traits_data_string,
                areas=taxon.geography_data_string,
                )

    out = sys.stdout
    trees.write(
        file=out,
        schema=args.output_format,
        unquoted_underscores=args.preserve_underscores)

if __name__ == "__main__":
    main()

