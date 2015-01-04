#! /usr/bin/env python

import sys
import os
import argparse

import sys
import os
import csv
import collections
import tempfile
import subprocess
try:
    import lagrange
    IS_LAGRANGE_AVAILABLE = True
except ImportError:
    IS_LAGRANGE_AVAILABLE = False
import dendropy
from dendropy.model import birthdeath
from dendropy.utility import processio
from archipelago import model

LAGRANGE_CONFIGURATION_TEMPLATE = """\
#!/usr/bin/env python
import os
import lagrange
data = \"\"\"\
### begin data
{{
 'area_adjacency': {area_adjacency},
 'area_dispersal': {area_dispersal},
 'area_labels': {area_labels},
 'base_rates': '__estimate__',
 'dispersal_durations': [10000.0],
 'dm_symmetric_entry': True,
 'excluded_ranges': [],
 'lagrange_version': '20130526',
 'max_range_size': {max_range_size},
 'model_name': '{model_name}',
 'newick_trees': [{{'included': [],
                   'name': 'Tree0',
                   'newick': '{newick_tree_string}',
                   'root_age': {root_age}}}],
 'ranges': {ranges},
 'taxa': {taxon_name_list},
 'taxon_range_data': {taxon_range_data},
}}
### end data
\"\"\"

i = 0
while 1:
    if not i:
        outfname = "{model_name}.results.txt"
    else:
        outfname = "{model_name}.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = open(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)

"""

class ArchipelagoProfiler(object):

    GEIGER_STATE_SYMBOLS = "123456789"

    def __init__(self,
            quiet=False,
            fail_on_estimation_error=True,
            debug_mode=False):
        self.quiet = quiet
        self.fail_on_estimation_error = fail_on_estimation_error
        self.debug_mode = debug_mode
        self.lagrange_estimation = IS_LAGRANGE_AVAILABLE
        if self.debug_mode:
            self.tree_file_name = "profiler.tree.nexus"
            self.traits_data_file_name = "profiler.traits.data.txt"
            self.geography_data_file_name = "profiler.geography.data.txt"
            self.commands_file_name = "profiler.commands.txt"
        else:
            self.tree_file = tempfile.NamedTemporaryFile()
            self.traits_data_file = tempfile.NamedTemporaryFile()
            self.geography_data_file = tempfile.NamedTemporaryFile()
            self.commands_file = tempfile.NamedTemporaryFile()
            self.output_file = tempfile.NamedTemporaryFile()
            self.tree_file_name = self.tree_file.name
            self.traits_data_file_name = self.traits_data_file.name
            self.geography_data_file_name = self.geography_data_file.name
            self.commands_file_name = self.commands_file.name

    def send_message(self, *args, **kwargs):
        if self.quiet:
            return
        sep = kwargs.pop("sep", " ")
        s = sep.join(str(a) for a in args)
        sys.stderr.write("-profiler- {}\n".format(s))

    def write_profiles(self,
            dest,
            profiles,
            suppress_headers=False):
        fieldnames = list(profiles[0].keys())
        writer = csv.DictWriter(dest, fieldnames=fieldnames)
        if not suppress_headers:
            writer.writeheader()
        writer.writerows(profiles)

    def profile_trees_from_path(self,
            trees_filepath,
            schema="newick",
            generating_model=None):
        trees = dendropy.TreeList.get_from_path(
                trees_filepath,
                schema=schema,
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        trees.tree_filepath = trees_filepath
        profiles = self.profile_trees(
                trees=trees,
                generating_model=generating_model,
                is_lineages_decoded=False,
                decode_lineages_from="node",)
        return profiles

    def profile_trees(self,
            trees,
            generating_model=None,
            is_lineages_decoded=False,
            decode_lineages_from="node"):
        profiles = []
        for tree_idx, tree in enumerate(trees):
            if hasattr(trees, "tree_filepath"):
                tree.tree_filepath = trees.tree_filepath
                tree.tree_offset = tree_idx
            r = self.profile_tree(
                    tree=tree,
                    generating_model=generating_model,
                    is_lineages_decoded=is_lineages_decoded)
            profiles.append(r)
        return profiles

    def profile_tree(self,
            tree,
            generating_model=None,
            is_lineages_decoded=False,
            decode_lineages_from="node"):
        if not is_lineages_decoded:
            model.ArchipelagoModel.decode_tree_lineages_from_labels(
                    tree=tree,
                    leaf_nodes_only=True,
                    encoded_source=decode_lineages_from)

        # intitialize profile_results
        profile_results = collections.OrderedDict()

        # prepopulate with params that are available
        if generating_model is not None:
            self.store_generating_model_parameters(generating_model, profile_results)

        # basic profile
        if hasattr(tree, "tree_filepath"):
            profile_results["tree.filepath"] = tree.tree_filepath
            if hasattr(tree, "tree_offset"):
                profile_results["tree.offset"] = tree.tree_offset
        profile_results["num.tips"] = len(list(nd for nd in tree.leaf_node_iter()))
        tree.calc_node_ages()
        root_age = tree.seed_node.age
        profile_results["root.age"] = root_age

        # estimate birth rate
        self.estimate_pure_birth_rate(tree, profile_results)

        # set up for external processing
        ## fix negative edge lengths
        self.preprocess_edge_lengths(tree)
        ## set up taxa
        self.preprocess_tree_taxa(tree)
        ## store tree
        self.create_working_tree_data(tree)

        ## create traits
        trait_names = self.create_traits_data(
                tree=tree,
                generating_model=generating_model)

        ## process traits
        self.estimate_trait_evolution_rates(
                tree=tree,
                profile_results=profile_results,
                trait_names=trait_names)

        # process areas
        self.estimate_pure_dispersal_rate(
                tree=tree,
                profile_results=profile_results)
        # self.estimate_lagrange_rates(
        #         tree=tree,
        #         profile_results=profile_results)

        # clean up
        self.restore_tree_taxa(tree)
        self.restore_edge_lengths(tree)

        # return
        return profile_results

    def store_generating_model_parameters(self, generating_model, profile_results):
        profile_results["num.areas"] = len(generating_model.geography.areas)
        profile_results["num.focal.areas"] = len(generating_model.geography.focal_area_indexes)
        profile_results["num.supplemental.areas"] = len(generating_model.geography.supplemental_area_indexes)
        profile_results["lineage.birth.rate.definition"] = generating_model.lineage_birth_rate_function.definition_content
        profile_results["lineage.birth.rate.description"] = generating_model.lineage_birth_rate_function.description
        profile_results["lineage.death.rate.definition"] = generating_model.lineage_death_rate_function.definition_content
        profile_results["lineage.death.rate.description"] = generating_model.lineage_death_rate_function.description
        profile_results["lineage.dispersal.rate.definition"] = generating_model.lineage_dispersal_rate_function.definition_content
        profile_results["lineage.dispersal.rate.description"] = generating_model.lineage_dispersal_rate_function.description
        for trait_idx, trait in enumerate(generating_model.trait_types):
            profile_results["trait.{}.true.transition.rate".format(trait.label)] = trait.transition_rate
        return profile_results

    def estimate_pure_birth_rate(self, tree, profile_results):
        try:
            bdfit = birthdeath.fit_pure_birth_model_to_tree(tree)
        except ZeroDivisionError:
            rate = 0.0
        else:
            rate = bdfit["birth_rate"]
        profile_results["pure.birth.rate"] = rate

    def estimate_trait_evolution_rates(self,
            tree,
            profile_results,
            trait_names):
        rcmds = []
        rcmds.append("library(parallel, quietly=T)")
        rcmds.append("library(ape, quietly=T)")
        rcmds.append("library(geiger, quietly=T)")
        rcmds.append("tree1 <- read.nexus('{}')".format(self.tree_file_name))
        rcmds.append("traits <- read.csv('{}', header=F, row.names=1)".format(self.traits_data_file_name))
        for trait_idx, trait_name in enumerate(trait_names):
            trait_var = "trait{}".format(trait_idx)
            rcmds.append("{} <- traits[,{}]".format(trait_var, trait_idx+1))
            rcmds.append("names({}) <- row.names(traits)".format(trait_var))
            rcmds.append("m = fitDiscrete(tree1, {})".format(trait_var))
            rcmds.append(r"cat(c(m$opt$q12), sep='\n')")
        rcmds = "\n".join(rcmds)
        rfile = open(self.commands_file_name, "w")
        rfile.write(rcmds + "\n")
        rfile.flush()
        rfile.close()
        shell_cmd = ["R",
                "--vanilla",
                "--no-save",
                "--slave",
                "--silent",
                "-f",
                self.commands_file_name]
        p = subprocess.Popen(
                shell_cmd,
                stdout=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p)
        if p.returncode != 0:
            if self.fail_on_estimation_error:
                raise Exception(p.returncode)
            else:
                rows = ["NA" for i in range(len(trait_names))]
        else:
            rows = [row.strip() for row in stdout.split("\n")]
            rows = [float(row) for row in rows if row]
            assert len(rows) == len(trait_names), rows
        for field_name, rate in zip(trait_names, rows):
            profile_results["trait.{}.est.transition.rate".format(field_name)] = rate

    def estimate_pure_dispersal_rate(self,
            tree,
            profile_results,):
        self.create_bayestraits_geography_file(tree, output_path=self.geography_data_file_name)
        bt_commands = []
        bt_commands.append("1") # multstate
        bt_commands.append("1") # ml; 2 == mcmc
        bt_commands.append("restrictall q01")
        bt_commands.append("run")
        bt_commands = "\n".join(bt_commands)
        p = subprocess.Popen(
                [
                    "BayesTraits",
                    self.tree_file_name,
                    self.geography_data_file_name,
                ],
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p, bt_commands)
        stdout = stdout.split("\n")
        result = dict(zip(stdout[-3].split("\t"), stdout[-2].split("\t")))
        rate = float(result['q01'])
        profile_results["geographical.transition.rate"] = rate

    def estimate_dec_rates(self,
            tree,
            profile_results):
        self.create_biogeobears_geography_file(tree=tree)

    def estimate_lagrange_rates(self,
            tree,
            profile_results):
        if not self.estimate_lagrange_rates:
            return
        lagrange_commands = self.compose_lagrange_template(tree=tree)
        commandsf = open(self.commands_file_name, "w")
        commandsf.write(lagrange_commands)
        commandsf.flush()
        commandsf.close()
        shell_cmd = ["python", self.commands_file_name]
        p = subprocess.Popen(
                shell_cmd,
                stdout=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p)

    def create_working_tree_data(self, tree):
        with open(self.tree_file_name, "w") as tf:
            tree.write_to_stream(
                    tf,
                    "nexus",
                    suppress_leaf_taxon_labels=False,
                    suppress_leaf_node_labels=True,
                    suppress_internal_node_labels=True,
                    suppress_internal_taxon_labels=False,
                    translate_tree_taxa=True,
                    )
            tf.flush()
            tf.close()

    def create_traits_data(self,
            tree,
            generating_model=None):
        if generating_model is None:
            num_trait_types = len(tree.taxon_namespace[0].traits_vector)
            trait_names = ["trait{}".format(i+1) for i in range(num_trait_types)]
        else:
            num_trait_types = len(generating_model.geography.trait_types)
            trait_names = [trait.label for trait in generating_model.geography.trait_types]
        assert len(trait_names) == num_trait_types
        lineage_trait_states = collections.OrderedDict()
        for trait_idx in range(num_trait_types):
            state_symbols = {}
            for taxon in tree.taxon_namespace:
                lineage_label = taxon.label
                state = taxon.traits_vector[trait_idx]
                try:
                    symbol = state_symbols[state]
                except KeyError:
                    symbol = ArchipelagoProfiler.GEIGER_STATE_SYMBOLS[len(state_symbols)]
                    state_symbols[state] = symbol
                try:
                    lineage_trait_states[lineage_label].append(symbol)
                except KeyError:
                    lineage_trait_states[lineage_label] = [ symbol ]
        with open(self.traits_data_file_name, "w") as dataf:
            for lineage_label in lineage_trait_states:
                traits = ",".join(lineage_trait_states[lineage_label])
                dataf.write("{},{}\n".format(lineage_label, traits))
            dataf.flush()
            dataf.close()
        return trait_names

    def create_bayestraits_geography_file(self, tree, output_path):
        sep = "\t"
        dataf = open(output_path, "w")
        for taxon in tree.taxon_namespace:
            incidences = [str(i) for i in taxon.distribution_vector]
            dataf.write("{}{}{}\n".format(taxon.label, sep, sep.join(incidences)))
        dataf.flush()
        dataf.close()

    def create_biogeobears_geography_file(self, tree, output_path):
        sep = "\t"
        area_labels = ["a{}".format(idx+1) for idx, a in tree.taxon_namespace[0].distribution_vector]
        dataf = open(output_path, "w")
        dataf.write("{num_lineages}\t{num_areas}\t({area_labels})\n".format(
            len(tree.taxon_namespace),
            len(area_labels),
            " ".join(area_labels),
            ))
        for taxon in tree.taxon_namespace:
            incidences = [str(i) for i in taxon.distribution_vector]
            if area_names:
                assert len(area_names) == len(incidences)
            dataf.write("{}{}{}\n".format(taxon.label, sep, "".join(incidences)))
        dataf.flush()
        dataf.close()

    def compose_lagrange_template(self, tree):
        area_names = sorted(self.reconstruct_areas(tree))
        kwargs = {}
        kwargs["area_adjacency"] = str([[1] * len(area_names)] * len(area_names))
        kwargs["area_dispersal"] = str([[1.0] * len(area_names)] * len(area_names))
        kwargs["area_labels"] = str(area_names)
        kwargs["max_range_size"] = len(area_names)
        kwargs["model_name"] = str(id(tree))
        kwargs["newick_tree_string"] = tree.as_string("newick").replace("\n", "")
        assert tree.seed_node.age
        kwargs["root_age"] = tree.seed_node.age

        kwargs["taxon_name_list"] = [taxon.label for taxon in tree.taxon_namespace]

        taxon_range_data = {}
        for taxon in tree.taxon_namespace:
            taxon_area_indexes = tuple([area_idx for area_idx, i in enumerate(taxon.distribution_vector) if str(i) == "1"])
            taxon_range_data[taxon.label] = taxon_area_indexes
        kwargs["taxon_range_data"] = taxon_range_data

        ## EVERY PERMUTATION OF AREAS
        # ranges = []
        # for i in range(len(area_names)):
        #     x = list(itertools.permutations(area_indexes, i))
        #     ranges.extend(x)
        # ranges = sorted(set(ranges))
        # kwargs["ranges"] = str(ranges)

        ranges = set()
        for a in taxon_range_data.values():
            ranges.add(a)
        ranges = sorted(ranges)
        kwargs["ranges"] = str(ranges)

        return LAGRANGE_CONFIGURATION_TEMPLATE.format(**kwargs)

    def get_root_age(self, tree):
        tree.calc_node_ages()
        root_age = tree.seed_node.age
        return root_age

    def reconstruct_areas(self, tree):
        # cannot rely on generating model for number of
        # areas because we do not know if supplemental
        # areas are included in node data
        num_areas = len(tree.taxon_namespace[0].distribution_vector)
        area_names = ["a{}".format(i+1) for i in range(num_areas)]
        return area_names

    def preprocess_tree_taxa(self, tree):
        tree.original_taxon_namespace = tree.taxon_namespace
        tree.taxon_namespace = dendropy.TaxonNamespace()
        for node_idx, node in enumerate(tree.leaf_node_iter()):
            node.original_taxon = node.taxon
            node.taxon = tree.taxon_namespace.new_taxon(label=node.label)
            node.taxon.traits_vector = node.traits_vector
            node.taxon.distribution_vector = node.distribution_vector

    def preprocess_edge_lengths(self, tree):
        for nd in tree:
            nd.edge.original_length = nd.edge.length
            if nd.edge.length is None:
                nd.edge.length = 0.0
            elif nd.edge.length < 0.0:
                self.send_message("Setting negative edge length of {} to 0.0".format(nd.edge.length))
                nd.edge.length = 0.0

    def restore_edge_lengths(self, tree):
        for nd in tree:
            nd.edge.length = nd.edge.original_length
            del nd.edge.original_length

    def restore_tree_taxa(self, tree):
        tree.taxon_namespace = tree.original_taxon_namespace
        del tree.original_taxon_namespace
        for node_idx, node in enumerate(tree.leaf_node_iter()):
            node.taxon = node.original_taxon
            del node.original_taxon

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
    profiler = ArchipelagoProfiler()
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


