#! /usr/bin/env python

import sys
import os
import csv
import collections
import tempfile
import subprocess
import re
import argparse
import dendropy
from dendropy.model import birthdeath
from dendropy.utility import processio
from archipelago import model
from archipelago.estimate_biogeobears import BiogeobearsEstimator

class ArchipelagoProfiler(object):

    GEIGER_STATE_SYMBOLS = "123456789"
    LAGRANGE_CPP_EXTRACT_PATTERN = re.compile(r".*^dis: ([0-9eE\-\.]+) ext: ([0-9eE\-\.]+).*", re.MULTILINE | re.DOTALL)
    DEFAULT_MINIMUM_BRANCH_LENGTH = 1e-6

    @staticmethod
    def get_profile_options_parser():
        parser = argparse.ArgumentParser(add_help=False)
        profile_metric_options = parser.add_argument_group("Profile Metrics")
        profile_metric_options.add_argument("--no-estimate-pure-birth",
                action="store_true",
                default=False,
                help="Do NOT estimate birth rate under a pure-birth model.")
        profile_metric_options.add_argument("--no-estimate-trait-transition",
                action="store_true",
                default=False,
                help="Do NOT estimate trait transition rate.")
        profile_metric_options.add_argument("--no-estimate-area-transition",
                action="store_true",
                default=False,
                help="Do NOT estimate area transition rate.")
        profile_metric_options.add_argument("--trait-transition-rate-estimation-method",
                default="bayestraits",
                help="Estimate trait transition rates using 'geiger' or 'bayestraits'")
        profile_metric_options.add_argument("--estimate-dec-biogeobears",
                action="store_true",
                default=False,
                help="Estimate parameters under Lagrange's DEC model (using BioGeoBears).")
        profile_metric_options.add_argument("--estimate-dec-lagrange",
                action="store_true",
                default=False,
                help="Estimate parameters under Lagrange's DEC model (using Lagrange).")
        profile_metric_options.add_argument("-b", "--minimum-branch-length",
                default=ArchipelagoProfiler.DEFAULT_MINIMUM_BRANCH_LENGTH,
                type=float,
                help="Minimum branch length (edges will be forced to this length).")

        profile_run_options = parser.add_argument_group("Profile Run Options")
        profile_run_options.add_argument("-v", "--verbosity",
                default=1,
                type=int,
                help="Progress noisiness")
        profile_run_options.add_argument("--debug-mode",
                action="store_true",
                default=False,
                help="Run in debug mode (work files will not be deleted).")
        profile_run_options.add_argument("--ignore-estimation-errors",
                action="store_true",
                default=False,
                help="Ignore errors raised by estimation internally or by external programs")
        return parser

    @staticmethod
    def from_option_args(args):
        profiler = ArchipelagoProfiler(
                minimum_branch_length = args.minimum_branch_length,
                is_estimate_pure_birth_rate=not args.no_estimate_pure_birth,
                is_estimate_trait_transition_rates=not args.no_estimate_trait_transition,
                trait_transition_rate_estimation_method=args.trait_transition_rate_estimation_method,
                is_estimate_area_transition_rates=not args.no_estimate_area_transition,
                is_estimate_dec_biogeobears=args.estimate_dec_biogeobears,
                is_estimate_dec_lagrange=args.estimate_dec_lagrange,
                verbosity=args.verbosity,
                fail_on_estimation_error=not args.ignore_estimation_errors,
                debug_mode=args.debug_mode,
                )
        return profiler

    def __init__(self,
            is_estimate_pure_birth_rate=True,
            is_estimate_trait_transition_rates=True,
            is_estimate_area_transition_rates=True,
            trait_transition_rate_estimation_method="bayestraits",
            is_estimate_dec_biogeobears=False,
            is_estimate_dec_lagrange=False,
            minimum_branch_length=DEFAULT_MINIMUM_BRANCH_LENGTH,
            verbosity=2,
            fail_on_estimation_error=True,
            debug_mode=False):
        self.minimum_branch_length = minimum_branch_length
        self.is_estimate_pure_birth_rate = is_estimate_pure_birth_rate
        self.is_estimate_trait_transition_rates = is_estimate_trait_transition_rates
        self.is_estimate_area_transition_rates = is_estimate_area_transition_rates
        self.trait_transition_rate_estimation_method = trait_transition_rate_estimation_method
        self.is_estimate_dec_biogeobears = is_estimate_dec_biogeobears
        self.is_estimate_dec_lagrange = is_estimate_dec_lagrange
        self.quiet = verbosity <= 1
        self.fail_on_estimation_error = fail_on_estimation_error
        self.debug_mode = debug_mode
        if self.debug_mode:
            self.tree_file_name = "profiler.tree.nexus"
            self.newick_tree_file_name = "profiler.tree.newick"
            self.traits_data_file_name = "profiler.traits.data.txt"
            self.geography_data_file_name = "profiler.geography.data.txt"
            self.commands_file_name = "profiler.commands.txt"
            self.output_file_name = "profiler.output.txt"
        else:
            self.tree_file = tempfile.NamedTemporaryFile()
            self.newick_tree_file = tempfile.NamedTemporaryFile()
            self.traits_data_file = tempfile.NamedTemporaryFile()
            self.geography_data_file = tempfile.NamedTemporaryFile()
            self.commands_file = tempfile.NamedTemporaryFile()
            self.output_file = tempfile.NamedTemporaryFile()
            self.tree_file_name = self.tree_file.name
            self.newick_tree_file_name = self.newick_tree_file.name
            self.traits_data_file_name = self.traits_data_file.name
            self.geography_data_file_name = self.geography_data_file.name
            self.commands_file_name = self.commands_file.name
            self.output_file_name = self.output_file.name
        self.biogeobears_estimator = BiogeobearsEstimator(
                commands_file_name=self.commands_file_name,
                results_file_name=self.output_file_name,
                debug_mode=self.debug_mode,
                )

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
        writer = csv.DictWriter(dest,
                fieldnames=fieldnames,
                lineterminator=os.linesep,
                )
        if not suppress_headers:
            writer.writeheader()
        writer.writerows(profiles)

    def profile_trees_from_path(self,
            trees_filepath,
            schema="newick",
            generating_model=None,
            is_lineages_decoded=False,
            lineage_data_source="node",
            traits_filepath=None,
            distribution_filepath=None,
            ):
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
                is_lineages_decoded=is_lineages_decoded,
                lineage_data_source=lineage_data_source,
                traits_filepath=traits_filepath,
                distribution_filepath=distribution_filepath,
                )
        return profiles

    def profile_trees(self,
            trees,
            generating_model=None,
            is_lineages_decoded=False,
            lineage_data_source="node",
            traits_filepath=None,
            distribution_filepath=None,
            ):
        profiles = []
        for tree_idx, tree in enumerate(trees):
            if hasattr(trees, "tree_filepath"):
                tree.tree_filepath = trees.tree_filepath
                tree.tree_offset = tree_idx
            r = self.profile_tree(
                    tree=tree,
                    generating_model=generating_model,
                    is_lineages_decoded=is_lineages_decoded,
                    lineage_data_source=lineage_data_source,
                    traits_filepath=traits_filepath,
                    distribution_filepath=distribution_filepath,
                    )
            profiles.append(r)
        return profiles

    def profile_tree(self,
            tree,
            generating_model=None,
            is_lineages_decoded=False,
            lineage_data_source="node",
            traits_filepath=None,
            distribution_filepath=None,
            ):
        if not is_lineages_decoded:
            model.ArchipelagoModel.set_lineage_data(
                    tree=tree,
                    leaf_nodes_only=True,
                    lineage_data_source=lineage_data_source,
                    traits_filepath=traits_filepath,
                    distribution_filepath=distribution_filepath,
                    )

        # intitialize profile_results
        profile_results = collections.OrderedDict()

        # basic profile
        if hasattr(tree, "tree_filepath"):
            profile_results["tree.filepath"] = tree.tree_filepath
            if hasattr(tree, "tree_offset"):
                profile_results["tree.offset"] = tree.tree_offset
        profile_results["num.tips"] = len(list(nd for nd in tree.leaf_node_iter()))
        tree.calc_node_ages()
        root_age = tree.seed_node.age
        profile_results["root.age"] = root_age

        # prepopulate with params that are available
        if generating_model is not None:
            self.store_generating_model_parameters(generating_model, profile_results)

        # estimate birth rate
        if self.is_estimate_pure_birth_rate:
            birth_rate = self.estimate_pure_birth_rate(tree, profile_results)

        # set up for external processing
        ## fix zero or negative edge lengths
        self.preprocess_edge_lengths(tree)
        ## set up taxa
        self.preprocess_tree_taxa(tree)
        ## store tree
        self.create_working_tree_data(tree)

        ## process traits
        self.estimate_trait_transition_rates(
                tree=tree,
                profile_results=profile_results,
                generating_model=generating_model)

        # process areas
        if len(tree.taxon_namespace[0].distribution_vector) > 1:
            if self.is_estimate_area_transition_rates:
                self.estimate_pure_dispersal_weight(
                        tree=tree,
                        profile_results=profile_results)
            if self.is_estimate_dec_biogeobears:
                self.estimate_dec_rates_biogeobears(
                        tree=tree,
                        profile_results=profile_results,)
            if self.is_estimate_dec_lagrange:
                self.estimate_dec_rates_lagrange(
                        tree=tree,
                        profile_results=profile_results,)

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
        profile_results["lineage.dispersal.weight.definition"] = generating_model.lineage_dispersal_weight_function.definition_content
        profile_results["lineage.dispersal.weight.description"] = generating_model.lineage_dispersal_weight_function.description
        profile_results["global.dispersal.rate"] = generating_model.global_dispersal_rate
        profile_results["mean.dispersal.rate"] = generating_model.mean_dispersal_rate
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
        return rate

    def estimate_trait_transition_rates(self,
            tree,
            profile_results,
            generating_model):

        if not self.is_estimate_trait_transition_rates:
            return None

        if generating_model is None:
            num_trait_types = len(tree.taxon_namespace[0].traits_vector)
            trait_names = ["trait{}".format(i+1) for i in range(num_trait_types)]
        else:
            num_trait_types = len(generating_model.trait_types)
            trait_names = [trait.label for trait in generating_model.trait_types]
        if num_trait_types == 0:
            return None

        if self.trait_transition_rate_estimation_method == "geiger":
            return self.estimate_trait_transition_rates_using_geiger(
                    tree=tree,
                    profile_results=profile_results,
                    trait_names=trait_names)
        elif self.trait_transition_rate_estimation_method == "bayestraits":
            return self.estimate_trait_transition_rates_using_bayestraits(
                    tree=tree,
                    profile_results=profile_results,
                    trait_names=trait_names)
        else:
            raise ValueError("Expecting 'geiger' or 'bayestraits', but instead found: '{}'".format(self.trait_transition_rate_estimation_method))

    def estimate_trait_transition_rates_using_bayestraits(self,
            tree,
            profile_results,
            trait_names):
        for trait_idx, trait_name in enumerate(trait_names):
            symbols = self.create_bayestraits_traits_data(
                    tree,
                    trait_idx,
                    output_path=self.traits_data_file_name)
            master_rate = "q{}{}".format(symbols[0],symbols[1])
            bt_commands = []
            bt_commands.append("1") # multstate
            bt_commands.append("1") # ml; 2 == mcmc
            bt_commands.append("restrictall {}".format(master_rate))
            bt_commands.append("run")
            bt_commands = "\n".join(bt_commands)
            p = subprocess.Popen(
                    [
                        "BayesTraits",
                        self.tree_file_name,
                        self.traits_data_file_name,
                    ],
                    stdout=subprocess.PIPE,
                    stdin=subprocess.PIPE,
                    )
            stdout, stderr = processio.communicate(p, bt_commands)
            stdout_rows = stdout.split("\n")
            targeted_row_idx = None
            for row_idx, row in enumerate(stdout_rows):
                # if "q01" in row and "q10" in row:
                if row.startswith("Tree No\tLh\tq"):
                    targeted_row_idx = row_idx + 1
                    break
            if targeted_row_idx is None:
                if self.fail_on_estimation_error:
                    raise Exception("Failed to extract results from BayesTraits estimation")
                else:
                    rate = "NA"
            else:
                result = dict(zip(stdout_rows[targeted_row_idx-1].split("\t"), stdout_rows[targeted_row_idx].split("\t")))
                rate = float(result[master_rate])
            profile_results["trait.{}.est.transition.rate".format(trait_name)] = rate

    def estimate_trait_transition_rates_using_geiger(self,
            tree,
            profile_results,
            trait_names):
        self.create_geiger_traits_data(tree=tree, num_trait_types=len(trait_names))
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

    def estimate_pure_dispersal_weight(self,
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
        stdout_rows = stdout.split("\n")
        targeted_row_idx = None
        for row_idx, row in enumerate(stdout_rows):
            # if "q01" in row and "q10" in row:
            if row.startswith("Tree No\tLh\tq"):
                targeted_row_idx = row_idx + 1
                break
        if targeted_row_idx is None:
            if self.fail_on_estimation_error:
                raise Exception("Failed to extract results from BayesTraits estimation")
            else:
                rate = "NA"
        else:
            result = dict(zip(stdout_rows[targeted_row_idx-1].split("\t"), stdout_rows[targeted_row_idx].split("\t")))
            rate = float(result['q01'])
        profile_results["area.est.transition.rate"] = rate

    def estimate_dec_rates_biogeobears(self,
            tree,
            profile_results,
            **kwargs):
        tree.write_to_path(self.newick_tree_file_name, "newick")
        self.create_biogeobears_geography_file(tree=tree, output_path=self.geography_data_file_name)
        dec_results = self.biogeobears_estimator.estimate_dec(
                newick_tree_filepath=self.newick_tree_file_name,
                geography_filepath=self.geography_data_file_name,
                max_range_size=len(tree.taxon_namespace[0].distribution_vector),
                **kwargs
                )
        profile_results["biogeobears.dec.dispersal.rate"] = dec_results["d"]
        profile_results["biogeobears.dec.extinction.rate"] = dec_results["e"]

    def estimate_dec_rates_lagrange(self,
            tree,
            profile_results,
            **kwargs):
        tree.write_to_path(
                self.newick_tree_file_name,
                "newick",
                suppress_rooting=True,
                )
        self.create_lagrangecpp_geography_file(tree=tree, output_path=self.geography_data_file_name)
        configf = open(self.commands_file_name, "w")
        configf.write("treefile = {}\n".format(self.newick_tree_file_name))
        configf.write("datafile = {}\n".format(self.geography_data_file_name))
        configf.flush()
        configf.close()
        shell_cmd = ["lagrange_cpp", self.commands_file_name]
        try:
            p = subprocess.Popen(
                    shell_cmd,
                    stdout=subprocess.PIPE,
                    )
        except OSError as e:
            raise OSError("Failed to execute command: {}".format(" ".join(shell_cmd)))
        stdout, stderr = processio.communicate(p)
        if p.returncode != 0:
            if self.fail_on_estimation_error:
                raise Exception(p.returncode)
            else:
                profile_results["lagrange.dec.dispersal.rate"] = "NA"
                profile_results["lagrange.dec.extinction.rate"] = "NA"
        else:
            match = ArchipelagoProfiler.LAGRANGE_CPP_EXTRACT_PATTERN.match(stdout)
            if not match:
                if self.fail_on_estimation_error:
                    raise Exception("Failed to extract results from Lagrange estimation")
                else:
                    profile_results["lagrange.dec.dispersal.rate"] = "NA"
                    profile_results["lagrange.dec.extinction.rate"] = "NA"
            else:
                results = match.groups(1)
                profile_results["lagrange.dec.dispersal.rate"] = float(results[0])
                profile_results["lagrange.dec.extinction.rate"] = float(results[1])

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

    def create_geiger_traits_data(self, tree, num_trait_types):
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

    def create_bayestraits_traits_data(self,
            tree,
            trait_idx,
            output_path):
        sep = "\t"
        dataf = open(output_path, "w")
        symbols = set()
        for taxon in tree.taxon_namespace:
            trait_state = str(taxon.traits_vector[trait_idx])
            symbols.add(trait_state)
            dataf.write("{}{}{}\n".format(taxon.label, sep, trait_state))
        dataf.flush()
        dataf.close()
        return sorted(symbols)

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
        area_labels = ["a{}".format(idx+1) for idx, a in enumerate(tree.taxon_namespace[0].distribution_vector)]
        dataf = open(output_path, "w")
        dataf.write("{num_lineages}\t{num_areas}\t({area_labels})\n".format(
            num_lineages=len(tree.taxon_namespace),
            num_areas=len(area_labels),
            area_labels=" ".join(area_labels),
            ))
        for taxon in tree.taxon_namespace:
            incidences = [str(i) for i in taxon.distribution_vector]
            dataf.write("{}{}{}\n".format(taxon.label, sep, "".join(incidences)))
        dataf.flush()
        dataf.close()

    def create_lagrangecpp_geography_file(self, tree, output_path):
        sep = "\t"
        dataf = open(output_path, "w")
        dataf.write("{num_lineages}\t{num_areas})\n".format(
            num_lineages=len(tree.taxon_namespace),
            num_areas=len(tree.taxon_namespace[0].distribution_vector),
            ))
        for taxon in tree.taxon_namespace:
            incidences = [str(i) for i in taxon.distribution_vector]
            dataf.write("{}{}{}\n".format(taxon.label, sep, "".join(incidences)))
        dataf.flush()
        dataf.close()

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
                nd.edge.length = self.minimum_branch_length
            elif nd.edge.length <= self.minimum_branch_length:
                self.send_message("Setting edge length of {} to {}".format(nd.edge.length, self.minimum_branch_length))
                nd.edge.length = self.minimum_branch_length

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




