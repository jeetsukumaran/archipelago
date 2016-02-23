#! /usr/bin/env python
from __future__ import division
import sys
import os
import collections
import dendropy
import subprocess
from dendropy.calculate import treemeasure
from dendropy.model import birthdeath
from dendropy.utility import processio
from dendropy.calculate import statistics
from archipelago import model
from archipelago.utility import USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE

class Rcalculator(object):

    RESULT_FLAG_LEADER = "[!!!]"

    class HomogenousCommunityException(Exception):
        pass

    def __init__(self):
        pass

    def execute_rscript(self, script, prefix_key="predictor."):
        # x = open("t.R", "w")
        # x.write(script)
        # x.flush()
        # x.close()
        cmd = []
        cmd.append("Rscript")
        cmd.append("--vanilla")
        cmd.append("-")
        p = subprocess.Popen(cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p, script)
        if p.returncode != 0:
            print(script)
            for row in stderr.split("\n"):
                print("# {}".format(row))
            sys.exit(p.returncode)
        results = {}
        num_lines_with_results = 0
        for line in stdout.split("\n"):
            if not line.startswith(Rcalculator.RESULT_FLAG_LEADER):
                continue
            parts = line[len(Rcalculator.RESULT_FLAG_LEADER):].split("=")
            assert len(parts) == 2
            key = parts[0].strip()
            try:
                value = float(parts[1].strip())
            except ValueError as e:
                value = "NA"
            results[prefix_key + key] = value
            num_lines_with_results += 1
        return results

    def _compose_cophenetic_matrix(
            self,
            dists,
            taxon_names,
            byrow=True,
            ):
        return "matrix(c({data}), nrow={nrow}, byrow={byrow}, dimnames=list(c({names}),c({names})))".format(
            data=",".join("{}".format(d) for d in dists),
            nrow=len(taxon_names),
            byrow="T" if byrow else "F",
            names=",".join(taxon_names))

    def _compose_community_matrix(
            self,
            data,
            comm_names,
            taxon_names,
            byrow=True,
            ):
        if len(set(data)) == 1:
            raise Rcalculator.HomogenousCommunityException()
        return "matrix(c({data}), nrow={nrow}, byrow={byrow}, dimnames=list(c({comm_names}),c({taxon_names})))".format(
            data=",".join("{}".format(d) for d in data),
            nrow=len(comm_names),
            byrow="T" if byrow else "F",
            comm_names=",".join(comm_names),
            taxon_names=",".join(taxon_names)
            )

    def calc_summary_stats(
            self,
            tree,
            patristic_distance_matrix,
            total_tree_length,
            total_tree_edges,
            area_taxa,
            trait_taxa,
            ):
        pdm = patristic_distance_matrix
        taxon_names = []
        weighted_dists = []
        unweighted_dists = []
        normalized_weighted_dists = []
        normalized_unweighted_dists = []
        for taxon1 in tree.taxon_namespace:
            taxon_names.append("'{}'".format(taxon1.label))
            # for taxon2 in tree.taxon_namespace:
            #     weighted_dist = pdm(taxon1, taxon2)
            #     unweighted_dist = pdm.path_edge_count(taxon1, taxon2)
            #     normalized_weighted_dist = weighted_dist / total_tree_length
            #     normalized_unweighted_dist = unweighted_dist / total_tree_edges
            #     weighted_dists.append(weighted_dist)
            #     unweighted_dists.append(unweighted_dist)
            #     normalized_weighted_dists.append(normalized_weighted_dist)
            #     normalized_unweighted_dists.append(normalized_unweighted_dist)
        rscript = []
        rscript.append("suppressMessages(library(picante))")
        community_regimes = []
        community_regimes.append(
            ("Area", area_taxa, "by_area"),
        )
        # community_regimes.append(
        #     ("Trait", trait_taxa, "by_area"),
        # )
        for trait_idx in trait_taxa:
            # sys.stderr.write("{}: {}\n".format(individual_trait_taxa, trait_taxa[individual_trait_taxa]))
            community_regimes.append(
                    ("Trait.{}.".format(trait_idx + USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE), trait_taxa[trait_idx], "by_trait_{}".format(trait_idx + USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE))
            )
            # for trait_state in trait_taxa[trait_idx]:
            #     community_regimes.append(
            #             ("Trait.{}.{}".format(trait_idx, trait_state), trait_taxa[trait_idx][trait_state], "by_trait_{}_{}".format(trait_idx, trait_state))
            #     )

        for dists, dists_desc in (
                    # (weighted_dists, "weighted"),
                    # (unweighted_dists, "unweighted"),
                    (normalized_weighted_dists, "normalized_weighted"),
                    (normalized_unweighted_dists, "normalized_unweighted"),
                ):
            cophenetic_dist_matrix_str = self._compose_cophenetic_matrix(
                    dists=dists,
                    taxon_names=taxon_names)
            cophenetic_dist_matrix_name = "{}_cophenetic_dist_matrix".format(dists_desc)
            rscript.append("{} <- {}".format(cophenetic_dist_matrix_name, cophenetic_dist_matrix_str))
            for comm_prefix, comm_data, comm_desc in community_regimes:
                comm_names = []
                pa_data = []
                for idx in comm_data:
                    comm_taxa = [taxon for taxon in comm_data[idx]]
                    # print("# {}: {}: {} of {}: {}\n".format(comm_desc, idx, len(comm_taxa), len(tree_taxa), ", ".join(t.label for t in comm_taxa)))
                    # sys.stderr.write("# {}: {}: {} of {}: {}\n".format(comm_desc, idx, len(comm_taxa), len(tree.taxon_namespace), ", ".join(t.label for t in comm_taxa)))
                    comm_names.append("'{}{}'".format(comm_prefix, idx))
                    # sys.stderr.write("comm_names={}\n".format(comm_names))
                    for taxon in tree.taxon_namespace:
                        if taxon in comm_taxa:
                            pa_data.append(1)
                        else:
                            pa_data.append(0)
                comm_pa_matrix_name = "community_{}".format(comm_desc)
                try:
                    comm_pa_matrix_str = self._compose_community_matrix(
                            data=pa_data,
                            comm_names=comm_names,
                            taxon_names=["'{}'".format(t.label) for t in tree.taxon_namespace])
                except Rcalculator.HomogenousCommunityException:
                    continue
                # sys.stderr.write("{} <- {}\n".format(comm_pa_matrix_name, comm_pa_matrix_str))
                rscript.append("{} <- {}".format(comm_pa_matrix_name, comm_pa_matrix_str))
                nruns = 100
                prefix = "{comm}.{dists}".format(comm=comm_pa_matrix_name,
                        dists=cophenetic_dist_matrix_name.replace("_cophenetic_dist_matrix", "")).replace("_", ".")
                out = "stdout()"
                # out = "'z.txt'"
                for stat_type in ("mpd", "mntd"):
                    rscript.append("result <- ses.{stat_type}({comm},{dists},null.model='taxa.labels',abundance.weighted=FALSE,runs={nruns})".format(
                        stat_type=stat_type,
                        comm=comm_pa_matrix_name,
                        dists=cophenetic_dist_matrix_name,
                        nruns=nruns))
                    rscript.append("result.df <- as.data.frame(result)")
                    if comm_desc.startswith("by_trait"):
                        rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.Z.', rownames(result.df), ' = ', result.df${stat_type}.obs.z, '\\n', sep=''), {out})".format(
                            stat_type=stat_type,
                            result_flag=Rcalculator.RESULT_FLAG_LEADER,
                            prefix=prefix,
                            out=out,
                            ))
                        rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.p.', rownames(result.df), ' = ', result.df${stat_type}.obs.p, '\\n', sep=''), {out})".format(
                            stat_type=stat_type,
                            result_flag=Rcalculator.RESULT_FLAG_LEADER,
                            prefix=prefix,
                            out=out,
                            ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.Z.mean', ' = ', mean(result.df${stat_type}.obs.z), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.p.mean', ' = ', mean(result.df${stat_type}.obs.p), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.Z.var', ' = ', var(result.df${stat_type}.obs.z), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.p.var', ' = ', var(result.df${stat_type}.obs.p), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
        rscript = "\n".join(rscript)
        results = self.execute_rscript(rscript)
        tree.stats.update(results)
        return results

class TreeSummarizer(object):

    class IncompleteRaditionException(Exception):
        pass

    class IncompleteAreaRadiationException(IncompleteRaditionException):
        pass

    class IncompleteTraitRaditionException(IncompleteRaditionException):
        pass

    def __init__(self,
            drop_trees_not_spanning_all_areas=True,
            drop_trees_not_spanning_multiple_traits=False,
            trait_indexes_to_exclude=None,
            trait_states_to_exclude=None,
            ):
        """
        Creates summaries for trees.

        Parameters
        ----------
        drop_trees_not_spanning_all_areas : bool
            Skip calculations for trees that do not span all areas.
        drop_trees_not_spanning_multiple_traits : bool
            Skip calculations for trees that do not span more than one trait.
        trait_indexes_to_exclude : iterable
            0-based indexes of traits to skip in calculations.
        trait_states_to_exclude : iterable of tuples
            Tuples in the form of (a,b), where 'a' is the 0-based index of the
            trait and 'b' is the state to skip in calculations.
        """
        self.drop_trees_not_spanning_all_areas = drop_trees_not_spanning_all_areas
        self.drop_trees_not_spanning_multiple_traits = drop_trees_not_spanning_multiple_traits
        if trait_indexes_to_exclude:
            self.trait_indexes_to_exclude = set(trait_indexes_to_exclude)
        else:
            self.trait_indexes_to_exclude = set([])
        if trait_states_to_exclude:
            self.trait_states_to_exclude = set(trait_states_to_exclude)
        else:
            self.trait_states_to_exclude = set([])
        self.rcalc = Rcalculator()
        self.stat_name_prefix = "predictor"
        self.stat_name_delimiter = "."
        self.num_randomization_replicates = 10

    def get_mean_patristic_distance(self, pdm, nodes):
        if len(nodes) <= 1:
            return "NA", "NA"
        weighted_dist = 0.0
        unweighted_dist = 0.0
        ncomps = 0
        for idx1, nd1 in enumerate(nodes[:-1]):
            for nd2 in nodes[idx1+1:]:
                assert nd1.taxon
                assert nd2.taxon
                weighted_dist += pdm(nd1.taxon, nd2.taxon)
                unweighted_dist += pdm.path_edge_count(nd1.taxon, nd2.taxon)
                ncomps += 1
        return weighted_dist/ncomps, unweighted_dist/ncomps

    def summarize_trees(self, trees):
        processed_trees = []
        summary_fieldnames = set()
        summary_results = []
        for tree in list(trees):
            try:
                self.summarize_tree(tree)
                processed_trees.append(tree)
                summary_fieldnames.update(tree.stats.keys())
                summary_results.append(collections.OrderedDict(tree.stats))
            except TreeSummarizer.IncompleteAreaRadiationException:
                pass
        return processed_trees, summary_fieldnames, summary_results

    def summarize_tree(self, tree):
        self.preprocess_tree_taxa(tree)
        area_taxa_map = collections.defaultdict(list)
        trait_taxa_map = collections.defaultdict(lambda: collections.defaultdict(list))
        for taxon in tree.taxon_namespace:
            for area_idx, area_presence in enumerate(taxon.distribution_vector):
                if area_presence:
                    area_taxa_map[area_idx].append(taxon)
            for trait_idx, trait_state in enumerate(taxon.traits_vector):
                # sys.stderr.write("==> {}   {}   {}\n".format(trait_idx, trait_state,taxon))
                if self.trait_indexes_to_exclude and trait_idx in self.trait_indexes_to_exclude:
                    continue
                if self.trait_states_to_exclude and (trait_idx, trait_state) in self.trait_states_to_exclude:
                    continue
                trait_taxa_map[trait_idx][trait_state].append(taxon)
        num_areas = len(tree.taxon_namespace[0].distribution_vector)
        if len(area_taxa_map) < num_areas and self.drop_trees_not_spanning_all_areas:
            raise TreeSummarizer.IncompleteAreaRadiationException()
        for trait_idx in trait_taxa_map:
            if len(trait_taxa_map[trait_idx]) < 2 and self.drop_trees_not_spanning_multiple_traits:
                raise TreeSummarizer.IncompleteTraitRaditionException()

        # print("---")
        # for a in area_taxa:
        #     print("{}: {}".format(a, [x.label for x in area_taxa[a]]))
        # print("---")
        # for t in trait_taxa[0]:
        #     print("{}: {}".format(t, [x.label for x in trait_taxa[0][t]]))

        tree.stats = collections.defaultdict(lambda:"NA")
        tree.stats["model_id"] = tree.annotations.get_value("model_id", "NA")
        total_tree_length = 0.0
        total_tree_edges = 0.0
        for nd in tree:
            total_tree_edges += 1.0
            if nd.edge.length is None:
                nd.edge.length = 0.0
            total_tree_length += nd.edge.length
        # rstats = self.rcalc.calc_summary_stats(
        #         tree=tree,
        #         patristic_distance_matrix=pdm,
        #         total_tree_length=total_tree_length,
        #         total_tree_edges=total_tree_edges,
        #         area_taxa=area_taxa,
        #         trait_taxa=trait_taxa,
        #         )
        self.calc_summary_stats(
                tree=tree,
                area_taxa_map=area_taxa_map,
                trait_taxa_map=trait_taxa_map)
        self.restore_tree_taxa(tree)

    def preprocess_tree_taxa(self, tree):
        model.ArchipelagoModel.set_lineage_data(
                tree=tree,
                leaf_nodes_only=True,
                lineage_data_source='node')
        tree.original_taxon_namespace = tree.taxon_namespace
        tree.taxon_namespace = dendropy.TaxonNamespace()
        for node_idx, node in enumerate(tree.leaf_node_iter()):
            node.original_taxon = node.taxon
            node.taxon = tree.taxon_namespace.new_taxon(label=node.label)
            node.taxon.traits_vector = node.traits_vector
            node.taxon.distribution_vector = node.distribution_vector

    def restore_tree_taxa(self, tree):
        tree.taxon_namespace = tree.original_taxon_namespace
        del tree.original_taxon_namespace
        for node_idx, node in enumerate(tree.leaf_node_iter()):
            node.taxon = node.original_taxon
            del node.original_taxon

    def calc_summary_stats(self,
            tree,
            area_taxa_map,
            trait_taxa_map,):
        stats = {}
        phylogenetic_distance_matrix = tree.phylogenetic_distance_matrix()
        stats.update(self._calc_trait_based_stats(
            phylogenetic_distance_matrix=phylogenetic_distance_matrix,
            trait_taxa_map=trait_taxa_map,
            ))
        tree.stats.update(stats)
        return stats

    def _calc_trait_based_stats(self,
            phylogenetic_distance_matrix,
            trait_taxa_map,
            ):
        trait_assemblage_memberships, trait_assemblage_descriptions = self._get_trait_community_regimes(trait_taxa_map)
        assert len(trait_assemblage_memberships) == len(trait_assemblage_descriptions)
        results = self._calc_community_ecology_stats(
                phylogenetic_distance_matrix=phylogenetic_distance_matrix,
                assemblage_memberships=trait_assemblage_memberships,
                assemblage_descriptions=trait_assemblage_descriptions,
                report_character_state_specific_results=True,
                report_character_class_wide_results=True,
                )
        return results
        # for assemblage_desc, results_desc in zip(results, trait_assemblage_descriptions):

    def _get_trait_community_regimes(self, trait_taxa_map):
        assemblage_descriptions = []
        assemblage_memberships = []
        for trait_idx in trait_taxa_map:
            for trait_state_idx in trait_taxa_map[trait_idx]:
                assemblage_memberships.append( trait_taxa_map[trait_idx][trait_state_idx] )
                regime = {
                    "assemblage_basis_class_id": "trait{}{}".format(self.stat_name_delimiter, trait_idx + USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE),
                    "assemblage_basis_state_id": "state{}".format(trait_state_idx),
                    "assemblage_name": "Trait{}".format(
                        trait_idx + USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE,
                        )
                }
                assemblage_descriptions.append(regime)
        return assemblage_memberships, assemblage_descriptions

    def _calc_community_ecology_stats(self,
            phylogenetic_distance_matrix,
            assemblage_memberships,
            assemblage_descriptions,
            report_character_state_specific_results=True,
            report_character_class_wide_results=True,
            ):
        summary_statistics_suite = {}
        assert len(assemblage_descriptions) == len(assemblage_memberships)
        for edge_weighted_desc in ("weighted", "unweighted"):
            if edge_weighted_desc:
                is_weighted_edge_distances = True
            else:
                is_weighted_edge_distances = False
            for underlying_statistic_type_desc in ("mpd", "mntd"):
                if underlying_statistic_type_desc == "mpd":
                    stat_fn_name = "standardized_effect_size_mean_pairwise_distance"
                else:
                    stat_fn_name = "standardized_effect_size_mean_nearest_taxon_distance"
                stat_fn = getattr(phylogenetic_distance_matrix, stat_fn_name)
                results_group = stat_fn(
                    assemblage_memberships=assemblage_memberships,
                    is_weighted_edge_distances=is_weighted_edge_distances,
                    is_normalize_by_tree_size=True,
                    num_randomization_replicates=self.num_randomization_replicates,
                    )
                assert len(results_group) == len(assemblage_memberships)
                results_by_character_class = {}
                results_by_character_class["z"] = collections.defaultdict(list)
                results_by_character_class["p"] = collections.defaultdict(list)
                for result, assemblage_desc in zip(results_group, assemblage_descriptions):
                    character_class_statistic_prefix = self.stat_name_delimiter.join([
                        self.stat_name_prefix,
                        "community",
                        "by",
                        assemblage_desc["assemblage_basis_class_id"],
                        edge_weighted_desc,
                        underlying_statistic_type_desc,
                        # assemblage_desc["assemblage_basis_state_id"],
                        ])
                    for ses_result_statistic in ("z", "p"): # z = score, p = p-value; turns out the latter is quite informative
                        ses_result_statistic_value = getattr(result, ses_result_statistic)
                        if report_character_state_specific_results:
                            character_state_statistic_name = self.stat_name_delimiter.join([
                                character_class_statistic_prefix,
                                ses_result_statistic,
                                assemblage_desc["assemblage_basis_state_id"],
                                ])
                            assert character_state_statistic_name not in summary_statistics_suite
                            summary_statistics_suite[character_state_statistic_name] = ses_result_statistic_value
                        if report_character_class_wide_results:
                            results_by_character_class[ses_result_statistic][character_class_statistic_prefix].append(ses_result_statistic_value)
        if report_character_class_wide_results:
            for ses_result_statistic in results_by_character_class:
                assert len(results_by_character_class[ses_result_statistic]) > 0
                for character_class_statistic_prefix in results_by_character_class[ses_result_statistic]:
                    sn_title = self.stat_name_delimiter.join([
                        character_class_statistic_prefix,
                        ses_result_statistic,
                        ])
                    mean, var = statistics.mean_and_sample_variance(results_by_character_class[ses_result_statistic][character_class_statistic_prefix])
                    key1 = "{}{}{}".format(sn_title, self.stat_name_delimiter, "mean")
                    assert key1 not in summary_statistics_suite
                    summary_statistics_suite[key1] = mean
                    key2 = "{}{}{}".format(sn_title, self.stat_name_delimiter, "var")
                    assert key2 not in summary_statistics_suite
                    summary_statistics_suite[key2] = var
        return summary_statistics_suite
