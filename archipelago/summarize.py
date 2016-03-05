#! /usr/bin/env python
from __future__ import division
import sys
import os
import collections
import dendropy
from dendropy.calculate import treemeasure
from dendropy.model import birthdeath
from dendropy.utility import processio
from dendropy.calculate import statistics
from archipelago import model
from archipelago import utility
from archipelago.utility import USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE

class TreeSummarizer(object):

    class SingleTaxonAssemblageException(Exception):
        pass

    class SingleTaxonAreaAssemblageException(SingleTaxonAssemblageException):
        def __init__(self, area_idx, current_tree_idx):
            TreeSummarizer.SingleTaxonAssemblageException.__init__(self)
            self.area_idx = area_idx

    class SingleTaxonTraitStateAssemblageException(SingleTaxonAssemblageException):
        def __init__(self, trait_idx, trait_state_idx, current_tree_idx):
            TreeSummarizer.SingleTaxonAssemblageException.__init__(self)
            self.trait_idx = trait_idx
            self.trait_state_idx = trait_state_idx

    class IncompleteRaditionException(Exception):
        pass

    class IncompleteAreaRadiationException(IncompleteRaditionException):
        pass

    class IncompleteTraitRaditionException(IncompleteRaditionException):
        pass

    def __init__(self,
            drop_trees_not_spanning_all_areas=True,
            drop_trees_not_spanning_multiple_traits=False,
            drop_trees_with_single_lineage_areas=False,
            drop_trees_with_single_lineage_trait_states=False,
            trait_indexes_to_exclude=None,
            trait_states_to_exclude=None,
            run_logger=None,
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
        drop_trees_with_single_lineage_areas : bool
            If False then communities defined by areas, in which there
            is only one lineage in such a "community" will be skipped, but the
            other communities defined on the tree will be processed. If True,
            then the ENTIRE tree will be skipped if even one community by area
            has this issue.
        drop_trees_with_single_lineage_trait_states : bool
            If False then communities defined by traits, in which there
            is only one lineage in such a "community" will be skipped, but the
            other communities defined on the tree will be processed. If True,
            then the ENTIRE tree will be skipped if even one community by trait
            has this issue.
        """
        self.drop_trees_not_spanning_all_areas = drop_trees_not_spanning_all_areas
        self.drop_trees_not_spanning_multiple_traits = drop_trees_not_spanning_multiple_traits
        self.drop_trees_with_single_lineage_areas = drop_trees_with_single_lineage_areas
        self.drop_trees_with_single_lineage_trait_states = drop_trees_with_single_lineage_trait_states
        if trait_indexes_to_exclude:
            self.trait_indexes_to_exclude = set(trait_indexes_to_exclude)
        else:
            self.trait_indexes_to_exclude = set([])
        if trait_states_to_exclude:
            self.trait_states_to_exclude = set(trait_states_to_exclude)
        else:
            self.trait_states_to_exclude = set([])
        self.run_logger = run_logger
        self.progress_message_frequency_percentage = 1
        self.stat_name_prefix = "predictor"
        self.stat_name_delimiter = "."
        self.num_randomization_replicates = 100
        self._current_tree_idx = None

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

    def summarize_trees(self, trees,):
        processed_trees = []
        summary_fieldnames = set()
        summary_results = []
        trees = list(trees)
        self._current_tree_idx = None
        for tree_idx, tree in enumerate(trees):
            self._current_tree_idx = tree_idx
            try:
                if self.run_logger and not (int(float(tree_idx)/len(trees) * 100) % self.progress_message_frequency_percentage):
                    self.run_logger.info("Tree {} of {}".format( tree_idx+1, len(trees)))
                self.summarize_tree(tree)
                processed_trees.append(tree)
                summary_fieldnames.update(tree.stats.keys())
                summary_results.append(collections.OrderedDict(tree.stats))
            except TreeSummarizer.IncompleteAreaRadiationException:
                self.run_logger.warning("Skipping tree {} (1-based index): Not all areas occupied".format(
                self._current_tree_idx+1))
            except TreeSummarizer.SingleTaxonAreaAssemblageException as e:
                self.run_logger.warning("Skipping tree {}: Area {} (0-based index)has only one lineage".format(
                self._current_tree_idx+1,
                e.area_idx))
            except TreeSummarizer.SingleTaxonTraitStateAssemblageException as e:
                self.run_logger.warning("Skipping tree {}: Trait {}, state {} has only one lineage".format(
                self._current_tree_idx+1,
                e.trait_idx,
                e.trait_state_idx,
                ))
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

        tree.stats = collections.defaultdict(lambda:"NA")
        tree.stats["model_id"] = tree.annotations.get_value("model_id", "NA")
        total_tree_length = 0.0
        total_tree_edges = 0.0
        for nd in tree:
            total_tree_edges += 1.0
            if nd.edge.length is None:
                nd.edge.length = 0.0
            total_tree_length += nd.edge.length
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
        phylogenetic_distance_matrix = tree.phylogenetic_distance_matrix()
        tree.stats.update(self._calc_area_based_stats(
            phylogenetic_distance_matrix=phylogenetic_distance_matrix,
            area_taxa_map=area_taxa_map,
            ))
        tree.stats.update(self._calc_trait_based_stats(
            phylogenetic_distance_matrix=phylogenetic_distance_matrix,
            trait_taxa_map=trait_taxa_map,
            ))
        return tree.stats

    def _calc_area_based_stats(self,
            phylogenetic_distance_matrix,
            area_taxa_map,
            ):
        area_assemblage_memberships, area_assemblage_descriptions = self._get_area_community_regimes(area_taxa_map)
        assert len(area_assemblage_memberships) == len(area_assemblage_descriptions)
        results = self._calc_community_ecology_stats(
                phylogenetic_distance_matrix=phylogenetic_distance_matrix,
                assemblage_memberships=area_assemblage_memberships,
                assemblage_descriptions=area_assemblage_descriptions,
                report_character_state_specific_results=False,
                report_character_class_wide_results=True,
                )
        return results

    def _get_area_community_regimes(self, area_taxa_map):
        assemblage_descriptions = []
        assemblage_memberships = []
        for area_idx in area_taxa_map:
            area_taxa = area_taxa_map[area_idx]
            if len(area_taxa) < 2:
                if not self.drop_trees_with_single_lineage_areas:
                    if self.run_logger:
                        self.run_logger.warning("Skipping statistics calculation for community-by-area for area {} (0-based index) of tree {} (1-based index): only one lineage in area".format(
                            area_idx, self._current_tree_idx+1))
                    continue
                else:
                    raise TreeSummarizer.SingleTaxonAreaAssemblageException(area_idx, self._current_tree_idx+1)
            assemblage_memberships.append( area_taxa )
            regime = {
                "assemblage_basis_class_id": "area",
                "assemblage_basis_state_id": "state{}{}".format(self.stat_name_delimiter, area_idx),
            }
            assemblage_descriptions.append(regime)
        if not assemblage_memberships:
            raise TreeSummarizer.IncompleteAreaRadiationException()
        return assemblage_memberships, assemblage_descriptions

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

    def _get_trait_community_regimes(self, trait_taxa_map):
        assemblage_descriptions = []
        assemblage_memberships = []
        for trait_idx in trait_taxa_map:
            for trait_state_idx in trait_taxa_map[trait_idx]:
                tt = trait_taxa_map[trait_idx][trait_state_idx]
                if len(tt) < 2:
                    if not self.drop_trees_with_single_lineage_trait_states:
                        if self.run_logger:
                            self.run_logger.warning("Skipping statistics calculation for community-by-trait-state for trait {}, state {} of tree {}: only one lineage with trait state".format(
                                trait_idx+USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE, trait_state_idx, self._current_tree_idx+1))
                        continue
                    else:
                        raise TreeSummarizer.SingleTaxonTraitStateAssemblageException(trait_idx+USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE, trait_state_idx, self._current_tree_idx+1)
                assemblage_memberships.append(tt)
                regime = {
                    "assemblage_basis_class_id": "trait{}{}".format(self.stat_name_delimiter, trait_idx + USER_SPECIFIED_TRAIT_TYPE_INDEX_START_VALUE),
                    "assemblage_basis_state_id": "state{}{}".format(self.stat_name_delimiter, trait_state_idx),
                }
                assemblage_descriptions.append(regime)
        if not assemblage_memberships:
            raise TreeSummarizer.IncompleteTraitRaditionException()
        return assemblage_memberships, assemblage_descriptions

    def _calc_community_ecology_stats(self,
            phylogenetic_distance_matrix,
            assemblage_memberships,
            assemblage_descriptions,
            report_character_state_specific_results=True,
            report_character_class_wide_results=True,
            ):

        assert len(assemblage_descriptions) == len(assemblage_memberships)

        summary_statistics_suite = {}
        results_by_character_class = {}
        stat_scores_to_be_harvested = ("obs", "z", "p",) # z = score, p = p-value (turns out this is quite informative)
        for sstbh in stat_scores_to_be_harvested:
            results_by_character_class[sstbh] = collections.defaultdict(list)

        for edge_weighted_desc in ("unweighted", "weighted"):
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
                for result, assemblage_desc in zip(results_group, assemblage_descriptions):
                    for ses_result_statistic in stat_scores_to_be_harvested:
                        character_class_statistic_prefix = self.stat_name_delimiter.join([
                            self.stat_name_prefix,
                            "community",
                            "by",
                            assemblage_desc["assemblage_basis_class_id"],
                            ])
                        statistic_subtype_desc = self.stat_name_delimiter.join([
                            edge_weighted_desc,
                            underlying_statistic_type_desc,
                            # assemblage_desc["assemblage_basis_state_id"],
                            ])
                        character_class_statistic_key = tuple([character_class_statistic_prefix, statistic_subtype_desc])
                        ses_result_statistic_value = getattr(result, ses_result_statistic)
                        if ses_result_statistic_value is None:
                            continue
                        if report_character_state_specific_results:
                            character_state_statistic_name = self.stat_name_delimiter.join([
                                character_class_statistic_prefix,
                                assemblage_desc["assemblage_basis_state_id"],
                                statistic_subtype_desc,
                                ses_result_statistic,
                                ])
                            assert character_state_statistic_name not in summary_statistics_suite
                            summary_statistics_suite[character_state_statistic_name] = ses_result_statistic_value
                        if report_character_class_wide_results:
                            results_by_character_class[ses_result_statistic][character_class_statistic_key].append(ses_result_statistic_value)
        if report_character_class_wide_results:
            for ses_result_statistic in results_by_character_class:
                if len(results_by_character_class[ses_result_statistic]) == 0:
                    continue
                for key in results_by_character_class[ses_result_statistic]:
                    character_class_statistic_prefix, statistic_subtype_desc = key
                    svalues = results_by_character_class[ses_result_statistic][key]
                    mean_var = statistics.mean_and_sample_variance(svalues)
                    for s, sdesc in zip( mean_var, ("mean", "var"), ):
                        sn_title = self.stat_name_delimiter.join([
                            character_class_statistic_prefix,
                            sdesc,
                            statistic_subtype_desc,
                            ses_result_statistic,
                            ])
                        assert sn_title not in summary_statistics_suite
                        summary_statistics_suite[sn_title] = s
        return summary_statistics_suite
