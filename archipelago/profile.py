#! /usr/bin/env python

import collections
from dendropy.model import birthdeath
import dendropy
from archipelago import estimate

class TreeProfiler(object):

    def __init__(self):
        self.trait_transition_rate_estimator = estimate.TraitEvolutionRateEstimator()

    def diagnose_num_trait_types(self,
            trees,
            is_trees_decoded=False,
            is_suppressed_taxa=False):
        if not is_trees_decoded:
            model.ArchipelagoModel.decode_tree_lineages_from_labels(
                    trees=trees,
                    is_suppressed_taxa=is_suppressed_taxa)
        sample_node = next(trees.leaf_node_iter())
        num_trait_types = len(sample_node.traits_vector)
        trees.num_trait_types = num_trait_types

    def estimate_pure_birth(self, trees, tree_results_map):
        for tree in trees:
            try:
                bdfit = birthdeath.fit_pure_birth_model_to_tree(tree)
            except ValueError:
                pass
            try:
                tree_results_map[tree]["pure.birth.rate"] = bdfit["birth_rate"]
            except KeyError:
                tree_results_map[tree] = {"pure.birth.rate": bdfit["birth_rate"]}
        return tree_results_map

    def estimate_trait_transition_rates(self,
            trees,
            tree_results_map,
            trait_estimated_transition_rate_field_names,
            is_trees_decoded=False,
            is_suppressed_taxa=False,
            ):
        self.trait_transition_rate_estimator.estimate_trait_evolution_rate(
                trees=trees,
                tree_results_map=tree_results_map,
                trait_estimated_transition_rate_field_names=trait_estimated_transition_rate_field_names,
                is_trees_decoded=is_trees_decoded,
                is_suppressed_taxa=is_suppressed_taxa,)




