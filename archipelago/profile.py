#! /usr/bin/env python

import collections
from dendropy.model import birthdeath
import dendropy

class TreeProfiler(object):

    def estimate_pure_birth(self, trees, result_map):
        for tree in trees:
            try:
                bdfit = birthdeath.fit_pure_birth_model_to_tree(tree)
            except ValueError:
                pass
            try:
                result_map[tree]["pure.birth.rate"] = bdfit["birth_rate"]
            except KeyError:
                result_map[tree] = {"pure.birth.rate": bdfit["birth_rate"]}
        return result_map



