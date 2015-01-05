#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import os
import re
import collections
import tempfile
import subprocess
import dendropy
from dendropy.utility import processio

R_TEMPLATE = """\
library(optimx)
library(FD)
# library(snow)
library(parallel)
library(BioGeoBEARS)
{patch_code}
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
{param_settings}
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$max_range_size = {max_range_size}
BioGeoBEARS_run_object$geogfn = "{geography_filepath}"
BioGeoBEARS_run_object$trfn = "{tree_filepath}"
res = bears_optim_run(BioGeoBEARS_run_object)
sink("{results_filepath}")
res$outputs@params_table
sink()
"""

PARAM_SETTING_TEMPLATE = """\
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["{param_name}","{param_aspect}"] = {value}
"""

class BiogeobearsEstimator(object):

    def __init__(self,
            commands_file_name=None,
            results_file_name=None,
            fail_on_estimation_error=True,
            debug_mode=False,
            ):
        if commands_file_name is None:
            self.commands_file = tempfile.NamedTemporaryFile()
            self.commands_file_name = self.commands_file.name
        else:
            self.commands_file_name = commands_file_name
        if results_file_name is None:
            self.results_file = tempfile.NamedTemporaryFile()
            self.results_file_name = self.results_file.name
        else:
            self.results_file_name = results_file_name
        # self.results_file_name = "debugbgb.txt"
        self.fail_on_estimation_error = fail_on_estimation_error
        self.debug_mode = debug_mode
        self.path_to_libexec = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "libexec")
        self.patch_filenames = [
                "BioGeoBEARS_basics_v1.R",
                "BioGeoBEARS_generics_v1.R",
                "BioGeoBEARS_classes_v1.R",
                "BioGeoBEARS_models_v1.R",
                "BioGeoBEARS_plots_v1.R",
                "BioGeoBEARS_readwrite_v1.R",
                "BioGeoBEARS_simulate_v1.R",
                "BioGeoBEARS_stratified_v1.R",
                "BioGeoBEARS_univ_model_v1.R",
                "calc_loglike_sp_v01.R",
                ]
        self.patch_filepaths = [os.path.join(self.path_to_libexec, f) for f in self.patch_filenames]
        self.patch_code = "\n".join(["source('{}')".format(f) for f in self.patch_filepaths])

    def estimate_dec(self,
            newick_tree_filepath,
            geography_filepath,
            max_range_size,
            **kwargs
            ):

        param_settings = []
        for param_name in ("b", "e", "d"):
            if "fixed_" + param_name in kwargs:
                for param_aspect in ("min", "max", "init", "est"):
                    param_settings.append(PARAM_SETTING_TEMPLATE.format(
                        param_name="b", param_aspect=param_aspect, value=kwargs["fixed_"+param_name]))
            else:
                for param_aspect in ("min_", "max_", "init_", "est_"):
                    if param_aspect + param_name in kwargs:
                        param_settings.append(PARAM_SETTING_TEMPLATE.format(
                            param_name="b", param_aspect=param_aspect[:-1], value=kwargs[param_aspect+param_name]))
        param_settings = "\n".join(param_settings)
        rcmds = R_TEMPLATE.format(
            patch_code=self.patch_code,
            param_settings=param_settings,
            tree_filepath=newick_tree_filepath,
            geography_filepath=geography_filepath,
            max_range_size=max_range_size,
            results_filepath=self.results_file_name,
            )
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
                stdout=subprocess.PIPE if not self.debug_mode else None,
                stderr=subprocess.PIPE if not self.debug_mode else None,
                )
        stdout, stderr = processio.communicate(p)
        if p.returncode != 0:
            if self.fail_on_estimation_error:
                raise Exception("Non-zero return code: {}\n{}\n{}".format(
                        p.returncode,
                        stdout,
                        stderr,
                        ))
            else:
                return None
        results_rows = open(self.results_file_name, "r").read().split("\n")
        results_table = collections.OrderedDict()
        for row in results_rows[1:21]:
            cols = row.split()
            if cols[0] == "desc" or cols[0] == "note":
                break
            try:
                results_table[cols[0]] = float(cols[5])
            except IndexError:
                raise IndexError(cols)
        return results_table

