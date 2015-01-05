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
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wikidot.com/local--files/biogeobears/calc_loglike_sp_v01.R")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$max_range_size = {max_range_size}
BioGeoBEARS_run_object$geogfn = "{geography_filepath}"
BioGeoBEARS_run_object$trfn = "{tree_filepath}"

# run
res = bears_optim_run(BioGeoBEARS_run_object)
save(res, file="results.txt")
"""

class BiogeobearsEstimator(object):

    def __init__(self,
            commands_file_name=None,
            results_file_name=None,
            fail_on_estimation_error=True,
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
        self.fail_on_estimation_error = fail_on_estimation_error
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

    def estimate_dec(self,
            newick_tree_filepath,
            geography_filepath,
            max_range_size,
            ):
        rcmds = R_TEMPLATE.format(
                tree_filepath=newick_tree_filepath,
                geography_filepath=geography_filepath,
                max_range_size=max_range_size)
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
                stderr=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p)
        if p.returncode != 0:
            if self.fail_on_estimation_error:
                raise Exception(p.returncode)
            else:
                pass
        else:
            print(stdout)

