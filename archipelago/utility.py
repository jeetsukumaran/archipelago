#! /usr/bin/env python

##############################################################################
##
##  Copyright 2010-2014 Jeet Sukumaran.
##  All rights reserved.
##
##  Reoccurrences and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met:
##
##      * Reoccurrencess of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
##      * Reoccurrencess in binary form must reproduce the above copyright
##        notice, this list of conditions and the following disclaimer in the
##        documentation and/or other materials provided with the occurrences.
##      * The names of its contributors may not be used to endorse or promote
##        products derived from this software without specific prior written
##        permission.
##
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
##  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
##  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
##  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN OR MARK T. HOLDER
##  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
##  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
##  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
##  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
##  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
##  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
##  POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import os
import logging
import inspect
import tempfile
import pprint
import collections
import dendropy
import json

_LOGGING_LEVEL_ENVAR = "ARCHIPELAGO_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "ARCHIPELAGO_LOGGING_FORMAT"

def dump_stack():
    for frame, filename, line_num, func, source_code, source_index in inspect.stack()[2:]:
        if source_code is None:
            print("{}: {}".format(filename, line_num))
        else:
            print("{}: {}: {}".format(filename, line_num, source_code[source_index].strip()))

def encode_lineage(node, exclude_areas=None):
    traits = "".join(str(i) for i in node.traits_vector)
    if exclude_areas is None:
        areas = "".join(str(i) for i in node.distribution_vector)
    else:
        areas = "".join(str(i) for idx, i in enumerate(node.distribution_vector) if idx not in exclude_areas)
    return "s{}.{}.{}".format(node.index, traits, areas)

def read_model_from_json_path(filepath):
    return json.load(open(filepath, "r"))

def read_model_from_python_path(filepath):
    src = open(filepath, "rb").read()
    return eval(src)

class IndexGenerator(object):

    def __init__(self, start=0):
        self.start = start
        self.index = start

    def __next__(self):
        c = self.index
        self.index += 1
        return c
    next = __next__

    def reset(self, start=None):
        if start is None:
            start = self.start
        self.index = start

class OutOfRegionError(Exception):
    pass

class ColorAssigner(object):

    COLORS = (
        "#ff0000", "#00ff00", "#0000ff", "#ffff00", "#ff00ff", "#00ffff", "#000000",
        "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#808080",
        "#c00000", "#00c000", "#0000c0", "#c0c000", "#c000c0", "#00c0c0", "#c0c0c0",
        "#400000", "#004000", "#000040", "#404000", "#400040", "#004040", "#404040",
        "#200000", "#002000", "#000020", "#202000", "#200020", "#002020", "#202020",
        "#600000", "#006000", "#000060", "#606000", "#600060", "#006060", "#606060",
        "#a00000", "#00a000", "#0000a0", "#a0a000", "#a000a0", "#00a0a0", "#a0a0a0",
        "#e00000", "#00e000", "#0000e0", "#e0e000", "#e000e0", "#00e0e0", "#e0e0e0",
        "#666666",
        )
    COLOR_INDEXES = {}
    for idx, c in enumerate(COLORS):
        COLOR_INDEXES[c] = idx

    def __init__(self, offset=0):
        self.assigned_colors = {}
        self.offset = offset

    def __getitem__(self, code):
        try:
            return self.assigned_colors[code]
        except KeyError:
            on_bits = []
            for idx, i in enumerate(code):
                if i == "1":
                    on_bits.append(idx)
            if len(on_bits) > 1:
                self.assigned_colors[code] = "#666666"
            elif len(on_bits) == 0:
                raise OutOfRegionError
            else:
                self.assigned_colors[code] = self.COLORS[on_bits[0]+self.offset]
            return self.assigned_colors[code]

class NameToSymbolMap(object):

    SYMBOLS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    SYMBOL_INDEXES = {}
    for idx, s in enumerate(SYMBOLS):
        SYMBOL_INDEXES[s] = idx

    def __init__(self):
        self.assigned_symbols = {}

    def __getitem__(self, name):
        try:
            return self.assigned_symbols[name]
        except KeyError:
            sidx = len(self.assigned_symbols)
            assert sidx < len(NameToSymbolMap.SYMBOLS)
            self.assigned_symbols[name] = NameToSymbolMap.SYMBOLS[sidx]
            return self.assigned_symbols[name]

    def __iter__(self):
        for idx in range(len(self.assigned_symbols)):
            return NameToSymbolMap.SYMBOLS[idx]

    def __len__(self):
        return len(self.assigned_symbols)

class TreeProcessor(object):

    def __init__(self,
            areas_to_exclude=None,
            drop_stunted_trees=False,
            ):
        """
        Parameters
        ----------
        areas_to_exclude : list of integers
            List of (0-based) indexes of areas to exclude from analysis.
            Incidences on these areas will be ignored, and reference to these
            areas will be removed from taxon labels before decoding. Any
            trees that do not include tips that occur on areass that are not
            excluded will be removed from trees.
        drop_stunted_trees : bool
            If `True`, trees with less than 2 tips will be removed.
        """
        self.area_colors = ColorAssigner()
        self.habitat_colors = ColorAssigner()
        self.drop_stunted_trees = drop_stunted_trees
        self.areas_to_exclude = areas_to_exclude

    def decode_labeled_item_biogeography(self, t):
        if hasattr(t, "is_encoded") and t.is_encoded:
            return
        if areas_to_exclude:
            p = t.label.split(".")
            t.label = ".".join([p[0], p[1][1:], p[2]])
        label_parts = t.label.split(".")
        t.area_code = label_parts[1]
        # if self.exclude_first_area_as_continental_source_outside_of_analysis:
        #     t.area_code = t.area_code[1:]
        try:
            t.area_color = self.area_colors[t.area_code]
            t.habitat_code = label_parts[2]
            t.habitat_color = self.habitat_colors[t.habitat_code]
            t.is_encoded = True
            t.out_of_region = False
        except OutOfRegionError:
            # taxon only occurs on first area, and needs to be removed
            t.out_of_region = True

    def decode_taxon_biogeography(self, taxa):
        if hasattr(taxa, "is_encoded") and taxa.is_encoded:
            return []
        taxa_to_exclude = []
        for t in taxa:
            self.decode_labeled_item_biogeography(t)
            if t.out_of_region:
                taxa_to_exclude.append(t)
        taxa.is_encoded = True
        return taxa_to_exclude

    def prune_trees(self, trees, taxa_to_exclude):
        if not taxa_to_exclude:
            return trees
        to_keep = dendropy.TreeList(taxon_namespace=trees.taxon_namespace)
        for tree in trees:
            try:
                tree.prune_taxa(taxa_to_exclude)
                if not self.drop_stunted_trees or tree.seed_node.num_child_nodes() > 1:
                    to_keep.append(tree)
            except AttributeError:
                # trying to prune root node
                pass
        for taxon in taxa_to_exclude:
            trees.taxon_namespace.remove_taxon(taxon)
        return to_keep

    def write_colorized_trees(self,
            outf,
            trees,
            schema,
            decode_trees):
        if schema == "by-area":
            taxon_color_attr = "area_color"
            branch_color_attr = "habitat_color"
        elif schema == "by-habitat":
            taxon_color_attr = "habitat_color"
            branch_color_attr = "area_color"
        else:
            raise ValueError("Unrecognized schema: {}".format(schema))
        if decode_trees:
            trees = self.decode_trees(trees)
        for taxon in trees.taxon_namespace:
            taxon.annotations["!color"] = getattr(taxon, taxon_color_attr)
        trees.write_to_stream(outf, "nexus")

    def decode_trees(self, trees):
        taxa_to_exclude = self.decode_taxon_biogeography(trees.taxon_namespace)
        trees = self.prune_trees(trees, taxa_to_exclude)
        return trees

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        self.handlers = []
        if kwargs.get("log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            handler1.setLevel(stderr_logging_level)
            handler1.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler1)
            self.handlers.append(handler1)
        if kwargs.get("log_to_file", True):
            if "log_stream" in kwargs:
                log_stream = kwargs.get("log_stream")
            else:
                log_stream = open(kwargs.get("log_path", self.name + ".log"), "w")
            handler2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            handler2.setLevel(file_logging_level)
            handler2.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler2)
            self.handlers.append(handler2)
        self._system = None

    def _get_system(self):
        return self._system

    def _set_system(self, system):
        self._system = system
        if self._system is None:
            for handler in self.handlers:
                handler.setFormatter(self.get_default_formatter())
        else:
            for handler in self.handlers:
                handler.setFormatter(self.get_simulation_generation_formatter())

    system = property(_get_system, _set_system)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simulation_generation_formatter(self):
        # f = logging.Formatter("[%(asctime)s] t = %(elapsed_time)10.6f: %(message)s")
        f = logging.Formatter("[%(asctime)s] %(simulation_time)s%(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def supplemental_info_d(self):
        if self._system is not None:
            # return {
            #   "simulation_time" : "[t = {:10.6f}] ".format(self._system.elapsed_time),
            # }
            if self._system.elapsed_time == 0:
                return {
                        "simulation_time" : "Setup: ",
                        }
            else:
                return {
                        "simulation_time" : "[t = {:13.6f}] ".format(self._system.elapsed_time),
                        }
        else:
            return None

    def debug(self, msg):
        self._log.debug("[DEBUG] {}".format(msg), extra=self.supplemental_info_d())

    def info(self, msg):
        self._log.info(msg, extra=self.supplemental_info_d())

    def warning(self, msg):
        self._log.warning(msg, extra=self.supplemental_info_d())

    def error(self, msg):
        self._log.error(msg, extra=self.supplemental_info_d())

    def critical(self, msg):
        self._log.critical(msg, extra=self.supplemental_info_d())

