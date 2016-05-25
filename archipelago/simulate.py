#! /usr/bin/env python

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import random
import collections
import argparse
import pprint
import copy
import json
from distutils.util import strtobool

import dendropy
from dendropy.utility import textprocessing

import archipelago
from archipelago import eventlog
from archipelago import model
from archipelago import utility
from archipelago import error

class Snapshot(object):

    def __init__(self, simulator):
        self.start_time = simulator.elapsed_time
        self.end_time = None
        self.duration = None
        all_areas_complete_tree_sio = StringIO()
        all_areas_extant_tree_sio = StringIO()
        focal_areas_tree_sio = StringIO()
        histories_sio = StringIO()
        simulator.store_sample(
                all_areas_complete_tree_out=all_areas_complete_tree_sio,
                all_areas_extant_tree_out=all_areas_extant_tree_sio,
                focal_areas_tree_out=focal_areas_tree_sio,
                histories_out=histories_sio,
                )
        self.all_areas_complete_tree_str = all_areas_complete_tree_sio.getvalue()
        self.all_areas_extant_tree_str = all_areas_extant_tree_sio.getvalue()
        self.focal_areas_tree_str = focal_areas_tree_sio.getvalue()
        if simulator.event_log:
            histories_str = histories_sio.getvalue()
            self.history_d = json.loads(histories_str)
        else:
            self.history_d = None

    def close(self, end_time):
        self.end_time = end_time
        self.duration = self.end_time - self.start_time

    def compose_snapshot(self, time_to_add):
        # This must be some of the ugliest code I've ever written
        # Basically, just to add some time to the frozen histories
        results = {}
        if self.history_d:
            taxon_namespace = dendropy.TaxonNamespace(self.history_d["leaf_labels"])
        else:
            taxon_namespace = dendropy.TaxoNamespace()
        all_areas_complete_tree = dendropy.Tree.get(
                data=self.all_areas_complete_tree_str,
                schema="newick",
                taxon_namespace=taxon_namespace,
                rooting="force-rooted",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        eventlog.EventLog.prepare_tree_for_event_serialization(
                all_areas_complete_tree,
                time_to_add_to_extant_tips=time_to_add,
                is_set_taxa=False)
        results["all-areas.complete"] = all_areas_complete_tree.as_string("newick",
                suppress_internal_node_labels=False,
                suppress_leaf_node_labels=False,
                )

        for tree_desc, tree_str in (
                ("all-areas.extant", self.all_areas_extant_tree_str,),
                ("focal-areas", self.focal_areas_tree_str,),
                ):
            tree = dendropy.Tree.get(
                    data=tree_str,
                    schema="newick",
                    rooting="force-rooted",
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True,
                    )
            for nd in tree.leaf_node_iter():
                nd.edge.length += time_to_add
            results[tree_desc] = tree.as_string("newick",
                    suppress_internal_node_labels=False,
                    suppress_leaf_node_labels=False,
                    )

        if not self.history_d:
            rvals.append(None)
        else:
            self.history_d["tree"] = eventlog.EventLog.compose_tree_data(
                    tree=all_areas_complete_tree,
                    max_event_time=self.history_d["tree"]["max_event_time"])
            # self.history_d["lineages"] = eventlog.EventLog.compose_lineage_definitions(tree=all_areas_complete_tree)
            # assert len(self.history_d["leaf_labels"]) == len(all_areas_complete_tree.taxon_namespace)
            for lineage_d in self.history_d["lineages"]:
                if lineage_d["is_leaf"]:
                    lineage_d["lineage_end_time"] += time_to_add
                    lineage_d["lineage_duration"] = lineage_d["lineage_end_time"] - lineage_d["lineage_start_time"]
            results["history"] = json.dumps(self.history_d, indent=4, separators=(',', ': '))
        return results

class ArchipelagoSimulator(object):

    @staticmethod
    def get_fixed_value_function(v, description):
        f = lambda x: v
        f.__doc__ = description
        return f

    @staticmethod
    def compose_all_areas_complete_trees_filepath(output_prefix):
        return output_prefix + ".all-areas.complete.trees"

    @staticmethod
    def compose_all_areas_extant_trees_filepath(output_prefix):
        return output_prefix + ".all-areas.extant.trees"

    @staticmethod
    def compose_focal_areas_trees_filepath(output_prefix):
        return output_prefix + ".focal-areas.trees"

    @staticmethod
    def compose_histories_filepath(output_prefix):
        return output_prefix + ".histories.json"

    @staticmethod
    def simple_node_label_function(node):
        return "s{}".format(node.index)

    def __init__(self,
            archipelago_model,
            config_d,
            is_verbose_setup):

        # configure
        self.elapsed_time = 0.0 # need to be here for logging
        config_d = dict(config_d) # make copy so we can pop items
        self.configure_simulator(config_d, verbose=is_verbose_setup)

        # set up model
        self.model = archipelago_model

        # set up geography
        self.geography = self.model.new_geography()

        ### Initialize event log
        if self.is_store_histories:
            self.event_log = eventlog.EventLog()
        else:
            self.event_log = None

        # initialize phylogeny
        self.phylogeny = model.Phylogeny(
                model=self.model,
                geography=self.geography,
                rng=self.rng,
                log_event=self.log_event if self.event_log is not None else None,
                debug_mode=self.debug_mode,
                run_logger=self.run_logger,
                )
        self.phylogeny.bootstrap()

        # begin logging generations
        self.run_logger.system = self

        # gsa snapshots
        self.snapshots = []

    def configure_simulator(self, config_d, verbose=True):

        self.name = config_d.pop("name", None)
        if self.name is None:
            self.name = str(id(self))
        self.output_prefix = config_d.pop("output_prefix", "archipelago-{}".format(self.name))

        self.run_logger = config_d.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(
                    name="archipelago",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_path=self.output_prefix + ".log",
                    file_logging_level=config_d.pop("file_logging_level", "info"),
                    )
        self.run_logger.system = self

        if verbose:
            self.run_logger.info("Configuring simulation '{}'".format(self.name))

        if config_d.pop("store_all_areas_complete_trees", True):
            self.all_areas_complete_trees_file = config_d.pop("all_areas_complete_trees_file", None)
            if self.all_areas_complete_trees_file is None:
                self.all_areas_complete_trees_file = open(ArchipelagoSimulator.compose_all_areas_complete_trees_filepath(self.output_prefix), "w")
            if verbose:
                self.run_logger.info("All areas complete trees filepath: {}".format(self.all_areas_complete_trees_file.name))
        else:
            self.all_areas_complete_trees_file = None
            if verbose:
                self.run_logger.info("All areas complete trees will not be stored")

        if config_d.pop("store_all_areas_extant_trees", True):
            self.all_areas_extant_trees_file = config_d.pop("all_areas_extant_trees_file", None)
            if self.all_areas_extant_trees_file is None:
                self.all_areas_extant_trees_file = open(ArchipelagoSimulator.compose_all_areas_extant_trees_filepath(self.output_prefix), "w")
            if verbose:
                self.run_logger.info("All areas extant lineage trees filepath: {}".format(self.all_areas_extant_trees_file.name))
        else:
            self.all_areas_extant_trees_file = None
            if verbose:
                self.run_logger.info("All areas extantlineage trees will not be stored")

        if config_d.pop("store_focal_areas_trees", True):
            self.focal_areas_trees_file = config_d.pop("focal_areas_trees_file", None)
            if self.focal_areas_trees_file is None:
                self.focal_areas_trees_file = open(ArchipelagoSimulator.compose_focal_areas_trees_filepath(self.output_prefix), "w")
            if verbose:
                self.run_logger.info("Focal area extant lineage trees filepath: {}".format(self.focal_areas_trees_file.name))
        else:
            self.focal_areas_trees_file = None
            if verbose:
                self.run_logger.info("Focal area extant lineage trees will not be stored")

        if config_d.pop("store_histories", False):
            self.is_store_histories = True
            self.histories_file = config_d.pop("histories_file", None)
            if self.histories_file is None:
                self.histories_file = open(ArchipelagoSimulator.compose_histories_filepath(self.output_prefix), "w")
            if verbose:
                self.run_logger.info("Event histories will be written out to: {}".format(self.histories_file.name))
        else:
            self.is_store_histories = False
            self.histories_file = None
            if verbose:
                self.run_logger.info("Event histories will not be written out")

        self.is_suppress_internal_node_labels = config_d.pop("suppress_internal_node_labels", False)
        if verbose:
            self.run_logger.info("Internal node labels will{} be written on trees".format(" not" if self.is_suppress_internal_node_labels else ""))

        self.is_encode_nodes = config_d.pop("encode_nodes", True)
        if verbose:
            if self.is_encode_nodes:
                self.run_logger.info("Trait states and ranges will be encoded on node labels")
            else:
                self.run_logger.info("Trait states and ranges will NOT be encoded on node labels")

        self.is_annotate_nodes = config_d.pop("annotate_nodes", False)
        if verbose:
            if self.is_annotate_nodes:
                self.run_logger.info("Trait states and ranges will be annotated on node labels")
            else:
                self.run_logger.info("Trait states and ranges will NOT be annotated on node labels")

        self.rng = config_d.pop("rng", None)
        if self.rng is None:
            self.random_seed = config_d.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            if verbose:
                self.run_logger.info("Initializing with random seed {}".format(self.random_seed))
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in config_d:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            if verbose:
                self.run_logger.info("Using existing random number generator")

        self.log_frequency = config_d.pop("log_frequency", None)

        self.debug_mode = config_d.pop("debug_mode", False)
        if verbose and self.debug_mode:
            self.run_logger.info("Running in DEBUG mode")

        if config_d.pop("store_model_description", True):
            self.model_description_file = config_d.pop("model_description_file", None)
            if self.model_description_file is None:
                self.model_description_file = open(self.output_prefix + ".model.log.json", "w")
            if verbose:
                self.run_logger.info("Model description filepath: {}".format(self.model_description_file.name))
        else:
            self.model_description_file = None
            if verbose:
                self.run_logger.info("Model description will not be stored")

        if config_d:
            raise TypeError("Unsupported configuration keywords: {}".format(config_d))

    def run(self):

        ### Save model
        if self.model_description_file is not None:
            self.model.write_model(self.model_description_file)

        ### Initialize time
        self.elapsed_time = 0.0

        ### Initialize logging
        ### None: default logging, 0: no logging
        if self.log_frequency is None:
            if self.model.target_focal_area_lineages:
                default_log_frequency = 1
            else:
                default_log_frequency = self.model.max_time/100
        if self.log_frequency:
            if self.model.target_focal_area_lineages:
                last_logged_num_tips = 0
            else:
                last_logged_time = 0.0

        ### Initialize debugging
        if self.debug_mode:
            num_events = 0

        ### Initialize termination conditiong checking
        ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
        ntips = len(self.phylogeny.current_lineages)

        ### Keep track of snapshot for GSA termination condition
        current_snapshot = None

        while True:

            ### DEBUG
            if self.debug_mode:
                num_events += 1
                # self.run_logger.debug("Pre-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.run_logger.debug("Post-event {}: debug check: current lineages = {}".format(num_events, len(self.phylogeny.current_lineages)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    assert lineage.areas

            ### LOGGING
            if self.log_frequency:
                if self.model.target_focal_area_lineages:
                    if ntips_in_focal_areas - last_logged_num_tips >= self.log_frequency:
                        last_logged_num_tips = ntips_in_focal_areas
                        self.run_logger.info("{} lineages occurring in focal areas, {} lineages across all areas".format(ntips_in_focal_areas, ntips))
                else:
                    if self.elapsed_time - last_logged_time >= self.log_frequency:
                        last_logged_time = self.elapsed_time
                        self.run_logger.info("{} lineages occurring in focal areas, {} lineages across all areas".format(ntips_in_focal_areas, ntips))

            ### EVENT SCHEDULING
            event_calls, event_rates = self.schedule_events()
            sum_of_event_rates = sum(event_rates)
            if self.debug_mode:
                if sum_of_event_rates == 0:
                    self.run_logger.debug("Sum of event rates is 0: {}".format(event_rates))

            time_till_event = self.rng.expovariate(sum_of_event_rates)
            self.elapsed_time += time_till_event
            if self.model.max_time and self.elapsed_time > self.model.max_time:
                time_to_add = self.model.max_time - self.elapsed_time
                self.elapsed_time = self.model.max_time
                for lineage in self.phylogeny.iterate_current_lineages():
                    lineage.edge.length += time_to_add
                self.run_logger.info("Termination condition of t = {} reached: storing results and terminating".format(self.elapsed_time))
                self.store_sample(
                    all_areas_complete_tree_out=self.all_areas_complete_trees_file,
                    all_areas_extant_tree_out=self.all_areas_extant_trees_file,
                    focal_areas_tree_out=self.focal_areas_trees_file,
                    histories_out=self.histories_file,
                    )
                break
            else:
                for lineage in self.phylogeny.iterate_current_lineages():
                    lineage.edge.length += time_till_event

            ### EVENT SELECTION AND EXECUTION
            event_idx = model.weighted_index_choice(
                    weights=event_rates,
                    sum_of_weights=sum_of_event_rates,
                    rng=self.rng)
            if self.debug_mode:
                self.run_logger.debug("Event {}: {}".format(num_events, event_calls[event_idx]))

            try:
                event_calls[event_idx][0](**event_calls[event_idx][1])
            except model.Lineage.NullDistributionException as e:
                self.phylogeny.extinguish_lineage(e.lineage, extinction_type="null area set")

            ### DEBUG
            if self.debug_mode:
                # self.run_logger.debug("Post-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.run_logger.debug("Post-event {}: debug check: current lineages = {}".format(num_events, len(self.phylogeny.current_lineages)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    assert lineage.areas

            ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
            ntips = len(self.phylogeny.current_lineages)
            if self.model.gsa_termination_focal_area_lineages:
                if ntips_in_focal_areas >= self.model.gsa_termination_focal_area_lineages:
                    # select/process one of the previously stored snapshots, write to final results file,
                    # and then break
                    if current_snapshot is not None:
                        current_snapshot.close(end_time=self.elapsed_time)
                    self.select_and_save_gsa_snapshot()
                    break
                elif ntips_in_focal_areas == self.model.target_focal_area_lineages:
                    if current_snapshot is not None:
                        current_snapshot.close(end_time=self.elapsed_time)
                    current_snapshot = Snapshot(self)
                    self.snapshots.append(current_snapshot)
                elif current_snapshot is not None and ntips_in_focal_areas != self.model.gsa_termination_focal_area_lineages:
                    current_snapshot.close(end_time=self.elapsed_time)
                    current_snapshot = None
            elif self.model.target_focal_area_lineages and ntips_in_focal_areas >= self.model.target_focal_area_lineages:
                self.run_logger.info("Termination condition of {} lineages in focal areas reached at t = {}: storing results and terminating".format(self.model.target_focal_area_lineages, self.elapsed_time))
                self.store_sample(
                    all_areas_complete_tree_out=self.all_areas_complete_trees_file,
                    all_areas_extant_tree_out=self.all_areas_extant_trees_file,
                    focal_areas_tree_out=self.focal_areas_trees_file,
                    histories_out=self.histories_file,
                    )
                break

    def log_event(self, **kwargs):
        self.event_log.register_event(
                event_time=self.elapsed_time,
                **kwargs)

    def schedule_events(self):
        master_event_calls = []
        master_event_rates = []
        event_fluxes = {}
        event_calls = {}
        event_weights = {}
        for event_type in ("birth", "death", "area_gain", "area_loss"):
            event_fluxes[event_type] = 0.0
            event_calls[event_type] = []
            event_weights[event_type] = []

        # if self.debug_mode:
        #     num_current_lineages = len(self.phylogeny.current_lineages)
        #     self.run_logger.debug("Scheduling events for {} current lineages".format(
        #         num_current_lineages))

        for lineage in self.phylogeny.iterate_current_lineages():
            if self.debug_mode:
                assert lineage.is_extant
                assert lineage.areas
                lineage.debug_check()

            # speciation
            event_fluxes["birth"] += self.model.mean_birth_rate
            for area in lineage.areas:
                birth_weight = self.model.lineage_birth_weight_function(lineage=lineage, area=area)
                if birth_weight:
                    event_calls["birth"].append( (self.phylogeny.split_lineage, {"lineage": lineage, "area": area}) )
                    event_weights["birth"].append(birth_weight)

            # global extinction
            event_fluxes["death"] += self.model.mean_death_rate
            if self.model.mean_death_rate:
                for area in lineage.areas:
                    death_weight = self.model.lineage_death_weight_function(lineage=lineage, area=area)
                    if death_weight:
                        event_calls["death"].append( (self.phylogeny.extinguish_lineage, {"lineage": lineage, "extinction_type": "death"}) )
                        event_weights["death"].append(death_weight)

            # trait evolution
            for trait_idx, current_state_idx in enumerate(lineage.traits_vector):
                for proposed_state_idx in range(self.model.trait_types[trait_idx].nstates):
                    if proposed_state_idx == current_state_idx:
                        continue
                    trait_transition_rate = self.model.trait_types[trait_idx].transition_rate_matrix[current_state_idx][proposed_state_idx]
                    if trait_transition_rate:
                        master_event_calls.append( (self.phylogeny.evolve_trait, {"lineage": lineage, "trait_idx": trait_idx, "state_idx": proposed_state_idx}) )
                        master_event_rates.append(trait_transition_rate)

            # Dispersal/Area Gain

            ## check ##
            # for area in self.geography.areas:
            #     if area not in lineage.areas:
            #         event_calls.append((lineage.add_area, {"area": area}))
            #         event_rates.append(self.model.global_area_gain_rate)

            area_gain_event_parameters, area_gain_event_rates, area_gain_rates_marginalized_by_destination_area = self.geography.calculate_raw_area_gain_events(
                    lineage=lineage,
                    lineage_area_gain_weight_fn=self.model.lineage_area_gain_weight_function,
                    simulation_elapsed_time=self.elapsed_time)
            num_source_areas = len(lineage.areas)
            num_dest_areas = len(self.geography.areas) - num_source_areas
            event_fluxes["area_gain"] += (self.model.global_area_gain_rate * (num_source_areas * num_dest_areas))
            if area_gain_event_rates:
                for area_idx, area_gain_rate in enumerate(area_gain_rates_marginalized_by_destination_area):
                    if area_gain_rate:
                        event_calls["area_gain"].append((lineage.add_area, {"area": self.geography.areas[area_idx], "is_log_event": True}))
                        event_weights["area_gain"].append(area_gain_rate)

            # DEC/local extinction
            event_fluxes["area_loss"] += (self.model.mean_area_loss_rate * len(lineage.areas))
            for area in lineage.areas:
                area_loss_weight = self.model.lineage_area_loss_weight_function(lineage=lineage, area=area)
                if area_loss_weight:
                    event_calls["area_loss"].append( (lineage.remove_area, {"area": area, "is_log_event": True}) )
                    event_weights["area_loss"].append(area_loss_weight)

        for event_type in event_fluxes:
            subevent_flux = event_fluxes[event_type]
            subevent_weights = event_weights[event_type]
            normalization_factor = float(sum(subevent_weights))
            subevent_rates = [ subevent_flux * (w/normalization_factor) for w in subevent_weights ]
            master_event_calls.extend( event_calls[event_type] )
            master_event_rates.extend( subevent_rates )

        return master_event_calls, master_event_rates

    def select_and_save_gsa_snapshot(self):
        if not self.snapshots:
            raise error.FailedToMeetFocalAreaLineageTargetException("Did not generate required number of target lineages in focal areas")
        weights = [s.duration for s in self.snapshots]
        sum_of_weights = sum(weights)
        rnd = self.rng.uniform(0, 1) * sum_of_weights
        diff = None
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                diff = rnd
                break
        selected_snapshot = self.snapshots[i]
        assert diff is not None
        time_to_add = (selected_snapshot.end_time + rnd) - selected_snapshot.start_time
        assert time_to_add >= 0
        self.run_logger.info("GSA Termination condition of {} lineages in focal areas reached at t = {}: selecting snapshot {} (of {}) with {} lineages in focal areas, spanning time from {} to {}, at time {} and terminating".format(
            self.model.gsa_termination_focal_area_lineages,
            self.elapsed_time,
            i+1,
            len(self.snapshots),
            self.model.target_focal_area_lineages,
            selected_snapshot.start_time,
            selected_snapshot.end_time,
            selected_snapshot.start_time + time_to_add,
            ))
        results = selected_snapshot.compose_snapshot(time_to_add=time_to_add)
        if self.all_areas_complete_trees_file is not None:
            self.all_areas_complete_trees_file.write(results["all-areas.complete"])
        if self.all_areas_extant_trees_file is not None:
            self.all_areas_extant_trees_file.write(results["all-areas.extant"])
        if self.focal_areas_trees_file is not None:
            self.focal_areas_trees_file.write(results["focal-areas"])
        if self.histories_file is not None:
            self.histories_file.write(results["history"])

    def store_sample(self,
            all_areas_complete_tree_out,
            all_areas_extant_tree_out,
            focal_areas_tree_out,
            histories_out,
            ):
        results = self.phylogeny.generate_tree_strings_for_serialization(
            is_encode_nodes=self.is_encode_nodes,
            is_annotate_nodes=self.is_annotate_nodes,
            is_suppress_internal_node_labels=self.is_suppress_internal_node_labels,
            )
        if all_areas_complete_tree_out is not None:
            all_areas_complete_tree_out.write(results["all-areas.complete"])
        if all_areas_extant_tree_out is not None:
            all_areas_extant_tree_out.write(results["all-areas.extant"])
        if focal_areas_tree_out is not None:
            focal_areas_tree_out.write(results["focal-areas"])
        if histories_out and self.event_log:
            self.event_log.write_histories(
                    out=histories_out,
                    tree=self.phylogeny,)
        return

#         if focal_areas_tree_out is not None:
#             focal_areas_tree = self.phylogeny.extract_focal_areas_tree()
#             n = len(focal_areas_tree.seed_node._child_nodes)
#             if n < 2:
#                 raise error.InsufficientFocalAreaLineagesSimulationException("Insufficient lineages in focal area: {}".format(n))
#             self.write_focal_areas_tree(
#                     out=focal_areas_tree_out,
#                     tree=focal_areas_tree,
#                     )
#         if all_areas_tree_out is not None:
#             self.write_all_areas_tree(
#                     out=all_areas_tree_out,
#                     tree=self.phylogeny,
#                     )
        # if histories_out is not None and self.event_log is not None:
        #     self.write_histories(out=histories_out,
        #             tree=self.phylogeny)

    def write_histories(self, out, tree):
            if self.is_encode_nodes:
                labelf = lambda x: x.encode_lineage(
                        set_label=False,
                        add_annotation=self.is_annotate_nodes,
                        exclude_supplemental_areas=True)
            else:
                labelf = ArchipelagoSimulator.simple_node_label_function
            self.event_log.write_histories(
                    out=out,
                    tree=self.phylogeny,
                    node_label_fn=labelf)

    def write_focal_areas_tree(self, out, tree):
        if self.is_encode_nodes:
            labelf = lambda x: x.encode_lineage(
                    set_label=False,
                    add_annotation=self.is_annotate_nodes,
                    exclude_supplemental_areas=True)
        else:
            labelf = ArchipelagoSimulator.simple_node_label_function
        tree.write_to_stream(
                out,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_fn=labelf,
                suppress_internal_node_labels=self.is_suppress_internal_node_labels,
                )

    def write_all_areas_tree(self, out, tree):
        if self.is_encode_nodes:
            labelf = lambda x: x.encode_lineage(
                    set_label=False,
                    add_annotation=self.is_annotate_nodes,
                    exclude_supplemental_areas=False)
        else:
            labelf = ArchipelagoSimulator.simple_node_label_function
        tree.write_to_stream(
                out,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_fn=labelf,
                suppress_internal_node_labels=self.is_suppress_internal_node_labels,
                )

    def debug_compose_tree(self, tree):
        labelf = lambda x: self.model.encode_lineage(x,
                set_label=False,
                add_annotation=False,
                exclude_supplemental_areas=False)
        s = tree.as_string(
                "newick",
                node_label_compose_fn=self.model.encode_all_areas_lineage,
                suppress_edge_lengths=True)
        return s.replace("\n", "")


def repeat_run(
        output_prefix,
        nreps,
        model_definition_source,
        model_definition_type="python-dict",
        config_d=None,
        interpolate_missing_model_values=False,
        random_seed=None,
        stderr_logging_level=None,
        file_logging_level=None,
        maximum_num_restarts_per_replicates=100,
        debug_mode=False):
    """
    Executes multiple runs of the Archipelago simulator under identical
    parameters to produce the specified number of replicates, discarding failed
    runs.

    Parameters
    ----------
    nreps : integer
        Number of replicates to produce.
    config_d : dict
        Simulator configuration parameters as keyword-value pairs. To be
        re-used for each replicate.
    model_definition_source : object
        See 'model_definition_type' argument for values this can take.
    model_definition_type : str
        Whether 'model_definition_source' is:

            - 'python-dict' : a Python dictionary defining the model.
            - 'python-dict-str' : a string providing a Python dictionary
                defining the model.
            - 'python-dict-filepath' : a path to a Python file to be evaluated;
                the file should be a valid Python script containing nothing but a
                dictionary defining the model.
            - 'json-filepath': a path to a JSON file containing a dictionary
                defining the model.

    interpolate_missing_model_values : bool
        Allow missing values in model to be populated by default values (inadvisable).
    random_seed : integer
        Random seed to be used (for single random number generator across all
        replicates).
    stderr_logging_level : string or None
        Message level threshold for screen logs; if 'none' or `None`, screen
        logs will be supprsed.
    file_logging_level : string or None
        Message level threshold for file logs; if 'none' or `None`, file
        logs will be supprsed.
    maximum_num_restarts_per_replicates : int
        A failed replicate (due to e.g., total extinction of all taxa) will be
        re-run. This limits the number of re-runs.
    """
    if output_prefix is None:
        output_prefix = config_d.pop("output_prefix", "archipelago")
    if config_d is None:
        config_d = {}
    config_d["output_prefix"] = output_prefix
    if stderr_logging_level is None:
        if debug_mode:
            log_to_stderr = True
            stderr_logging_level = "debug"
        else:
            log_to_stderr = True
            stderr_logging_level = "info"
    elif stderr_logging_level.lower() == "none":
        log_to_stderr = False
    else:
        log_to_stderr = True
    if file_logging_level is None :
        if debug_mode:
            log_to_file = True
            file_logging_level = "debug"
        else:
            log_to_file = True
            file_logging_level = "info"
    elif file_logging_level.lower() == "none":
        log_to_file = False
    else:
        log_to_file = True
    config_d["debug_mode"] = debug_mode
    if "run_logger" not in config_d:
        config_d["run_logger"] = utility.RunLogger(
                name="archipelago",
                log_to_stderr=log_to_stderr,
                stderr_logging_level=stderr_logging_level,
                log_to_file=log_to_file,
                log_path=output_prefix + ".log",
                file_logging_level=file_logging_level,
                )
    run_logger = config_d["run_logger"]
    run_logger.info("-archipelago- Starting: {}".format(archipelago.description()))
    if "rng" not in config_d:
        if random_seed is None:
            random_seed = config_d.pop("random_seed", None)
            if random_seed is None:
                random_seed = random.randint(0, sys.maxsize)
        else:
            random_seed = int(random_seed)
        run_logger.info("-archipelago- Initializing with random seed: {}".format(random_seed))
        config_d["rng"] = random.Random(random_seed)
    else:
        run_logger.info("-archipelago- Using existing RNG: {}".format(config_d["rng"]))

    if config_d.get("store_all_areas_complete_trees", True) and "all_areas_complete_trees_file" not in config_d:
        config_d["all_areas_complete_trees_file"] = open(ArchipelagoSimulator.compose_all_areas_complete_trees_filepath(output_prefix), "w")
    if config_d.get("store_all_areas_extant_trees", True) and "all_areas_extant_trees_file" not in config_d:
        config_d["all_areas_extant_trees_file"] = open(ArchipelagoSimulator.compose_all_areas_extant_trees_filepath(output_prefix), "w")
    if config_d.get("store_focal_areas_trees", True) and "focal_areas_trees_file" not in config_d:
        config_d["focal_areas_trees_file"] = open(ArchipelagoSimulator.compose_focal_areas_trees_filepath(output_prefix), "w")

    if config_d.get("store_histories", False) and "areas_histories_file" not in config_d:
        config_d["histories_file"] = open(ArchipelagoSimulator.compose_histories_filepath(output_prefix), "w")
    if "histories_file" in config_d:
        config_d["histories_file"].write("[\n")

    try:
        current_rep = 0
        while current_rep < nreps:
            if current_rep > 0 and "histories_file" in config_d:
                config_d["histories_file"].write(",\n")
            simulation_name="Run{}".format((current_rep+1))
            run_output_prefix = "{}.R{:04d}".format(output_prefix, current_rep+1)
            run_logger.info("-archipelago- Replicate {} of {}: Starting".format(current_rep+1, nreps))
            num_restarts = 0
            while True:
                if num_restarts == 0 and current_rep == 0:
                    is_verbose_setup = True
                    model_setup_logger = run_logger
                else:
                    is_verbose_setup = False
                    model_setup_logger = None
                archipelago_model = model.ArchipelagoModel.create(
                        model_definition_source=model_definition_source,
                        model_definition_type=model_definition_type,
                        interpolate_missing_model_values=interpolate_missing_model_values,
                        run_logger=model_setup_logger,
                        )
                archipelago_simulator = ArchipelagoSimulator(
                    archipelago_model=archipelago_model,
                    config_d=config_d,
                    is_verbose_setup=is_verbose_setup)
                try:
                    archipelago_simulator.run()
                    run_logger.system = None
                except error.ArchipelagoException as e:
                    run_logger.system = None
                    run_logger.info("-archipelago- Replicate {} of {}: Simulation failure before termination condition at t = {}: {}".format(current_rep+1, nreps, archipelago_simulator.elapsed_time, e))
                    num_restarts += 1
                    if num_restarts > maximum_num_restarts_per_replicates:
                        run_logger.info("-archipelago- Replicate {} of {}: Maximum number of restarts exceeded: aborting".format(current_rep+1, nreps))
                        break
                    else:
                        run_logger.info("-archipelago- Replicate {} of {}: Restarting replicate (number of restarts: {})".format(current_rep+1, nreps, num_restarts))
                else:
                    run_logger.system = None
                    run_logger.info("-archipelago- Replicate {} of {}: Completed to termination condition at t = {}".format(current_rep+1, nreps, archipelago_simulator.elapsed_time))
                    num_restarts = 0
                    break
            if "histories_file" in config_d:
                config_d["histories_file"].flush()
            current_rep += 1
    except KeyboardInterrupt:
        pass
    if "histories_file" in config_d:
        config_d["histories_file"].write("\n]\n")
        config_d["histories_file"].flush()
        config_d["histories_file"].close()
