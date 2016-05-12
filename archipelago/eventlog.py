#! /usr/bin/env python

class EventLog(object):

    def __init__(self):
        self.edge_events = {}

    def log_event(self,
            event_edge,
            event_type,
            event_subtype,
            area_idx,
            child0_edge,
            child1_edge,
            ):
        ev = {
            "event_edge": event_edge,
            "event_type": event_type,
            "event_subtype": event_subtype,
            "area_idx": area_idx,
            "child0_edge": child0_edge,
            "child1_edge": child1_edge,
            }
        self.edge_events[event_edge] = ev

    def log_lineage_extinction(self, event_edge):
        del self.edge_events[event_edge]
