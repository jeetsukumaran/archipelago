#! /usr/bin/env python

class EventLog(object):

    def __init__(self):
        self.lineage_events = {}

    def register_event(self,
            lineage,
            event_type,
            event_subtype,
            area_idx,
            child0_lineage,
            child1_lineage,
            ):
        ev = {
            "lineage": lineage,
            "event_type": event_type,
            "event_subtype": event_subtype,
            "area_idx": area_idx,
            "child0_lineage": child0_lineage,
            "child1_lineage": child1_lineage,
            }
        self.lineage_events[lineage] = ev

    def log_lineage_extinction(self, lineage):
        del self.lineage_events[lineage]
