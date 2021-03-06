#! /usr/bin/env python

import sys
import random
import archipelago

## functions to return rates
def birth_weight1(**kwargs):
    """
    Speciation rate is proportional to the number of areas occupied.
    """
    lineage = kwargs["lineage"]
    num_areas = len(lineage.areas)
    r = float(num_areas) / (num_areas + 1)
    return r

def area_gain_weight1(**kwargs):
    """
    Something makes dispersal out of area s1 very high relative to the others
    if the lineage's trait state is 0.
    """
    lineage = kwargs["lineage"]
    from_area = kwargs["from_area"]
    to_area = kwargs["to_area"]
    if from_area.index == 0 and lineage.trait_state("q1") == 0:
        return 100.0
    else:
        return 1.0

def area_loss_weight1(**kwargs):
    """
    Another very convoluted condition: if the current lineage has trait state
    1, and there are lineages with trait state 0 in the area, then loss weight
    is really high.
    """
    lineage = kwargs["lineage"]
    area = kwargs["area"]
    if lineage.trait_state("q1") == 1:
        for other_lineage in area.lineages:
            if other_lineage.trait_state("q1") == 0:
                return 100.0
    return 1.0

## model definition
model_definition_source = {
    "areas": [
        {'is_supplemental': True, 'label': 's1'},
        {'is_supplemental': False, 'label': 'a1'},
        {'is_supplemental': False, 'label': 'a2'},
        {'is_supplemental': False, 'label': 'a3'},
        {'is_supplemental': False, 'label': 'a4'},
        ],
    "traits" : [
        {"label": "q1", "nstates": 3, "transition_rate": 0.01, },
    ],
    "diversification": {
        "mean_birth_rate" : 1.0,
        "lineage_birth_weight": {
            "definition_type": "function_object",
            "definition": birth_weight1,
        },
        "mean_death_rate" : 0.0,
        "lineage_death_weight": {
            "definition_type": "fixed_value",
            "definition": 1.0
        },
    },
    "anagenetic_range_evolution": {
        "global_area_gain_rate": 0.01,
        "lineage_area_gain_weight": {
            "definition_type": "function_object",
            "definition": area_gain_weight1,
        },
        "mean_area_loss_rate": 1.0,
        "lineage_area_loss_weight": {
            "definition_type": "function_object",
            "definition": area_loss_weight1,
        },
    },
    "cladogenetic_range_evolution": {
        "sympatric_subset_speciation_weight": 1.0,
        "single_area_vicariance_speciation_weight": 1.0,
        "widespread_vicariance_speciation_weight": 1.0,
        "founder_event_speciation_weight": 0.0
    },
    "termination_conditions": {
        "target_focal_area_lineages": 20,
    }
}

# execute simulations
archipelago.run(
        output_prefix="results/m1",
        nreps=10,
        model_definition_source=model_definition_source,
        model_definition_type="python-dict",
        random_seed=random.randint(0, sys.maxsize))

