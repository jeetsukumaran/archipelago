#! /usr/bin/env python

import sys
import random
import archipelago

## functions to return rates
def birth_rate1(lineage):
    """
    Speciation rate is proportional to the number of areas occupied.
    """
    num_areas = sum(lineage.distribution_vector)
    r = float(num_areas) / (len(lineage.distribution_vector) + 1)
    return r

def dispersal_rate1(lineage):
    """
    Something makes dispersal out of area s1 very high relative to the others.
    """
    if lineage.distribution_vector[0] == 1:
        return 100.0
    else:
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
        "lineage_birth_rate": {
            "definition_type": "function_object",
            "definition": birth_rate1,
        },
        "lineage_death_rate": {
            "definition_type": "fixed_value",
            "definition": 0,
        },
    },
    "dispersal": {
        "global_dispersal_rate": 0.01,
        "lineage_dispersal_weight": {
            "definition_type": "function_object",
            "definition": dispersal_rate1,
        },
    },
    "cladogenesis": {
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

