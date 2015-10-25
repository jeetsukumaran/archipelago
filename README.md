# Archipelago

## Introduction

This is a suite of programs to generate and analyze data under the Archipelago model.
It consists of four programs:

-   ``archipelago-profile-trees.py``
    Calculates summary statistics for trees. Typically, this is used to derive
    the calibration process parameter values for the simulator.

-   ``archipelago-simulate.py``
    Simulates data under the "Archipelago" model. Typically, this is used to
    generate phylogenies to be used to derive the training data sets.

-   ``archipelago-summarize.py``
    Calculates summary statistics for one or more phylogenies. This will be used to derive
    the target data from an empirical phylogeny as well as the training data from simulated phylogenies.

-   ``archipelago-classify.py``
    Given a set of training data (summary statistics) and target data,
    classifies the target data using a Discriminant Analysis of Principal
    Components (DAPC).

## Installation

### Pre-requisites

-   Python 2.7 or higher
-   [DendroPy Phylogenetic Computing Library](http://dendropy.org/) ([Version 4 or above](https://github.com/jeetsukumaran/DendroPy/tree/DendroPy4))
-   [R](http://www.r-project.org/)
-   The following *R* packages:
    -   adegenet
    -   picante
    -   BioGeoBears
    -   GEIGER

### Installing from Source

Run the following in the top-level directory of the project:

    $ python setup.py install

## Usage

### Preparation of Data

The primary input is an ultrametic (and rooted) phylogeny in which each of the
tip lineages are associated with presence/absence data over a set of
geographical areas, as well as one or more character traits. This information
currently has to be encoded into the labels of the tips. Future versions of
archipelago will provide programs that encode this information using data
supplied by the user. Currently, however, users will have to manually encode
this information themselves.

The tip label for each lineage consists of three components separated by a
carat ('``^``') character:

    s68^1.2.0^0011

The first component ("s68" in the above example) is an arbitrary species index
that uniquely identifies each lineage.

The second component ("1.2.0" in the above example) is the vector of trait
states for the species. The states of each trait for the lineage is represented
by an integer index separated by a period. In this example, there are three
trait types, with the value of the first trait type "1", the second "2", and
the third "0". Note that the first trait is indexed by 0, not 1. That is, if
there are ``n`` states, then the set of valid trait state indexes is {0, 1, 2,
..., n-1}. All lineages in the system must have the same number of trait types
(three, in this example), and the trait state values are assumed to be given in
sequential order.

The third component ("0011" in the above example) is the geographical
presence/absence vector for the lineage, where "0" indicates the absence of a
lineage from the area, and "1" indicates the presence. In this example, the
lineage is absent from the first and second areas, but present in the third and
fourth. All lineages must have their presences or absences specified in all
areas.

Consider a system given by the following phylogeny:

    [&R] ((A:1,B:1):3,(C:3,(D:2,E:2):1):1);

where the traits are represented by the following character matrix:

    A   001
    B   011
    C   201
    D   100
    E   211

and the distribution over six areas is represented by the following incidence matrix:

    A   110001
    B   010011
    C   101101
    D   001100
    E   011111

We would then encode this information in the phylogeny as follows:

    [&R] ((A^0.0.1^110001:1,B^0.1.1^010011:1):3,(C^2.0.1^101101:3,(D^1.0.0^001100:2,E^2.1.1^011111:2):1):1);

### Calculation of Calibration Parameters

At the very least, we need to provide the simulation program,
"``archipelago-simulate.py``" the following parameters:

    (1) the birth rate
    (2) the extinction rate
    (3) the trait transition rate (for each of the trait types)
    (4) the global dispersal rate

This is done conveniently for us by the program: "``archipelago-profile-trees.py``".
Full help on running the program is available by typing:

    $ archipelago-profile-trees.py --help

This program takes as its input one or more trees, in either Newick, Nexus or
NeXML format, with trait and geographical data encoded in the labels as
described above.
Options are available to specify the parameters that will be estimated.
The output of the program will be a CSV (comma-separated value) file, with a
single row for each tree passed in as output, and the columns the data fields.

### Simulation of Training Data

The simulation program is "``archipelago-simulate.py``".
Full help is available by typing:

    $ archipelago-simulate.py --help

This program requires a *model* file, which is a Python script containing a Python dictionary.
For e.g.:

    {
        "areas": [
            {'is_supplemental': True, 'label': 's1'},
            {'is_supplemental': False, 'label': 'a1'},
            {'is_supplemental': False, 'label': 'a2'},
            {'is_supplemental': False, 'label': 'a3'},
            {'is_supplemental': False, 'label': 'a4'},
            {'is_supplemental': False, 'label': 'a5'},
            {'is_supplemental': False, 'label': 'a6'},
            ],
        "traits" : [
            {"label": "q1", "nstates": 3, "transition_rate": 0.01, },
            {"label": "q2", "nstates": 2, "transition_rate": 0.01, },
            {"label": "q3", "nstates": 2, "transition_rate": 0.01, },
        ],
        "diversification": {
            "lineage_birth_rate": {
                "definition_type": "lambda_definition",
                "definition": "lambda lineage: 0.01",
                "description": "fixed: 0.01"
            }
            "lineage_death_rate": {
                "definition_type": "lambda_definition",
                "definition": "lambda lineage: 0.00",
                "description": "fixed: 0.00"
            }
        },
        "anagenetic_range_evolution": {
            "global_area_gain_rate": 0.01,
            "lineage_area_gain_weight": {
                "definition_type": "lambda_definition",
                "definition": "lambda lineage: 0.01",
                "description": "fixed: 0.01"
            },
            "lineage_area_loss_rate": {
                "definition_type": "lambda_definition",
                "definition": "lambda lineage: 0.0",
                "description": "fixed: 0.0"
            }
        },
        "cladogenetic_range_evolution": {
            "sympatric_subset_speciation_weight": 1.0,
            "single_area_vicariance_speciation_weight": 1.0,
            "widespread_vicariance_speciation_weight": 1.0,
            "founder_event_speciation_weight": 0.0
        },
        "termination_conditions": {
            "target_focal_area_lineages": 50,
        }
    }

The above example sets up a geography consisting of 6 focal areas, "a1" through
"a6", and one supplemental area, "s1".
As the connection weights for the areas are not specified (i.e. giving the relative weight of dispersal from area "a1" to area "a2" for example), then dispersal between all areas are assumed to be equal by default.
Similarly, as the relative diversity weights of the areas are not specified,
then the areas are all assumed to have equal diversity by default.

Three traits are defined: "q1", with 3 states; "q2", with 2 states; and "q3", with 2 states.
The trait evolution rate is 0.01 for all three trait types.
As the optional trait evolution transition weight matrix is not specified, an
equal-rate model is assumed by default

The diversification regime gives constant/fixed birth rate of 0.01.

The global area gain rate (= dispersal rate in the DEC model) is set to 0.01.
Here, the lineage area gain weight is also set to a fixed-value function that
always returns 0.01.
If we wanted to model a trait-dependent lineage-specific area gain weight, where the weight of the dispersal depended on the trait state of trait type "q1", we might do something like:

        "anagenetic_range_evolution": {
            "global_area_gain_rate": 0.01,
            "lineage_area_gain_weight": {
                "definition_type": "trait_state_index_map:q1",
                "definition": [
                    2.0,
                    1.0,
                    0.0
                ]
            },
        },

Here, if the state of a lineage's "q1" trait was "0", the lineage-specific dispersal weight would be 2.0.
If the state of a lineage's "q1" trait was "1", the lineage-specific dispersal weight would be 1.0.
And if the state of a lineage's "q1" trait was "2", the lineage-specific dispersal weight would be 0.0.

The rate of area loss of a lineage (= "extinction" in the DEC model) is fixed to a constant value of 0.0.

Finally, the example specifies a termination condition of 50 (extant or tip) lineages occurring in the focal areas.

The simulation produces as output a single ultrametric phylogeny for each replicate, with
the trait and geographical data encoded in the labels.

## Calculation of Summary Statistics

The program "``archipelago-summarize.py``" calculates the suite of summary
statistics used in the reported studies. It takes as its input one or more
trees with the trait and geographical information encoded in the tip labels as
described above.
The output of the program is a table of summary statistics in CSV format, with
one row per tree in input data, and in the same order.
If the input data is the (encoded) target phylogeny, then the output will be the *target dataset*.
If the input data are the results of the simulations, then the output will be the *training dataset*.
The training data set needs to have a column with the field name
"model.category" identifying the generating model for each row of data.
This can be added by using the "``-l``" or "``--labels``" option of "``archipelago-summarize.py``".
For example, if calculating summary statistics for data generated under a
"constrained" model, "simulation1.focal-areas.trees":

    $ archipelago-summarize.py -l model.category:constrained simulation1.focal-areas.trees

## Classification of Target Data

The program "``archipelago-classify.py``" takes a target dataset and one or
more training data sets as input. The target dataset and training dataset(s)
need to have the same set of summary statistics. The target dataset(s) need to
have a column identifying the generating model for each element of data,
"model.category". The number of principal component axes to retain for the DAPC
function also needs to be specified. The recommended option is to specify
"``--optimize-npca 0``".

The output of the program will be a CSV table showing the posterior
probabilities and model assignments of each of the items in the target dataset.
