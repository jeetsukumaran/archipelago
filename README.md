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
-   [DendroPy Phylogenetic Computing Library, Version 4 or above](http://dendropy.org/)
-   [R](http://www.r-project.org/)
-   The following *R* packages:
    -   adegenet
    -   picante
    -   BioGeoBears
    -   GEIGER

### Installing from Source

Run the following in the top-level directory of the project:

    $ python setup.py install

## Documentation

Full documentation is available online here:

   http://jeetsukumaran.github.io/archipelago

or in the ``archipelago.pdf`` document found in the root directory of the
project.
