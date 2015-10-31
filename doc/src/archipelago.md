
# ARCHIPELAGO

## Introduction

"Archipelago" is the name of a generative phylogeny-based model that simultaneously incorporates the diversification processes of speciation and extinction nad the biogeographical processes of area gain ("dispersal") and area loss ("extirpation"), with these processes being differentially regulated by ecological or other traits that are themselves co-evolving on the phylogeny.
The theory and background to this model and its usage is described in the following paper:

    [details]


This software project, "archipelago" presents a suite of programs to generate and analyze data under the Archipelago model.
The primary objective of the analysis is to *classify a dataset with respect to the model that generated it*.
This is, thus, a computational biogeographical model selection analysis program that exploits the power and flexibility of the Archipelago model to allow you to ask and answer historical biogeographical questions of a nature and complexity that are not possible under any other approach.
In particular, instead of asking questions about ancestral area patterns, you can ask questions about processes, and how some processes (e.g., ecology) affect other processes (e.g., dispersal or speciation).


## Installation

### Pre-requisites

-   Python 2.7 or higher
-   [DendroPy Phylogenetic Computing Library](http://dendropy.org/) ([Version 4 or above](https://github.com/jeetsukumaran/DendroPy/tree/DendroPy4))
-   [R](http://www.r-project.org/)
-   The following *R* packages are required to calculate the summary statistics and carry out the model classification:
    -   adegenet
    -   picante
-   The following *R* packages are required by the profiler:
    -   BioGeoBears
    -   GEIGER

### Installing from Source

Run the following in the top-level directory of the project:

    $ python setup.py install

## General Usage


### Workflow

![archipelago analysis workflow](figs/flowchart1.png)
