# Analysis code for Ingham et al. (2020)

This repository contains analysis code for the paper:

Ingham, VA, Elg, S, Nagi, SC & Dondelinger, F, "The transcriptional regulation of sub-lethal insecticide exposure in the major malaria vector Anopheles coluzzii", In Preparation, 2020.

## Network Inference

The network inference code can be run using the main function:

    source('NetworkInference/apply_EDISON.R')
 
 Note that we are using some features from the development version of the EDISON package. This version can be installed as follows:
 
     devtools::install_github(repo='FrankD/EDISON/Package/EDISON/', ref='specify-predictors')
     
All other packages used can be installed directly from CRAN.

## Simulation Study for Validation Experiment

For the simulation study that informed the design of the validation experiment, the entire analysis can be reproduced by running the following RMarkdown file:

    knitr::knit('SimulationStudy/simulation_test_reliability.Rmd')
