# fiji-zikv-model: a model of Zika virus transmission fit to surveillance, serological and molecular data from Fiji between 2013 and 2017 

## Summary
This repository contains code and data to accompany analysis of Zika virus transmission in Fiji between 2013 and 2017 and a comparison of flavivirus dynamics in the country. 

## Guide to scripts
Open `fiji-zikv-model.Rproj` 

Optionally adjust objects in `Rscripts/preamble_zikvfiji.R`. Here you can change: 
* the length of the MCMC run
* the number of chains used in fitting 
* to include/exclude serological data in the model fitting  
* to include/exclude cross protection during 2013--14 DENV-3 outbreak 
* length of burn in of posteriors
* give the model run a unique name

## Guide to modelling in order
1. `denv3-fiji-mcmc.R` Model fitting to 2013-14 DENV-3 epidemic and temperature data
2. `beast/mono_2_beast.xml` Model sequence data to extract tMRCA for Fiji
3. `zika-fiji-mcmc.R` Main ZIKV transmission model fitting

`plot-zika-fiji.R` Analysis of model fit. Figures 1-3, table 1. 

## Reference
To be added

## DOI
[![DOI](https://zenodo.org/badge/241614932.svg)](https://zenodo.org/badge/latestdoi/241614932)