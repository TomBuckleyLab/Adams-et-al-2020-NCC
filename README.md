# Adams-et-al-2020
Code and data for Adams et al 2020, Nature Climate Change, "Diminishing CO2-driven gains in water use efficiency of global forests"

https://doi.org/10.5281/zenodo.3693240

This repository contains R code and associated data used to generate output presented in  the paper referenced above. Measurements of stable carbon isotope discrimination in tree rings, or quantities derived therefrom, were collated from the literature and used to calculate intrinsic water use efficiency (W, the ratio of net CO2 assimilation rate to stomatal conductance) and its rate of change with respect to atmospheric CO2 concentration (ca), or dW/dca, over the 20th century. 

isotopes.R: R code to process tree ring isotope data series reported in different forms (e.g., d13C, WUE, ci,  Delta) and calculate rates of change of imputed intrinsic water use efficiency with respect to atmospheric CO2 concentration. 

Two input files are also included: 

input.csv: isotope data and associated ancillary data collated from the literature. Columns are as follows:
 - ID: identifier for specific location or specier
 - year: year of reported measurement
 - y: data value reported by authors (e.g., d13C)
 - study: current authors' index of previously published study
 - datatype: form of "y"; see R code for more explanation
 - tree: identifier unique to each isotope data series (combination of study and ID)
 - site: identifier unique to each site (may not be unique to each "tree")
 - author: shorthand notation for author of original study
 - spp: species name
 - fixer: y/n nitrogen fixing species
 - lifeform: d/e deciduous/evergreen
 - group: a/g angiosperm/gymnosperm
 - lat: latitude in degrees
 - lon: longitude in degrees
 - ca: atmospheric CO2 in ppm (micromoles per mole)
 - d13Catm: delta-13C (isotopic composition relative to standard) of atmospheric CO2 (permille)

processed.csv: input file described above, augmented with calculated outputs using the first part of the code in isotopes.R:
 - Delta: isotopic discrimination
 - A_ca: imputed ratio of net CO2 assimilation rate to ca based on an assumed doubling ratio of A wrt ca
 - ci_ca: imputed ratio of intercellular to ambient CO2 concentrations based on A_ca and Delta
 - iWUE: imputed intrinsic water use efficiency (referred to as W in the paper) based on ca and ci_ca
 - ci_ca_simple: ci/ca calculated ignoring effects of mesophyll conductance and photorespiration
 - iWUE_simple: iWUE calculated from ca and ci_ca_simple

Climate data referred to in the study are not included, as they are from a variety of sources, some of which do not allow redistribution; they can be reconstructed from those sources (see main text of the paper) given the latitudes and longitudes of each study site in input.csv. 
