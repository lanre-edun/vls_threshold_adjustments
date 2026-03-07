[![DOI](https://zenodo.org/badge/1170620682.svg)](https://doi.org/10.5281/zenodo.18905261)
# Title
Validating HIV viral suppression threshold adjustments for comparable estimates using data from nationally representative household surveys in sub-Saharan Africa

# Description
Using nationally representative Population-based HIV Impact Assessment surveys (PHIAs) from sub-Saharan Africa this analyses aimed to (1) evaluate the models estimated by Johnson et al. used for adjusting VLS estimates reported at alternative thresholds (<50, <200, or <400 copies/mL) to the UNAIDS recommended threshold (<1000 copies/mL), in African populations; (2) compare alternative statistical models for describing observed VL distributions; and (3) assess sex and age differences in VL distribution parameters among PLHIV on ART in sub-Saharan Africa. 

# Setup 
Analysis were conducted in R version 4.3.1 and packages required for analysis are specified in individual R scripts

# Data 
The PHIA data used in this analysis are available upon request from the PHIA data portal (https://phia-data.icap.columbia.edu/datasets). The Botswana 2021 survey data can be accessed via the Statistics Botswana microdata portal (https://microdata.statsbots.org.bw:4443/index.php/catalog/26), and the Nigeria 2018 survey data via the National Bureau of Statistics microdata portal (https://microdata.nigerianstat.gov.ng/index.php/catalog/65). 

UNAIDS Estimates 2024 Spectrum files were sources from (https://hivtools.unaids.org/spectrum-file-request/)  

# Usage and structure
1. Data preparation: data preparation code are in data_preparation.R and include code to read in surveys, functions required to clean data and define core variables
2. Data exploration: data exploration code and code for plotting core figure are in exploration and figures script
3. Compute survey weighted viral load suppression estimates for all surveys using survey_vls_estimates and by age and sex using survey_vls_estimates_age_sex
4. Model fitting using brms: fit Weibull, reverse Weibull, Frechet, gamma and lognormal models to PHIA survey data to obtain pooled shape estimates (brms_regression) and sex and age specific shape estimates (brms_regression_sex)
5. Compute adjusted viral load suppression estimate at recommended threshold (<1000) from <50, <200 or <400 copies/mL using the explore distributions and parameters
Using: 
  - adjusted_new_params: using parameters from Johnson et al and PHIA calibrated paramaters
  - adjusted_age_phia_params: using age-specific parameters from PHIA calibration
  - adjusted_sex_phia_params: using sex-specific parameters from PHIA calibration
  - adjusted_gamma_frechet: using PHIA calibrated parameters for the Frechet, gamma and lognormal models
  - adjusted_age_gamma_frechet: using age-specific parameters from PHIA calibration
  - adjusted_sex_gamm_frechet: using sex-specific parameters from PHIA calibration
  - adjusted_pareto: using plausible paramaters for the Pareto distribution
6. Pareto simulations: using plausible parameters for the Pareto distribution we explore the lower bound of the location parameter required to achieve different levels of VLS
7. Sensitivity analyses: repeat analyses using subset of people with HIV experienced on ART (>12 months) and from earlier PHIAs. 

# Plots
Code for scatter plots comparing adjusted vs. observed VLS estimates are in scripts for step 5, while code for other plots and tables are in exploration and figures script.

# License
  Creative Commons Attribution’ (CC BY)
