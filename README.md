# comparability_mobility_data
 Code and data for "Assessing comparability between mobility data sources" manuscript

## Functions code
Functions source code: `functions_sims_aim3.R`
This includes the similation functions, as well as adapted functions for fitting mobility models using informative priors.  

## Source code
Generating mobility matrices: `combining_datasets_2022_02_12.Rmd`
Run first to set up initial conditions for measles simulations: `set_up_ICs_aim3.Rmd`
Rmd with simulations: `simulations_aim3.Rmd`
Make figures and tables using created output: `figures and tables.Rmd`

## Data sets
Data included:  
*Note: individual-level data not provided due to privacy and confidentiality of participants. 
We provide district-level origin-destination matrices for three data sets:  
1) Facebook: `data/M_fb_all.RDS`
2) Mobile phone data: `data/OD_data_diag.rds`
3) Travel survey: `data/mobility_mat_travelsurvey.RDS` 
Note that DHS data are publically available from https://dhsprogram.com/methodology/survey/survey-display-542.cfm   
These origin-destination matrices were used to fit the departure-diffusion models, compare probabilities of departure, and combine in the Bayesian framework.

## Predicted origin-destination matrices using different combinations of departure and diffusion data sets
Each object is a set of 100 predicted origin-destination matrices. These serve as inputs in the measles simulations
1) mobile phone departure, mobile phone diffusion: `data/pred/M_CDR_full.rds`
2) Facebook departure, Facebook diffusion: `data/pred/M_FB_full.rds`
3) Travel survey departure, travel survey diffusion: `data/pred/M_TS_full.rds`
4) DHS departure, travel survey diffusion: `data/pred/M_DHS_TS_full.rds`
5) DHS departure, Facebook diffusion: `data/pred/M_DHS_FB_full.rds`
6) Mobile phone adjusted using DHS departure, mobile phone diffusion: `data/pred/M_CDRwDHS_CDR_full.rds`
7) Facebook adjusted using DHS departure, Facebook diffusion: `data/pred/M_FBwDHS_FB_full.rds`
8) DHS departure, mobile phone diffusion: `data/pred/M_DHS_CDR_full.rds`
9) Pooled raw departure, mobile phone diffusion: `data/pred/M_pool_CDR_full.rds`
10) Pooled weighted departure, mobile phone diffusion: `data/pred/M_pool_wt_CDR_full.rds`
11) Pooled raw departure, Facebook diffusion: `data/pred/M_pool_FB_full.rds`
12) Pooled raw departure, travel survey diffusion: `data/pred/M_pool_TS_full.rds`
13) Facebook departure, Facebook with travel survey prior diffusion: `data/pred/M_FBwTSpr_FB_full.rds`