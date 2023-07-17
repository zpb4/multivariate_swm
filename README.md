# Multivariate-SWM
Code repository to support WRR manuscript 'A multisite Stochastic Watershed Modeling (SWM) approach with intermittency for regional low flow and flood risk analysis'   
Submitted XXX

## Description
The code below supports the data processing, model fitting, generation, and plotting routines to support the aforementioned manuscript.
## Getting started
### Dependencies
Raw data to support this code in in the following Zenodo repository: https://doi.org/10.5281/zenodo.8155751
### Installing
Requires following R packages:
* rmGarch
* fGarch
* doParallel

### Executing program
The workflow below is configured to run from the file configuration when the Zenodo repository is unzipped and stored in a repository named 'data'
#### Multivariate SWM fitting
Numbering indicates order in which scripts must be run  
Runtimes (in parentheses at end) are estimated with parallelization where applicable on an HPC resource 

1) mvswm_data-process.R: Processes hydrologic simulation (SAC-SMA) and CDEC full natural flow data from .txt files in data repository (<1 min)
2) mvswm_fit-model_intermittent_correlated.R: Fits and generates from logistic regression intermittency model to raw multisite data (15 min)
3) mvswm_fit-model_corr.R: Fits correlated multivariate error model to multisite SAC-SMA model errors (5 min)
4) mvswm_fit-model_ind.R: Fits independent error models to multisite SAC-SMA model errors (5 min)

#### Multivariate SWM generation

1) mvswm_synthetic-gen_corr.R: Generates specified number of correlated multisite SWM samples; saved in 'out' repository (30 min per 1000 samples)
2) mvswm_synthetic-gen_ind.R: Generates specified number of independent multisite SWM samples; saved in 'out' repository (30 min per 1000 samples)

#### Multivariate SWM postprocessing

1) mvswm_intermittency_post-process.R: Applies intermittency postprocessing to generate multisite SWM samples; saved in 'out' repository (10 min per 1000 samples)
2) obs_sim_gev_design_process.R: Calculates observation and simulation based GEV estimates for design events (5 min)
3) mvswm_syn_gev_design_process.R: Calculates GEV estimates for each SWM sample's design events (20 min per 1000 samples)
4) mvswm_cal_trivar-stats_hi-flow.R: Calculates observed and SWM high flow design statistics (20 min per 1000 samples)
5) mvswm_cal_trivar-stats_lo-flow.R: Calculates observed and SWM low flow design statistics (20 min per 1000 samples)

#### Plotting routines {plot subdirectory}

1) plot/plot_all-site-correlation_fig1.R: Plot Spearman error correlations across all 14 sites
2) plot/plot_fit-distributions_intermittency_corr_fig3_ggplot.R: Plots verification statistics for Figure 3
3) plot/plot_ts_density_bar_fig4_ggplot.R: Plots validation statistics for Figure 4

#### Miscellaneous

#### Contact
Zach Brodeur, zpb4@cornell.edu

#### Plotting routines

- main_plot.R: Main plotting script for manuscript figures
