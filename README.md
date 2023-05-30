# markov_model

Code and data for running the analysis from the paper: Transient increased risk of influenza infection following RSV infection in South Africa; findings from the PHIRST study, South Africa, 2016-2018.

Files are: 
- descriptive_analysis.R: This loads and combines the datasets, and runs a descriptive analysis 
- crude_analyses.R: This runs the statistical/crude analyses
- prep_data.R: This prepares the data for use in the markov model
- model_run_covariates.R: This runs the markov modelling
- evaluate_output.R: This evaluates the output from the markov modelling and creates output plots

A summarised version of the data is available in Data/rsv_flu_data.csv

Analysis was conducted using R version 4.2.2
