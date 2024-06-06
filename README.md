# Robust Knowledge-Guided Biclustering for Multi-omics Data (Briefings in Bioinformatics, 2024)
Qiyiwen Zhang, Changgee Chang & Qi Long

## Abstract 
We propose a novel Bayesian biclustering method called Bayesian graph-guided biclustering (BGB). Specifically, we introduce a new hierarchical sparsity-inducing prior to effectively incorporate biological graph information, and establish a unified framework to model multi-view data. We develop an efficient Markov chain Monte Carlo (MCMC) algorithm to conduct posterior sampling and inference.

## Code Description
- BGB_mcmcfunction.R contains the main function for the model BGB.
- analyze_function.R contains the functions for biclustering result evaluation.  
- analyze_res.R contains the sample code to analyze the bicluster result based on MCMC chain.  
- data_generation.R contains the R code to generate the synthetic data.
- sim_gau.R contains the sample code for the model implementation and model selection. 
- gbcmetric.R contains the another types of measurement function for biclustering results. 

## Link
[Original paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10701104/)
