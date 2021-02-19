# complex_coev_sccs

This repo contains data and code to reproduce the analyses and figures in Ringen, Martin, and Jaeggi (in prep) *Novel phylogenetic methods reveal that resource-use intensification
drives the evolution of “complex” societies.*

processing.R creates the dataset and phylogeny used in all analyses. model_comparison.R contains code for fitting all of the static models M0-MS, presented in part 1 of the paper. dynamic_model.R contains code for fitting the dynamical ornstein-uhlenbeck model of co-evolution presented in part 2 of the paper. Stan code supporting these scripts is located in the "stan_models" folder.

