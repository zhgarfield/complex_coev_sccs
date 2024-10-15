# Model Fitting and Comparison 

# NOTE: This script relies on "rstan" and convenience functions from the "rethinking" package. For installation instructions see https://github.com/rmcelreath/rethinking)
if (!require("tidyverse")) install.packages("tidyverse")
# if (!require("phytools")) install.packages("phytools")
if (!require("qgraph")) install.packages("qgraph") # Correlation network plot
if (!require("patchwork")) install.packages("patchwork") # Correlation network plot

library(rstan)
library(rethinking)
library(tidyverse)

# Read in data
d <- read.csv("data_analysis.csv", stringsAsFactors = F)
#sccs_tree <- ape::read.tree("SCCS_supertree.tre")

#sum(d$socname == sccs_tree$tip.label)/length(d$socname) # checking for discrepancies between phylo tree names and dataframe names, should return 1 if all match

#### Data dictionary ####
## id: sccs id number
## socname: society name
## v149: Writing and Records
## v150: Fixity of Residence
## v151: Agriculture
## v152: Urbanization
## v153: Technological Specialization
## v154: Land Transport
## v155: Money
## v156: Density of Population
## v157: Political Integration
## v158: Social Stratification
## v204: Dependence on Hunting
## v20: Food Storage
## v820: Subsistence Economy (primary mode)
## v179: Latitude
## v181: Longitude
## v200: Region

#### Create phylogenetic distance matrix ##
#dist_mat <- cophenetic.phylo(sccs_tree) # pairwise distance matrix using branch lengths
#dist_mat <- dist_mat / max(dist_mat) # scaling matrix to [0,1]

#### Organizing data ##
storage <- ifelse(d$v20 > 1, 1, 0) # dichtomizing food storage
storage <- ifelse(is.na(storage), -99, storage) # flagging missing storage observations
hunting <- d$v204 + 1 # making the lowest value 1
K_hunt <- max(hunting) # number of ordinal levels

J <- 12 # number of variables
K <- max(d$v149) # number of ordinal levels
N <- nrow(d) # num observations

# Reverse coding hunting to be high -> low
hunting_rev <- as.factor(hunting)
levels(hunting_rev) <- rev(levels(hunting_rev))
hunting_rev <- as.numeric( as.character(hunting_rev) )

# Organize data into list for stan, now adding a third factor for Community Size (CS)
data_list <- list(
  y = cbind(d$v151, d$v150, d$v156, d$v152, d$v149, d$v153, d$v154, d$v155, d$v157, d$v158),  # all variables
  hunting = hunting_rev,
  storage = storage,
  N = N,
  J = J,
  K = K,
  K_hunt = K_hunt,
  zeros = matrix(0, N, J) # we'll use these within stan models
)

# Fit stan model with three latent factors
n_iter <- 1000  # number of samples
n_chains <- 10  # number of chains

fit_m3 <- stan(file = "stan_models/m3_three_factors.stan", 
               data = data_list, 
               iter = n_iter, 
               chains = n_chains, 
               cores = n_chains, 
               init = 0,  
               control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_m3, "fit_m3_three_factors.rds")

# Extract the samples
post_samples <- as.array(fit_m3)

# Extract parameters related to res_v using regex matching
res_v_samples <- post_samples[, , grepl("^res_v", dimnames(post_samples)$parameters)]

# Check the dimensions of the extracted res_v_samples
dim(res_v_samples)  # Should reflect three factors now

# Extract parameter names
param_names <- dimnames(post_samples)$parameters
res_v_param_names <- param_names[grepl("^res_v", param_names)]

# Extract the factors
RI_factors_samples <- res_v_samples[, , grepl("res_v\\[.*,1\\]", res_v_param_names)]  # RI factor
TSD_factors_samples <- res_v_samples[, , grepl("res_v\\[.*,2\\]", res_v_param_names)]  # TSD factor
CS_factors_samples <- res_v_samples[, , grepl("res_v\\[.*,3\\]", res_v_param_names)]   # CS factor

# Calculate the posterior mean for RI, TSD, and CS factors across all iterations and chains
RI_factors_mean <- apply(RI_factors_samples, 3, mean)  # Mean RI values for each society
TSD_factors_mean <- apply(TSD_factors_samples, 3, mean)  # Mean TSD values for each society
CS_factors_mean <- apply(CS_factors_samples, 3, mean)  # Mean CS values for each society

# Combine RI, TSD, and CS factor means into a data frame with society names
factor_values <- data.frame(
  society = d$socname,
  sccs_id = d$id,
  RI_factor = RI_factors_mean,
  TSD_factor = TSD_factors_mean,
  CS_factor = CS_factors_mean
)

# Factor Loadings Extraction

lambda_TSD_samples <- post_samples[, , grepl("^lambda_TSD", param_names)]  # TSD factor loadings (6 elements)
lambda_RI_samples <- post_samples[, , grepl("^lambda_RI", param_names)]  # RI factor loadings (6 elements)
lambda_CS_samples <- post_samples[, , grepl("^lambda_CS", param_names)]  # CS factor loadings (2 elements for Urbanization and Density of Population)

# Calculate the posterior mean for each loading
lambda_TSD_mean <- apply(lambda_TSD_samples, 3, mean)  # Mean TSD loadings (6 variables)
lambda_RI_mean <- apply(lambda_RI_samples, 3, mean)  # Mean RI loadings (6 variables)
lambda_CS_mean <- apply(lambda_CS_samples, 3, mean)  # Mean CS loadings (2 variables)

# Update variables for each factor
TSD_variables <- c("Writing and Records (v149)", "Land Transport (v154)", "Technological Specialization (v153)", "Money (v155)", "Political Integration (v157)", "Social Stratification (v158)")
RI_variables <- c("Agriculture (v151)", "Dependence on Hunting (v204)", "Food Storage (v20)", "Density of Population (v156)", "Urbanization (v152)", "Fixity of Residence (v150)")
CS_variables <- c("Urbanization (v152)", "Density of Population (v156)")  # These load on the new CS factor

# Create data frames for each factor
loadings_df_RI <- data.frame(variable = RI_variables, RI_loading = c(1, lambda_RI_mean))  # Fix first factor loading to 1 for RI
loadings_df_TSD <- data.frame(variable = TSD_variables, TSD_loading = lambda_TSD_mean)
loadings_df_CS <- data.frame(variable = CS_variables, CS_loading = lambda_CS_mean)

# Combine into one final data frame
final_loadings_df <- rbind(loadings_df_RI, loadings_df_TSD, loadings_df_CS)


