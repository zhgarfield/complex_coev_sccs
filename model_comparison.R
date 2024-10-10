# Model Fitting and Comparison 

# NOTE: This script relies on "rstan" and convenience functions from the "rethinking" package. For installation instructions see https://github.com/rmcelreath/rethinking)
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("phytools")) install.packages("phytools")
if (!require("qgraph")) install.packages("qgraph") # Correlation network plot
if (!require("patchwork")) install.packages("patchwork") # Correlation network plot
library(rethinking)

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

# Organize data into list for stan
data_list <- list(
  y = cbind(d$v151, d$v150, d$v156, d$v152, d$v149, d$v153, d$v154, d$v155, d$v157, d$v158), # all the Murdock complexity variables
  hunting = hunting_rev,
  storage = storage,
  N = N,
  J = J,
  K = K,
  K_hunt = K_hunt,
  #dist_mat = dist_mat,
  zeros = matrix(0, N, J) # we'll use these within stan models
)

##### Fit stan models, takes multiple days to fit them all! ####
n_iter <- 500 # how many samples
n_chains <- 10 # how many markov chains

## Uncomment models below to fit once, then save for future use to avoid refitting

######################################################
#### Saturated Model (mS) ############################
#fit_mS <- stan( file="stan_models/mS.stan", data=data_list, iter=n_iter, chains=n_chains, cores=n_chains, init="0", control=list(adapt_delta=0.9) )
#saveRDS(fit_mS, "fit_mS.rds")

#### Model Zero ########################################
#fit_m0 <- stan( file="stan_models/m0.stan", data=data_list, iter=n_iter, chains=n_chains, cores=n_chains, init="0", control=list(adapt_delta=0.9) )
#saveRDS(fit_m0, "fit_m0.rds")

#######################################################
#### Model One ########################################
#fit_m1 <- stan( file="stan_models/m1.stan", data=data_list, iter=n_iter, chains=n_chains, cores=n_chains, init="0", control=list(adapt_delta=0.9) )
#saveRDS(fit_m1, "fit_m1.rds")

#### Model Two ########################################
#fit_m2 <- stan( file="stan_models/m2.stan", data=data_list, iter=n_iter, chains=n_chains, cores=n_chains, init="0", control=list(adapt_delta=0.9) )
#saveRDS(fit_m2, "fit_m2.rds")

#######################################################
#### Model Three ######################################
fit_m3 <- stan(file = "stan_models/m3.stan", data = data_list, iter = n_iter, chains = n_chains, cores = n_chains, init = 0, control = list(adapt_delta = 0.9))
saveRDS(fit_m3, "fit_m3.rds")

#### Model comparison with PSIS LOOCV #################

## Bring in previously fit model objects and estimate PSIS-LOO
# # fit_m0 <- readRDS("fit_m0.rds"); loo_m0 <- loo(fit_m0)
# # fit_m1 <- readRDS("fit_m1.rds"); loo_m1 <- loo(fit_m1)
# # fit_m2 <- readRDS("fit_m2.rds"); loo_m2 <- loo(fit_m2)
# # fit_m3 <- readRDS("fit_m3.rds"); loo_m3 <- loo(fit_m3)
# # fit_mS <- readRDS("fit_mS.rds"); loo_mS <- loo(fit_mS)
# 
# ## Difference in elpd from each model to the saturated model
# elpd_df <- data.frame(model = paste0("M",0:3), elpd_diff = rep(NA,4), se_diff = rep(NA,4))
# models <- c("m0","m1", "m2", "m3")
# 
# for (m in 1:length(models)) {
#   
#   m_compare <- paste0("loo_", models[m])
#   
#   # Difference in elpd, make sure we get the right sign of the difference depending on the comparison
#   elpd_df$elpd_diff[m] <- ifelse( rownames(loo::loo_compare(loo_mS, eval(as.symbol(m_compare))))[2] == "model2", loo::loo_compare(loo_mS, eval(as.symbol(m_compare)))[2,][1], loo::loo_compare(loo_mS, eval(as.symbol(m_compare)))[2,][1] * -1 )
#   
#   # SE of difference
#   elpd_df$se_diff[m] <- loo::loo_compare(loo_mS, eval(as.symbol(m_compare)))[2,][2]
# }
# 
# 
# ## Difference in elpd (M2, M3) - MS
# elpd_df2 <- data.frame(model = c("m2"), elpd_diff = NA, se_diff = NA)
# 
# # Difference in elpd, make sure we get the right sign of the difference depending on the comparison
# elpd_df2$elpd_diff[1] <- -1*loo::loo_compare(loo_m3, loo_m2)[2,][1]
# 
# # SE of difference
# elpd_df2$se_diff[1] <- loo::loo_compare(loo_m3, loo_m2)[2,][2]
# 
# ## Plot both difference plots  
# svg("elpd_diff.svg", width=10, height=7, pointsize=12)
# 
# diff1 <- ggplot(elpd_df, aes(x=as.numeric(elpd_diff), y=model)) + geom_point(size=4) + geom_errorbarh(aes(xmin=elpd_diff - se_diff*2, xmax=elpd_diff + se_diff*2), height=0, lwd=2) + geom_vline(aes(xintercept=0), linetype="dashed") + ylab("") + xlab(expression(paste(Delta, "(ELPD - MS ELPD)"))) + ggtitle("Model Performance") + theme_bw(base_size=24) + theme(panel.grid.minor = element_blank())
# 
# diff2 <- ggplot(elpd_df2, aes(x=as.numeric(elpd_diff), y=model)) + geom_point(size=4) + geom_errorbarh(aes(xmin=elpd_diff - se_diff*2, xmax=elpd_diff + se_diff*2), height=0, lwd=2) + geom_vline(aes(xintercept=0), linetype="dashed") + ylab("") + xlab(expression(paste(Delta, "(M3 ELPD - M2 ELPD)"))) + ggtitle("") + theme_bw(base_size=24) + theme(panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(limits=c(-15,120))
# 
# diff1 + diff2
# 
# dev.off()

## Compare candidate models using different weighting methods
#stacking_weights <- loo::loo_model_weights(list(loo_m1, loo_m2, loo_m3, loo_m4), method="stacking")
#PBMA_weights <- loo::loo_model_weights(list(loo_m1, loo_m2, loo_m3, loo_m4), method="pseudobma")
#WAIC_weights <- rethinking::compare(fit_m1, fit_m2, fit_m3, fit_m4)

#####################################################
#### Phylogentic correlation network ################
# fit_mS <- readRDS("fit_mS.rds")
# post <- extract.samples(fit_mS)
# n_samps <- length(post$lp__)
# 
# cor_phy <- apply(post$Rho_phy, 2:3, median)
# rownames(cor_phy) <- c('Wrt','FoR','Agr','Urb','L.S.','Lnt','Mny','DoP','PlI','ScS','Hnt','FdS')
# colnames(cor_phy) <- c('Wrt','FoR','Agr','Urb','L.S.','Lnt','Mny','DoP','PlI','ScS','Hnt','FdS')
# 
# svg("cor_network.svg", width=7, height=4, pointsize=12)
# 
# qgraph(cor_phy,
#        graph="cor",
#        layout="spring",
#        color = c("#87ceeb", "slategray"),
#        posCol = "#70a494",
#        negCol = "#de8a5a",
#        groups = c("TSD",rep("RI",3), rep("TSD",3), "RI", "TSD", "TSD", "RI", "RI"),
#        vsize=9,
#        borders = T,
#        minimum = 0.1,
#        fade = T )
# 
# dev.off()
# 
# ###################################################
# #### Phylogenetic signal in m3 ####################
# post <- extract.samples(fit_m3)
# 
# var_res_RI <- apply(post$res_v[,,1], 1 ,var)
# var_phy_RI <- apply(post$phy_v[,,1], 1 ,var)
# 
# var_res_TSD <- apply(post$res_v[,,2], 1 ,var)
# var_phy_TSD <- apply(post$phy_v[,,2], 1 ,var)
# 
# r2_phy_RI <- var_phy_RI / (var_phy_RI + var_res_RI)
# r2_phy_TSD <- var_phy_TSD / (var_phy_TSD + var_res_TSD)
# 
# round(median(r2_phy_RI),2)
# round(HPDI(r2_phy_RI, prob=0.9),2)
# 
# round(median(r2_phy_TSD),2)
# round(HPDI(r2_phy_TSD, prob=0.9),2)

# ADDITIONS TO THE ORIGINAL SCRIPT

# Extract the samples
# Convert the stanfit object to an array
post_samples <- as.array(fit_m3)

# Extract parameters related to res_v using regex matching
res_v_samples <- post_samples[, , grepl("^res_v", dimnames(post_samples)$parameters)]

# Check the dimensions of the extracted res_v_samples
dim(res_v_samples)  # Check the dimensions

# Extract parameter names
param_names <- dimnames(post_samples)$parameters
res_v_param_names <- param_names[grepl("^res_v", param_names)]

# Extract the RI factors (where second index = 1)
RI_factors_samples <- res_v_samples[, , grepl("res_v\\[.*,1\\]", res_v_param_names)]

# Extract the TSD factors (where second index = 2)
TSD_factors_samples <- res_v_samples[, , grepl("res_v\\[.*,2\\]", res_v_param_names)]

# Calculate the posterior mean for RI and TSD factors across all iterations and chains
RI_factors_mean <- apply(RI_factors_samples, 3, mean)  # Mean RI values for each society
TSD_factors_mean <- apply(TSD_factors_samples, 3, mean)  # Mean TSD values for each society

# Combine RI and TSD factor means into a data frame with society names
factor_values <- data.frame(
  society = d$socname,
  sccs_id = d$id,
  RI_factor = RI_factors_mean,
  TSD_factor = TSD_factors_mean
)




# Extract the factor loadings for RI and TSD using grepl
lambda_TSD_samples <- post_samples[, , grepl("^lambda_TSD", param_names)]  # TSD factor (6 elements)
lambda_RI_samples <- post_samples[, , grepl("^lambda_RI", param_names)]  # RI factor (5 elements + 1 fixed to 1)

# Calculate the posterior mean for each loading
lambda_TSD_mean <- apply(lambda_TSD_samples, 3, mean)  # Mean TSD loadings (6 variables)
lambda_RI_mean <- apply(lambda_RI_samples, 3, mean)  # Mean RI loadings (6 variables)

# Correct the assignments for each factor:
TSD_variables <- c("Writing and Records (v149)", "Land Transport (v154)", 
                   "Technological Specialization (v153)", "Money (v155)", 
                   "Political Integration (v157)", "Social Stratification (v158)")

RI_variables <- c("Agriculture (v151)", "Dependence on Hunting (v204)", 
                  "Food Storage (v20)", "Density of Population (v156)", 
                  "Urbanization (v152)", "Fixity of Residence (v150)")

# Create the correct loadings data frame for RI (6 variables total, with Agriculture fixed to 1)
loadings_df_RI <- data.frame(
  variable = RI_variables,
  RI_loading = c(1, lambda_RI_mean)  # Fixed loading of 1 for Agriculture (v151)
)

# Create the correct loadings data frame for TSD (6 variables total)
loadings_df_TSD <- data.frame(
  variable = TSD_variables,
  TSD_loading = lambda_TSD_mean  # All TSD loadings
)

# Combine the two data frames into one
final_loadings_df <- rbind(loadings_df_RI, loadings_df_TSD)

