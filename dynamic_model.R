library(rethinking)
library(phangorn)
library(phytools)
library(tidyverse)
library(patchwork)
library(phaseR)
library(deSolve)

# Read in data
d <- read.csv("data_analysis.csv", stringsAsFactors = F)
sccs_tree <- ape::read.tree("SCCS_supertree.tre")

setdiff(sccs_tree$tip.label, d$socname) # checking for discrepancies between phylo tree names and dataframe names, should return "character(0)" if all is well

N <- nrow(d) # number of tips

#### Cut up phylogenetic tree into segments ########
tree <- sccs_tree
plot(tree)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree))

times <- node.depth.edgelength(tree) # date of each node, measured in thousands of years

split_points <- sort(unique(times))
node_time <- match(times, split_points)

node_seq <- seq(from=1, to=length(node_time))
node_seq <- node_seq[order(node_time)]
parent <- Ancestors(tree, node_seq, type="parent") 

parent_time <- rep(NA, length(node_seq))
parent_time[1] <- -99
for (i in 2:length(parent_time)) parent_time[i] <- (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]]) / max(node.depth.edgelength(tree))

N_seg <- length(node_seq)
#######################
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

y <- cbind(d$v151, d$v150, d$v156, d$v152, d$v149, d$v153, d$v154, d$v155, d$v157, d$v158) # all the Murdock complexity variables

#### Center priors for observed variables based on their sample means #############
prior_y <- matrix(NA, nrow = 10, ncol = K-1)

for (j in 1:(J-2)) {
  cumu_logit <- logit( cumsum( table(y[,j]) ) / length(y[,j][y[,j] >= 0]) )
  
  prior_y[j,] <- cumu_logit[-K] # omitting the last category, which is always +Inf
}

prior_storage <- logit( mean(storage[storage >= 0]))

prior_hunting <- (logit( cumsum( table(hunting_rev) ) / length(hunting_rev[hunting_rev >= 0]) ))[-K_hunt]

## Organize data into a list for Stan
data_list <- list(
  N = N,
  J = J,
  N_seg = N_seg,
  node_seq = node_seq,
  parent = parent,
  ts = cbind(0.0001, parent_time),
  date = times[order(node_time)],
  y = y,
  hunting = hunting_rev,
  storage = storage,
  K = max(d$v151),
  K_hunt = K_hunt,
  prior_y = prior_y,
  prior_storage = prior_storage,
  prior_hunting = prior_hunting
)

n_chains <- 8
n_cores <- 8
n_iter <- 1000

mod <- stan_model(file="stan_models/mcoev_OU.stan")
fit_mcoev <- sampling(mod, data=data_list, chains=n_chains, cores=n_cores, iter=n_iter, init="0", control=list(adapt_delta = 0.99))

#saveRDS(fit_mcoev, "fit_mcoev.rds")
#fit_mcoev <- readRDS("fit_mcoev.rds")

post <- extract.samples(fit_mcoev)

###################################
##### Scatter plot of posterior median trait values for sample societies
svg("coev_scatter.svg", width=6, height=6, pointsize=12)
par(cex=1.5, pty="s")

plot( scale(apply(post$z[,1:data_list$N,2,2], 2, median)) ~ scale(apply(post$z[,1:data_list$N,2,1], 2, median)), pch=16, col=col.alpha("black",0.6),xlab="RI (z-score)", ylab="TPC (z-score)")
points( scale(apply(post$z[,1:data_list$N,2,2], 2, median)) ~ scale(apply(post$z[,1:data_list$N,2,1], 2, median)), col="black")

dev.off()

#### Get parameter values and distribution of latent variables among study societes ###
theta_med <- apply(post$theta, 2, median) # posterior median parameter values for OU model

mean_RI <- mean(apply(post$z[,1:data_list$N,2,1], 2, median))
mean_TPC <- mean(apply(post$z[,1:data_list$N,2,2], 2, median))

sd_RI <- median(apply(post$z[,1:data_list$N,2,1], 1, sd))
sd_TPC <- median(apply(post$z[,1:data_list$N,2,2], 1, sd))

low_RI <- mean_RI - sd_RI*2 # -2SD
high_RI <- mean_RI + sd_RI*2 # + 2SD

low_TPC <- mean_TPC - sd_TPC*2 # -2SD
high_TPC <- mean_TPC + sd_TPC*2 # + 2SD

##### Flow field diagram ################################
OU <- function(t, y, parameters) {
  y1 <- y[1]
  y2 <- y[2]
  
  theta <- c()
  for (t in 1:length(parameters)) theta[t] <- parameters[t]
  
  dy <- numeric(2)
  dy[1] <- (theta[13]*(theta[1] + theta[2]*y1 + theta[3]*y1 + theta[4]*y2 + theta[5]*y2^2 + theta[6]*y1*y2 - y1))
  dy[2] <- (theta[14]*(theta[7] + theta[8]*y2 + theta[9]*y2 + theta[10]*y1 + theta[11]*y1^2 + theta[12]*y1*y2 - y2))
  
  list(dy)
}

svg("phaseplane.svg", width=6, height=6, pointsize=12)
par(cex=1.5, pty="s")

## Plot phase plane
OU.flowField <- flowField(OU, xlim = c(low_RI, high_RI), ylim = c(low_TPC, high_TPC), parameters=apply(post$theta, 2, median), add = FALSE, xlab="RI (z-score)", ylab="TPC (z-score)", axes=T, arrow.type="proportional", col="black", bg="black", frac=1.1, xaxt='n', yaxt='n', xaxs="i", yaxs='i', points=10)

## Add nullclines to phase plane
nullclines(OU, xlim =  c(low_RI, high_RI), ylim = c(low_TPC, high_TPC), apply(post$theta, 2, median), xlab="Resource Use Intensification", ylab="Techno-political Complexity", points=30, axes=F, col=c( "skyblue", "slategray"), add.legend=F, lwd=3)

# Add axes
axis(1, at=c(low_RI, mean_RI, high_RI), labels=(c(low_RI, mean_RI,high_RI)-mean_RI)/sd_RI)
axis(2, at=c(low_TPC, mean_TPC, high_TPC), labels=(c(low_TPC, mean_TPC,high_TPC)-mean_TPC)/sd_TPC)

dev.off()

##############################################################
##### Make predictions for difference in trait values ########
OU_ode <- function( time, y, parms ) {
  with(as.list(c(y,parms)), {
    
    dy1 <- parms[13]*(parms[1] + parms[2]*y[1] + parms[3]*y[1]^2 + parms[4]*y[2] + parms[5]*y[2]^2 + parms[6]*y[1]*y[2] - y[1])
    
    dy2 <- parms[14]*(parms[7] + parms[8]*y[2] + parms[9]*y[2]^2 + parms[10]*y[1] + parms[11]*y[1]^2 + parms[12]*y[1]*y[2] - y[2])
    
    list(c(dy1,dy2))
  })
}

n_samps <- length(post$lp__) # num posterior samples

mean_TPC <- median(mean_TPC)
sd_TPC <- median(sd_TPC)

mean_RI <- median(mean_RI)
sd_RI <- median(sd_RI)

# Over what time scale do we want to predict change? 500 years.
time_scale <- 500 / (max(node.depth.edgelength(tree)) * 100) 

times <- c(0,time_scale)

delta_predict <- function( RI = mean_RI, TPC = mean_TPC, resp = NA, parms, samp_id ) {
  
  y <- c(RI, TPC)
  out <- ode(y, times, OU_ode, parms)
  
  if (resp == "TPC") {
    
    std_dev <- sd(post$z[samp_id,1:data_list$N,2,2])
    
    delta_z <- (out[2,3] - out[1,3]) / std_dev
  }
  
  else if (resp == "RI") {
    
    std_dev <- sd(post$z[samp_id,1:data_list$N,2,1])
    
    delta_z <- (out[2,2] - out[1,2]) / std_dev
  }
  
  return(delta_z)
}

## Now, predict optimal trait values of RI and TPC across a range of parameter values
n_preds <- 30

# Sequences of latent trait values to predict along
RI_seq <- seq(from = mean_RI - sd_RI*2, to = mean_RI + sd_RI*2, length.out = n_preds)
TPC_seq <- seq(from = mean_TPC - sd_TPC*2, to = mean_TPC + sd_TPC*2, length.out = n_preds)

# Matrices to hold model predictions
pred_RI_low <- matrix(NA, nrow = n_samps, ncol = n_preds); pred_RI_mid <- pred_RI_low; pred_RI_high <- pred_RI_low

pred_TPC_low <- matrix(NA, nrow = n_samps, ncol = n_preds); pred_TPC_mid <- pred_TPC_low; pred_TPC_high <- pred_TPC_low

for (i in 1:n_samps) {
  
  theta <- post$theta[i,]
  
for (j in 1:n_preds) {
  pred_RI_low[i,j] <- delta_predict(TPC = TPC_seq[j], RI = mean_RI - sd_RI*2, resp = "RI", parms = theta, samp_id = i)
  
  pred_RI_mid[i,j] <- delta_predict(TPC = TPC_seq[j], RI = mean_RI, resp = "RI", parms = theta, samp_id = i)
  
  pred_RI_high[i,j] <- delta_predict(TPC = TPC_seq[j], RI = mean_RI + sd_RI*2, resp = "RI", parms = theta, samp_id = i)
  
  pred_TPC_low[i,j] <- delta_predict(RI = RI_seq[j], TPC = mean_TPC - sd_TPC*2, resp = "TPC", parms = theta, samp_id = i)
  
  pred_TPC_mid[i,j] <- delta_predict(RI = RI_seq[j], TPC = mean_TPC, resp = "TPC", parms = theta, samp_id = i)
  
  pred_TPC_high[i,j] <- delta_predict(RI = RI_seq[j], TPC = mean_TPC + sd_TPC*2, resp = "TPC", parms = theta, samp_id = i)
}
}

# Bring all of the predictions together
pred_all_long <- as.data.frame(
  rbind(pred_RI_low, pred_RI_mid, pred_RI_high, pred_TPC_low, pred_TPC_mid, pred_TPC_high)
) %>% 
  mutate(samp_id = rep(1:n_samps, 6)) %>% 
  pivot_longer(-samp_id) %>% 
  mutate(x = as.numeric(substr(name, 2, nchar(name))),
         y_fixed = rep( rep( c("Low (-2 SD)", "Average (0)", "High (+2 SD)"), each = n_samps*n_preds ), 2 ),
         resp = rep( c("RI", "TPC"), each = n()/2)
  )

pred_all_long$facet_label <- ifelse( pred_all_long$resp == "TPC", paste("TPC =", pred_all_long$y_fixed), paste("RI =", pred_all_long$y_fixed))

# Summarize posterior
pred_summary <- pred_all_long %>% 
  group_by(facet_label, x) %>% 
  summarise( med = median(value), resp = unique(resp),
             lower_90 = PI(value, prob=0.9)[1], upper_90 = PI(value, prob=0.9)[2],
             lower_60 = PI(value, prob=0.6)[1], upper_60 = PI(value, prob=0.6)[2],
             lower_30 = PI(value, prob=0.3)[1], upper_30 = PI(value, prob=0.3)[2] )


# Set facet labels for each response
dTPC <- ggplot(filter(pred_summary, resp == "TPC")) + 
  facet_wrap(~fct_relevel(facet_label, "TPC = Low (-2 SD)", "TPC = Average (0)", "High (+2 SD)")) +
  geom_line(aes(x = x, y = med), lwd=1) +
  geom_ribbon(aes(x = x, ymin=lower_90, ymax=upper_90), alpha=0.15, fill = "slategray", color=NA) +
  geom_ribbon(aes(x = x, ymin=lower_60, ymax=upper_60), alpha=0.15, fill = "slategray", color=NA) +
  geom_ribbon(aes(x = x, ymin=lower_30, ymax=upper_30), alpha=0.15, fill = "slategray", color=NA) +
  geom_hline(yintercept = 0, lty="dashed", col="red", lwd=1) +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)*n_preds, labels = c("-2", "-1", "0", "1","2"), expand=c(0,0)) +
  scale_y_continuous(limits=c(-2,2.5)) +
  theme_bw(base_size = 12) + 
  theme(panel.spacing = unit(2, "lines"), strip.background = element_rect(colour=NA, fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(paste(Delta, "TPC")[" (500 years)"])) + 
  xlab("Resource-Use Intensification (RI) z-score")

dRI <- ggplot(filter(pred_summary, resp == "RI")) + 
  facet_wrap(~fct_relevel(facet_label, "RI = Low (-2 SD)", "RI = Average (0)", "RI = High (+2 SD)")) +
  geom_line(aes(x = x, y = med), lwd=1) +
  geom_ribbon(aes(x = x, ymin=lower_90, ymax=upper_90), alpha=0.15, fill = "skyblue", color=NA) +
  geom_ribbon(aes(x = x, ymin=lower_60, ymax=upper_60), alpha=0.15, fill = "skyblue", color=NA) +
  geom_ribbon(aes(x = x, ymin=lower_30, ymax=upper_30), alpha=0.15, fill = "skyblue", color=NA) +
  geom_hline(yintercept = 0, lty="dashed", col="red", lwd=1) +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)*n_preds, labels = c("-2", "-1", "0", "1","2"), expand=c(0,0)) +
  scale_y_continuous(limits=c(-2,2.5)) +
  theme_bw(base_size = 12) + 
  theme(panel.spacing = unit(2, "lines"), strip.background = element_rect(colour=NA, fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(paste(Delta, "RI")[" (500 years)"])) + 
  xlab("Technopolitical Complexity (TPC) z-score")

## Compose plots together and save ####
dTPC / dRI
ggsave(filename = "fig3.pdf", width = 6, height = 5, dpi=600)

#############################################################
##### EVERTTHING BELOW THIS IS A WIP!! ######################

#### R^2 plots ###########################
r2 <- matrix(NA, nrow = n_samps, ncol = J)

for (j in 1:(J-2)) {
  
  if (j == 1) r2[,j] <- apply(post$z[,1:data_list$N,2,1], 1, var) / (apply(post$z[,1:data_list$N,2,1], 1, var) + pi^2/3)
  
  else if (j <= 4) r2[,j] <- apply(post$z[,1:data_list$N,2,1]*post$lambda[,j-1], 1, var) / (apply(post$z[,1:data_list$N,2,1]*post$lambda[,j-1], 1, var) + pi^2/3)
  
  else r2[,j] <- apply(post$z[,1:data_list$N,2,2]*post$lambda[,j-1], 1, var) / (apply(post$z[,1:data_list$N,2,2]*post$lambda[,j-1], 1, var) + pi^2/3)
}

r2[,J-1] <- apply(post$z[,1:data_list$N,2,1]*post$lambda[,J-2], 1, var) / (apply(post$z[,1:data_list$N,2,1]*post$lambda[,J-2], 1, var) + pi^2/3)
r2[,J] <- apply(post$z[,1:data_list$N,2,1]*post$lambda[,J-1], 1, var) / (apply(post$z[,1:data_list$N,2,1]*post$lambda[,J-1], 1, var) + pi^2/3)

r2 <- as.data.frame(r2)
names(r2) <- c("Agriculture", "Fixity of Residence", "Density of Population", "Urbanization", "Writing & Records", "Tech. Specialization", "Land Transportation", "Money", "Political Integration", "Social Stratification")

r2$samp_id <- 1:n_samps

r2_long <- r2 %>% pivot_longer(-samp_id)
##############################################

##### Counterfactual plots for latent vars ###
n_samps <- length(post$lp__)

# add posterior estimates
post_c <- array( post$c, dim=c(n_samps, 10, 4), dimnames = list(post = 1:n_samps, var = c("Agriculture", "Fixity of Residence", "Density of Population", "Urbanization", "Writing and Records", "Technological Specialization", "Land Transport", "Money", "Political Integration", "Social Stratification"),cutpoint = as.character(1:4)))

# reshape
post_long <- post_c %>% as.tbl_cube(met_name = "est") %>% as_tibble()
post_long$cutpoint <- as.character(post_long$cutpoint)

# add the final cutpoint, with log-odds of +Inf
post_5 <- post_long %>% filter(cutpoint == "1")
post_5$cutpoint <- "5"
post_5$est <- Inf

post_long <- bind_rows(post_long, post_5)

# get latent factor loadings
lambda <- as.data.frame( cbind(rep(1, n_samps), post$lambda[,1:9]) )
colnames(lambda) <- c("Agriculture", "Fixity of Residence", "Density of Population", "Urbanization", "Writing and Records", "Technological Specialization", "Land Transport", "Money", "Political Integration", "Social Stratification")
lambda$post <- 1:n_samps

lambda_long <- lambda %>% pivot_longer(-post, names_to = "var", values_to="lambda")

# Combine cutpoints with lambdas
cut_lam <- left_join(post_long, lambda_long, by=c("post", "var"))

# Create prediction sequence
cut_lam <- cut_lam %>% slice(rep(1:n(), each = 20)) # repeat the data frame 20 times
cut_lam$pred_seq <- rep(seq(from=0, to=5, length.out=20), times=nrow(cut_lam)/20) # populate latent factor values

cut_lam$cum_prob <- inv_logit( cut_lam$est - cut_lam$lambda*cut_lam$pred_seq ) # cumulative probability

cut_lam2 <- cut_lam %>% group_by(post, var, pred_seq) %>% arrange(cutpoint) %>% mutate(prob = cum_prob - lag(cum_prob, n=1, default=0)) # convert from cumulative prob to prob by substracting the previous cutpoint cumulative prob

# Now summarise with median and HPDI
cut_lam_sum <- cut_lam2 %>% group_by(var, cutpoint, pred_seq) %>% summarise(med = median(prob),
                                                                            lower1=HPDI(prob,0.96/4)[1],
                                                                            upper1=HPDI(prob,0.96/4)[2],
                                                                            lower2=HPDI(prob,(0.96/4)*2)[1],
                                                                            upper2=HPDI(prob,(0.96/4)*2)[2],
                                                                            lower3=HPDI(prob,(0.96/4)*3)[1],
                                                                            upper3=HPDI(prob,(0.96/4)*3)[2],
                                                                            lower4=HPDI(prob,0.96)[1],
                                                                            upper4=HPDI(prob,0.96)[2])

cut_lam_sum$group <- ifelse(cut_lam_sum$var %in% c("Agriculture", "Fixity of Residence", "Density of Population", "Urbanization"), "Resource Intensification", "Techno-political Complexity")

svg( file="dynamic_cutpoints.svg", width=9, height=8, pointsize=12 )

ggplot(cut_lam_sum, aes(x=pred_seq, y=med, color=cutpoint, fill=cutpoint)) + 
  facet_wrap(group ~ var, scales="free_x") + 
  geom_line() + 
  geom_ribbon(aes(ymin=lower4, ymax=upper4), alpha=0.1, color=NA) +
  geom_ribbon(aes(ymin=lower3, ymax=upper3), alpha=0.1, color=NA) +
  geom_ribbon(aes(ymin=lower2, ymax=upper2), alpha=0.1, color=NA) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1), alpha=0.1, color=NA) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  theme_bw() + 
  xlab("Latent Factor Score") + 
  ylab("Probability") + 
  theme(strip.background =element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill = "Ordinal Value", color= "Ordinal Value")

dev.off()

########################################
##### Hunting prediction ###############
lam_hunt <- data.frame(
  lam = post$lambda[,10],
  samp = 1:length(post$lp__)
)

c_hunt <- as.data.frame(post$c_hunt)
names(c_hunt) <- 1:ncol(c_hunt)
c_hunt$"10" <- Inf
c_hunt$samp <- 1:nrow(c_hunt)

c_hunt_long <- c_hunt %>% pivot_longer(-samp)
hunt_long <- left_join(c_hunt_long, lam_hunt)

hunt_long$name <- as.numeric(hunt_long$name) - 1

# Create prediction sequence
hunt_long <- hunt_long %>% slice(rep(1:n(), each = 20)) # repeat the data frame 20 times
hunt_long$pred_seq <- rep(seq(from=0, to=5, length.out=20), times=nrow(hunt_long)/20) # populate latent factor values

hunt_long$cum_prob <- inv_logit( hunt_long$value - hunt_long$lam*hunt_long$pred_seq ) # cumulative probability

hunt_long2 <- hunt_long %>% group_by(samp, pred_seq) %>% arrange(name) %>% mutate(prob = cum_prob - lag(cum_prob, n=1, default=0)) # convert from cumulative prob to prob by substracting the previous cutpoint cumulative prob

# Now summarise with median and HPDI
hunt_long_sum <- hunt_long2 %>% group_by(name, pred_seq) %>% summarise(med = median(prob),
                                                                            lower1=HPDI(prob,0.96/4)[1],
                                                                            upper1=HPDI(prob,0.96/4)[2],
                                                                            lower2=HPDI(prob,(0.96/4)*2)[1],
                                                                            upper2=HPDI(prob,(0.96/4)*2)[2],
                                                                            lower3=HPDI(prob,(0.96/4)*3)[1],
                                                                            upper3=HPDI(prob,(0.96/4)*3)[2],
                                                                            lower4=HPDI(prob,0.96)[1],
                                                                            upper4=HPDI(prob,0.96)[2])


hunt_long_sum$hunting <- factor(hunt_long_sum$name, labels=c("86-100%","76-85%","66-75%","56-65%","46-55%","36-45%","26-35%","16-25%","6-15%","0-5%"))
#levels(hunt_long_sum$hunting) <- rev(levels(hunt_long_sum$hunting))

ggplot(hunt_long_sum, aes(x=pred_seq, y=med, color=hunting, fill=hunting)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=lower4, ymax=upper4), alpha=0.1, color=NA) +
  geom_ribbon(aes(ymin=lower3, ymax=upper3), alpha=0.1, color=NA) +
  geom_ribbon(aes(ymin=lower2, ymax=upper2), alpha=0.1, color=NA) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1), alpha=0.1, color=NA) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  theme_bw() + 
  xlab("Resource-Use Intensification") + 
  ylab("Probability") + 
  theme(strip.background =element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(fill = "% Dependence on Hunting", color= "% Dependence on Hunting")

########################################
##### Food storage pred ################
lam_fs <- data.frame(
  lam = post$lambda[,11],
  alpha = post$a_storage,
  samp = 1:length(post$lp__)
)

fs_long <- lam_fs %>% slice(rep(1:n(), each = 20)) # repeat the data frame 20 times
fs_long$pred_seq <- rep(seq(from=0, to=5, length.out=20), times=nrow(fs_long)/20) # populate latent factor values

fs_long$prob <- inv_logit( fs_long$alpha + fs_long$lam*fs_long$pred_seq )

# Now summarise with median and HPDI
fs_long_sum <- fs_long %>% group_by(pred_seq) %>% summarise(med = median(prob),
                                                                       lower1=HPDI(prob,0.96/4)[1],
                                                                       upper1=HPDI(prob,0.96/4)[2],
                                                                       lower2=HPDI(prob,(0.96/4)*2)[1],
                                                                       upper2=HPDI(prob,(0.96/4)*2)[2],
                                                                       lower3=HPDI(prob,(0.96/4)*3)[1],
                                                                       upper3=HPDI(prob,(0.96/4)*3)[2],
                                                                       lower4=HPDI(prob,0.96)[1],
                                                                       upper4=HPDI(prob,0.96)[2])



ggplot(fs_long_sum, aes(x=pred_seq, y=med)) + 
  geom_line(color="#00A08A") + 
  geom_ribbon(aes(ymin=lower4, ymax=upper4), alpha=0.1, color=NA, fill="#00A08A") +
  geom_ribbon(aes(ymin=lower3, ymax=upper3), alpha=0.1, color=NA, fill="#00A08A") +
  geom_ribbon(aes(ymin=lower2, ymax=upper2), alpha=0.1, color=NA, fill="#00A08A") +
  geom_ribbon(aes(ymin=lower1, ymax=upper1), alpha=0.1, color=NA, fill="#00A08A") + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  theme_bw() + 
  xlab("Resource-Use Intensification") + 
  ylab("Pr (Food Storage)") + 
  theme(strip.background =element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##################################################
##### Distribution of latent vars ################

latent_df <- data.frame(
  est = c(apply(post$z[,1:N,2,1], 2, median), apply(post$z[,1:N,2,2], 2, median)),
  factor = rep(c("Resource-Use Intensification", "Technopolitical Complexity"), each=N)
)

ggplot(latent_df, aes(x=est)) + geom_density(aes(fill=factor),color=NA, alpha=0.6) + theme_bw(base_size=15) + xlab("Latent Factor Values") + ylab("") + scale_fill_manual(values=c("skyblue","slategray")) + scale_y_continuous(expand=c(0,0)) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.title = element_blank())


##################################################
##### Alpha/sigma plot ###########################

alpha_sig <- data.frame(
  ratio = c(post$alpha[,1]/post$sigma[,1], post$alpha[,2]/post$sigma[,2]),
  factor = rep(c("Resource-Use Intensification", "Technopolitical Complexity"), each=n_samps)
)

ggplot(alpha_sig, aes(x=ratio)) + geom_density(aes(fill=factor),color=NA) + theme_bw(base_size=15) + xlab(expression(alpha/sigma)) + ylab("") + scale_fill_manual(values=c("skyblue","slategray")) + scale_y_continuous(expand=c(0,0)) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.title = element_blank())


##########################################
#### Generate table for manuscript #######
par_names <- c("$\\theta_{[\\textrm{RI},\\textrm{RI}]}^1$", "$\\theta_{[\\textrm{RI},\\textrm{RI}]}^2$", "$\\theta_{[\\textrm{RI},\\textrm{TPC}]}^1$", "$\\theta_{[RI,TPC]}^2$", "$\\theta_{[\\textrm{RI},\\textrm{RI*TPC}]}$",
               "$\\theta_{[\\textrm{TPC},\\textrm{TPC}]}^1$", "$\\theta_{[\\textrm{TPC},\\textrm{TPC}]}^2$", "$\\theta_{[\\textrm{TPC},\\textrm{RI}]}^1$", "$\\theta_{[\\textrm{TPC},\\textrm{RI}]}^2$", "$\\theta_{[\\textrm{TPC},\\textrm{TPC*RI}]}$",
               "$\\alpha_{[\\textrm{RI}]}$", "$\\alpha_{[\\textrm{TPC}]}$",
               "$\\sigma_{[\\textrm{RI}]}$", "$\\sigma_{[\\textrm{TPC}]}$")

median_values <- c( apply(post$theta, 2, median), median(post$sigma[,1]), median(post$sigma[,2]) )
SE_values <- c( apply(post$theta, 2, sd), sd(post$sigma[,1]), sd(post$sigma[,2]) )
lower <- round( c( apply(post$theta, 2, HPDI, 0.9)[1,], HPDI(post$sigma[,1], 0.9)[1], HPDI(post$sigma[,2], 0.9)[1] ), 2)
upper <- round( c( apply(post$theta, 2, HPDI, 0.9)[2,], HPDI(post$sigma[,1], 0.9)[2], HPDI(post$sigma[,2], 0.9)[2] ), 2)

table_df <- data.frame(
  Parameter = par_names,
  Median = median_values,
  SE = SE_values,
  HPDI = paste0( paste0( paste0("[",lower), paste0("," ,upper) ), "]" ) 
)

names(table_df)[4] <- "HPDI (90%)"

## This codes gets plugged into markdown document
library(xtable)
print(xtable(table_df), sanitize.text.function = identity, include.rownames = F, sanitize.colnames.function = NULL)

