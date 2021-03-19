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

# Line up date of each node with the split points in the tree
split_points <- sort(unique(times))
node_time <- match(times, split_points)

# Create a sequence of nodes, respecting temporal order
node_seq <- seq(from=1, to=length(node_time))
node_seq <- node_seq[order(node_time)]
parent <- Ancestors(tree, node_seq, type="parent") 

# Parent time indicates amount of time since the parent node, scaled by the total depth of the tree
parent_time <- rep(NA, length(node_seq))
parent_time[1] <- -99
for (i in 2:length(parent_time)) parent_time[i] <- (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]]) / max(node.depth.edgelength(tree))

N_seg <- length(node_seq)  # total num segments in the tree\

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

## Organize data into a list for Stan
data_list <- list(
  N = N,
  J = J,
  N_seg = N_seg,
  node_seq = node_seq,
  parent = parent,
  ts = parent_time,
  y = y,
  hunting = hunting_rev,
  storage = storage,
  K = max(d$v151),
  K_hunt = K_hunt
)

mod <- stan_model("stan_models/mcoev_SDE.stan")

fit <- sampling(mod, data=data_list, iter=500, init="0", chains=10, cores=10)

#saveRDS(fit_mcoev, "fit_mcoev.rds")
#fit_mcoev <- readRDS("fit_mcoev.rds")

post <- extract.samples(fit)

##### Scatter plot of posterior median trait values for sample societies
#svg("coev_scatter.svg", width=6, height=6, pointsize=12)
#par(cex=1.5, pty="s")

plot( scale(post_eta_median$TPC) ~ scale(post_eta_median$RI), pch=16, col=col.alpha("black",0.6),xlab="RI (z-score)", ylab="TPC (z-score)")
points( scale(post_eta_median$TPC) ~ scale(post_eta_median$RI), col="black")

dev.off()

#### Get parameter values and distribution of latent variables among study societes ###
mean_RI <- mean(post_eta_median$RI)
mean_TPC <- mean(post_eta_median$TPC)

sd_RI <- sd(post_eta_median$RI)
sd_TPC <- sd(post_eta_median$TPC)

low_RI <- mean_RI - sd_RI*2 # -2SD
high_RI <- mean_RI + sd_RI*2# + 2SD

low_TPC <- mean_TPC - sd_TPC*2 # -2SD
high_TPC <- mean_TPC + sd_TPC*2 # + 2SD

A <- matrix( apply( cbind(post$`A[1,1]`,post$`A[1,2]`, post$`A[2,1]`,post$`A[2,2]`), 2, median),
                    nrow = 2, ncol = 2, byrow = T) 

b <- c( median(post$`b0[1]`), median(post$`b0[2]`) )

##### Flow field diagram ################################
OU <- function(t, y, parameters) {

  parms <- parameters
  dy <- numeric(2)
  dy[1] <- y[1]*A[1,1] + y[2]*A[2,1] + b[1]
  dy[2] <- y[2]*A[2,2] + y[1]*A[1,2] + b[2]
  
  list(dy)
}

#svg("phaseplane.svg", width=6, height=6, pointsize=12)
#par(cex=1.5, pty="s")

## Plot phase plane
OU.flowField <- flowField(OU, xlim = c(low_RI, high_RI), ylim = c(low_TPC, high_TPC), parameters = NA, add = FALSE, xlab="RI (z-score)", ylab="TPC (z-score)", points=20, col="cornflowerblue")

## Add nullclines to phase plane
nullclines(OU, xlim =  c(low_RI, high_RI), ylim = c(low_TPC, high_TPC), parameters = NA, xlab="Resource Use Intensification", ylab="Techno-political Complexity", points=10, axes=F, col=c( "skyblue", "slategray"), add.legend=F, lwd=3)

# Add axes
axis(1, at=c(low_RI, mean_RI, high_RI), labels=(c(low_RI, mean_RI,high_RI)-mean_RI)/sd_RI)
axis(2, at=c(low_TPC, mean_TPC, high_TPC), labels=(c(low_TPC, mean_TPC,high_TPC)-mean_TPC)/sd_TPC)

dev.off()

##############################################################
##### Make predictions for difference in trait values ########
OU_ode <- function( time, y, parms ) {
  with(as.list(c(y,parms)), {
    
    dy1 <- 
    
    dy2 <- 
    
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
