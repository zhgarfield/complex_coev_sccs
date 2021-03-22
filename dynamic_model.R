library(rethinking)
library(phangorn)
library(phytools)
library(tidyverse)
library(patchwork)
library(colorspace)
library(viridis)
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

y <- cbind(d$v151, d$v150, d$v156, d$v152, d$v149, d$v153, d$v154, d$v155, d$v157, d$v158) # all the Murdock complexity variables: 151 = agriculture, 150 = fixity of residence, 156 = density of pop, 152 = urbanization, 149 = writing, 153 = technological/labor specialization, 154 = land transport technology, 155 = money, 157 = political integration, 158 = social stratification

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

fit <- sampling(mod, data=data_list, iter=2000, init="0", chains=10, cores=10)

#saveRDS(fit, "fit_mcoev.rds")
#fit <- readRDS("fit_mcoev.rds")

post <- extract.samples(fit)
n_samps <- length(post$lp__)

#### Extract median eta values for sample societies
RI_median <- apply(post$eta[,1:N,1], 2, median)
TSD_median <- apply(post$eta[,1:N,2], 2, median)

##### Scatter plot of posterior median trait values for sample societies
svg("coev_scatter.svg", width=6, height=6, pointsize=12)
par(cex=1.5, pty="s")

plot( scale(TSD_median) ~ scale(RI_median), pch=16, col=col.alpha("black",0.6),xlab="RI (z-score)", ylab="TSD (z-score)")
points( scale(TSD_median) ~ scale(RI_median), col="black")

dev.off()

#### Get parameter values and distribution of latent variables among study societes ###
median_RI <- median(RI_median)
median_TSD <- median(TSD_median)

mad_RI <- mad(RI_median)
mad_TSD <- mad(TSD_median)

low_RI <- median_RI - mad_RI*2 # -2SD
high_RI <- median_RI + mad_RI*2# + 2SD

low_TSD <- median_TSD - mad_TSD*2 # -2SD
high_TSD <- median_TSD + mad_TSD*2 # + 2SD

A <- apply(post$A, 2:3, median)
b <- apply(post$b, 2, median)

sigma <- apply(post$sigma, 2, median)

##### Flow field diagram ################################
OU <- function(t, y, parameters) {

  dy <- numeric(2)
  dy[1] <- y[1]*A[1,1] + y[2]*A[1,2] + b[1]
  dy[2] <- y[2]*A[2,2] + y[1]*A[2,1] + b[2]
  
  list(dy)
}

svg("phaseplane.svg", width=6, height=6, pointsize=12)
par(pty="s")

## Plot phase plane
OU.flowField <- flowField(OU, xlim = c(low_RI, high_RI), ylim = c(low_TSD, high_TSD), parameters = NA, add = FALSE, xlab="", ylab="", points=10, col="grey", xaxt='n', yaxt='n', arrow.type = "proportional", frac=1.5, xaxs="i", yaxs="i", axes=F , lwd=2)

mtext(side=1, "RI (z-score)", at=median_RI, line=2.5, cex=1.3)
mtext(side=2, "TSD (z-score)", at=median_TSD, line=2.5, cex=1.3)

## Add nullclines to phase plane
nc <- nullclines(OU, xlim =  c(low_RI, high_RI), ylim = c(low_TSD, high_TSD), parameters = NA, points=20, axes=F, col=c("gray10","#60C4EB"), add.legend=F, lwd=4)

# Add axes
axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)
axis(2, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

dev.off()

##############################################################
##### Delta theta (how do the traits change eq values) #######
delta_theta <- function(RI, TSD, A, b, resp) {
  
  med_RI <- median(RI)
  diff_RI <- mad(RI) 
  
  med_TSD <- median(TSD)
  diff_TSD <- mad(TSD) 
  
  if (resp == "RI") {
  med_theta <- -(med_TSD*A[1,2] + b[1])/A[1,1]
  diff_theta <- -((med_TSD + diff_TSD)*A[1,2] + b[1])/A[1,1]
  
  return( (diff_theta - med_theta) / diff_RI )
  }
  
  else if (resp == "TSD") {
    med_theta <- -(med_RI*A[2,1] + b[2])/A[2,2]
    diff_theta <- -((med_RI + diff_RI)*A[2,1] + b[2])/A[2,2]
    
    return( (diff_theta - med_theta) / diff_TSD )
  }
}

delta_theta_RI <- c()
delta_theta_TSD <- c()

for (i in 1:n_samps) {
  delta_theta_RI[i] <- delta_theta(RI = post$eta[i,1:N,1], TSD = post$eta[i,1:N,2], A = post$A[i,,], b = post$b[i,], resp="RI")
  delta_theta_TSD[i] <- delta_theta(RI = post$eta[i,1:N,1], TSD = post$eta[i,1:N,2], A = post$A[i,,], b = post$b[i,], resp="TSD")
}

delta_theta_df <- data.frame(
  delta_theta = c(delta_theta_RI, delta_theta_TSD),
  resp = rep(c("RI", "TSD"), each = n_samps)
)

round(sum(delta_theta_RI > 0)/n_samps, 2)
round(sum(delta_theta_TSD > 0)/n_samps, 2)


p_delta_theta <- ggplot(delta_theta_df, aes(x=delta_theta)) + 
  geom_density(aes(fill=resp, color=resp), alpha=0.8) + 
  geom_vline(xintercept = 0, color="indianred", linetype='dashed', lwd=1) +
  annotate("text",label=expression("RI" %->% "TSD"), x=1.15, y=1.5, color="deepskyblue4", size=3.5) +
  annotate("text",label=expression(paste("    PP"[">0"]," > 0.99")), x=1.15, y=1.35, color="deepskyblue4", size=3.5) +
  annotate("text",label=expression("TSD" %->% "RI"), x=0.4, y=2.2, color="gray10", size=3.5) +
  annotate("text",label=expression(paste("    PP"[">0"]," = 0.77")), x=0.4, y=2.05, color="gray10", size=3.5) +
  theme_minimal() + 
  scale_x_continuous(limits=c(-1.2,2)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,2.6)) +
  scale_color_manual(values=c("gray10", "#60C4EB")) +
  scale_fill_manual(values=c("gray10", "#60C4EB")) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position = "none") +
  xlab(expression(paste(Delta, theta)["z"])) +
  ylab("")

ggsave(filename = "fig3.pdf", plot=p_delta_theta, width = 6, height = 3.5, dpi=600)

##############################################################
##### Make predictions for difference in trait values ########
OU_sde <- function( RI, TSD, resp ) {
  
  if (resp == "RI") {
  delta_sigma <- ((A[1,1]*RI + A[1,2]*TSD + b[1])/mad_RI) / (sigma[1]/mad_RI)
  }
  
  if (resp == "TSD") {
    delta_sigma <- ((A[2,2]*TSD + A[2,1]*RI + b[2])/mad_TSD) / (sigma[2]/mad_TSD)
  }
 
  return(delta_sigma)
}

RI_seq <- seq(from=low_RI, to=high_RI, length.out=20)
TSD_seq <- seq(from=low_TSD, to=high_TSD, length.out=20)

preds <- expand.grid(x=RI_seq, y=TSD_seq)
preds$RI <- NA # placeholder
preds$TSD <- NA

for (i in 1:nrow(preds)) {
  preds$RI[i] <- OU_sde( RI = preds$x[i], TSD = preds$y[i], resp="RI" )
  preds$TSD[i] <- OU_sde( RI = preds$x[i], TSD = preds$y[i], resp="TSD" )
}

preds_long <- preds %>% pivot_longer(-c(x,y))

preds_long$name <- factor( preds_long$name, labels=c(expression(paste(Delta,"RI") ), expression(paste(Delta,"TSD"))) )


delta_sigma <- ggplot(preds_long, aes(x=(x-median_RI)/mad_RI, y=(y - median_TSD)/mad_TSD,z=value,fill = value)) +
  facet_wrap(~name,labeller = label_parsed) +
  geom_raster() +
  geom_contour(colour = "white", breaks=c(-1,1)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_continuous_divergingx(palette = "Geyser", trans = "reverse",guide = guide_legend(reverse = TRUE) ) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing = unit(1.5, "lines"),strip.background = element_blank()) + 
  labs(fill = expression( frac(paste(Delta,alpha),sigma))) +
  xlab("RI (z-score)") +
  ylab("TSD (z-score)") + 
  coord_fixed()

## Compose plots together and save ####
ggsave(filename = "fig4.pdf", plot=delta_sigma, width = 6, height = 3.5, dpi=600)

#############################################################
#### Correlation of manifest variables with latent vars #####

manifest_cors <- matrix( NA, nrow=n_samps, ncol=J )

for (j in 1:J) {
  
  if (j < 5) eta_var <- apply(post$eta[,1:N,1]*post$lambda[,j], 1, var)
  if (j > 4 & j < J-1) eta_var <- apply(post$eta[,1:N,2]*post$lambda[,j], 1, var)
  if (j > J-2) eta_var <- apply(post$eta[,1:N,1]*post$lambda[,j], 1, var)
  
  manifest_cors[,j] <- sqrt(eta_var / (eta_var + pi^2/3))
}

round(apply(manifest_cors, 2, median),2)
round(apply(manifest_cors, 2, HPDI, prob=0.9),2)


#############################################################
#### Make predictions on original data scales ###############
#############################################################
dev.off()
{ svg(file="glm_preds.svg",width=8.5, height=11)

par(mfrow=c(4,3), pty="s")

n_preds <- 30
RI <- seq(from=low_RI, to=high_RI, length.out=n_preds)
TSD <- seq(from=low_TSD, to=high_RI, length.out=n_preds)

ordinal_cols <- viridis::viridis(n = K)

#151 = agriculture, 150 = fixity of residence, 156 = density of pop, 152 = urbanization, 149 = writing, 153 = technological/labor specialization, 154 = land transport technology, 155 = money, 157 = political integration, 158 = social stratification

####### Agriculture #######
plot(NULL, xlim=c(low_RI, high_RI), ylim=c(0,1), ylab="Pr (Agriculture)", xlab="RI (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Agriculture", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,1,l] - post$lambda[,1]*RI[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,1,l] - post$lambda[,1]*RI[j] ) - inv_logit( post$c[,1,l-1] - post$lambda[,1]*RI[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,1,l-1] - post$lambda[,1]*RI[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)

legend(x=-2, y=1, legend=(c("None", "<10% food supply", ">10%; secondary", "Primary; not intensive","Primary; intensive")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)
  
####### Fixity of Residence #######
plot(NULL, xlim=c(low_RI, high_RI), ylim=c(0,1), ylab="Pr (Fixity of Residence)", xlab="RI (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Fixity of Residence", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,2,l] - post$lambda[,2]*RI[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,2,l] - post$lambda[,2]*RI[j] ) - inv_logit( post$c[,2,l-1] - post$lambda[,2]*RI[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,2,l-1] - post$lambda[,2]*RI[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)

legend(x=-2, y=1, legend=(c("Nomadic", "Seminomadic", "Semisedentary", "Sedentary; impermanent", "Sedentary")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Density of Population #######
plot(NULL, xlim=c(low_RI, high_RI), ylim=c(0,1), ylab="Pr (Density of Pop.)", xlab="RI (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Density of Population", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,3,l] - post$lambda[,3]*RI[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,3,l] - post$lambda[,3]*RI[j] ) - inv_logit( post$c[,3,l-1] - post$lambda[,3]*RI[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,3,l-1] - post$lambda[,3]*RI[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)

legend(x=0, y=1, legend=(c("<1 per sq mi", "1-5 per sq mi", "5-25 per sq mi", "26 - 100 per sq mi", ">100 per sq mi")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Urbanization #######
plot(NULL, xlim=c(low_RI, high_RI), ylim=c(0,1), ylab="Pr (Urbanization)", xlab="RI (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Urbanization", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,4,l] - post$lambda[,4]*RI[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,4,l] - post$lambda[,4]*RI[j] ) - inv_logit( post$c[,4,l-1] - post$lambda[,4]*RI[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,4,l-1] - post$lambda[,4]*RI[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)

legend(x=-2, y=1, legend=(c("<100 persons", "100-200 persons", "200-400 persons", "400 - 1000 persons", ">1000 perons")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Writing #######
plot(NULL, xlim=c(low_TSD, high_TSD), ylim=c(0,1), ylab="Pr (Writing)", xlab="TSD (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Writing", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,5,l] - post$lambda[,5]*TSD[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,5,l] - post$lambda[,5]*TSD[j] ) - inv_logit( post$c[,5,l-1] - post$lambda[,5]*TSD[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,5,l-1] - post$lambda[,5]*TSD[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

legend(x=-1.5, y=1, legend=(c("None", "Mnemonic devices", "Nonwritten records", "True writing; no records", "True writing; records")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Tech/Labor Specialization #######
plot(NULL, xlim=c(low_TSD, high_TSD), ylim=c(0,1), ylab="Pr (Tech Specialization)", xlab="TSD (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Technological/Labor Specialization", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,6,l] - post$lambda[,6]*TSD[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,6,l] - post$lambda[,6]*TSD[j] ) - inv_logit( post$c[,6,l-1] - post$lambda[,6]*TSD[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,6,l-1] - post$lambda[,6]*TSD[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

legend(x=-1.5, y=1, legend=(c("None", "Pottery only", "Loom weaving", "Metalwork", "Smiths, weavers, potters")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)

####### Land Transport Technology #######
plot(NULL, xlim=c(low_TSD, high_TSD), ylim=c(0,1), ylab="Pr (Land Transport)", xlab="TSD (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Land Transport Technology", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,8,l] - post$lambda[,8]*TSD[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,8,l] - post$lambda[,8]*TSD[j] ) - inv_logit( post$c[,8,l-1] - post$lambda[,8]*TSD[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,8,l-1] - post$lambda[,8]*TSD[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

legend(x=-1.2, y=1, legend=(c("Human only", "Pack animals", "Draft animals", "Animal-drawn vehicles", "Automotive vehicles")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Money #######
plot(NULL, xlim=c(low_TSD, high_TSD), ylim=c(0,1), ylab="Pr (Money)", xlab="TSD (z-score)", xaxt='n',xaxs="i", yaxs="i", cex=0.8)
mtext("Money")

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,8,l] - post$lambda[,8]*TSD[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,8,l] - post$lambda[,8]*TSD[j] ) - inv_logit( post$c[,8,l-1] - post$lambda[,8]*TSD[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,8,l-1] - post$lambda[,8]*TSD[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

legend(x=-1.5, y=1, legend=(c("None","Domestiacally usable articles", "Alien currency", "Elementary forms", "True money")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Political Integration #######
plot(NULL, xlim=c(low_TSD, high_TSD), ylim=c(0,1), ylab="Pr (Political Integration)", xlab="TSD (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Political Integration", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,9,l] - post$lambda[,9]*TSD[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,9,l] - post$lambda[,9]*TSD[j] ) - inv_logit( post$c[,9,l-1] - post$lambda[,9]*TSD[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,9,l-1] - post$lambda[,9]*TSD[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

legend(x=-1.5, y=1, legend=(c("None","Autonomous local communities", "1 level", "2 levels","3 levels")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Political Integration #######
plot(NULL, xlim=c(low_TSD, high_TSD), ylim=c(0,1), ylab="Pr (Political Integration)", xlab="TSD (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Political Integration", cex=0.8)

for (l in 1:K) {
  
  pred <- matrix(NA, nrow=n_samps, ncol=n_preds)
  
  for (j in 1:n_preds) {
    
    if (l == 1) pred[,j] = inv_logit( post$c[,10,l] - post$lambda[,10]*TSD[j] ) 
    if (l > 1 & l < K) pred[,j] = inv_logit( post$c[,10,l] - post$lambda[,10]*TSD[j] ) - inv_logit( post$c[,10,l-1] - post$lambda[,10]*TSD[j] )
    if (l == K) pred[,j] = 1 - inv_logit( post$c[,10,l-1] - post$lambda[,10]*TSD[j] )
  }
  
  lines(apply(pred, 2, median), x=RI, col=ordinal_cols[l], lwd=2)
  shade(apply(pred, 2, PI, prob=0.9), RI, col=col.alpha(ordinal_cols[l],0.15))
}

axis(1, at=c(low_TSD, median_TSD, high_TSD), labels=(c(low_TSD, median_TSD,high_TSD)-median_TSD)/mad_TSD)

legend(x=-1.5, y=1, legend=(c("Egalitarian","Hereditary slavery","2 classes, no castes", "2 classes, castes", "3 social classes, castes")), lwd=2, lty="solid", col=ordinal_cols, bty='n', cex=0.8)


####### Food Storage #######
pred_store <- matrix(NA, nrow=n_samps, ncol=n_preds)

for (j in 1:n_preds) {
  pred_store[,j] = inv_logit(post$a_storage + post$lambda[,J]*RI[j])
}

plot(y=apply(pred_store, 2, median), x=RI, ylim=c(0,1), type="l", ylab="Pr (Food Storage)", xlab="RI (z-score)", xaxt='n', col="cornflowerblue",lwd=2, xaxs="i", yaxs="i")
mtext("Food Storage", cex=0.8)

axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)
shade(apply(pred_store, 2, PI, prob=0.9), RI, col=col.alpha("cornflowerblue"))

####### Hunting #######
hunt_cols <- viridis::viridis(n = K_hunt)

plot(NULL, xlim=c(low_RI, high_RI), ylim=c(0,1), ylab="Pr (Dependence on Hunting)", xlab="RI (z-score)", xaxt='n',xaxs="i", yaxs="i")
mtext("Dependence on Hunting", cex=0.8)

for (l in 1:K_hunt) {
pred_hunt <- matrix(NA, nrow=n_samps, ncol=n_preds)

for (j in 1:n_preds) {

 if (l == 1) pred_hunt[,j] = inv_logit( post$c_hunt[,l] - post$lambda[,J-1]*RI[j] ) 
 if (l > 1 & l < K_hunt) pred_hunt[,j] = inv_logit( post$c_hunt[,l] - post$lambda[,J-1]*RI[j] ) - inv_logit( post$c_hunt[,l-1] - post$lambda[,J-1]*RI[j] )
 if (l == K_hunt) pred_hunt[,j] = 1 - inv_logit( post$c_hunt[,l-1] - post$lambda[,J-1]*RI[j] )
}

lines(apply(pred_hunt, 2, median), x=RI, col=hunt_cols[l], lwd=2)
shade(apply(pred_hunt, 2, PI, prob=0.9), RI, col=col.alpha(hunt_cols[l],0.15))
}

axis(1, at=c(low_RI, median_RI, high_RI), labels=(c(low_RI, median_RI,high_RI)-median_RI)/mad_RI)

legend(x=-2, y=1, legend=rev(c("0-5%", "6-15%", "16-25%", "26-35%","36-45%","46-55%","56-65%","66-75%", "76-85%", "86-100%")), lwd=2, lty="solid", col=hunt_cols, bty='n', cex=0.8)

dev.off()
} ## end plot


