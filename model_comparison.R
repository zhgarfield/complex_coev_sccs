# Model Fitting and Comparison 

# NOTE: This script relies on "rstan" and convience functions from the "rethinking" package. For installation instructions see https://github.com/rmcelreath/rethinking)
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("phytools")) install.packages("phytools")
if (!require("qgraph")) install.packages("qgraph") # Correlation network plot
library(rethinking)

# Read in data
d <- read.csv("data_analysis.csv", stringsAsFactors = F)
sccs_tree <- ape::read.tree("SCCS_supertree.tre")

setdiff(sccs_tree$tip.label, d$socname) # checking for discrepancies between phylo tree names and dataframe names, should return "character(0)" if all is well

#### Data dictionary ####
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
dist_mat <- cophenetic.phylo(sccs_tree) # pairwise distance matrix using branch lengths
dist_mat <- dist_mat / max(dist_mat) # scaling matrix to [0,1]

#### Organizing data ##
storage <- ifelse(d$v20 > 1, 1, 0) # dichtomizing food storage
storage <- ifelse(is.na(storage), -99, storage) # flagging missing storage observations
hunting <- d$v204 + 1 # making the lowest value 1
K_hunt <- max(hunting) # number of ordinal levels

J <- 12 # number of variables
K <- max(d$v149) # number of ordinal levels
N <- nrow(d) # num observations

# Organize data into list for stan
data_list <- list(
  y = cbind(d$v151, d$v150, d$v156, d$v152, d$v149, d$v153, d$v154, d$v155, d$v157, d$v158), # all the Murdock complexity variables
  hunting = hunting,
  storage = storage,
  N = N,
  J = J,
  K = K,
  K_hunt = K_hunt,
  dist_mat = dist_mat,
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
#fit_m3 <- stan( file="stan_models/m3.stan", data=data_list, iter=n_iter, chains=n_chains, cores=n_chains, init="0", control=list(adapt_delta=0.9) )
#saveRDS(fit_m3, "fit_m3.rds")

#### Model comparison with PSIS LOOCV #################

## Bring in previously fit model objects and estimate PSIS-LOO
fit_m0 <- readRDS("fit_m0.rds"); loo_m0 <- loo(fit_m0)
fit_m1 <- readRDS("fit_m1.rds"); loo_m1 <- loo(fit_m1)
fit_m2 <- readRDS("fit_m2.rds"); loo_m2 <- loo(fit_m2)
fit_m3 <- readRDS("fit_m3.rds"); loo_m3 <- loo(fit_m3)
fit_mS <- readRDS("fit_mS.rds"); loo_mS <- loo(fit_mS)

## Difference in elpd from each model to the saturated model
elpd_df <- data.frame(model = paste0("M",0:3), elpd_diff = rep(NA,4), se_diff = rep(NA,4))
models <- c("m0","m1", "m2", "m3")

for (m in 1:length(models)) {
  
  m_compare <- paste0("loo_", models[m])
  
  # Difference in elpd, make sure we get the right sign of the difference depending on the comparison
  elpd_df$elpd_diff[m] <- ifelse( rownames(loo::loo_compare(loo_mS, eval(as.symbol(m_compare))))[2] == "model2", loo::loo_compare(loo_mS, eval(as.symbol(m_compare)))[2,][1], loo::loo_compare(loo_mS, eval(as.symbol(m_compare)))[2,][1] * -1 )
  
  # SE of difference
  elpd_df$se_diff[m] <- loo::loo_compare(loo_mS, eval(as.symbol(m_compare)))[2,][2]
}


## Difference in elpd (M2, M3) - M4
elpd_df2 <- data.frame(model = c("m2"), elpd_diff = NA, se_diff = NA)

# Difference in elpd, make sure we get the right sign of the difference depending on the comparison
elpd_df2$elpd_diff[1] <- ifelse( rownames(loo::loo_compare(loo_m2, loo_m3))[2] == "model2", loo::loo_compare(loo_m2, loo_m3)[2,][1], loo::loo_compare(loo_m2, loo_m3)[2,][1] * -1 ) 

# SE of difference
elpd_df2$se_diff[1] <- loo::loo_compare(loo_m2, loo_m3)[2,][2]

## Plot both difference plots  
library(patchwork)  

svg("elpd_diff.svg", width=10, height=7, pointsize=12)

diff1 <- ggplot(elpd_df, aes(x=as.numeric(elpd_diff), y=model)) + geom_point(size=4) + geom_errorbarh(aes(xmin=elpd_diff - se_diff*2, xmax=elpd_diff + se_diff*2), height=0, lwd=2) + geom_vline(aes(xintercept=0), linetype="dashed") + ylab("") + xlab(expression(paste(Delta, "(ELPD - MS ELPD)"))) + ggtitle("Model Performance") + theme_bw(base_size=24) + theme(panel.grid.minor = element_blank())

diff2 <- ggplot(elpd_df2, aes(x=as.numeric(elpd_diff), y=model)) + geom_point(size=4) + geom_errorbarh(aes(xmin=elpd_diff - se_diff*2, xmax=elpd_diff + se_diff*2), height=0, lwd=2) + geom_vline(aes(xintercept=0), linetype="dashed") + ylab("") + xlab(expression(paste(Delta, "(M3 ELPD - M2 ELPD)"))) + ggtitle("") + theme_bw(base_size=24) + theme(panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(limits=c(-15,55))

diff1 + diff2

dev.off()

## Compare candidate models using different weighting methods
#stacking_weights <- loo::loo_model_weights(list(loo_m1, loo_m2, loo_m3, loo_m4), method="stacking")
#PBMA_weights <- loo::loo_model_weights(list(loo_m1, loo_m2, loo_m3, loo_m4), method="pseudobma")
#WAIC_weights <- rethinking::compare(fit_m1, fit_m2, fit_m3, fit_m4)

#####################################################
#### Phylogentic correlation network ################
fit_mS <- readRDS("fit_mS.rds")
post <- extract.samples(fit_mS)
n_samps <- length(post$lp__)

cor_phy <- apply(post$Rho_phy, 2:3, median)
rownames(cor_phy) <- c('Wrt','FoR','Agr','Urb','L.S.','Lnt','Mny','DoP','PlI','ScS','Hnt','FdS')
colnames(cor_phy) <- c('Wrt','FoR','Agr','Urb','L.S.','Lnt','Mny','DoP','PlI','ScS','Hnt','FdS')

svg("cor_network.svg", width=7, height=4, pointsize=12)

qgraph(cor_phy,
       graph="cor",
       layout="spring",
       color = c("#87ceeb", "slategray"),
       posCol = "#70a494",
       negCol = "#de8a5a",
       groups = c("TPC",rep("RI",3), rep("TPC",3), "RI", "TPC", "TPC", "RI", "RI"),
       vsize=9,
       borders = T,
       minimum = 0.1,
       fade = T )

dev.off()


####### NOTE: EVERYTHING BELOW THIS IS A WIP!! ################################

#008080,#70a494,#b4c8a,#f6edbd,#edbb8a,#de8a5a,#ca562c
###############################################################################
##### ELPD breakdown by region/subsistence mode ###############################
library(ggridges)

loo_m3 <- loo_list[[4]]
m3_pw <- as.data.frame( matrix( loo_m3$pointwise[,1] , nrow=N, ncol=9 ) )

names(m3_pw) <- c("FoR","DoP","Urb","Wrt","T.S.","LnT", "Mny","PlI", "ScS")
m3_pw$sum_elpd <- rowSums(m3_pw)
m3_pw$socname <- d$socname

m3_pw <- left_join( select(m3_pw,socname, sum_elpd),  select(d, socname, v820, v200, v201) )

### Summary score of complexity variables (and agr, following original)
m3_pw$sum_score <- d$v149 + d$v150 + d$v151 + d$v152 + d$v153 + d$v154 + d$v155 + d$v156 + d$v157 + d$v158

### Label Murdock's world region codes
m3_pw$subtype <- factor(m3_pw$v820, labels=c("Gathering", "Hunting", "Fishing", "Incipient Agriculture", "Domestic Animals", "Extensive Agriculture", "Intensive Agriculture", "Other"))

m3_pw$Region <- factor(m3_pw$v200, labels=c("Africa","Circum-Mediterranean", "East Eurasia",  "Insular Pacific", "North America", "South America"))

m3_pw$Area <- factor( paste0(m3_pw$v200, m3_pw$v201), labels=c(
  "Upper Nile",
  "African Hunters",
  "S.Afr. Bantu",
  "C. Bantu",
  "N.E. Bantu",
  "Equit. Bantu",
  "Guinea Coast",
  "W. Sudan",
  "Nigerian Pt.",
  "E. Sudan",
  "Sem.Near E",
  "Ethiop-Horn",
  "Mosl.  Sudan",
  "Sahara",
  "N. Africa",
  "S. Europe",
  "N.W. Europe",
  "E. Europe",
  "Turk-Caus.",
  "S.E. Asia",
  "Middle East",
  "Cntrl. Asia",
  "Arctic Asia",
  "East Asia",
  "Himalayas",
  "N.-C. India",
  "South India",
  "Indian Ocn.",
  "Assam-Burma",
  "E. Polynesia",
  "Phl-Formosa",
  "W. Indonesia",
  "E. Indonesia",
  "New Guinea",
  "Australia",
  "Micronesia",
  "W. Melanesia",
  "E. Melanesia",
  "W. Polynesia",
  "C.Mexico",
  "Arctic Amer",
  "N.W. Coast",
  "California",
  "Gr.Basin-Pl",
  "Plains",
  "Prarie",
  "E. Woodlands",
  "Southwest",
  "N.W. Mexico",
  "E. Brazil",
  "C.America",
  "Caribbean",
  "Guiana",
  "Lower Amaz",
  "Inner Amaz",
  "Andes",
  "Chile-Pata",
  "Gran Chago",
  "Mato Grosso"
))

ggplot(m3_pw, aes(x=sum_elpd, y=fct_reorder(Region, sum_elpd))) + geom_density_ridges(scale=0.8, color=NA) + geom_jitter(height=0.05, width=0, alpha=0.6) + xlab("ELPD") + ylab("") + theme_bw(base_size=15)

ggplot(filter(m3_pw, subtype != "Other"), aes(x=sum_elpd, y=fct_reorder(subtype, sum_elpd))) + geom_density_ridges(scale=0.8, color=NA) + geom_jitter(height=0.05, width=0, alpha=0.6) + xlab("ELPD") + ylab("") + theme_bw(base_size=15)

ggplot(m3_pw, aes(x=sum_elpd, y=fct_reorder(Area, sum_elpd))) + facet_wrap(~Region, scales="free_y") + geom_jitter(height=0.05, width=0, alpha=0.6) + xlab("ELPD") + ylab("") + theme_bw(base_size=15)

ggplot(m3_pw, aes(y=sum_elpd, x=sum_score)) + geom_point(alpha=0.6) + theme_bw(base_size=15) + xlab("Complexity Summary Score") + ylab("ELPD")

###############################################################################
###############################################################################
##### GGM comparison & exploratory factor analysis ############################
if (!require("psych")) install.packages("psych")

# Extract the median correlation matrix from m4
m4_fit <- readRDS("fit_m4.rds")

med_cor <- apply(extract.samples(m4_fit)$Rho, 2:3, median)
rownames(med_cor) <- c("Agr","FoR","DoP","Urb","Wrt","T.S.","LnT", "Mny","PlI", "ScS", "Hnt", "FdS"); colnames(med_cor) <- rownames(med_cor)

# Exploratory factor analysis suggests two factor solution
fa.parallel(med_cor, n.obs = 186, fa = "fa", n.iter = 1000, error.bars = TRUE, fm="gls")

# Identify the two-factor solution
efa <- fa(med_cor, n.obs=186, nfactors=2, fm="gls", rotate="geominQ")
efa

# Breaking down elpd, where do the two models perform best?
# first, extract the posterior complexity latent score
post_m2 <- extract.samples( readRDS("fit_m2.rds") )
latent_C <- apply(post_m2$eta_C, 2, median)

d$submode <- ifelse(d$v820 == 8, 7, d$v820) # recode 2 "trade" cases as intensive agriculture


m3_loo_pw <- loo_list[[4]]$pointwise[,1]
m3_loo_pw <- data.frame( elpd = m3_loo_pw, socname = rep(d$socname,9), latent_C = rep(latent_C,9), submode = rep(d$submode, 9), region = rep(d$v200,9)  )

m3_loo_pw <- loo_list[[5]]$pointwise[,1]
m3_loo_pw <- data.frame( elpd = m3_loo_pw, socname = rep(d$socname,9), latent_C = rep(latent_C,9), submode = rep(d$submode, 9), region = rep(d$v200,9)  )


# Average across complexity measures
m3_loo_pw_sum <- m3_loo_pw %>% group_by(socname) %>% summarise(elpd = sum(elpd), latent_C = mean(latent_C), submode = mean(submode), region = mean(region))
m3_loo_pw_sum$sub2 <- ifelse(m3_loo_pw_sum$submode %in% c(1:3), "Hunting, Gathering, Fishing", m3_loo_pw_sum$submode)
m3_loo_pw_sum$sub2 <- ifelse(m3_loo_pw_sum$submode %in% c(4,6), "Horticulture/Incipient Agriculture", m3_loo_pw_sum$sub2)
m3_loo_pw_sum$sub2 <- ifelse(m3_loo_pw_sum$submode == 5, "Pastoralism", m3_loo_pw_sum$sub2)
m3_loo_pw_sum$sub2 <- ifelse(m3_loo_pw_sum$submode == 7, "Intensive Agriculture", m3_loo_pw_sum$sub2)

ggplot(m3_loo_pw_sum, aes(x=fct_reorder(sub2, elpd), sub2, y=elpd)) + geom_point(aes(color=sub2),alpha=0.6) + geom_violinhalf(aes(color=sub2, fill=sub2), alpha=0.2) + coord_flip() + xlab("") + theme_bw() + theme(legend.position = "none")

m3_loo_pw_sum$region <- as.factor(m3_loo_pw_sum$region)

ggplot(m3_loo_pw_sum, aes(x=fct_reorder(region, elpd), region, y=elpd)) + geom_point(aes(color=region),alpha=0.6) + geom_violinhalf(aes(color=region, fill=region), alpha=0.2) + coord_flip() + xlab("") + theme_bw() + theme(legend.position = "none")

#### Exploratory comparison: m0 - m7 ################
# Difference in elpd from each model to the top model
elpd_df <- data.frame(model = rep(NA,8), elpd_diff = rep(NA,8), se_diff = rep(NA,8))

for (m in 1:8) {
  elpd_df$model[m] = paste0("m",m-1)
  elpd_df$elpd_diff[m] = loo::loo_compare( loo_list[[5]], loo_list[[m]] )[2,][1]
  elpd_df$se_diff[m] = loo::loo_compare( loo_list[[5]], loo_list[[m]] )[2,][2]
}

ggplot(elpd_df[-c(5),], aes(x=as.numeric(elpd_diff), y=model)) + geom_point() + geom_errorbarh(aes(xmin=elpd_diff - se_diff*2, xmax=elpd_diff + se_diff*2), height=0) + geom_vline(aes(xintercept=0)) + theme_bw(base_size=18) + ylab("") + xlab(expression(paste(Delta, "(elpd - m4 elpd)"))) + ggtitle("Model Performance") + theme()






p <- ggtree(sccs_tree)


d_summary <- data.frame(socname = d$socname)

p3 <- facet_plot(p, panel='bar', data=d2, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y), size=3, color='blue4') 

## Descriptive dot plot for Murdock codes
d_codes <- as.data.frame(data_list$y)
colnames(d_codes) <- c("Wrt", "FoR", "Agr", "Urb", "T.S.", "LnT", "Mny", "DoP", "PlI", "ScS")

d_codes$socname <- d$socname
d_codes$socname <- ifelse(substr(d_codes$socname, 1, 7) == "MISSING", substr(d_codes$socname, 9, nchar(d_codes$socname)), d_codes$socname) # remove MISSING labels from phylogeny
d_codes$socname <- sub("_"," ", d_codes$socname) # replace underscores with spaces

code_long <- d_codes %>% pivot_longer(-socname) %>% group_by(socname) %>% mutate(sum_score = sum(value))

code_long$group <- ifelse(code_long$name %in% c("FoR", "Agr", "Urb", "DoP"), "Resource Intensification", "Techno-political Complexity")

pdf(file="data_long.pdf", height=11, width=8.5)
ggplot(code_long, aes(y=fct_reorder(socname, sum_score), x=name)) + geom_tile(aes(fill=as.character(value))) + scale_fill_viridis_d() + theme_minimal(base_size=5.5) + theme(axis.ticks = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_blank(),panel.spacing.x = unit(0,"line"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab("") + xlab("") + coord_equal() + geom_vline(xintercept = 4.5, col="red")
dev.off()




