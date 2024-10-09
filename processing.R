library(tidyverse)
library(phytools)

######################################
d <- read.csv("data-raw/SCCS-var1-2000.csv", stringsAsFactors = F) # reading in SCCS dataset, should just give reduced data for publication
#d <- d[-187,] # dropping extraneous row

# Selecting study variables from SCCS codebook
d <- d[,c("sccs.","socname", "v149", "v150", "v151", "v152", "v153", "v154", "v155", "v156", "v157", "v158", "v204", "v20", "v820", "v200", "v201")]

##### Bring in location data from DPLACE ################
d_loc <- read.csv("data-raw/location.csv", stringsAsFactors = F)

d_loc <- d_loc %>%
  filter(Source == "Standard cross-cultural sample") %>%
  select(Society.id, Revised.latitude, Revised.longitude) %>% 
  mutate(sccsn = substr(Society.id, 5, nchar(Society.id))) %>%
  mutate(sccsn = as.numeric(sccsn)) %>%
  select(-Society.id)

d <- left_join(d, d_loc, by=c("sccs." = "sccsn"))
names(d)[1] <- "id"

# Manually adding lat/lon for one society where data missing from DPLACE
d$Revised.latitude[d$socname == "Iban"] <- 2
d$Revised.longitude[d$socname == "Iban"] <- 110



sccs_tree <- read.nexus("data-raw/sccs_supertree.nex")
setdiff(sccs_tree$tip.label, d$socname) # checking for discrepancies between phylo tree names and dataframe names

# "MISSING" indicates that genetic data missing for the target pop, position based on linguistic data
d[c(1,2,14,25,27,29,30,33,39,46,48,53,59,63,66,82,93,97,102,106,124,133,134,135,136,137,138,139,140,143,144,146,158,164,167,174,176,178), "socname"] <-  c("Nama_Hottentot", "Kung_Bushmen", "Nkundo_Mongo", "Pastoral_Fulani", "Massa_Masa", "Fur_Darfur", "Otoro_Nuba", "MISSING_Kafa_Kaffa", "Kenuzi_Nubians","Rwala_Bedouin", "Gheg_Albanians", "Yurak_Samoyed", "Punjabi_West", "Uttar_Pradesh", "Khalka_Mongols", "Negri_Sembilan", "MISSING_Kimam", "New_Ireland", "Mbau_Fijians", "Western_Samoans", "Copper_Eskimo", "MISSING_Twana", "MISSING_Yurok", "MISSING_Pomo_Eastern", "Yokuts_Lake", "Paiute_North", "MISSING_Klamath","MISSING_Kutenai", "Gros_Ventre", "MISSING_Omaha", "MISSING_Huron", "MISSING_Natchez", "Cuna_Tule", "Carib_Barama", "Cubeo_Tucano", "MISSING_Nambicuara", "MISSING_Timbira", "MISSING_Botocudo")

setdiff(sccs_tree$tip.label, d$socname) # should be 0

sccs_tree <- drop.tip(sccs_tree, subset(sccs_tree$tip.label, !(sccs_tree$tip.label %in% d$socname)))

# Putting the phylogeny and the data frame in the same order, this is very important!
d <- d[match(sccs_tree$tip.label, d$socname),]

# Now we can remove the "MISSING" from the labels
d$socname <- ifelse(substr(d$socname, 1, nchar("MISSING")) == "MISSING", substr(d$socname, nchar("MISSING")+2, nchar(d$socname)), d$socname)

# Do the same for phylogeny
sccs_tree$tip.label <- ifelse(substr(sccs_tree$tip.label, 1, nchar("MISSING")) == "MISSING", substr(sccs_tree$tip.label, nchar("MISSING")+2, nchar(sccs_tree$tip.label)), sccs_tree$tip.label)

# Export processed files for analysis
write.csv(d, "data_analysis.csv", row.names = F)
write.tree(sccs_tree, "SCCS_supertree.tre")
