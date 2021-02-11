#### F. distichus common garden indicator taxa analysis ###
# Get indicator ASVs between sites and over time in common environment
#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(indicspecies)
library(ggpubr)

#### Read in dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/Combined/CG_FT_rock_fucus_combined_phyloseq_r1500_Jan2021.RDS")
dimnames(tax_table(cg.ft.fr))

#### Compare all communities at time 1 to all communities at time 5 ####
# Get otu table and meta data for sample comparisions 
cgf.1.5 <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "CG" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$sample.number %in% c(1,5)) #subset phyloseq object
cgf.meta <- as.data.frame(unclass(sample_data(cgf.1.5))) %>% mutate_if(is.factor, as.character) # get otu table (community data matrix) where sites/samples are rows and species/ASVs are columns
cgf.otu <- as.data.frame(unclass(otu_table(cgf.1.5)))
cgf.otu.mat <- as.matrix(cgf.otu)
ind.1.5 <- cgf.meta$sampleid
cgf.otu.mat <- cgf.otu.mat[ind.1.5, ] #order matrix to match metadata

cg.sn <- cgf.meta$sample.number #get vector of indicator variable

#### Indispecies analysis 1 vs 5 ####
cg.sn.indval <- multipatt(cgf.otu.mat, cg.sn, control = how(nperm=999))
summary(cg.sn.indval, indvalcomp = T) #view results

#### Import results ###
cgf15.indval <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indvall_all_day1vs5_Feb2021.csv")

# Add taxonomy data
cg.ft.tax <- as.data.frame(unclass(tax_table(cg.ft.fr)))
cgf15.indval.tax <- left_join(cgf15.indval, cg.ft.tax)
cgf15.indval.tax$Sequence <- toupper(cgf15.indval.tax$Sequence)
view(cgf15.indval.tax)
# Save
write.csv(cgf15.indval.tax, file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indvall_all_day1vs5_with_taxonomy_Feb2021.csv")
############################################################################################################################################


############################################################################################################################################
#### Compare all communities at time 1 with all possible comparsions ####
# Get otu table and meta data for sample comparisions 
cgf.1 <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "CG" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$sample.number ==1) #subset phyloseq object
cgf.meta1 <- as.data.frame(unclass(sample_data(cgf.1))) %>% mutate_if(is.factor, as.character) # get otu table (community data matrix) where sites/samples are rows and species/ASVs are columns
cgf.otu1 <- as.data.frame(unclass(otu_table(cgf.1)))
cgf.otu.mat1 <- as.matrix(cgf.otu1)
ind.1 <- cgf.meta1$sampleid
cgf.otu.mat1 <- cgf.otu.mat1[ind.1, ] #order matrix to match metadata

cg.site <- cgf.meta1$origin #get vector of indicator variable

#### Indispecies analysis ####
cg.site.indval <- multipatt(cgf.otu.mat1, cg.site, duleg = T, control = how(nperm=999))
summary(cg.site.indval, indvalcomp = T) #view results

# all comparisons
cg.site.indval.all <- multipatt(cgf.otu.mat1, cg.site, control = how(nperm=999))
summary(cg.site.indval.all, indvalcomp = T) 


#### Import results ###
cg.site.indval <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indvall_site_Feb2021.csv")
cg.site.indval.all <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indvall_all_day1_Feb2021.csv")


# Add taxonomy data
# cg.ft.tax <- as.data.frame(unclass(tax_table(cg.ft.fr)))
cgf.site.indval.tax <- left_join(cg.site.indval, cg.ft.tax)
cgf.site.indval.tax$Sequence <- toupper(cgf.site.indval.tax$Sequence)

cg.site.indval.all.tax <- left_join(cg.site.indval.all, cg.ft.tax)
cg.site.indval.all.tax$Sequence <- toupper(cg.site.indval.all.tax$Sequence)


# Save
write.csv(cgf.site.indval.tax, file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indvall_site_with_taxonomy_Feb2021.csv")
write.csv(cg.site.indval.all.tax, file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indval_all_day1_Feb2021_with_taxonomy.csv")
############################################################################################################################################


############################################################################################################################################
#### Compare all communities at time 5 with all possible comparisons ####
cgf.5 <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "CG" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$sample.number ==5) #subset phyloseq object
cgf.meta5 <- as.data.frame(unclass(sample_data(cgf.5))) %>% mutate_if(is.factor, as.character) # get otu table (community data matrix) where sites/samples are rows and species/ASVs are columns
cgf.otu5 <- as.data.frame(unclass(otu_table(cgf.5)))
cgf.otu.mat5 <- as.matrix(cgf.otu5)
ind.5 <- cgf.meta5$sampleid
cgf.otu.mat5 <- cgf.otu.mat5[ind.5, ] #order matrix to match metadata

cg.end <- cgf.meta5$origin #get vector of indicator variable

#### Indispecies analysis ####
cg.end.indval <- multipatt(cgf.otu.mat5, cg.end, control = how(nperm=999))
summary(cg.end.indval, indvalcomp = T) #view results

#### Import results ###
cg.end.indval <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indvall_end_Feb2021.csv")

# Add taxonomy data
# cg.ft.tax <- as.data.frame(unclass(tax_table(cg.ft.fr)))
cgf.end.indval.tax <- left_join(cg.end.indval, cg.ft.tax)
$Sequence <- toupper(cgf.end.indval.tax$Sequence)

write.csv(cgf.end.indval.tax,file="~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indval_end_Feb2021_with_taxonomy.csv")
############################################################################################################################################


############################################################################################################################################
##### Combine beginning, end, and multiple group comparisons ####
# Filter to greatest indval statistics (> 0.7)
all.1.v.5 <- cgf15.indval.tax %>% filter(indval.stat >= 0.7) # all samples at beginning vs. end of common garden experiment
all.1.v.5$group <- paste("all_hosts_day", all.1.v.5$group, sep=".")  # modify group to indicate comparison

all.origin.comps.1 <- cg.site.indval.all.tax %>% filter(indval.stat >= 0.7) # comparing hosts from different origin sites at start of common garden experiment
all.origin.comps.1$group <- paste(all.origin.comps.1$group, "day.1", sep="_") # modify group to indicate day 1

all.origin.comps.5 <- cgf.end.indval.tax %>% filter(indval.stat >= 0.7) # comparing hosts from different origin sites at end of common garden experiment
all.origin.comps.5$group <- paste(all.origin.comps.5$group, "day.5", sep="_") # modify group to indicate day 5


# Combine results from multiple indval tests
cg.indval.tableS1 <- rbind(all.1.v.5,all.origin.comps.1,all.origin.comps.5)

#### Save final subsetted results table ####
write.csv(cg.indval.tableS1, file = "~/Desktop/Desktop2020/CG_FT/Common_Garden/Data/CG_Indval_TableS1_Feb2021.csv") 
############################################################################################################################################




