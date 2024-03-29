#### Common Garden taxa summary barplots and within-between comparison of Jaccard dissimilarities ####

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(gridExtra)
library(rstatix)

#### Read combined dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/CG_FT_combined_phyloseq_r1500_Jan2021.RDS")
cg.ft.tax <- as.data.frame(unclass(tax_table(cg.ft.fr)))

# Subset to common garden and fucus samples and create dataframe
cg.ps <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "CG")
#cg.df <- psmelt(cg.ps) # create dataframe 
#cg.df <- cg.df %>% mutate_if(is.factor, as.character)
#write.csv(cg.df, file="~/Desktop/Desktop2020/CG_FT/Data/CG-FT_CG_dataframe_r1500_Jan2021.RDS")
cg.df <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Data/CG-FT_CG_dataframe_r1500_Jan2021.RDS")
cg.df <- cg.df %>% mutate_if(is.factor, as.character)
cg.f <- cg.df %>% filter(type == "fucus")

#### Get relative abundance ####
cg.f <- cg.f %>%  mutate(relabund = Abundance/1500 * 100)

#### Find most abundant ASVS ####
cg.f.4 <- cg.f %>% filter(relabund > 5)
length(unique(cg.f.4$OTU))
length(unique(cg.f.4$Family))  #10
length(unique(cg.f.4$Genus)) #25
genera.4 <- as.character(unique(cg.f.4$Genus)) #Group by genera
sort(genera.4)

#### 26 colors ordered by Family ####
c26.ordered <- c("snow4",
                 "#6185e8",
                 "#d47aa3",
                 "#dc375f",
                 "darkred",
                 "#e28e24",
                 "#e05621",
                 "#ba4334",
                 "midnightblue",
                 "#4ca578",
                 "#a0b92c",
                 "#ab88d3",
                 "#b865af",
                 "#905ed9",
                 "#3a8b37",
                 "#57b946",
                 "#c68337",
                 "#ddb82d",
                 "lightgoldenrod",
                 "lightskyblue",
                 "#5e56ad",
                 "#3ebab5",
                 "#5f95cc",
                 "#c954b8",
                 "#9A4C4C",
                 "#898b48",
                 "#925218")

#Look at colors
scales::show_col(c26.ordered)

# Update Family/Genus names 
cg.f$Genus2 <- ifelse(cg.f$Genus %in% genera.4, cg.f$Genus, "Other")
unique(cg.f$Genus2)
cg.f$Genus3 <- ifelse(cg.f$Genus2 %in% genera.4, paste(cg.f$Family, cg.f$Genus2, sep="_"), "Other")
sort(unique(cg.f$Genus3))


#### Fix awkward taxonomy assignments
unique(gsub("_unknown_.*", "_Unknown",cg.f$Genus3))
cg.f$Genus3 <- gsub("_unknown_.*", "_Unknown",cg.f$Genus3)
cg.f$Genus3[cg.f$Genus3 == "Flavobacteriaceae_Nonlabens"] <- "Flavobacteriaceae_Unknown"
cg.f$Genus3[cg.f$Genus3 == "unknown_Family_Unknown"] <- "Bacteroidia_Unknown"
cg.f$Genus3[cg.f$Genus3 == "unknown_Bacteroidia_Unknown"] <- "Bacteroidia_Unknown"
sort(unique(cg.f$Genus3))
sort(unique(cg.f$sampleid))

# Add column for sample day and host individual
cg.f.id <- cg.f %>% select(individual, origin, sample.number) %>% distinct()
cg.f.id <- cg.f.id %>% ungroup() %>% group_by(origin, sample.number)  %>% mutate(id2 = seq_along(individual))
cg.f2 <- left_join(cg.f, cg.f.id)
unique(cg.f2$id2)
cg.f2$id3 <- paste(cg.f2$sample.number, cg.f2$id2, sep="-")

#### Make stacked barplots ####
cg.f2$origin <- factor(cg.f2$origin, levels=c("PB", "NB", "WB west wall", "WB high", "WB low"))
cg.f.1.5 <- cg.f2 %>% filter(sample.number %in% c(1,2,5))


#### Arrange taxonomic assignments in alphabetical order with "Other" at the beginning ####
paste0(sort(unique(cg.f.1.5$Genus3)), collapse="','") 
cg.f.1.5$Genus3 <- factor(cg.f.1.5$Genus3, levels=c('Other','Bacteroidia_Unknown','Flavobacteriaceae_Algitalea','Flavobacteriaceae_Croceitalea','Flavobacteriaceae_Dokdonia','Flavobacteriaceae_Maribacter','Flavobacteriaceae_Psychroserpens','Flavobacteriaceae_Unknown','Hyphomonadaceae_Litorimonas','Pirellulaceae_Blastopirellula',
                                                    'Pseudoalteromonadaceae_Pseudoalteromonas','Rhodobacteraceae_Octadecabacter','Rhodobacteraceae_Sulfitobacter','Rhodobacteraceae_Unknown','Rubritaleaceae_Luteolibacter','Rubritaleaceae_Rubritalea','Saprospiraceae_Lewinella','Saprospiraceae_Membranicola','Saprospiraceae_Unknown',
                                                    'Shewanellaceae_Shewanella','Syntrophomonadaceae_Syntrophomonas','Thiohalorhabdaceae_Granulosicoccus'))

sort(unique(cg.f.1.5$Genus3))


# Plot and save
pdf("~/Desktop/Desktop2020/CG_FT/Common_Garden/Figures/CG_barplot_Family_Genus_1-5_ordered_8Feb2021.pdf", 
    width = 13, # define plot width and height
    height = 6)
c <- ggplot(data=cg.f.1.5, aes(x=id3, y=relabund, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  facet_grid(~origin, scales="free", space="free")+
  theme_classic(base_size = 12)+
  theme(panel.spacing.x=unit(0, "lines")) +
  theme(strip.text.x = element_text(size = 12)) +
  scale_fill_manual(values = c26.ordered) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size =10)) +
  theme(axis.text.y=element_text(size =12)) +
  theme(legend.position = "bottom", legend.text=element_text(size=12)) +
  ylab("Relative abundance") +
  xlab("Sample day_Host individual") +
  guides(fill=guide_legend(ncol=3)) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) 
c
dev.off()
##################################################################

##################################################################
#### Compare Jaccard dissimilarity of samples from the same site 
# or between different sites of origin over time in the common garden

#### Distance matrix ####
cg.f.ps <- subset_samples(cg.ps, sample_data(cg.ps)$type == 'fucus')
cg.f.ps <- subset_samples(cg.f.ps, sample_data(cg.f.ps)$sample.number %in% c(1,5))
cg.f.dist <- phyloseq::distance(cg.f.ps, method='jaccard', type='samples') #make distance object
cg.f.distmat <- as.matrix(cg.f.dist) #convert to matrix

#### Jaccard dissimilarity ####
# 0 means both samples share exact the same species
# 1 means both samples have no species in common

cg.f.distmeta <- as.data.frame(unclass(sample_data(cg.f.ps)))
colnames(cg.f.distmeta)
cg.f.distmeta2 <- cg.f.distmeta %>% select(sampleid, type, date, origin, sample.number, individual)

#### Create dataframe of pairwise comparisons ####
pairwise.cg <-subset(reshape2::melt(cg.f.distmat, na.rm = TRUE))
colnames(pairwise.cg) <- c('s.1','s.2','jacc')
pairwise.cg <- pairwise.cg %>% filter(! s.1 == s.2) # remove pairwise comparisons of the same sample

colnames(cg.f.distmeta2)
names(cg.f.distmeta2)[names(cg.f.distmeta2) == 'sampleid'] <- 's'
cg.1 <- cg.f.distmeta2
colnames(cg.1) <- paste(colnames(cg.1), "1", sep = ".")
cg.2 <- cg.f.distmeta2
colnames(cg.2) <- paste(colnames(cg.2), "2", sep = ".")

pairwise.cg.m <- left_join(pairwise.cg, cg.1)
pairwise.cg.m <- left_join(pairwise.cg.m, cg.2)

pairwise.cg.m.filt2 <- pairwise.cg.m
pairwise.cg.m.filt2$comp <- paste(pairwise.cg.m.filt2$sample.number.1, pairwise.cg.m.filt2$sample.number.2, sep="-")
pairwise.cg.m.filt2$within_btwn <- ifelse(pairwise.cg.m.filt2$origin.1 == pairwise.cg.m.filt2$origin.2, "within site", "between sites")
pairwise.cg.m.filt2$origin.1 <- factor(pairwise.cg.m.filt2$origin.1, levels=c("PB", "NB", "WB west wall", "WB high", "WB low"))

pairwise.select <- pairwise.cg.m.filt2 %>% filter(comp %in% c("1-1","5-5"))
pairwise.select$comp2 <- ifelse(pairwise.select$comp == "1-1", "1", "5")
pairwise.select$within_btwn2 <- factor(pairwise.select$within_btwn, levels=c("within site","between sites"))
pairwise.select$origin.1 <- factor(pairwise.select$origin.1, levels=c("PB", "NB", "WB west wall", "WB high", "WB low"))

#### Calculate paired t-test for each combo ####
pwt <- pairwise.select %>% ungroup() %>% group_by(within_btwn2, origin.1) %>% pairwise_t_test(jacc ~ comp2)

#### Save pairwise-t-test results ####
pwt <- pwt %>% select(-c(p:p.signif))
colnames(pwt) <- c("origin","comparison","beta_metric","group1","group2", "n1","n2", "p.adjust", "p.adjust.signif")
write.csv(pwt, file="~/Desktop/Manuscripts_2021/CG_FT_2021/CG-FT_Supplementary_Info_2021/CG_Jaccard_comparison_TableS1_Feb2021.csv")


# Creat dataframe of paired T-test significance values to annotate plot
paste(shQuote(unique(pairwise.select$origin.1)), collapse=", ")
cg.origin <- c('PB', 'NB', 'WB west wall', 'WB high', 'WB low', 'PB', 'NB', 'WB west wall', 'WB high', 'WB low')
strrep(c("'within site', ","'between sites', "), 5)
cg.compare <- c('within site', 'within site', 'within site', 'within site', 'within site', 'between sites', 'between sites', 'between sites', 'between sites', 'between sites')
sig <- c(NA,NA,"*",NA,"*","*",NA,NA,NA,NA )
sig.df <- data.frame(cg.compare, cg.origin, sig)
colnames(sig.df) <- c("within_btwn2", "origin.1", "sig")

#### Plot Jaccard within and between sites ####
pdf("~/Desktop/Desktop2020/CG_FT/Common_Garden/Figures/CG_within_between_Jaccard_Feb2021.pdf", 
    width = 13, # define plot width and height
    height = 4)
wbp <- ggplot(pairwise.select, aes(x=comp2, y = jacc)) +
  geom_jitter(width=0.2) + 
  geom_boxplot(fill='#A4A4A4', alpha= 0.2, color="black") +
  facet_grid(~within_btwn2+origin.1) +
  theme_classic()+
  theme(panel.spacing.x=unit(0, "lines")) +
  xlab("Sample day") +
  ylab("Jaccard dissimilarity index") +
  theme(strip.text.x = element_text(size = 12)) +
  theme(axis.text=element_text(size =12), axis.title=element_text(size=12)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) 
wbp + geom_text(x = 1.5, y = 0.965, aes(label = sig), data = sig.df, size = 9)
dev.off()

# Add statistics annoation to plot
wbp2 <- wbp + geom_text(x = 1.5, y = 0.965, aes(label = sig), data = sig.df, size = 9)
wbp2

#### Combine Stacked Barplot and Jaccard Within-Between comparisons ###
library(ggpubr)

# Save with figure labels
pdf(file="~/Desktop/Manuscripts_2021/CG_FT_2021/CG-FT_Figures_Resubmission_2021/CG-FT_Figure3_CG_Barplots_Jaccard_Feb2021.pdf",
    width =12, height = 9)
ggarrange(c, wbp2, ncol = 1, labels = c("A","B"), heights=c(1.75,1))
dev.off()

