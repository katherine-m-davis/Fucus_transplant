#### Common Garden differences over time ####

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(viridis)
library(gridExtra)

#### Read combined dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/Combined/CG_FT_rock_fucus_combined_phyloseq_r1500.RDS")

# Only Fucus
cg.ft.f <-subset_samples(cg.ft.fr, ! sample_data(cg.ft.fr)$type %in% c("rock","seawater"))

# Only samples taken on day 1 and 5
cg.ft.f15 <-subset_samples(cg.ft.f, sample_data(cg.ft.f)$sample.number %in% c(1,5))
cg.ft.f15 <-subset_samples(cg.ft.f15, sample_data(cg.ft.f15)$study=="CG")


#### Ordinate and plot #####
cg.ft.f.pcoa15 <- ordinate(cg.ft.f15, method="PCoA", distance="bray") # Initial timepoint
sample_data(cg.ft.f15)$sample.number<- as.factor(sample_data(cg.ft.f15)$sample.number)


#### Figure 2. PCoA plot of CG samples by site on day 1 and 5 ####
# Colors + Shapes
cg.palette <- c("#999999", "#E69F00", "#56B4E9", "darkblue", "cyan4")
day.shape <- c(16, 8)

# Plot and save
pdf(file="~/Desktop/Desktop2020/Manuscripts/CG_FT_Manuscript_2020/CG_FT_Manuscript_Figures/CG-FT_site-morphotype_PCoA_Jan21.pdf",
    width = 5.5, height = 3.5)
sm <- plot_ordination(cg.ft.f15, cg.ft.f.pcoa15, type = "samples", color="origin_morph", shape= "sample.number") +
  theme_classic(base_size = 12) + 
  scale_color_manual(values=cg.palette) + 
  geom_point(size=3) +
  scale_shape_manual(values = day.shape) +
  labs(color='Site (Morphotype)', shape = 'Sample day') 
sm
dev.off()
##################################################################

##################################################################
#### Are samples in Common Garden at day 1 differentiated by site of origin ####
# Only controls from transplant experiment
cg.f <-subset_samples(cg.ft.f, sample_data(cg.ft.fr)$study =="CG")
cg.f.1 <- subset_samples(cg.f, sample_data(cg.f)$sample.number %in% c(0,1))

cg.f.distmeta.1 <- as.data.frame(unclass(sample_data(cg.f.1)))
cg.f.distmeta.1 %>% group_by(origin, sample.number) %>% summarise(n())

cg1 <- as.character(cg.f.distmeta.1$sampleid)

# Make distance matrix
cg.f.dist.1 <- distance(cg.f.1, method='bray', type='samples') #make distance object
cg.f.distmat.1 <- as.matrix(cg.f.dist.1) #convert to matrix
# reorder distance matrix to match meta data
cg.f.distmat.1 <- cg.f.distmat.1[cg1,cg1]

#### PERMANOVA site differences at start of common garden ####
permanova.cg1.site <- adonis(cg.f.distmat.1 ~ origin , data = cg.f.distmeta.1)
permanova.cg1.site
#----------------------------------------------------------------# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origin     4    5.5339 1.38348  8.6352 0.48281  0.001 ***
# Residuals 37    5.9279 0.16021         0.51719           
# Total     41   11.4618                 1.00000         
#----------------------------------------------------------------# 
##################################################################

##################################################################
#### Are samples in Common Garden at day 5 differentiated by site of origin ####
cg.f.5 <- subset_samples(cg.f, sample_data(cg.f)$sample.number == 5)

cg.f.distmeta.5 <- as.data.frame(unclass(sample_data(cg.f.5)))
cg.f.distmeta.5 %>% group_by(origin, sample.number) %>% summarise(n())

# Make distance matrix
cg.f.dist.5 <- distance(cg.f.5, method='bray', type='samples') #make distance object
cg.f.distmat.5 <- as.matrix(cg.f.dist.5) #convert to matrix
# reorder distance matrix to match meta data
cg5 <- as.character(cg.f.distmeta.5$sampleid)
cg.f.distmat.5 <- cg.f.distmat.5[cg5,cg5]

#### PERMANOVA site differences at end of common garden ####
permanova.cg5.site <- adonis(cg.f.distmat.5 ~ origin , data = cg.f.distmeta.5)
permanova.cg5.site
#----------------------------------------------------------------# 
#           Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# origin     4    3.4676 0.86691  5.7874 0.5365  0.001 ***
# Residuals 20    2.9958 0.14979         0.4635           
# Total     24    6.4634                 1.0000 
#----------------------------------------------------------------# 
##################################################################

##################################################################

