#### Examining F. distichus microbiota differences by site and morphology ####

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(viridis)
library(gridExtra)

#### Read combined dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/Combined/CG_FT_rock_fucus_combined_phyloseq_r1500.RDS")

# Subset to only Fucus
cg.ft.f <-subset_samples(cg.ft.fr, ! sample_data(cg.ft.fr)$type %in% c("rock","seawater"))

# Subset to only controls from transplant experiment
cg.ft.f <-subset_samples(cg.ft.f, ! sample_data(cg.ft.f)$individual =="transplant")
cg.ft.f.meta <- as.data.frame(unclass(sample_data(cg.ft.f)))

# Subset to only samples taken initially
cg.ft.f0 <-subset_samples(cg.ft.f, sample_data(cg.ft.f)$sample.number %in% c(0,1))

##################################################################
#### Ordinate and plot #####
cg.ft.f.pcoa0 <- ordinate(cg.ft.f0, method="PCoA", distance="bray") # Initial timepoint

# Colors + Shapes
cg.palette <- c("#999999", "#E69F00", "#56B4E9", "darkblue", "cyan4")
day.shape <- c(15, 16, 8)

sample_data(cg.ft.f0)$sample.number <- as.factor(sample_data(cg.ft.f0)$sample.number)
cg.ft.f0.meta <- as.data.frame(unclass(sample_data(cg.ft.f0)))
cg.ft.f0.meta %>% group_by(origin) %>% summarise(n())

# Plot
plot_ordination(cg.ft.f0, cg.ft.f.pcoa0, type = "samples", color="origin", shape= "sample.number") +
  theme_classic(base_size = 12) + 
  scale_color_manual(values=cg.palette) + 
  geom_point(size=4) +
  scale_shape_manual(values = day.shape)
##################################################################

##################################################################
######## PERMANOVA by site and morphology at initial timepoint ####
cg.f.dist0 <- distance(cg.ft.f0, method='bray', type='samples') #make distance object
cg.f.distmat0 <- as.matrix(cg.f.dist0) #convert to matrix
sample.order <- rownames(cg.f.distmat0) #save sampleid order from rownames of distance matrix

#### Filter metadata to match distance matrices ####
cg.f.distmeta0 <- as.data.frame(unclass(sample_data(cg.ft.f0))) #get metadata
cg.f.distmeta0 %>% group_by(origin) %>% summarise(n()) # get sample sizes per site

# reorder distance matrix to match meta data
m <- as.character(cg.f.distmeta0$sampleid) #get list of sample ids
cg.f.distmat0 <- cg.f.distmat0[m,m]

# Are communities shaped by morphology?
permanova.morph <- adonis(cg.f.distmat0 ~ morph , data = cg.f.distmeta0)
permanova.morph
#----------------------------------------------------------------# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# morph      2    5.8162 2.90810  16.237 0.39858  0.001 ***
# Residuals 49    8.7762 0.17911         0.60142           
# Total     51   14.5924                 1.00000            
#----------------------------------------------------------------# 

# Are communities shaped by site of origin?
permanova.site <- adonis(cg.f.distmat0 ~ origin , data = cg.f.distmeta0)
permanova.site
#----------------------------------------------------------------# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origin     4    7.2403 1.81008  11.571 0.49617  0.001 ***
# Residuals 47    7.3521 0.15643         0.50383           
# Total     51   14.5924                 1.00000   
#----------------------------------------------------------------# 

# Are communities shaped by site of origin nested within morphology?
permanova.morph.site <- adonis(cg.f.distmat0 ~ morph/origin, strata = cg.f.distmeta0$morph, data = cg.f.distmeta0)
permanova.morph.site
#----------------------------------------------------------------# 
#              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# morph         2    5.8162 2.90810  18.591 0.39858  0.001 ***
# morph:origin  2    1.4241 0.71206   4.552 0.09759  0.001 ***
# Residuals    47    7.3521 0.15643         0.50383           
# Total        51   14.5924                 1.00000
#----------------------------------------------------------------# 
##################################################################

##################################################################
#### Test for site effect within Morphotype A ####
cg.ft.fA <- subset_samples(cg.ft.f0, sample_data(cg.ft.f0)$morph == 'A')
cg.f.distA <- distance(cg.ft.fA, method='bray', type='samples') #make distance object
cg.f.distmatA <- as.matrix(cg.f.distA) #convert to matrix
sample.orderA <- rownames(cg.f.distmatA) #save sampleid order from rownames of distance matrix

#### Filter metadata to match distance matrices ####
cg.f.distmetaA <- as.data.frame(unclass(sample_data(cg.ft.fA))) #get metadata
cg.f.distmetaA %>% group_by(origin) %>% summarise(n()) # get sample sizes per site

# reorder distance matrix to match meta data
mA <- as.character(cg.f.distmetaA$sampleid) #get list of sample ids
cg.f.distmatA <- cg.f.distmatA[mA,mA]

# Are communities shaped by morphology?
permanova.Asite <- adonis(cg.f.distmatA ~ origin , data = cg.f.distmetaA)
permanova.Asite
#----------------------------------------------------------------# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origin     2    1.4241 0.71206  4.0864 0.30077  0.001 ***
# Residuals 19    3.3107 0.17425         0.69923           
# Total     21    4.7349                 1.00000  
#----------------------------------------------------------------# 
