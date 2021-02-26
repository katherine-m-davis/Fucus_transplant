#### Fucus reciprocal transplant: change over time in transplant vs control samples #####
# Do control microbiota change over time?
# Do transplant microbiota over time?
# Do transplant microbiota become different from site of origin controls?
# Do transplant microbiota become similar to destination site controls?

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(gridExtra)

#### Read in dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/CG_FT_combined_phyloseq_r1500_Jan2021.RDS")

#### Subset ####
ft.ps <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT")
ft.ps.meta <- as.data.frame(unclass(sample_data(ft.ps)))
ft.ps.meta %>% group_by(type) %>% summarise(n())

##### WB low -> PB compared to WB low ->WB low AND PB -> PB controls ####
ft.wp <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$treatment %in% c('WB low->PB', 'WB low->WB low', 'PB->PB'))
ft.wp.meta <- as.data.frame(unclass(sample_data(ft.wp)))

# Subset to only Fucus samples
ft.wpf <- subset_samples(ft.wp, sample_data(ft.wp)$type == 'fucus')
ft.wp.meta %>% group_by(type, treatment, sample.number) %>% summarise(n()) #Summarise sample numbers per group
pcoa.ft.wpf <- ordinate(ft.wpf, method = "PCoA", distance = "bray") #Ordinate


#### Plot PCoA with WB low->PB and both control sites ####

# Day shapes and trmt colors
ft.day.shapes <- c(16,1,2,6,5,8)
#ft.trmt.colors <- c("#E69F00","#463806", "#7ACE54", "darkblue")
ft.trmt.colors <- c("#E69F00", "#7ACE54", "darkblue")
# Make and save plot
sample_data(ft.wpf)$sample.number <- as.factor(sample_data(ft.wpf)$sample.number)
wfpp1 <- plot_ordination(ft.wpf, pcoa.ft.wpf, type = "samples", color="treatment", shape ="sample.number") +
  theme_classic(base_size = 12) +
  geom_point(size=3) +
  scale_shape_manual(values= ft.day.shapes) +
  scale_color_manual(values = ft.trmt.colors) +
  labs(color='Treatment', shape = 'Sample day') 
wfpp1$layers <- wfpp1$layers[-1]
wfpp1
##################################################################

##################################################################
#### Is there change in control samples over time?----------------# 

#### WB low controls over time #### 
ft.wb <- subset_samples(ft, sample_data(ft)$treatment =='WB low->WB low')
ft.wb.meta <- as.data.frame(unclass(sample_data(ft.wb))) #Get metadata
ft.wb.meta %>% group_by(type, treatment, sample.number) %>% summarise(n()) # samples per timepoint and habitat type

# Fucus
ft.wb.f <- subset_samples(ft.wb, sample_data(ft.wb)$type == "fucus")

# Distance matrices
ft.dist.wbf <- distance(ft.wb.f, method='bray', type='samples') #make distance object
ft.distmat.wbf <- as.matrix(ft.dist.wbf) #convert to matrix

# reorder distance matrix to match meta data
ft.distmeta.wbf <- as.data.frame(unclass(sample_data(ft.wb.f)))
wb <- as.character(ft.distmeta.wbf$sampleid) #order by sampleid
ft.distmat.wbf <- ft.distmat.wbf[wb,wb] #order by sampleid
ft.distmeta.wbf %>% group_by(sample.number) %>% summarise(n()) #Get summary of sample numbers per group

wb #check that rows of matrix match
rownames(ft.distmat.wbf)

# PERMANOVA for change over time in  in WB low -> WB low fucus controls (No significant change)
permanova.ft.wbf <- adonis(ft.distmat.wbf ~  sample.number,  data = ft.distmeta.wbf)
permanova.ft.wbf
#----------------------------------------------------------------# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# sample.number  1   0.22115 0.22115   1.535 0.07128  0.094 .
# Residuals     20   2.88151 0.14408         0.92872         
# Total         21   3.10266                 1.00000  
#----------------------------------------------------------------# 
##################################################################

#### PB controls over time ####
ft.pb <- subset_samples(ft, sample_data(ft)$treatment =='PB->PB')
ft.pb.meta <- as.data.frame(unclass(sample_data(ft.pb)))
ft.pb.meta %>% group_by(type, treatment, sample.number) %>% summarise(n()) # samples per timepoint and habitat type

# Fucus
ft.pb.f <- subset_samples(ft.pb, sample_data(ft.pb)$type == "fucus")

# Distance matrices
ft.dist.pbf <- distance(ft.pb.f, method='bray', type='samples') #make distance object
ft.distmat.pbf <- as.matrix(ft.dist.pbf) #convert to matrix

# reorder distance matrix to match meta data
ft.distmeta.pbf <- as.data.frame(unclass(sample_data(ft.pb.f)))
pb <- as.character(ft.distmeta.pbf$sampleid)
ft.distmat.pbf <- ft.distmat.pbf[pb,pb]
ft.distmeta.pbf %>% group_by(sample.number) %>% summarise(n()) #Get summary of sample numbers per group

####PERMANOVA for change over time in PB->PB fucus controls (No significant change)
permanova.ft.pbf <- adonis(ft.distmat.pbf ~  sample.number,  data = ft.distmeta.pbf)
permanova.ft.pbf
#----------------------------------------------------------------# 
#               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# sample.number  1   0.07074 0.070739 0.59644 0.03044   0.92
# Residuals     19   2.25343 0.118602         0.96956       
# Total         20   2.32417                  1.00000 
#----------------------------------------------------------------# 
##################################################################

##################################################################
#### Is there change in experimental samples over time?----------# 

#### PB -> WB low transplants over time ####
ft.pw <- subset_samples(ft, sample_data(ft)$treatment =='PB->WB low')
ft.pw.meta <- as.data.frame(unclass(sample_data(ft.pw)))
ft.pw.meta %>% group_by(type, treatment, sample.number) %>% summarise(n()) # samples per timepoint and habitat type

# Fucus
ft.pw.f <- subset_samples(ft.pw, sample_data(ft.pw)$type == "fucus")

# Distance matrices
ft.dist.pwf <- distance(ft.pw.f, method='bray', type='samples') #make distance object
ft.distmat.pwf <- as.matrix(ft.dist.pwf) #convert to matrix

# reorder distance matrix to match meta data
ft.distmeta.pwf <- as.data.frame(unclass(sample_data(ft.pw.f)))
pw <- as.character(ft.distmeta.pwf$sampleid)
ft.distmat.pwf <- ft.distmat.pwf[pw,pw]
ft.distmeta.pwf %>% group_by(sample.number) %>% summarise(n()) #Get summary of sample numbers per group

#### PERMANOVA for change over time in PB->WB low fucus transplants (No significant change)
permanova.ft.pwf <- adonis(ft.distmat.pwf ~  sample.number,  data = ft.distmeta.pwf)
permanova.ft.pwf
#----------------------------------------------------------------# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# sample.number  1   0.14507 0.14507  1.3908 0.11224  0.164
# Residuals     11   1.14736 0.10430         0.88776       
# Total         12   1.29242                 1.00000 
#----------------------------------------------------------------# 
##################################################################


#### WB low -> PB transplants over time ####
ft.wp <- subset_samples(ft, sample_data(ft)$treatment =='WB low->PB')
ft.wp.meta <- as.data.frame(unclass(sample_data(ft.wp)))
ft.wp.meta %>% group_by(type, treatment, sample.number) %>% summarise(n()) # samples per timepoint and habitat type

# Fucus
ft.wp.f <- subset_samples(ft.wp, sample_data(ft.wp)$type == "fucus")

# Distance matrices
ft.dist.wpf <- distance(ft.wp.f, method='bray', type='samples') #make distance object
ft.distmat.wpf <- as.matrix(ft.dist.wpf) #convert to matrix

# reorder distance matrix to match meta data
ft.distmeta.wpf <- as.data.frame(unclass(sample_data(ft.wp.f)))
wp <- as.character(ft.distmeta.wpf$sampleid)
ft.distmat.wpf <- ft.distmat.pwf[wp,wp]
ft.distmeta.wpf %>% group_by(sample.number) %>% summarise(n()) #Get summary of sample numbers per group

#### PERMANOVA for change over time in WB low-> PB fucus transplants (Significant change)
permanova.ft.wpf <- adonis(ft.distmat.wpf ~  sample.number,  data = ft.distmeta.wpf)
permanova.ft.wpf
#----------------------------------------------------------------# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sample.number  1    0.9295 0.92948  5.1018 0.10183  0.001 ***
# Residuals     45    8.1984 0.18219         0.89817           
# Total         46    9.1279                 1.00000  
#----------------------------------------------------------------# 
##################################################################


##################################################################
# Do Rock transplants change overtime?---------------------------# 

#### PB -> WB low rocks over time ####
ft.pw.r <- subset_samples(ft.pw, sample_data(ft.pw)$type == "rock")

# Distance matrices
ft.dist.pwr <- distance(ft.pw.r, method='bray', type='samples') #make distance object
ft.distmat.pwr <- as.matrix(ft.dist.pwr) #convert to matrix

# reorder distance matrix to match meta data
ft.distmeta.pwr <- as.data.frame(unclass(sample_data(ft.pw.r)))
pwr <- as.character(ft.distmeta.pwr$sampleid)
ft.distmat.pwr <- ft.distmat.pwr[pwr,pwr]
ft.distmeta.pwr %>% group_by(sample.number) %>% summarise(n()) #Get summary of sample numbers per group

#### PERMANOVA for change over time in PB->WB low rock transplants (No significant change)
permanova.ft.pwr <- adonis(ft.distmat.pwr ~  sample.number,  data = ft.distmeta.pwr)
permanova.ft.pwr
#----------------------------------------------------------------# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# sample.number  1    0.2744 0.27440 0.91263 0.11534  0.543
# Residuals      7    2.1047 0.30067         0.88466       
# Total          8    2.3791                 1.00000       
#----------------------------------------------------------------# 
##################################################################

##################################################################
#### WB low -> PB rocks over time ####
ft.wp.r <- subset_samples(ft.wp, sample_data(ft.wp)$type == "rock")

# Distance matrices
ft.dist.wpr <- distance(ft.wp.r, method='bray', type='samples') #make distance object
ft.distmat.wpr <- as.matrix(ft.dist.wpr) #convert to matrix

# reorder distance matrix to match meta data
ft.distmeta.wpr <- as.data.frame(unclass(sample_data(ft.wp.r)))
wpr <- as.character(ft.distmeta.wpr$sampleid)
ft.distmat.wpr <- ft.distmat.pwr[wpr,wpr]
ft.distmeta.wpr %>% group_by(sample.number) %>% summarise(n()) #Get summary of sample numbers per group

#### PERMANOVA for change over time in WB low -> PB rock transplants (Significant change)
permanova.ft.wpr <- adonis(ft.distmat.wpr ~  sample.number,  data = ft.distmeta.wpr)
permanova.ft.wpr
#----------------------------------------------------------------# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# sample.number  1    0.6125 0.61249  2.4145 0.15664  0.016 *
# Residuals     13    3.2977 0.25367         0.84336         
# Total         14    3.9102                 1.00000
#----------------------------------------------------------------# 
##################################################################


