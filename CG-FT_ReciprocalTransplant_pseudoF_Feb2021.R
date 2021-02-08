library(phyloseq)
library(vegan)
library(tidyverse)
library(gridExtra)

#### Read in dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/CG_FT_combined_phyloseq_r1500_Jan2021.RDS")

##################################################################
# Is there a change in how differentiated transplant communities-#
# are over time relative to controls?----------------------------# 

##################################################################
#### Compare WB low->PB transplants to WB low controls ####

# Subset by origin and habitat type for WB 
ft.f.wb <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$origin == "WB low")
ft.f.wb.meta <- as.data.frame(unclass(sample_data(ft.f.wb)))

# Further subset by sample number (day of study)
ft.f.wb0 <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == 0)
ft.f.wb1 <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == 1)
ft.f.wb3 <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == 3)
ft.f.wb5 <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == 5)

# Comparsion T0
ft.f.wb0.dist <- distance(ft.f.wb0, method='bray', type='samples') #make distance object
ft.f.wb0.distmat <- as.matrix(ft.f.wb0.dist)
ft.f.wb.meta0 <- as.data.frame(unclass(sample_data(ft.f.wb0)))
ft.f.wb.meta0 <- ft.f.wb.meta0 %>% mutate_if(is.factor, as.character)
meta0 <- unique(ft.f.wb.meta0$sampleid)
ft.f.wb0.distmat <- ft.f.wb0.distmat[meta0,meta0]
permanova.ft.wb0 <- adonis(ft.f.wb0.distmat ~ treatment, data = ft.f.wb.meta0)
permanova.ft.wb0
#PERMANOVA WB low -> PBtransplant vs WB control at time 0-----#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1   0.38379 0.38379  2.3256 0.13423  0.012 *
# Residuals 15   2.47543 0.16503         0.86577         
# Total     16   2.85922                 1.00000  
#-------------------------------------------------------------#

# Comparsion T1
ft.f.wb1.dist <- distance(ft.f.wb1, method='bray', type='samples') #make distance object
ft.f.wb1.distmat <- as.matrix(ft.f.wb1.dist)
ft.f.wb.meta1 <- as.data.frame(unclass(sample_data(ft.f.wb1)))
ft.f.wb.meta1 <- ft.f.wb.meta1 %>% mutate_if(is.factor, as.character)
meta1 <- unique(ft.f.wb.meta1$sampleid)
ft.f.wb1.distmat <- ft.f.wb1.distmat[meta1,meta1]
permanova.ft.wb1 <- adonis(ft.f.wb1.distmat ~ treatment, data = ft.f.wb.meta1)
permanova.ft.wb1
#PERMANOVA WB low -> PBtransplant vs WB control at time 1------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1    0.4290 0.42900  2.7891 0.21808  0.018 *
# Residuals 10    1.5382 0.15382         0.78192         
# Total     11    1.9672                 1.00000        
#--------------------------------------------------------------#

# Comparsion T3
ft.f.wb3.dist <- distance(ft.f.wb3, method='bray', type='samples') #make distance object
ft.f.wb3.distmat <- as.matrix(ft.f.wb3.dist)
ft.f.wb.meta3 <- as.data.frame(unclass(sample_data(ft.f.wb3)))
ft.f.wb.meta3 <- ft.f.wb.meta3 %>% mutate_if(is.factor, as.character)
meta3 <- unique(ft.f.wb.meta3$sampleid)
ft.f.wb3.distmat <- ft.f.wb3.distmat[meta3,meta3]
permanova.ft.wb3 <- adonis(ft.f.wb3.distmat ~ treatment, data = ft.f.wb.meta3)
permanova.ft.wb3
#PERMANOVA WB low -> PBtransplant vs WB control at time 3-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1    0.4608 0.46084  3.2051 0.14434  0.003 **
# Residuals 19    2.7319 0.14378         0.85566          
# Total     20    3.1927                 1.00000                  
#------------------------------------------------------------#

# Comparsion T5
ft.f.wb5.dist <- distance(ft.f.wb5, method='bray', type='samples') #make distance object
ft.f.wb5.distmat <- as.matrix(ft.f.wb5.dist)
ft.f.wb.meta5 <- as.data.frame(unclass(sample_data(ft.f.wb5)))
ft.f.wb.meta5 <- ft.f.wb.meta5 %>% mutate_if(is.factor, as.character)
meta5 <- unique(ft.f.wb.meta5$sampleid)
ft.f.wb5.distmat <- ft.f.wb5.distmat[meta5,meta5]
permanova.ft.wb5 <- adonis(ft.f.wb5.distmat ~ treatment, data = ft.f.wb.meta5)
permanova.ft.wb5
#PERMANOVA transplant vs control at time 5-------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1    0.4846 0.48460  2.6249 0.14893  0.004 **
# Residuals 15    2.7693 0.18462         0.85107          
# Total     16    3.2539                 1.00000       
#------------------------------------------------------------#

#### Plot WB low->PB transplant vs. WB low control pseudo-F over time ####
# How differentiated are experimental treatments at each time
# compared to controls at origin site?
wb.sample.number <- c(0,1,3,5) # sample numbers 
wbw.pseudo_f <- c(permanova.ft.wb0[["aov.tab"]]$F.Model[1], 
                  permanova.ft.wb1[["aov.tab"]]$F.Model[1],
                  permanova.ft.wb3[["aov.tab"]]$F.Model[1],
                  permanova.ft.wb5[["aov.tab"]]$F.Model[1]) #pseudo-F values

wbw.pseudo_f.df <- as.data.frame(cbind(wb.sample.number, wbw.pseudo_f)) # combine
wbw.pseudo_f.df$permanova.comparison <- "WB low->PB transplant to WB low control"
wbw.pseudo_f.df$transplant <- "WB low->PB transplants"
wbw.pseudo_f.df$compared.to <- "origin site"
colnames(wbw.pseudo_f.df) <- c("sample.number", "pseudo_F", "permanova.comparsion", "transplant", "compared.to")
ggplot(wbw.pseudo_f.df, aes(x = sample.number, y = pseudo_F)) + geom_point() + theme_classic()
##################################################################

##################################################################
#### Compare WB low->PB transplants to PB controls ####
# Subset by treatment and habitat type 
ft.f.wbp <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$treatment %in% c("WB low->PB", "PB->PB"))
ft.f.wbp.meta <- as.data.frame(unclass(sample_data(ft.f.wbp)))

# Further subset by sample number (day of study)
ft.f.wbp0 <- subset_samples(ft.f.wbp, sample_data(ft.f.wbp)$sample.number == 0)
ft.f.wbp1 <- subset_samples(ft.f.wbp, sample_data(ft.f.wbp)$sample.number == 1)
ft.f.wbp3 <- subset_samples(ft.f.wbp, sample_data(ft.f.wbp)$sample.number == 3)
ft.f.wbp5 <- subset_samples(ft.f.wbp, sample_data(ft.f.wbp)$sample.number == 5)

# Comparsion T0
ft.f.wbp0.dist <- distance(ft.f.wbp0, method='bray', type='samples') #make distance object
ft.f.wbp0.distmat <- as.matrix(ft.f.wbp0.dist)
ft.f.wbp.meta0 <- as.data.frame(unclass(sample_data(ft.f.wbp0)))
ft.f.wbp.meta0 <- ft.f.wbp.meta0 %>% mutate_if(is.factor, as.character)
meta0 <- unique(ft.f.wbp.meta0$sampleid)
ft.f.wbp0.distmat <- ft.f.wbp0.distmat[meta0,meta0]
permanova.ft.wbp0 <- adonis(ft.f.wbp0.distmat ~ treatment, data = ft.f.wbp.meta0)
permanova.ft.wbp0
#PERMANOVA WB low -> PB transplant vs PB control at time 0-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1    1.7043 1.70430  11.019 0.47869  0.002 **
# Residuals 12    1.8560 0.15467         0.52131          
# Total     13    3.5603                 1.00000   
#-------------------------------------------------------------#

# Comparsion T1
ft.f.wbp1.dist <- distance(ft.f.wbp1, method='bray', type='samples') #make distance object
ft.f.wbp1.distmat <- as.matrix(ft.f.wbp1.dist)
ft.f.wbp.meta1 <- as.data.frame(unclass(sample_data(ft.f.wbp1)))
ft.f.wbp.meta1 <- ft.f.wbp.meta1 %>% mutate_if(is.factor, as.character)
meta1 <- unique(ft.f.wbp.meta1$sampleid)
ft.f.wbp1.distmat <- ft.f.wbp1.distmat[meta1,meta1]
permanova.ft.wbp1 <- adonis(ft.f.wbp1.distmat ~ treatment, data = ft.f.wbp.meta1)
permanova.ft.wbp1
#PERMANOVA WB low -> PBtransplant vs PB control at time 1------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  1    1.6635 1.66347  11.502 0.46944  0.001 ***
# Residuals 13    1.8800 0.14462         0.53056           
# Total     14    3.5435                 1.00000            
#--------------------------------------------------------------#

# Comparsion T3
ft.f.wbp3.dist <- distance(ft.f.wbp3, method='bray', type='samples') #make distance object
ft.f.wbp3.distmat <- as.matrix(ft.f.wbp3.dist)
ft.f.wbp.meta3 <- as.data.frame(unclass(sample_data(ft.f.wbp3)))
ft.f.wbp.meta3 <- ft.f.wbp.meta3 %>% mutate_if(is.factor, as.character)
meta3 <- unique(ft.f.wbp.meta3$sampleid)
ft.f.wbp3.distmat <- ft.f.wbp3.distmat[meta3,meta3]
permanova.ft.wbp3 <- adonis(ft.f.wbp3.distmat ~ treatment, data = ft.f.wbp.meta3)
permanova.ft.wbp3
#PERMANOVA WB low -> PBtransplant vs PB control at time 3-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1    1.4125 1.41252  9.3702 0.34235  0.002 **
# Residuals 18    2.7134 0.15075         0.65765          
# Total     19    4.1259                 1.00000                  
#------------------------------------------------------------#

# Comparsion T5
ft.f.wbp5.dist <- distance(ft.f.wbp5, method='bray', type='samples') #make distance object
ft.f.wbp5.distmat <- as.matrix(ft.f.wbp5.dist)
ft.f.wbp.meta5 <- as.data.frame(unclass(sample_data(ft.f.wbp5)))
ft.f.wbp.meta5 <- ft.f.wbp.meta5 %>% mutate_if(is.factor, as.character)
meta5 <- unique(ft.f.wbp.meta5$sampleid)
ft.f.wbp5.distmat <- ft.f.wbp5.distmat[meta5,meta5]
permanova.ft.wbp5 <- adonis(ft.f.wbp5.distmat ~ treatment, data = ft.f.wbp.meta5)
permanova.ft.wbp5
#PERMANOVA WB low -> PB transplant vs PB control at time 5------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  1    1.5356 1.53562  8.3295 0.35704  0.001 ***
# Residuals 15    2.7654 0.18436         0.64296           
# Total     16    4.3010                 1.00000         
#------------------------------------------------------------#

#### Plot WB low->PB transplant vs. PB control pseudo-F over time ####
# How differentiated are experimental treatments at each time
# compared to controls at destination site?
wb.sample.number <- c(0,1,3,5) # sample numbers 
wbp.pseudo_f <- c(permanova.ft.wbp0[["aov.tab"]]$F.Model[1], 
                  permanova.ft.wbp1[["aov.tab"]]$F.Model[1],
                  permanova.ft.wbp3[["aov.tab"]]$F.Model[1],
                  permanova.ft.wbp5[["aov.tab"]]$F.Model[1]) #pseudo-F values

wbp.pseudo_f.df <- as.data.frame(cbind(wb.sample.number, wbp.pseudo_f)) # combine
wbp.pseudo_f.df$permanova.comparison <- "WB low->PB transplant to PB control"
wbp.pseudo_f.df$transplant <- "WB low->PB transplants"
wbp.pseudo_f.df$compared.to <- "destination site"
colnames(wbp.pseudo_f.df) <- c("sample.number", "pseudo_F", "permanova.comparsion", "transplant", "compared.to")
ggplot(wbp.pseudo_f.df, aes(x = sample.number, y = pseudo_F)) + geom_point() + theme_classic()
##################################################################


##################################################################
#### Compare PB->WB low transplants to PB controls ####
# Subset by treatment and habitat type 
ft.f.pb <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$origin == "PB")
ft.f.pb.meta <- as.data.frame(unclass(sample_data(ft.f.pb)))
unique(ft.f.pb.meta$date)

# Further subset by sample number (day of study)
ft.f.pb0 <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$date == "15-Jun-18")
ft.f.pb1 <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$date == "16-Jun-18")
ft.f.pb2 <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$date == "17-Jun-18")

# Comparsion T0
ft.f.pb0.dist <- distance(ft.f.pb0, method='bray', type='samples') #make distance object
ft.f.pb0.distmat <- as.matrix(ft.f.pb0.dist)
ft.f.pb.meta0 <- as.data.frame(unclass(sample_data(ft.f.pb0)))
ft.f.pb.meta0 <- ft.f.pb.meta0 %>% mutate_if(is.factor, as.character)
meta0 <- unique(ft.f.pb.meta0$sampleid)
ft.f.pb0.distmat <- ft.f.pb0.distmat[meta0,meta0]
permanova.ft.pb0 <- adonis(ft.f.pb0.distmat ~ treatment, data = ft.f.pb.meta0)
permanova.ft.pb0
#PERMANOVA PB->WB low transplant vs PB control at time 0-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treatment  1   0.12048 0.12048  1.0916 0.21439 0.5333
# Residuals  4   0.44148 0.11037         0.78561       
# Total      5   0.56196                 1.00000    
#-------------------------------------------------------------#

# Comparsion T1
ft.f.pb1.dist <- distance(ft.f.pb1, method='bray', type='samples') #make distance object
ft.f.pb1.distmat <- as.matrix(ft.f.pb1.dist)
ft.f.pb.meta1 <- as.data.frame(unclass(sample_data(ft.f.pb1)))
ft.f.pb.meta1 <- ft.f.pb.meta1 %>% mutate_if(is.factor, as.character)
meta1 <- unique(ft.f.pb.meta1$sampleid)
ft.f.pb1.distmat <- ft.f.pb1.distmat[meta1,meta1]
permanova.ft.pb1 <- adonis(ft.f.pb1.distmat ~ treatment, data = ft.f.pb.meta1)
permanova.ft.pb1
#PERMANOVA PB->WB low transplant vs PB control at time 1-----#
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# treatment  1   0.07609 0.076091 0.67099 0.18278    0.6
# Residuals  3   0.34020 0.113402         0.81722       
# Total      4   0.41630                  1.00000           
#--------------------------------------------------------------#

# Comparsion T2
ft.f.pb2.dist <- distance(ft.f.pb2, method='bray', type='samples') #make distance object
ft.f.pb2.distmat <- as.matrix(ft.f.pb2.dist)
ft.f.pb.meta2 <- as.data.frame(unclass(sample_data(ft.f.pb2)))
ft.f.pb.meta2 <- ft.f.pb.meta2 %>% mutate_if(is.factor, as.character)
meta2 <- unique(ft.f.pb.meta2$sampleid)
ft.f.pb2.distmat <- ft.f.pb2.distmat[meta2,meta2]
permanova.ft.pb2 <- adonis(ft.f.pb2.distmat ~ treatment, data = ft.f.pb.meta2)
permanova.ft.pb2
#PERMANOVA PB->WB low transplant vs PB control at time 2-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1   0.18105 0.18105  1.7112 0.13462  0.045 *
# Residuals 11   1.16383 0.10580         0.86538         
# Total     12   1.34488                 1.00000                 
#------------------------------------------------------------#

#### Plot PB->WB low transplant vs. PB control pseudo-F over time ####
# How differentiated are experimental treatments at each time
# compared to controls at origin site?
pb.sample.number <- c(0,1,2) # sample numbers 
pb.pseudo_f <- c(permanova.ft.pb0[["aov.tab"]]$F.Model[1], 
                  permanova.ft.pb1[["aov.tab"]]$F.Model[1],
                  permanova.ft.pb2[["aov.tab"]]$F.Model[1]) #pseudo-F values

pb.pseudo_f.df <- as.data.frame(cbind(pb.sample.number, pb.pseudo_f)) # combine
pb.pseudo_f.df$permanova.comparison <- "PB->WB low transplant to PB control"
pb.pseudo_f.df$transplant <- "PB->WB low transplants"
pb.pseudo_f.df$compared.to <- "origin site"
colnames(pb.pseudo_f.df) <- c("sample.number", "pseudo_F", "permanova.comparsion", "transplant", "compared.to")
ggplot(pb.pseudo_f.df, aes(x = sample.number, y = pseudo_F)) + geom_point() + theme_classic()
##################################################################

##################################################################
#### Compare PB->WB low transplants to WB low controls ####
# Subset by treatment and habitat type 
ft.f.pbw <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT" & sample_data(cg.ft.fr)$type == "fucus" & sample_data(cg.ft.fr)$treatment %in% c("PB->WB low","WB low->WB low"))
ft.f.pbw.meta <- as.data.frame(unclass(sample_data(ft.f.pbw)))
unique(ft.f.pbw.meta$date)

# ft.f.pbw.filt <- subset_samples(ft.f.pbw,sample_data(ft.f.pbw)$date %in% c("15-Jun-18","16-Jun-18","17-Jun-18"))
# ft.pbw.ord <- ordinate(ft.f.pbw.filt, method="PCoA", distance = "bray")
# plot_ordination(ft.f.pbw.filt, ft.pbw.ord, type = "samples", color = "treatment", shape = "date") + theme_classic()

# Further subset by sample number (day of study)
ft.f.pbw0 <- subset_samples(ft.f.pbw, sample_data(ft.f.pbw)$date == "15-Jun-18")
ft.f.pbw1 <- subset_samples(ft.f.pbw, sample_data(ft.f.pbw)$date == "16-Jun-18")
ft.f.pbw2 <- subset_samples(ft.f.pbw, sample_data(ft.f.pbw)$date == "17-Jun-18")

# Comparsion T0
ft.f.pbw0.dist <- distance(ft.f.pbw0, method='bray', type='samples') #make distance object
ft.f.pbw0.distmat <- as.matrix(ft.f.pbw0.dist)
ft.f.pbw.meta0 <- as.data.frame(unclass(sample_data(ft.f.pbw0)))
ft.f.pbw.meta0 <- ft.f.pbw.meta0 %>% mutate_if(is.factor, as.character)
meta0 <- unique(ft.f.pbw.meta0$sampleid)
ft.f.pbw0.distmat <- ft.f.pbw0.distmat[meta0,meta0]
permanova.ft.pbw0 <- adonis(ft.f.pbw0.distmat ~ treatment, data = ft.f.pbw.meta0)
permanova.ft.pbw0
#PERMANOVA PB->WB low transplant vs WB low control at time 0-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1   0.97080 0.97080  10.553 0.67853  0.046 *
# Residuals  5   0.45995 0.09199         0.32147         
# Total      6   1.43075                 1.00000    
#-------------------------------------------------------------#

# Comparsion T1
ft.f.pbw1.dist <- distance(ft.f.pbw1, method='bray', type='samples') #make distance object
ft.f.pbw1.distmat <- as.matrix(ft.f.pbw1.dist)
ft.f.pbw.meta1 <- as.data.frame(unclass(sample_data(ft.f.pbw1)))
ft.f.pbw.meta1 <- ft.f.pbw.meta1 %>% mutate_if(is.factor, as.character)
meta1 <- unique(ft.f.pbw.meta1$sampleid)
ft.f.pbw1.distmat <- ft.f.pbw1.distmat[meta1,meta1]
permanova.ft.pbw1 <- adonis(ft.f.pbw1.distmat ~ treatment, data = ft.f.pbw.meta1)
permanova.ft.pbw1
#PERMANOVA PB->WB low transplant vs WB low control at time 1-----#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treatment  1   0.71835 0.71835  4.7958 0.61518    0.1
# Residuals  3   0.44936 0.14979         0.38482       
# Total      4   1.16772                 1.00000           
#--------------------------------------------------------------#

# Comparsion T2
ft.f.pbw2.dist <- distance(ft.f.pbw2, method='bray', type='samples') #make distance object
ft.f.pbw2.distmat <- as.matrix(ft.f.pbw2.dist)
ft.f.pbw.meta2 <- as.data.frame(unclass(sample_data(ft.f.pbw2)))
ft.f.pbw.meta2 <- ft.f.pbw.meta2 %>% mutate_if(is.factor, as.character)
meta2 <- unique(ft.f.pbw.meta2$sampleid)
ft.f.pbw2.distmat <- ft.f.pbw2.distmat[meta2,meta2]
permanova.ft.pbw2 <- adonis(ft.f.pbw2.distmat ~ treatment, data = ft.f.pbw.meta2)
permanova.ft.pbw2
#PERMANOVA PB->WB low transplant vs WB low control at time 2-----#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1    1.5569 1.55689  14.666 0.57142  0.002 **
# Residuals 11    1.1677 0.10615         0.42858          
# Total     12    2.7246                 1.00000                    
#------------------------------------------------------------#

#### Plot PB->WB low transplant vs. WB low control pseudo-F over time ####
# How differentiated are experimental treatments at each time
# compared to controls at origin site?
pbw.sample.number <- c(0,1,2) # sample numbers 
pbw.pseudo_f <- c(permanova.ft.pbw0[["aov.tab"]]$F.Model[1], 
                 permanova.ft.pbw1[["aov.tab"]]$F.Model[1],
                 permanova.ft.pbw2[["aov.tab"]]$F.Model[1]) #pseudo-F values

pbw.pseudo_f.df <- as.data.frame(cbind(pbw.sample.number, pbw.pseudo_f)) # combine
pbw.pseudo_f.df$permanova.comparison <- "PB->WB low transplant to WB low control"
pbw.pseudo_f.df$transplant <- "PB->WB low transplants"
pbw.pseudo_f.df$compared.to <- "destination site"
colnames(pbw.pseudo_f.df) <- c("sample.number", "pseudo_F", "permanova.comparsion", "transplant", "compared.to")
ggplot(pbw.pseudo_f.df, aes(x = sample.number, y = pseudo_F)) + geom_point() + theme_classic()
##################################################################

##################################################################
#### Combine results of pseudo-F comparisions--------------- ####
all.pseudo_F.df <- rbind(wbw.pseudo_f.df, wbp.pseudo_f.df, pbw.pseudo_f.df, pb.pseudo_f.df)

# Plot pseudo-F values over time
pdf(file="~/Desktop/CG-FT_revised_figures/CG-FT_pseudo_F_SI.pdf", width = 6, height = 3)
ggplot(all.pseudo_F.df, aes(x =sample.number, y = pseudo_F)) + 
  theme_classic(base_size = 11) +
  geom_point(aes(colour = factor(compared.to)), size =2) + 
  facet_grid(~transplant)+ 
    scale_color_manual(values = c("grey0","grey63")) +
  xlab("Sample day") +
  ylab("pseudo-F value") +
  labs(color="Compared to controls at") 
dev.off()

# Check for significant change in pseudo-F value over time
all.pseudo_F.df %>% group_by(compared.to, transplant) %>% 
  group_modify(~ broom::tidy(aov(pseudo_F ~ sample.number, data = .x)))
##################################################################

