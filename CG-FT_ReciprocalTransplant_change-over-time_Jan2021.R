#### Reciprocal Transplant differentiation between treatment and controls over time ####
# Calculate pseudo-F from PERMANOVA of transplant vs control on each sampling day
# Compare pseudo-F value change over time

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(viridis)
library(gridExtra)

#### Read in dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/Combined/CG_FT_rock_fucus_combined_phyloseq_r1500.RDS")

# Subset to only Fucus samples
ft.fr <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == 'FT')
ft.fr.meta <- as.data.frame(unclass(sample_data(ft.fr)))

#### PERMANOVA for each time point ####
sample_data(ft.fr)$trt.time <- paste(sample_data(ft.fr)$treatment, sample_data(ft.fr)$sample.number, sep="_")
sample_data(ft.fr)$trt.time

# Separate by Fucus and Rock 
ft.f.wb <- subset_samples(ft.fr, sample_data(ft.fr)$origin == "WB low" & sample_data(ft.fr)$type == "fucus")
ft.r.wb <- subset_samples(ft.fr, sample_data(ft.fr)$origin == "WB low" & sample_data(ft.fr)$type == "rock")

#### Subset to controls and transplants at each timepoint ####
#### Compare tranpslants at T0 to WB low controls 
ft.f.wb.0c <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == 0)
wb.f.0c.meta <- data.frame(sample_data(ft.f.wb.0c))
wb.f.0c.meta %>% group_by(location) %>% summarise(n())
pcoa.ft.wb0c<- ordinate(ft.f.wb.0c, method = "PCoA", distance = "bray")

sample_data(ft.f.wb.0c)$sample.number <- as.character(sample_data(ft.f.wb.0c)$sample.number)
plot_ordination(ft.f.wb.0c, pcoa.ft.wb0c, type = "samples", color="treatment", shape ="sample.number") +
  theme_classic(base_size = 12) +
  geom_point(size=3) +
  scale_shape_manual(values= ft.day.shapes) +
  scale_color_manual(values = ft.trmt.colors) +
  labs(color='Treatment', shape = 'Sample day') 



#### Make distance object T0
wb.f.0c.dist <- distance(ft.f.wb.0c, method='bray', type='samples') #make distance object
wb.f.0c.distmat <- as.matrix(wb.f.0c.dist) #convert to matrix
wb0c <- wb.f.0c.meta$sampleid
wb.f.0c.distmat <- wb.f.0c.distmat[wb0c,wb0c] # match metadata order

# PERMANOVA T0 transplants compared to controls
pa.f.wb.0c <- adonis(wb.f.0c.distmat ~ location, data = wb.f.0c.meta)
pa.f.wb.0c
# ----------------------------------------------------------------#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# trt.time   1   0.38379 0.38379  2.3256 0.13423  0.014 *
# Residuals 15   2.47543 0.16503         0.86577         
# Total     16   2.85922                 1.00000 
# ----------------------------------------------------------------#
####################################################################

####################################################################
#### Compare tranpslants at T1 to WB low controls 
ft.f.wb.1c <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$individual == "neighbor"|c(sample_data(ft.f.wb)$individual == "transplant" & sample_data(ft.f.wb)$sample.number == 1))
wb.f.1c.meta <- data.frame(sample_data(ft.f.wb.1c))
wb.f.1c.meta %>% group_by(location) %>% summarise(n())

#### Make distance object T1
wb.f.1c.dist <- distance(ft.f.wb.1c, method='bray', type='samples') #make distance object
wb.f.1c.distmat <- as.matrix(wb.f.1c.dist) #convert to matrix
wb1c <- wb.f.1c.meta$sampleid
wb.f.1c.distmat <- wb.f.1c.distmat[wb1c,wb1c] # match metadata order

pcoa.ft.wb1c<- ordinate(ft.f.wb.1c, method = "PCoA", distance = "bray")
sample_data(ft.f.wb.1c)$sample.number <- as.character(sample_data(ft.f.wb.1c)$sample.number)
plot_ordination(ft.f.wb.1c, pcoa.ft.wb1c, type = "samples", color="treatment", shape ="sample.number") +
  theme_classic(base_size = 12) +
  geom_point(size=3) +
  scale_shape_manual(values= ft.day.shapes) +
  scale_color_manual(values = ft.trmt.colors) +
  labs(color='Treatment', shape = 'Sample day') 


# PERMANOVA T1 transplants compared to controls
pa.f.wb.1c <- adonis(wb.f.1c.distmat ~ location, data = wb.f.1c.meta)
pa.f.wb.1c
# ----------------------------------------------------------------#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# trt.time   1   0.38379 0.38379  2.3256 0.13423  0.014 *
# Residuals 15   2.47543 0.16503         0.86577         
# Total     16   2.85922                 1.00000 
# ----------------------------------------------------------------#
####################################################################

####################################################################
#### Compare tranpslants at T3 to WB low controls 
ft.f.wb.3c <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$individual == "neighbor"|c(sample_data(ft.f.wb)$individual == "transplant" & sample_data(ft.f.wb)$sample.number == 3))
wb.f.3c.meta <- data.frame(sample_data(ft.f.wb.3c))
wb.f.3c.meta %>% group_by(location) %>% summarise(n())

#### Make distance object T3
wb.f.3c.dist <- distance(ft.f.wb.3c, method='bray', type='samples') #make distance object
wb.f.3c.distmat <- as.matrix(wb.f.3c.dist) #convert to matrix
wb3c <- wb.f.3c.meta$sampleid
wb.f.3c.distmat <- wb.f.3c.distmat[wb3c,wb3c] # match metadata order

# PERMANOVA T3 transplants compared to controls
pa.f.wb.3c <- adonis(wb.f.3c.distmat ~ location, data = wb.f.3c.meta)
pa.f.wb.3c
#----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# location   1    0.7394 0.73939   5.048 0.14004  0.001 ***
# Residuals 31    4.5406 0.14647         0.85996           
# Total     32    5.2800                 1.00000 
# ----------------------------------------------------------------#
####################################################################

####################################################################
#### Compare tranpslants at T5 to WB low controls 
ft.f.wb.5c <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$individual == "neighbor"|c(sample_data(ft.f.wb)$individual == "transplant" & sample_data(ft.f.wb)$sample.number == 5))
wb.f.5c.meta <- data.frame(sample_data(ft.f.wb.5c))
wb.f.5c.meta %>% group_by(location) %>% summarise(n())

#### Make distance object T5
wb.f.5c.dist <- distance(ft.f.wb.5c, method='bray', type='samples') #make distance object
wb.f.5c.distmat <- as.matrix(wb.f.5c.dist) #convert to matrix
wb5c <- wb.f.5c.meta$sampleid
wb.f.5c.distmat <- wb.f.5c.distmat[wb5c,wb5c] # match metadata order

# PERMANOVA T3 transplants compared to controls
pa.f.wb.5c <- adonis(wb.f.5c.distmat ~ location, data = wb.f.5c.meta)
pa.f.wb.5c
#----------------------------------------------------------------#
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# location   1    0.9991 0.99913  5.9748 0.18119  0.001 ***
# Residuals 27    4.5151 0.16723         0.81881           
# Total     28    5.5142                 1.00000 
# ----------------------------------------------------------------#
####################################################################

####################################################################

#### Compare to original time 0 WB controls ####
ft.f.wb.0.all <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$trt.time %in% c("WB->WB_0", "WB->PB_0"))
wb.f.0.all.dist <- distance(ft.f.wb.0.all, method='bray', type='samples') #make distance object
ft.f.wb.0.all.meta <- data.frame(sample_data(ft.f.wb.0.all))
ft.f.wb.0.all.meta$trt.time

pa.ft.f.wb.0.all <- adonis(wb.f.0.all.dist ~ trt.time, data = ft.f.wb.0.all.meta)
pa.ft.f.wb.0.all


#### Compare to original time 1 WB controls and trt ####
ft.f.wb.1.all <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$trt.time %in% c("WB low->WB low_0", "WB low->PB_1", "WB low->WB low_1"))
wb.f.1.all.dist <- distance(ft.f.wb.1.all, method='bray', type='samples') #make distance object
ft.f.wb.1.all.meta <- data.frame(sample_data(ft.f.wb.1.all))
ft.f.wb.1.all.meta$trt.time

pa.ft.f.wb.1.all <- adonis(wb.f.1.all.dist ~ sample.number, data = ft.f.wb.1.all.meta)
pa.ft.f.wb.1.all

# "WB->PB_0", "WB->PB_1"
# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# location   1   0.40127 0.40127  2.4302 0.13942  0.012 *
# Residuals 15   2.47678 0.16512         0.86058         
# Total     16   2.87805                 1.00000    
# ----------------------------------------------------------------#

# "WB->PB_0", "WB->PB_0", "WB->PB_1"
# ----------------------------------------------------------------#
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# sample.number  1    0.2952 0.29515  1.7134 0.06414  0.066 .
# Residuals     25    4.3066 0.17227         0.93586         
# Total         26    4.6018                 1.00000     
# ----------------------------------------------------------------#


#### Time 3 ####

ft.f.wb.3.0 <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == "3" & sample_data(ft.f.wb)$treatment == "WB->PB"|sample_data(ft.f.wb)$sample.number == "0" & sample_data(ft.f.wb)$treatment == "WB->WB")
wb.f.3.0.dist <- distance(ft.f.wb.3.0, method='bray', type='samples') #make distance object
wb.f.3.0.meta <- data.frame(sample_data(ft.f.wb.3.0))

pa.f.wb.3.0 <- adonis(wb.f.3.0.dist ~ treatment, data = wb.f.3.0.meta)
pa.f.wb.3.0

# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  1    0.5675 0.56748   3.703 0.14407  0.001 ***
# Residuals 22    3.3714 0.15325         0.85593           
# Total     23    3.9389                 1.00000    
# ----------------------------------------------------------------#

ft.f.wb.3.all <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$trt.time %in% c("WB->WB_0", "WB->PB_0", "WB->PB_3"))
wb.f.3.all.dist <- distance(ft.f.wb.3.all, method='bray', type='samples') #make distance object
ft.f.wb.3.all.meta <- data.frame(sample_data(ft.f.wb.3.all))
ft.f.wb.3.all.meta$trt.time

pa.ft.f.wb.3.all <- adonis(wb.f.3.all.dist ~ sample.number, data = ft.f.wb.3.all.meta)
pa.ft.f.wb.3.all

# ----------------------------------------------------------------#
#               Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample.number  1    0.8562 0.85620  5.0634 0.1404  0.001 ***
# Residuals     31    5.2419 0.16909         0.8596           
# Total         32    6.0981                 1.0000       
# ----------------------------------------------------------------#



ft.f.wb.5.0 <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$sample.number == "5" & sample_data(ft.f.wb)$treatment == "WB->PB"|sample_data(ft.f.wb)$sample.number == "0" & sample_data(ft.f.wb)$treatment == "WB->WB")
wb.f.5.0.dist <- distance(ft.f.wb.5.0, method='bray', type='samples') #make distance object
wb.f.5.0.meta <- data.frame(sample_data(ft.f.wb.5.0))

pa.f.wb.5.0 <- adonis(wb.f.5.0.dist ~ treatment, data = wb.f.5.0.meta)
pa.f.wb.5.0

# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  1    0.7552 0.75522  3.9471 0.16483  0.001 ***
# Residuals 20    3.8267 0.19133         0.83517           
# Total     21    4.5819                 1.00000      
# ----------------------------------------------------------------#


ft.f.wb.5.all <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$trt.time %in% c("WB->WB_0", "WB->PB_0", "WB->PB_5"))
wb.f.5.all.dist <- distance(ft.f.wb.3.all, method='bray', type='samples') #make distance object
ft.f.wb.5.all.meta <- data.frame(sample_data(ft.f.wb.5.all))
ft.f.wb.5.all.meta$trt.time

pa.ft.f.wb.5.all <- adonis(wb.f.5.all.dist ~ sample.number, data = ft.f.wb.5.all.meta)
pa.ft.f.wb.5.all

# ----------------------------------------------------------------#
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sample.number  1    1.0780 1.07804  5.4875 0.15912  0.001 ***
# Residuals     29    5.6972 0.19645         0.84088           
# Total         30    6.7752                 1.00000  
# ----------------------------------------------------------------#




f.stat <- c(0, 2.4302, 1.7134, 3.703, 5.0634,3.9471, 5.4875)
s.num <- c(0, 1,1, 3,3, 5,5)
wb.wb <- cbind(f.stat,s.num)
wb.wb <- data.frame(wb.wb)
wb.wb$transplant <- "WB->PB"
wb.wb$comparison <- "WB"
wb.wb$type <- "F.distichus"
str(wb.wb)



wb.lm <- lm(f.stat ~ s.num)
summary(wb.lm)
plot(f.stat ~ s.num)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   1.0063     0.6763   1.488   0.1969  
# s.num         0.8500     0.2139   3.975   0.0106 *
# Residual standard error: 1.041 on 5 degrees of freedom
# Multiple R-squared:  0.7596,	Adjusted R-squared:  0.7115 
# F-statistic:  15.8 on 1 and 5 DF,  p-value: 0.01058


#### WB vs PB ####
ft.f.pb <- subset_samples(ft.fr, sample_data(ft.fr)$origin == "PB" & sample_data(ft.fr)$type == "fucus")
ft.r.pb <- subset_samples(ft.fr, sample_data(ft.fr)$origin == "PB" & sample_data(ft.fr)$type == "rock")


#### PERMANOVA by time PB fucus trans ####
ft.f.pb.trans <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$treatment == "PB->WB"|sample_data(ft.f.pb)$trt.time == "PB->PB_0")
pb.f.tran.dist <- distance(ft.f.pb.trans, method='bray', type='samples') #make distance object
ft.f.pb.trans.meta <- data.frame(sample_data(ft.f.pb.trans))
ft.f.pb.trans.meta$sample.number
pb.trans <-  adonis(pb.f.tran.dist ~ sample.number, data = ft.f.pb.trans.meta)
pb.trans

# ----------------------------------------------------------------#
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# sample.number  2   0.30146 0.15073  1.4461 0.18199  0.105
# Residuals     13   1.35502 0.10423         0.81801       
# Total         15   1.65648                 1.00000 
# ----------------------------------------------------------------#


#### PERMANOVA by time WB fucus trans ####
ft.f.wb.trans <- subset_samples(ft.f.wb, sample_data(ft.f.wb)$treatment == "WB->PB"|sample_data(ft.f.wb)$trt.time == "WB->WB_0")
wb.f.tran.dist <- distance(ft.f.wb.trans, method='bray', type='samples') #make distance object
ft.f.wb.trans.meta <- data.frame(sample_data(ft.f.wb.trans))
ft.f.wb.trans.meta$sample.number
wb.trans <-  adonis(wb.f.tran.dist ~ sample.number, data = ft.f.wb.trans.meta)
wb.trans

# ----------------------------------------------------------------#
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sample.number  3    1.6938 0.56461  3.1534 0.15146  0.001 ***
# Residuals     53    9.4896 0.17905         0.84854           
# Total         56   11.1835                 1.00000 
# ----------------------------------------------------------------#


#### PERMANOVA by time PB rock trans ####
ft.r.pb <- subset_samples(ft.fr, sample_data(ft.fr)$origin == "PB" & sample_data(ft.fr)$type == "rock")
ft.r.pb.trans <- subset_samples(ft.r.pb, sample_data(ft.r.pb)$treatment == "PB->WB"|sample_data(ft.r.pb)$trt.time == "PB->PB_0")
pb.r.tran.dist <- distance(ft.r.pb.trans, method='bray', type='samples') #make distance object
ft.r.pb.trans.meta <- data.frame(sample_data(ft.r.pb.trans))
ft.r.pb.trans.meta$sample.number
pb.r.trans <-  adonis(pb.r.tran.dist ~ sample.number, data = ft.r.pb.trans.meta)
pb.r.trans
# ----------------------------------------------------------------#
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# sample.number  2   0.53827 0.26914 0.86342 0.25671   0.64
# Residuals      5   1.55855 0.31171         0.74329       
# Total          7   2.09682                 1.00000   
# ----------------------------------------------------------------#

ft.r.wb <- subset_samples(ft.fr, sample_data(ft.fr)$origin == "WB" & sample_data(ft.fr)$type == "rock")
ft.r.wb.trans <- subset_samples(ft.r.wb, sample_data(ft.r.wb)$treatment == "WB->PB"|sample_data(ft.r.wb)$trt.time == "WB->WB_0")
wb.r.tran.dist <- distance(ft.r.wb.trans, method='bray', type='samples') #make distance object
ft.r.wb.trans.meta <- data.frame(sample_data(ft.r.wb.trans))
ft.r.wb.trans.meta$sample.number
wb.r.trans <-  adonis(wb.r.tran.dist ~ sample.number, data = ft.r.wb.trans.meta)
wb.r.trans









#### Compare to original time 0 PB controls ####
ft.f.pb.1.all <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$trt.time %in% c("PB->WB_1", "PB->PB_3"))
pb.f.1.all.dist <- distance(ft.f.pb.1.all, method='bray', type='samples') #make distance object
ft.f.pb.1.all.meta <- data.frame(sample_data(ft.f.pb.1.all))
ft.f.pb.1.all.meta$trt.time

pa.ft.f.pb.1.all <- adonis(pb.f.1.all.dist ~ trt.time, data = ft.f.pb.1.all.meta)
pa.ft.f.pb.1.all

# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# trt.time   1   0.17123 0.17123  1.3189 0.14153  0.213
# Residuals  8   1.03859 0.12982         0.85847       
# Total      9   1.20981                 1.00000 
# ----------------------------------------------------------------#

ft.f.pb.1.all2 <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$trt.time %in% c("PB->WB_1", "PB->PB_3", "PB->PB_4"))
pb.f.1.all.dist2 <- distance(ft.f.pb.1.all2, method='bray', type='samples') #make distance object
ft.f.pb.1.all.meta2 <- data.frame(sample_data(ft.f.pb.1.all2))
ft.f.pb.1.all.meta2$location

pa.ft.f.pb.1.all2 <- adonis(pb.f.1.all.dist2 ~ location, data = ft.f.pb.1.all.meta2)
pa.ft.f.pb.1.all2

# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# location   1   0.13021 0.13021  1.0328 0.09361  0.401
# Residuals 10   1.26074 0.12607         0.90639       
# Total     11   1.39095                 1.00000  
# ----------------------------------------------------------------#



ft.f.pb.2.all <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$trt.time %in% c("PB->WB_2", "PB->PB_3"))
pb.f.2.all.dist <- distance(ft.f.pb.2.all, method='bray', type='samples') #make distance object
ft.f.pb.2.all.meta <- data.frame(sample_data(ft.f.pb.2.all))
ft.f.pb.2.all.meta$trt.time

pa.ft.f.wb.2.all <- adonis(pb.f.2.all.dist ~ trt.time, data = ft.f.pb.2.all.meta)
pa.ft.f.wb.2.all

# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# trt.time   1   0.21233 0.21233  1.8918 0.12704  0.012 *
# Residuals 13   1.45909 0.11224         0.87296         
# Total     14   1.67142                 1.00000         
# ----------------------------------------------------------------#


ft.f.pb.2.all2 <- subset_samples(ft.f.pb, sample_data(ft.f.pb)$trt.time %in% c("PB->WB_2", "PB->PB_3", "PB->PB_5"))
pb.f.2.all.dist2 <- distance(ft.f.pb.2.all2, method='bray', type='samples') #make distance object
ft.f.pb.2.all.meta2 <- data.frame(sample_data(ft.f.pb.2.all2))
ft.f.pb.2.all.meta2$location

pa.ft.f.wb.2.all2 <- adonis(pb.f.2.all.dist2 ~ location, data = ft.f.pb.2.all.meta2)
pa.ft.f.wb.2.all2

# ----------------------------------------------------------------#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# location   1   0.19546 0.19546  1.6427 0.08363  0.056 .
# Residuals 18   2.14177 0.11899         0.91637         
# Total     19   2.33723                 1.00000 
# ----------------------------------------------------------------#

f.stat <- c(0, 1.3189, 1.0328, 1.6427,1.8918)
s.num <- c(0, 1,1, 2,2)
pb.pb <- cbind(f.stat, s.num)
pb.pb <- data.frame(pb.pb)
pb.pb$transplant <- "PB->WB"
pb.pb$comparison <- "PB"
pb.pb$type <- "F.distichus"

plot(s.num.pb,f.stat.pb)
pb.lm <- lm(f.stat.pb ~ s.num.pb)
summary(pb.lm)

anova(pb.lm)
anova(wb.lm)
#               Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   0.1670     0.2009   0.831   0.4669   
# s.num.pb      0.8419     0.1421   5.926   0.0096 **
# 
# Residual standard error: 0.2377 on 3 degrees of freedom
# Multiple R-squared:  0.9213,	Adjusted R-squared:  0.8951 
# F-statistic: 35.11 on 1 and 3 DF,  p-value: 0.009603




##### Combine plots ####
f.comps <- rbind(wb.wb, pb.pb)

ft.palette <- c("#E69F00", "chartreuse3","springgreen4", "darkblue")

ggplot(f.comps, aes(x= s.num, y=f.stat, color = transplant)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=F) + 
  theme_classic() +
  scale_color_manual(values = c("chartreuse3","springgreen4"))


