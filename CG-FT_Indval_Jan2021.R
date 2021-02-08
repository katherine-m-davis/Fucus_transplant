#### F. distichus reciprocal transplant indicator taxa analysis between control sites and transplants over time ####
# Get indicator ASVs between WB low and PB controls
#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(indicspecies)
library(ggpubr)

#### Read in dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/Combined/CG_FT_rock_fucus_combined_phyloseq_r1500_Jan2021.RDS")
dimnames(tax_table(cg.ft.fr))

#### Get otu table and meta data for sample comparisions ####
ft.controls <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT" & sample_data(cg.ft.fr)$individual == "neighbor" & sample_data(cg.ft.fr)$type == "fucus") #subset phyloseq object
ft.c.meta <- as.data.frame(unclass(sample_data(ft.controls))) %>% mutate_if(is.factor, as.character) # get otu table (community data matrix) where sites/samples are rows and species/ASVs are columns
ft.c.otu <- as.data.frame(unclass(otu_table(ft.controls)))
ft.c.otu.mat <- as.matrix(ft.c.otu)
ind <- ft.c.meta$sampleid
ft.c.otu.mat <- ft.c.otu.mat[ind, ] #order matrix to match metadata
origin <- ft.c.meta$origin #get vector of indicator variable
  
#### Indispecies analysis ####
ft.ctrl <- multipatt(ft.c.otu.mat, origin, control = how(nperm=999))
summary(ft.ctrl, indvalcomp = T)

#### Import results ###
ft.indval <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Data/FT_control_INDVAL_results_Jan2021.csv")

# Add taxonomy data
cg.ft.tax <- as.data.frame(unclass(tax_table(cg.ft.fr)))
ft.indval.tax <- left_join(ft.indval, cg.ft.tax)


#### Get relative abundance of top Indicator ASVs across samples ####
ft.indval.9 <- ft.indval.tax %>% filter(indval.stat >= 0.9) # select ASVs with indval stat greater than 0.9
length(unique(ft.indval.9$ASV)) #47
# Data frame of phyloseq object
# ft.f <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT" & sample_data(cg.ft.fr)$type == "fucus")
# ftdf <- psmelt(ft.f)
# saveRDS(ftdf, file="~/Desktop/Desktop2020/CG_FT/Data/Combined/FT_df_r1500_Jan2021.RDS")
ftdf <- readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/Combined/FT_df_r1500_Jan2021.RDS")

#### Calculate prevalence and mean relative abundance ####
ftdf <- ftdf %>% group_by(treatment, OTU) %>%  mutate(n_samples_total=n()) %>% mutate(n_samples_present=n())
ftdf <- ftdf %>% mutate(perc.samples = n_samples_present/n_samples_total * 100)
ftdf <- ftdf %>% mutate(relabund = Abundance/1500 * 100)
ftdf <- ftdf %>% group_by(treatment, OTU) %>% mutate(mean.relabund = mean(relabund), sd.relabund= sd(relabund))
ftdf$key2 <- parse_number(ftdf$key)
ftdf$sample.ind <- paste(ftdf$sample.number, ftdf$key2, sep="_")

#### Filter to indicator taxa with indval stat > 0.8 ####
ftdf.ind <- ftdf %>% ungroup() %>% filter(OTU %in% ft.indval.9$ASV)
ftdf.ind <- ftdf.ind %>% mutate_if(is.factor, as.character)

# Format for plotting  
ftdf.ind$Genus[ftdf.ind$Genus == "Nonlabens"] <- "unknown_Flavobacteriaceae"
ftdf.ind$OTU2 <- paste(ftdf.ind$Genus, ftdf.ind$OTU, sep="_")
ftdf.ind$OTU3 <- paste(ftdf.ind$Family, ftdf.ind$OTU2, sep="_")
sort(unique(ftdf.ind$OTU3))
ftdf.ind$trt.type <- ifelse(ftdf.ind$treatment %in% c("PB->PB", "WB low->WB low"), "control","transplant")
ftdf.ind.f <- ftdf.ind %>% filter(type == "fucus")
length(unique(ftdf.ind$OTU2)) #76

ftdf.ind.f$color <- ifelse(ftdf.ind.f$origin == "PB", "#E69F00", "darkblue")  #add column for color corresponding to indicator site
# Create legend title with square root symbol
ft.heatmap.legend <- "\u221A ASV \nrelative \nabundance"  


#### Plot heatmap ####
ft.hm <- ggplot(ftdf.ind.f, aes(sample.ind, OTU3, fill= sqrt(relabund))) + 
  theme_classic(base_size = 5) +
  scale_fill_gradient(ft.heatmap.legend, low="#FFFFFF", high = "grey20") +
  facet_grid(~ treatment+trt.type, scales = "free_x") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = sort(unique(ftdf.ind.f$OTU3))) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.spacing.x=unit(0, "lines")) +
  xlab(expression(paste("Sample day_", italic("F. distichus")," individual")))
ft.hm


#### Make color code figure for heatmap ####
ft.indval.9$Genus[ft.indval.9$Genus == "Nonlabens"] <- "unknown_Flavobacteriaceae"
ft.indval.9$OTU2 <- paste(ft.indval.9$Genus, ft.indval.9$ASV, sep="_")
ft.indval.9$OTU3 <- paste(ft.indval.9$Family, ft.indval.9$OTU2, sep="_")
ft.indval.9$OTU3 <- gsub("\\_unknown_.*_", "_unknown_", ft.indval.9$OTU3)
ft.indval.9$xlab <- "0-0" #add random labels for plot alignment
ft.indval.9$facet <- "Treatment" #add random labels for faceted plot alignment
ft.indval.9$facet2 <- "." #add random labels for faceted plot alignment

# Make legend title for color guide
color_legend_title <- "Indicator of site"

ft.hm.col <- ggplot(ft.indval.9, aes(x=xlab, y=OTU3, fill= group)) + 
  geom_tile() +
  theme_classic(base_size =5) +
  facet_grid(~facet + facet2) +
  scale_fill_manual(color_legend_title, values = c("#E69F00", "darkblue"))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), plot.margin = rep(unit(0,"null"),4))+
  theme(axis.ticks.x= element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(expand = expansion(add=0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = sort(ft.indval.9$OTU3))  +
    theme(strip.text.x = element_blank())
ft.hm.col 


# Save heatmap color legend for adding manually
ft.heatmap.colors <- get_legend(ft.hm.col)
pdf(file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/CG-FT_FT_Indicator_Heatmap_color_legend_2021.pdf", height =2, width =2)
ft.hm.color.legend.plot <- as_ggplot(ft.heatmap.colors)
ft.hm.color.legend.plot
dev.off()

# Save heatmap relative abundance legend for adding manually
ft.heatmap.abund <- get_legend(ft.hm) 

png(file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/CG-FT_FT_Indicator_Heatmap_abundance_legend_2021.png", height =2, width =2)
ft.hm.abund.legend.plot <- as_ggplot(ft.heatmap.abund)
ft.hm.abund.legend.plot
dev.off()

ft.hm.abund.legend.plot

# Combine plots
ftheatmap <- ggarrange(ft.hm.col + rremove("legend") + rremove("xlab") + rremove("x.text"), ft.hm, nrow =1,  align="h", widths = c(1,3.5))
ftheatmap

#### Save figure
png(file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/CG-FT_FT_Indicator_Heatmap_Feb2021.png", width = 180, height = 100, 
      units = "mm", res = 400)
ggarrange(ft.hm.col + rremove("legend") + rremove("xlab") + rremove("x.text"), ft.hm, nrow =1,  align="h", widths = c(1,3.5))
dev.off()


ftheatmapplot <-ggarrange(ft.hm.col + rremove("legend") + rremove("xlab") + rremove("x.text"), ft.hm, nrow =1,  align="h", widths = c(1,3.5))
ggsave(ftheatmap,filename = "~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/ggsave_CG-FT_FT_Indicator_Heatmap_Jan2021.pdf",device = "pdf",width = 180,height = 100,unit="mm")


