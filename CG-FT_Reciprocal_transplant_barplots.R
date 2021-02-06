#### Fucus transplant taxa summary barplots ####

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(gridExtra)

###################################################################################################################
#### Read combined dataset ########################################################################################
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/CG_FT_combined_phyloseq_r1500_Jan2021.RDS")
cg.ft.tax <- as.data.frame(unclass(tax_table(cg.ft.fr)))

# Create or read in dataframe
# ft.ps <- subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$study == "FT")
# ft.df <- psmelt(ft.ps) # create dataframe 
# write.csv(ft.df, file = "~/Desktop/Desktop2020/CG_FT/Data/CG-FT_FT_dataframe_r1500_Jan2021.csv")

ft.df <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Data/CG-FT_FT_dataframe_r1500_Jan2021.csv")
ft.df <- ft.df %>% mutate_if(is.factor, as.character)

#Add column for individual.number
ft.df$individual.number <- sub('.*\\-', '', ft.df$sampleid_old)

# Extract number from mixed numeric & character string
regexp <- "[[:digit:]]+"
# process string
ft.df$host.individual <- str_extract(ft.df$key, regexp)
ft.df$host.individual <- as.numeric(ft.df$host.individual)


# Combine variables
ft.df$sample.day_host.number <- paste(ft.df$sample.number, ft.df$host.individual, sep="_")
ft.df$sample.day_host.number <- paste(ft.df$sample.day_host.number, ft.df$individual.number, sep=".")
####################################################################################################################

####################################################################################################################
##### Make barplots for F. distichus samples #######################################################################
ft.f <- ft.df %>% filter(type == "fucus")

#### Get relative abundance ####
ft.f <- ft.f %>%  mutate(relabund = Abundance/1500 * 100)

#### Find most abundant ASVS ####
ft.f.5 <- ft.f %>% filter(relabund > 5)
length(unique(ft.f.5$OTU)) #47
length(unique(ft.f.5$Family))  #11
length(unique(ft.f.5$Genus)) #23
genera.5 <- as.character(unique(ft.f.5$Genus)) #Group by genera
sort(genera.5)

# Update Family/Genus names 
ft.f$Genus2 <- ifelse(ft.f$Genus %in% genera.5, ft.f$Genus, "A_Other")
unique(ft.f$Genus2)
ft.f$Genus3 <- ifelse(ft.f$Genus2 %in% genera.5, paste(ft.f$Family, ft.f$Genus2, sep="_"), "A_Other")
sort(unique(ft.f$Genus3))

#### Fix awkward taxonomy assignments
unique(gsub("_unknown_.*", "_Unknown",ft.f$Genus3))
ft.f$Genus3 <- gsub("_unknown_.*", "_Unknown",ft.f$Genus3)
ft.f$Genus3[ft.f$Genus3 == "Flavobacteriaceae_Nonlabens"] <- "Flavobacteriaceae_Unknown"
ft.f$Genus3[ft.f$Genus3 == "unknown_Family_Unknown"] <- "Bacteroidia_Unknown"
ft.f$Genus3[ft.f$Genus3 == "unknown_Bacteroidia_Unknown"] <- "Bacteroidia_Unknown"
sort(unique(ft.f$Genus3))
sort(unique(ft.f$sampleid))

# Add column for sample day and host individual
# Get host individual
# regexp <- "[[:digit:]]+"
# 
# # process string
# ft.f$key2 <- str_extract(ft.f$key, regexp)
# ft.f$host.individual <- as.numeric(ft.f$key2)
# 
# # Combine variables
# ft.f$sample.day_host.number <- paste(ft.f$sample.number, ft.f$individual.number, sep="_")
# ft.f$sample.day_host.number <- paste(ft.f$sample.day_host.number, ft.f$host.individual, sep="_")
# 
# # Find replicated sample day and host individual in WB low -> WB low control
# ft.f2 <- ft.f %>% filter(! c(treatment == "WB low->WB low" & sample.number == 0 & individual == "transplant"))

ft.f.g <- ft.f %>% group_by(Genus3, treatment, sample.day_host.number, sampleid_old) %>% summarise(relabund = sum(relabund))
####################################################################################################################


####################################################################################################################
##### Make barplots for Rock samples ###############################################################################
ft.r <- ft.df %>% filter(type == "rock")

#### Get relative abundance ####
ft.r <- ft.r %>%  mutate(relabund = Abundance/1500 * 100)

#### Find most abundant ASVS ####
ft.r.5 <- ft.r %>% filter(relabund > 5)
length(unique(ft.r.5$OTU)) #46
length(unique(ft.r.5$Family))  #9
length(unique(ft.r.5$Genus)) #18
genera.5r <- as.character(unique(ft.r.5$Genus)) #Group by genera
sort(genera.5r)

# Update Family/Genus names 
ft.r$Genus2 <- ifelse(ft.r$Genus %in% genera.5r, ft.r$Genus, "A_Other")
unique(ft.r$Genus2)
ft.r$Genus3 <- ifelse(ft.r$Genus2 %in% genera.5r, paste(ft.r$Family, ft.r$Genus2, sep="_"), "A_Other")
sort(unique(ft.r$Genus3))


#### Fix awkward taxonomy assignments
unique(gsub("_unknown_.*", "_Unknown",ft.r$Genus3))
ft.r$Genus3 <- gsub("_unknown_.*", "_Unknown",ft.r$Genus3)
ft.r$Genus3[ft.r$Genus3 == "Flavobacteriaceae_Nonlabens"] <- "Flavobacteriaceae_Unknown"
ft.r$Genus3[ft.r$Genus3 == "unknown_Family_Unknown"] <- "Bacteroidia_Unknown"
ft.r$Genus3[ft.r$Genus3 == "unknown_Bacteroidia_Unknown"] <- "Bacteroidia_Unknown"
sort(unique(ft.r$Genus3))
sort(unique(ft.r$sampleid))


# Add column for sample day and host individual
# Get host individual
# regexp <- "[[:digit:]]+"
# 
# # process string
# ft.r$key2 <- str_extract(ft.r$key, regexp)
# ft.r$host.individual <- as.numeric(ft.r$key2)
# 
# # Combine variables
# ft.r$sample.day_host.number <- paste(ft.r$sample.number, ft.r$host.individual, sep="_")
# 
# # Find replicated sample day and host individual in WB low -> WB low control
# ft.r2 <- ft.r %>% filter(treatment == "WB low->WB low" & sample.number == 0) %>% select(OTU:treatment, key2) %>% distinct()
# 
# ft.r3 <- ft.r %>% filter(!c(treatment == "WB low->WB low" & date == "13-Jun-18"))

ft.r.g <- ft.r %>% ungroup() %>% group_by(Genus3, treatment, sample.day_host.number, sampleid_old) %>% summarise(relabund = sum(relabund))
####################################################################################################################

####################################################################################################################
# Make a color lookup table for plotting #-------------------------------------------------------------------------#

#### Make combined list of rock and F. distichus top genera ####
g.all <- c(unique(ft.f.g$Genus3), unique(ft.r.g$Genus3))
g5_unique <- unique(g.all)


c31.ordered <- c("snow4",
                 "#6185e8",
                 "#925218",
                 "#FF69B4",
                 "#dc375f",
                 "darkred",
                 "#DC143C",
                 "#e05621",
                 "#ba4334",
                 "#FF7F50",
                 "#4ca578",
                 "darkslategray",
                 "midnightblue",
                 "#a0b92c",
                 "#898b48",
                 "#e28e24",
                 "#c68337",
                 "#ab88d3",
                 "#b865af",
                 "hotpink4",
                 "#905ed9",
                 "#5e56ad",
                 "#3a8b37",
                 "#57b946",
                 "lightskyblue",
                 "#3ebab5",
                 "#5f95cc",
                 "slategray3",
                 "lightgoldenrod",
                 "#c954b8",
                 "#0000CD")


scales::show_col(c31.ordered) #plot colors
sort(g5_unique) # display sorted list of all taxa

# Combine taxa and colors in dataframe
FT.color.lookup <- as.data.frame(cbind(c31.ordered, sort(g5_unique))) %>% mutate_if(is.factor, as.character)
colnames(FT.color.lookup) <- c("FTcolor", "Genus3")

# Save lookup table
#write.csv(FT.color.lookup, file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Data/FT_color_lookup_table.csv")

# Read color lookup table
FT.color.lookup <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Data/FT_color_lookup_table.csv")
FT.color.lookup <- FT.color.lookup %>% select(- 'X') %>% mutate_if(is.factor, as.character)

# Assign taxa names to colors
FTcolors <- FT.color.lookup$FTcolor
names(FTcolors) <- FT.color.lookup$Genus3
head(FTcolors)
####################################################################################################################


####################################################################################################################
#### Make stacked barplots for Fucus ------------------------------------------------------------------------------#

# Plot and save
pdf("~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/FT_barplot_Fucus_Family_Genus_ordered_6Feb2021.pdf", 
    width = 13, # define plot width and height
    height = 6)
f <- ggplot(data=ft.f.g, aes(x=sample.day_host.number, y=relabund, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  facet_grid(~treatment, scales="free", space="free")+
  theme_classic(base_size = 12)+
  theme(panel.spacing.x=unit(0, "lines")) +
  theme(strip.text.x = element_text(size = 12)) +
  scale_fill_manual(values = FTcolors) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size =10)) +
  theme(axis.text.y=element_text(size =12)) +
  theme(legend.position = "bottom", legend.text=element_text(size=12)) +
  ylab("Relative abundance") +
  xlab("Sample day_Host individual") +
  guides(fill=guide_legend(ncol=3)) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "F. distichus") +
  theme(plot.title = element_text(face = "italic"))
f
dev.off()
####################################################################################################################


####################################################################################################################
#### Make stacked barplots for Rocks ------------------------------------------------------------------------------#

# Plot and save
pdf("~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/FT_barplot_Rocks_Family_Genus_ordered_5Feb2021.pdf", 
    width = 13, # define plot width and height
    height = 6)
r <- ggplot(data=ft.r.g, aes(x=sample.day_host.number, y=relabund, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  facet_grid(~treatment, scales="free", space="free")+
  theme_classic(base_size = 12)+
  theme(panel.spacing.x=unit(0, "lines")) +
  theme(strip.text.x = element_text(size = 12)) +
  scale_fill_manual(values = FTcolors) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size =10)) +
  theme(axis.text.y=element_text(size =12)) +
  theme(legend.position = "bottom", legend.text=element_text(size=12)) +
  ylab("Relative abundance") +
  xlab("Sample day_Host individual") +
  guides(fill=guide_legend(ncol=3)) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Rock substrate") 
r
dev.off()
####################################################################################################################




####################################################################################################################
#### Get combined legend for Supplementary Figure -----------------------------------------------------------------#
FT.color.lookup$axis.x <- seq_along(1:length(FT.color.lookup$Genus3))
FT.color.lookup$axis.y <- 1

lp <- ggplot(FT.color.lookup, aes(x=axis.x, y=axis.y, fill= Genus3)) + 
  theme_classic(base_size = 12)+
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values =FTcolors) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 4))

# Extract the legend. Returns a gtable
leg <- get_legend(lp)

# Convert to a ggplot and print
gg.lp <- as_ggplot(leg) 
gg.lp <- gg.lp + theme(plot.margin = unit(c(0, 1, 0, 4), "cm"))
####################################################################################################################
#### Combine plots with Figure A & B labels for Supplementary Figure ----------------------------------------------#
library(ggpubr)

fp <- f + theme(legend.position = "none")
rp <- r + theme(legend.position = "none")


pdf("~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/FT_barplot_combined_fucus-rock_ordered_6Feb2021.pdf", 
    width = 14, # define plot width and height,
    height = 12)
grid.arrange(fp, rp, gg.lp, nrow = 3, heights=c(1,1,0.5))
dev.off()


