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
ft.f$Genus2 <- ifelse(ft.f$Genus %in% genera.5, ft.f$Genus, "Other")
unique(ft.f$Genus2)
ft.f$Genus3 <- ifelse(ft.f$Genus2 %in% genera.5, paste(ft.f$Family, ft.f$Genus2, sep="_"), "Other")
sort(unique(ft.f$Genus3))

#### Fix awkward taxonomy assignments
unique(gsub("_unknown_.*", "_Unknown",ft.f$Genus3))
ft.f$Genus3 <- gsub("_unknown_.*", "_Unknown",ft.f$Genus3)
ft.f$Genus3[ft.f$Genus3 == "Flavobacteriaceae_Nonlabens"] <- "Flavobacteriaceae_Unknown"
ft.f$Genus3[ft.f$Genus3 == "unknown_Family_Unknown"] <- "Bacteroidia_Unknown"
ft.f$Genus3[ft.f$Genus3 == "unknown_Bacteroidia_Unknown"] <- "Bacteroidia_Unknown"
sort(unique(ft.f$Genus3))

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

#### Arrange taxonomic assignments in alphabetical order with "Other" at the beginning ####
paste0(sort(unique(ft.f.g$Genus3)), collapse="','") 
ft.f.g$Genus3 <- factor(ft.f.g$Genus3, levels=c('Other','Arenicellaceae_Arenicella','Flavobacteriaceae_Croceitalea','Flavobacteriaceae_Dokdonia','Flavobacteriaceae_Maribacter','Flavobacteriaceae_Olleya','Flavobacteriaceae_Unknown','Hyphomonadaceae_Litorimonas','Pirellulaceae_Blastopirellula','Pseudoalteromonadaceae_Pseudoalteromonas','Rhodobacteraceae_Litoreibacter','Rhodobacteraceae_Octadecabacter','Rhodobacteraceae_Planktomarina','Rhodobacteraceae_Sulfitobacter','Rhodobacteraceae_Unknown','Rubritaleaceae_Luteolibacter','Rubritaleaceae_Rubritalea','Saprospiraceae_Lewinella','Saprospiraceae_Membranicola','Saprospiraceae_Unknown','Sphingomonadaceae_Unknown','Thiohalorhabdaceae_Granulosicoccus','Vibrionaceae_Vibrio'))
sort(unique(ft.f.g$Genus3))
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
ft.r$Genus2 <- ifelse(ft.r$Genus %in% genera.5r, ft.r$Genus, "Other")
unique(ft.r$Genus2)
ft.r$Genus3 <- ifelse(ft.r$Genus2 %in% genera.5r, paste(ft.r$Family, ft.r$Genus2, sep="_"), "Other")
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

#### Arrange taxonomic assignments in alphabetical order with "Other" at the beginning ####
paste0(sort(unique(ft.r.g$Genus3)), collapse="','") 
ft.r.g$Genus3 <- factor(ft.r.g$Genus3, levels=c('Other','Cyclobacteriaceae_Tunicatimonas','Flavobacteriaceae_NS3a_marine_group','Flavobacteriaceae_Tenacibaculum','Flavobacteriaceae_Unknown',
                                                'Hyphomonadaceae_Fretibacter','Pseudoalteromonadaceae_Pseudoalteromonas','Pseudoalteromonadaceae_Psychrosphaera',
                                                'Rhizobiaceae_Pseudahrensia','Rhizobiaceae_Unknown','Rhodobacteraceae_Litoreibacter','Rhodobacteraceae_Planktomarina','Rhodobacteraceae_Unknown',
                                                'Saprospiraceae_Lewinella','Saprospiraceae_Portibacter','Saprospiraceae_Unknown','Thiohalorhabdaceae_Granulosicoccus','Vibrionaceae_Vibrio'))
sort(unique(ft.r.g$Genus3))
####################################################################################################################

####################################################################################################################
# Make a color lookup table for plotting #-------------------------------------------------------------------------#

#### Make combined list of rock and F. distichus top genera ####
g.all <- c(unique(ft.f.g$Genus3), unique(ft.r.g$Genus3))
g5_unique <- unique(g.all)

# c31.ordered <- c('#0000CD','#0b7e8d','#3a8b37','#3ebab5','#3f0482','#4ca578','#57b946','#5e56ad','#5f95cc','#898b48','#8B4513','#905ed9','#a0b92c','#ab88d3','#b865af','#ba4334','#c68337','#c954b8','#DC143C','#dc375f','#e05621','#e28e24','#FF69B4','#FF7F50','darkred','darkslategray','hotpink4','lightgoldenrod','lightskyblue','midnightblue','snow4')
# # scales::show_col(c31.ordered) #plot colors
# sort(g5_unique) # display sorted list of all taxa
# 
# # Combine taxa and colors in dataframe
# FT.color.lookup <- as.data.frame(cbind(c31.ordered, sort(g5_unique))) %>% mutate_if(is.factor, as.character)
# colnames(FT.color.lookup) <- c("FTcolor", "Genus3")
#
# Save lookup table
#write.csv(FT.color.lookup, file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Data/FT_color_lookup_table.csv")

# Read color lookup table
FT.color.lookup <- read.csv(file="~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Data/FT_color_lookup_table.csv")
FT.color.lookup <- FT.color.lookup %>% select(- 'X') %>% mutate_if(is.factor, as.character)

paste0(sort(unique(FT.color.lookup$FTcolor)), collapse="','") 

# Assign taxa names to colors
FTcolors <- FT.color.lookup$FTcolor
names(FTcolors) <- FT.color.lookup$Genus3
head(FTcolors)
scales::show_col(FTcolors)
####################################################################################################################


####################################################################################################################
#### Make stacked barplots for Fucus ------------------------------------------------------------------------------#

# Plot and save
pdf("~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/FT_barplot_Fucus_Family_Genus_ordered_9Feb2021.pdf", 
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
pdf("~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/FT_barplot_Rocks_Family_Genus_ordered_9Feb2021.pdf", 
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
sort(FT.color.lookup$Genus3)
paste0(sort(unique(FT.color.lookup$Genus3)), collapse="','") 
FT.color.lookup$Genus3 <- factor(FT.color.lookup$Genus3, levels=c('Other', 'Arenicellaceae_Arenicella','Cyclobacteriaceae_Tunicatimonas','Flavobacteriaceae_Croceitalea','Flavobacteriaceae_Dokdonia','Flavobacteriaceae_Maribacter',
                                                         'Flavobacteriaceae_NS3a_marine_group','Flavobacteriaceae_Olleya','Flavobacteriaceae_Tenacibaculum','Flavobacteriaceae_Unknown','Hyphomonadaceae_Fretibacter','Hyphomonadaceae_Litorimonas',
                                                         'Pirellulaceae_Blastopirellula','Pseudoalteromonadaceae_Pseudoalteromonas','Pseudoalteromonadaceae_Psychrosphaera','Rhizobiaceae_Pseudahrensia','Rhizobiaceae_Unknown','Rhodobacteraceae_Litoreibacter',
                                                         'Rhodobacteraceae_Octadecabacter','Rhodobacteraceae_Planktomarina','Rhodobacteraceae_Sulfitobacter','Rhodobacteraceae_Unknown','Rubritaleaceae_Luteolibacter','Rubritaleaceae_Rubritalea','Saprospiraceae_Lewinella',
                                                         'Saprospiraceae_Membranicola','Saprospiraceae_Portibacter','Saprospiraceae_Unknown','Sphingomonadaceae_Unknown','Thiohalorhabdaceae_Granulosicoccus','Vibrionaceae_Vibrio'))


lp <- ggplot(FT.color.lookup, aes(x=axis.x, y=axis.y, fill= Genus3)) + 
  theme_classic(base_size = 11)+
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values =FTcolors) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 4))

# Extract the legend. Returns a gtable
leg <- get_legend(lp)

# Convert to a ggplot and print
gg.lp <- as_ggplot(leg) 
gg.lp <- gg.lp + theme(plot.margin = unit(c(0, 6, 0, 6), "cm"))
gg.lp
####################################################################################################################
#### Combine plots with Figure A & B labels for Supplementary Figure ----------------------------------------------#
library(cowplot)

fp <- f + theme(legend.position = "none")
rp <- r + theme(legend.position = "none")


pdf("~/Desktop/Desktop2020/CG_FT/Fucus_Transplant/Figures/FT_barplot_combined_fucus-rock_ordered_9Feb2021.pdf", 
    width = 14, # define plot width and height,
    height = 12)
cowplot::plot_grid(fp, rp, gg.lp, nrow = 3, labels = c('A', 'B', ''))
dev.off()


#grid.arrange(grobs=list(fp, rp, gg.lp), labels = c("A","B",NA), nrow = 3, heights=c(1,1,0.5))



