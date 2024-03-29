#### Fucus reciprocal transplant PCoA plots ####

#### Libraries ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(gridExtra)
library(rstatix)

#### Read combined dataset ####
cg.ft.fr <-readRDS(file="~/Desktop/Desktop2020/CG_FT/Data/CG_FT_combined_phyloseq_r1500_Jan2021.RDS")


#### Figure 4 plots reciprocal transplant by habitat type and sample day ####
# Filter dataset for transplant experiment 
unique(sample_data(cg.ft.fr)$project)
ft <-subset_samples(cg.ft.fr, sample_data(cg.ft.fr)$project == "fucus transplant")

# Update type from fucus to macroalgal host
sample_data(ft)$type <- ifelse(sample_data(ft)$type == "fucus", "macroalgal host", sample_data(ft)$type)
unique(sample_data(ft)$type)

# Only control samples
ft.c <-subset_samples(ft, ! sample_data(ft)$individual =="transplant")
sample_data(ft.c)$individual

# Fucus and Rock ordiation
pcoa.ft.c <- ordinate(ft.c, method = "PCoA", distance = "bray")

# Colors + Shapes
ft.palette <- c("#E69F00", "darkblue")

#### Plot and save Figure 4A ####
pdf(file="~/Desktop/Desktop2020/Manuscripts/CG_FT_Manuscript_2020/CG_FT_Figures_FINAL_Feb2021/CG-FT_ft_controls_PCoA_Jan21.pdf",
    width = 5, height = 3.5)
ftc <- plot_ordination(ft.c, pcoa.ft.c, type = "samples", color="origin", shape ="type") +
  theme_classic(base_size = 12) + 
  scale_color_manual(values=ft.palette) + 
  geom_point(size=3) +
  labs(color='Site', shape = 'Sample type')
ftc
dev.off()


# Sample numbers for controls
ft.c.df <- as.data.frame(unclass(sample_data(ft.c)))
ft.c.df %>%  group_by(origin, type) %>% summarise(n())
# 1 PB     fucus    21
# 2 PB     rock     13
# 3 WB low fucus    17
# 4 WB low rock     17


#### Fucus Transplant by sample day and treatment, Figure 4B ####
# Subset to only fucus host samples
ft.f <- subset_samples(ft, sample_data(ft)$type == "macroalgal host")  
sample_data(ft.f)$sample.number <- as.character(sample_data(ft.f)$sample.number)

# Ordinate
pcoa.ft.f <- ordinate(ft.f, method = "PCoA", distance = "bray")

# Day shapes and trmt colors
ft.day.shapes <- c(16,1,2,6,5,8)
ft.trmt.colors <- c("#E69F00","#463806", "#7ACE54", "darkblue")

#### Plot and save Figure 4B ####
pdf(file="~/Desktop/Desktop2020/Manuscripts/CG_FT_Manuscript_2020/CG_FT_Figures_FINAL_Feb2021/CG-FT_ft_treatment_day_shapes_PCoA_Jan21.pdf",
    width = 5, height = 3.5)
p2 <- plot_ordination(ft.f, pcoa.ft.f, type = "samples", color="treatment", shape ="sample.number") +
  theme_classic(base_size = 12) +
  geom_point(size=3) +
  scale_shape_manual(values= ft.day.shapes) +
  scale_color_manual(values = ft.trmt.colors) +
  labs(color='Treatment', shape = 'Sample day') +
  labs(title = "F. distichus") +
  theme(plot.title = element_text(face = "italic"))
p2$layers <- p2$layers[-1]
p2
dev.off()


#### Rock Transplant by sample day and treatment, Figure 4B ####
# Subset to only rock substrate samples
ft.r <- subset_samples(ft, sample_data(ft)$type == "rock") 
sample_data(ft.r)$sample.number <- as.character(sample_data(ft.r)$sample.number)

# Ordinate
pcoa.ft.r <- ordinate(ft.r, method = "PCoA", distance = "bray")

# Day shapes and trmt colors
ft.day.shapes <- c(16,1,2,6,5,8)
ft.trmt.colors <- c("#E69F00","#463806", "#7ACE54", "darkblue")

#### Plot and save Figure 4B ####
pdf(file="~/Desktop/Desktop2020/Manuscripts/CG_FT_Manuscript_2020/CG_FT_Figures_FINAL_Feb2021/CG-FT_ft_treatment_day_shapes_PCoA_Jan21.pdf",
    width = 5, height = 3.5)
p3 <- plot_ordination(ft.r, pcoa.ft.r, type = "samples", color="treatment", shape ="sample.number") +
  theme_classic(base_size = 12) +
  geom_point(size=3) +
  scale_shape_manual(values= ft.day.shapes) +
  scale_color_manual(values = ft.trmt.colors) +
  labs(color='Treatment', shape = 'Sample day') +
  labs(title = "Rock substrate") +
  theme(legend.position = "none")
p3$layers <- p3$layers[-1]
p3
dev.off()


# Extract the legends. Returns a gtable
legA <- get_legend(ftc)
legB <- get_legend(p2)
# Convert to a ggplot and print
gg.A <- as_ggplot(legA) 
gg.B <- as_ggplot(legB) 
gg.B 


gg.A <-gg.A + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
####################################################################################################################
#### Combine plots with Figure A & B labels for Supplementary Figure ----------------------------------------------#
library(cowplot)

ftcl <- ftc + theme(legend.position = "none")
p2l <- p2 + theme(legend.position = "none")
p3l <- p3 + theme(legend.position = "none")


# Save combined Figure 4 with labels
pdf(file="~/Desktop/Manuscripts_2021/CG_FT_2021/CG-FT_Figures_Resubmission_2021/CG-FT_Figure4_FT_PCoA_Feb2021.pdf",
    width =7, height = 11)
ggarrange(ggarrange(ftcl, NULL, p2l, NULL, p3l, ncol = 1, labels = c("A","","B", "", "C"), heights=c(1,0.05,1,0.05,1)),
  ggarrange(gg.A, NULL, gg.B, NULL, NULL, ncol = 1, heights=c(1,0.25,1,0.05,0.7)), ncol = 2, widths=c(1, 0.5))
dev.off()



          
