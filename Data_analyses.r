##### Loading packages #####
library(ggtree) 
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggtreeExtra)
library(ggstance)
library(ggsignif)
library(ggpmisc)
library(broom)
library(cowplot)
library(plyr)
library(psych)
library(RVAideMemoire)
library(vegan)
library(ggcorrplot)
library(ape)
library(tAI)
library(dplyr)
library(tidyr)
library(phytools)
library(tidyverse)
library(ggrepel)
library(Biostrings)
library(threadr)
library(ggh4x)





##### Loading data #####

# Loading phylogeny
tree = read.tree("Source_files/Source_fig1.tree")

# Loading species data
species_data = read.delim("Source_files/Source_species_metadata.txt")
species_data = as.data.frame(species_data)

# Loading genome stats
genome_stats = read.delim("Source_files/Source_genome_stats.txt")

# Loading growth results
growth_parameters = read.delim("Source_files/Source_growth_parameters.txt")

# Loading BUSCO results
busco = read.delim("Source_files/Source_BUSCO_basidiomycota.txt")

# Loading tRNA gene copy number
tRNA_copy_number = read.delim("Source_files/Source_tRNA_copy_number.txt")

# Loading tRNA genetic distance
tRNA_genetic_distance = read.delim("Source_files/Source_tRNA_genetic_distance.txt")

# Loading growth temperatures
growth_temp = read.delim("Source_files/Source_temp_growth.txt")





##### ——— MAIN FIGURES ——— #####

#### Figure 1 ####


## Process genome stats for plotting ##

# For plotting in bars, it is required to subtract carbohydrate and lipid proteins from the total number of proteins
genome_stats$other_proteins = genome_stats$genes - genome_stats$carbo_genes - genome_stats$lipid_genes

# For plotting in bars, it is required to subtract transposable elements and repetitive DNA from the genome size
# Transposable elements were not possible to calculate for some genomes and is therefore missing (NA)
# Replace NA with 0 for calculations
genome_stats$TEcoverage0 = genome_stats$TEcoverage
genome_stats$TEcoverage0[is.na(genome_stats$TEcoverage0)] <- 0
genome_stats$non_repetitive_DNA = genome_stats$genome_size - genome_stats$masked_repeats - genome_stats$TEcoverage0 

## Test calculations
genome_stats$genome_size == genome_stats$non_repetitive_DNA + genome_stats$masked_repeats + genome_stats$TEcoverage0 

# Reshape data.frame into tabular 
genome_stats_melt = melt(genome_stats)



# Draws tree
Fig1 <- ggtree(tree, branch.length='branch.length', size = 1) + 
  geom_treescale(x=0, y=0, fontsize = 2, linesize = 1)


# Adds the tip colour by substrate
Fig1 <- Fig1 %<+% species_data + 
  geom_tiplab(aes(label = factor(Name)), size=3.5) +
  geom_nodelab(size=2.5,
               hjust = 1.5,
               vjust=-0.5)

# Adds genome size panel
Fig1 = Fig1 + geom_facet(panel ='Genome size (Mbp)', 
                         data = subset(genome_stats_melt, variable == "non_repetitive_DNA" | 
                                                          variable == "masked_repeats" | 
                                                          variable == "TEcoverage"), 
                         geom = geom_bar, 
                         stat = "identity", 
                         position = 'stack',
                         aes(x=value/1000000, fill = factor(variable, levels=c("non_repetitive_DNA", 
                                                                               "masked_repeats", 
                                                                               "TEcoverage"))), 
                         colour="black",
                         size = 0.6,
                         orientation = 'y', 
                         width = .6)

# Adds proteome panel
Fig1 = Fig1 + geom_facet(panel ='Proteome (N)', 
                         data = subset(genome_stats_melt, variable == "other_proteins" | 
                                                          variable == "carbo_genes" | 
                                                          variable == "lipid_genes"),
                         geom = geom_bar,
                         stat = "identity", 
                         position = 'stack',
                         aes(x = value, fill = factor(variable, levels=c("other_proteins",
                                                                         "lipid_genes", 
                                                                         "carbo_genes"))),
                         colour = "black",
                         size = 0.6,
                         orientation = 'y', 
                         width = .6) +
  theme_tree2()+
  theme(legend.position = "right", 
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.background = element_blank())


# Adds tRNA panel
Fig1 = Fig1 + geom_facet(panel ='tRNA content (N)', 
                         data = subset(genome_stats_melt, variable == "tRNAs"),
                         geom = geom_col, 
                         aes(x = value, fill = variable), 
                         colour="black",
                         size = 0.6,
                         orientation = 'y', 
                         width = .6) +
  theme_tree2() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.background = element_blank())+
  scale_fill_manual(name = "", 
                    labels = c("non_repetitive_DNA" = "Non-repetitive DNA", 
                               "masked_repeats" = "Repetitive DNA", 
                               "TEcoverage" = "Transposable Elements",
                               "other_proteins" = "Other proteins",
                               "lipid_genes" = "Lipid metabolism",
                               "carbo_genes" = "Carbohydrate metabolism",
                               "tRNAs" = "tRNA genes"),
                    values = c("tRNAs" = "#8C649B",
                               "carbo_genes" = "#ff3333",
                               "lipid_genes" = "#ffff00",
                               "other_proteins" = "#9C7A70",
                               "non_repetitive_DNA" = "#619B6D", 
                               "masked_repeats" = "#CCCCCC",  
                               "TEcoverage" = "black")) +
  xlim_tree(1)


# Scales panels
Fig1 = facet_widths(Fig1, c(Tree = 2, "Genome size (Mbp)" = 1, "Proteome (N)" = 1, "tRNA content (N)" = 1))
Fig1

# Export Fig1, before editing
pdf("Figures/Fig1_raw.pdf", height = 10, width = 15)
Fig1
dev.off()





#### Figure 2 ####

Fig2 = ggarrange(
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="genome_size") , aes(x=Lifestyle, y=value/1000000, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Genome size (Mbps)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0.1,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(15,35), breaks = seq(15, 35, 5), expand = expansion(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="GC_genome") , aes(x=Lifestyle, y=value*100, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("GC content (%)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0.1,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(44,65), breaks = seq(45, 65, 5), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="tRNAs") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("tRNA genes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0.1,0,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(0, 1550), breaks = seq(0, 1550, 250), expand = expand_scale(mult = c(0, .1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="repeat_content") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Repeat content (%)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="TE_content") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("TE content (%)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.6, 0.5), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=6)", "Environmental" = "Saprotroph\n(N=13)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1")) 
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="genes") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Proteins (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(6000, 10000), breaks = seq(6000, 10000, 1000), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="secreted_genes") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Secreted proteins (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(200, 800), breaks = seq(200, 800, 100), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="CAZymes") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("CAZymes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(150, 310), breaks = seq(150, 310, 50), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="CAZymes_secreted") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Secreted CAZymes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(50, 125), breaks = seq(50, 125, 25), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="carbo_genes") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Carbohydrate genes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(300, 505), breaks = seq(300, 505, 50), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1")) 
  
  ,
  
  ggplot(subset(genome_stats_melt, Order == "Trichosporonales" & One_strain_per_species == "Yes" &  variable=="lipid_genes") , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Lipid genes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0.2,0,0), 'lines')) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.3, map_signif_level=T)+
    scale_y_continuous(limits = c(150, 310), breaks = seq(150, 310, 50), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=10)", "Environmental" = "Saprotroph\n(N=20)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1")) 
  
  ,
  
  legend = "none",
  align = "hv",
  ncol = 3,
  nrow = 4
)


pdf("Figures/Fig2.pdf", height = 8, width = 6.2, useDingbats = FALSE)
Fig2
dev.off()

tiff("Figures/Fig2.tiff", height = 8, width = 6.2, units = "in", compression = "lzw+p", res = 360)
Fig2
dev.off()



#### Figure 3 ####

S = read.delim("Source_files/Source_S.txt")

S$Lifestyle_2[ S$Genus == "Apiotrichum"] <- "Emerging pathogens\n(Apiotrichum)"
S$Lifestyle_2[ S$Genus == "Cutaneotrichosporon"] <- "Common pathogens\n(Cutaneotrichosporon, Trichosporon, Cryptococcus)"
S$Lifestyle_2[ S$Genus == "Trichosporon"] <- "Common pathogens\n(Cutaneotrichosporon, Trichosporon, Cryptococcus)"
S$Lifestyle_2[ S$Genus == "Cryptococcus"] <- "Common pathogens\n(Cutaneotrichosporon, Trichosporon, Cryptococcus)"



Fig3 = ggplot(subset(S, !is.na(Lifestyle_2) ), aes(x=Lifestyle, y=as.numeric(lipids)/as.numeric(carbs), fill = Lifestyle))+
  geom_boxplot(colour = "black")+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), alpha = 0.7, colour = "black") +
  ylab(expression(paste(frac(S ["lipid transport and metabolism"], S ["carbohydrate transport and metabolism"]))))+
  xlab("") +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = 7),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(colour = "black"))+
  geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.2, margin_top = 0.07, 
              map_signif_level=function(p) sprintf("italic(P) == %.2g", p), parse = T)+
  facet_wrap(Lifestyle_2 ~ ., nrow = 1)+
  scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
  scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=XX)", 
                            "Environmental" = "Saprotrophic\n(N=XX)"))+
  scale_y_continuous(breaks = seq(0.75, 1.20, 0.05), limits = c(0.75, 1.20), expand = expansion(mult = c(0, 0)))


pdf("Figures/Fig3_raw.pdf", width = 6.2, height = 4)
Fig3
dev.off()




#### Figure 4 ####

Fig4 = ggarrange(
  
  ggplot(growth_parameters, 
         aes(x = S, y = max_gr, colour = medium)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method="lm", se=FALSE, fullrange=F, level=0.95, linewidth = 1.3) +
    stat_fit_glance(method = 'lm',
                    geom = "text_npc",
                    method.args = list(formula = y ~ x),
                    aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value))),
                    label.x = 'left', label.y = c(0.95, 0.88), size = 3.18, fontface = "bold") +
    geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, linewidth = 1, colour = "black", linetype = "dashed") +
    stat_fit_glance(method = 'lm', 
                    geom = "text_npc",
                    method.args = list(formula = y ~ x), colour = "black",
                    aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value))),
                    label.x = 'left', label.y = 0.81, size = 3.18, fontface = "bold") +
    xlab("S") +
    ylab(expression(paste("Maximum growth rate ( ", h^-1, ")")))  +
    theme_bw() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank()) +
    labs(colour = "") +
    scale_colour_manual(name = "Substrate", labels=c('YNBoleic' = "Lipids (Oleic acid)", 'PD' = "Carbohydrates (Potato dextrose)"),
                        values = c('PD' = "dodgerblue",
                                   'YNBoleic' = "darkorange1")) +  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0.01,0.01)) +
    scale_x_continuous(limits = c(0.5,0.86), breaks = seq(0.5,0.86, 0.05), expand = c(0.01,0.01)) +
    labs(tag = expression(bold("A"))) 
  
  ,
  
  ggplot(growth_parameters, 
         aes(x = S, y = lag_time, colour = medium)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method="lm", se=FALSE, fullrange=F, level=0.95, linewidth = 1.3) +
    stat_fit_glance(method = 'lm',
                    geom = "text_npc",
                    method.args = list(formula = y ~ x),
                    aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value))),
                    label.x = 'left', label.y = c(0.19, 0.12), size = 3.18, fontface = "bold") +
    geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, linewidth = 1, colour = "black", linetype = "dashed") +
    stat_fit_glance(method = 'lm', 
                    geom = "text_npc",
                    method.args = list(formula = y ~ x), colour = "black",
                    aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value))),
                    label.x = 'left', label.y = 0.05, size = 3.18, fontface = "bold") +
    xlab("S") +
    ylab("Lag time (h)")  +
    theme_bw() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank()) +
    labs(colour = "") +
    scale_colour_manual(name = "Substrate", labels=c('YNBoleic' = "Lipids (Oleic acid)", 'PD' = "Carbohydrates (Potato dextrose)"),
                        values = c('PD' = "dodgerblue",
                                   'YNBoleic' = "darkorange1")) +  scale_y_continuous(limits = c(0,35), breaks = seq(0, 35, 5), expand = c(0.01,0.01)) +
    scale_x_continuous(limits = c(0.5,0.86), breaks = seq(0.5,0.86, 0.05), expand = c(0.01,0.01)) +
    labs(tag = expression(bold("B"))) 
  
  ,
  
  common.legend = T, 
  legend = "bottom", 
  align = "h", 
  nrow = 1)


# Save plot
pdf("Figures/Fig4.pdf", width = 6.29, height = 3.5)
Fig4
dev.off()


tiff("Figures/Fig4.tiff", width = 6.29, height = 3.5, units = "in", compression = "lzw+p", res = 360)
Fig4
dev.off()



#### Figure 5 ####

growth_OD_substrates = read.delim("Source_files/Source_growth_OD_substrates.txt")

Fig5A = ggplot(subset(growth_OD_substrates, Lifestyle == "Environmental"), aes(x = Minutes/1440, y = OD_avg, colour = Temperature)) +
  geom_line(linewidth = 1.5, alpha = 0.6) +
  facet_nested(. ~ Lifestyle + Genus + sample, labeller = labeller(Lifestyle = c("Clinical" = "Opportunistic pathogens", "Environmental" = "Saprotrophs")), scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.background.x = element_rect(colour = "black"),
        strip.text = element_text(colour = "black", face = "bold", size = 8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.margin=margin(0,0,0,0),
        panel.spacing =  unit(0.1, "lines")) +
  scale_x_continuous(limits = c(0,7), breaks = c(1, 3, 5, 7), expand = c(0,0)) +
  xlab("Day") +
  ylab(expression(paste(OD [" 600 nm"]))) +
  scale_colour_manual(values = c("33 °C" = "#440154" ,"37 °C" = "#22AA88")) +
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0, 1.2, 0.2), expand = c(0,0)) 


Fig5B = ggplot(subset(growth_OD_substrates, Lifestyle == "Clinical"), aes(x = Minutes/1440, y = OD_avg, colour = Temperature)) +
  geom_line(linewidth = 1.5, alpha = 0.6) +
  facet_nested(. ~ Lifestyle + Genus + sample, labeller = labeller(Lifestyle = c("Clinical" = "Opportunistic pathogens", "Environmental" = "Saprotrophs")), scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.background.x = element_rect(colour = "black"),
        strip.text = element_text(colour = "black", face = "bold", size = 8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.margin=margin(0,0,0,0),
        panel.spacing =  unit(0.1, "lines")) +
  scale_x_continuous(limits = c(0,7), breaks = c(1, 3, 5, 7), expand = c(0,0)) +
  xlab("Day") +
  ylab(expression(paste(OD [" 600 nm"]))) +
  scale_colour_manual(values = c("33 °C" = "#440154" ,"37 °C" = "#22AA88")) +
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0, 1.2, 0.2), expand = c(0,0)) 


pdf("Figures/Fig5.pdf", width = 6.29, height = 4)
ggarrange(Fig5A, Fig5B, common.legend = T, legend = "bottom", align = "hv", nrow = 2, ncol=1)
dev.off()

tiff("Figures/Fig5.tiff", width = 6.29, height = 4, units = "in", compression = "lzw+p", res = 360)
ggarrange(Fig5A, Fig5B, common.legend = T, legend = "bottom", align = "hv", nrow = 2, ncol=1)
dev.off()





##### ——— EXTENDED DATA FIGURES ——— #####

#### Extended Data Figure 1 ####

# For plotting purposes, Cryptococcus has been labeled ZCryptococcus

busco = subset(busco, select = -c(Lineage, Complete))
busco = melt(busco)

Extended_Data_Fig1 = ggplot(busco, aes(y=Species, x=value, fill=variable)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A"), labels = c("Complete and single-copy BUSCOs (S)", "Complete and duplicated BUSCOs (D)", "Fragmented BUSCOs (F)", "Missing BUSCOs (M)"))+
  scale_y_discrete(limits=rev, labels = c("ZCryptococcus amylolentus CBS 6039" = "Cryptococcus amylolentus CBS 6039",
                                          "ZCryptococcus deneoformans JEC21" = "Cryptococcus deneoformans JEC21",
                                          "ZCryptococcus floricola DSM 27421" = "Cryptococcus floricola DSM 27421",
                                          "ZCryptococcus gattii WM276" = "Cryptococcus gattii WM276")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.position = "bottom",
        legend.justification = c(1.5,0),
        legend.margin = margin(t = -0.2, r = 0, b = 1, l = 0, unit = "pt"),
        legend.title = element_blank(),
        legend.key.size=unit(8,"point"),
        legend.text = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank()) +
  xlab("BUSCO") +
  ylab("") +
  guides(fill=guide_legend(ncol=2))+
  scale_x_continuous(expand = c(0, 0, 0, 0.05)) +
  annotate("text", label = "*", y = 29, x = 1.05, size = 5, vjust = 0.75, family = "serif") +
  annotate("text", label = "*", y = 28, x = 1.05, size = 5, vjust = 0.75, family = "serif") +
  annotate("text", label = "*", y = 14, x = 1.05, size = 5, vjust = 0.75, family = "serif") +
  annotate("text", label = "*", y = 11, x = 1.05, size = 5, vjust = 0.75, family = "serif") +
  annotate("text", label = "*", y = 10, x = 1.05, size = 5, vjust = 0.75, family = "serif")


# Save plot
pdf("Figures/Extended_Data_Fig1.pdf", height = 7.8, width = 5.8)
Extended_Data_Fig1
dev.off()

tiff("Figures/Extended_Data_Fig1.tiff", height = 7.8, width = 5.8, units = "in", compression = "lzw+p", res = 360)
Extended_Data_Fig1
dev.off()





#### Extended Data Figure 2 ####

Extended_Data_Fig2 = ggarrange(
  
  ggplot(subset(genome_stats, Order == "Trichosporonales" & One_strain_per_species == "Yes"), aes(x=genome_size/1000000, y=repeat_content))+
  geom_point()+
  theme_bw()+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"))+
  xlab("Genome size (Mbps)")+
  ylab("Repeat content (%)")+
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, colour = 'red', linewidth = 1.3) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', after_stat(r.squared), after_stat(p.value))),
                  label.x = 'left', label.y = 'top', size = 3) +
  scale_x_continuous(limits = c(15,35), breaks = seq(15, 35, by = 5), expand = c(0.03, 0)) +
  scale_y_continuous(limits = c(5, 25), breaks = seq(5, 25, by = 5), expand = c(0.03, 0)) +
  labs(tag = expression(bold("A")))
  
  ,
  
  ggplot(subset(genome_stats, Order == "Trichosporonales" & One_strain_per_species == "Yes"), aes(x=genome_size/1000000, y=TE_content))+
  geom_point()+
  theme_bw()+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"))+
  xlab("Genome size (Mbps)")+
  ylab("TE content (%)")+
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, colour = 'red', linewidth = 1.3) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', after_stat(r.squared), after_stat(p.value))),
                  label.x = 'left', label.y = 'top', size = 3) +
  scale_x_continuous(limits = c(15,35), breaks = seq(15, 35, by = 5), expand = c(0.03, 0)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.6, by = 0.2), expand = c(0.03, 0)) +
  labs(tag = expression(bold("B"))) 

  ,
  
  align = "hv", legend = "none")



# Save plot
pdf("Figures/Extended_Data_Fig2.pdf", height = 3, width = 7)
Extended_Data_Fig2
dev.off()

tiff("Figures/Extended_Data_Fig2.tiff", height = 3, width = 7, units = "in", compression = "lzw+p", res = 360)
Extended_Data_Fig2
dev.off()


#### Extended Data Figure 3 ####


Extended_Data_Fig3 = ggarrange(
  
  ggplot(subset(genome_stats_melt, variable=="genome_size" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value/1000000, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Genome size (Mbps)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(15,35), breaks = seq(15, 35, 5), expand = expansion(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "A") +
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="GC_genome" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value*100, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("GC content (%)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(44,65), breaks = seq(45, 65, 5), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "B") +
    facet_grid(. ~ Genus, scales = "free", space = "free")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="tRNAs" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("tRNA genes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(0, 1550), breaks = seq(0, 1550, 250), expand = expand_scale(mult = c(0, .1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "C") +
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="repeat_content" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Repeat content (%)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "D") +
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="TE_content" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("TE content (%)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.6, 0.5), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "E") + 
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="genes" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Proteins (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(6000, 10000), breaks = seq(6000, 10000, 1000), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "F") + 
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="secreted_genes" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Secreted proteins (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(200, 800), breaks = seq(200, 800, 100), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "G") + 
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="CAZymes" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")), aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("CAZymes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(150, 310), breaks = seq(150, 310, 50), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "H") + 
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable=="CAZymes_secreted" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum") & value != is.na(value)) , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Secreted CAZymes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(50, 125), breaks = seq(50, 125, 25), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
    labs(tag = "I") +
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable == "carbo_genes" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")) , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Carbohydrate genes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(300, 505), breaks = seq(300, 505, 50), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1")) +
    labs(tag = "J") +
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  ggplot(subset(genome_stats_melt, variable == "lipid_genes" & One_strain_per_species == "Yes" & (Genus == "Trichosporon" | Genus == "Cryptococcus" | Genus == "Cutaneotrichosporon" | Genus == "Apiotrichum")) , aes(x=Lifestyle, y=value, fill = Lifestyle))+
    geom_boxplot(colour="black", outlier.colour = "black")+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4), alpha = 0.6, size = 1, colour = "black") +
    ylab("Lipid genes (N)")+
    xlab("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'lines')) +
    scale_y_continuous(limits = c(150, 310), breaks = seq(150, 310, 50), expand = expand_scale(mult = c(0, 0.1)))+
    labs(fill = "Lifestyle") +
    scale_x_discrete(labels=c("Clinical" = "Opportunistic\npathogen\n(N=#)", "Environmental" = "Saprotroph\n(N=#)"))+
    scale_fill_manual(values = c("darkorange1", "dodgerblue1")) +
    labs(tag = "K") +
    facet_grid(. ~ Genus, scales = "free", space = "free_x")
  
  ,
  
  align = "hv",
  ncol = 1)


pdf("Figures/Extended_Data_Fig3_raw.pdf", height = 22, width = 6)
Extended_Data_Fig3
dev.off()




#### Extended Data Figure 4 ####

Extended_Data_Fig4 = ggarrange( 
  
  ggplot(subset(genome_stats, Order == "Trichosporonales" & One_strain_per_species == "Yes"), aes(x=genome_size/1000000, y=tRNAs)) +
  geom_point() +
  theme_bw() +
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "bottom") +
  xlab("Genome size (Mbps)") +
  ylab("tRNA genes (N)") +
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, colour = 'black', linetype = "dashed", linewidth = 1.3) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value))),
                  label.x = 'left', label.y = 0.83, size = 3) +
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, size = 1.3, mapping = aes(colour = Lifestyle)) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value)), colour = Lifestyle),
                  label.x = 'left', label.y = c(0.95,0.89), size = 3) +
  scale_colour_manual(name = "Lifestyle",
                      values = c("darkorange1", "dodgerblue1"), 
                      labels=c("Clinical" = "Opportunistic pathogen",
                               "Environmental" = "Saprotrophic"))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  labs(tag = expression(bold("A"))) 

,

  ggplot(subset(genome_stats, Order == "Trichosporonales" & One_strain_per_species == "Yes"), aes(x=TE_content, y=tRNAs))+
  geom_point() +
  theme_bw() +
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "bottom") +
  xlab("TE content (%)") +
  ylab("tRNA genes (N)") +
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.85, colour = 'black', linetype = "dashed", linewidth = 1.3) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value))),
                  label.x = 'left', label.y = 0.83, size = 3)+
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, size = 1.3, mapping = aes(colour = Lifestyle)) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', stat(r.squared), stat(p.value)), colour = Lifestyle),
                  label.x = 'left', label.y = c(0.95,0.89), size = 3) +
  scale_colour_manual(name = "Lifestyle",
                      values = c("darkorange1", "dodgerblue1"), 
                      labels=c("Clinical" = "Opportunistic pathogen",
                               "Environmental" = "Saprotrophic")) +
  labs(tag = expression(bold("B")))

,

align = "hv", common.legend = TRUE, legend = "bottom")

# Save plot
pdf("Figures/Extended_Data_Fig4.pdf", height = 4, width = 7)
Extended_Data_Fig4
dev.off()

tiff("Figures/Extended_Data_Fig4.tiff", height = 4, width = 7, units = "in", compression = "lzw+p", res = 360)
Extended_Data_Fig4
dev.off()




#### Extended Data Figure 5 ####

Extended_Data_Fig5 = ggplot(S, aes(x = Lifestyle, y=as.numeric(lipids)/as.numeric(carbs), fill = Lifestyle))+
  geom_boxplot(colour = "black") +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), alpha = 0.7) +
  ylab(expression(paste(frac(S ["lipid transport and metabolism"], S ["carbohydrate transport and metabolism"]))))+
  xlab("Lifestyle") +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = 5),
        axis.text = element_text(colour = "black", size = 6),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", size = 8))+
  geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.2, margin_top = 0.07, 
              map_signif_level=T, textsize = 2)+
  facet_grid(. ~ Genus, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = c("darkorange1", "dodgerblue1"))+
  scale_x_discrete(labels=c("Clinical" = "OP", 
                            "Environmental" = "S"))+
  scale_y_continuous(breaks = seq(0.6, 1.2, 0.05), limits = c(0.6, 1.2))


pdf("Figures/Extended_Data_Fig5.pdf", width = 6.29, height = 3)
Extended_Data_Fig5
dev.off()

tiff("Figures/Extended_Data_Fig5.tiff", width = 6.29, height = 3, units = "in", compression = "lzw+p", res = 360)
Extended_Data_Fig5
dev.off()




##### ——— SUPPLEMENTARY FIGURES ——— #####


#### Supplementary Figure 1 ####

Supplementary_Fig1 = ggplot(subset(genome_stats, Order == "Trichosporonales"), aes(x=tRNAs, fill = Lifestyle)) +  
  geom_bar(colour = "black", position="stack") +
  scale_fill_manual(name = "Lifestyle", labels = c("Opportunistic pathogens", "Saprotrophic"), values=c("Clinical" = "darkorange1","Environmental" = "dodgerblue"))+
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        legend.key.size = unit(0.4,"cm"),
        legend.position.inside = c(0.85, 0.8)) +
  scale_x_binned(limits = c(0, 1500), breaks = seq(0, 1500, by = 100), expand = c(0.01, 0.1)) +
  scale_y_continuous(breaks = seq(0, 11, by = 1), expand = c(0, 0.1)) +
  xlab("Total number of tRNA genes") +
  ylab("Number of species") +
  geom_vline(aes(xintercept = median(tRNAs), 
                 linetype="median"), colour = "black", linewidth=1, key_glyph = "path") +
  geom_vline(aes(xintercept = mean(tRNAs), 
                 linetype="mean"), colour = "black", linewidth=1, key_glyph = "path") +
  scale_linetype_manual(name = "Statistics", values = c(median = "dashed", mean = "dotted")) +
  guides(linetype = guide_legend(order = 2, position = "inside", keywidth = unit(1.2,"cm")), 
         fill = guide_legend(order = 1, position = "inside"))



pdf("Figures/Supplementary_Fig1.pdf", width = 7, height = 4)
Supplementary_Fig1
dev.off()

tiff("Figures/Supplementary_Fig1.tiff", width = 7, height = 4, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig1
dev.off()



#### Supplementary Figure 2 ####

unique_anticodons = read.delim("Source_files/Source_unique_anticodons.txt")

Supplementary_Fig2 = ggplot(unique_anticodons, aes(x = unique_anticodons)) + 
  geom_histogram(bins=8, colour = "white", fill = "black")+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank())+
  scale_x_continuous(breaks = seq(39, 46, by = 1), expand = c(0, 0.05))+
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0.1))+
  xlab("Number of unique tRNA anticodon types")+
  ylab("Number of species")+
  geom_vline(xintercept = median(unique_anticodons$unique_anticodons), linetype="dashed", 
             color = "darkgrey", linewidth=1)


pdf("Figures/Supplementary_Fig2.pdf", width = 3.5, height = 3)
Supplementary_Fig2
dev.off()

tiff("Figures/Supplementary_Fig2.tiff", width = 3.5, height = 3, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig2
dev.off()


#### Supplementary Figure 3 ####

Supplementary_Fig3_tree <- ggtree(tree, branch.length='branch.length', size = 1) + 
  geom_treescale(x=0, y=13, fontsize = 2, linesize = 1)

Supplementary_Fig3_tree <- Supplementary_Fig3_tree %<+% species_data + 
  geom_tiplab(aes(label = factor(Name)), size=3.5)+
  geom_nodelab(size=2.5,
               hjust = 1.5,
               vjust=-0.5) +
  theme(legend.position = "")


Supplementary_Fig3 = Supplementary_Fig3_tree + xlim(0,2) +
  ggplot(subset(tRNA_copy_number, codon_counts > 0 & Aminoacid != "Ter"), aes(x=AAAnticodon, y=species)) +
  geom_point(alpha=0.75, shape=21, aes(size=codon_counts, fill=codon_counts))+
  scale_fill_viridis_b(option = "H", begin = 0.15, end = 1) +
  scale_size(range = c(1,8))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  theme(legend.position = "right")+
  scale_y_discrete(limits = rev(c("Apiotrichum_mycotoxinovorans_ACCC_20271",
                                  "Apiotrichum_mycotoxinovorans_CICC_1454",
                                  "Apiotrichum_mycotoxinivorans_GMU1709",
                                  "Apiotrichum_veenhuisii_JCM_10691",
                                  "Apiotrichum_akiyoshidainum_HP2023",
                                  "Apiotrichum_laibachii_JCM_2947",
                                  "Apiotrichum_gracile_JCM_10018",
                                  "Apiotrichum_domesticum_JCM_9580",
                                  "Apiotrichum_montevideense_JCM_9937",
                                  "Apiotrichum_brassicae_JCM_1599",
                                  "Apiotrichum_siamense_L8in5",
                                  "Apiotrichum_porosum_JCM_1458",
                                  "Apiotrichum_porosum_DSM_27194",
                                  "Apiotrichum_gamsii_JCM_9941",
                                  "Pascua_guehoae_Phaff_60_59",
                                  "Pascua_guehoae_JCM_10690",
                                  "Prillingera_fragicola_JCM_1530",
                                  "Cutaneotrichosporon_oleaginosus_IBC0246",
                                  "Cutaneotrichosporon_oleaginosum_ATCC_20509_reseq",
                                  "Cutaneotrichosporon_oleaginosus_ATCC_20508",
                                  "Cutaneotrichosporon_cutaneum_JCM_1462",
                                  "Cutaneotrichosporon_dermatis_JCM_11170",
                                  "Cutaneotrichosporon_dermatis_ATCC_204094",
                                  "Cutaneotrichosporon_arboriformis_JCM_14201",
                                  "Cutaneotrichosporon_cyanovorans_JCM_31833",
                                  "Cutaneotrichosporon_curvatus_JCM_1532",
                                  "Cutaneotrichosporon_daszewskae_JCM_11166",
                                  "Trichosporon_asahii_N5_275_008G1",
                                  "Trichosporon_asahii_CBS_8904",
                                  "Trichosporon_asahii_ATCC_201110",
                                  "Trichosporon_asahii_JCM_2466_CBS_2479_reseq",
                                  "Trichosporon_faecale_JCM_2941",
                                  "Trichosporon_inkin_JCM_9195_ATCC_18020_reseq",
                                  "Vanrija_humicola_JCM_1457",
                                  "Vanrija_humicola_CBS_4282",
                                  "Vanrija_humicola_UJ1",
                                  "Vanrija_humicola_ATCC_9949",
                                  "Vanrija_pseudolonga_DUCC4014",
                                  "Haglerozyma_chiarellii_ATCC_MYA-4694",
                                  "Takashimella_koratensis_JCM_12878",
                                  "Takashimella_tepidaria_JCM_11965",
                                  "Cryptococcus_floricola_DSM_27421",
                                  "Cryptococcus_amylolentus_CBS_6039",
                                  "Cryptococcus_deneoformans_JEC21",
                                  "Cryptococcus_gattii_WM276")))


pdf("Figures/Supplementary_Fig3.pdf", width = 15, height = 8)
Supplementary_Fig3
dev.off()

tiff("Figures/Supplementary_Fig3.tiff", width = 15, height = 8, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig3
dev.off()



#### Supplementary Figure 4 ####

Blomk = data.frame(matrix(ncol = 4, nrow = 0))
colnames(Blomk) <- c("Anticodon","K", "P", "dataset")


BlomK_data = subset(tRNA_copy_number, AACodon != "Asn(AAT)" &
                      AACodon != "Asp(GAT)" &
                      AACodon != "Cys(TGT)" &
                      AACodon != "His(CAT)" &
                      AACodon != "Pro(CCC)" &
                      AACodon != "Ser(TCC)" &
                      AACodon != "Tyr(TAT)" &
                      AACodon != "Val(GTC)" &
                      Aminoacid != "Ter" &
                      Order == "Trichosporonales")


for (i in unique(BlomK_data$AAAnticodon)) {
  anticodon=i
  test1 = subset(BlomK_data, AAAnticodon == i)
  rownames(test1) = test1$species
  test1 = subset(test1, select = c(species, codon_counts ))
  names(test1)=c("species", "codon_counts")
  
  K = phylosig(tree, setNames(test1[,"codon_counts"], rownames(test1)), method="K", test=TRUE, nsim=100)
  Blomk = rbind(Blomk, paste(c(i,K$K, K$P, "tRNA counts")))
}

colnames(Blomk) <- c("Anticodon","K", "P", "dataset")

Supplementary_Fig4 = ggplot(subset(Blomk, Blomk$P <=0.05 & Blomk$dataset=="tRNA counts"), aes(x=Anticodon, y=as.numeric(K)))+
  geom_point(size=2, aes(colour = cut(as.numeric(K), c(-Inf, 0.98, 1.02, Inf))))+
  geom_hline(yintercept=1)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        panel.grid = element_line(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_blank())+
  ylab("Blomberg's K") +
  xlab("tRNA gene") +
  scale_color_manual(name = "BM evolution",
                     values = c("(-Inf,0.98]" = "blue",
                                "(0.98,1.02]" = "black",
                                "(1.02, Inf]" = "red"),
                     breaks = c("(1.02, Inf]", "(0.98,1.02]", "(-Inf,0.98]"),
                     labels = c("More than expected", "As expected", "Less than expected"))



pdf("Figures/Supplementary_Fig4.pdf", height = 4, width = 7)
Supplementary_Fig4
dev.off()

tiff("Figures/Supplementary_Fig4.tiff", height = 4, width = 7, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig4
dev.off()



#### Supplementary Figure 5 ####

# Phylogenetic distance
pairwise_phylodistance = cophenetic.phylo(tree)
pairwise_phylodistance = as.data.frame(pairwise_phylodistance)

# Sort rows and columns alphabetically
pairwise_phylodistance = pairwise_phylodistance[order(row.names(pairwise_phylodistance)), ]
pairwise_phylodistance = pairwise_phylodistance[ , order(colnames(pairwise_phylodistance))]

# Extract triangles
pairwise_phylodistance[upper.tri(pairwise_phylodistance, diag = TRUE)] <- NA

# Melt
pairwise_phylodistance_melt <- melt(as.matrix(pairwise_phylodistance), na.rm = TRUE)
colnames(pairwise_phylodistance_melt) = c("species1", "species2", "phylodistance")


# tRNA structure correlation
tRNA_matrix = t(acast(subset(tRNA_copy_number, Order == "Trichosporonales"), species~AACodon, value.var="codon_counts", sum))
tRNA_matrix = as.data.frame(tRNA_matrix)


# Sort rows and columns alphabetically
tRNA_matrix = tRNA_matrix[order(row.names(tRNA_matrix)), ]
tRNA_matrix = tRNA_matrix[ , order(colnames(tRNA_matrix))]


# Pairwise correlations 
corr_tRNA = corr.test(tRNA_matrix)

# Retrieve R and P
corr_tRNA_r = corr_tRNA$r
corr_tRNA_p = corr_tRNA$p


# Extract triangles
corr_tRNA_r[upper.tri(corr_tRNA_r, diag = TRUE)] <- NA
corr_tRNA_p[upper.tri(corr_tRNA_p, diag = TRUE)] <- NA

corr_tRNA_melt_r <- melt(corr_tRNA_r, na.rm = TRUE)
corr_tRNA_melt_p <- melt(corr_tRNA_p, na.rm = TRUE)

corr_tRNA_melt = merge(corr_tRNA_melt_r, corr_tRNA_melt_p, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2")) 
colnames(corr_tRNA_melt) = c("species1", "species2", "R", "P")


# Merge data
corr_tRNA_phylodistance = merge(corr_tRNA_melt, pairwise_phylodistance_melt, by.x = c("species1", "species2") , by.y = c("species1", "species2"))

# Remove non.significant correlations
corr_tRNA_phylodistance = subset(corr_tRNA_phylodistance, P < 0.05)


# Calculate Bray-Curtis dissimilarity 
dissim = vegdist(t(tRNA_matrix), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
dissim = as.matrix(dissim)

dissim[upper.tri(dissim, diag = TRUE)] <- NA

dissim_melt <- melt(as.matrix(dissim), na.rm = TRUE)
colnames(dissim_melt) = c("species1", "species2", "dissimilarity")

dissim_melt$similarity = 1 - dissim_melt$dissimilarity


# Merge data
corr_tRNA_phylodistance = merge(corr_tRNA_phylodistance, dissim_melt[, c("species1", "species2", "similarity")], by.x = c("species1", "species2") , by.y = c("species1", "species2"))

corr_tRNA_phylodistance = merge(corr_tRNA_phylodistance, species_data[, c("id", "Genus", "Lifestyle")], by.x = "species1", by.y = "id")

corr_tRNA_phylodistance = merge(corr_tRNA_phylodistance, species_data[, c("id", "Genus",  "Lifestyle")], by.x = "species2", by.y = "id")

colnames(corr_tRNA_phylodistance) = c("species2", "species1", "R", "P", "phylodistance", "similarity", "Genus1", "Lifestyle.1", "Genus2", "Lifestyle.2")


# Create categories based on genus or lifestyles
corr_tRNA_phylodistance$category = ifelse(
  corr_tRNA_phylodistance$Genus1 == corr_tRNA_phylodistance$Genus2, "1. Same genus", 
  ifelse(
    corr_tRNA_phylodistance$Genus1 != corr_tRNA_phylodistance$Genus2, "2. Different genera", "error"))


corr_tRNA_phylodistance$lifestyle = ifelse(
  corr_tRNA_phylodistance$Lifestyle.1 == "Clinical" & corr_tRNA_phylodistance$Lifestyle.2 == "Clinical", "1. Clinical", 
  ifelse(
    corr_tRNA_phylodistance$Lifestyle.1 == "Environmental" & 
      corr_tRNA_phylodistance$Lifestyle.2 == "Environmental", "2. Environmental", 
    "3. Different lifestyles"))


# Remove Takashimella from dataset
corr_tRNA_phylodistance_subset = corr_tRNA_phylodistance
corr_tRNA_phylodistance_subset = corr_tRNA_phylodistance_subset[!grepl("Takashimella", corr_tRNA_phylodistance_subset$species1),]  
corr_tRNA_phylodistance_subset = corr_tRNA_phylodistance_subset[!grepl("Takashimella", corr_tRNA_phylodistance_subset$species2),] 
corr_tRNA_phylodistance_subset = subset(corr_tRNA_phylodistance_subset, R >= 0.8)


Supplementary_Fig5_R = ggplot(corr_tRNA_phylodistance_subset, aes(x=phylodistance, y=R, colour = category))+
  geom_point(alpha = 0.7, size = 1)+
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, aes(colour = category), size = 1.3) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', after_stat(r.squared), after_stat(p.value))),
                  label.x = 'left', label.y = 'bottom', size = 3)+
  xlab("Phylogenetic distance")+
  ylab("Correlation coeficient (R)")+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom") +
  scale_colour_manual(values = c("#9467BD", "#2CA02C"))


Supplementary_Fig5_sim = ggplot(corr_tRNA_phylodistance_subset, aes(x=phylodistance, y=similarity, colour = category))+
  geom_point(alpha = 0.7, size = 1)+
  geom_smooth(method="lm", se=FALSE, fullrange=T, level=0.95, aes(colour = category), size = 1.3) +
  stat_fit_glance(method = 'lm', 
                  geom = "text_npc",
                  method.args = list(formula = y ~ x),
                  aes(label = sprintf('R² = %.4f, P = %.2g', after_stat(r.squared), after_stat(p.value))),
                  label.x = 'left', label.y = 'bottom', size = 3)+
  xlab("Phylogenetic distance")+
  ylab("Bray-Curtis similarity (%)")+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom") +
  scale_colour_manual(values = c("#9467BD", "#2CA02C"))

pdf("Figures/Supplementary_Fig5_raw.pdf", height = 4, width = 8)
ggarrange(Supplementary_Fig5_R, Supplementary_Fig5_sim, align = "hv", common.legend = TRUE, legend = "bottom")
dev.off()

tiff("Figures/Supplementary_Fig5_raw.tiff", height = 4, width = 8, units = "in", compression = "lzw+p", res = 360)
ggarrange(Supplementary_Fig5_R, Supplementary_Fig5_sim, align = "hv", common.legend = TRUE, legend = "bottom")
dev.off()



#### Supplementary Figure 6 ####

# Format tRNA dataset
tRNA_copy_number$codon_counts_group = ifelse(tRNA_copy_number$codon_counts>101, "101-105", 
                                      ifelse(tRNA_copy_number$codon_counts<=100 & tRNA_copy_number$codon_counts>=96, "96-100", 
                                      ifelse(tRNA_copy_number$codon_counts<=95 & tRNA_copy_number$codon_counts>=91, "91-95",
                                      ifelse(tRNA_copy_number$codon_counts<=90 & tRNA_copy_number$codon_counts>=86, "86-90", 
                                      ifelse(tRNA_copy_number$codon_counts<=85 & tRNA_copy_number$codon_counts>=81, "81-85",
                                      ifelse(tRNA_copy_number$codon_counts<=80 & tRNA_copy_number$codon_counts>=76, "76-80",
                                      ifelse(tRNA_copy_number$codon_counts<=75 & tRNA_copy_number$codon_counts>=71, "71-75",
                                      ifelse(tRNA_copy_number$codon_counts<=70 & tRNA_copy_number$codon_counts>=66, "66-70",
                                      ifelse(tRNA_copy_number$codon_counts<=65 & tRNA_copy_number$codon_counts>=61, "61-65",
                                      ifelse(tRNA_copy_number$codon_counts<=60 & tRNA_copy_number$codon_counts>=56, "56-60",
                                      ifelse(tRNA_copy_number$codon_counts<=55 & tRNA_copy_number$codon_counts>=51, "51-55",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=50 & tRNA_copy_number$codon_counts>=46, "46-50",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=45 & tRNA_copy_number$codon_counts>=41, "41-45",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=40 & tRNA_copy_number$codon_counts>=36, "36-40",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=35 & tRNA_copy_number$codon_counts>=31, "31-35",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=30 & tRNA_copy_number$codon_counts>=26, "26-30",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=25 & tRNA_copy_number$codon_counts>=21, "21-25",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=20 & tRNA_copy_number$codon_counts>=16, "16-20",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=15 & tRNA_copy_number$codon_counts>=11, "11-15",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=10 & tRNA_copy_number$codon_counts>=6, "6-10",                                         
                                      ifelse(tRNA_copy_number$codon_counts<=5 & tRNA_copy_number$codon_counts>=2, "2-5", 
                                      ifelse(tRNA_copy_number$codon_counts<=1 &tRNA_copy_number$codon_counts>=0, NA, "error")
                                      )))))))))))))))))))))


Supplementary_Fig6_A = ggplot(tRNA_copy_number, aes(x=codon_counts, y=Mean_Distance))+
  geom_point(size = 0.5)+
  xlab("tRNA gene copy number (N)")+
  ylab("Intragenomic mean genetic distance") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(colour="black"),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black"))+
  scale_x_continuous(limits = c(0,105), breaks = seq(0, 105, 5), expand = c(0.01,0.01)) +
  labs(tag = "A")

Supplementary_Fig6_B = ggplot(tRNA_copy_number, aes(x=codon_counts_group, y=Mean_Distance))+
  geom_boxplot(fill="grey", outlier.size = 0.5, colour = "black", lwd=0.4)+
  xlab("tRNA gene copy number (N)")+
  ylab("Intragenomic mean genetic distance") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(colour="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black")) +
  scale_x_discrete(limits = c("2-5",
                              "6-10",
                              "11-15",
                              "16-20",
                              "21-25",
                              "26-30",
                              "31-35",
                              "36-40",
                              "41-45",
                              "46-50",
                              "51-55",
                              "56-60",
                              "61-65",
                              "66-70",
                              "71-75",
                              "76-80",
                              "81-85",
                              "86-90",
                              "91-95",
                              "96-100",
                              "101-105")) +
  labs(tag = "B")


pdf("Figures/Supplementary_Fig6.pdf", height = 7, width = 6.8)
plot_grid(Supplementary_Fig6_A, Supplementary_Fig6_B, align = "v", rel_widths = c(1, 1), ncol = 1)
dev.off()

tiff("Figures/Supplementary_Fig6.tiff", height = 7, width = 6.8, units = "in", compression = "lzw+p", res = 360)
plot_grid(Supplementary_Fig6_A, Supplementary_Fig6_B, align = "v", rel_widths = c(1, 1), ncol = 1)
dev.off()



#### Supplementary Figure 7 ####

Supplementary_Fig7 = ggplot(tRNA_genetic_distance, aes(x = Anticodon, y = Distance))+
  geom_violin(data = tRNA_genetic_distance, aes(fill = Aminoacid), scale = "width", colour = "black", linewidth = 0.2)+
  ylab("Genetic distance")+
  xlab("tRNA gene")+
  theme_bw()+
  theme(legend.position = "",
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(colour = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 7),
        strip.background = element_rect(colour = "black"),
        strip.text = element_text(colour="black", size = 10)) +
  facet_wrap(. ~ Aminoacid, scales = "free_x", ncol = 5)


pdf("Figures/Supplementary_Fig7.pdf", height = 7, width = 8)
Supplementary_Fig7
dev.off()

tiff("Figures/Supplementary_Fig7.tiff", height = 7, width = 8, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig7
dev.off()



#### Supplementary Figure 8 ####

Supplementary_Fig8 = ggplot(na.omit(subset(tRNA_copy_number, Aminoacid != "Ter" & Order == "Trichosporonales")), aes(x = reorder(AAAnticodon, Mean_Distance, median, na.rm = TRUE), y = Mean_Distance, drop = TRUE))+
  geom_boxplot(fill = "darkgrey", colour = "black")+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  xlab("tRNA gene") +
  ylab("Mean distance")


pdf("Figures/Supplementary_Fig8.pdf", height = 5, width = 7)
Supplementary_Fig8
dev.off()

tiff("Figures/Supplementary_Fig8.tiff", height = 5, width = 7, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig8
dev.off()




#### Supplementary Figure 9 ####

Supplementary_Fig9 = ggplot(na.omit(subset(tRNA_copy_number, Aminoacid != "Ter" & Order == "Trichosporonales")), aes(x = reorder(species, Mean_Distance, FUN = median), y = Mean_Distance))+
  geom_boxplot(fill = "darkgrey", colour = "black")+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  xlab("Species") +
  ylab("Mean distance")


pdf("Figures/Supplementary_Fig9.pdf", height = 7, width = 7)
Supplementary_Fig9
dev.off()

tiff("Figures/Supplementary_Fig9.tiff", height = 7, width = 7, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig9
dev.off()



#### Supplementary Figure 10 ####

rscu = read.delim("Source_files/Source_Cryptococcus_RSCU.txt")

rscu$species = factor(rscu$species, levels = c("Cryptococcus amylolentus CBS 6039", "Cryptococcus floricola DSM 27421",
                                               "Cryptococcus deneoformans JEC21", "Cryptococcus gattii WM276") ,
                      labels = c("Cryptococcus amylolentus CBS 6039" = "C. amylolentus CBS 6039", 
                                 "Cryptococcus floricola DSM 27421" = "C. floricola DSM 27421",
                                 "Cryptococcus deneoformans JEC21" = "C. deneoformans JEC21",
                                 "Cryptococcus gattii WM276" = "C. gattii WM276"))  

base_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")

num_codons <- rscu %>% pull(Codon) %>% unique() %>% length()
extended_palette <- rep(base_palette, length.out = num_codons)


Supplementary_Fig10 = ggplot(rscu, aes(x = value, y = fct_rev(dataset), fill = Codon)) +
  geom_col(position = "fill", colour = "black") +
  facet_grid(aa ~ species, 
             scales = "free") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.background = element_rect(fill = "#F0F0F0"),
        panel.spacing.x = unit(0.7, "lines"),
        strip.background.y = element_rect(colour = "black", fill = "#E0E0E0"),
        strip.text.x = element_text(face = "bold", size = 6),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8)) +
  scale_fill_manual(values = extended_palette) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_discrete(labels = c("RSCU genome" = "RSCU",
                              "tRNAs" = "tRNAs")) +
  ylab("") +
  xlab("Proportion")


pdf("Figures/Supplementary_Fig10.pdf", height = 9.44, width = 6.29)
Supplementary_Fig10
dev.off()


tiff("Figures/Supplementary_Fig10.tiff", height = 9.44, width = 6.29, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig10
dev.off()



#### Supplementary Figure 11 ####

all_functions_tAI = read.delim("Source_files/Source_all_functions_tAI.txt")


Supplementary_Fig11 = ggarrange(

  ggplot(subset(all_functions_tAI, S_norm > 0 & Process == "METABOLISM" ), aes(x = Lifestyle, y = S_norm, fill = Lifestyle)) +
    geom_boxplot() +
    geom_jitter(position=position_jitterdodge(jitter.width = 0.2), alpha = 0.7, aes(colour = Genus), size = 1) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.2, margin_top = 0.05,
                map_signif_level=function(p) sprintf("italic(P) == %.2g", p), parse = T, tip_length = 0.01, size = 0.4, textsize = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), 
                       breaks = seq(0.85, 1.53, 0.05), limits = c(0.85, 1.53)) +
    scale_colour_manual(values = c("yellowgreen", "darkred", "purple")) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_rect(fill = NA, colour = "black"),
          axis.text = element_text(colour = "black", size = 8),
          axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(size = 10)) +
    facet_grid(Process ~ Function, labeller = labeller(Function = label_wrap_gen(width = 22), Process = label_wrap_gen(width = 20))) +
    xlab("") +
    scale_fill_manual(name = "Lifestyle", labels = c("Opportunistic pathogens (OP)", "Saprotrophic (S)"),
                      values=c("Clinical" = "darkorange1","Environmental" = "dodgerblue")) +
    scale_x_discrete(labels=c("Clinical" = "OP",
                              "Environmental" = "S")) +
    ylab("Normalized S")
  
  ,
  
  ggplot(subset(all_functions_tAI, S_norm > 0 & Process == "CELLULAR PROCESSES AND SIGNALING" ), aes(x = Lifestyle, y = S_norm, fill = Lifestyle)) +
    geom_boxplot() +
    geom_jitter(position=position_jitterdodge(jitter.width = 0.2), alpha = 0.7, aes(colour = Genus), size = 1) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.2, margin_top = 0.05,
                map_signif_level=function(p) sprintf("italic(P) == %.2g", p), parse = T, tip_length = 0.01, size = 0.4, textsize = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), 
                       breaks = seq(0.2, 1.82, 0.1), limits = c(0.2, 1.82)) +
    scale_colour_manual(values = c("yellowgreen", "darkred", "purple")) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_rect(fill = NA, colour = "black"),
          axis.text = element_text(colour = "black", size = 8),
          axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(size = 10)) +
    facet_grid(Process ~ Function, labeller = labeller(Function = label_wrap_gen(width = 20), Process = label_wrap_gen(width = 20))) +
    xlab("") +
    scale_fill_manual(name = "Lifestyle", labels = c("Opportunistic pathogens (OP)", "Saprotrophic (S)"),
                      values=c("Clinical" = "darkorange1","Environmental" = "dodgerblue")) +
    scale_x_discrete(labels=c("Clinical" = "OP",
                              "Environmental" = "S")) +
    ylab("Normalized S")
  
  ,
  
  ggplot(subset(all_functions_tAI, S_norm > 0 & Process == "INFORMATION STORAGE AND PROCESSING" ), aes(x = Lifestyle, y = S_norm, fill = Lifestyle)) +
    geom_boxplot() +
    geom_jitter(position=position_jitterdodge(jitter.width = 0.2), alpha = 0.7, aes(colour = Genus), size = 1) +
    geom_signif(test="wilcox.test", comparisons = list(c("Clinical", "Environmental")), step_increase = 0.2, margin_top = 0.05,
                map_signif_level=function(p) sprintf("italic(P) == %.2g", p), parse = T, tip_length = 0.01, size = 0.4, textsize = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), 
                       breaks = seq(0.3, 1.78, 0.1), limits = c(0.3, 1.78)) +
    scale_colour_manual(values = c("yellowgreen", "darkred", "purple")) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_rect(fill = NA, colour = "black"),
          axis.text = element_text(colour = "black", size = 8),
          axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(size = 10)) +
    facet_grid(Process ~ Function, labeller = labeller(Function = label_wrap_gen(width = 25), Process = label_wrap_gen(width = 20))) +
    xlab("") +
    scale_fill_manual(name = "Lifestyle", labels = c("Opportunistic pathogens (OP)", "Saprotrophic (S)"),
                      values=c("Clinical" = "darkorange1","Environmental" = "dodgerblue")) +
    scale_x_discrete(labels=c("Clinical" = "OP",
                              "Environmental" = "S")) +
    ylab("Normalized S")
  
  ,
  
  legend = "bottom",
  common.legend = TRUE,
  nrow = 3
)



pdf("Figures/Supplementary_Fig11.pdf", width = 12, height = 10)
Supplementary_Fig11
dev.off()

tiff("Figures/Supplementary_Fig11.tiff", width = 12, height = 10, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig11
dev.off()



#### Supplementary Figure 12 ####

S_norm = read.delim("Source_files/Source_S_norm.txt")

Supplementary_Fig12 = ggplot(S_norm, aes(x = Dataset, y = value, fill = Dataset)) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), alpha = 0.6) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5)) +
  ylab("S value") +
  scale_fill_manual(values = c("carbs" = "dodgerblue",
                               "S_carbs_norm" = "#82CFFF",
                               "lipids" = "darkorange1",
                               "S_lipids_norm" = "#FF9E45",
                               "S_genome" = "#85B32D")) +
  scale_x_discrete(labels = c("lipids" = "S\nlipids",
                              "carbs" = "S\ncarbohydrates",
                              "S_genome" = "S\ngenome",
                              "S_carbs_norm" = "S norm\ncarbohydrates",
                              "S_lipids_norm" = "S norm\nlipids"),
                   limits = c("lipids",
                              "carbs",
                              "S_genome",
                              "S_lipids_norm",
                              "S_carbs_norm")) +
  geom_signif(test="wilcox.test", comparisons = list(c("lipids", "carbs"), c("S_carbs_norm", "S_lipids_norm")), step_increase = 0, margin_top = 0.1, 
              map_signif_level=function(p) sprintf("italic(P)-value == %.2g", p), parse = T) +
  facet_grid(. ~ factor(Lifestyle, levels = c("Environmental", "Clinical")), scales = "free_x", labeller = as_labeller(c("Clinical" = "Opportunistic pathogen", 
                                                                         "Environmental" = "Saprotrophic"))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), 
                     breaks = seq(0.40, 1.50, 0.1), limits = c(0.40, 1.50))



pdf("Figures/Supplementary_Fig12.pdf", height = 4, width = 7.5, useDingbats = FALSE)
Supplementary_Fig12
dev.off()

tiff("Figures/Supplementary_Fig12.tiff", height = 4, width = 7.5, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig12
dev.off()



#### Supplementary Figure 13 ####

Supplementary_Fig13 = ggplot(growth_temp, aes(x = Day, y = OD, fill = Lifestyle)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_nested(Temp ~ Genus + sample, scales = "free_x") +
  theme(axis.text.x = element_text(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text.x.top = element_text(angle = 0),
        strip.text.y = element_text(angle = 0),
        strip.background.x = element_rect(fill = NA, colour = "black"),
        strip.background.y = element_rect(fill = "grey95", colour = "black"),
        panel.background = element_rect(fill = NA, colour = "black")) +
  xlab("Day") +
  ylab(expression(paste(OD [" 600 nm"]))) +
  scale_fill_manual(name = "Lifestyle", values = c("Clinical" = "darkorange1", "Environmental" = "dodgerblue1"), 
                    labels=c("Clinical" = "Opportunistic pathogen", 
                             "Environmental" = "Saprotrophic"))


pdf("Figures/Supplementary_Fig13.pdf",  width = 9, height = 6)
Supplementary_Fig13
dev.off()

tiff("Figures/Supplementary_Fig13.tiff",  width = 9, height = 6, units = "in", compression = "lzw+p", res = 360)
Supplementary_Fig13
dev.off()