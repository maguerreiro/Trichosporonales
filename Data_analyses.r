
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

# Loading colour palette with 25 colours
colours25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown")
colours25b <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gray50", "steelblue4", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "skyblue2", "blue1", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown")

# Loading phylogeny
tree = read.tree("Source_files/Source_fig1.tree")

# Loading species data
species_data = read.delim("Source_files/Source_species_metadata.txt")
species_data = as.data.frame(species_data)

# Loading genome stats
genome_stats = read.delim("Source_files/Source_genome_stats.txt")

# Loading EggNOG annotations and functional categories
eggnog = read.delim2("Source_files/Source_EggNOG.annotations", header = 1)
functions = read.delim2("Source_files/Source_EggNOG_categories.txt", header = 1)

# Loading gene length
gene_length = read.delim("Source_files/Source_gene_length.txt", header = F)
colnames(gene_length) = c("gene", "gene.length")

# Loading tAI and S calculation for the whole genome
tai = read.delim()
S = read.delim()

# Loading growth results
growth_substrates = read.delim("Source_files/Source_growth_substrates.txt")

# Loading BUSCO results
busco = read.delim("Source_files/Source_BUSCO_basidiomycota.txt")











##### Format EggNOG results #####

# subsets eggnog to include only the OG category (narr_og_cat)
eggnog_cat = subset(eggnog, select = c(query_name, narr_og_cat))

# remove rows with no eggnog categories (-)
eggnog_cat <- subset(eggnog_cat, eggnog_cat$narr_og_cat!="-")







##### Process genome stats for plotting #####

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







##### ——— MAIN FIGURES ——— #####

#### Figure 1 ####

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



#### Figure 4 ####

Fig4 = ggarrange(
  
  ggplot(growth_substrates, 
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
  
  ggplot(growth_substrates, 
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



##### ——— EXTENDED DATA FIGURES ——— #####

#### Extended Data Figure 1 ####

# For plotting purposes, Cryptococcus has been labeled ZCryptococcus

busco = subset(busco, select = -c(Lineage, Complete))
busco = melt(busco)

Extended_Data_Fig1 = ggplot(busco, aes(y=Species, x=value, fill=variable)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = colours25, labels = c("Complete and single-copy BUSCOs (S)", "Complete and duplicated BUSCOs (D)", "Fragmented BUSCOs (F)", "Missing BUSCOs (M)"))+
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


