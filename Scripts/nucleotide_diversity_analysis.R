library(rstatix)
library(gghighlight)
library(tidyverse)


### Reading in pi files ####
# Pi 5 Kb window for whole diverse pop, clusters (not reported)
all_snps_5k_window_one_pop_pi <- read.table("all_snps_5k_window_one_pop/all_snps_5k_window_one_pop_pi.txt", header = T)
all_snps_5k_window_clusters_pi <- read.table("all_snps_5k_window_clusters/all_snps_5k_window_clusters_pi.txt", header = T)

# Pi scaffolds for whole pop, clusters
all_snps_scaffolds_one_pop_pi <- read.table("all_snps_scaffolds_one_pop/all_snps_scaffolds_one_pop_pi.txt", header = T)
all_snps_scaffolds_clusters_pi <- read.table("all_snps_scaffolds_clusters/all_snps_scaffolds_clusters_pi.txt", header = T)

# Pi SCO genes for whole pop, clusters
all_snps_primary_genes_sco_one_pop_pi <- read.table("all_snps_primary_genes_sco_one_pop/all_snps_primary_genes_sco_one_pop_pi.txt", header = T)
all_snps_primary_genes_sco_clusters_pi <- read.table("all_snps_primary_genes_sco_clusters/all_snps_primary_genes_sco_clusters_pi.txt", header = T)

# Pi SCO genes for zero- and four-fold degenerate sites
all_snps_4_fold_gff_one_pop_pi <- read.table("all_snps_4-fold_gff_one_pop/all_snps_4-fold_gff_one_pop_pi.txt", header = T)
all_snps_0_fold_gff_one_pop_pi <- read.table("all_snps_0-fold_gff_one_pop/all_snps_0-fold_gff_one_pop_pi.txt", header = T) 

# DXY SCO genes clusters
all_snps_primary_genes_sco_clusters_dxy <- read.table("all_snps_primary_genes_sco_clusters/all_snps_primary_genes_sco_clusters_dxy.txt",header = T)


### Data summaries, filter for windows with less than 10 observations ####
all_snps_5k_window_one_pop_pi <- all_snps_5k_window_one_pop_pi %>% 
  filter(no_sites >= 10)

all_snps_5k_window_one_pop_pi %>% 
  group_by(pop) %>%
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_5k_window_clusters_pi <- all_snps_5k_window_clusters_pi %>% 
  filter(no_sites >= 10) 

all_snps_5k_window_clusters_pi %>% 
  group_by(pop) %>%
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_scaffolds_one_pop_pi <- all_snps_scaffolds_one_pop_pi %>% 
  filter(no_sites >= 10)
  
all_snps_scaffolds_one_pop_pi %>% 
  group_by(pop) %>%
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_scaffolds_clusters_pi <- all_snps_scaffolds_clusters_pi %>% 
  filter(no_sites >= 10)

all_snps_scaffolds_clusters_pi %>% 
  group_by(pop) %>%
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_primary_genes_sco_one_pop_pi <- all_snps_primary_genes_sco_one_pop_pi %>% 
  filter(no_sites >= 10)

all_snps_primary_genes_sco_one_pop_pi %>% 
  group_by(pop) %>%
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_primary_genes_sco_clusters_pi <- all_snps_primary_genes_sco_clusters_pi %>% 
  filter(no_sites >= 10)

all_snps_primary_genes_sco_clusters_pi %>% 
  group_by(pop) %>%
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_primary_genes_sco_clusters_dxy <- all_snps_primary_genes_sco_clusters_dxy %>% 
  filter(no_sites >= 10) %>% 
  unite(pop_comparison, c("pop1", "pop2"))

all_snps_primary_genes_sco_clusters_dxy %>% 
  group_by(pop_comparison) %>% 
  summarise(across(avg_dxy, list(mean = mean, median = median, sd = sd)))

all_snps_0_fold_gff_one_pop_pi <- all_snps_0_fold_gff_one_pop_pi %>% 
  filter(no_sites >= 10)

all_snps_0_fold_gff_one_pop_pi %>% 
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))

all_snps_4_fold_gff_one_pop_pi <- all_snps_4_fold_gff_one_pop_pi %>% 
  filter(no_sites >= 10)

all_snps_4_fold_gff_one_pop_pi %>% 
  summarise(across(avg_pi, list(mean = mean, median = median, sd = sd)))


### Plots #####

# Make scaffolds a factor
all_snps_primary_genes_sco_one_pop_pi$chromosome <- factor(all_snps_primary_genes_sco_one_pop_pi$chromosome)

# Top 1% of pi estimates
top_1 <- all_snps_primary_genes_sco_one_pop_pi[all_snps_primary_genes_sco_one_pop_pi$avg_pi > quantile(all_snps_primary_genes_sco_one_pop_pi$avg_pi,prob=1-1/100),]

# Figure 3A
ggplot(all_snps_primary_genes_sco_one_pop_pi, aes(x = chromosome, y = avg_pi, color = chromosome)) +
  geom_point(color = "#BC3C29FF") +
  # scale_color_manual(values = rep(c("black", "gray53"), 2227)) +
  theme_Publication() +
  theme(legend.position = "none") +
  xlab("Scaffold") +
  ylab(expression(bold(pi))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  gghighlight(avg_pi > quantile(avg_pi,prob=1-1/100),
              unhighlighted_params = list(colour = "gray53", alpha = 0.5))

ggsave("all_pi_dotplot.tiff", dpi = 300, height = 6, width = 19)

# Figure 3B
ggplot(all_snps_primary_genes_sco_one_pop_pi, aes(x = pop, y = avg_pi)) +
  geom_violin() +
  geom_boxplot(outlier.alpha = 0.3, width = 0.15) +
  theme_Publication() +
  xlab("Diverse population") +
  ylab(expression(bold(pi))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

ggsave("all_pi_violinplot.tiff", dpi = 300, height = 6, width = 6)

# Figure 3C
ggplot(all_snps_primary_genes_sco_clusters_dxy, aes(x = pop_comparison, y = avg_dxy)) +
  geom_violin() +
  geom_boxplot(outlier.alpha = 0.3, width = 0.15) +
  theme_Publication() +
  ylab(expression(bold(italic(d)[XY]))) +
  scale_x_discrete(name="Cluster",
                   breaks=c("Northern_Coastal_Central", "Central_Southern_Interior", "Northern_Coastal_Southern_Interior"),
                   labels=c("Northern-Coastal/\nCentral", "Central/\nSouthern-Interior", "Northern-Coastal/\nSouthern-Interior"))

ggsave("all_dxy_violinplot.tiff", dpi = 300, height = 6, width = 7)

# Combine zero- and four-fold estimates
all_snps_0_4_gff <- rbind(data.frame(pi = all_snps_0_fold_gff_one_pop_pi$avg_pi, degen = rep("zero", 2772)), data.frame(pi = all_snps_4_fold_gff_one_pop_pi$avg_pi, degen = rep("four", 2179)))

# Figure 3D
ggplot(all_snps_0_4_gff, aes(x = degen, y = pi)) +
  geom_violin() +
  geom_boxplot(outlier.alpha = 0.3, width = 0.15) +
  theme_Publication() +
  ylab(expression(bold(pi))) +
  scale_x_discrete(name="Sites",
                   breaks=c("four", "zero"),
                   labels=c(expression(paste(bold(pi[4]))), expression(paste(bold(pi[0])))))
ggsave("0-4-fold_pi.tiff", height = 6, width = 6)

### Annotations for top 1% ####  
# Load annotations
gene_annotations <- read.delim("gene.functions_dashes_removed.txt", skip = 1, header = F,
                               sep = "\t")

# Load GFF
gff <- read.table("Tplicatav3.1c.primaryTrs_dashes_removed_sorted.gff3")

# Create variable for "window"
gff_genes <- gff %>% 
  filter(V3 == "gene") %>% 
  unite(window, c(V4, V5))

# Simplify gene names
gff_genes[, "V9"] <- str_extract_all(gff_genes[, "V9"], "Name=Thupl.*", simplify = T)[,1]
gff_genes[, "V9"] <- str_extract_all(gff_genes[, "V9"], "Thupl.*", simplify = TRUE)[,1]

# Match GFF windows to pi windows
top_1_bed_window_adjusted <- top_1 %>% 
  mutate(window_pos_1_adjusted = (window_pos_1 + 1)) %>% 
  unite(window, c(window_pos_1_adjusted, window_pos_2))
  
top_1_gff <- gff_genes %>% 
  filter(window %in% top_1_bed_window_adjusted$window) %>% 
  rename(Gene = V9)

top_1_gene_annotations <- gene_annotations %>% 
  rename(Gene = V1) %>% 
  filter(Gene %in% top_1_gff$Gene)
  
top_1_annotations <- merge(top_1_gff, top_1_gene_annotations, by = "Gene", all = T) 

write.csv(top_1_annotations, "top_1_annotations.csv")

### Cluster differences - statistical tests ####

# Pi
NC_cluster_pi <- all_snps_primary_genes_sco_clusters_pi %>% 
  filter(pop == "Northern_Coastal")

C_cluster_pi <- all_snps_primary_genes_sco_clusters_pi %>% 
  filter(pop == "Central")

SI_cluster_pi <- all_snps_primary_genes_sco_clusters_pi %>% 
  filter(pop == "Southern_Interior")

mean(NC_cluster_pi$avg_pi, na.rm = T)
mean(C_cluster_pi$avg_pi, na.rm = T)
mean(SI_cluster_pi$avg_pi, na.rm = T)

clusters_df <- data.frame(NC_cluster_pi$avg_pi, C_cluster_pi$avg_pi, SI_cluster_pi$avg_pi)

wilcox.test(NC_cluster_pi$avg_pi, C_cluster_pi$avg_pi)
wilcox.test(NC_cluster_pi$avg_pi, SI_cluster_pi$avg_pi)
wilcox.test(SI_cluster_pi$avg_pi, C_cluster_pi$avg_pi)

kruskal.test(clusters_df)

# dxy
NC_C_cluster_dxy <- all_snps_primary_genes_sco_clusters_dxy %>% 
  filter(pop_comparison == "Northern_Coastal_Central")

C_SI_cluster_dxy <- all_snps_primary_genes_sco_clusters_dxy %>% 
  filter(pop_comparison == "Central_Southern_Interior")

NC_SI_cluster_dxy <- all_snps_primary_genes_sco_clusters_dxy %>% 
  filter(pop_comparison == "Northern_Coastal_Southern_Interior")


mean(NC_C_cluster_dxy$avg_dxy, na.rm = T)
mean(C_SI_cluster_dxy$avg_dxy, na.rm = T)
mean(NC_SI_cluster_dxy$avg_dxy, na.rm = T)

clusters_df_dxy <- data.frame(NC_C_cluster_dxy$avg_dxy, C_SI_cluster_dxy$avg_dxy, NC_SI_cluster_dxy$avg_dxy)

kruskal.test(clusters_df_dxy)

# Pi
NC_cluster_pi <- all_snps_all_primary_genes_clusters_pi %>% 
  filter(pop == "Northern_Coastal")

C_cluster_pi <- all_snps_all_primary_genes_clusters_pi %>% 
  filter(pop == "Central")

SI_cluster_pi <- all_snps_all_primary_genes_clusters_pi %>% 
  filter(pop == "Southern_Interior")

mean(NC_cluster_pi$avg_pi, na.rm = T)
mean(C_cluster_pi$avg_pi, na.rm = T)
mean(SI_cluster_pi$avg_pi, na.rm = T)

clusters_df <- data.frame(NC_cluster_pi$avg_pi, C_cluster_pi$avg_pi, SI_cluster_pi$avg_pi)

wilcox.test(NC_cluster_pi$avg_pi, C_cluster_pi$avg_pi)
wilcox.test(NC_cluster_pi$avg_pi, SI_cluster_pi$avg_pi)
wilcox.test(SI_cluster_pi$avg_pi, C_cluster_pi$avg_pi)

kruskal.test(clusters_df)

# dxy
NC_C_cluster_dxy <- all_snps_all_primary_genes_clusters_dxy %>% 
  filter(pop_comparison == "Northern_Coastal_Central")

C_SI_cluster_dxy <- all_snps_all_primary_genes_clusters_dxy %>% 
  filter(pop_comparison == "Central_Southern_Interior")

NC_SI_cluster_dxy <- all_snps_all_primary_genes_clusters_dxy %>% 
  filter(pop_comparison == "Northern_Coastal_Southern_Interior")


mean(NC_C_cluster_dxy$avg_dxy, na.rm = T)
mean(C_SI_cluster_dxy$avg_dxy, na.rm = T)
mean(NC_SI_cluster_dxy$avg_dxy, na.rm = T)

clusters_df_dxy <- data.frame(NC_C_cluster_dxy$avg_dxy, C_SI_cluster_dxy$avg_dxy, NC_SI_cluster_dxy$avg_dxy)

kruskal.test(clusters_df_dxy)
