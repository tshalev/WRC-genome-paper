rm(list=ls())
library(reshape2)
library(ggsci)
library(tidyverse)
source("~/UBC/GSAT/PhD/WRC/r_scripts/publication_theme.r")

library(extrafont)
library(extrafontdb)

setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/S_lines_pop_gen/")

parent_vep <- read.table("../vep_1.3_annotation_final_compound_effects_removed.txt", header = T)
significant_snps <- read.table("significant_snps_chisq_holm_1e-5.txt")
diverse_maf <- read.table("../diverse_pop_snp_maf.txt", header = T)

# Getting outlier SNPs and adding MAF
S_selected_snps_vep <- parent_vep %>%
  filter(SNP %in% significant_snps$V1) %>% 
  arrange(SNP)

maf_selected_snps <- diverse_maf %>% 
  filter(SNP %in% S_selected_snps_vep$SNP) %>% 
  arrange(SNP)

S_selected_snps_vep_w_maf <- S_selected_snps_vep %>% 
  mutate(maf = maf_selected_snps$maf, minor_allele = maf_selected_snps$minor_allele)

# SNP consequences from vep file
SNP_consequences_parents <- parent_vep %>%
  # filter(!Uploaded_variation %in% S_selected_snps_vep$Uploaded_variation) %>% 
  select(Consequence) %>%
  group_by(Consequence) %>% 
  summarise(count_parents =n()) %>% 
  mutate(freq_parents = count_parents / sum(count_parents))

SNP_consequences_selected <- S_selected_snps_vep %>% 
  select(Consequence) %>%
  group_by(Consequence) %>% 
  summarise(count_selected =n()) %>% 
  mutate(freq_selected = count_selected / sum(count_selected))

SNP_consequences_all <- merge(SNP_consequences_parents, SNP_consequences_selected, all = T)
SNP_consequences_all[is.na(SNP_consequences_all)] <- 0

SNP_consequences_all_freq <- SNP_consequences_all %>% 
  select(Consequence, freq_parents, freq_selected) %>% 
  arrange(desc(freq_selected))

SNP_consequences_all_freq$Consequence <- factor(SNP_consequences_all_freq$Consequence, 
                                              levels = c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
                                                         "start_lost", "missense_variant", "splice_region_variant", "stop_retained_variant",
                                                         "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
                                                         "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"))

parent_vep$Consequence <- factor(parent_vep$Consequence, 
                                                levels = c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
                                                           "start_lost", "missense_variant", "splice_region_variant", "stop_retained_variant",
                                                           "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
                                                           "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"))

S_selected_snps_vep$Consequence <- factor(S_selected_snps_vep$Consequence, 
                                 levels = c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
                                            "start_lost", "missense_variant", "splice_region_variant", "stop_retained_variant",
                                            "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
                                            "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"))

melted <- melt(SNP_consequences_all_freq)
# melted$Consequence <- factor(melted$Consequence, levels = melted$Consequence[order(melted$value)])

vep <- ggplot(parent_vep, aes(x = Consequence)) +
  geom_histogram(stat = "count") +
  ylab("Count") +
  coord_flip() +
  theme_Publication() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  # scale_y_continuous(name = "Frequency", breaks = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), limits = c(0,0.3)) +
  scale_x_discrete(limits = rev(levels(melted$Consequence)), breaks = c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
                                                                        "start_lost", "missense_variant", "splice_region_variant", "stop_retained_variant",
                                                                        "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
                                                                        "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"),
                   labels = c("Splice acceptor variant", "Splice donor variant", "Stop gained", "Stop lost",
                                                                        "Start lost", "Missense variant", "Splice region variant", "Stop retained variant",
                                                                        "Synonymous variant", "5' UTR variant", "3' UTR variant", "Intron variant",
                                                                        "Upstream gene variant", "Downstream gene variant", "Intergenic variant")) 
  

ggplot(SNP_consequences_parents, aes(x = Consequence, y = freq_parents)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("vep.svg", vep, height = 6, width = 11)
ggsave("S_vep.svg", S_vep, height = 6, width = 11)

# Plot for Figure 5
vep_outlier <- ggplot(melted, aes (x = Consequence, y =value)) +
  geom_bar(aes(fill = variable), stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_Publication() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  scale_y_continuous(name = "Frequency", breaks = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), limits = c(0,0.3)) +
  scale_x_discrete(limits = rev(levels(melted$Consequence)), breaks = c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
                                                                         "start_lost", "missense_variant", "splice_region_variant", "stop_retained_variant",
                                                                         "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
                                                                         "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"),
                   labels = c("Splice acceptor variant", "Splice donor variant", "Stop gained", "Stop lost",
                                                                        "Start lost", "Missense variant", "Splice region variant", "Stop retained variant",
                                                                        "Synonymous variant", "5' UTR variant", "3' UTR variant", "Intron variant",
                                                                        "Upstream gene variant", "Downstream gene variant", "Intergenic variant")) +
  # geom_signif(comparisons = list(c("Missense variant", "5' UTR variant", "3' UTR variant", "Downstream gene variant", "Intergenic variant")),
  #             map_signif_level = T) +
  scale_fill_manual(breaks=c("freq_parents", "freq_selected"),
                    labels=c("All SNPs", "Outlier SNPs"),
                    values = pal_nejm()(8)) +
  guides(y_continuous = FALSE) +
  guides(x_discrete = FALSE) +
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.95))
  # geom_bracket(
  #   xmin = group1, xmax = group2,
  #   y.position = 0.05, label = "*",
  #   tip.length = 0.01, coord.flip = T)
  # geom_signif(stat="identity",
  #             data=data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),
  #                             y=c(0.25, 0.15), annotation=c("**", "NS")),
  #             aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
  
vep_outlier


ggsave("vep_outlier.svg", vep_outlier, height = 7, width = 10)

# Fisher Exact Test for overrepresentation of SNP consequences
table_maker_func <- function(x){
  matrix(c(SNP_consequences_all[x, 2], 
           (sum(SNP_consequences_all$count_parents) - SNP_consequences_all[x, 2]),
           SNP_consequences_all[x, 4],
           (sum(SNP_consequences_all$count_selected) - SNP_consequences_all[x, 4])),
         nrow = 2,
         dimnames = list(c(paste0(SNP_consequences_all[x, 1]), "Other"), 
                         c("All", "Selected")))
}  

for (i in 1:24){
  x <- table_maker_func(i)
  write.table(x, "contingency_tables.txt", append = T, col.names = T, row.names = T, quote = F, sep = "\t")
}

count_table <- SNP_consequences_all %>% 
  select(count_parents, count_selected)

row.names(count_table) = SNP_consequences_all$Consequence

chisq.test(count_table, simulate.p.value = T)

for(i in 1:15){
  x <- table_maker_func(i)
  y <- fisher.test(x)
  lapply(y, write, "FET_results_compound_effects_removed.txt", append = T)
}

f <- fisher.test(x)

S4_fate_table <- read.table("S_lines_pop_gen/S4_fate_table_significant.txt", header = T)
significant_snps <- S4_fate_table %>% 
  filter(p.adj.holm < 1e-05)

significant_fixed <- significant_snps %>% 
  filter(ref_fixed > not_fixed | alt_fixed > not_fixed)
  
significant_fixed_vep <- S_selected_snps_vep_w_maf %>% 
  filter(SNP %in% significant_fixed$SNP)

stop_start_splice_gained_lost_snps <- S_selected_snps_vep_w_maf %>% 
  filter(Consequence =="stop_gained" | Consequence == "stop_lost" | Consequence == "start_lost" | 
           Consequence == "splice_acceptor_variant" | Consequence == "splice_donor_variant") %>% 
  arrange(SNP)

S4_fate_start_stop_splice_gained_lost <- S4_fate_table_props_for_chisq %>% 
  filter(SNP %in% stop_start_splice_gained_lost_snps$SNP) %>% 
  arrange(SNP) %>% 
  mutate(maf = stop_start_splice_gained_lost_snps$maf, minor_allele = stop_start_splice_gained_lost_snps$minor_allele)


### JGI annotations ##########################

gene_annotations <- read.table("../gene.functions.txt", skip = 1, header = F,
                               sep = "\t")
gff <- read.table("../Tplicatav3.1c.gene_exons_dashes_removed_final.gff3")

## significant snps ###
significant_snp_annotations <- gene_annotations %>%
  filter(V1 %in% S_selected_snps_vep$Gene) %>% 
  rename(Gene = V1)

significant_snp_annotations <- merge(S_selected_snps_vep, significant_snp_annotations, by = "Gene", all = T)  

# write.table(significant_snp_annotations, "significant_snp_JGI_annotations.txt", quote = F, col.names = T, row.names = F)

### GO annotations ######
significant_snp_GO <- significant_snp_annotations %>% 
  filter(V3 == "GO")

table(significant_snp_GO$V4)
summary(factor(significant_snp_GO$V4))

table(significant_snp_GO$V2)

significant_snp_GO_small <- significant_snp_GO %>% 
  select(Consequence, V4)

tbl_summary(significant_snp_GO_small) %>% 
  as_flex_table()

ggplot(significant_snp_GO_small, aes(x = V4)) +
  geom_histogram(stat = "count")

# write_csv(significant_snp_GO, "significant_snps_GO.csv")

significant_snp_GO_one_snp_per_gene <- read.table("significant_snps_GO_one_snp_per_gene.csv", sep =",", header = T)
significant_snp_GO_one_snp_per_gene_small <- significant_snp_GO_one_snp_per_gene %>% 
  select(Consequence, V4)

tbl_summary(significant_snp_GO_one_snp_per_gene_small)

## snps het in all lines ###
het_snp_vep <- parent_vep %>%
  filter(SNP %in% rownames(het_snps)) %>% 
  arrange(SNP)

het_snp_annotations <- gene_annotations %>%
  filter(V1 %in% het_snp_vep$Gene) %>% 
  rename(Gene = V1)

het_snp_annotations <- merge(het_snp_vep, het_snp_annotations, by = "Gene", all = T)  

het_snp_GO <- het_snp_annotations %>% 
  filter(V3 == "GO")

table(het_snp_GO$V4)
ggplot(het_snp_GO, aes(x = V4)) +
  geom_histogram(stat = "count")

## all snps

all_snp_annotations <- gene_annotations %>%
  filter(V1 %in% parent_vep$Gene) %>% 
  rename(Gene = V1)

all_snp_annotations <- merge(parent_vep, all_snp_annotations, by = "Gene", all = T)  

all_snp_GO <- all_snp_annotations %>% 
  filter(V3 == "GO")

# write.table(all_snp_GO, "all_snps_GO.txt", quote = F, col.names = T, row.names = F)

all_snp_GO_one_snp_per_gene <- read.table("all_snps_GO_one_snp_per_gene.txt", sep = "\t", header = T)
write.csv(table(all_snp_GO_one_snp_per_gene$V2), "all_snp_GO_ID_table.csv")
write.csv(table(all_snp_GO_one_snp_per_gene$V4), "all_snp_GO_term_table.csv")

select_GO_counts <- read.table("GO_table_code.csv", sep = ",", header = T)
all_GO_counts <- read.table("all_snp_GO_ID_table.csv", sep = ",", header = T)

merged_GO_counts <- merge(all_GO_counts, select_GO_counts, by = "Var1", all = T)
merged_GO_counts[is.na(merged_CO_counts)] <- 0

table_maker_func <- function(x){
  matrix(c(merged_GO_counts[x, 4], 
           (sum(merged_GO_counts$Freq.x) - merged_GO_counts[x, 4]),
           merged_GO_counts[x, 7],
           (sum(merged_GO_counts$Freq.y) - merged_GO_counts[x, 7])),
         nrow = 2,
         dimnames = list(c(paste0(merged_GO_counts[x, 1],"_",merged_GO_counts[x,3]), "Other"), 
                         c("All", "Selected")))
}  

for(i in 1:933){
  x <- table_maker_func(i)
  y <- fisher.test(x, alternative = "less")
  z <- data.frame(Name = row.names(x)[1], p = y$p.value, estimate = y$estimate)
  lapply(z, write, "FET_results_GO_categories.txt", append = T)
}
