rm(list = ls())

## Script to obtain population genomics statistics using the PopGenome package for R ##

## Preliminary steps
# Load library PopGenome
library(PopGenome)
library(reshape2)
library(tidyverse)

setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/")

path <- "vcf_split_snps_in_genes_final_annotation/"

# Get list of scaffold lengths in VCF
v3_scaffold_lengths <- read.table("genome_scaffold_lengths_dashes_removed.txt")
# scaffolds_with_genes <- read.table("scaffolds_with_genes.txt")

# ggplot(v3_scaffold_lengths_in_vcf, aes(x = V2)) +
  # geom_histogram(bins = 100)

# write.table(v3_scaffold_lengths_in_vcf, "scaffold_lengths_scaffolds_with_genes.txt", col.names = F, row.names = F, quote = F)

# Read in scaffold names
scaffold_names <- scan("scaffolds_with_genes_final_annotation.txt", what = "factor")

# scaffold_names <- factor(scaffold_names)

v3_scaffold_lengths_in_vcf <- v3_scaffold_lengths %>% 
  filter(V1 %in% scaffold_names)

# write.table(v3_scaffold_lengths_in_vcf, "scaffold_lengths_scaffolds_in_gff_and_vcf_v3.1c.txt", col.names = F, row.names = F, quote = F)

# scaffold_names_new <- read.table("scaffolds_in_gff_and_vcf_v3.1c.txt")
# scaffold_names_old <- read.table("scaffolds_in_gff_and_vcf.txt")
# 
# filter(scaffold_names_old,!(V1 %in% scaffold_names_new$V1))
# filter(scaffold_names_new,!(V1 %in% scaffold_names_old$V1))

# read in GFF as table
gff <- read.table("Tplicatav3.1c.gene_exons_dashes_removed_final.gff3")

# Filtering GFF file to match VCF
gff_scaffolds_with_genes <- gff %>% 
  filter(V1 %in% scaffold_names)

# write.table(gff_scaffolds_with_genes, "gff_nonsynonymous_scaffolds_v3.1c_final.gff3", quote = F, col.names = F, row.names = F, sep = "\t")
# GFF_split_into_scaffolds("gff_nonsynonymous_scaffolds_v3.1c_final.gff3", "gff_nonsynonymous_scaffolds_v3.1c_final")

# # Create list of matching scaffolds between VCF and GFF
# gff_scaffold_names_with_genes <- gff_scaffolds_with_genes %>% 
#   select(V1) %>% 
#   filter(V1 %in% scaffold_names)
# 
# gff_genes <- gff_scaffolds_with_genes %>% 
#   filter(V3 == "gene")
# 
# gff_scaffolds_with_genes_unique <- unique(gff_scaffolds_with_genes$V1, na.rm = T)
# gff_scaffolds_with_genes_unique <- as.vector(gff_scaffolds_with_genes_unique)
# gff_scaffolds_with_genes_unique <- as.factor(gff_scaffolds_with_genes_unique)
# 
# 
# # write.table(gff_scaffolds_with_genes_unique, "unique_scaffold_names_with_genes_v3.1c.txt", col.names = F, row.names = F, quote = F)
# gff_scaffold_names_correct_order <- scan("unique_scaffold_names_with_genes_v3.1c.txt", what = "factor")
# # gff_scaffold_names_correct_order <- factor(gff_scaffold_names_correct_order)

# Scaffold length
# Read in scaffold lengths
# scaffold_lengths <- read.table("scaffold_lengths_scaffolds_in_gff_and_vcf_v3.1c.txt")
scaffold_lengths_with_genes <- v3_scaffold_lengths %>% 
  filter(V1 %in% scaffold_names) %>% 
  arrange(V1)

# write.table(scaffold_lengths_gff_vcf, "v3_scaffold_lengths_in_vcf_gff_1456_mislabeled_S_lines_removed.txt", quote = F, col.names = F, row.names = F)
# # 
# write.table(, "scaffold_lengths_scaffolds_with_genes_just_lengths_v3.1c.txt", quote = F, col.names = F, row.names = F)

# scaffold_lengths <- scan("scaffold_lengths_scaffolds_in_gff_and_vcf_v3.1c_just_lengths.txt")

scaffold_lengths <- as.numeric(scaffold_lengths_with_genes$V2)
length = c(scaffold_lengths)

# Define subpopulations
parents <- read.table("parent_samples.txt", header = F)
parents$V1 <- as.character(parents$V1)
FS <- read.table("F_samples_corrected.txt", header = F)
FS$V1 <- as.character(FS$V1)
S1 <- read.table("S1_samples_corrected.txt", header = F)
S1$V1 <- as.character(S1$V1)
S2 <- read.table("S2_samples_corrected.txt", header = F)
S2$V1 <- as.character(S2$V1)
S3 <- read.table("S3_samples_corrected.txt", header = F)
S3$V1 <- as.character(S3$V1)
S4 <- read.table("S4_samples_corrected.txt", header = F)
S4$V1 <- as.character(S4$V1)
S5 <- read.table("S5_samples_corrected.txt", header = F)
S5$V1 <- as.character(S5$V1)

# Define line subpopulations
Line_1 <- read.table("../pop_gen_v3_snps/Line_1_samples.012.indv", header = F, stringsAsFactors = F)
Line_6 <- read.table("../pop_gen_v3_snps/Line_6_samples.012.indv", header = F, stringsAsFactors = F)
Line_7 <- read.table("../pop_gen_v3_snps/Line_7_samples.012.indv", header = F, stringsAsFactors = F)
Line_8 <- read.table("../pop_gen_v3_snps/Line_8_samples.012.indv", header = F, stringsAsFactors = F)
Line_10 <- read.table("../pop_gen_v3_snps/Line_10_samples.012.indv", header = F, stringsAsFactors = F)
Line_13 <- read.table("../pop_gen_v3_snps/Line_13_samples.012.indv", header = F, stringsAsFactors = F)
Line_16 <- read.table("../pop_gen_v3_snps/Line_16_samples.012.indv", header = F, stringsAsFactors = F)
Line_17 <- read.table("../pop_gen_v3_snps/Line_17_samples.012.indv", header = F, stringsAsFactors = F)
Line_19 <- read.table("../pop_gen_v3_snps/Line_19_samples.012.indv", header = F, stringsAsFactors = F)
Line_20 <- read.table("../pop_gen_v3_snps/Line_20_samples.012.indv", header = F, stringsAsFactors = F)
Line_21 <- read.table("../pop_gen_v3_snps/Line_21_samples.012.indv", header = F, stringsAsFactors = F)
Line_23 <- read.table("../pop_gen_v3_snps/Line_23_samples.012.indv", header = F, stringsAsFactors = F)
Line_26 <- read.table("../pop_gen_v3_snps/Line_26_samples.012.indv", header = F, stringsAsFactors = F)
Line_27 <- read.table("../pop_gen_v3_snps/Line_27_samples.012.indv", header = F, stringsAsFactors = F)
Line_29 <- read.table("../pop_gen_v3_snps/Line_29_samples.012.indv", header = F, stringsAsFactors = F)

# Define generation in line subpopulations
Line_6_FS <- read.table("Line_6_FS_samples.txt", header = F, stringsAsFactors = F)
Line_6_S1 <- read.table("Line_6_S1_samples.txt", header = F, stringsAsFactors = F)
Line_6_S4 <- read.table("Line_6_S4_samples.txt", header = F, stringsAsFactors = F)
Line_7_FS <- read.table("Line_7_FS_samples.txt", header = F, stringsAsFactors = F)
Line_7_S1 <- read.table("Line_7_S1_samples.txt", header = F, stringsAsFactors = F)
Line_7_S4 <- read.table("Line_7_S4_samples.txt", header = F, stringsAsFactors = F)
Line_10_FS <- read.table("Line_10_FS_samples.txt", header = F, stringsAsFactors = F)
Line_10_S1 <- read.table("Line_10_S1_samples.txt", header = F, stringsAsFactors = F)
Line_10_S4 <- read.table("Line_10_S4_samples.txt", header = F, stringsAsFactors = F)
Line_13_FS <- read.table("Line_13_FS_samples.txt", header = F, stringsAsFactors = F)
Line_13_S1 <- read.table("Line_13_S1_samples.txt", header = F, stringsAsFactors = F)
Line_13_S4 <- read.table("Line_13_S4_samples.txt", header = F, stringsAsFactors = F)
Line_16_FS <- read.table("Line_16_FS_samples.txt", header = F, stringsAsFactors = F)
Line_16_S1 <- read.table("Line_16_S1_samples.txt", header = F, stringsAsFactors = F)
Line_16_S4 <- read.table("Line_16_S4_samples.txt", header = F, stringsAsFactors = F)
Line_17_FS <- read.table("Line_17_FS_samples.txt", header = F, stringsAsFactors = F)
Line_17_S1 <- read.table("Line_17_S1_samples.txt", header = F, stringsAsFactors = F)
Line_17_S4 <- read.table("Line_17_S4_samples.txt", header = F, stringsAsFactors = F)
Line_19_FS <- read.table("Line_19_FS_samples.txt", header = F, stringsAsFactors = F)
Line_19_S1 <- read.table("Line_19_S1_samples.txt", header = F, stringsAsFactors = F)
Line_19_S4 <- read.table("Line_19_S4_samples.txt", header = F, stringsAsFactors = F)

# Test with just one line
FS_1_1 <- "1_1-F"
S1_1_1 <- "111-S1"
S2_1_1 <- "111-2-S2"
S3_1_1 <- "111-22-S3"
S4_1_1 <- "111-221-S4"
S5_1_1 <- "111-221-S5"

# Split parents by location
VI_BC <- read.table("VI_BC.txt", header = F, stringsAsFactors = F)
COAST_BC <- read.table("COAST_BC.txt", header = F, stringsAsFactors = F)
HG_BC <- read.table("HG_BC.txt", header = F, stringsAsFactors = F)
INT_BC <- read.table("INT_BC.txt", header = F, stringsAsFactors = F)
CA_OR <- read.table("CA_OR.txt", header = F, stringsAsFactors = F)

## Analysis of diverse population altogether - synonymous non-synonymous
# Load VCF and GFF files for each scaffold, set subpopulations, split data by gene and calculate neutrality and linkage statistics 
for (i in 1:length(scaffold_names)){
  GENOME.class <- readVCF(paste(path, "VCF_", scaffold_names[i], ".recode.vcf.gz", sep = ""), numcols = 100000, 
                          tid = scaffold_names[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE, gffpath = paste(path, scaffold_names[i], sep = ""))
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, scaffold_names[i], sep = ""),
                                        # chr = scaffold_names[i], feature = "gene")
  gff_info <- get_gff_info(gff.file = paste(path, scaffold_names[i], sep = ""),
                           chr= scaffold_names[i], feature = "gene")
                           if (class(gff_info) == "matrix"){
                           gff_info <- list(gff_info)
                           }
  # Get synonymous/nonsynonymous sites
  GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste(path, scaffold_names[i], ".fa", sep = ""))
  # Set populations
  # GENOME.class <- set.populations(GENOME.class, list(parents$V1, FS$V1, S1$V1, S2$V1, S3$V1, S4$V1, S5$V1), diploid = T)
  # GENOME.class <- set.populations(GENOME.class, list(VI_BC$V1, COAST_BC$V1, HG_BC$V1, INT_BC$V1, CA_OR$V1), diploid = T)
  genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)
  # Split by gene
  # genes <- splitting.data(GENOME.class, subsites="gene")

  # Neutrality stats for non-synonymous sites
  genes <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_parents_sc",
                         scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[2]],
  #             file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_FS_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[3]],
  #             file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_S1_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[4]],
  #             file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_S2_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[5]],
  #             file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_S3_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[6]],
  #             file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_S4_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[7]],
  #             file=paste("Neutrality_stats_gene_nonsyn/neutrality_stats_nonsyn_S5_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # Neutrality stats for synonymous sites
  genes <- neutrality.stats(genes, subsites="syn", FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_parents_sc",
                         scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[2]],
  #             file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_FS_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[3]],
  #             file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_S1_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[4]],
  #             file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_S2_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[5]],
  #             file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_S3_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[6]],
  #             file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_S4_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[7]],
  #             file=paste("Neutrality_stats_gene_syn/neutrality_stats_syn_S5_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # Linkage stats for non-synonymous sites
  # genes <- linkage.stats(genes, subsites="nonsyn")
  # write.table(get.linkage(genes)[[1]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_parents_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[2]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_FS_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[3]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_S1_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[4]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_S2_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[5]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_S3_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[6]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_S4_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[7]],
  #             file=paste("Linkage_stats_gene_nonsyn/linkage_stats_nonsyn_S5_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # 
  # # Linkage stats for synonymous sites
  # genes <- linkage.stats(genes, subsites="syn")
  # write.table(get.linkage(genes)[[1]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_parents_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[2]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_FS_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[3]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_S1_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[4]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_S2_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[5]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_S3_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[6]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_S4_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[7]],
  #             file=paste("Linkage_stats_gene_syn/linkage_stats_syn_S5_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # FST stats non-synonymous
  # genes <- F_ST.stats(genes, mode = "nucleotide", subsites = "nonsyn")
  # # write.table(get.F_ST(genes), file=paste("FST_stats_gene/fst_stats_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_gene_nonsyn/pairwise_FST_nonsyn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nucleotide.F_ST, file=paste("FST_stats_gene_nonsyn/fixation_index_based_on_MAF_nonsyn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nuc.diversity.within, file=paste("FST_stats_gene_nonsyn/nucleotide_diversity_within_populations_nonsyn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_gene_nonsyn/nucleotide_diversity_one_vc_all_nonsyn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_gene_nonsyn/nucleotide_diversity_between_populations_nonsyn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # 
  # # FST stats synonymous
  # genes <- F_ST.stats(genes, mode = "nucleotide", subsites = "syn")
  # # write.table(get.F_ST(genes), file=paste("FST_stats_gene/fst_stats_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_gene_syn/pairwise_FST_syn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nucleotide.F_ST, file=paste("FST_stats_gene_syn/fixation_index_based_on_MAF_syn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nuc.diversity.within, file=paste("FST_stats_gene_syn/nucleotide_diversity_within_populations_syn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_gene_syn/nucleotide_diversity_one_vc_all_syn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_gene_syn/nucleotide_diversity_between_populations_syn_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # Diversity stats non-synonymous
  genes <- diversity.stats(genes, subsites="nonsyn", pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_nonsyn/pi_nucleotide_diversity_within_populations_nonsyn_sc", scaffold_names[i],
                                          "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # Diversity stats synonymous
  genes <- diversity.stats(genes, subsites="syn", pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_syn/pi_nucleotide_diversity_within_populations_syn_sc", scaffold_names[i],
                                          "_by_gene.txt", sep=""), quote=F, sep='\t')
}
  # genes <- MKT(genes, do.fisher.test = T)
  # for (x in 1:genes@genelength){
  #   write.table(get.MKT(genes)[[x]], file=paste("MKT_stats_gene_with_snps_diploid/mkt_stats_gene",x,"_sc",
  #                                               gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # }
  # 
  # gff_info <- get_gff_info(gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
  # chr= gff_scaffold_names_correct_order[i], feature = "gene")
  # if (class(gff_info) == "matrix"){
  # gff_info <- list(gff_info)
 # }

# Analysis for diverse population - by gene  
for (i in 1:length(scaffold_names)){
  GENOME.class <- readVCF(paste(path, "VCF_", scaffold_names[i], ".recode.vcf.gz", sep = ""), numcols = 100000, 
                          tid = scaffold_names[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE)
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, scaffold_names[i], sep = ""),
  #                                       chr = scaffold_names[i], feature = "gene")
  gff_info <- get_gff_info(gff.file = paste(path, scaffold_names[i], sep = ""),
                           chr= scaffold_names[i], feature = "gene")
  if (class(gff_info) == "matrix"){
    gff_info <- list(gff_info)
  }
 
  # Set populations
  # GENOME.class <- set.populations(GENOME.class, list(parents$V1, FS$V1, S1$V1, S2$V1, S3$V1, S4$V1, S5$V1), diploid = T)
  # GENOME.class <- set.populations(GENOME.class, list(VI_BC$V1, COAST_BC$V1, HG_BC$V1, INT_BC$V1, CA_OR$V1), diploid = T)
  genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)
  
  # Neutrality stats
  genes <- neutrality.stats(genes, FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_new/neutrality_stats_parents_sc",
                         scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[2]],
  #             file=paste("Neutrality_stats_genes/neutrality_stats_FS_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[3]],
  #             file=paste("Neutrality_stats_genes/neutrality_stats_S1_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[4]],
  #             file=paste("Neutrality_stats_genes/neutrality_stats_S2_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[5]],
  #             file=paste("Neutrality_stats_genes/neutrality_stats_S3_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[6]],
  #             file=paste("Neutrality_stats_genes/neutrality_stats_S4_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[7]],
  #             file=paste("Neutrality_stats_genes/neutrality_stats_S5_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # # Linkage stats
  genes <- linkage.stats(genes)
  write.table(get.linkage(genes)[[1]], file=paste("Linkage_stats_genes/linkage_stats_parents_sc",
                                                  scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[2]], file=paste("Linkage_stats_genes/linkage_stats_FS_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[3]], file=paste("Linkage_stats_genes/linkage_stats_S1_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[4]], file=paste("Linkage_stats_genes/linkage_stats_S2_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[5]], file=paste("Linkage_stats_genes/linkage_stats_S3_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[6]], file=paste("Linkage_stats_genes/linkage_stats_S4_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[7]], file=paste("Linkage_stats_genes/linkage_stats_S5_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # 
  # # FST stats
  # genes <- F_ST.stats(genes, mode = "nucleotide")
  #    # write.table(get.F_ST(genes), file=paste("FST_stats_genes/fst_stats_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  #    write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_genes/pairwise_FST_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  #    write.table(genes@nucleotide.F_ST, file=paste("FST_stats_genes/fixation_index_based_on_MAF_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  #    write.table(genes@nuc.diversity.within, file=paste("FST_stats_genes/nucleotide_diversity_within_populations_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  #    write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_genes/nucleotide_diversity_one_vc_all_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  #    write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_genes/nucleotide_diversity_between_populations_sc", scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # Diversity stats
  genes <- diversity.stats(genes, pi=TRUE)
     write.table(genes@Pi, file=paste("FST_stats_gene_new/pi_nucleotide_diversity_within_populations_sc", scaffold_names[i],
                                      "_by_gene.txt", sep=""), quote=F, sep='\t')
}
## Analysis of S lines - not used
# Load VCF and GFF files for each scaffold, set subpopulations, split data by gene and calculate neutrality and linkage statistics 
for (i in 1:length(gff_scaffold_names_correct_order)){
  GENOME.class <- readVCF(paste(path, "VCF_", gff_scaffold_names_correct_order[i], ".recode.vcf.gz", sep = ""), numcols = 100000, 
                          tid = gff_scaffold_names_correct_order[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE, gffpath = paste(path, gff_scaffold_names_correct_order[i], sep = ""))
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
  #                                       chr = gff_scaffold_names_correct_order[i], feature = "gene")
  gff_info <- get_gff_info(gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
                           chr= gff_scaffold_names_correct_order[i], feature = "gene")
  if (class(gff_info) == "matrix"){
    gff_info <- list(gff_info)
  }
  # Get synonymous/nonsynonymous sites
  GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste(path, gff_scaffold_names_correct_order[i], ".fa", sep = ""))
  # Set populations
  GENOME.class <- set.populations(GENOME.class, list(parents$V1, FS$V1, S1$V1, S2$V1, S3$V1, S4$V1, S5$V1), diploid = T)
  # GENOME.class <- set.populations(GENOME.class, list(parents$V1, FS$V1, S1$V1, S2$V1, S3$V1), diploid = T)
  # genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)
  # Split by gene
  genes <- splitting.data(GENOME.class, subsites="gene")
  
  # # Neutrality stats for non-synonymous sites
  genes <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_parents_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[2]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_FS_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[3]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_S1_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[4]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_S2_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[5]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_S3_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[6]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_S4_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[7]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_/neutrality_stats_nonsyn_S5_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # Neutrality stats for synonymous sites
  genes <- neutrality.stats(genes, subsites="syn", FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_parents_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[2]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_FS_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[3]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_S1_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[4]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_S2_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[5]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_S3_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[6]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_S4_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[7]],
              file=paste("Neutrality_stats_gene_with_snps_syn_/neutrality_stats_syn_S5_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # # Linkage stats
  genes <- linkage.stats(genes, subsites = "nonsyn")
  write.table(get.linkage(genes)[[1]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__parents_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[2]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__FS_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[3]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__S1_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[4]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__S2_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[5]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__S3_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[6]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__S4_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[7]], file=paste("Linkage_stats_gene_with_snps_nonsyn_/Linkage_stats_gene_with_snps_nonsyn__S5_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # # Linkage stats
  genes <- linkage.stats(genes, subsites = "syn")
  write.table(get.linkage(genes)[[1]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__parents_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[2]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__FS_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[3]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__S1_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[4]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__S2_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[5]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__S3_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[6]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__S4_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[7]], file=paste("Linkage_stats_gene_with_snps_syn_/Linkage_stats_gene_with_snps_syn__S5_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # FST stats non-synonymous
  genes <- F_ST.stats(genes, subsites = "nonsyn", mode = "nucleotide")
  # write.table(get.F_ST(genes), file=paste("FST_stats_gene_with_snps_nonsyn_/FST_stats_gene_with_snps_nonsyn__sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_gene_with_snps_nonsyn_/pairwise_FST_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nucleotide.F_ST, file=paste("FST_stats_gene_with_snps_nonsyn_/fixation_index_based_on_MAF_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.diversity.within, file=paste("FST_stats_gene_with_snps_nonsyn_/nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_gene_with_snps_nonsyn_/nucleotide_diversity_one_vc_all_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_gene_with_snps_nonsyn_/nucleotide_diversity_between_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # Diversity stats non-synonymous
  genes <- diversity.stats(genes, subsites="nonsyn", pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_nonsyn_/pi_nucleotide_diversity_within_populations_nonsyn_sc", gff_scaffold_names_correct_order[i],
                                   "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # FST stats synonymous
  genes <- F_ST.stats(genes, subsites = "syn", mode = "nucleotide")
  # write.table(get.F_ST(genes), file=paste("FST_stats_gene_with_snps_syn_/FST_stats_gene_with_snps_syn__sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_gene_with_snps_syn_/pairwise_FST_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nucleotide.F_ST, file=paste("FST_stats_gene_with_snps_syn_/fixation_index_based_on_MAF_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.diversity.within, file=paste("FST_stats_gene_with_snps_syn_/nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_gene_with_snps_syn_/nucleotide_diversity_one_vc_all_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_gene_with_snps_syn_/nucleotide_diversity_between_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # Diversity stats synonymous
  genes <- diversity.stats(genes, subsites="syn", pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_syn_/pi_nucleotide_diversity_within_populations_syn_sc", gff_scaffold_names_correct_order[i],
                                   "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # genes <- MKT(genes, do.fisher.test = T)
  # for (x in 1:genes@genelength){
  #   write.table(get.MKT(genes)[[x]], file=paste("MKT_stats_gene_with_snps_diploid/mkt_stats_gene",x,"_sc",
  #                                               gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
}

## Analysis of S lines - not used
# Load VCF and GFF files for each scaffold, set subpopulations, split data by gene and calculate neutrality and linkage statistics 
for (i in 41:length(gff_scaffold_names_correct_order)){
  GENOME.class <- readVCF(paste(path, "VCF_", gff_scaffold_names_correct_order[i], ".recode.vcf.gz", sep = ""), numcols = 100000, 
                          tid = gff_scaffold_names_correct_order[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE, gffpath = paste(path, gff_scaffold_names_correct_order[i], sep = ""))
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
  #                                       chr = gff_scaffold_names_correct_order[i], feature = "gene")
  # Get synonymous/nonsynonymous sites
  GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste(path, gff_scaffold_names_correct_order[i], ".fa", sep = ""))
  # Set populations
  GENOME.class <- set.populations(GENOME.class, list(Line_1$V1, Line_6$V1, Line_7$V1, Line_8$V1, Line_10$V1, Line_13$V1, Line_16$V1,
                                                     Line_17$V1, Line_19$V1, Line_20$V1, Line_21$V1, Line_23$V1, Line_26$V1, Line_27$V1, 
                                                     Line_29$V1), diploid = T)
  # genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)
  # Split by gene
  genes <- splitting.data(GENOME.class, subsites="gene")
  
  # Neutrality stats for non-synonymous sites
  genes <- neutrality.stats(genes, subsites = "nonsyn", FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_1_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[2]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_6_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[3]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_7_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[4]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_8_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[5]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_10_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[6]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_13_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[7]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_16_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[8]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_17_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[9]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_19_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[10]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_20_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[11]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_21_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[12]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_23_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[13]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_26_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[14]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_27_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[15]],
              file=paste("Neutrality_stats_gene_with_snps_nonsyn_lines/neutrality_stats_nonsyn_Line_29_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  
  # Neutrality stats for synonymous sites
  genes <- neutrality.stats(genes, subsites = "syn", FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_1_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[2]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_6_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[3]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_7_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[4]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_8_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[5]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_10_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[6]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_13_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[7]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_16_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[8]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_17_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[9]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_19_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[10]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_20_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[11]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_21_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[12]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_23_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[13]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_26_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[14]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_27_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[15]],
              file=paste("Neutrality_stats_gene_with_snps_syn_lines/neutrality_stats_syn_Line_29_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  
  # Diversity stats non-synonymous
  genes <- diversity.stats(genes, subsites="nonsyn", pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_nonsyn_lines/pi_nucleotide_diversity_within_populations_nonsyn_sc", gff_scaffold_names_correct_order[i],
                                   "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # Diversity stats synonymous
  genes <- diversity.stats(genes, subsites="syn", pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_syn_lines/pi_nucleotide_diversity_within_populations_syn_sc", gff_scaffold_names_correct_order[i],
                                   "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  genes <- MKT(genes, do.fisher.test = T)
  for (x in 1:genes@genelength){
    write.table(get.MKT(genes)[[x]], file=paste("MKT_stats_gene_with_snps_lines/mkt_stats_gene",x,"_sc",
                                                gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  }
  
  # genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)
  # # Neutrality stats
  # genes <- neutrality.stats(genes, FAST=TRUE)
  # write.table(get.neutrality(genes, theta=TRUE)[[1]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_1_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[2]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_6_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[3]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_7_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[4]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_8_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[5]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_10_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[6]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_13_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[7]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_16_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[8]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_17_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[9]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_19_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[10]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_20_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[11]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_21_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[12]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_23_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[13]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_26_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[14]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_27_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[15]],
  #             file=paste("Neutrality_stats_gene_with_snps_lines/neutrality_stats_Line_29_sc", 
  #                        gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # 
  # # Linkage stats
  # genes <- linkage.stats(genes)
  # write.table(get.linkage(genes)[[1]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_1_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[2]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_6_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[3]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_7_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[4]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_8_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[5]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_10_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[6]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_13_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[7]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_16_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[8]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_17_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[9]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_19_sc", 
  #                                                 gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[10]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_20_sc", 
  #                                                  gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[11]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_21_sc", 
  #                                                  gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[12]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_23_sc", 
  #                                                  gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[13]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_26_sc", 
  #                                                  gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[14]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_27_sc", 
  #                                                  gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(get.linkage(genes)[[15]], file=paste("Linkage_stats_gene_with_snps_lines/linkage_stats_Line_29_sc", 
  #                                                  gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # 
  # # FST stats
  # genes <- F_ST.stats(genes)
  # write.table(get.F_ST(genes), file=paste("FST_stats_gene_with_snps_lines/fst_stats_sc", gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_gene_with_snps_lines/pairwise_FST_sc", gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nucleotide.F_ST, file=paste("FST_stats_gene_with_snps_lines/fixation_index_based_on_MAF_sc", gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nuc.diversity.within, file=paste("FST_stats_gene_with_snps_lines/nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_gene_with_snps_lines/nucleotide_diversity_one_vc_all_sc", gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_gene_with_snps_lines/nucleotide_diversity_between_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
  # # Diversity stats
  # genes <- diversity.stats(genes, pi=TRUE)
  # write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_lines/pi_nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i],
  #                                  "_by_gene_with_snps_lines.txt", sep=""), quote=F, sep='\t')
}

## Analysis of S lines (not used)
# Load VCF and GFF files for each scaffold, set subpopulations, split data by gene and calculate neutrality and linkage statistics 
for (i in 2498:length(gff_scaffold_names_correct_order)){
  GENOME.class <- readVCF(paste(path, "VCF_", gff_scaffold_names_correct_order[i], ".recode.vcf.gz", sep = ""), numcols = 100000, 
                          tid = gff_scaffold_names_correct_order[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE)
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
  #                                       chr = gff_scaffold_names_correct_order[i], feature = "gene")
  gff_info <- get_gff_info(gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
                           chr= gff_scaffold_names_correct_order[i], feature = "gene")
  if (class(gff_info) == "matrix"){
    gff_info <- list(gff_info)
  }
  # Get synonymous/nonsynonymous sites
  # GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste(path, gff_scaffold_names_correct_order[i], ".fa", sep = ""))
  # Set populations

  GENOME.class <- set.populations(GENOME.class, list(Line_6_FS$V1, Line_6_S1$V1, Line_6_S4$V1, 
                                                     Line_7_FS$V1, Line_7_S1$V1, Line_7_S4$V1, 
                                                     Line_10_FS$V1, Line_10_S1$V1, Line_10_S4$V1, Line_13_FS$V1, Line_13_S1$V1, 
                                                     Line_13_S4$V1, Line_16_FS$V1, Line_16_S1$V1, Line_16_S4$V1, Line_17_FS$V1, 
                                                     Line_17_S1$V1, Line_17_S4$V1, Line_19_FS$V1, Line_19_S1$V1, Line_19_S4$V1), diploid = T)
  genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)

  # Neutrality stats
  genes <- neutrality.stats(genes, FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_6_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[2]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_6_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[3]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_6_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[4]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_7_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[5]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_7_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[6]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_7_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[7]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_10_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[8]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_10_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[9]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_10_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[10]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_13_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[11]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_13_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[12]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_13_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[13]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_16_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[14]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_16_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[15]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_16_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[16]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_17_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[17]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_17_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[18]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_17_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[19]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_19_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[20]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_19_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[21]],
              file=paste("Neutrality_stats_gene_generation_line_1456_mislabeled_removed/neutrality_stats_Line_19_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  
  # Linkage stats
  genes <- linkage.stats(genes)
  write.table(get.linkage(genes)[[1]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_6_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[2]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_6_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[3]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_6_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[4]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_7_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[5]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_7_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[6]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_7_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[7]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_10_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[8]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_10_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[9]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_10_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[10]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_13_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[11]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_13_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[12]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_13_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[13]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_16_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[14]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_16_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[15]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_16_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[16]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_17_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[17]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_17_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[18]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_17_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[19]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_19_FS_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[20]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_19_S1_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[21]],
              file=paste("Linkage_stats_gene_generation_line_1456_mislabeled_removed/linkage_stats_Line_19_S4_sc", 
                         gff_scaffold_names_correct_order[i], "_by_gene_lines.txt", sep=""), quote=F, sep='\t')
  
  # FST stats
  genes <- F_ST.stats(genes)
  write.table(get.F_ST(genes), file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/fst_stats_sc", gff_scaffold_names_correct_order[i], "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/pairwise_FST_sc", gff_scaffold_names_correct_order[i], "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nucleotide.F_ST, file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/fixation_index_based_on_MAF_sc", gff_scaffold_names_correct_order[i], "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.diversity.within, file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/nucleotide_diversity_one_vs_all_sc", gff_scaffold_names_correct_order[i], "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/nucleotide_diversity_between_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
  # Diversity stats
  genes <- diversity.stats(genes, pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_gene_generation_line_1456_mislabeled_removed/pi_nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i],
                                   "_by_gene_generation_line.txt", sep=""), quote=F, sep='\t')
}

## Parents split by location
# Load VCF and GFF files for each scaffold, set subpopulations, split data by gene and calculate neutrality and linkage statistics 
for (i in 1:length(gff_scaffold_names_correct_order)){
  GENOME.class <- readVCF(paste(path, "VCF_", gff_scaffold_names_correct_order[i], ".recode.vcf.gz", sep = ""), numcols = 100000, 
                          tid = gff_scaffold_names_correct_order[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE)
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
  #                                       chr = gff_scaffold_names_correct_order[i], feature = "gene")
  gff_info <- get_gff_info(gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
                           chr= gff_scaffold_names_correct_order[i], feature = "gene")
  if (class(gff_info) == "matrix"){
    gff_info <- list(gff_info)
  }
  # Get synonymous/nonsynonymous sites
  # GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste(path, gff_scaffold_names_correct_order[i], ".fa", sep = ""))
  # Set populations
  # GENOME.class <- set.populations(GENOME.class, list(parents$V1, FS$V1, S1$V1, S2$V1, S3$V1, S4$V1, S5$V1), diploid = T)
  GENOME.class <- set.populations(GENOME.class, list(VI_BC$V1, COAST_BC$V1, HG_BC$V1, INT_BC$V1, CA_OR$V1), diploid = T)
  genes <- splitting.data(GENOME.class, positions = gff_info, type = 2)
  # Split by gene
  # genes <- splitting.data(GENOME.class, subsites="gene")
  
  # # Neutrality stats for non-synonymous sites
  # genes <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)
  # write.table(get.neutrality(genes, theta=TRUE)[[1]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_parents_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[2]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_FS_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[3]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_S1_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[4]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_S2_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[5]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_S3_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[6]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_S4_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[7]],
  #             file=paste("Neutrality_stats_gene_with_snps_nonsyn/neutrality_stats_nonsyn_S5_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # 
  # # Neutrality stats for synonymous sites
  # genes <- neutrality.stats(genes, subsites="syn", FAST=TRUE)
  # write.table(get.neutrality(genes, theta=TRUE)[[1]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_parents_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[2]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_FS_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[3]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_S1_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[4]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_S2_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[5]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_S3_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[6]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_S4_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # write.table(get.neutrality(genes, theta=TRUE)[[7]],
  #             file=paste("Neutrality_stats_gene_with_snps_syn/neutrality_stats_syn_S5_sc",
  #                        gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # 
  # # Diversity stats non-synonymous
  # genes <- diversity.stats(genes, subsites="nonsyn", pi=TRUE)
  # write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_nonsyn/pi_nucleotide_diversity_within_populations_nonsyn_sc", gff_scaffold_names_correct_order[i],
  #                                  "_by_gene.txt", sep=""), quote=F, sep='\t')
  # 
  # # Diversity stats synonymous
  # genes <- diversity.stats(genes, subsites="syn", pi=TRUE)
  # write.table(genes@Pi, file=paste("FST_stats_gene_with_snps_syn/pi_nucleotide_diversity_within_populations_syn_sc", gff_scaffold_names_correct_order[i],
  #                                  "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # genes <- MKT(genes, do.fisher.test = T)
  # for (x in 1:genes@genelength){
  #   write.table(get.MKT(genes)[[x]], file=paste("MKT_stats_gene_with_snps_diploid/mkt_stats_gene",x,"_sc",
  #                                               gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # }
  # 
  # gff_info <- get_gff_info(gff.file = paste(path, gff_scaffold_names_correct_order[i], sep = ""),
  # chr= gff_scaffold_names_correct_order[i], feature = "gene")
  # if (class(gff_info) == "matrix"){
  # gff_info <- list(gff_info)
  # }
  
  # Neutrality stats
  genes <- neutrality.stats(genes, FAST=TRUE)
  write.table(get.neutrality(genes, theta=TRUE)[[1]],
              file=paste("Neutrality_stats_parents_by_location/neutrality_stats_VI_BC_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[2]],
              file=paste("Neutrality_stats_parents_by_location/neutrality_stats_COAST_BC_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[3]],
              file=paste("Neutrality_stats_parents_by_location/neutrality_stats_HG_BC_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[4]],
              file=paste("Neutrality_stats_parents_by_location/neutrality_stats_INT_BC_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.neutrality(genes, theta=TRUE)[[5]],
              file=paste("Neutrality_stats_parents_by_location/neutrality_stats_CA_OR_sc",
                         gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # # Linkage stats
  genes <- linkage.stats(genes)
  write.table(get.linkage(genes)[[1]], file=paste("Linkage_stats_parents_by_location/linkage_stats_VI_BC_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[2]], file=paste("Linkage_stats_parents_by_location/linkage_stats_COAST_BC_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[3]], file=paste("Linkage_stats_parents_by_location/linkage_stats_HG_BC_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[4]], file=paste("Linkage_stats_parents_by_location/linkage_stats_INT_BC_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(get.linkage(genes)[[5]], file=paste("Linkage_stats_parents_by_location/linkage_stats_CA_OR_sc",
                                                  gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  
  # FST stats
  genes <- F_ST.stats(genes, mode = "nucleotide")
  # write.table(get.F_ST(genes), file=paste("FST_stats_parents_by_location/fst_stats_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.F_ST.pairwise), file=paste("FST_stats_parents_by_location/pairwise_FST_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nucleotide.F_ST, file=paste("FST_stats_parents_by_location/fixation_index_based_on_MAF_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.diversity.within, file=paste("FST_stats_parents_by_location/nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(genes@nuc.F_ST.vs.all, file=paste("FST_stats_parents_by_location/nucleotide_diversity_one_vc_all_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  write.table(t(genes@nuc.diversity.between), file=paste("FST_stats_parents_by_location/nucleotide_diversity_between_populations_sc", gff_scaffold_names_correct_order[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
  # Diversity stats
  genes <- diversity.stats(genes, pi=TRUE)
  write.table(genes@Pi, file=paste("FST_stats_parents_by_location/pi_nucleotide_diversity_within_populations_sc", gff_scaffold_names_correct_order[i],
                                   "_by_gene.txt", sep=""), quote=F, sep='\t')
}


#### Normalising pop gen stats, statistics for testing differences
setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/PopGenome_Results_final_annotation/")

# Read in diversity stats
# all_pairwise_FST_parent_locations_genes <- read.table("FST_stats_parents_by_location/pairwise_fst_parents_fixed.txt", header = T)
# all_pairwise_FST_parent_locations_scaffolds <- read.table("FST_stats_parents_by_location_scaffolds/all_pairwise_fst_corrected.txt", header = T)
# 
# all_pairwise_FST_parent_locations_syn <- read.table("FST_stats_gene_with_snps_syn_parents_by_location/pairwise_fst_parents_syn_fixed.txt", header = T)
# all_pairwise_FST_parent_locations_nonsyn <- read.table("FST_stats_gene_with_snps_nonsyn_parents_by_location/pairwise_fst_parents_nonsyn_fixed.txt", header = T)

# all_nucl_between_parent_locations_genes <- read.table("FST_stats_parents_by_location/all_nucleotide_diversity_between_stats_fixed.txt", header = T) 
# all_nucl_between_parent_locations_scaffolds <- read.table("FST_stats_parents_by_location_scaffolds/all_nucleotide_diversity_between_stats_fixed.txt", header = T) 

all_pi <- read.table("FST_stats_genes/all_pi_stats_fixed.txt")
all_pi_parents_by_location <- read.table("FST_stats_parents_by_location/all_pi_stats_fixed.txt")

all_pi_syn <- read.table("FST_stats_gene_syn/all_pi_stats_syn_fixed.txt")
all_pi_nonsyn <- read.table("FST_stats_gene_nonsyn/all_pi_stats_nonsyn_fixed.txt")

all_pi_syn_new <- read.table("FST_stats_gene_syn_new/all_pi_stats_fixed.txt")
all_pi_nonsyn_new <- read.table("FST_stats_gene_nonsyn_new/all_pi_stats_fixed.txt")

all_pi_parents_by_location_syn <- read.table("FST_stats_gene_with_snps_syn_parents_by_location/all_pi_stats_syn_fixed.txt")
all_pi_parents_by_location_nonsyn <- read.table("FST_stats_gene_with_snps_nonsyn_parents_by_location/all_pi_stats_nonsyn_fixed.txt")

# Read in neutrality stats
all_neutrality_gene_parents_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_parents_fixed.txt", sep = "\t")
all_neutrality_gene_parents_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_parents_fixed.txt", sep = "\t")
all_neutrality_gene_FS_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_FS_fixed.txt", sep = "\t")
all_neutrality_gene_FS_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_FS_fixed.txt", sep = "\t")
all_neutrality_gene_S1_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_S1_fixed.txt", sep = "\t")
all_neutrality_gene_S1_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_S1_fixed.txt", sep = "\t")
all_neutrality_gene_S2_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_S2_fixed.txt", sep = "\t")
all_neutrality_gene_S2_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_S2_fixed.txt", sep = "\t")
all_neutrality_gene_S3_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_S3_fixed.txt", sep = "\t")
all_neutrality_gene_S3_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_S3_fixed.txt", sep = "\t")
all_neutrality_gene_S4_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_S4_fixed.txt", sep = "\t")
all_neutrality_gene_S4_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_S4_fixed.txt", sep = "\t")
all_neutrality_gene_S5_syn <- read.table("Neutrality_stats_gene_syn/all_neutrality_stats_syn_S5_fixed.txt", sep = "\t")
all_neutrality_gene_S5_nonsyn <- read.table("Neutrality_stats_gene_nonsyn/all_neutrality_stats_nonsyn_S5_fixed.txt", sep = "\t")

all_neutrality_gene_parents <- read.table("Neutrality_stats_genes/all_neutrality_stats_parents_fixed.txt", sep = "\t")
all_neutrality_gene_FS <- read.table("Neutrality_stats_genes/all_neutrality_stats_FS_fixed.txt", sep = "\t")
all_neutrality_gene_S1 <- read.table("Neutrality_stats_genes/all_neutrality_stats_S1_fixed.txt", sep = "\t")
all_neutrality_gene_S2 <- read.table("Neutrality_stats_genes/all_neutrality_stats_S2_fixed.txt", sep = "\t")
all_neutrality_gene_S3 <- read.table("Neutrality_stats_genes/all_neutrality_stats_S3_fixed.txt", sep = "\t")
all_neutrality_gene_S4 <- read.table("Neutrality_stats_genes/all_neutrality_stats_S4_fixed.txt", sep = "\t")
all_neutrality_gene_S5 <- read.table("Neutrality_stats_genes/all_neutrality_stats_S5_fixed.txt", sep = "\t")

all_neutrality_gene_VI_BC <- read.table("Neutrality_stats_parents_by_location/all_neutrality_stats_VI_BC_fixed.txt", sep = "\t")
all_neutrality_gene_COAST_BC <- read.table("Neutrality_stats_parents_by_location/all_neutrality_stats_COAST_BC_fixed.txt", sep = "\t")
all_neutrality_gene_HG_BC <- read.table("Neutrality_stats_parents_by_location/all_neutrality_stats_HG_BC_fixed.txt", sep = "\t")
all_neutrality_gene_INT_BC <- read.table("Neutrality_stats_parents_by_location/all_neutrality_stats_INT_BC_fixed.txt", sep = "\t")
all_neutrality_gene_CA_OR <- read.table("Neutrality_stats_parents_by_location/all_neutrality_stats_CA_OR_fixed.txt", sep = "\t")

all_neutrality_gene_syn_VI_BC <- read.table("Neutrality_stats_gene_with_snps_syn_parents_by_location/all_neutrality_stats_syn_VI_BC_fixed.txt", sep = "\t")
all_neutrality_gene_syn_COAST_BC <- read.table("Neutrality_stats_gene_with_snps_syn_parents_by_location/all_neutrality_stats_syn_COAST_BC_fixed.txt", sep = "\t")
all_neutrality_gene_syn_HG_BC <- read.table("Neutrality_stats_gene_with_snps_syn_parents_by_location/all_neutrality_stats_syn_HG_BC_fixed.txt", sep = "\t")
all_neutrality_gene_syn_INT_BC <- read.table("Neutrality_stats_gene_with_snps_syn_parents_by_location/all_neutrality_stats_syn_INT_BC_fixed.txt", sep = "\t")
all_neutrality_gene_syn_CA_OR <- read.table("Neutrality_stats_gene_with_snps_syn_parents_by_location/all_neutrality_stats_syn_CA_OR_fixed.txt", sep = "\t")

all_neutrality_gene_nonsyn_VI_BC <- read.table("Neutrality_stats_gene_with_snps_nonsyn_parents_by_location/all_neutrality_stats_nonsyn_VI_BC_fixed.txt", sep = "\t")
all_neutrality_gene_nonsyn_COAST_BC <- read.table("Neutrality_stats_gene_with_snps_nonsyn_parents_by_location/all_neutrality_stats_nonsyn_COAST_BC_fixed.txt", sep = "\t")
all_neutrality_gene_nonsyn_HG_BC <- read.table("Neutrality_stats_gene_with_snps_nonsyn_parents_by_location/all_neutrality_stats_nonsyn_HG_BC_fixed.txt", sep = "\t")
all_neutrality_gene_nonsyn_INT_BC <- read.table("Neutrality_stats_gene_with_snps_nonsyn_parents_by_location/all_neutrality_stats_nonsyn_INT_BC_fixed.txt", sep = "\t")
all_neutrality_gene_nonsyn_CA_OR <- read.table("Neutrality_stats_gene_with_snps_nonsyn_parents_by_location/all_neutrality_stats_nonsyn_CA_OR_fixed.txt", sep = "\t")


# Read in linkage stats
all_linkage_gene_parents <- read.table("Linkage_stats_genes/all_linkage_stats_parents_fixed.txt", sep = "\t")
all_linkage_gene_FS <- read.table("Linkage_stats_genes/all_linkage_stats_FS_fixed.txt", sep = "\t")
all_linkage_gene_S1 <- read.table("Linkage_stats_genes/all_linkage_stats_S1_fixed.txt", sep = "\t")
all_linkage_gene_S2 <- read.table("Linkage_stats_genes/all_linkage_stats_S2_fixed.txt", sep = "\t")
all_linkage_gene_S3 <- read.table("Linkage_stats_genes/all_linkage_stats_S3_fixed.txt", sep = "\t")
all_linkage_gene_S4 <- read.table("Linkage_stats_genes/all_linkage_stats_S4_fixed.txt", sep = "\t")
all_linkage_gene_S5 <- read.table("Linkage_stats_genes/all_linkage_stats_S5_fixed.txt", sep = "\t")

all_linkage_gene_parents_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_parents_fixed.txt", sep = "\t")
all_linkage_gene_parents_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_parents_fixed.txt", sep = "\t")
all_linkage_gene_FS_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_FS_fixed.txt", sep = "\t")
all_linkage_gene_FS_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_FS_fixed.txt", sep = "\t")
all_linkage_gene_S1_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_S1_fixed.txt", sep = "\t")
all_linkage_gene_S1_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_S1_fixed.txt", sep = "\t")
all_linkage_gene_S2_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_S2_fixed.txt", sep = "\t")
all_linkage_gene_S2_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_S2_fixed.txt", sep = "\t")
all_linkage_gene_S3_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_S3_fixed.txt", sep = "\t")
all_linkage_gene_S3_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_S3_fixed.txt", sep = "\t")
all_linkage_gene_S4_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_S4_fixed.txt", sep = "\t")
all_linkage_gene_S4_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_S4_fixed.txt", sep = "\t")
all_linkage_gene_S5_nonsyn <- read.table("Linkage_stats_gene_nonsyn/all_linkage_stats_nonsyn_S5_fixed.txt", sep = "\t")
all_linkage_gene_S5_syn <- read.table("Linkage_stats_gene_syn/all_linkage_stats_syn_S5_fixed.txt", sep = "\t")

all_linkage_gene_VI_BC <- read.table("Linkage_stats_parents_by_location/all_linkage_stats_VI_BC_fixed.txt", sep = "\t")
all_linkage_gene_COAST_BC <- read.table("Linkage_stats_parents_by_location/all_linkage_stats_COAST_BC_fixed.txt", sep = "\t")
all_linkage_gene_HG_BC <- read.table("Linkage_stats_parents_by_location/all_linkage_stats_HG_BC_fixed.txt", sep = "\t")
all_linkage_gene_INT_BC <- read.table("Linkage_stats_parents_by_location/all_linkage_stats_INT_BC_fixed.txt", sep = "\t")
all_linkage_gene_CA_OR <- read.table("Linkage_stats_parents_by_location/all_linkage_stats_CA_OR_fixed.txt", sep = "\t")

all_linkage_gene_syn_VI_BC <- read.table("Linkage_stats_gene_with_snps_syn_parents_by_location/all_linkage_stats_syn_VI_BC_fixed.txt", sep = "\t")
all_linkage_gene_syn_COAST_BC <- read.table("Linkage_stats_gene_with_snps_syn_parents_by_location/all_linkage_stats_syn_COAST_BC_fixed.txt", sep = "\t")
all_linkage_gene_syn_HG_BC <- read.table("Linkage_stats_gene_with_snps_syn_parents_by_location/all_linkage_stats_syn_HG_BC_fixed.txt", sep = "\t")
all_linkage_gene_syn_INT_BC <- read.table("Linkage_stats_gene_with_snps_syn_parents_by_location/all_linkage_stats_syn_INT_BC_fixed.txt", sep = "\t")
all_linkage_gene_syn_CA_OR <- read.table("Linkage_stats_gene_with_snps_syn_parents_by_location/all_linkage_stats_syn_CA_OR_fixed.txt", sep = "\t")

all_linkage_gene_nonsyn_VI_BC <- read.table("Linkage_stats_gene_with_snps_nonsyn_parents_by_location/all_linkage_stats_nonsyn_VI_BC_fixed.txt", sep = "\t")
all_linkage_gene_nonsyn_COAST_BC <- read.table("Linkage_stats_gene_with_snps_nonsyn_parents_by_location/all_linkage_stats_nonsyn_COAST_BC_fixed.txt", sep = "\t")
all_linkage_gene_nonsyn_HG_BC <- read.table("Linkage_stats_gene_with_snps_nonsyn_parents_by_location/all_linkage_stats_nonsyn_HG_BC_fixed.txt", sep = "\t")
all_linkage_gene_nonsyn_INT_BC <- read.table("Linkage_stats_gene_with_snps_nonsyn_parents_by_location/all_linkage_stats_nonsyn_INT_BC_fixed.txt", sep = "\t")
all_linkage_gene_nonsyn_CA_OR <- read.table("Linkage_stats_gene_with_snps_nonsyn_parents_by_location/all_linkage_stats_nonsyn_CA_OR_fixed.txt", sep = "\t")

# read in GFF as table
gff <- read.table("../Tplicatav3.1c.gene_exons_dashes_removed_final.gff3")

# Read in scaffold names
scaffold_names <- scan("../synonymous_scaffolds.txt", what = "factor")
# scaffold_names_genes <- scan("../scaffolds_with_genes.txt", what = "factor")

scaffold_names <- factor(scaffold_names)
# scaffold_names <- factor(scaffold_names_genes)

# Scaffold lengths
v3_scaffold_lengths <- read.table("../genome_scaffold_lengths_dashes_removed.txt")

v3_scaffold_lengths_in_vcf <- v3_scaffold_lengths %>% 
  filter(V1 %in% scaffold_names) %>% 
  arrange(V1)

# Filtering GFF file to match VCF
gff_genes_correct_order <- gff %>% 
  filter(V1 %in% scaffold_names) %>% 
  filter(V3 == "gene") %>%
  arrange(V1)



gff_genes_per_scaffold <- gff_genes_correct_order %>%
  group_by(V1) %>% 
  summarise(n_genes_on_scaffold = n())


# 
#### Pi estimates ###################################
### Pi Parents, S generations ########
pi_with_gene_names <- all_pi %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  # filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>%
  mutate(norm_pi_parents = V2 / gene_length, norm_pi_FS = V3 / gene_length, norm_pi_S1 = V4 / gene_length,
         norm_pi_S2 = V5 / gene_length, norm_pi_S3 = V6 / gene_length, norm_pi_S4 = V7 / gene_length, norm_pi_S5 = V8 / gene_length)

mean(pi_with_gene_names$norm_pi_parents)
mean(pi_with_gene_names$norm_pi_FS)
mean(pi_with_gene_names$norm_pi_S1)
mean(pi_with_gene_names$norm_pi_S2)
mean(pi_with_gene_names$norm_pi_S3)
mean(pi_with_gene_names$norm_pi_S4)
mean(pi_with_gene_names$norm_pi_S5)

wilcox.test(pi_with_gene_names$norm_pi_FS, pi_with_gene_names$norm_pi_S1)

x <- data.frame(pi_with_gene_names[,13:17])
kruskal.test(x)

#nonsyn
pi_nonsyn <- all_pi_nonsyn_new %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  # filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>%
  mutate(norm_pi_parents = V2 / gene_length)#, norm_pi_FS = V3 / gene_length, norm_pi_S1 = V4 / gene_length,
# norm_pi_S2 = V5 / gene_length, norm_pi_S3 = V6 / gene_length, norm_pi_S4 = V7 / gene_length, norm_pi_S5 = V8 / gene_length)

mean(pi_nonsyn$norm_pi_parents)
mean(pi_nonsyn$norm_pi_FS)
mean(pi_nonsyn$norm_pi_S1)
mean(pi_nonsyn$norm_pi_S2)
mean(pi_nonsyn$norm_pi_S3)
mean(pi_nonsyn$norm_pi_S4)
mean(pi_nonsyn$norm_pi_S5)

# syn
pi_syn <- all_pi_syn_new %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  # filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>%
  mutate(norm_pi_parents = V2 / gene_length)#, norm_pi_FS = V3 / gene_length, norm_pi_S1 = V4 / gene_length,
# norm_pi_S2 = V5 / gene_length, norm_pi_S3 = V6 / gene_length, norm_pi_S4 = V7 / gene_length, norm_pi_S5 = V8 / gene_length)

mean(pi_syn$norm_pi_parents)
pi_nonsyn_pi_syn_ratio <- data.frame(pi_nonsyn = pi_nonsyn$norm_pi_parents, pi_syn = pi_syn$norm_pi_parents, ratio = pi_nonsyn$norm_pi_parents / pi_syn$norm_pi_parents)
mean(pi_syn$norm_pi_FS)
mean(pi_syn$norm_pi_S1)
mean(pi_syn$norm_pi_S2)
mean(pi_syn$norm_pi_S3)
mean(pi_syn$norm_pi_S4)
mean(pi_syn$norm_pi_S5)

wilcox.test(pi_nonsyn$norm_pi_parents, pi_syn$norm_pi_parents, paired = T, alternative = "greater")
wilcox.test(pi_nonsyn$norm_pi_FS, pi_syn$norm_pi_FS, paired = T, alternative = "greater")
wilcox.test(pi_nonsyn$norm_pi_S1, pi_syn$norm_pi_S1, paired = T, alternative = "greater")
wilcox.test(pi_nonsyn$norm_pi_S2, pi_syn$norm_pi_S2, paired = T, alternative = "greater")
wilcox.test(pi_nonsyn$norm_pi_S3, pi_syn$norm_pi_S3, paired = T, alternative = "greater")
wilcox.test(pi_nonsyn$norm_pi_S4, pi_syn$norm_pi_S4, paired = T, alternative = "greater")
wilcox.test(pi_nonsyn$norm_pi_S5, pi_syn$norm_pi_S5, paired = T, alternative = "greater")

pi_nonsyn_syn_parents <- data.frame(pi_syn = pi_syn$norm_pi_parents, pi_nonsyn = pi_nonsyn$norm_pi_parents)

pi_nonsyn_syn_parents <- pi_nonsyn_syn_parents %>% 
  mutate(piA_piS = pi_nonsyn/pi_syn)

ggplot(pi_nonsyn_syn_parents, aes(x = piA_piS)) +
  geom_boxplot()

### Pi Parents by location ########
pi_parents_by_location <- all_pi_parents_by_location %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  # filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>%
  mutate(norm_pi_VI_BC = V2 / gene_length, norm_pi_COAST_BC = V3 / gene_length, norm_pi_HG_BC = V4 / gene_length,
         norm_pi_INT_BC = V5 / gene_length, norm_pi_CA_OR = V6 / gene_length)

mean(pi_parents_by_location$norm_pi_VI_BC)
mean(pi_parents_by_location$norm_pi_COAST_BC)
mean(pi_parents_by_location$norm_pi_HG_BC)
mean(pi_parents_by_location$norm_pi_INT_BC)
mean(pi_parents_by_location$norm_pi_CA_OR)

wilcox.test(pi_parents_by_location$norm_pi_VI_BC, pi_parents_by_location$norm_pi_COAST_BC)
wilcox.test(pi_parents_by_location$norm_pi_VI_BC, pi_parents_by_location$norm_pi_HG_BC)
wilcox.test(pi_parents_by_location$norm_pi_VI_BC, pi_parents_by_location$norm_pi_INT_BC)
wilcox.test(pi_parents_by_location$norm_pi_CA_OR, pi_parents_by_location$norm_pi_VI_BC)

wilcox.test(pi_parents_by_location$norm_pi_COAST_BC, pi_parents_by_location$norm_pi_HG_BC)
wilcox.test(pi_parents_by_location$norm_pi_COAST_BC, pi_parents_by_location$norm_pi_INT_BC)
wilcox.test(pi_parents_by_location$norm_pi_COAST_BC, pi_parents_by_location$norm_pi_CA_OR)

wilcox.test(pi_parents_by_location$norm_pi_HG_BC, pi_parents_by_location$norm_pi_INT_BC)
wilcox.test(pi_parents_by_location$norm_pi_HG_BC, pi_parents_by_location$norm_pi_CA_OR, alternative = "less")

wilcox.test(pi_with_gene_names$norm_pi_parents, pi_parents_by_location$norm_pi_VI_BC)
wilcox.test(pi_with_gene_names$norm_pi_parents, pi_parents_by_location$norm_pi_COAST_BC, alternative = "greater")
wilcox.test(pi_with_gene_names$norm_pi_parents, pi_parents_by_location$norm_pi_HG_BC, alternative = "greater")
wilcox.test(pi_with_gene_names$norm_pi_parents, pi_parents_by_location$norm_pi_INT_BC, alternative = "greater")
wilcox.test(pi_with_gene_names$norm_pi_parents, pi_parents_by_location$norm_pi_CA_OR, alternative = "greater")

x <- data.frame(pi_parents_by_location[10:14])
kruskal.test(x)


# nonsyn
pi_parents_by_location_nonsyn <- all_pi_parents_by_location_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  # filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>%
  mutate(norm_pi_VI_BC = V2 / gene_length, norm_pi_COAST_BC = V3 / gene_length, norm_pi_HG_BC = V4 / gene_length,
         norm_pi_INT_BC = V5 / gene_length, norm_pi_CA_OR = V6 / gene_length)

mean(pi_parents_by_location_nonsyn$norm_pi_VI_BC)
mean(pi_parents_by_location_nonsyn$norm_pi_COAST_BC)
mean(pi_parents_by_location_nonsyn$norm_pi_HG_BC)
mean(pi_parents_by_location_nonsyn$norm_pi_INT_BC)
mean(pi_parents_by_location_nonsyn$norm_pi_CA_OR)

# syn
pi_parents_by_location_syn <- all_pi_parents_by_location_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  # filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>%
  mutate(norm_pi_VI_BC = V2 / gene_length, norm_pi_COAST_BC = V3 / gene_length, norm_pi_HG_BC = V4 / gene_length,
         norm_pi_INT_BC = V5 / gene_length, norm_pi_CA_OR = V6 / gene_length)

mean(pi_parents_by_location_syn$norm_pi_VI_BC)
mean(pi_parents_by_location_syn$norm_pi_COAST_BC)
mean(pi_parents_by_location_syn$norm_pi_HG_BC)
mean(pi_parents_by_location_syn$norm_pi_INT_BC)
mean(pi_parents_by_location_syn$norm_pi_CA_OR)


wilcox.test(pi_parents_by_location_syn$norm_pi_VI_BC, pi_parents_by_location_nonsyn$norm_pi_VI_BC)

# 
# write.table(pi_with_gene_names$gene, "genes_in_scaffolds.txt", row.names = F, col.names = F, quote = F)
# 
# pi_syn <- all_pi_syn %>% 
#   filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>% 
#   mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
#   filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0 & V7 == 0 & V8 == 0)) %>% 
#   mutate(norm_pi_parents = V2 / gene_length, norm_pi_FS = V3 / gene_length, norm_pi_S1 = V4 / gene_length, 
#          norm_pi_S2 = V5 / gene_length, norm_pi_S3 = V6 / gene_length, norm_pi_S4 = V7 / gene_length, norm_pi_S5 = V8 / gene_length)
# pi_with_scaffold_names <- all_pi_parent_locations_scaffolds %>% 
#   mutate(scaffold = scaffold_names, scaffold_length = v3_scaffold_lengths_in_vcf$V2) %>%
#   filter(!(V2 == 0 & V3 == 0 & V4 == 0 & V5 == 0 & V6 == 0)) %>% 
#   mutate(norm_pi_VI = V2 / scaffold_length, norm_pi_COAST = V3 / scaffold_length, norm_pi_HG = V4 / scaffold_length, 
#          norm_pi_INT = V5 / scaffold_length, norm_pi_CA_OR = V6 / scaffold_length)
# 
# mean(pi_with_scaffold_names$norm_pi_VI)
# mean(pi_with_scaffold_names$norm_pi_COAST)
# mean(pi_with_scaffold_names$norm_pi_HG)
# mean(pi_with_scaffold_names$norm_pi_INT)
# mean(pi_with_scaffold_names$norm_pi_CA_OR)
# 
# pi_with_scaffold_names_melted <- reshape2::melt(pi_with_scaffold_names[,-c(1,8:13)])
# 
# ggplot(pi_with_scaffold_names_melted, aes(value, fill = variable, colour = variable)) +
#   geom_density(alpha = 0.1) #+
#   # xlim(0, 0.00001)
# 
# nucl_between_with_gene_names <- all_nucl_between_parent_locations_genes %>%
#   mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
#   filter(gene %in% genes_with_snps$V1) %>% 
#   mutate(norm_VI_BC.COAST_BC = VI_BC.COAST_BC / gene_length, norm_VI_BC.HG_BC = VI_BC.HG_BC / gene_length, 
#          norm_VI_BC.INT_BC = VI_BC.INT_BC / gene_length, norm_VI_BC.CA_OR = VI_BC.CA_OR / gene_length, 
#          norm_COAST_BC.HG_BC = COAST_BC.HG_BC / gene_length, norm_COAST_BC.INT_BC = COAST_BC.INT_BC / gene_length,
#          norm_COAST_BC.CA_OR = COAST_BC.CA_OR / gene_length, norm_HG_BC.INT_BC = HG_BC.INT_BC / gene_length,
#          norm_HG_BC.CA_OR = HG_BC.CA_OR / gene_length, norm_INT_BC.CA_OR = INT_BC.CA_OR / gene_length)
# 
# clipr::write_clip(data.frame(mean(nucl_between_with_gene_names$norm_VI_BC.COAST_BC, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_VI_BC.HG_BC, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_VI_BC.INT_BC, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_VI_BC.CA_OR, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_COAST_BC.HG_BC, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_COAST_BC.INT_BC, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_COAST_BC.CA_OR, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_HG_BC.INT_BC, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_HG_BC.CA_OR, na.rm = T),
#                              mean(nucl_between_with_gene_names$norm_INT_BC.CA_OR, na.rm = T)))

#### Neutrality stats ###############################
# #Parents, S generations
neutrality_parents <- all_neutrality_gene_parents %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_parents = V11 / gene_length)

neutrality_parents$norm_theta_parents[is.na(neutrality_parents$norm_theta_parents)] <- 0

neutrality_FS <- all_neutrality_gene_FS %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_FS = V11 / gene_length)

neutrality_FS$norm_theta_FS[is.na(neutrality_FS$norm_theta_FS)] <- 0

neutrality_S1 <- all_neutrality_gene_S1 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S1 = V11 / gene_length)

neutrality_S1$norm_theta_S1[is.na(neutrality_S1$norm_theta_S1)] <- 0

neutrality_S2 <- all_neutrality_gene_S2 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S2 = V11 / gene_length)

neutrality_S2$norm_theta_S2[is.na(neutrality_S2$norm_theta_S2)] <- 0

neutrality_S3 <- all_neutrality_gene_S3 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S3 = V11 / gene_length)

neutrality_S3$norm_theta_S3[is.na(neutrality_S3$norm_theta_S3)] <- 0

neutrality_S4 <- all_neutrality_gene_S4 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S4 = V11 / gene_length)

neutrality_S4$norm_theta_S4[is.na(neutrality_S4$norm_theta_S4)] <- 0

neutrality_S5 <- all_neutrality_gene_S5 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S5 = V11 / gene_length)

neutrality_S5$norm_theta_S5[is.na(neutrality_S5$norm_theta_S5)] <- 0

mean(neutrality_parents$norm_theta_parents, na.rm = T)
mean(neutrality_FS$norm_theta_FS, na.rm = T)
mean(neutrality_S1$norm_theta_S1, na.rm = T)
mean(neutrality_S2$norm_theta_S2, na.rm = T)
mean(neutrality_S3$norm_theta_S3, na.rm = T)
mean(neutrality_S4$norm_theta_S4, na.rm = T)
mean(neutrality_S5$norm_theta_S5, na.rm = T)

mean(neutrality_parents$V2, na.rm = T)
mean(neutrality_FS$V2, na.rm = T)
mean(neutrality_S1$V2, na.rm = T)
mean(neutrality_S2$V2, na.rm = T)
mean(neutrality_S3$V2, na.rm = T)
mean(neutrality_S4$V2, na.rm = T)
mean(neutrality_S5$V2, na.rm = T)

wilcox.test(neutrality_FS$norm_theta_FS, neutrality_S4$norm_theta_S4)

wilcox.test(neutrality_FS_syn$V2, neutrality_FS_nonsyn$V2, paired = T)
wilcox.test(neutrality_S1_syn$V2, neutrality_S1_nonsyn$V2, paired = T)
wilcox.test(neutrality_S2_syn$V2, neutrality_S2_nonsyn$V2, paired = T)
wilcox.test(neutrality_S3_syn$V2, neutrality_S3_nonsyn$V2, paired = T)
wilcox.test(neutrality_S4_syn$V2, neutrality_S4_nonsyn$V2, paired = T)
wilcox.test(neutrality_S5_syn$V2, neutrality_S5_nonsyn$V2, paired = T)


x <- data.frame(cbind(neutrality_FS$norm_theta_FS, neutrality_S1$norm_theta_S1,
                      neutrality_S2$norm_theta_S2, neutrality_S3$norm_theta_S3, neutrality_S4$norm_theta_S4))
kruskal.test(x)


outlier_D_parents <- neutrality_parents %>% 
  arrange(desc(V2)) %>% 
  top_frac(-0.01, V2)

write.table(outlier_D_parents, "outlier_D_parents.txt", quote = F, col.names = T, row.names = F)

FS_S5_D_values <- neutrality_FS %>% 
  select(V2, scaffold, gene_length, gene) %>% 
  mutate(S1_D = neutrality_S1$V2, S2_D = neutrality_S2$V2, S3_D = neutrality_S3$V2, S4_D = neutrality_S4$V2, S5_D = neutrality_S5$V2)

outlier_D_FS <- FS_S5_D_values %>% 
  arrange(desc(V2)) %>% 
  top_frac(-0.01, V2)

outlier_D_S4 <- FS_S5_D_values %>% 
  arrange(desc(S4_D)) %>% 
  top_frac(-0.01, S4_D)

write.table( outlier_D_S4, "outlier_D_S4.txt", col.names = T, row.names = F, quote = F)

outlier_D_S5 <- FS_S5_D_values %>% 
  arrange(desc(S5_D)) %>% 
  top_frac(-0.01, S5_D)

FS_S4_D_outlier_overlap <- outlier_D_FS %>% 
  filter(gene %in% outlier_D_S4$gene)

FS_S4_outlier_no_overlap <- outlier_D_S4 %>% 
  filter(!(gene %in% outlier_D_FS$gene))

FS_S5_decreasing <- FS_S5_D_values %>% 
  filter(V2 > S1_D & V2 > S2_D & V2 > S3_D & V2 > S4_D & S1_D > S2_D & S1_D > S3_D & S1_D > S4_D & S2_D > S3_D & S2_D > S4_D & S3_D > S4_D &
           S4_D < 0 & V2 > 0)

write.table(FS_S5_decreasing, "FS_S4_D_devreasing.txt", col.names = T, row.names = F, quote = F)
write.table(FS_S4_outlier_no_overlap, "FS_S4_outlier_no_overlap.txt", col.names = T, row.names = F, quote = F)

FS_S5_D <- data.frame(FS = neutrality_FS$V2, S1 = neutrality_S1$V2, S2 = neutrality_S2$V2, S3 = neutrality_S3$V2,
                      S4 = neutrality_S4$V2, S5 = neutrality_S5$V2)

FS_S5_D_syn_nonsyn <- data.frame(FS = c(neutrality_FS_syn$V2, neutrality_FS_nonsyn$V2), 
                                 S1 = c(neutrality_S1_syn$V2, neutrality_S1_nonsyn$V2), 
                                 S2 = c(neutrality_S2_syn$V2, neutrality_S2_nonsyn$V2),
                                 S3 = c(neutrality_S3_syn$V2, neutrality_S3_nonsyn$V2), 
                                 S4 = c(neutrality_S4_syn$V2, neutrality_S4_nonsyn$V2),
                                 S5 = c(neutrality_S5_syn$V2, neutrality_S5_nonsyn$V2),
                                 site = c(rep("syn", 33360), rep("nonsyn", 33360)))

melted_FS_S5_D <- melt(FS_S5_D)
melted_FS_S5_D_syn_nonsyn <- melt(FS_S5_D_syn_nonsyn)

ggplot(melted_FS_S5_D_syn_nonsyn, aes(x = variable, y = value, fill = site)) +
  geom_boxplot()

# nonsyn neutrality
neutrality_parents_nonsyn <- all_neutrality_gene_parents_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_parents_nonsyn = V11 / gene_length)

neutrality_parents_nonsyn$norm_theta_parents_nonsyn[is.na(neutrality_parents_nonsyn$norm_theta_parents_nonsyn)] <- 0

neutrality_FS_nonsyn <- all_neutrality_gene_FS_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_FS_nonsyn = V11 / gene_length)

neutrality_FS_nonsyn$norm_theta_FS_nonsyn[is.na(neutrality_FS_nonsyn$norm_theta_FS_nonsyn)] <- 0

neutrality_S1_nonsyn <- all_neutrality_gene_S1_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S1_nonsyn = V11 / gene_length)

neutrality_S1_nonsyn$norm_theta_S1_nonsyn[is.na(neutrality_S1_nonsyn$norm_theta_S1_nonsyn)] <- 0

neutrality_S2_nonsyn <- all_neutrality_gene_S2_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S2_nonsyn = V11 / gene_length)

neutrality_S2_nonsyn$norm_theta_S2_nonsyn[is.na(neutrality_S2_nonsyn$norm_theta_S2_nonsyn)] <- 0

neutrality_S3_nonsyn <- all_neutrality_gene_S3_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S3_nonsyn = V11 / gene_length)

neutrality_S3_nonsyn$norm_theta_S3_nonsyn[is.na(neutrality_S3_nonsyn$norm_theta_S3_nonsyn)] <- 0

neutrality_S4_nonsyn <- all_neutrality_gene_S4_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S4_nonsyn = V11 / gene_length)

neutrality_S4_nonsyn$norm_theta_S4_nonsyn[is.na(neutrality_S4_nonsyn$norm_theta_S4_nonsyn)] <- 0

neutrality_S5_nonsyn <- all_neutrality_gene_S5_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S5_nonsyn = V11 / gene_length)

neutrality_S5_nonsyn$norm_theta_S5_nonsyn[is.na(neutrality_S5_nonsyn$norm_theta_S5_nonsyn)] <- 0

mean(neutrality_parents_nonsyn$norm_theta_parents_nonsyn, na.rm = T)
mean(neutrality_FS_nonsyn$norm_theta_FS_nonsyn, na.rm = T)
mean(neutrality_S1_nonsyn$norm_theta_S1_nonsyn, na.rm = T)
mean(neutrality_S2_nonsyn$norm_theta_S2_nonsyn, na.rm = T)
mean(neutrality_S3_nonsyn$norm_theta_S3_nonsyn, na.rm = T)
mean(neutrality_S4_nonsyn$norm_theta_S4_nonsyn, na.rm = T)
mean(neutrality_S5_nonsyn$norm_theta_S5_nonsyn, na.rm = T)

mean(neutrality_parents_nonsyn$V2, na.rm = T)
mean(neutrality_FS_nonsyn$V2, na.rm = T)
mean(neutrality_S1_nonsyn$V2, na.rm = T)
mean(neutrality_S2_nonsyn$V2, na.rm = T)
mean(neutrality_S3_nonsyn$V2, na.rm = T)
mean(neutrality_S4_nonsyn$V2, na.rm = T)
mean(neutrality_S5_nonsyn$V2, na.rm = T)

# syn neutrality
neutrality_parents_syn <- all_neutrality_gene_parents_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_parents_syn = V11 / gene_length)

neutrality_parents_syn$norm_theta_parents_syn[is.na(neutrality_parents_syn$norm_theta_parents_syn)] <- 0

neutrality_FS_syn <- all_neutrality_gene_FS_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_FS_syn = V11 / gene_length)

neutrality_FS_syn$norm_theta_FS_syn[is.na(neutrality_FS_syn$norm_theta_FS_syn)] <- 0

neutrality_S1_syn <- all_neutrality_gene_S1_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S1_syn = V11 / gene_length)

neutrality_S1_syn$norm_theta_S1_syn[is.na(neutrality_S1_syn$norm_theta_S1_syn)] <- 0

neutrality_S2_syn <- all_neutrality_gene_S2_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S2_syn = V11 / gene_length)

neutrality_S2_syn$norm_theta_S2_syn[is.na(neutrality_S2_syn$norm_theta_S2_syn)] <- 0

neutrality_S3_syn <- all_neutrality_gene_S3_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S3_syn = V11 / gene_length)

neutrality_S3_syn$norm_theta_S3_syn[is.na(neutrality_S3_syn$norm_theta_S3_syn)] <- 0

neutrality_S4_syn <- all_neutrality_gene_S4_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S4_syn = V11 / gene_length)

neutrality_S4_syn$norm_theta_S4_syn[is.na(neutrality_S4_syn$norm_theta_S4_syn)] <- 0

neutrality_S5_syn <- all_neutrality_gene_S5_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_S5_syn = V11 / gene_length)

neutrality_S5_syn$norm_theta_S5_syn[is.na(neutrality_S5_syn$norm_theta_S5_syn)] <- 0

mean(neutrality_parents_syn$norm_theta_parents_syn, na.rm = T)
mean(neutrality_FS_syn$norm_theta_FS_syn, na.rm = T)
mean(neutrality_S1_syn$norm_theta_S1_syn, na.rm = T)
mean(neutrality_S2_syn$norm_theta_S2_syn, na.rm = T)
mean(neutrality_S3_syn$norm_theta_S3_syn, na.rm = T)
mean(neutrality_S4_syn$norm_theta_S4_syn, na.rm = T)
mean(neutrality_S5_syn$norm_theta_S5_syn, na.rm = T)

mean(neutrality_parents_syn$V2, na.rm = T)
mean(neutrality_FS_syn$V2, na.rm = T)
mean(neutrality_S1_syn$V2, na.rm = T)
mean(neutrality_S2_syn$V2, na.rm = T)
mean(neutrality_S3_syn$V2, na.rm = T)
mean(neutrality_S4_syn$V2, na.rm = T)
mean(neutrality_S5_syn$V2, na.rm = T)

wilcox.test(neutrality_parents_nonsyn$norm_theta_parents_nonsyn, neutrality_parents_syn$norm_theta_parents_syn, paired = T, alternative = "greater")

### Neutrality - parents locations #######
neutrality_VI_BC <- all_neutrality_gene_VI_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_VI_BC = V11 / gene_length)

neutrality_VI_BC$norm_theta_VI_BC[is.na(neutrality_VI_BC$norm_theta_VI_BC)] <- 0

neutrality_COAST_BC <- all_neutrality_gene_COAST_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_COAST_BC = V11 / gene_length)

neutrality_COAST_BC$norm_theta_COAST_BC[is.na(neutrality_COAST_BC$norm_theta_COAST_BC)] <- 0

neutrality_HG_BC <- all_neutrality_gene_HG_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_HG_BC = V11 / gene_length)

neutrality_HG_BC$norm_theta_HG_BC[is.na(neutrality_HG_BC$norm_theta_HG_BC)] <- 0

neutrality_INT_BC <- all_neutrality_gene_INT_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_INT_BC = V11 / gene_length)

neutrality_INT_BC$norm_theta_INT_BC[is.na(neutrality_INT_BC$norm_theta_INT_BC)] <- 0

neutrality_CA_OR <- all_neutrality_gene_CA_OR %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_CA_OR = V11 / gene_length)

neutrality_CA_OR$norm_theta_CA_OR[is.na(neutrality_CA_OR$norm_theta_CA_OR)] <- 0

mean(neutrality_VI_BC$norm_theta_VI_BC)
mean(neutrality_COAST_BC$norm_theta_COAST_BC)
mean(neutrality_HG_BC$norm_theta_HG_BC)
mean(neutrality_INT_BC$norm_theta_INT_BC)
mean(neutrality_CA_OR$norm_theta_CA_OR)

mean(neutrality_VI_BC$V2, na.rm = T)
mean(neutrality_COAST_BC$V2, na.rm = T)
mean(neutrality_HG_BC$V2, na.rm = T)
mean(neutrality_INT_BC$V2, na.rm = T)
mean(neutrality_CA_OR$V2, na.rm = T)

# nonsyn neutrality - parents by location
neutrality_VI_BC_nonsyn <- all_neutrality_gene_nonsyn_VI_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_VI_BC_nonsyn = V11 / gene_length)

neutrality_VI_BC_nonsyn$norm_theta_VI_BC_nonsyn[is.na(neutrality_VI_BC_nonsyn$norm_theta_VI_BC_nonsyn)] <- 0

neutrality_COAST_BC_nonsyn <- all_neutrality_gene_nonsyn_COAST_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_COAST_BC_nonsyn = V11 / gene_length)

neutrality_COAST_BC_nonsyn$norm_theta_COAST_BC_nonsyn[is.na(neutrality_COAST_BC_nonsyn$norm_theta_COAST_BC_nonsyn)] <- 0

neutrality_HG_BC_nonsyn <- all_neutrality_gene_nonsyn_HG_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_HG_BC_nonsyn = V11 / gene_length)

neutrality_HG_BC_nonsyn$norm_theta_HG_BC_nonsyn[is.na(neutrality_HG_BC_nonsyn$norm_theta_HG_BC_nonsyn)] <- 0

neutrality_INT_BC_nonsyn <- all_neutrality_gene_nonsyn_INT_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_INT_BC_nonsyn = V11 / gene_length)

neutrality_INT_BC_nonsyn$norm_theta_INT_BC_nonsyn[is.na(neutrality_INT_BC_nonsyn$norm_theta_INT_BC_nonsyn)] <- 0

neutrality_CA_OR_nonsyn <- all_neutrality_gene_nonsyn_CA_OR %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_CA_OR_nonsyn = V11 / gene_length)

neutrality_CA_OR_nonsyn$norm_theta_CA_OR_nonsyn[is.na(neutrality_CA_OR_nonsyn$norm_theta_CA_OR_nonsyn)] <- 0

mean(neutrality_VI_BC_nonsyn$norm_theta_VI_BC_nonsyn)
mean(neutrality_COAST_BC_nonsyn$norm_theta_COAST_BC_nonsyn)
mean(neutrality_HG_BC_nonsyn$norm_theta_HG_BC_nonsyn)
mean(neutrality_INT_BC_nonsyn$norm_theta_INT_BC_nonsyn)
mean(neutrality_CA_OR_nonsyn$norm_theta_CA_OR_nonsyn)

mean(neutrality_VI_BC_nonsyn$V2, na.rm = T)
mean(neutrality_COAST_BC_nonsyn$V2, na.rm = T)
mean(neutrality_HG_BC_nonsyn$V2, na.rm = T)
mean(neutrality_INT_BC_nonsyn$V2, na.rm = T)
mean(neutrality_CA_OR_nonsyn$V2, na.rm = T)

# syn neutrality - parents by locaation
neutrality_VI_BC_syn <- all_neutrality_gene_syn_VI_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_VI_BC_syn = V11 / gene_length)

neutrality_VI_BC_syn$norm_theta_VI_BC_syn[is.na(neutrality_VI_BC_syn$norm_theta_VI_BC_syn)] <- 0

neutrality_COAST_BC_syn <- all_neutrality_gene_syn_COAST_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_COAST_BC_syn = V11 / gene_length)

neutrality_COAST_BC_syn$norm_theta_COAST_BC_syn[is.na(neutrality_COAST_BC_syn$norm_theta_COAST_BC_syn)] <- 0

neutrality_HG_BC_syn <- all_neutrality_gene_syn_HG_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_HG_BC_syn = V11 / gene_length)

neutrality_HG_BC_syn$norm_theta_HG_BC_syn[is.na(neutrality_HG_BC_syn$norm_theta_HG_BC_syn)] <- 0

neutrality_INT_BC_syn <- all_neutrality_gene_syn_INT_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_INT_BC_syn = V11 / gene_length)

neutrality_INT_BC_syn$norm_theta_INT_BC_syn[is.na(neutrality_INT_BC_syn$norm_theta_INT_BC_syn)] <- 0

neutrality_CA_OR_syn <- all_neutrality_gene_syn_CA_OR %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_theta_CA_OR_syn = V11 / gene_length)

neutrality_CA_OR_syn$norm_theta_CA_OR_syn[is.na(neutrality_CA_OR_syn$norm_theta_CA_OR_syn)] <- 0

mean(neutrality_VI_BC_syn$norm_theta_VI_BC_syn)
mean(neutrality_COAST_BC_syn$norm_theta_COAST_BC_syn)
mean(neutrality_HG_BC_syn$norm_theta_HG_BC_syn)
mean(neutrality_INT_BC_syn$norm_theta_INT_BC_syn)
mean(neutrality_CA_OR_syn$norm_theta_CA_OR_syn)

mean(neutrality_VI_BC_syn$V2, na.rm = T)
mean(neutrality_COAST_BC_syn$V2, na.rm = T)
mean(neutrality_HG_BC_syn$V2, na.rm = T)
mean(neutrality_INT_BC_syn$V2, na.rm = T)
mean(neutrality_CA_OR_syn$V2, na.rm = T)

### Linkage all gens #####################
linkage_parents <- all_linkage_gene_parents %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_parents = V2 / gene_length)

# linkage_parents$norm_wall_b_parents[is.na(linkage_parents$norm_wall_b_parents)] <- 0

linkage_FS <- all_linkage_gene_FS %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_FS = V2 / gene_length)

# linkage_FS$norm_wall_b_FS[is.na(linkage_FS$norm_wall_b_FS)] <- 0

linkage_S1 <- all_linkage_gene_S1 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S1 = V2 / gene_length)

# linkage_S1$norm_wall_b_S1[is.na(linkage_S1$norm_wall_b_S1)] <- 0

linkage_S2 <- all_linkage_gene_S2 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S2 = V2 / gene_length)

# linkage_S2$norm_wall_b_S2[is.na(linkage_S2$norm_wall_b_S2)] <- 0

linkage_S3 <- all_linkage_gene_S3 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S3 = V2 / gene_length)

# linkage_S3$norm_wall_b_S3[is.na(linkage_S3$norm_wall_b_S3)] <- 0

linkage_S4 <- all_linkage_gene_S4 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S4 = V2 / gene_length)

# linkage_S4$norm_wall_b_S4[is.na(linkage_S4$norm_wall_b_S4)] <- 0

linkage_S5 <- all_linkage_gene_S5 %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S5 = V2 / gene_length)

# linkage_S5$norm_wall_b_S5[is.na(linkage_S5$norm_wall_b_S5)] <- 0

mean(linkage_parents$norm_wall_b_parents, na.rm = T)
mean(linkage_FS$norm_wall_b_FS, na.rm = T)
mean(linkage_S1$norm_wall_b_S1, na.rm = T)
mean(linkage_S2$norm_wall_b_S2, na.rm = T)
mean(linkage_S3$norm_wall_b_S3, na.rm = T)
mean(linkage_S4$norm_wall_b_S4, na.rm = T)
mean(linkage_S5$norm_wall_b_S5, na.rm = T)


wilcox.test(linkage_FS$norm_wall_b_FS, linkage_S1$norm_wall_b_S1, alternative = "less")
wilcox.test(linkage_FS$norm_wall_b_FS, linkage_S2$norm_wall_b_S2)
wilcox.test(linkage_FS$norm_wall_b_FS, linkage_S3$norm_wall_b_S3)
wilcox.test(linkage_FS$norm_wall_b_FS, linkage_S4$norm_wall_b_S4, alternative = "less")
wilcox.test(linkage_S1$norm_wall_b_S1, linkage_S2$norm_wall_b_S2)
wilcox.test(linkage_S1$norm_wall_b_S1, linkage_S3$norm_wall_b_S3)


wilcox.test(linkage_parents$norm_wall_b_parents, linkage_FS$norm_wall_b_FS)

x <- data.frame(cbind(linkage_FS$norm_wall_b_FS, linkage_S1$norm_wall_b_S1,
                      linkage_S2$norm_wall_b_S2, linkage_S3$norm_wall_b_S3, linkage_S4$norm_wall_b_S4, linkage_S5$norm_wall_b_S5))
kruskal.test(x)

# nonsyn linkage
linkage_parents_nonsyn <- all_linkage_gene_parents_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_parents_nonsyn = V2 / gene_length)

# linkage_parents_nonsyn$norm_wall_b_parents_nonsyn[is.na(linkage_parents_nonsyn$norm_wall_b_parents_nonsyn)] <- 0

linkage_FS_nonsyn <- all_linkage_gene_FS_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_FS_nonsyn = V2 / gene_length)

# linkage_FS_nonsyn$norm_wall_b_FS_nonsyn[is.na(linkage_FS_nonsyn$norm_wall_b_FS_nonsyn)] <- 0

linkage_S1_nonsyn <- all_linkage_gene_S1_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S1_nonsyn = V2 / gene_length)

# linkage_S1_nonsyn$norm_wall_b_S1_nonsyn[is.na(linkage_S1_nonsyn$norm_wall_b_S1_nonsyn)] <- 0

linkage_S2_nonsyn <- all_linkage_gene_S2_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S2_nonsyn = V2 / gene_length)

# linkage_S2_nonsyn$norm_wall_b_S2_nonsyn[is.na(linkage_S2_nonsyn$norm_wall_b_S2_nonsyn)] <- 0

linkage_S3_nonsyn <- all_linkage_gene_S3_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S3_nonsyn = V2 / gene_length)

# linkage_S3_nonsyn$norm_wall_b_S3_nonsyn[is.na(linkage_S3_nonsyn$norm_wall_b_S3_nonsyn)] <- 0

linkage_S4_nonsyn <- all_linkage_gene_S4_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S4_nonsyn = V2 / gene_length)

# linkage_S4_nonsyn$norm_wall_b_S4_nonsyn[is.na(linkage_S4_nonsyn$norm_wall_b_S4_nonsyn)] <- 0

linkage_S5_nonsyn <- all_linkage_gene_S5_nonsyn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S5_nonsyn = V2 / gene_length)

# linkage_S5_nonsyn$norm_wall_b_S5_nonsyn[is.na(linkage_S5_nonsyn$norm_wall_b_S5_nonsyn)] <- 0

mean(linkage_parents_nonsyn$norm_wall_b_parents_nonsyn, na.rm = T)
mean(linkage_FS_nonsyn$norm_wall_b_FS_nonsyn, na.rm = T)
mean(linkage_S1_nonsyn$norm_wall_b_S1_nonsyn, na.rm = T)
mean(linkage_S2_nonsyn$norm_wall_b_S2_nonsyn, na.rm = T)
mean(linkage_S3_nonsyn$norm_wall_b_S3_nonsyn, na.rm = T)
mean(linkage_S4_nonsyn$norm_wall_b_S4_nonsyn, na.rm = T)
mean(linkage_S5_nonsyn$norm_wall_b_S5_nonsyn, na.rm = T)

# syn linkage
linkage_parents_syn <- all_linkage_gene_parents_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_parents_syn = V2 / gene_length)

# linkage_parents_syn$norm_wall_b_parents_syn[is.na(linkage_parents_syn$norm_wall_b_parents_syn)] <- 0

linkage_FS_syn <- all_linkage_gene_FS_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_FS_syn = V2 / gene_length)

# linkage_FS_syn$norm_wall_b_FS_syn[is.na(linkage_FS_syn$norm_wall_b_FS_syn)] <- 0

linkage_S1_syn <- all_linkage_gene_S1_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S1_syn = V2 / gene_length)

# linkage_S1_syn$norm_wall_b_S1_syn[is.na(linkage_S1_syn$norm_wall_b_S1_syn)] <- 0

linkage_S2_syn <- all_linkage_gene_S2_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S2_syn = V2 / gene_length)

# linkage_S2_syn$norm_wall_b_S2_syn[is.na(linkage_S2_syn$norm_wall_b_S2_syn)] <- 0

linkage_S3_syn <- all_linkage_gene_S3_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S3_syn = V2 / gene_length)

# linkage_S3_syn$norm_wall_b_S3_syn[is.na(linkage_S3_syn$norm_wall_b_S3_syn)] <- 0

linkage_S4_syn <- all_linkage_gene_S4_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S4_syn = V2 / gene_length)

# linkage_S4_syn$norm_wall_b_S4_syn[is.na(linkage_S4_syn$norm_wall_b_S4_syn)] <- 0

linkage_S5_syn <- all_linkage_gene_S5_syn %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_S5_syn = V2 / gene_length)

# linkage_S5_syn$norm_wall_b_S5_syn[is.na(linkage_S5_syn$norm_wall_b_S5_syn)] <- 0

mean(linkage_parents_syn$norm_wall_b_parents_syn, na.rm = T)
mean(linkage_FS_syn$norm_wall_b_FS_syn, na.rm = T)
mean(linkage_S1_syn$norm_wall_b_S1_syn, na.rm = T)
mean(linkage_S2_syn$norm_wall_b_S2_syn, na.rm = T)
mean(linkage_S3_syn$norm_wall_b_S3_syn, na.rm = T)
mean(linkage_S4_syn$norm_wall_b_S4_syn, na.rm = T)
mean(linkage_S5_syn$norm_wall_b_S5_syn, na.rm = T)

mean(linkage_parents_syn$V2, na.rm = T)
mean(linkage_FS_syn$V2, na.rm = T)
mean(linkage_S1_syn$V2, na.rm = T)
mean(linkage_S2_syn$V2, na.rm = T)
mean(linkage_S3_syn$V2, na.rm = T)
mean(linkage_S4_syn$V2, na.rm = T)
mean(linkage_S5_syn$V2, na.rm = T)

wilcox.test(linkage_parents_syn$norm_wall_b_parents_syn, linkage_parents_nonsyn$norm_wall_b_parents_nonsyn, alternative = "less", paired = T)

### linkage - parents locations #######
linkage_VI_BC <- all_linkage_gene_VI_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_VI_BC = V2 / gene_length)

linkage_VI_BC$norm_wall_b_VI_BC[is.na(linkage_VI_BC$norm_wall_b_VI_BC)] <- 0

linkage_COAST_BC <- all_linkage_gene_COAST_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_COAST_BC = V2 / gene_length)

linkage_COAST_BC$norm_wall_b_COAST_BC[is.na(linkage_COAST_BC$norm_wall_b_COAST_BC)] <- 0

linkage_HG_BC <- all_linkage_gene_HG_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_HG_BC = V2 / gene_length)

linkage_HG_BC$norm_wall_b_HG_BC[is.na(linkage_HG_BC$norm_wall_b_HG_BC)] <- 0

linkage_INT_BC <- all_linkage_gene_INT_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_INT_BC = V2 / gene_length)

linkage_INT_BC$norm_wall_b_INT_BC[is.na(linkage_INT_BC$norm_wall_b_INT_BC)] <- 0

linkage_CA_OR <- all_linkage_gene_CA_OR %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_CA_OR = V2 / gene_length)

linkage_CA_OR$norm_wall_b_CA_OR[is.na(linkage_CA_OR$norm_wall_b_CA_OR)] <- 0

mean(linkage_VI_BC$norm_wall_b_VI_BC)
mean(linkage_COAST_BC$norm_wall_b_COAST_BC)
mean(linkage_HG_BC$norm_wall_b_HG_BC)
mean(linkage_INT_BC$norm_wall_b_INT_BC)
mean(linkage_CA_OR$norm_wall_b_CA_OR)

# nonsyn linkage - parents by location
linkage_VI_BC_nonsyn <- all_linkage_gene_nonsyn_VI_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_VI_BC_nonsyn = V2 / gene_length)

linkage_VI_BC_nonsyn$norm_wall_b_VI_BC_nonsyn[is.na(linkage_VI_BC_nonsyn$norm_wall_b_VI_BC_nonsyn)] <- 0

linkage_COAST_BC_nonsyn <- all_linkage_gene_nonsyn_COAST_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_COAST_BC_nonsyn = V2 / gene_length)

linkage_COAST_BC_nonsyn$norm_wall_b_COAST_BC_nonsyn[is.na(linkage_COAST_BC_nonsyn$norm_wall_b_COAST_BC_nonsyn)] <- 0

linkage_HG_BC_nonsyn <- all_linkage_gene_nonsyn_HG_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_HG_BC_nonsyn = V2 / gene_length)

linkage_HG_BC_nonsyn$norm_wall_b_HG_BC_nonsyn[is.na(linkage_HG_BC_nonsyn$norm_wall_b_HG_BC_nonsyn)] <- 0

linkage_INT_BC_nonsyn <- all_linkage_gene_nonsyn_INT_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_INT_BC_nonsyn = V2 / gene_length)

linkage_INT_BC_nonsyn$norm_wall_b_INT_BC_nonsyn[is.na(linkage_INT_BC_nonsyn$norm_wall_b_INT_BC_nonsyn)] <- 0

linkage_CA_OR_nonsyn <- all_linkage_gene_nonsyn_CA_OR %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_CA_OR_nonsyn = V2 / gene_length)

linkage_CA_OR_nonsyn$norm_wall_b_CA_OR_nonsyn[is.na(linkage_CA_OR_nonsyn$norm_wall_b_CA_OR_nonsyn)] <- 0

mean(linkage_VI_BC_nonsyn$norm_wall_b_VI_BC_nonsyn)
mean(linkage_COAST_BC_nonsyn$norm_wall_b_COAST_BC_nonsyn)
mean(linkage_HG_BC_nonsyn$norm_wall_b_HG_BC_nonsyn)
mean(linkage_INT_BC_nonsyn$norm_wall_b_INT_BC_nonsyn)
mean(linkage_CA_OR_nonsyn$norm_wall_b_CA_OR_nonsyn)

# syn linkage - parents by locaation
linkage_VI_BC_syn <- all_linkage_gene_syn_VI_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_VI_BC_syn = V2 / gene_length)

linkage_VI_BC_syn$norm_wall_b_VI_BC_syn[is.na(linkage_VI_BC_syn$norm_wall_b_VI_BC_syn)] <- 0

linkage_COAST_BC_syn <- all_linkage_gene_syn_COAST_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_COAST_BC_syn = V2 / gene_length)

linkage_COAST_BC_syn$norm_wall_b_COAST_BC_syn[is.na(linkage_COAST_BC_syn$norm_wall_b_COAST_BC_syn)] <- 0

linkage_HG_BC_syn <- all_linkage_gene_syn_HG_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_HG_BC_syn = V2 / gene_length)

linkage_HG_BC_syn$norm_wall_b_HG_BC_syn[is.na(linkage_HG_BC_syn$norm_wall_b_HG_BC_syn)] <- 0

linkage_INT_BC_syn <- all_linkage_gene_syn_INT_BC %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_INT_BC_syn = V2 / gene_length)

linkage_INT_BC_syn$norm_wall_b_INT_BC_syn[is.na(linkage_INT_BC_syn$norm_wall_b_INT_BC_syn)] <- 0

linkage_CA_OR_syn <- all_linkage_gene_syn_CA_OR %>%
  mutate(scaffold = gff_genes_correct_order$V1, gene_length = gff_genes_correct_order$V5 - gff_genes_correct_order$V4, gene = gff_genes_correct_order$V9) %>%
  mutate(norm_wall_b_CA_OR_syn = V2 / gene_length)

linkage_CA_OR_syn$norm_wall_b_CA_OR_syn[is.na(linkage_CA_OR_syn$norm_wall_b_CA_OR_syn)] <- 0

mean(linkage_VI_BC_syn$norm_wall_b_VI_BC_syn)
mean(linkage_COAST_BC_syn$norm_wall_b_COAST_BC_syn)
mean(linkage_HG_BC_syn$norm_wall_b_HG_BC_syn)
mean(linkage_INT_BC_syn$norm_wall_b_INT_BC_syn)
mean(linkage_CA_OR_syn$norm_wall_b_CA_OR_syn)

