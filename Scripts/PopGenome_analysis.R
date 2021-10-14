rm(list = ls())

## Script to obtain population genomics statistics using the PopGenome package for R ##

## Preliminary steps
# Load library PopGenome
library(PopGenome)
library(reshape2)
library(tidyverse)

path <- "sco_gene_scaffolds/"

# Get list of scaffold lengths in VCF
v3_scaffold_lengths <- read.table("genome_scaffold_lengths_dashes_removed.txt")

# Read in scaffold names
scaffold_names <- scan("sco_gene_scaffolds.txt", what = "factor")

# scaffold_names <- factor(scaffold_names)

v3_scaffold_lengths_in_vcf <- v3_scaffold_lengths %>% 
  filter(V1 %in% scaffold_names)



# read in GFF as table
gff <- read.table("Tplicatav3.1c.primaryTrs_dashes_removed_sorted.gff3")

# Filtering GFF file to match VCF
gff_scaffolds_with_genes <- gff %>% 
  filter(V1 %in% scaffold_names)

write.table(gff_scaffolds_with_genes, "gff_sco_genes.gff3", quote = F, col.names = F, row.names = F, sep = "\t")
GFF_split_into_scaffolds("gff_sco_genes.gff3", "gff_sco_genes")

# Scaffold length
# Read in scaffold lengths
# scaffold_lengths <- read.table("scaffold_lengths_scaffolds_in_gff_and_vcf_v3.1c.txt")
scaffold_lengths_with_genes <- v3_scaffold_lengths %>% 
  filter(V1 %in% scaffold_names) %>% 
  arrange(V1)

# write.table(scaffold_lengths_gff_vcf, "v3_scaffold_lengths_in_vcf_gff_1456_mislabeled_S_lines_removed.txt", quote = F, col.names = F, row.names = F)

scaffold_lengths <- as.numeric(scaffold_lengths_with_genes$V2)
length = c(scaffold_lengths)



# Analysis for diverse population - by gene  
for (i in 1:length(scaffold_names)){
  GENOME.class <- readVCF(paste(path, "VCF_", scaffold_names[i], ".recode.vcf.gz", sep = ""), numcols = 1000000, 
                          tid = scaffold_names[i], frompos=1, topos=length[i], approx=FALSE, out="", parallel=FALSE, 
                          include.unknown=TRUE, gffpath = paste(path, scaffold_names[i], sep = ""))
  # genes <- split_data_into_GFF_features(GENOME.class, gff.file = paste(path, scaffold_names[i], sep = ""),
  # chr = scaffold_names[i], feature = "gene")
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
  # genes <- neutrality.stats(genes, FAST=TRUE)
  # write.table(get.neutrality(genes, theta=TRUE)[[1]],
  #             file=paste("Neutrality_stats_gene_new/neutrality_stats_parents_sc",
  #                        scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
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
  # genes <- linkage.stats(genes)
  # write.table(get.linkage(genes)[[1]], file=paste("Linkage_stats_genes/linkage_stats_parents_sc",
  #                                                 scaffold_names[i], "_by_gene.txt", sep=""), quote=F, sep='\t')
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
     write.table(genes@Pi, file=paste("FST_stats_sco/pi_nucleotide_diversity_within_populations_sc", scaffold_names[i],
                                      "_by_gene.txt", sep=""), quote=F, sep='\t')
}
