rm(list=ls())

library(bigsnpr)
library(tidyverse)

setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/")

bedfile <- "LD_pruning_data/diverse_pop.g95minmeanDP15maxmeanDP60minQ30AB28.bed"
plink2 <- download_plink2("LD_pruning_data")

# file <- snp_fastImpute(bedfile)
# 
# G <- bedfile$genotypes
# CHR <- bedfile$map$chromosome
# POS <- bedfile$map$physical.pos
# NCORES <- nb_cores()
# G <- snp_fastImputeSimple(Gna = G, ncores = NCORES, method = "mean2")
# 
# 
# big_counts(G, ind.col = 1:12)

rel <- snp_plinkKINGQC(
  plink2.path = plink2,
  bedfile.in = bedfile,
  thr.king = 2^-3.5,
  make.bed = FALSE,
  ncores = nb_cores(),
  extra.options = "--allow-extra-chr"
)
str(rel)

# svd1 <- big_randomSVD(G, snp_scaleBinom(), ncores = NCORES)

(obj.bed <- bed(bedfile))
# (obj.bed2 <- bed(bedfile2))

ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)
# ind.rel2 <- match(c(rel2$ID1, rel2$ID2), obj.bed2$fam$sample.ID)

ind.norel <- rows_along(obj.bed)[-ind.rel]
# ind.norel2 <- rows_along(obj.bed2)[-ind.rel2]

obj.clumping <- bed_clumping(obj.bed, ind.row = ind.rel, 
                       ncores = nb_cores(), thr.r2 = 0.2)
# obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.norel, k = 20, 
                       # ncores = nb_cores())

bimfile <- read.table("LD_pruning_data/diverse_pop.g95minmeanDP15maxmeanDP60minQ30AB28.bim")
pruned <- bimfile$V2[obj.clumping]

write.table(pruned, "LD_pruned_snp_IDs_diverse_r2_0.2.txt", quote = F, row.names = F, col.names = F)

