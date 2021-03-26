library(adegenet)
library(hierfstat)
library(vcfR)
library(pegas)
library(reshape2)
library(ggsci)
library(apex)
library(ggpubr)
library(mmod)
library(factoextra)
library(FactoMineR)
library(poppr)
library(ggpubr)
library(rstatix)
library(tidyverse)

source("~/UBC/GSAT/PhD/WRC/r_scripts/publication_theme.r")

setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/")

### Corrected allele matrices for pop gen, heterozygosity #####
F_lines_allele_matrix <- t(read.table("S_lines_pop_gen/F_lines_alleles.txt", check.names = F))
S1_lines_allele_matrix <- t(read.table("S_lines_pop_gen/S1_lines_alleles.txt", check.names = F))
S2_lines_allele_matrix <- t(read.table("S_lines_pop_gen/S2_lines_alleles.txt", check.names = F))
S3_lines_allele_matrix <- t(read.table("S_lines_pop_gen/S3_lines_alleles.txt", check.names = F))
S4_lines_allele_matrix <- t(read.table("S_lines_pop_gen/S4_lines_alleles.txt", check.names = F))
S5_lines_allele_matrix <- t(read.table("S_lines_pop_gen/S5_lines_alleles.txt", check.names = F))

all_allele_matrix <- rbind(F_lines_allele_matrix, S1_lines_allele_matrix, S2_lines_allele_matrix, S3_lines_allele_matrix, S4_lines_allele_matrix, S5_lines_allele_matrix)
FS_S4_allele_matrix <- rbind(F_lines_allele_matrix, S1_lines_allele_matrix, S2_lines_allele_matrix, S3_lines_allele_matrix, S4_lines_allele_matrix)

# Popmap files
popmap_generations <- read.table("../popmap_generations.txt")
popmap_FS_lines <- read.table("popmap_FS_lines_new.txt")
popmap_parental_lines <- read.table("popmap_parental_lines.txt")

# Convert to genind
FS_S5_genind <- df2genind(all_allele_matrix, sep = "/", NA.char = "NA")
FS_S4_genind <- df2genind(FS_S4_allele_matrix, sep = "/", NA.char = "NA")
F_lines_genind <- FS_S5_genind[1:28]
S1_lines_genind <- FS_S5_genind[29:56]
S2_lines_genind <- FS_S5_genind[57:84]
S3_lines_genind <- FS_S5_genind[85:112]
S4_lines_genind <- FS_S5_genind[113:140]
S5_lines_genind <- FS_S5_genind[141:151]

Line_1_1_genind <- FS_S5_genind[c(1, 21, 41, 61, 81, 101)]

Line_1_1_F_genind <- FS_S5_genind[1]
Line_1_1_S1_genind <- FS_S5_genind[29]
Line_1_1_S2_genind <- FS_S5_genind[57]
Line_1_1_S3_genind <- FS_S5_genind[85]
Line_1_1_S4_genind <- FS_S5_genind[113]
Line_1_1_S5_genind <- FS_S5_genind[141]

Line_6_1_F_genind <- FS_S5_genind[2]
Line_6_1_S1_genind <- FS_S5_genind[30]
Line_6_1_S2_genind <- FS_S5_genind[58]
Line_6_1_S3_genind <- FS_S5_genind[86]
Line_6_1_S4_genind <- FS_S5_genind[114]

Line_6_4_F_genind <- FS_S5_genind[3]
Line_6_4_S1_genind <- FS_S5_genind[31]
Line_6_4_S2_genind <- FS_S5_genind[59]
Line_6_4_S3_genind <- FS_S5_genind[87]
Line_6_4_S4_genind <- FS_S5_genind[115]
Line_6_4_S5_genind <- FS_S5_genind[142]

Line_7_2_F_genind <- FS_S5_genind[4]
Line_7_2_S1_genind <- FS_S5_genind[32]
Line_7_2_S2_genind <- FS_S5_genind[60]
Line_7_2_S3_genind <- FS_S5_genind[88]
Line_7_2_S4_genind <- FS_S5_genind[116]

Line_7_4_F_genind <- FS_S5_genind[5]
Line_7_4_S1_genind <- FS_S5_genind[33]
Line_7_4_S2_genind <- FS_S5_genind[61]
Line_7_4_S3_genind <- FS_S5_genind[89]
Line_7_4_S4_genind <- FS_S5_genind[117]

Line_8_2_F_genind <- FS_S5_genind[6]
Line_8_2_S1_genind <- FS_S5_genind[34]
Line_8_2_S2_genind <- FS_S5_genind[62]
Line_8_2_S3_genind <- FS_S5_genind[90]
Line_8_2_S4_genind <- FS_S5_genind[118]

Line_8_4_F_genind <- FS_S5_genind[7]
Line_8_4_S1_genind <- FS_S5_genind[35]
Line_8_4_S2_genind <- FS_S5_genind[63]
Line_8_4_S3_genind <- FS_S5_genind[91]
Line_8_4_S4_genind <- FS_S5_genind[119]

Line_10_2_F_genind <- FS_S5_genind[8]
Line_10_2_S1_genind <- FS_S5_genind[36]
Line_10_2_S2_genind <- FS_S5_genind[64]
Line_10_2_S3_genind <- FS_S5_genind[92]
Line_10_2_S4_genind <- FS_S5_genind[120]
Line_10_2_S5_genind <- FS_S5_genind[143]

Line_10_5_F_genind <- FS_S5_genind[9]
Line_10_5_S1_genind <- FS_S5_genind[37]
Line_10_5_S2_genind <- FS_S5_genind[65]
Line_10_5_S3_genind <- FS_S5_genind[93]
Line_10_5_S4_genind <- FS_S5_genind[121]

Line_13_3_F_genind <- FS_S5_genind[10]
Line_13_3_S1_genind <- FS_S5_genind[38]
Line_13_3_S2_genind <- FS_S5_genind[66]
Line_13_3_S3_genind <- FS_S5_genind[94]
Line_13_3_S4_genind <- FS_S5_genind[122]

Line_13_4_F_genind <- FS_S5_genind[11]
Line_13_4_S1_genind <- FS_S5_genind[39]
Line_13_4_S2_genind <- FS_S5_genind[67]
Line_13_4_S3_genind <- FS_S5_genind[95]
Line_13_4_S4_genind <- FS_S5_genind[123]

Line_16_1_F_genind <- FS_S5_genind[12]
Line_16_1_S1_genind <- FS_S5_genind[40]
Line_16_1_S2_genind <- FS_S5_genind[68]
Line_16_1_S3_genind <- FS_S5_genind[96]
Line_16_1_S4_genind <- FS_S5_genind[124]

Line_16_5_F_genind <- FS_S5_genind[13]
Line_16_5_S1_genind <- FS_S5_genind[41]
Line_16_5_S2_genind <- FS_S5_genind[69]
Line_16_5_S3_genind <- FS_S5_genind[97]
Line_16_5_S4_genind <- FS_S5_genind[125]
Line_16_5_S5_genind <- FS_S5_genind[144]

Line_17_2_F_genind <- FS_S5_genind[14]
Line_17_2_S1_genind <- FS_S5_genind[42]
Line_17_2_S2_genind <- FS_S5_genind[70]
Line_17_2_S3_genind <- FS_S5_genind[98]
Line_17_2_S4_genind <- FS_S5_genind[126]
Line_17_2_S5_genind <- FS_S5_genind[145]

Line_17_5_F_genind <- FS_S5_genind[15]
Line_17_5_S1_genind <- FS_S5_genind[43]
Line_17_5_S2_genind <- FS_S5_genind[71]
Line_17_5_S3_genind <- FS_S5_genind[99]
Line_17_5_S4_genind <- FS_S5_genind[127]
Line_17_5_S5_genind <- FS_S5_genind[146]

Line_19_2_F_genind <- FS_S5_genind[16]
Line_19_2_S1_genind <- FS_S5_genind[44]
Line_19_2_S2_genind <- FS_S5_genind[72]
Line_19_2_S3_genind <- FS_S5_genind[100]
Line_19_2_S4_genind <- FS_S5_genind[128]

Line_19_5_F_genind <- FS_S5_genind[17]
Line_19_5_S1_genind <- FS_S5_genind[45]
Line_19_5_S2_genind <- FS_S5_genind[73]
Line_19_5_S3_genind <- FS_S5_genind[101]
Line_19_5_S4_genind <- FS_S5_genind[129]

Line_20_1_F_genind <- FS_S5_genind[18]
Line_20_1_S1_genind <- FS_S5_genind[46]
Line_20_1_S2_genind <- FS_S5_genind[74]
Line_20_1_S3_genind <- FS_S5_genind[102]
Line_20_1_S4_genind <- FS_S5_genind[130]

Line_20_4_F_genind <- FS_S5_genind[19]
Line_20_4_S1_genind <- FS_S5_genind[47]
Line_20_4_S2_genind <- FS_S5_genind[75]
Line_20_4_S3_genind <- FS_S5_genind[103]
Line_20_4_S4_genind <- FS_S5_genind[131]

Line_21_2_F_genind <- FS_S5_genind[20]
Line_21_2_S1_genind <- FS_S5_genind[48]
Line_21_2_S2_genind <- FS_S5_genind[76]
Line_21_2_S3_genind <- FS_S5_genind[104]
Line_21_2_S4_genind <- FS_S5_genind[132]
Line_21_2_S5_genind <- FS_S5_genind[147]

Line_21_6_F_genind <- FS_S5_genind[21]
Line_21_6_S1_genind <- FS_S5_genind[49]
Line_21_6_S2_genind <- FS_S5_genind[77]
Line_21_6_S3_genind <- FS_S5_genind[105]
Line_21_6_S4_genind <- FS_S5_genind[133]
Line_21_6_S5_genind <- FS_S5_genind[148]

Line_23_2_F_genind <- FS_S5_genind[22]
Line_23_2_S1_genind <- FS_S5_genind[50]
Line_23_2_S2_genind <- FS_S5_genind[78]
Line_23_2_S3_genind <- FS_S5_genind[106]
Line_23_2_S4_genind <- FS_S5_genind[134]
Line_23_2_S5_genind <- FS_S5_genind[149]

Line_23_4_F_genind <- FS_S5_genind[23]
Line_23_4_S1_genind <- FS_S5_genind[51]
Line_23_4_S2_genind <- FS_S5_genind[79]
Line_23_4_S3_genind <- FS_S5_genind[107]
Line_23_4_S4_genind <- FS_S5_genind[135]

Line_26_1_F_genind <- FS_S5_genind[24]
Line_26_1_S1_genind <- FS_S5_genind[52]
Line_26_1_S2_genind <- FS_S5_genind[80]
Line_26_1_S3_genind <- FS_S5_genind[108]
Line_26_1_S4_genind <- FS_S5_genind[136]

Line_26_4_F_genind <- FS_S5_genind[25]
Line_26_4_S1_genind <- FS_S5_genind[53]
Line_26_4_S2_genind <- FS_S5_genind[81]
Line_26_4_S3_genind <- FS_S5_genind[109]
Line_26_4_S4_genind <- FS_S5_genind[137]

Line_27_2_F_genind <- FS_S5_genind[26]
Line_27_2_S1_genind <- FS_S5_genind[54]
Line_27_2_S2_genind <- FS_S5_genind[82]
Line_27_2_S3_genind <- FS_S5_genind[110]
Line_27_2_S4_genind <- FS_S5_genind[138]

Line_29_2_F_genind <- FS_S5_genind[27]
Line_29_2_S1_genind <- FS_S5_genind[55]
Line_29_2_S2_genind <- FS_S5_genind[83]
Line_29_2_S3_genind <- FS_S5_genind[111]
Line_29_2_S4_genind <- FS_S5_genind[139]
Line_29_2_S5_genind <- FS_S5_genind[150]

Line_29_4_F_genind <- FS_S5_genind[28]
Line_29_4_S1_genind <- FS_S5_genind[56]
Line_29_4_S2_genind <- FS_S5_genind[84]
Line_29_4_S3_genind <- FS_S5_genind[112]
Line_29_4_S4_genind <- FS_S5_genind[140]
Line_29_4_S5_genind <- FS_S5_genind[151]


# FS_S4_clusters <- find.clusters(FS_S4_genind)
# FS_S4_dapc <- dapc(FS_S4_genind, FS_S4_clusters$grp)
# scatter(FS_S4_dapc)
# 
# FS_S5_stats <- genind2hierfstat(FS_S5_genind, pop = popmap_generations$V2)
# FS_S4_stats_FS_lines <- genind2hierfstat(FS_S4_genind, pop = popmap_FS_lines$V2)
# FS_S5_stats_parental_lines <- genind2hierfstat(FS_S4_genind, pop = popmap_parental_lines$V2)
# 
# # Genotype PCA
# x <- indpca(FS_S5_stats)
# 
# pdf(file = "FS-S5_generations_PCA_small.pdf")
# plot(x, cex = 0.2)
# dev.off()
# 
# y <- indpca(FS_S5_stats_FS_lines)
# 
# pdf(file = "FS-S5_lines_PCA_small.pdf")
# plot(y, cex = 0.2)
# dev.off()
# 
# z <- indpca(FS_S5_stats_parental_lines)
# 
# pdf(file = "FS-S5_parental_lines_PCA_small.pdf")
# plot(z, cex = 0.2)
# dev.off()

# Get summary for heterozygosity
FS_S4_div <- summary(FS_S4_genind) 
F_lines_div <- summary(F_lines_genind)
S1_lines_div <- summary(S1_lines_genind)
S2_lines_div <- summary(S2_lines_genind)
S3_lines_div <- summary(S3_lines_genind)
S4_lines_div <- summary(S4_lines_genind)
S5_lines_div <- summary(S5_lines_genind)

Line_1_1_F_div <- summary(Line_1_1_F_genind)
Line_1_1_S1_div <- summary(Line_1_1_S1_genind)
Line_1_1_S2_div <- summary(Line_1_1_S2_genind)
Line_1_1_S3_div <- summary(Line_1_1_S3_genind)
Line_1_1_S4_div <- summary(Line_1_1_S4_genind)
Line_1_1_S5_div <- summary(Line_1_1_S5_genind)

Line_6_1_F_div <- summary(Line_6_1_F_genind)
Line_6_1_S1_div <- summary(Line_6_1_S1_genind)
Line_6_1_S2_div <- summary(Line_6_1_S2_genind)
Line_6_1_S3_div <- summary(Line_6_1_S3_genind)
Line_6_1_S4_div <- summary(Line_6_1_S4_genind)

Line_6_4_F_div <- summary(Line_6_4_F_genind)
Line_6_4_S1_div <- summary(Line_6_4_S1_genind)
Line_6_4_S2_div <- summary(Line_6_4_S2_genind)
Line_6_4_S3_div <- summary(Line_6_4_S3_genind)
Line_6_4_S4_div <- summary(Line_6_4_S4_genind)
Line_6_4_S5_div <- summary(Line_6_4_S5_genind)

Line_7_2_F_div <- summary(Line_7_2_F_genind)
Line_7_2_S1_div <- summary(Line_7_2_S1_genind)
Line_7_2_S2_div <- summary(Line_7_2_S2_genind)
Line_7_2_S3_div <- summary(Line_7_2_S3_genind)
Line_7_2_S4_div <- summary(Line_7_2_S4_genind)

Line_7_4_F_div <- summary(Line_7_4_F_genind)
Line_7_4_S1_div <- summary(Line_7_4_S1_genind)
Line_7_4_S2_div <- summary(Line_7_4_S2_genind)
Line_7_4_S3_div <- summary(Line_7_4_S3_genind)
Line_7_4_S4_div <- summary(Line_7_4_S4_genind)

Line_8_2_F_div <- summary(Line_8_2_F_genind)
Line_8_2_S1_div <- summary(Line_8_2_S1_genind)
Line_8_2_S2_div <- summary(Line_8_2_S2_genind)
Line_8_2_S3_div <- summary(Line_8_2_S3_genind)
Line_8_2_S4_div <- summary(Line_8_2_S4_genind)

Line_8_4_F_div <- summary(Line_8_4_F_genind)
Line_8_4_S1_div <- summary(Line_8_4_S1_genind)
Line_8_4_S2_div <- summary(Line_8_4_S2_genind)
Line_8_4_S3_div <- summary(Line_8_4_S3_genind)
Line_8_4_S4_div <- summary(Line_8_4_S4_genind)

Line_10_2_F_div <- summary(Line_10_2_F_genind)
Line_10_2_S1_div <- summary(Line_10_2_S1_genind)
Line_10_2_S2_div <- summary(Line_10_2_S2_genind)
Line_10_2_S3_div <- summary(Line_10_2_S3_genind)
Line_10_2_S4_div <- summary(Line_10_2_S4_genind)
Line_10_2_S5_div <- summary(Line_10_2_S5_genind)

Line_10_5_F_div <- summary(Line_10_5_F_genind)
Line_10_5_S1_div <- summary(Line_10_5_S1_genind)
Line_10_5_S2_div <- summary(Line_10_5_S2_genind)
Line_10_5_S3_div <- summary(Line_10_5_S3_genind)
Line_10_5_S4_div <- summary(Line_10_5_S4_genind)

Line_13_3_F_div <- summary(Line_13_3_F_genind)
Line_13_3_S1_div <- summary(Line_13_3_S1_genind)
Line_13_3_S2_div <- summary(Line_13_3_S2_genind)
Line_13_3_S3_div <- summary(Line_13_3_S3_genind)
Line_13_3_S4_div <- summary(Line_13_3_S4_genind)

Line_13_4_F_div <- summary(Line_13_4_F_genind)
Line_13_4_S1_div <- summary(Line_13_4_S1_genind)
Line_13_4_S2_div <- summary(Line_13_4_S2_genind)
Line_13_4_S3_div <- summary(Line_13_4_S3_genind)
Line_13_4_S4_div <- summary(Line_13_4_S4_genind)

Line_16_1_F_div <- summary(Line_16_1_F_genind)
Line_16_1_S1_div <- summary(Line_16_1_S1_genind)
Line_16_1_S2_div <- summary(Line_16_1_S2_genind)
Line_16_1_S3_div <- summary(Line_16_1_S3_genind)
Line_16_1_S4_div <- summary(Line_16_1_S4_genind)

Line_16_5_F_div <- summary(Line_16_5_F_genind)
Line_16_5_S1_div <- summary(Line_16_5_S1_genind)
Line_16_5_S2_div <- summary(Line_16_5_S2_genind)
Line_16_5_S3_div <- summary(Line_16_5_S3_genind)
Line_16_5_S4_div <- summary(Line_16_5_S4_genind)
Line_16_5_S5_div <- summary(Line_16_5_S5_genind)

Line_17_2_F_div <- summary(Line_17_2_F_genind)
Line_17_2_S1_div <- summary(Line_17_2_S1_genind)
Line_17_2_S2_div <- summary(Line_17_2_S2_genind)
Line_17_2_S3_div <- summary(Line_17_2_S3_genind)
Line_17_2_S4_div <- summary(Line_17_2_S4_genind)
Line_17_2_S5_div <- summary(Line_17_2_S5_genind)

Line_17_5_F_div <- summary(Line_17_5_F_genind)
Line_17_5_S1_div <- summary(Line_17_5_S1_genind)
Line_17_5_S2_div <- summary(Line_17_5_S2_genind)
Line_17_5_S3_div <- summary(Line_17_5_S3_genind)
Line_17_5_S4_div <- summary(Line_17_5_S4_genind)
Line_17_5_S5_div <- summary(Line_17_5_S5_genind)

Line_19_2_F_div <- summary(Line_19_2_F_genind)
Line_19_2_S1_div <- summary(Line_19_2_S1_genind)
Line_19_2_S2_div <- summary(Line_19_2_S2_genind)
Line_19_2_S3_div <- summary(Line_19_2_S3_genind)
Line_19_2_S4_div <- summary(Line_19_2_S4_genind)

Line_19_5_F_div <- summary(Line_19_5_F_genind)
Line_19_5_S1_div <- summary(Line_19_5_S1_genind)
Line_19_5_S2_div <- summary(Line_19_5_S2_genind)
Line_19_5_S3_div <- summary(Line_19_5_S3_genind)
Line_19_5_S4_div <- summary(Line_19_5_S4_genind)

Line_20_1_F_div <- summary(Line_20_1_F_genind)
Line_20_1_S1_div <- summary(Line_20_1_S1_genind)
Line_20_1_S2_div <- summary(Line_20_1_S2_genind)
Line_20_1_S3_div <- summary(Line_20_1_S3_genind)
Line_20_1_S4_div <- summary(Line_20_1_S4_genind)

Line_20_4_F_div <- summary(Line_20_4_F_genind)
Line_20_4_S1_div <- summary(Line_20_4_S1_genind)
Line_20_4_S2_div <- summary(Line_20_4_S2_genind)
Line_20_4_S3_div <- summary(Line_20_4_S3_genind)
Line_20_4_S4_div <- summary(Line_20_4_S4_genind)

Line_21_2_F_div <- summary(Line_21_2_F_genind)
Line_21_2_S1_div <- summary(Line_21_2_S1_genind)
Line_21_2_S2_div <- summary(Line_21_2_S2_genind)
Line_21_2_S3_div <- summary(Line_21_2_S3_genind)
Line_21_2_S4_div <- summary(Line_21_2_S4_genind)
Line_21_2_S5_div <- summary(Line_21_2_S5_genind)

Line_21_6_F_div <- summary(Line_21_6_F_genind)
Line_21_6_S1_div <- summary(Line_21_6_S1_genind)
Line_21_6_S2_div <- summary(Line_21_6_S2_genind)
Line_21_6_S3_div <- summary(Line_21_6_S3_genind)
Line_21_6_S4_div <- summary(Line_21_6_S4_genind)
Line_21_6_S5_div <- summary(Line_21_6_S5_genind)

Line_23_2_F_div <- summary(Line_23_2_F_genind)
Line_23_2_S1_div <- summary(Line_23_2_S1_genind)
Line_23_2_S2_div <- summary(Line_23_2_S2_genind)
Line_23_2_S3_div <- summary(Line_23_2_S3_genind)
Line_23_2_S4_div <- summary(Line_23_2_S4_genind)
Line_23_2_S5_div <- summary(Line_23_2_S5_genind)

Line_23_4_F_div <- summary(Line_23_4_F_genind)
Line_23_4_S1_div <- summary(Line_23_4_S1_genind)
Line_23_4_S2_div <- summary(Line_23_4_S2_genind)
Line_23_4_S3_div <- summary(Line_23_4_S3_genind)
Line_23_4_S4_div <- summary(Line_23_4_S4_genind)

Line_26_1_F_div <- summary(Line_26_1_F_genind)
Line_26_1_S1_div <- summary(Line_26_1_S1_genind)
Line_26_1_S2_div <- summary(Line_26_1_S2_genind)
Line_26_1_S3_div <- summary(Line_26_1_S3_genind)
Line_26_1_S4_div <- summary(Line_26_1_S4_genind)

Line_26_4_F_div <- summary(Line_26_4_F_genind)
Line_26_4_S1_div <- summary(Line_26_4_S1_genind)
Line_26_4_S2_div <- summary(Line_26_4_S2_genind)
Line_26_4_S3_div <- summary(Line_26_4_S3_genind)
Line_26_4_S4_div <- summary(Line_26_4_S4_genind)

Line_27_2_F_div <- summary(Line_27_2_F_genind)
Line_27_2_S1_div <- summary(Line_27_2_S1_genind)
Line_27_2_S2_div <- summary(Line_27_2_S2_genind)
Line_27_2_S3_div <- summary(Line_27_2_S3_genind)
Line_27_2_S4_div <- summary(Line_27_2_S4_genind)

Line_29_2_F_div <- summary(Line_29_2_F_genind)
Line_29_2_S1_div <- summary(Line_29_2_S1_genind)
Line_29_2_S2_div <- summary(Line_29_2_S2_genind)
Line_29_2_S3_div <- summary(Line_29_2_S3_genind)
Line_29_2_S4_div <- summary(Line_29_2_S4_genind)
Line_29_2_S5_div <- summary(Line_29_2_S5_genind)

Line_29_4_F_div <- summary(Line_29_4_F_genind)
Line_29_4_S1_div <- summary(Line_29_4_S1_genind)
Line_29_4_S2_div <- summary(Line_29_4_S2_genind)
Line_29_4_S3_div <- summary(Line_29_4_S3_genind)
Line_29_4_S4_div <- summary(Line_29_4_S4_genind)
Line_29_4_S5_div <- summary(Line_29_4_S5_genind)

mean(F_lines_div$Hobs)
mean(F_lines_div$Hexp)

mean(S1_lines_div$Hobs)
mean(S1_lines_div$Hexp)

mean(S2_lines_div$Hobs)
mean(S2_lines_div$Hexp)

mean(S3_lines_div$Hobs)
mean(S3_lines_div$Hexp)

mean(S4_lines_div$Hobs)
mean(S4_lines_div$Hexp)

mean(S5_lines_div$Hobs)
mean(S5_lines_div$Hexp)

# write.table(data.frame(F_lines_div$Hobs, F_lines_div$Hexp), "F_lines_Het.txt", quote = F)
# write.table(data.frame(S1_lines_div$Hobs, S1_lines_div$Hexp), "S1_lines_Het.txt", quote = F)
# write.table(data.frame(S2_lines_div$Hobs, S2_lines_div$Hexp), "S2_lines_Het.txt", quote = F)
# write.table(data.frame(S3_lines_div$Hobs, S3_lines_div$Hexp), "S3_lines_Het.txt", quote = F)
# write.table(data.frame(S4_lines_div$Hobs, S4_lines_div$Hexp), "S4_lines_Het.txt", quote = F)
# write.table(data.frame(S5_lines_div$Hobs, S5_lines_div$Hexp), "S5_lines_Het.txt", quote = F)

# Combine generations for figure
# Line 1_1
generation_hobs_Line_1_1 <- data.frame(F_Hobs = Line_1_1_F_div$Hobs,
                                       S1_Hobs = Line_1_1_S1_div$Hobs,
                                       S2_Hobs = Line_1_1_S2_div$Hobs,
                                       S3_Hobs = Line_1_1_S3_div$Hobs,
                                       S4_Hobs = Line_1_1_S4_div$Hobs,
                                       S5_Hobs = Line_1_1_S5_div$Hobs)

colnames(generation_hobs_Line_1_1) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_1_1)
melted_generation_hobs_Line_1_1 <- melt(generation_hobs_Line_1_1)

generation_hobs_Line_1_1_summary <- data.frame(FS = mean(generation_hobs_Line_1_1$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_1_1$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_1_1$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_1_1$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_1_1$S4, na.rm = T),
                                               S5 = mean(generation_hobs_Line_1_1$S5, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_1_1_summary <- melt(generation_hobs_Line_1_1_summary)

ggplot(melted_generation_hobs_Line_1_1, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_1_1_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_1_1 <- melted_generation_hobs_Line_1_1 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_1_1_summary <- data.frame(FS = mean(expectations_halving_Line_1_1$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_1_1$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_1_1$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_1_1$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_1_1$S4, na.rm = T),
                                                    S5 = mean(expectations_halving_Line_1_1$S5, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_1_1_summary_melted <- melt(expectations_halving_Line_1_1_summary)

obs_exp_1_1 <- rbind(melted_generation_hobs_Line_1_1_summary, expectations_halving_Line_1_1_summary_melted)

obs_exp_1_1_plot <- ggplot(obs_exp_1_1, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 1_1")


wilcox.test(generation_hobs_Line_1_1$S1, expectations_halving_Line_1_1$S1, paired = T)
wilcox.test(generation_hobs_Line_1_1$S2, expectations_halving_Line_1_1$S2, paired = T)
wilcox.test(generation_hobs_Line_1_1$S3, expectations_halving_Line_1_1$S3, paired = T)
wilcox.test(generation_hobs_Line_1_1$S4, expectations_halving_Line_1_1$S4, paired = T)
wilcox.test(generation_hobs_Line_1_1$S5, expectations_halving_Line_1_1$S5, paired = T)

write.table(generation_hobs_Line_1_1, "Line_1_1_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_1_1, "Line_1_1_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 6_1
generation_hobs_Line_6_1 <- data.frame(F_Hobs = Line_6_1_F_div$Hobs,
                                       S1_Hobs = Line_6_1_S1_div$Hobs,
                                       S2_Hobs = Line_6_1_S2_div$Hobs,
                                       S3_Hobs = Line_6_1_S3_div$Hobs,
                                       S4_Hobs = Line_6_1_S4_div$Hobs)

colnames(generation_hobs_Line_6_1) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_6_1)
melted_generation_hobs_Line_6_1 <- melt(generation_hobs_Line_6_1)

generation_hobs_Line_6_1_summary <- data.frame(FS = mean(generation_hobs_Line_6_1$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_6_1$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_6_1$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_6_1$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_6_1$S4, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_6_1_summary <- melt(generation_hobs_Line_6_1_summary)

ggplot(melted_generation_hobs_Line_6_1, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_6_1_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_6_1 <- melted_generation_hobs_Line_6_1 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_6_1_summary <- data.frame(FS = mean(expectations_halving_Line_6_1$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_6_1$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_6_1$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_6_1$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_6_1$S4, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_6_1_summary_melted <- melt(expectations_halving_Line_6_1_summary)

obs_exp_6_1 <- rbind(expectations_halving_Line_6_1_summary_melted, melted_generation_hobs_Line_6_1_summary)

obs_exp_6_1_plot <- ggplot(obs_exp_6_1, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 6_1")

wilcox.test(generation_hobs_Line_6_1$S1, expectations_halving_Line_6_1$S1, paired = T)
wilcox.test(generation_hobs_Line_6_1$S2, expectations_halving_Line_6_1$S2, paired = T)
wilcox.test(generation_hobs_Line_6_1$S3, expectations_halving_Line_6_1$S3, paired = T)
wilcox.test(generation_hobs_Line_6_1$S4, expectations_halving_Line_6_1$S4, paired = T)

write.table(generation_hobs_Line_6_1, "Line_6_1_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_6_1, "Line_6_1_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 6_4
generation_hobs_Line_6_4 <- data.frame(F_Hobs = Line_6_4_F_div$Hobs,
                                       S1_Hobs = Line_6_4_S1_div$Hobs,
                                       S2_Hobs = Line_6_4_S2_div$Hobs,
                                       S3_Hobs = Line_6_4_S3_div$Hobs,
                                       S4_Hobs = Line_6_4_S4_div$Hobs,
                                       S5_Hobs = Line_6_4_S5_div$Hobs)

colnames(generation_hobs_Line_6_4) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_6_4)
melted_generation_hobs_Line_6_4 <- melt(generation_hobs_Line_6_4)

generation_hobs_Line_6_4_summary <- data.frame(FS = mean(generation_hobs_Line_6_4$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_6_4$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_6_4$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_6_4$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_6_4$S4, na.rm = T),
                                               S5 = mean(generation_hobs_Line_6_4$S5, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_6_4_summary <- melt(generation_hobs_Line_6_4_summary)

ggplot(melted_generation_hobs_Line_6_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_6_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_6_4 <- melted_generation_hobs_Line_6_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_6_4_summary <- data.frame(FS = mean(expectations_halving_Line_6_4$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_6_4$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_6_4$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_6_4$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_6_4$S4, na.rm = T),
                                                    S5 = mean(expectations_halving_Line_6_4$S5, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_6_4_summary_melted <- melt(expectations_halving_Line_6_4_summary)

obs_exp_6_4 <- rbind(expectations_halving_Line_6_4_summary_melted, melted_generation_hobs_Line_6_4_summary)

obs_exp_6_4_plot <- ggplot(obs_exp_6_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 6_4")

wilcox.test(generation_hobs_Line_6_4$S1, expectations_halving_Line_6_4$S1, paired = T)
wilcox.test(generation_hobs_Line_6_4$S2, expectations_halving_Line_6_4$S2, paired = T)
wilcox.test(generation_hobs_Line_6_4$S3, expectations_halving_Line_6_4$S3, paired = T)
wilcox.test(generation_hobs_Line_6_4$S4, expectations_halving_Line_6_4$S4, paired = T)
wilcox.test(generation_hobs_Line_6_4$S5, expectations_halving_Line_6_4$S5, paired = T)

write.table(generation_hobs_Line_6_4, "Line_6_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_6_4, "Line_6_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 7_2
generation_hobs_Line_7_2 <- data.frame(F_Hobs = Line_7_2_F_div$Hobs,
                                       S1_Hobs = Line_7_2_S1_div$Hobs,
                                       S2_Hobs = Line_7_2_S2_div$Hobs,
                                       S3_Hobs = Line_7_2_S3_div$Hobs,
                                       S4_Hobs = Line_7_2_S4_div$Hobs)

colnames(generation_hobs_Line_7_2) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_7_2)
melted_generation_hobs_Line_7_2 <- melt(generation_hobs_Line_7_2)

generation_hobs_Line_7_2_summary <- data.frame(FS = mean(generation_hobs_Line_7_2$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_7_2$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_7_2$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_7_2$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_7_2$S4, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_7_2_summary <- melt(generation_hobs_Line_7_2_summary)

ggplot(melted_generation_hobs_Line_7_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_7_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_7_2 <- melted_generation_hobs_Line_7_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_7_2_summary <- data.frame(FS = mean(expectations_halving_Line_7_2$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_7_2$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_7_2$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_7_2$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_7_2$S4, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_7_2_summary_melted <- melt(expectations_halving_Line_7_2_summary)

obs_exp_7_2 <- rbind(expectations_halving_Line_7_2_summary_melted, melted_generation_hobs_Line_7_2_summary)

obs_exp_7_2_plot <- ggplot(obs_exp_7_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 7_2")

wilcox.test(generation_hobs_Line_7_2$S1, expectations_halving_Line_7_2$S1, paired = T)
wilcox.test(generation_hobs_Line_7_2$S2, expectations_halving_Line_7_2$S2, paired = T)
wilcox.test(generation_hobs_Line_7_2$S3, expectations_halving_Line_7_2$S3, paired = T)
wilcox.test(generation_hobs_Line_7_2$S4, expectations_halving_Line_7_2$S4, paired = T)

write.table(generation_hobs_Line_7_2, "Line_7_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_7_2, "Line_7_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 7_4
generation_hobs_Line_7_4 <- data.frame(F_Hobs = Line_7_4_F_div$Hobs,
                                       S1_Hobs = Line_7_4_S1_div$Hobs,
                                       S2_Hobs = Line_7_4_S2_div$Hobs,
                                       S3_Hobs = Line_7_4_S3_div$Hobs,
                                       S4_Hobs = Line_7_4_S4_div$Hobs)

colnames(generation_hobs_Line_7_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_7_4)
melted_generation_hobs_Line_7_4 <- melt(generation_hobs_Line_7_4)

generation_hobs_Line_7_4_summary <- data.frame(FS = mean(generation_hobs_Line_7_4$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_7_4$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_7_4$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_7_4$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_7_4$S4, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_7_4_summary <- melt(generation_hobs_Line_7_4_summary)

ggplot(melted_generation_hobs_Line_7_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_7_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_7_4 <- melted_generation_hobs_Line_7_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_7_4_summary <- data.frame(FS = mean(expectations_halving_Line_7_4$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_7_4$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_7_4$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_7_4$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_7_4$S4, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_7_4_summary_melted <- melt(expectations_halving_Line_7_4_summary)

obs_exp_7_4 <- rbind(expectations_halving_Line_7_4_summary_melted, melted_generation_hobs_Line_7_4_summary)

obs_exp_7_4_plot <- ggplot(obs_exp_7_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 7_4")

wilcox.test(generation_hobs_Line_7_4$S1, expectations_halving_Line_7_4$S1, paired = T)
wilcox.test(generation_hobs_Line_7_4$S2, expectations_halving_Line_7_4$S2, paired = T)
wilcox.test(generation_hobs_Line_7_4$S3, expectations_halving_Line_7_4$S3, paired = T)
wilcox.test(generation_hobs_Line_7_4$S4, expectations_halving_Line_7_4$S4, paired = T)

write.table(generation_hobs_Line_7_4, "Line_7_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_7_4, "Line_7_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 8_2
generation_hobs_Line_8_2 <- data.frame(F_Hobs = Line_8_2_F_div$Hobs,
                                       S1_Hobs = Line_8_2_S1_div$Hobs,
                                       S2_Hobs = Line_8_2_S2_div$Hobs,
                                       S3_Hobs = Line_8_2_S3_div$Hobs,
                                       S4_Hobs = Line_8_2_S4_div$Hobs)

colnames(generation_hobs_Line_8_2) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_8_2)
melted_generation_hobs_Line_8_2 <- melt(generation_hobs_Line_8_2)

generation_hobs_Line_8_2_summary <- data.frame(FS = mean(generation_hobs_Line_8_2$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_8_2$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_8_2$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_8_2$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_8_2$S4, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_8_2_summary <- melt(generation_hobs_Line_8_2_summary)

ggplot(melted_generation_hobs_Line_8_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_8_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_8_2 <- melted_generation_hobs_Line_8_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_8_2_summary <- data.frame(FS = mean(expectations_halving_Line_8_2$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_8_2$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_8_2$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_8_2$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_8_2$S4, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_8_2_summary_melted <- melt(expectations_halving_Line_8_2_summary)

obs_exp_8_2 <- rbind(expectations_halving_Line_8_2_summary_melted, melted_generation_hobs_Line_8_2_summary)

obs_exp_8_2_plot <- ggplot(obs_exp_8_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 8_2")

wilcox.test(generation_hobs_Line_8_2$S1, expectations_halving_Line_8_2$S1, paired = T)
wilcox.test(generation_hobs_Line_8_2$S2, expectations_halving_Line_8_2$S2, paired = T)
wilcox.test(generation_hobs_Line_8_2$S3, expectations_halving_Line_8_2$S3, paired = T)
wilcox.test(generation_hobs_Line_8_2$S4, expectations_halving_Line_8_2$S4, paired = T)

write.table(generation_hobs_Line_8_2, "Line_8_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_8_2, "Line_8_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 8_4
generation_hobs_Line_8_4 <- data.frame(F_Hobs = Line_8_4_F_div$Hobs,
                                       S1_Hobs = Line_8_4_S1_div$Hobs,
                                       S2_Hobs = Line_8_4_S2_div$Hobs,
                                       S3_Hobs = Line_8_4_S3_div$Hobs,
                                       S4_Hobs = Line_8_4_S4_div$Hobs)

colnames(generation_hobs_Line_8_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_8_4)
melted_generation_hobs_Line_8_4 <- melt(generation_hobs_Line_8_4)

generation_hobs_Line_8_4_summary <- data.frame(FS = mean(generation_hobs_Line_8_4$F, na.rm = T), 
                                               S1 = mean(generation_hobs_Line_8_4$S1, na.rm = T),
                                               S2 = mean(generation_hobs_Line_8_4$S2, na.rm = T),
                                               S3 = mean(generation_hobs_Line_8_4$S3, na.rm = T),
                                               S4 = mean(generation_hobs_Line_8_4$S4, na.rm = T),
                                               H = "Observed")

melted_generation_hobs_Line_8_4_summary <- melt(generation_hobs_Line_8_4_summary)

ggplot(melted_generation_hobs_Line_8_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_8_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_8_4 <- melted_generation_hobs_Line_8_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_8_4_summary <- data.frame(FS = mean(expectations_halving_Line_8_4$F, na.rm = T), 
                                                    S1 = mean(expectations_halving_Line_8_4$S1, na.rm = T),
                                                    S2 = mean(expectations_halving_Line_8_4$S2, na.rm = T),
                                                    S3 = mean(expectations_halving_Line_8_4$S3, na.rm = T),
                                                    S4 = mean(expectations_halving_Line_8_4$S4, na.rm = T),
                                                    H = "Expected")
expectations_halving_Line_8_4_summary_melted <- melt(expectations_halving_Line_8_4_summary)

obs_exp_8_4 <- rbind(expectations_halving_Line_8_4_summary_melted, melted_generation_hobs_Line_8_4_summary)

obs_exp_8_4_plot <- ggplot(obs_exp_8_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 8_4")

wilcox.test(generation_hobs_Line_8_4$S1, expectations_halving_Line_8_4$S1, paired = T)
wilcox.test(generation_hobs_Line_8_4$S2, expectations_halving_Line_8_4$S2, paired = T)
wilcox.test(generation_hobs_Line_8_4$S3, expectations_halving_Line_8_4$S3, paired = T)
wilcox.test(generation_hobs_Line_8_4$S4, expectations_halving_Line_8_4$S4, paired = T)

write.table(generation_hobs_Line_8_4, "Line_8_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_8_4, "Line_8_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 10_2
generation_hobs_Line_10_2 <- data.frame(F_Hobs = Line_10_2_F_div$Hobs,
                                        S1_Hobs = Line_10_2_S1_div$Hobs,
                                        S2_Hobs = Line_10_2_S2_div$Hobs,
                                        S3_Hobs = Line_10_2_S3_div$Hobs,
                                        S4_Hobs = Line_10_2_S4_div$Hobs,
                                        S5_Hobs = Line_10_2_S5_div$Hobs)

colnames(generation_hobs_Line_10_2) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_10_2)
melted_generation_hobs_Line_10_2 <- melt(generation_hobs_Line_10_2)

generation_hobs_Line_10_2_summary <- data.frame(FS = mean(generation_hobs_Line_10_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_10_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_10_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_10_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_10_2$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_10_2$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_10_2_summary <- melt(generation_hobs_Line_10_2_summary)

ggplot(melted_generation_hobs_Line_10_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_10_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_10_2 <- melted_generation_hobs_Line_10_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_10_2_summary <- data.frame(FS = mean(expectations_halving_Line_10_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_10_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_10_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_10_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_10_2$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_10_2$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_10_2_summary_melted <- melt(expectations_halving_Line_10_2_summary)

obs_exp_10_2 <- rbind(expectations_halving_Line_10_2_summary_melted, melted_generation_hobs_Line_10_2_summary)

obs_exp_10_2_plot <- ggplot(obs_exp_10_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 10_2")

wilcox.test(generation_hobs_Line_10_2$S1, expectations_halving_Line_10_2$S1, paired = T)
wilcox.test(generation_hobs_Line_10_2$S2, expectations_halving_Line_10_2$S2, paired = T)
wilcox.test(generation_hobs_Line_10_2$S3, expectations_halving_Line_10_2$S3, paired = T)
wilcox.test(generation_hobs_Line_10_2$S4, expectations_halving_Line_10_2$S4, paired = T)
wilcox.test(generation_hobs_Line_10_2$S5, expectations_halving_Line_10_2$S5, paired = T)

write.table(generation_hobs_Line_10_2, "Line_10_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_10_2, "Line_10_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 10_5
generation_hobs_Line_10_5 <- data.frame(F_Hobs = Line_10_5_F_div$Hobs,
                                        S1_Hobs = Line_10_5_S1_div$Hobs,
                                        S2_Hobs = Line_10_5_S2_div$Hobs,
                                        S3_Hobs = Line_10_5_S3_div$Hobs,
                                        S4_Hobs = Line_10_5_S4_div$Hobs)

colnames(generation_hobs_Line_10_5) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_10_5)
melted_generation_hobs_Line_10_5 <- melt(generation_hobs_Line_10_5)

generation_hobs_Line_10_5_summary <- data.frame(FS = mean(generation_hobs_Line_10_5$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_10_5$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_10_5$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_10_5$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_10_5$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_10_5_summary <- melt(generation_hobs_Line_10_5_summary)

ggplot(melted_generation_hobs_Line_10_5, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_10_5_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_10_5 <- melted_generation_hobs_Line_10_5 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_10_5_summary <- data.frame(FS = mean(expectations_halving_Line_10_5$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_10_5$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_10_5$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_10_5$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_10_5$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_10_5_summary_melted <- melt(expectations_halving_Line_10_5_summary)

obs_exp_10_5 <- rbind(expectations_halving_Line_10_5_summary_melted, melted_generation_hobs_Line_10_5_summary)

obs_exp_10_5_plot <- ggplot(obs_exp_10_5, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 10_5")

wilcox.test(generation_hobs_Line_10_5$S1, expectations_halving_Line_10_5$S1, paired = T)
wilcox.test(generation_hobs_Line_10_5$S2, expectations_halving_Line_10_5$S2, paired = T)
wilcox.test(generation_hobs_Line_10_5$S3, expectations_halving_Line_10_5$S3, paired = T)
wilcox.test(generation_hobs_Line_10_5$S4, expectations_halving_Line_10_5$S4, paired = T)

write.table(generation_hobs_Line_10_5, "Line_10_5_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_10_5, "Line_10_5_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 13_3
generation_hobs_Line_13_3 <- data.frame(F_Hobs = Line_13_3_F_div$Hobs,
                                        S1_Hobs = Line_13_3_S1_div$Hobs,
                                        S2_Hobs = Line_13_3_S2_div$Hobs,
                                        S3_Hobs = Line_13_3_S3_div$Hobs,
                                        S4_Hobs = Line_13_3_S4_div$Hobs)

colnames(generation_hobs_Line_13_3) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_13_3)
melted_generation_hobs_Line_13_3 <- melt(generation_hobs_Line_13_3)

generation_hobs_Line_13_3_summary <- data.frame(FS = mean(generation_hobs_Line_13_3$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_13_3$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_13_3$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_13_3$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_13_3$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_13_3_summary <- melt(generation_hobs_Line_13_3_summary)

ggplot(melted_generation_hobs_Line_13_3, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_13_3_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_13_3 <- melted_generation_hobs_Line_13_3 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_13_3_summary <- data.frame(FS = mean(expectations_halving_Line_13_3$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_13_3$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_13_3$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_13_3$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_13_3$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_13_3_summary_melted <- melt(expectations_halving_Line_13_3_summary)

obs_exp_13_3 <- rbind(expectations_halving_Line_13_3_summary_melted, melted_generation_hobs_Line_13_3_summary)

obs_exp_13_3_plot <- ggplot(obs_exp_13_3, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 13_3")

wilcox.test(generation_hobs_Line_13_3$S1, expectations_halving_Line_13_3$S1, paired = T)
wilcox.test(generation_hobs_Line_13_3$S2, expectations_halving_Line_13_3$S2, paired = T)
wilcox.test(generation_hobs_Line_13_3$S3, expectations_halving_Line_13_3$S3, paired = T)
wilcox.test(generation_hobs_Line_13_3$S4, expectations_halving_Line_13_3$S4, paired = T)

write.table(generation_hobs_Line_13_3, "Line_13_3_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_13_3, "Line_13_3_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 13_4
generation_hobs_Line_13_4 <- data.frame(F_Hobs = Line_13_4_F_div$Hobs,
                                        S1_Hobs = Line_13_4_S1_div$Hobs,
                                        S2_Hobs = Line_13_4_S2_div$Hobs,
                                        S3_Hobs = Line_13_4_S3_div$Hobs,
                                        S4_Hobs = Line_13_4_S4_div$Hobs)

colnames(generation_hobs_Line_13_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_13_4)
melted_generation_hobs_Line_13_4 <- melt(generation_hobs_Line_13_4)

generation_hobs_Line_13_4_summary <- data.frame(FS = mean(generation_hobs_Line_13_4$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_13_4$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_13_4$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_13_4$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_13_4$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_13_4_summary <- melt(generation_hobs_Line_13_4_summary)

ggplot(melted_generation_hobs_Line_13_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_13_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_13_4 <- melted_generation_hobs_Line_13_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_13_4_summary <- data.frame(FS = mean(expectations_halving_Line_13_4$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_13_4$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_13_4$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_13_4$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_13_4$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_13_4_summary_melted <- melt(expectations_halving_Line_13_4_summary)

obs_exp_13_4 <- rbind(expectations_halving_Line_13_4_summary_melted, melted_generation_hobs_Line_13_4_summary)

obs_exp_13_4_plot <- ggplot(obs_exp_13_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 13_4")

wilcox.test(generation_hobs_Line_13_4$S1, expectations_halving_Line_13_4$S1, paired = T)
wilcox.test(generation_hobs_Line_13_4$S2, expectations_halving_Line_13_4$S2, paired = T)
wilcox.test(generation_hobs_Line_13_4$S3, expectations_halving_Line_13_4$S3, paired = T)
wilcox.test(generation_hobs_Line_13_4$S4, expectations_halving_Line_13_4$S4, paired = T)

write.table(generation_hobs_Line_13_4, "Line_13_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_13_4, "Line_13_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 16_1
generation_hobs_Line_16_1 <- data.frame(F_Hobs = Line_16_1_F_div$Hobs,
                                        S1_Hobs = Line_16_1_S1_div$Hobs,
                                        S2_Hobs = Line_16_1_S2_div$Hobs,
                                        S3_Hobs = Line_16_1_S3_div$Hobs,
                                        S4_Hobs = Line_16_1_S4_div$Hobs)

colnames(generation_hobs_Line_16_1) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_16_1)
melted_generation_hobs_Line_16_1 <- melt(generation_hobs_Line_16_1)

generation_hobs_Line_16_1_summary <- data.frame(FS = mean(generation_hobs_Line_16_1$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_16_1$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_16_1$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_16_1$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_16_1$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_16_1_summary <- melt(generation_hobs_Line_16_1_summary)

ggplot(melted_generation_hobs_Line_16_1, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_16_1_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_16_1 <- melted_generation_hobs_Line_16_1 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_16_1_summary <- data.frame(FS = mean(expectations_halving_Line_16_1$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_16_1$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_16_1$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_16_1$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_16_1$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_16_1_summary_melted <- melt(expectations_halving_Line_16_1_summary)

obs_exp_16_1 <- rbind(expectations_halving_Line_16_1_summary_melted, melted_generation_hobs_Line_16_1_summary)

obs_exp_16_1_plot <- ggplot(obs_exp_16_1, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 16_1")

wilcox.test(generation_hobs_Line_16_1$S1, expectations_halving_Line_16_1$S1, paired = T)
wilcox.test(generation_hobs_Line_16_1$S2, expectations_halving_Line_16_1$S2, paired = T)
wilcox.test(generation_hobs_Line_16_1$S3, expectations_halving_Line_16_1$S3, paired = T)
wilcox.test(generation_hobs_Line_16_1$S4, expectations_halving_Line_16_1$S4, paired = T)

write.table(generation_hobs_Line_16_1, "Line_16_1_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_16_1, "Line_16_1_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 16_5
generation_hobs_Line_16_5 <- data.frame(F_Hobs = Line_16_5_F_div$Hobs,
                                        S1_Hobs = Line_16_5_S1_div$Hobs,
                                        S2_Hobs = Line_16_5_S2_div$Hobs,
                                        S3_Hobs = Line_16_5_S3_div$Hobs,
                                        S4_Hobs = Line_16_5_S4_div$Hobs,
                                        S5_Hobs = Line_16_5_S5_div$Hobs)

colnames(generation_hobs_Line_16_5) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_16_5)
melted_generation_hobs_Line_16_5 <- melt(generation_hobs_Line_16_5)

generation_hobs_Line_16_5_summary <- data.frame(FS = mean(generation_hobs_Line_16_5$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_16_5$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_16_5$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_16_5$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_16_5$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_16_5$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_16_5_summary <- melt(generation_hobs_Line_16_5_summary)

ggplot(melted_generation_hobs_Line_16_5, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_16_5_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_16_5 <- melted_generation_hobs_Line_16_5 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_16_5_summary <- data.frame(FS = mean(expectations_halving_Line_16_5$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_16_5$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_16_5$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_16_5$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_16_5$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_16_5$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_16_5_summary_melted <- melt(expectations_halving_Line_16_5_summary)

obs_exp_16_5 <- rbind(expectations_halving_Line_16_5_summary_melted, melted_generation_hobs_Line_16_5_summary)

obs_exp_16_5_plot <- ggplot(obs_exp_16_5, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 16_5")

wilcox.test(generation_hobs_Line_16_5$S1, expectations_halving_Line_16_5$S1, paired = T)
wilcox.test(generation_hobs_Line_16_5$S2, expectations_halving_Line_16_5$S2, paired = T)
wilcox.test(generation_hobs_Line_16_5$S3, expectations_halving_Line_16_5$S3, paired = T)
wilcox.test(generation_hobs_Line_16_5$S4, expectations_halving_Line_16_5$S4, paired = T)
wilcox.test(generation_hobs_Line_16_5$S5, expectations_halving_Line_16_5$S5, paired = T)

write.table(generation_hobs_Line_16_5, "Line_16_5_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_16_5, "Line_16_5_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 17_2
generation_hobs_Line_17_2 <- data.frame(F_Hobs = Line_17_2_F_div$Hobs,
                                        S1_Hobs = Line_17_2_S1_div$Hobs,
                                        S2_Hobs = Line_17_2_S2_div$Hobs,
                                        S3_Hobs = Line_17_2_S3_div$Hobs,
                                        S4_Hobs = Line_17_2_S4_div$Hobs,
                                        S5_Hobs = Line_17_2_S5_div$Hobs)

colnames(generation_hobs_Line_17_2) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_17_2)
melted_generation_hobs_Line_17_2 <- melt(generation_hobs_Line_17_2)

generation_hobs_Line_17_2_summary <- data.frame(FS = mean(generation_hobs_Line_17_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_17_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_17_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_17_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_17_2$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_17_2$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_17_2_summary <- melt(generation_hobs_Line_17_2_summary)

ggplot(melted_generation_hobs_Line_17_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_17_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_17_2 <- melted_generation_hobs_Line_17_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_17_2_summary <- data.frame(FS = mean(expectations_halving_Line_17_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_17_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_17_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_17_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_17_2$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_17_2$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_17_2_summary_melted <- melt(expectations_halving_Line_17_2_summary)

obs_exp_17_2 <- rbind(expectations_halving_Line_17_2_summary_melted, melted_generation_hobs_Line_17_2_summary)

obs_exp_17_2_plot <- ggplot(obs_exp_17_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 17_2")

wilcox.test(generation_hobs_Line_17_2$S1, expectations_halving_Line_17_2$S1, paired = T)
wilcox.test(generation_hobs_Line_17_2$S2, expectations_halving_Line_17_2$S2, paired = T)
wilcox.test(generation_hobs_Line_17_2$S3, expectations_halving_Line_17_2$S3, paired = T)
wilcox.test(generation_hobs_Line_17_2$S4, expectations_halving_Line_17_2$S4, paired = T)
wilcox.test(generation_hobs_Line_17_2$S5, expectations_halving_Line_17_2$S5, paired = T)

write.table(generation_hobs_Line_17_2, "Line_17_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_17_2, "Line_17_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 17_5
generation_hobs_Line_17_5 <- data.frame(F_Hobs = Line_17_5_F_div$Hobs,
                                        S1_Hobs = Line_17_5_S1_div$Hobs,
                                        S2_Hobs = Line_17_5_S2_div$Hobs,
                                        S3_Hobs = Line_17_5_S3_div$Hobs,
                                        S4_Hobs = Line_17_5_S4_div$Hobs,
                                        S5_Hobs = Line_17_5_S5_div$Hobs)

colnames(generation_hobs_Line_17_5) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_17_5)
melted_generation_hobs_Line_17_5 <- melt(generation_hobs_Line_17_5)

generation_hobs_Line_17_5_summary <- data.frame(FS = mean(generation_hobs_Line_17_5$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_17_5$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_17_5$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_17_5$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_17_5$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_17_5$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_17_5_summary <- melt(generation_hobs_Line_17_5_summary)

ggplot(melted_generation_hobs_Line_17_5, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_17_5_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_17_5 <- melted_generation_hobs_Line_17_5 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_17_5_summary <- data.frame(FS = mean(expectations_halving_Line_17_5$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_17_5$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_17_5$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_17_5$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_17_5$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_17_5$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_17_5_summary_melted <- melt(expectations_halving_Line_17_5_summary)

obs_exp_17_5 <- rbind(expectations_halving_Line_17_5_summary_melted, melted_generation_hobs_Line_17_5_summary)

obs_exp_17_5_plot <- ggplot(obs_exp_17_5, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 17_5")

wilcox.test(generation_hobs_Line_17_5$S1, expectations_halving_Line_17_5$S1, paired = T)
wilcox.test(generation_hobs_Line_17_5$S2, expectations_halving_Line_17_5$S2, paired = T)
wilcox.test(generation_hobs_Line_17_5$S3, expectations_halving_Line_17_5$S3, paired = T)
wilcox.test(generation_hobs_Line_17_5$S4, expectations_halving_Line_17_5$S4, paired = T)
wilcox.test(generation_hobs_Line_17_5$S5, expectations_halving_Line_17_5$S5, paired = T)

write.table(generation_hobs_Line_17_5, "Line_17_5_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_17_5, "Line_17_5_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 19_2
generation_hobs_Line_19_2 <- data.frame(F_Hobs = Line_19_2_F_div$Hobs,
                                        S1_Hobs = Line_19_2_S1_div$Hobs,
                                        S2_Hobs = Line_19_2_S2_div$Hobs,
                                        S3_Hobs = Line_19_2_S3_div$Hobs,
                                        S4_Hobs = Line_19_2_S4_div$Hobs)

colnames(generation_hobs_Line_19_2) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_19_2)
melted_generation_hobs_Line_19_2 <- melt(generation_hobs_Line_19_2)

generation_hobs_Line_19_2_summary <- data.frame(FS = mean(generation_hobs_Line_19_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_19_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_19_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_19_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_19_2$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_19_2_summary <- melt(generation_hobs_Line_19_2_summary)

ggplot(melted_generation_hobs_Line_19_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_19_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_19_2 <- melted_generation_hobs_Line_19_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_19_2_summary <- data.frame(FS = mean(expectations_halving_Line_19_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_19_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_19_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_19_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_19_2$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_19_2_summary_melted <- melt(expectations_halving_Line_19_2_summary)

obs_exp_19_2 <- rbind(expectations_halving_Line_19_2_summary_melted, melted_generation_hobs_Line_19_2_summary)

obs_exp_19_2_plot <- ggplot(obs_exp_19_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 19_2")

wilcox.test(generation_hobs_Line_19_2$S1, expectations_halving_Line_19_2$S1, paired = T)
wilcox.test(generation_hobs_Line_19_2$S2, expectations_halving_Line_19_2$S2, paired = T)
wilcox.test(generation_hobs_Line_19_2$S3, expectations_halving_Line_19_2$S3, paired = T)
wilcox.test(generation_hobs_Line_19_2$S4, expectations_halving_Line_19_2$S4, paired = T)

write.table(generation_hobs_Line_19_2, "Line_19_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_19_2, "Line_19_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 19_5
generation_hobs_Line_19_5 <- data.frame(F_Hobs = Line_19_5_F_div$Hobs,
                                        S1_Hobs = Line_19_5_S1_div$Hobs,
                                        S2_Hobs = Line_19_5_S2_div$Hobs,
                                        S3_Hobs = Line_19_5_S3_div$Hobs,
                                        S4_Hobs = Line_19_5_S4_div$Hobs)

colnames(generation_hobs_Line_19_5) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_19_5)
melted_generation_hobs_Line_19_5 <- melt(generation_hobs_Line_19_5)

generation_hobs_Line_19_5_summary <- data.frame(FS = mean(generation_hobs_Line_19_5$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_19_5$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_19_5$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_19_5$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_19_5$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_19_5_summary <- melt(generation_hobs_Line_19_5_summary)

ggplot(melted_generation_hobs_Line_19_5, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_19_5_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_19_5 <- melted_generation_hobs_Line_19_5 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_19_5_summary <- data.frame(FS = mean(expectations_halving_Line_19_5$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_19_5$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_19_5$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_19_5$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_19_5$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_19_5_summary_melted <- melt(expectations_halving_Line_19_5_summary)

obs_exp_19_5 <- rbind(expectations_halving_Line_19_5_summary_melted, melted_generation_hobs_Line_19_5_summary)

obs_exp_19_5_plot <- ggplot(obs_exp_19_5, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 19_5")

wilcox.test(generation_hobs_Line_19_5$S1, expectations_halving_Line_19_5$S1, paired = T)
wilcox.test(generation_hobs_Line_19_5$S2, expectations_halving_Line_19_5$S2, paired = T)
wilcox.test(generation_hobs_Line_19_5$S3, expectations_halving_Line_19_5$S3, paired = T)
wilcox.test(generation_hobs_Line_19_5$S4, expectations_halving_Line_19_5$S4, paired = T)

write.table(generation_hobs_Line_19_5, "Line_19_5_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_19_5, "Line_19_5_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 20_1
generation_hobs_Line_20_1 <- data.frame(F_Hobs = Line_20_1_F_div$Hobs,
                                        S1_Hobs = Line_20_1_S1_div$Hobs,
                                        S2_Hobs = Line_20_1_S2_div$Hobs,
                                        S3_Hobs = Line_20_1_S3_div$Hobs,
                                        S4_Hobs = Line_20_1_S4_div$Hobs)

colnames(generation_hobs_Line_20_1) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_20_1)
melted_generation_hobs_Line_20_1 <- melt(generation_hobs_Line_20_1)

generation_hobs_Line_20_1_summary <- data.frame(FS = mean(generation_hobs_Line_20_1$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_20_1$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_20_1$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_20_1$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_20_1$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_20_1_summary <- melt(generation_hobs_Line_20_1_summary)

ggplot(melted_generation_hobs_Line_20_1, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_20_1_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_20_1 <- melted_generation_hobs_Line_20_1 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_20_1_summary <- data.frame(FS = mean(expectations_halving_Line_20_1$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_20_1$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_20_1$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_20_1$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_20_1$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_20_1_summary_melted <- melt(expectations_halving_Line_20_1_summary)

obs_exp_20_1 <- rbind(expectations_halving_Line_20_1_summary_melted, melted_generation_hobs_Line_20_1_summary)

obs_exp_20_1_plot <- ggplot(obs_exp_20_1, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 20_1")

wilcox.test(generation_hobs_Line_20_1$S1, expectations_halving_Line_20_1$S1, paired = T)
wilcox.test(generation_hobs_Line_20_1$S2, expectations_halving_Line_20_1$S2, paired = T)
wilcox.test(generation_hobs_Line_20_1$S3, expectations_halving_Line_20_1$S3, paired = T)
wilcox.test(generation_hobs_Line_20_1$S4, expectations_halving_Line_20_1$S4, paired = T)

write.table(generation_hobs_Line_20_1, "Line_20_1_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_20_1, "Line_20_1_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 20_4
generation_hobs_Line_20_4 <- data.frame(F_Hobs = Line_20_4_F_div$Hobs,
                                        S1_Hobs = Line_20_4_S1_div$Hobs,
                                        S2_Hobs = Line_20_4_S2_div$Hobs,
                                        S3_Hobs = Line_20_4_S3_div$Hobs,
                                        S4_Hobs = Line_20_4_S4_div$Hobs)

colnames(generation_hobs_Line_20_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_20_4)
melted_generation_hobs_Line_20_4 <- melt(generation_hobs_Line_20_4)

generation_hobs_Line_20_4_summary <- data.frame(FS = mean(generation_hobs_Line_20_4$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_20_4$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_20_4$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_20_4$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_20_4$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_20_4_summary <- melt(generation_hobs_Line_20_4_summary)

ggplot(melted_generation_hobs_Line_20_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_20_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_20_4 <- melted_generation_hobs_Line_20_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_20_4_summary <- data.frame(FS = mean(expectations_halving_Line_20_4$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_20_4$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_20_4$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_20_4$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_20_4$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_20_4_summary_melted <- melt(expectations_halving_Line_20_4_summary)

obs_exp_20_4 <- rbind(expectations_halving_Line_20_4_summary_melted, melted_generation_hobs_Line_20_4_summary)

obs_exp_20_4_plot <- ggplot(obs_exp_20_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 20_4")

wilcox.test(generation_hobs_Line_20_4$S1, expectations_halving_Line_20_4$S1, paired = T)
wilcox.test(generation_hobs_Line_20_4$S2, expectations_halving_Line_20_4$S2, paired = T)
wilcox.test(generation_hobs_Line_20_4$S3, expectations_halving_Line_20_4$S3, paired = T)
wilcox.test(generation_hobs_Line_20_4$S4, expectations_halving_Line_20_4$S4, paired = T)

write.table(generation_hobs_Line_20_4, "Line_20_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_20_4, "Line_20_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 21_2
generation_hobs_Line_21_2 <- data.frame(F_Hobs = Line_21_2_F_div$Hobs,
                                        S1_Hobs = Line_21_2_S1_div$Hobs,
                                        S2_Hobs = Line_21_2_S2_div$Hobs,
                                        S3_Hobs = Line_21_2_S3_div$Hobs,
                                        S4_Hobs = Line_21_2_S4_div$Hobs,
                                        S5_Hobs = Line_21_2_S5_div$Hobs)

colnames(generation_hobs_Line_21_2) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_21_2)
melted_generation_hobs_Line_21_2 <- melt(generation_hobs_Line_21_2)

generation_hobs_Line_21_2_summary <- data.frame(FS = mean(generation_hobs_Line_21_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_21_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_21_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_21_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_21_2$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_21_2$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_21_2_summary <- melt(generation_hobs_Line_21_2_summary)

ggplot(melted_generation_hobs_Line_21_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_21_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_21_2 <- melted_generation_hobs_Line_21_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_21_2_summary <- data.frame(FS = mean(expectations_halving_Line_21_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_21_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_21_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_21_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_21_2$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_21_2$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_21_2_summary_melted <- melt(expectations_halving_Line_21_2_summary)

obs_exp_21_2 <- rbind(expectations_halving_Line_21_2_summary_melted, melted_generation_hobs_Line_21_2_summary)

obs_exp_21_2_plot <- ggplot(obs_exp_21_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 21_2")

wilcox.test(generation_hobs_Line_21_2$S1, expectations_halving_Line_21_2$S1, paired = T)
wilcox.test(generation_hobs_Line_21_2$S2, expectations_halving_Line_21_2$S2, paired = T)
wilcox.test(generation_hobs_Line_21_2$S3, expectations_halving_Line_21_2$S3, paired = T)
wilcox.test(generation_hobs_Line_21_2$S4, expectations_halving_Line_21_2$S4, paired = T)
wilcox.test(generation_hobs_Line_21_2$S5, expectations_halving_Line_21_2$S5, paired = T)

write.table(generation_hobs_Line_21_2, "Line_21_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_21_2, "Line_21_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 21_6
generation_hobs_Line_21_6 <- data.frame(F_Hobs = Line_21_6_F_div$Hobs,
                                        S1_Hobs = Line_21_6_S1_div$Hobs,
                                        S2_Hobs = Line_21_6_S2_div$Hobs,
                                        S3_Hobs = Line_21_6_S3_div$Hobs,
                                        S4_Hobs = Line_21_6_S4_div$Hobs,
                                        S5_Hobs = Line_21_6_S5_div$Hobs)

colnames(generation_hobs_Line_21_6) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_21_6)
melted_generation_hobs_Line_21_6 <- melt(generation_hobs_Line_21_6)

generation_hobs_Line_21_6_summary <- data.frame(FS = mean(generation_hobs_Line_21_6$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_21_6$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_21_6$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_21_6$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_21_6$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_21_6$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_21_6_summary <- melt(generation_hobs_Line_21_6_summary)

ggplot(melted_generation_hobs_Line_21_6, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_21_6_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_21_6 <- melted_generation_hobs_Line_21_6 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_21_6_summary <- data.frame(FS = mean(expectations_halving_Line_21_6$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_21_6$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_21_6$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_21_6$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_21_6$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_21_6$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_21_6_summary_melted <- melt(expectations_halving_Line_21_6_summary)

obs_exp_21_6 <- rbind(expectations_halving_Line_21_6_summary_melted, melted_generation_hobs_Line_21_6_summary)

obs_exp_21_6_plot <- ggplot(obs_exp_21_6, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 21_6")

wilcox.test(generation_hobs_Line_21_6$S1, expectations_halving_Line_21_6$S1, paired = T)
wilcox.test(generation_hobs_Line_21_6$S2, expectations_halving_Line_21_6$S2, paired = T)
wilcox.test(generation_hobs_Line_21_6$S3, expectations_halving_Line_21_6$S3, paired = T)
wilcox.test(generation_hobs_Line_21_6$S4, expectations_halving_Line_21_6$S4, paired = T)
wilcox.test(generation_hobs_Line_21_6$S5, expectations_halving_Line_21_6$S5, paired = T)

write.table(generation_hobs_Line_21_6, "Line_21_6_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_21_6, "Line_21_6_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 23_2
generation_hobs_Line_23_2 <- data.frame(F_Hobs = Line_23_2_F_div$Hobs,
                                        S1_Hobs = Line_23_2_S1_div$Hobs,
                                        S2_Hobs = Line_23_2_S2_div$Hobs,
                                        S3_Hobs = Line_23_2_S3_div$Hobs,
                                        S4_Hobs = Line_23_2_S4_div$Hobs,
                                        S5_Hobs = Line_23_2_S5_div$Hobs)

colnames(generation_hobs_Line_23_2) <- c("F", "S1", "S2", "S3", "S4", "S5")
# generation_hobs <- rownames_to_column(generation_hobs_Line_23_2)
melted_generation_hobs_Line_23_2 <- melt(generation_hobs_Line_23_2)

generation_hobs_Line_23_2_summary <- data.frame(FS = mean(generation_hobs_Line_23_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_23_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_23_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_23_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_23_2$S4, na.rm = T),
                                                S5 = mean(generation_hobs_Line_23_2$S5, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_23_2_summary <- melt(generation_hobs_Line_23_2_summary)

ggplot(melted_generation_hobs_Line_23_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_23_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_23_2 <- melted_generation_hobs_Line_23_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4, S5)

expectations_halving_Line_23_2_summary <- data.frame(FS = mean(expectations_halving_Line_23_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_23_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_23_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_23_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_23_2$S4, na.rm = T),
                                                     S5 = mean(expectations_halving_Line_23_2$S5, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_23_2_summary_melted <- melt(expectations_halving_Line_23_2_summary)

obs_exp_23_2 <- rbind(expectations_halving_Line_23_2_summary_melted, melted_generation_hobs_Line_23_2_summary)

obs_exp_23_2_plot <- ggplot(obs_exp_23_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 23_2")

wilcox.test(generation_hobs_Line_23_2$S1, expectations_halving_Line_23_2$S1, paired = T)
wilcox.test(generation_hobs_Line_23_2$S2, expectations_halving_Line_23_2$S2, paired = T)
wilcox.test(generation_hobs_Line_23_2$S3, expectations_halving_Line_23_2$S3, paired = T)
wilcox.test(generation_hobs_Line_23_2$S4, expectations_halving_Line_23_2$S4, paired = T)
wilcox.test(generation_hobs_Line_23_2$S5, expectations_halving_Line_23_2$S5, paired = T)

write.table(generation_hobs_Line_23_2, "Line_23_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_23_2, "Line_23_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 23_4
generation_hobs_Line_23_4 <- data.frame(F_Hobs = Line_23_4_F_div$Hobs,
                                        S1_Hobs = Line_23_4_S1_div$Hobs,
                                        S2_Hobs = Line_23_4_S2_div$Hobs,
                                        S3_Hobs = Line_23_4_S3_div$Hobs,
                                        S4_Hobs = Line_23_4_S4_div$Hobs)

colnames(generation_hobs_Line_23_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_23_4)
melted_generation_hobs_Line_23_4 <- melt(generation_hobs_Line_23_4)

generation_hobs_Line_23_4_summary <- data.frame(FS = mean(generation_hobs_Line_23_4$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_23_4$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_23_4$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_23_4$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_23_4$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_23_4_summary <- melt(generation_hobs_Line_23_4_summary)

ggplot(melted_generation_hobs_Line_23_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_23_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_23_4 <- melted_generation_hobs_Line_23_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_23_4_summary <- data.frame(FS = mean(expectations_halving_Line_23_4$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_23_4$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_23_4$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_23_4$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_23_4$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_23_4_summary_melted <- melt(expectations_halving_Line_23_4_summary)

obs_exp_23_4 <- rbind(expectations_halving_Line_23_4_summary_melted, melted_generation_hobs_Line_23_4_summary)

obs_exp_23_4_plot <- ggplot(obs_exp_23_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 23_4")

wilcox.test(generation_hobs_Line_23_4$S1, expectations_halving_Line_23_4$S1, paired = T)
wilcox.test(generation_hobs_Line_23_4$S2, expectations_halving_Line_23_4$S2, paired = T)
wilcox.test(generation_hobs_Line_23_4$S3, expectations_halving_Line_23_4$S3, paired = T)
wilcox.test(generation_hobs_Line_23_4$S4, expectations_halving_Line_23_4$S4, paired = T)

write.table(generation_hobs_Line_23_4, "Line_23_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_23_4, "Line_23_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 26_1
generation_hobs_Line_26_1 <- data.frame(F_Hobs = Line_26_1_F_div$Hobs,
                                        S1_Hobs = Line_26_1_S1_div$Hobs,
                                        S2_Hobs = Line_26_1_S2_div$Hobs,
                                        S3_Hobs = Line_26_1_S3_div$Hobs,
                                        S4_Hobs = Line_26_1_S4_div$Hobs)

colnames(generation_hobs_Line_26_1) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_26_1)
melted_generation_hobs_Line_26_1 <- melt(generation_hobs_Line_26_1)

generation_hobs_Line_26_1_summary <- data.frame(FS = mean(generation_hobs_Line_26_1$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_26_1$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_26_1$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_26_1$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_26_1$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_26_1_summary <- melt(generation_hobs_Line_26_1_summary)

ggplot(melted_generation_hobs_Line_26_1, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_26_1_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_26_1 <- melted_generation_hobs_Line_26_1 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_26_1_summary <- data.frame(FS = mean(expectations_halving_Line_26_1$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_26_1$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_26_1$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_26_1$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_26_1$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_26_1_summary_melted <- melt(expectations_halving_Line_26_1_summary)

obs_exp_26_1 <- rbind(expectations_halving_Line_26_1_summary_melted, melted_generation_hobs_Line_26_1_summary)

obs_exp_26_1_plot <- ggplot(obs_exp_26_1, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 26_1")

wilcox.test(generation_hobs_Line_26_1$S1, expectations_halving_Line_26_1$S1, paired = T)
wilcox.test(generation_hobs_Line_26_1$S2, expectations_halving_Line_26_1$S2, paired = T)
wilcox.test(generation_hobs_Line_26_1$S3, expectations_halving_Line_26_1$S3, paired = T)
wilcox.test(generation_hobs_Line_26_1$S4, expectations_halving_Line_26_1$S4, paired = T)

write.table(generation_hobs_Line_26_1, "Line_26_1_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_26_1, "Line_26_1_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 26_4
generation_hobs_Line_26_4 <- data.frame(F_Hobs = Line_26_4_F_div$Hobs,
                                        S1_Hobs = Line_26_4_S1_div$Hobs,
                                        S2_Hobs = Line_26_4_S2_div$Hobs,
                                        S3_Hobs = Line_26_4_S3_div$Hobs,
                                        S4_Hobs = Line_26_4_S4_div$Hobs)

colnames(generation_hobs_Line_26_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_26_4)
melted_generation_hobs_Line_26_4 <- melt(generation_hobs_Line_26_4)

generation_hobs_Line_26_4_summary <- data.frame(FS = mean(generation_hobs_Line_26_4$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_26_4$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_26_4$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_26_4$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_26_4$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_26_4_summary <- melt(generation_hobs_Line_26_4_summary)

ggplot(melted_generation_hobs_Line_26_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_26_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_26_4 <- melted_generation_hobs_Line_26_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_26_4_summary <- data.frame(FS = mean(expectations_halving_Line_26_4$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_26_4$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_26_4$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_26_4$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_26_4$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_26_4_summary_melted <- melt(expectations_halving_Line_26_4_summary)

obs_exp_26_4 <- rbind(expectations_halving_Line_26_4_summary_melted, melted_generation_hobs_Line_26_4_summary)

obs_exp_26_4_plot <- ggplot(obs_exp_26_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 26_4")

wilcox.test(generation_hobs_Line_26_4$S1, expectations_halving_Line_26_4$S1, paired = T)
wilcox.test(generation_hobs_Line_26_4$S2, expectations_halving_Line_26_4$S2, paired = T)
wilcox.test(generation_hobs_Line_26_4$S3, expectations_halving_Line_26_4$S3, paired = T)
wilcox.test(generation_hobs_Line_26_4$S4, expectations_halving_Line_26_4$S4, paired = T)

write.table(generation_hobs_Line_26_4, "Line_26_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_26_4, "Line_26_4_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 27_2
generation_hobs_Line_27_2 <- data.frame(F_Hobs = Line_27_2_F_div$Hobs,
                                        S1_Hobs = Line_27_2_S1_div$Hobs,
                                        S2_Hobs = Line_27_2_S2_div$Hobs,
                                        S3_Hobs = Line_27_2_S3_div$Hobs,
                                        S4_Hobs = Line_27_2_S4_div$Hobs)

colnames(generation_hobs_Line_27_2) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_27_2)
melted_generation_hobs_Line_27_2 <- melt(generation_hobs_Line_27_2)

generation_hobs_Line_27_2_summary <- data.frame(FS = mean(generation_hobs_Line_27_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_27_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_27_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_27_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_27_2$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_27_2_summary <- melt(generation_hobs_Line_27_2_summary)

ggplot(melted_generation_hobs_Line_27_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_27_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_27_2 <- melted_generation_hobs_Line_27_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_27_2_summary <- data.frame(FS = mean(expectations_halving_Line_27_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_27_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_27_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_27_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_27_2$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_27_2_summary_melted <- melt(expectations_halving_Line_27_2_summary)

obs_exp_27_2 <- rbind(expectations_halving_Line_27_2_summary_melted, melted_generation_hobs_Line_27_2_summary)

obs_exp_27_2_plot <- ggplot(obs_exp_27_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 27_2")

wilcox.test(generation_hobs_Line_27_2$S1, expectations_halving_Line_27_2$S1, paired = T)
wilcox.test(generation_hobs_Line_27_2$S2, expectations_halving_Line_27_2$S2, paired = T)
wilcox.test(generation_hobs_Line_27_2$S3, expectations_halving_Line_27_2$S3, paired = T)
wilcox.test(generation_hobs_Line_27_2$S4, expectations_halving_Line_27_2$S4, paired = T)

write.table(generation_hobs_Line_27_2, "Line_27_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_27_2, "Line_27_2_hexp.txt", col.names = T, row.names = F, quote = F)


# Line 29_2
generation_hobs_Line_29_2 <- data.frame(F_Hobs = Line_29_2_F_div$Hobs,
                                        S1_Hobs = Line_29_2_S1_div$Hobs,
                                        S2_Hobs = Line_29_2_S2_div$Hobs,
                                        S3_Hobs = Line_29_2_S3_div$Hobs,
                                        S4_Hobs = Line_29_2_S4_div$Hobs)

colnames(generation_hobs_Line_29_2) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_29_2)
melted_generation_hobs_Line_29_2 <- melt(generation_hobs_Line_29_2)

generation_hobs_Line_29_2_summary <- data.frame(FS = mean(generation_hobs_Line_29_2$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_29_2$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_29_2$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_29_2$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_29_2$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_29_2_summary <- melt(generation_hobs_Line_29_2_summary)

ggplot(melted_generation_hobs_Line_29_2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_29_2_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_29_2 <- melted_generation_hobs_Line_29_2 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_29_2_summary <- data.frame(FS = mean(expectations_halving_Line_29_2$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_29_2$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_29_2$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_29_2$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_29_2$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_29_2_summary_melted <- melt(expectations_halving_Line_29_2_summary)

obs_exp_29_2 <- rbind(expectations_halving_Line_29_2_summary_melted, melted_generation_hobs_Line_29_2_summary)

obs_exp_29_2_plot <- ggplot(obs_exp_29_2, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  theme(legend.position = "none") +
  ggtitle("Line 29_2")

wilcox.test(generation_hobs_Line_29_2$S1, expectations_halving_Line_29_2$S1, paired = T)
wilcox.test(generation_hobs_Line_29_2$S2, expectations_halving_Line_29_2$S2, paired = T)
wilcox.test(generation_hobs_Line_29_2$S3, expectations_halving_Line_29_2$S3, paired = T)
wilcox.test(generation_hobs_Line_29_2$S4, expectations_halving_Line_29_2$S4, paired = T)

write.table(generation_hobs_Line_29_2, "Line_29_2_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_29_2, "Line_29_2_hexp.txt", col.names = T, row.names = F, quote = F)

# Line 29_4
generation_hobs_Line_29_4 <- data.frame(F_Hobs = Line_29_4_F_div$Hobs,
                                        S1_Hobs = Line_29_4_S1_div$Hobs,
                                        S2_Hobs = Line_29_4_S2_div$Hobs,
                                        S3_Hobs = Line_29_4_S3_div$Hobs,
                                        S4_Hobs = Line_29_4_S4_div$Hobs)

colnames(generation_hobs_Line_29_4) <- c("F", "S1", "S2", "S3", "S4")
# generation_hobs <- rownames_to_column(generation_hobs_Line_29_4)
melted_generation_hobs_Line_29_4 <- melt(generation_hobs_Line_29_4)

generation_hobs_Line_29_4_summary <- data.frame(FS = mean(generation_hobs_Line_29_4$F, na.rm = T), 
                                                S1 = mean(generation_hobs_Line_29_4$S1, na.rm = T),
                                                S2 = mean(generation_hobs_Line_29_4$S2, na.rm = T),
                                                S3 = mean(generation_hobs_Line_29_4$S3, na.rm = T),
                                                S4 = mean(generation_hobs_Line_29_4$S4, na.rm = T),
                                                H = "Observed")

melted_generation_hobs_Line_29_4_summary <- melt(generation_hobs_Line_29_4_summary)

ggplot(melted_generation_hobs_Line_29_4, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_Publication()

ggplot(melted_generation_hobs_Line_29_4_summary, aes(x = variable, y = value)) +
  geom_point() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1))

expectations_halving_Line_29_4 <- melted_generation_hobs_Line_29_4 %>% 
  filter(variable == "F") %>% 
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16)) %>% 
  rename(F = value) %>% 
  select(F, S1, S2, S3, S4)

expectations_halving_Line_29_4_summary <- data.frame(FS = mean(expectations_halving_Line_29_4$F, na.rm = T), 
                                                     S1 = mean(expectations_halving_Line_29_4$S1, na.rm = T),
                                                     S2 = mean(expectations_halving_Line_29_4$S2, na.rm = T),
                                                     S3 = mean(expectations_halving_Line_29_4$S3, na.rm = T),
                                                     S4 = mean(expectations_halving_Line_29_4$S4, na.rm = T),
                                                     H = "Expected")
expectations_halving_Line_29_4_summary_melted <- melt(expectations_halving_Line_29_4_summary)

obs_exp_29_4 <- rbind(expectations_halving_Line_29_4_summary_melted, melted_generation_hobs_Line_29_4_summary)

obs_exp_29_4_plot <- ggplot(obs_exp_29_4, aes(x = variable, y = value, colour = H)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1)) +
  ylim(c(0,1)) +
  scale_colour_nejm() +
  ylab("Heterozygosity") +
  xlab("Generation") +
  theme_Publication() +
  ggtitle("Line 29_4")

wilcox.test(generation_hobs_Line_29_4$S1, expectations_halving_Line_29_4$S1, paired = T)
wilcox.test(generation_hobs_Line_29_4$S2, expectations_halving_Line_29_4$S2, paired = T)
wilcox.test(generation_hobs_Line_29_4$S3, expectations_halving_Line_29_4$S3, paired = T)
wilcox.test(generation_hobs_Line_29_4$S4, expectations_halving_Line_29_4$S4, paired = T)

write.table(generation_hobs_Line_29_4, "Line_29_4_hobs.txt", col.names = T, row.names = T, quote = F)
write.table(expectations_halving_Line_29_4, "Line_29_4_hexp.txt", col.names = T, row.names = F, quote = F)

all_line_heterozygosity_plots <- ggarrange(obs_exp_1_1_plot, obs_exp_6_1_plot, obs_exp_6_4_plot, obs_exp_7_2_plot, obs_exp_7_4_plot, 
                                           obs_exp_8_2_plot, obs_exp_8_4_plot, obs_exp_13_3_plot, obs_exp_13_4_plot, obs_exp_16_1_plot, 
                                           obs_exp_16_5_plot, obs_exp_17_2_plot, obs_exp_17_5_plot, obs_exp_19_2_plot, obs_exp_19_5_plot, 
                                           obs_exp_20_1_plot, obs_exp_20_4_plot, obs_exp_21_2_plot, obs_exp_21_6_plot, obs_exp_23_2_plot, 
                                           obs_exp_23_4_plot, obs_exp_26_1_plot, obs_exp_26_4_plot, obs_exp_27_2_plot, obs_exp_29_2_plot, 
                                           obs_exp_29_4_plot)

ggsave("all_line_heterozygosity_plots.svg", all_line_heterozygosity_plots, width = 25, height = 20)

### Combining lines for figure ###########

FS_S5_mean_het <- data.frame(rbind(generation_hobs_Line_1_1_summary, generation_hobs_Line_6_4_summary, generation_hobs_Line_17_2_summary,
                                   generation_hobs_Line_17_5_summary, generation_hobs_Line_21_2_summary, generation_hobs_Line_21_6_summary, 
                                   generation_hobs_Line_23_2_summary, generation_hobs_Line_29_2_summary, generation_hobs_Line_29_4_summary)) %>% 
  mutate(Line = factor(c("1_1", "6_4", "17_2", "17_5", "21_2", "21_6", "23_2", "29_2", "29_4")))

melted_FS_S5_mean_het <- melt(FS_S5_mean_het)

ggplot(melted_FS_S5_mean_het, aes(x = variable, y = value, fill = Line, colour = Line, group = Line)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(0, 1) +
  scale_fill_nejm() +
  theme_Publication()

### Getting expected het from FS observed het ################
generation_FS_S5_hobs <- data.frame(F_hobs = F_lines_div$Hobs,
                                    S1_hobs = S1_lines_div$Hobs,
                                    S2_hobs = S2_lines_div$Hobs,
                                    S3_hobs = S3_lines_div$Hobs,
                                    S4_hobs = S4_lines_div$Hobs,
                                    S5_hobs = S5_lines_div$Hobs)


het_snps <- generation_FS_S5_hobs %>%
  filter(S4_hobs == 1) 


clipr::write_clip(rownames(het_snps))

# plot het snps on LG
WRC_LG <- read.table("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/cedar/FINAL_DATA_SET/filtered_normalised/traits_snps_from_genome_paper/bayesR_maf_filter/WRC_chromosomes.txt", header = T)
# WRC_LG$chr <- factor(WRC_LG$chr)

het_snps_in_LG <- het_snps %>% 
  rownames_to_column() %>% 
  filter(rowname %in% WRC_LG$SNP) %>% 
  rename(SNP = rowname)



het_snps_in_LG <- merge(het_snps_in_LG,WRC_LG, by = "SNP")

het_snps_in_LG <- het_snps_in_LG %>% 
  arrange(chr,accurate_pos)

het_snps_in_LG$SNP <- factor(het_snps_in_LG$SNP, levels=unique(het_snps_in_LG$SNP))

het_snps_in_LG <- het_snps_in_LG %>% 
  select(SNP, chr, accurate_pos)

don_het_snps_in_LG <- het_snps_in_LG %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(accurate_pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(het_snps_in_LG, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, accurate_pos) %>%
  mutate(BPcum=accurate_pos+tot)

axisdf = don_het_snps_in_LG %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum)) / 2 )

het_snps_in_LG_plot <- ggplot(don_het_snps_in_LG, aes(x=BPcum, y = 1)) +
  
  # Show all points
  geom_point(alpha = 0.8) +
  scale_color_nejm() +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  # geom_vline(xintercept = c(888232888, 1683028619, 2539001849, 3225909214, 3875153489, 4496253278, 5139632346, 5771913333, 6413980409, 7022253572), linetype = "dashed") +
  xlab("Linkage Group") +
  
  # Custom the theme:
  theme_Publication() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

het_snps_in_LG_plot

colnames(generation_FS_S5_hobs) <- c("F", "S1", "S2", "S3", "S4", "S5")
generation_FS_S5_hobs <- rownames_to_column(generation_FS_S5_hobs)
melted_generation_FS_S5_hobs <- melt(generation_FS_S5_hobs)
# 
# ggplot(melted_generation_FS_S5_hobs, aes(x = value, fill = variable)) +
#   geom_density(alpha = 0.5) +
#   theme_Publication()
# 
# # Expected heterozygosities (Maybe not needed)
# generation_hexp <- data.frame(F_Hexp = F_lines_div$Hexp,
#                               S1_Hexp = S1_lines_div$Hexp,
#                               S2_Hexp = S2_lines_div$Hexp,
#                               S3_Hexp = S3_lines_div$Hexp,
#                               S4_Hexp = S4_lines_div$Hexp,
#                               S5_Hexp = S5_lines_div$Hexp)
# colnames(generation_hexp) <- c("F", "S1", "S2", "S3", "S4", "S5") 
# 
# # melted_hobs <- melt(generation_hobs)  
# melted_hexp <- melt(generation_hexp)  
# 
# write.table(melted_generation_FS_S5_hobs, "melted_hobs.txt", row.names = F, quote = F)
# write.table(melted_hexp, "melted_hexp.txt", row.names = F, quote = F)
# 
# # generation_heterozygosities <- data.frame(generation = rbind(melted_hobs$variable, melted_hexp$variable),
# #                                            het = c(rep("Obs", 229578), rep("Exp", 229578)))
# 
# generation_heterozygosities <- read.table("melted_hobs_hexp.txt")
# 
# ggplot(generation_heterozygosities, aes(x = V1, y = V2, fill = V3)) +
#   geom_boxplot() +
#   theme_Publication() +
#   scale_fill_nejm()
# 
expectations_halving <- melted_generation_FS_S5_hobs %>%
  filter(variable == "F") %>%
  mutate(S1 = (value / 2), S2 = (value / 4), S3 = (value / 8), S4 = (value/ 16), S5 = (value / 32)) %>%
  rename(F = value) %>%
  select(F, S1, S2, S3, S4, S5, variable)
# 
melted_expectations_halving <- melt(expectations_halving)

generation_heterozygosities_exp_halving <- read.table("S_lines_pop_gen/melted_expectations_halving_obs_exp_from_obs_FS_exp_removed.txt", header = F)
FS_S5_hets <- generation_heterozygosities_exp_halving %>% 
  rename(Heterozygosity = V4)


FS_S5_hets %>% 
  filter(Heterozygosity == "Observed") %>% 
  group_by(V2) %>% 
  summarise(mean = mean(V3), median = median(V3), sd = sd(V3))
# summarise_all(list(~mean(.),~sd(.), ~list(c(summary(.))))) %>% 
# unnest_wider(list)

FS_het <- FS_S5_hets %>% 
  filter(Heterozygosity == "Observed", V2 == "FS")

S1_het <- FS_S5_hets %>% 
  filter(Heterozygosity == "Observed", V2 == "S1")

S2_het <- FS_S5_hets %>% 
  filter(Heterozygosity == "Observed", V2 == "S2")

S3_het <- FS_S5_hets %>% 
  filter(Heterozygosity == "Observed", V2 == "S3")

S4_het <- FS_S5_hets %>% 
  filter(Heterozygosity == "Observed", V2 == "S4")

S5_het <- FS_S5_hets %>% 
  filter(Heterozygosity == "Observed", V2 == "S5")

## Plots for Figure 4
het_boxplot_exp_halving <- ggplot(FS_S5_hets, aes(x = V2, y = V3, fill = Heterozygosity)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_Publication() +
  scale_fill_nejm() +
  scale_y_continuous(name = "Heterozygosity", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("Generation") +
  scale_x_discrete("Generation", breaks = unique(FS_S5_hets$V2),
                   labels = c("FS",
                              "S1",
                              "S2",
                              "S3",
                              "S4",
                              "S5")) +
  theme(legend.title = element_blank(), legend.position = "top")
# theme(axis.text.x = element_text(vjust = -4))

het_violinplot_exp_halving <- ggplot(generation_heterozygosities_exp_halving, aes(x = V2, y = V3)) +
  geom_violin(aes(color = V4), trim = F, position = position_dodge(0.8), scale = "count", width = 1.1) +
  geom_boxplot(aes(color = V4), width = 0.05, position = position_dodge(0.8), outlier.alpha = 0.3) +
  # stat_summary(fun = mean, geom="point", shape=23, size=2, position = position_dodge(0.9)) +
  theme_Publication() +
  scale_color_nejm() +
  scale_y_continuous(name = "Heterozygosity", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("Generation")

het_violinplot_exp_halving <- ggplot(FS_S5_hets, aes(x = V2, y = V3)) +
  geom_violin(aes(color = Heterozygosity), trim = T, position = position_dodge(1), scale = "width") +
  geom_boxplot(aes(color = Heterozygosity), width = 0.15, position = position_dodge(1), outlier.alpha = 0.3) +
  # stat_summary(fun = mean, geom="point", shape=23, size=2, position = position_dodge(0.9)) +
  theme_Publication() +
  scale_color_nejm() +
  scale_y_continuous(name = "Heterozygosity", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("Generation") +
  scale_x_discrete("Generation", breaks = unique(FS_S5_hets$V2),
                   labels = c("FS (n = 28)",
                              "S1 (n = 28)",
                              "S2 (n = 28)",
                              "S3 (n = 28)",
                              "S4 (n = 28)",
                              "S5 (n = 11)"))

het_boxplot_exp_halving
het_violinplot_exp_halving

ggsave("S_lines_pop_gen/het_boxplot_exp_halving_exp_from_obs_legend_top.svg", het_boxplot_exp_halving, width = 10, height = 7)
ggsave("S_lines_pop_gen/het_violinplot_scale_width_S5.svg", het_violinplot_exp_halving, width = 10, height = 7)


generation_heterozygosities_obs <- generation_heterozygosities_exp_halving %>% 
  filter(V4 == "Observed")


generation_heterozygosities_exp <- generation_heterozygosities_exp_halving %>% 
  filter(V4 == "Expected")

summary(generation_heterozygosities_exp$V2)

het_two_cols <- read.table("S_lines_pop_gen/melted_expectations_halving_obs_exp_from_obs_two_cols.txt", header = T)

FS <- het_two_cols %>% 
  filter(Generation == "F")

mean(FS$Observed)
mean(FS$Expected)

S1 <- het_two_cols %>% 
  filter(Generation == "S1")

summary(S1$Observed)
summary(S1$Expected)

wilcox.test(S1$Obs, S1$Exp, paired = T, alternative = "greater")
mood.test(S1$Obs, S1$Exp, paired = T, alternative = "greater")

S2 <- het_two_cols %>% 
  filter(Generation == "S2")

mean(S2$Observed)
mean(S2$Expected)
summary(S2$Observed)
summary(S2$Expected)

wilcox.test(S2$Obs, S2$Exp, paired = T, alternative = "greater")
mood.test(S2$Obs, S2$Exp, paired = T, alternative = "greater")

S3 <- het_two_cols %>% 
  filter(Generation == "S3")

mean(S3$Observed)
mean(S3$Expected)

wilcox.test(S3$Obs, S3$Exp, paired = T, alternative = "greater")

S4 <- het_two_cols %>% 
  filter(Generation == "S4")

mean(S4$Observed)
mean(S4$Expected)

wilcox.test(S4$Obs, S4$Exp, paired = T, alternative = "greater")

S5 <- het_two_cols %>% 
  filter(Generation == "S5")

mean(S5$Observed)
mean(S5$Expected)

wilcox.test(S5$Obs, S5$Exp, paired = T, alternative = "less")

### Sign test ####
#FS
FS_het <- generation_heterozygosities_exp_halving %>% 
  filter(V2 =="F")

FS_het %>% 
  group_by(V4) %>%
  get_summary_stats(V3, type = "median_iqr")

bxp <- ggpaired(FS_het, x = "V4", y = "V3", 
                order = c("Observed", "Expected"),
                ylab = "Het", xlab = "group")
bxp

FS_test <- FS_het %>% 
  sign_test(V3 ~ V4) %>% 
  add_significance()
FS_test

#S1
S1_het <- generation_heterozygosities_exp_halving %>% 
  filter(V2 =="S1")

S1_het %>% 
  group_by(V4) %>%
  summary(V3)

bxp <- ggpaired(S1_het, x = "V4", y = "V3", 
                order = c("Observed", "Expected"),
                ylab = "Het", xlab = "group")
bxp

S1_test <- S1_het %>% 
  sign_test(V3 ~ 1, mu = 0.161, alternative = "less") %>% 
  add_significance()
S1_test

#S2
S2_het <- generation_heterozygosities_exp_halving %>% 
  filter(V2 =="S2", V4 == "Observed")

S2_het %>% 
  group_by(V4) %>%
  get_summary_stats(V3, type = "median_iqr")

bxp <- ggpaired(S2_het, x = "V4", y = "V3", 
                order = c("Observed", "Expected"),
                ylab = "Het", xlab = "group")
bxp

S2_test <- S2_het %>% 
  sign_test(V3 ~ 1, mu = 0.0803) %>% 
  add_significance()
S2_test

#S3
S3_het <- generation_heterozygosities_exp_halving %>% 
  filter(V2 =="S3", V4 == "Observed")

S3_het %>% 
  group_by(V4) %>%
  summary(V3)

bxp <- ggpaired(S3_het, x = "V4", y = "V3", 
                order = c("Observed", "Expected"),
                ylab = "Het", xlab = "group")
bxp

S3_test <- S3_het %>% 
  sign_test(V3 ~ 1, mu = 0.0401) %>% 
  add_significance()
S3_test

#S4
S4_het <- generation_heterozygosities_exp_halving %>% 
  filter(V2 =="S4", V4 == "Observed")

S4_het %>% 
  group_by(V4) %>%
  summary(V3)

bxp <- ggpaired(S4_het, x = "V4", y = "V3", 
                order = c("Observed", "Expected"),
                ylab = "Het", xlab = "group")
bxp

S4_test <- S4_het %>% 
  sign_test(V3 ~ 1, mu = 0.0201) %>% 
  add_significance()
S4_test

#S5
S5_het <- generation_heterozygosities_exp_halving %>% 
  filter(V2 =="S5", V4 == "Observed")

S5_het %>% 
  group_by(V4) %>%
  get_summary_stats(V3, type = "median_iqr")

bxp <- ggpaired(S5_het, x = "V4", y = "V3", 
                order = c("Observed", "Expected"),
                ylab = "Het", xlab = "group")
bxp

S5_test <- S5_het %>% 
  sign_test(V3 ~ 1, mu = 0.0100) %>% 
  add_significance()
S5_test


### Heterozygosities in diverse pop #####
parents_div <- summary(parents_genind)

diverse_pops_sep <- seppop(parents_genind)
VI_div <- summary(diverse_pops_sep$Vancouver_Island)
CB_div <- summary(diverse_pops_sep$Coastal_BC)
HG_div <- summary(diverse_pops_sep$Haida_Gwaii)
IB_div <- summary(diverse_pops_sep$Interior_BC)
US_div <- summary(diverse_pops_sep$US)


diverse_hobs <- data.frame(diverse_Hobs = parents_div$Hobs,
                           VI_Hobs = VI_div$Hobs,
                           CB_Hobs = CB_div$Hobs,
                           HG_Hobs = HG_div$Hobs,
                           IB_Hobs = IB_div$Hobs,                              
                           US_Hobs = US_div$Hobs)

melted_diverse_hobs <- melt(diverse_hobs)

diverse_selfing_het <- read.table("S_lines_pop_gen/melted_expectations_halving_obs_exp_from_obs_FS_exp_removed_w_diverse.txt", header = F)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


het_violinplot_diverse_pop <- ggplot(melted_diverse_hobs, aes(x = variable, y = value)) +
  geom_violin(trim = T, scale = "area") +
  geom_boxplot(width = 0.15, outlier.alpha = 0.3, position = "jitter") +
  # stat_summary(fun = mean, geom="point", shape=23, size=2, position = position_dodge(0.9)) +
  theme_Publication() +
  scale_color_nejm() +
  scale_y_continuous(name = "Heterozygosity", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("Population") +
  scale_x_discrete("population", breaks = unique(melted_diverse_hobs$variable),
                   labels = c(expression(paste("Diverse Population\n(",italic("n"),"= 112)")),
                              expression(paste("Vancouver Island\n(",italic("n")," = 63)")),
                              expression(paste("Coastal BC\n(",italic("n")," = 26)")),
                              expression(paste("Haida Gwaii\n(",italic("n")," = 16)")),
                              expression(paste("Interior BC\n(",italic("n")," = 3)")),
                              expression(paste("Coastal NW US\n(",italic("n")," = 4)"))))

het_boxplot_diverse_pop <- ggplot(melted_diverse_hobs, aes(x = variable, y = value)) +
  geom_boxplot(outlier.alpha = 0.3) +
  # geom_boxplot_jitter(outlier.jitter.width = 0.1, outlier.jitter.height = 0.05) +
  scale_x_discrete("Population", breaks = unique(melted_diverse_hobs$variable),
                   labels = c("Diverse Population (n = 112)",
                              "Vancouver Island (n = 63)",
                              "Coastal BC (n = 26)",
                              "Haida Gwaii (n = 16)",
                              "Interior BC (n = 3)",
                              ("Coastal NW US (n = 4)"))) +
  # stat_summary(fun = mean, geom="point", shape=23, size=2, position = position_dodge(0.9)) +
  theme_Publication() +
  scale_color_nejm() +
  scale_y_continuous(name = "Heterozygosity", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(axis.text.x = element_text(vjust = -4)) +
  xlab("Population")

het_boxplot_selfing_diverse <- ggplot(diverse_selfing_het, aes(x = V2, y = V3), fill = V4) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_Publication() +
  scale_fill_nejm() +
  scale_y_continuous(name = "Heterozygosity", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("Generation") +
  scale_x_discrete("Generation", breaks = unique(diverse_selfing_het$V2),
                   labels = c("Diverse Population (n = 112)",
                              "FS (n = 28)",
                              "S1 (n = 28)",
                              "S2 (n = 28)",
                              "S3 (n = 28)",
                              "S4 (n = 28)",
                              "S5 (n = 11)"))


het_boxplot_diverse_pop
het_boxplot_selfing_diverse
het_violinplot_diverse_pop

ggsave("het_boxplot_diverse_pop.svg", het_boxplot_diverse_pop, width = 10, height = 7)
ggsave("S_lines_pop_gen/het_violinplot_scale_width_S5.svg", het_violinplot_exp_halving, width = 10, height = 7)

