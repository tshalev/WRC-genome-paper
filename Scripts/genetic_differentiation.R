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

# Read in VCFs
parents_vcf <- read.vcfR("diverse_pop.g95minmeanDP15maxmeanDP60minQ30AB28.recode.vcf")
parents_intergenic_vcf <- read.vcfR("diverse_pop_all_intergenic_snps.recode.vcf")
parents_intergenic_r2_2_vcf <- read.vcfR("diverse_pop_all_intergenic_snps.r2_2.recode.vcf")

# Read in strata
parents_popmap <- read.table("parents_strata.txt", sep = "\t")

# Convert VCF to genind
parents_genind <- vcfR2genind(parents_vcf, pop = parents_popmap$V2)
parents_genind_subpop <- vcfR2genind(parents_vcf, pop = parents_popmap$V3)
parents_intergenic_genind <- vcfR2genind(parents_intergenic_vcf, pop = parents_popmap$V2)
parents_intergenic_genind_subpop <- vcfR2genind(parents_intergenic_vcf, pop = parents_popmap$V3)
parents_intergenic_r2_2_genind <- vcfR2genind(parents_intergenic_r2_2_vcf, pop = parents_popmap$V2)

# Stats
parents_stats <- genind2hierfstat(parents_genind)
parents_intergenic_stats <- genind2hierfstat(parents_intergenic_genind)

# FST tests
population <- parents_popmap$V2
subpop <- parents_popmap$V3
basic.stats(parents_genind)
wc(parents_stats)
parents_div <- summary(parents_genind)
loci <- parents_stats[, -1]
global_varcomp <- varcomp.glob(levels = data.frame(population, subpop), loci, diploid = T)
g_test <- test.g(loci, level = population)
between_test <- test.between(loci, test.lev = population, rand.unit = subpop, nperm = 100)
genet.dist(parents_genind, method = "WC84")
genet.dist(parents_genind_subpop, method = "WC84")

basic.stats(parents_maf05_r2_2_genind)
wc(parents_maf05_r2_0_stats)
parents_div_all <- summary(parents_genind_all)
loci <- parents_intergenic_stats[, -1]
global_varcomp_intergenic <- varcomp.glob(levels = data.frame(population, subpop), loci, diploid = T)
g_test_intergenic <- test.g(loci, level = population)
between_test_intergenic <- test.between(loci, test.lev = population, rand.unit = subpop, nperm = 100)
genet.dist(parents_intergenic_genind, method = "WC84")
genet.dist(parents_intergenic_genind_subpop, method = "WC84")

population <- as.character(HG_INT_popmap$V2)
basic.stats(parents_genind_all)
wc(parents_stats_all)
parents_div_all <- summary(parents_genind_all)
loci <- HG_INT_stats[, -1]
global_varcomp <- varcomp.glob(levels = data.frame(population), loci, diploid = T)
g_test <- test.g(loci, level = population)
genet.dist(parents_genind_all, method = "WC84")

# Some PCAs for genetic distance
a <- indpca(all_uncorrected_generation_stats)

pdf(file = "all_uncorrected_generations_PCA.pdf")
plot(a, cex = 0.2)
dev.off()

b <- indpca(all_uncorrected_FS_line_stats)

pdf(file = "all_uncorrected_FS_line_PCA.pdf")
plot(b, cex = 0.2)
dev.off()

c <- indpca(all_uncorrected_parental_line_stats)

pdf(file = "all_uncorrected_parental_line_PCA.pdf")
plot(c, cex = 0.2)
dev.off()

d <- indpca(parents_LD_pruned_stats)
pdf(file = "all_uncorrected_parents_PCA.pdf")
plot(d, cex = 0.7)
dev.off()

d2 <- indpca(parents_uncorrected_stats)
pdf(file = "all_uncorrected_parents_PCA.pdf", width = 4, height = 5)
plot(d2, cex = 0.5)
dev.off()

e <- indpca(S_uncorrected_parental_line_stats)
pdf(file = "S_uncorrected_parental_line_PCA_coloured.pdf")
plot(e, cex = 0.5, col = S_uncorrected_parental_line_popmap$V2, )
dev.off()

f <- indpca(S_uncorrected_FS_line_stats)
pdf(file = "S_uncorrected_FS_lines_PCA.pdf")
plot(f, cex = 0.5)
dev.off()

g <- indpca(S_uncorrected_generation_stats)
pdf(file = "S_uncorrected_generations_PCA.pdf")
plot(g, cex = 0.5)
dev.off()

# Generating PCAs for Figure 2
parents_scaled <- scaleGen(parents_genind_subpop, NA.method = "mean")
pca_parents_subpop <- dudi.pca(parents_scaled,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca_parents_intergenic$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

parents_intergenic_scaled <- scaleGen(parents_intergenic_genind, NA.method = "mean")
pca_parents_intergenic <- dudi.pca(parents_intergenic_scaled,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca_parents_intergenic$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

parents_intergenic_r2_2_scaled <- scaleGen(parents_intergenic_r2_2_genind, NA.method = "mean")
pca_parents_intergenic_r2_2 <- dudi.pca(parents_intergenic_r2_2_scaled,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca_parents_intergenic_r2_2$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

col <- RColorBrewer::brewer.pal(10, "Paired")
col <- c(pal_nejm()(8), "#6A3D9A", "#A6CEE3")

get_eigenvalue(pca_parents)
fviz_eig(pca_parents, addlabels = T)
p <- fviz_pca_ind(pca_parents,
                  axes = c(1, 2),
             col.ind = parents_popmap$V2, # color by groups
             pointsize = 4,
             palette = col,
             # addEllipses = TRUE, # Concentration ellipse
             # ellipse.type = "confidence",
             # legend.title = "Population",
             # repel = T,
             label = "none",
             mean.point = FALSE,
             pointshape = 19
)

pca_plot <- p + 
  theme_Publication() +
  labs(title = NULL) +
  xlab("PC 1 (3.2%)") + 
  ylab("PC 2 (1.7%)") +
  scale_color_manual(name="Population",
                    breaks=c("Coastal_BC", "Haida_Gwaii", "Interior_BC", "US", "Vancouver_Island"),
                    labels=c("Coastal BC", "Haida Gwaii", "Interior BC", "Coastal NW US", "Vancouver Island"),
                    values = col) +
  theme(legend.position = c(0.9, 0.92), legend.background = element_rect(colour = "black"))
  
pca_plot
ggsave("pca_plot_pops_all_snps_pops.svg", pca_plot, width = 10, height = 9)

p <- fviz_pca_ind(pca_parents_subpop,
                  axes = c(1, 2),
                  col.ind = parents_popmap$V3, # color by groups
                  pointsize = 4,
                  palette = col,
                  # addEllipses = TRUE, # Concentration ellipse
                  # ellipse.type = "confidence",
                  # legend.title = "Population",
                  # repel = T,
                  label = "none",
                  mean.point = FALSE,
                  pointshape = 19
)

pca_plot <- p + 
  theme_Publication() +
  labs(title = NULL) +
  xlab("PC 1 (3.2%)") + 
  ylab("PC 2 (1.7%)") +
  scale_color_manual(name="Subpopulation",
                    breaks=c("Central_Coast", "Discovery_Islands", "Haida_Gwaii", "Interior_BC", "Lower_Mainland", "Sunshine_Coast", "US",
                             "Vancouver_Island_Centre", "Vancouver_Island_North", "Vancouver_Island_South"),
                    labels=c("Central Coast", "Discovery Islands", "Haida Gwaii", "Interior BC", "Lower Mainland", "Sunshine Coast", "Coastal NW US",
                             "Vancouver Island (Centre)", "Vancouver Island (North)", "Vancouver Island (South)"),
                    values = col) +
  theme(legend.position = c(0.86, 0.87), legend.background = element_rect(colour = "black"))

pca_plot
ggsave("pca_plot_pops_all_snps_subpops.svg", pca_plot, width = 10, height = 9)

# DAPC ###########################################
grp <- find.clusters(parents_genind, max.n.clust = 20, n.pca = 111, choose.n.clust = F)
grp$grp

dapc1 <- dapc(parents_intergenic_genind, grp$grp, n.pca = 96)
scatter(dapc1)

parents_clusters <- find.clusters(parents_genind_all)

all_uncorrected_clusters <- find.clusters(all_uncorrected_genind)

parents_k_cluster_dapc <- dapc(parents_genind_all, pop = parents_clusters$grp)
parents_pops_dapc <- dapc(parents_uncorrected_genind, pop = parents_uncorrected_popmap$V2)
scatter(parents_k_cluster_dapc)
scatter(parents_pops_dapc)

all_uncorrected_dapc <- dapc(all_uncorrected_genind, pop = all_uncorrected_FS_line_popmap$V2)
scatter(all_uncorrected_dapc)

S_uncorrected_clusters <- find.clusters(S_uncorrected_genind)
table.value(table(S_uncorrected_parental_line_popmap$V2, S_uncorrected_clusters$grp), col.lab=paste("inf", 1:18),row.lab=paste("ori", 1:15))

S_uncorrected_clusters_dapc <- dapc(S_uncorrected_genind, pop = S_uncorrected_clusters$grp)
scatter(S_uncorrected_clusters_dapc)

S_uncorrected_parental_line_dapc <- dapc(S_uncorrected_genind, pop = S_uncorrected_parental_line_popmap$V2)

scatter(S_uncorrected_parental_line_dapc)

scatter(S_uncorrected_clusters_dapc, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:12))
scatter(all_uncorrected_dapc, posi.da="bottomleft", bg="white", pch=1:12)

myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(S_uncorrected_clusters_dapc, posi.da="bottomright",  bg="white",pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,posi.pca="bottomleft")


pop(S_uncorrected_genind) <- S_uncorrected_clusters$grp
dapc.S_line_uncorrected <- dapc(S_uncorrected_genind)
myPal <- colorRampPalette(c("blue","gold","red"))
scatter(dapc.S_line_uncorrected, col=transp(myPal(18)), scree.da=FALSE,cell=1.5, cex=2, bg="white",cstar=0)
set.seed(4)
contrib <- loadingplot(dapc.S_line_uncorrected$var.contr, axis=2,thres=.07, lab.jitter=1)

compoplot(dapc.S_line_uncorrected, posi="bottomright",txt.leg=paste("Cluster", 1:18), lab="",ncol=1, xlab="individuals", col=funky(18))
compoplot(dapc.S_line_uncorrected, subset = 1:50, posi="bottomright",txt.leg=paste("Cluster", 1:18), lab="",ncol=1, xlab="individuals", col=funky(18))

temp <- which(apply(dapc.S_line_uncorrected$posterior,1, function(e) all(e<0.9)))
temp

compoplot(dapc.S_line_uncorrected, subset=temp, posi="bottomright",txt.leg=paste("Cluster", 1:15),ncol=2, col=funky(15))
    
write.table(dapc.S_line_uncorrected$grp, "cluster_groups.txt", quote = F, row.names = F, col.names = F)      

parents_uncorrected_clusters <- find.clusters(parents_uncorrected_genind)
pop(parents_uncorrected_genind) <- parents_uncorrected_popmap$V2
dapc.parents_uncorrected <- dapc(parents_uncorrected_genind)
myPal <- colorRampPalette(c("blue","gold","red"))
scatter(dapc.parents_uncorrected, col=transp(myPal(2)), scree.da=FALSE,cell=1.5, cex=2, bg="white",cstar=0)
set.seed(4)
contrib <- loadingplot(dapc.parents_uncorrected$var.contr, axis=2,thres=.07, lab.jitter=1)

compoplot(dapc.parents_uncorrected, posi="bottomright",txt.leg=paste("Cluster", 1:2), lab="",ncol=1, xlab="individuals", col=funky(2))
compoplot(dapc.parents_uncorrected, subset = 1:50, posi="bottomright",txt.leg=paste("Cluster", 1:2), lab="",ncol=1, xlab="individuals", col=funky(2))

temp <- which(apply(dapc.parents_uncorrected$posterior,1, function(e) all(e<0.9)))
temp

compoplot(dapc.S_line_uncorrected, subset=temp, posi="bottomright",txt.leg=paste("Cluster", 1:15),ncol=2, col=funky(15))

write.table(dapc.parents_uncorrected$grp, "cluster_groups_parents.txt", quote = F, row.names = F, col.names = F)      


### Test for significant differentiation #####
parents_pop_genind <- vcfR2genind(parents_vcf, pop = parents_uncorrected_popmap$V2)
BC_parents_genind <- vcfR2genind(BC_vcf, pop = BC_parents_popmap$V2)

### AMOVA ######
strata <- read.table("strata_for_poppr.txt", sep = "\t", header = T)
parents_genind <- vcfR2genind(parents_vcf, return.alleles = T, pop = parents_popmap$V2)
parents_intergenic_genind <- vcfR2genind(parents_intergenic_vcf, return.alleles = T, pop = parents_popmap$V2)

strata(parents_genind) <- data.frame(strata)
strata(parents_intergenic_genind) <- data.frame(strata)

parents_genclone <- as.genclone(parents_genind)
parents_intergenic_genclone <- as.genclone(parents_intergenic_genind)
parents_genclone

table(strata(parents_genind, ~Pop))
table(strata(parents_genind, ~Pop/Subpop, combine = FALSE))

table(strata(parents_intergenic_genind, ~Pop))
table(strata(parents_intergenic_genind, ~Pop/Subpop, combine = FALSE))

# AMOVA
parents_amova <- poppr.amova(parents_genclone, ~Pop/Subpop)
parents_intergenic_amova <- poppr.amova(parents_intergenic_genclone, ~Pop/Subpop)

# parents_amovacc <- poppr.amova(parents_genclone, ~Pop/Subpop, clonecorrect = T)

parents_amova
parents_intergenic_amova
parents_amovacc

parents_signif <- randtest(parents_amova, nrepet = 1000, p.adjust.method = "holm")
parents_intergenic_signif <- randtest(parents_intergenic_amova, nrepet = 1000, p.adjust.method = "holm")

# parents_signif_cc <- randtest(parents_amovacc, nrepet = 1000, p.adjust.method = "holm")
parents_signif
parents_intergenic_signif

parents_signif_cc
plot(parents_signif)


# Separating pops
sep_pop <- seppop(parents_genind)
names(sep_pop)
VI_BC_COAST_BC <- repool(sep_pop$VI_BC, sep_pop$COAST_BC)
VI_BC_HG_BC <- repool(sep_pop$VI_BC, sep_pop$HG_BC)
VI_BC_INT_BC <- repool(sep_pop$VI_BC, sep_pop$INT_BC)
VI_BC_CA_OR <- repool(sep_pop$VI_BC, sep_pop$CA_OR)
COAST_BC_HG_BC <- repool(sep_pop$COAST_BC, sep_pop$HG_BC)
COAST_BC_INT_BC <- repool(sep_pop$COAST_BC, sep_pop$INT_BC)
COAST_BC_CA_OR <- repool(sep_pop$COAST_BC, sep_pop$CA_OR)
HG_BC_INT_BC <- repool(sep_pop$HG_BC, sep_pop$INT_BC)
HG_BC_CA_OR <- repool(sep_pop$HG_BC, sep_pop$CA_OR)
INT_BC_CA_OR <- repool(sep_pop$INT_BC, sep_pop$CA_OR)

VI_BC_COAST_BC_genclone <- as.genclone(VI_BC_COAST_BC)
VI_BC_HG_BC_genclone <- as.genclone(VI_BC_HG_BC)
VI_BC_INT_BC_genclone <- as.genclone(VI_BC_INT_BC)
VI_BC_CA_OR_genclone <- as.genclone(VI_BC_CA_OR)
COAST_BC_HG_BC_genclone <- as.genclone(COAST_BC_HG_BC)
COAST_BC_INT_BC_genclone <- as.genclone(COAST_BC_INT_BC)
COAST_BC_CA_OR_genclone <- as.genclone(COAST_BC_CA_OR)
HG_BC_INT_BC_genclone <- as.genclone(HG_BC_INT_BC)
HG_BC_CA_OR_genclone <- as.genclone(HG_BC_INT_BC)
INT_BC_CA_OR_genclone <- as.genclone(INT_BC_CA_OR) 

VI_BC_COAST_BC_amova <- poppr.amova(VI_BC_COAST_BC_genclone, ~parents_genind.pop)
VI_BC_HG_BC_amova <- poppr.amova(VI_BC_HG_BC_genclone, ~parents_genind.pop)
VI_BC_INT_BC_amova <- poppr.amova(VI_BC_INT_BC_genclone, ~parents_genind.pop)
VI_BC_CA_OR_amova <- poppr.amova(VI_BC_CA_OR_genclone, ~parents_genind.pop)
COAST_BC_HG_BC_amova <- poppr.amova(COAST_BC_HG_BC_genclone, ~parents_genind.pop)
COAST_BC_INT_BC_amova <- poppr.amova(COAST_BC_INT_BC_genclone, ~parents_genind.pop)
COAST_BC_CA_OR_amova <- poppr.amova(COAST_BC_CA_OR_genclone, ~parents_genind.pop)
HG_BC_INT_BC_amova <- poppr.amova(HG_BC_INT_BC_genclone, ~parents_genind.pop)
HG_BC_CA_OR_amova <- poppr.amova(HG_BC_CA_OR_genclone, ~parents_genind.pop)
INT_BC_CA_OR_amova <- poppr.amova(INT_BC_CA_OR_genclone, ~parents_genind.pop)

VI_BC_COAST_BC_amova
VI_BC_COAST_BC_signif <- randtest(VI_BC_COAST_BC_amova, nrepet = 1000, p.adjust.method = "holm")
VI_BC_COAST_BC_signif
plot(VI_BC_COAST_BC_signif)

VI_BC_HG_BC_amova
VI_BC_HG_BC_signif <- randtest(VI_BC_HG_BC_amova, nrepet = 1000, p.adjust.method = "holm")
VI_BC_HG_BC_signif
plot(VI_BC_HG_BC_signif)

VI_BC_INT_BC_amova
VI_BC_INT_BC_signif <- randtest(VI_BC_INT_BC_amova, nrepet = 1000, p.adjust.method = "holm")
VI_BC_INT_BC_signif
plot(VI_BC_INT_BC_signif)

VI_BC_CA_OR_amova
VI_BC_CA_OR_signif <- randtest(VI_BC_CA_OR_amova, nrepet = 1000, p.adjust.method = "holm")
VI_BC_CA_OR_signif
plot(VI_BC_CA_OR_BC_signif)

COAST_BC_HG_BC_amova
COAST_BC_HG_BC_signif <- randtest(COAST_BC_HG_BC_amova, nrepet = 1000, p.adjust.method = "holm")
COAST_BC_HG_BC_signif
plot(COAST_BC_HG_BC_signif)

COAST_BC_INT_BC_amova
COAST_BC_INT_BC_signif <- randtest(COAST_BC_INT_BC_amova, nrepet = 1000, p.adjust.method = "holm")
COAST_BC_INT_BC_signif
plot(COAST_BC_INT_BC_signif)

COAST_BC_CA_OR_amova
COAST_BC_CA_OR_signif <- randtest(COAST_BC_CA_OR_amova, nrepet = 1000, p.adjust.method = "holm")
COAST_BC_CA_OR_signif
plot(COAST_BC_CA_OR_signif)

HG_BC_INT_BC_amova
HG_BC_INT_BC_signif <- randtest(HG_BC_INT_BC_amova, nrepet = 1000, p.adjust.method = "holm")
HG_BC_INT_BC_signif
plot(HG_BC_INT_BC_signif)

HG_BC_CA_OR_amova
HG_BC_CA_OR_signif <- randtest(HG_BC_CA_OR_amova, nrepet = 1000, p.adjust.method = "holm")
HG_BC_CA_OR_signif
plot(HG_BC_CA_OR_signif)

INT_BC_CA_OR_amova
INT_BC_CA_OR_signif <- randtest(INT_BC_CA_OR_amova, nrepet = 1000, p.adjust.method = "holm")
INT_BC_CA_OR_signif
plot(INT_BC_CA_OR_signif)
