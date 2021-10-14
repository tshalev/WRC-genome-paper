library(adegenet)
library(geodist)
library(hierfstat)
library(vcfR)
library(pegas)
library(reshape2)
library(ggsci)
library(apex)
library(mmod)
library(factoextra)
library(FactoMineR)
library(poppr)
library(ggpubr)
library(rstatix)
library(tidyverse)

source("~/UBC/GSAT/PhD/WRC/r_scripts/publication_theme.r")

# Read in VCF

parents_mac_3_r2_01_vcf <- read.vcfR("diverse_pop.final.g95minQ30minmeanDP15maxmeanDP70AB28.HWE_het.mac3.1490k_r2_01.samples_sorted.recode.vcf")

# Read in popmap
clusters_popmap <- read.table("popmap_diverse_pop_clusters.txt")

# Convert VCF to genind
parents_mac_3_r2_01_genind <- vcfR2genind(parents_mac_3_r2_01_vcf, pop = clusters_popmap$V2)
parents_mac_3_r2_01_no_pop_genind <- vcfR2genind(parents_mac_3_r2_01_vcf, pop = clusters_popmap$V1)

# Geographic locations for Mantel test
all_geo_dist <- read.table("diverse_pop_clusters_geog_locations.txt", header = T, sep = "\t")
clusters_geo_dist <- read.table("clusters_average_coord.txt", header = T, row.names = 1)

# Stats
parents_mac_3_r2_01_stats <- genind2hierfstat(parents_mac_3_r2_01_genind)
parents_mac_3_r2_01_no_pop_stats <- genind2hierfstat(parents_mac_3_r2_01_no_pop_genind)


### FST tests ####
clusters <- clusters_popmap$V2

basic.stats(parents_mac_3_r2_01_genind)
wc(parents_mac_3_r2_01_stats)
parents_div <- summary(parents_mac_3_r2_01_genind)
loci <- parents_mac_3_r2_01_stats[, -1]
# global_varcomp_parents_mac_3_r2_01 <- varcomp.glob(levels = data.frame(clusters, pops), loci, diploid = T)
g_test_parents_mac_3_r2_01 <- test.g(loci, level = clusters, nperm = 1000)
# between_test_parents_mac_3_r2_01 <- test.between(loci, test.lev = clusters, rand.unit = pops, nperm = 999)
genet.dist(parents_mac_3_r2_01_genind, method = "WC84")


### Mantel tests ####
all_geo_dist <- all_geo_dist %>% 
  select(Parent, Longitude, Latitude) %>% 
  column_to_rownames(var = "Parent")

# Exact geo dist
all_geo_dist_mat <- geodist(all_geo_dist, measure = "geodesic", )
rownames(all_geo_dist_mat) <- row.names(all_geo_dist)
colnames(all_geo_dist_mat) <- row.names(all_geo_dist)

all_geo_dist_mat[upper.tri(all_geo_dist_mat)] <- ""
all_geo_dist_mat[upper.tri(all_geo_dist_mat,diag=TRUE)] <- ""

# Euclidean gen dist
all_gen_dist_mat <- dist(parents_mac_3_r2_01_no_pop_genind, method = "euclidean", diag = F, upper = F, p)

all_geo_dist_vector <- c(all_geo_dist_mat)
all_gen_dist_vector <- c(Dgen)

all_geo_dist_vector <- as.numeric(all_geo_dist_vector[all_geo_dist_vector != ""])

geo_gen_dist_df <- data.frame(geo = all_geo_dist_vector/1000, gen = all_gen_dist_vector)

# Plot individual correlation ge gen dist
ggplot(geo_gen_dist_df, aes(x = geo, y = gen)) +
  theme_Publication() +
  geom_point() +
  xlab("Geographic distance (km)") +
  ylab("Genetic distance") +
  geom_smooth(method = "lm")
  geom_abline(slope = 0.2408384, intercept = 0.417928)

ggsave("mantel_indiv_plot.svg", dpi = 300, width = 10, height = 9)  

m <- lm(gen ~ geo, geo_gen_dist_df)
summary(m)
anova(m)

# Mantel test clusters
parents_mac_3_r2_01_genpop <- genind2genpop(parents_mac_3_r2_01_genind)

Dgen <- dist.genpop(parents_mac_3_r2_01_genpop, method = 2)
Dgeo <- dist(cbind(clusters_geo_dist$Longitude, clusters_geo_dist$Latitude))
as.matrix(Dgeo)[1:5,1:5]

ibd <- mantel.randtest(Dgen, Dgeo, nrepet = 9999)

ibd
plot(ibd)

mcor <- vegan::mantel.correlog(Dgen, Dgeo)
mcor
plot(mcor)

Dgen <- as.numeric(Dgen)
Dgeo <- as.numeric(Dgeo)

par(mar=c(4,4,0,0))
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=1,  
     xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ Dgeo))

lines(loess.smooth(log(Dgeo), Dgen), col="red")


# Generating PCA for Figure 1C

parents_mac_3_r2_01_scaled <- scaleGen(parents_mac_3_r2_01_genind, NA.method = "mean")
pca_parents_mac_3_r2_01 <- dudi.pca(parents_mac_3_r2_01_scaled,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca_parents_mac_3_r2_01$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

col <- c(pal_nejm()(8))

get_eigenvalue(pca_parents_mac_3_r2_01)
fviz_eig(pca_parents_mac_3_r2_01, addlabels = T)
p <- fviz_pca_ind(pca_parents_mac_3_r2_01,
                  axes = c(1, 2),
                  col.ind = clusters_popmap$V2, # color by groups
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
  xlab("PC 1 (2.98%)") + 
  ylab("PC 2 (1.49%)") +
  scale_color_manual(name="Cluster",
                     breaks=c("Northern_Coastal", "Central", "Southern_Interior"),
                     labels=c("Northern-Coastal", "Central", "Southern-Interior"),
                     values = col) +
  theme(legend.position = c(0.16, 0.92), legend.background = element_rect(colour = "black"))

pca_plot
ggsave("pca_plot_pops_all_snps_clusters.tiff", pca_plot, width = 6, height = 6)


### DAPC ###########################################

# Adapted from https://grunwaldlab.github.io/Population_Genetics_in_R/
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(parents_mac_3_r2_01_genind, n.pca = 200, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_pubr()
p1 <- p1 + xlab("Number of groups (K)")
p1

my_k <- 2:5

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(200)
  grp_l[[i]] <- find.clusters(parents_mac_3_r2_01_genind, n.pca = 200, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(parents_mac_3_r2_01_genind, pop = grp_l[[i]]$grp, n.pca = 25, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

col <- c(pal_nejm()(3), "#7876B1FF", "#6F99ADFF")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_pubr()
p2 <- p2 + scale_color_manual(values=col)
p2 <- p2 + scale_fill_manual(values=col)
p2

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- clusters_popmap$V2
my_df <- tmp

for(i in 2:length(dapc_l)){
  set.seed(200)
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- clusters_popmap$V2
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_pubr()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Dark2")
p3 <- p3 + scale_fill_manual(values=col)
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3 <- p3 + xlab("Tree")
p3

library(ggpubr)
ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)

ggsave("dapc_multiplot.svg", dpi = 300, width = 15, height = 10)

## DAPC cross-validation ####

parents_mac_3_r2_01_cluster_xdapc <- xvalDapc(tab(parents_mac_3_r2_01_genind, NA.method = "mean"), clusters_popmap$V2, n.pca = 1:200)
parents_mac_3_r2_01_cluster_xdapc <- xvalDapc(tab(parents_mac_3_r2_01_genind, NA.method = "mean"), pop(parents_mac_3_r2_01_genind), n.pca = 1:70, n.rep = 1000, parallel = "snow", ncpus = 6L)


scatter(parents_mac_3_r2_01_cluster_xdapc$DAPC, cex = 2, legend = TRUE, 
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,)

scatter(parents_mac_3_r2_01_clusters_dapc, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

par(new=TRUE)

df <- data.frame(x = parents_mac_3_r2_01_cluster_xdapc$DAPC$ind.coord[,1], y = parents_mac_3_r2_01_cluster_xdapc$DAPC$ind.coord[,2])

s.label(dfxy = df, xax=1, yax=2, label=as.character(clusters_popmap$V1),
        clabel=0.7, # change the size of the labels
        boxes=TRUE, # if points are spaced wide enough, can use TRUE to add boxes around the labels
        grid=FALSE, addaxes=FALSE) # do not draw lines or axes in addition to the labels

mac3_r2_01_dapc <- parents_mac_3_r2_01_cluster_xdapc$DAPC

dapc.results <- as.data.frame(mac3_r2_01_dapc$posterior)
dapc.results$pop <- pop(parents_mac_3_r2_01_genind)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

col <- c(pal_nejm()(8), "#6A3D9A", "#A6CEE3")
  
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = col) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
