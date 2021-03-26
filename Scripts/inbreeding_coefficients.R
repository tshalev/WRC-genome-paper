rm(list = ls())

library(asreml)
library(asremlPlus)
library(nadiv)
library(synbreed)
library(sommer)
library(ggsci)
library(tidyverse)

source("~/UBC/GSAT/PhD/WRC/r_scripts/publication_theme.r")

library(extrafont)
library(extrafontdb)

setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/")

diverse_ibc <- read.table("diverse_pop.g95minmeanDP15maxmeanDP60minQ30AB28.ibc", header = T)
diverse_r202_ibc <- read.table("diverse_pop.g95minmeanDP15maxmeanDP60minQ30AB28r202.ibc", header = T)

diverse_r202_ibc$pop_simple <- factor(diverse_r202_ibc$pop_simple, levels = c("Diverse_Population", "Vancouver_Island", "Coastal_BC", "Haida_Gwaii", 
                                                                              "Interior_BC", "Coastal_NW_US"))

### IBC from corrected allele matrices - S lines #################
FS_lines_ibc <- read.table("S_lines_pop_gen/F_lines.ibc", header = T)
FS_lines_ibc <- FS_lines_ibc %>% 
  mutate(Generation = rep("FS", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                             "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S1_lines_ibc <- read.table("S_lines_pop_gen/S1_lines.ibc", header = T)
S1_lines_ibc <- S1_lines_ibc %>% 
  mutate(Generation = rep("S1", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                             "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S2_lines_ibc <- read.table("S_lines_pop_gen/S2_lines.ibc", header = T)
S2_lines_ibc <- S2_lines_ibc %>% 
  mutate(Generation = rep("S2", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S3_lines_ibc <- read.table("S_lines_pop_gen/S3_lines.ibc", header = T)
S3_lines_ibc <- S3_lines_ibc %>% 
  mutate(Generation = rep("S3", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S4_lines_ibc <- read.table("S_lines_pop_gen/S4_lines.ibc", header = T)
S4_lines_ibc <- S4_lines_ibc %>% 
  mutate(Generation = rep("S4", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S5_lines_ibc <- read.table("S_lines_pop_gen/S5_lines.ibc", header = T)
S5_lines_ibc <- S5_lines_ibc %>% 
  mutate(Generation = rep("S5", 11), Line = c("1_1", "6_4", "10_2", "16_5", "17_2", "17_5", "21_2", "21_6", "23_2", "29_2", "29_4"))

S_lines_ibc <- rbind(FS_lines_ibc, S1_lines_ibc, S2_lines_ibc, S3_lines_ibc, S4_lines_ibc, S5_lines_ibc)
S_lines_ibc <- S_lines_ibc %>% 
  mutate(Fhat3_noneg = if_else(Fhat3 < 0, 0, Fhat3))

### IBC from corrected allele matrices - S lines - r2 < 0.2 #################
FS_lines_r202_ibc <- read.table("S_lines_pop_gen/F_lines_r202.ibc", header = T)
FS_lines_r202_ibc <- FS_lines_r202_ibc %>% 
  mutate(Generation = rep("FS", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                             "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S1_lines_r202_ibc <- read.table("S_lines_pop_gen/S1_lines_r202.ibc", header = T)
S1_lines_r202_ibc <- S1_lines_r202_ibc %>% 
  mutate(Generation = rep("S1", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S2_lines_r202_ibc <- read.table("S_lines_pop_gen/S2_lines_r202.ibc", header = T)
S2_lines_r202_ibc <- S2_lines_r202_ibc %>% 
  mutate(Generation = rep("S2", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S3_lines_r202_ibc <- read.table("S_lines_pop_gen/S3_lines_r202.ibc", header = T)
S3_lines_r202_ibc <- S3_lines_r202_ibc %>% 
  mutate(Generation = rep("S3", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S4_lines_r202_ibc <- read.table("S_lines_pop_gen/S4_lines_r202.ibc", header = T)
S4_lines_r202_ibc <- S4_lines_r202_ibc %>% 
  mutate(Generation = rep("S4", 28), Line = c("1_1", "6_1", "6_4", "7_2", "7_4", "8_2", "8_4", "10_2", "10_5", "13_3", "13_4", "16_1", "16_5", "17_2", 
                                              "17_5", "19_2", "19_5", "20_1", "20_4", "21_2", "21_6", "23_2", "23_4", "26_1", "26_4", "27_2", "29_2", "29_4"))

S5_lines_r202_ibc <- read.table("S_lines_pop_gen/S5_lines_r202.ibc", header = T)
S5_lines_r202_ibc <- S5_lines_r202_ibc %>% 
  mutate(Generation = rep("S5", 11), Line = c("1_1", "6_4", "10_2", "16_5", "17_2", "17_5", "21_2", "21_6", "23_2", "29_2", "29_4"))

S_lines_r202_ibc <- rbind(FS_lines_r202_ibc, S1_lines_r202_ibc, S2_lines_r202_ibc, S3_lines_r202_ibc, S4_lines_r202_ibc, S5_lines_r202_ibc)
S_lines_r202_ibc <- S_lines_r202_ibc %>% 
  mutate(Fhat3_noneg = if_else(Fhat3 < 0, 0, Fhat3), exp = c(rep(0, 28), rep(0.5, 28), rep(0.75, 28), rep(0.875, 28), rep(0.9375, 28), rep(0.96875, 11)))

## Summary stats
S_lines_ibc %>% 
  group_by(Generation) %>% 
  summary(Fhat3_noneg = median(Fhat3_noneg, na.rm = T))

S_lines_r202_ibc %>% 
  group_by(Generation) %>% 
  summarise(mean = mean(Fhat3_noneg), median = median(Fhat3_noneg))

diverse_ibc %>% 
  group_by(pop_simple) %>% 
  summarise(Fhat3_noneg = mean(Fhat3_noneg, na.rm = T))

diverse_r202_ibc %>% 
  # group_by(pop_simple) %>% 
  summarise(Fhat3_noneg = mean(Fhat3_noneg, na.rm = T), median = median(Fhat3_noneg, na.rm = T))

cor(S_lines_ibc$Fhat3_noneg, S_lines_r202_ibc$Fhat3_noneg)

# Plots for Figure 4
fhat_gen_violin <- ggplot(S_lines_r202_ibc, aes(x = Generation, y = Fhat3_noneg)) +
  geom_violin(scale = "area") +
  geom_point(aes(x = Generation, y = exp), shape = 23, fill = "red") +
  geom_boxplot(width = 0.15) +
  theme_Publication() +
  scale_color_nejm() +
  scale_x_discrete("Generation", breaks = unique(S_lines_r202_ibc$Generation),
                   labels = c("FS",
                              "S1",
                              "S2",
                              "S3",
                              "S4",
                              "S5")) +
  scale_y_continuous(name = expression(paste(bolditalic("F"))), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),limits = c(0, 1))

fhat_gen_violin

fhat_pop_violin <- ggplot(diverse_r202_ibc, aes(x = pop_simple, y = Fhat3_noneg)) +
  geom_violin(scale = "area") +
  geom_boxplot(width = 0.15) +
  theme_Publication() +
  scale_color_nejm() +
  scale_x_discrete("Population",
                   labels = c("Diverse Population (n = 112)",
                              "Vancouver Island (n = 63)",
                              "Coastal BC (n = 26)",
                              "Haida Gwaii (n = 16)",
                              "Interior BC (n = 3)",
                              "Coastal NW US (n = 4)")) +
  scale_y_continuous(name = expression(paste(bolditalic("F"))), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0,0.6))

fhat_pop_violin

ggsave("S_lines_pop_gen/fhat3_S_violin.svg", fhat_gen_violin, dpi = 300, width = 10, height = 7)
ggsave("fhat3_pop_violin.svg", fhat_pop_violin, dpi = 300, width = 10, height = 7)
