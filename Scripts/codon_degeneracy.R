library(tidyverse)

codon_changes <- read.table("codon_changes.txt", header = T)

codon_changes <- codon_changes %>% 
  mutate(a = str_replace(Codon1, "[ATCG]", "A"), t = str_replace(Codon1, "[ATCG]", "T"), c = str_replace(Codon1, "[ATCG]", "C"), g = str_replace(Codon1, "[ATCG]", "G"))

codon_changes$a <- toupper(codon_changes$a)
codon_changes$t <- toupper(codon_changes$t)
codon_changes$c <- toupper(codon_changes$c)
codon_changes$g <- toupper(codon_changes$g)

codon_changes <- codon_changes %>% 
  mutate(aa_a = if_else(a == "GCT" | a == "GCC" | a == "GCA" | a == "GCG", "A",  
                        if_else(a == "CGT" | a == "CGC" | a =="CGA" | a == "CGG" | a == "AGA" | a == "AGG", "R",
                                if_else(a == "AAT" | a == "AAC", "N", 
                                        if_else(a == "GAT" | a == "GAC", "D", 
                                                if_else(a == "TGT" | a == "TGC", "C",
                                                        if_else(a == "CAA" | a == "CAG", "Q", 
                                                                if_else(a == "GAA" | a == "GAG", "E",
                                                                        if_else(a == "GGT" | a == "GGC" | a =="GGA" | a == "GGG", "G",
                                                                                if_else(a == "CAT" | a == "CAC", "H",
                                                                                        if_else(a == "ATT" | a == "ATC" | a =="ATA", "I",
                                                                                                if_else(a == "CTT" | a == "CTC" | a =="CTA" | a == "CTG" | a == "TTA" | a == "TTG", "L",
                                                                                                        if_else(a == "AAA" | a == "AAG", "K",
                                                                                                                if_else(a == "ATG", "start",
                                                                                                                        if_else(a == "TTT" | a == "TTC", "F",
                                                                                                                                if_else(a == "CCT" | a == "CCC" | a =="CCA" | a == "CCG", "P",
                                                                                                                                        if_else(a == "TCT" | a == "TCC" | a =="TCA" | a == "TCG" | a == "AGT" | a == "AGC", "S",
                                                                                                                                                if_else(a == "ACT" | a == "ACC" | a =="ACA" | a == "ACG", "T",
                                                                                                                                                        if_else(a == "TGG", "W",
                                                                                                                                                                if_else(a == "TAT" | a == "TAC", "Y",
                                                                                                                                                                        if_else(a == "GTT" | a == "GTC" | a =="GTA" | a == "GTG", "V",
                                                                                                                                                                                if_else(a == "TAA" | a == "TGA" | a =="TAG", "stop", "unknown")))))))))))))))))))))) %>% 
  mutate(aa_g = if_else(g == "GCT" | g == "GCC" | g == "GCA" | g == "GCG", "A",  
                        if_else(g == "CGT" | g == "CGC" | g =="CGA" | g == "CGG" | g == "AGA" | g == "AGG", "R",
                                if_else(g == "AAT" | g == "AAC", "N", 
                                        if_else(g == "GAT" | g == "GAC", "D", 
                                                if_else(g == "TGT" | g == "TGC", "C",
                                                        if_else(g == "CAA" | g == "CAG", "Q", 
                                                                if_else(g == "GAA" | g == "GAG", "E",
                                                                        if_else(g == "GGT" | g == "GGC" | g =="GGA" | g == "GGG", "G",
                                                                                if_else(g == "CAT" | g == "CAC", "H",
                                                                                        if_else(g == "ATT" | g == "ATC" | g =="ATA", "I",
                                                                                                if_else(g == "CTT" | g == "CTC" | g =="CTA" | g == "CTG" | g == "TTA" | g == "TTG", "L",
                                                                                                        if_else(g == "AAA" | g == "AAG", "K",
                                                                                                                if_else(g == "ATG", "start",
                                                                                                                        if_else(g == "TTT" | g == "TTC", "F",
                                                                                                                                if_else(g == "CCT" | g == "CCC" | g =="CCA" | g == "CCG", "P",
                                                                                                                                        if_else(g == "TCT" | g == "TCC" | g =="TCA" | g == "TCG" | g == "AGT" | g == "AGC", "S",
                                                                                                                                                if_else(g == "ACT" | g == "ACC" | g =="ACA" | g == "ACG", "T",
                                                                                                                                                        if_else(g == "TGG", "W",
                                                                                                                                                                if_else(g == "TAT" | g == "TAC", "Y",
                                                                                                                                                                        if_else(g == "GTT" | g == "GTC" | g =="GTA" | g == "GTG", "V",
                                                                                                                                                                                if_else(g == "TAA" | g == "TGA" | g =="TAG", "stop", "unknown")))))))))))))))))))))) %>% 
  mutate(aa_c = if_else(c == "GCT" | c == "GCC" | c == "GCA" | c == "GCG", "A",  
                        if_else(c == "CGT" | c == "CGC" | c =="CGA" | c == "CGG" | c == "AGA" | c == "AGG", "R",
                                if_else(c == "AAT" | c == "AAC", "N", 
                                        if_else(c == "GAT" | c == "GAC", "D", 
                                                if_else(c == "TGT" | c == "TGC", "C",
                                                        if_else(c == "CAA" | c == "CAG", "Q", 
                                                                if_else(c == "GAA" | c == "GAG", "E",
                                                                        if_else(c == "GGT" | c == "GGC" | c =="GGA" | c == "GGG", "G",
                                                                                if_else(c == "CAT" | c == "CAC", "H",
                                                                                        if_else(c == "ATT" | c == "ATC" | c =="ATA", "I",
                                                                                                if_else(c == "CTT" | c == "CTC" | c =="CTA" | c == "CTG" | c == "TTA" | c == "TTG", "L",
                                                                                                        if_else(c == "AAA" | c == "AAG", "K",
                                                                                                                if_else(c == "ATG", "start",
                                                                                                                        if_else(c == "TTT" | c == "TTC", "F",
                                                                                                                                if_else(c == "CCT" | c == "CCC" | c =="CCA" | c == "CCG", "P",
                                                                                                                                        if_else(c == "TCT" | c == "TCC" | c =="TCA" | c == "TCG" | c == "AGT" | c == "AGC", "S",
                                                                                                                                                if_else(c == "ACT" | c == "ACC" | c =="ACA" | c == "ACG", "T",
                                                                                                                                                        if_else(c == "TGG", "W",
                                                                                                                                                                if_else(c == "TAT" | c == "TAC", "Y",
                                                                                                                                                                        if_else(c == "GTT" | c == "GTC" | c =="GTA" | c == "GTG", "V",
                                                                                                                                                                                if_else(c == "TAA" | c == "TGA" | c =="TAG", "stop", "unknown")))))))))))))))))))))) %>% 
  mutate(aa_t = if_else(t == "GCT" | t == "GCC" | t == "GCA" | t == "GCG", "A",  
                        if_else(t == "CGT" | t == "CGC" | t =="CGA" |t == "CGG" | t == "AGA" | t == "AGG", "R",
                                if_else(t == "AAT" | t == "AAC", "N", 
                                        if_else(t == "GAT" | t == "GAC", "D", 
                                                if_else(t == "TGT" | t == "TGC", "C",
                                                        if_else(t == "CAA" | t == "CAG", "Q", 
                                                                if_else(t == "GAA" | t == "GAG", "E",
                                                                        if_else(t == "GGT" | t == "GGC" | t =="GGA" | t == "GGG", "G",
                                                                                if_else(t == "CAT" | t == "CAC", "H",
                                                                                        if_else(t == "ATT" | t == "ATC" | t =="ATA", "I",
                                                                                                if_else(t == "CTT" | t == "CTC" | t =="CTA" | t == "CTG" | t == "TTA" | t == "TTG", "L",
                                                                                                        if_else(t == "AAA" | t == "AAG", "K",
                                                                                                                if_else(t == "ATG", "start",
                                                                                                                        if_else(t == "TTT" | t == "TTC", "F",
                                                                                                                                if_else(t == "CCT" | t == "CCC" | t =="CCA" | t == "CCG", "P",
                                                                                                                                        if_else(t == "TCT" | t == "TCC" | t =="TCA" | t == "TCG" | t == "AGT" | t == "AGC", "S",
                                                                                                                                                if_else(t == "ACT" | t == "ACC" | t =="ACA" | t == "ACG", "T",
                                                                                                                                                        if_else(t == "TGG", "W",
                                                                                                                                                                if_else(t == "TAT" | t == "TAC", "Y",
                                                                                                                                                                        if_else(t == "GTT" | t == "GTC" | t =="GTA" | t == "GTG", "V",
                                                                                                                                                                                if_else(t == "TAA" | t == "TGA" | t =="TAG", "stop", "unknown"))))))))))))))))))))))

codon_changes <- codon_changes %>% 
  mutate(degeneracy = if_else(aa_a == aa_g & aa_a == aa_c & aa_a == aa_t & aa_g == aa_c & aa_g == aa_t & aa_c == aa_t, "4-fold",
                              if_else(aa_a != aa_g & aa_a != aa_c & aa_a != aa_t & aa_g != aa_c & aa_g != aa_t & aa_c != aa_t, "0-fold",
                                      if_else(aa_a == aa_g & aa_a != aa_c & aa_a != aa_t & aa_g != aa_c & aa_g != aa_t & aa_c != aa_t, "2-fold",
                                              if_else(aa_a != aa_g & aa_a == aa_c & aa_a != aa_t & aa_g != aa_c & aa_g != aa_t & aa_c != aa_t, "2-fold",
                                                      if_else(aa_a != aa_g & aa_a != aa_c & aa_a == aa_t & aa_g != aa_c & aa_g != aa_t & aa_c != aa_t, "2-fold",
                                                              if_else(aa_a != aa_g & aa_a != aa_c & aa_a != aa_t & aa_g == aa_c & aa_g != aa_t & aa_c != aa_t, "2-fold",
                                                                      if_else(aa_a != aa_g & aa_a != aa_c & aa_a != aa_t & aa_g != aa_c & aa_g == aa_t & aa_c != aa_t, "2-fold",
                                                                              if_else(aa_a != aa_g & aa_a != aa_c & aa_a != aa_t & aa_g != aa_c & aa_g != aa_t & aa_c == aa_t, "2-fold", "unknown")))))))))
  

four_fold_sites <- codon_changes %>% 
   filter(degeneracy == "4-fold") %>% 
   select(SNP)

zero_fold_sites <- codon_changes %>% 
   filter(degeneracy == "0-fold") %>% 
   select(SNP)

write.table(four_fold_sites, "4-fold_variant_SNPs.txt", quote = F, col.names = F, row.names = F) 
write.table(zero_fold_sites, "0-fold_variant_SNPs.txt", quote = F, col.names = F, row.names = F) 
