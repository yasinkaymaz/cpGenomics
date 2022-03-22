
library(coRdon)

library(Biostrings)



# dnaLD94 <- readSet(file="~/Downloads/LD94.fasta")
# LD94 <- codonTable(dnaLD94)
# dnaHD59 <- readSet(file="~/Downloads/HD59.fasta")
# HD59 <- codonTable(dnaHD59)
# 
# cc <- codonCounts(HD59)
# cc
# l <- getlen(HD59)#sequence lengths
# head(l)
# length(HD59) #number of sequences
# 
# ko <- getKO(HD59)
# head(ko)

dnaErvo <- readSet(file="/Users/yasinkaymaz/Documents/ChlP_GS/Nurbanu/Interspecies_Cp_Phylo/dNdS-2/Lens_ervoides.CDS.fasta")
ervo <- codonTable(dnaErvo)
cce <- codonCounts(ervo)
dnaCuli <- readSet(file="/Users/yasinkaymaz/Documents/ChlP_GS/Nurbanu/Interspecies_Cp_Phylo/dNdS-2/Lens_culinaris.CDS.fasta")
culi <- codonTable(dnaCuli)
ccc <- codonCounts(culi)


ribp <- c("rpl14","rpl16","rpl2","rpl20","rpl23","rpl32","rpl33","rpl36","rps11","rps12","rps14","rps15","rps19","rps2","rps3","rps4","rps7","rps8")

milc <- MILC(ervo, subsets = list(Ribs = c(names(dnaErvo) %in% ribp)))
head(milc)

xlab <- "MILC distance from sample centroid"
ylab <- "MILC distance from ribosomal genes"

milc_ervo <- MILC(ervo, subsets = list(Ribs = c(names(dnaErvo) %in% ribp)))
Bplot(x = "Ribs", y = "self", data = milc_ervo, size = 8) +
  labs(x = xlab, y = ylab)


#Predict gene expression relative to ribosomal genes
melp_ervo <- MELP(ervo, subsets = list(Ribs = c(names(dnaErvo) %in% ribp)))
melp_culi <- MELP(culi, subsets = list(Ribs = c(names(dnaCuli) %in% ribp)))

data.frame(names(dnaErvo), melp_ervo, melp_culi) %>% View

GENETIC_CODE
length(GENETIC_CODE)

gc <- data.frame(Codon = names(GENETIC_CODE), AminoAcid = GENETIC_CODE)
cce.df <- data.frame(Codon = names(colSums(cce)), CU_ervoides=colSums(cce))

ccc.df <- data.frame(Codon = names(colSums(ccc)), CuliCU=colSums(ccc))


ervo <- codonTable(dnaErvo)
cce <- codonCounts(ervo)
cce <- as.data.frame(cce)
rownames(cce) <- names(dnaErvo)
cce <- as.data.frame(t(cce))
cce$Cod <- rownames(cce)

gc %>% 
  full_join(cce.df, by="Codon") %>% 
  as_tibble() %>%
  group_by(AminoAcid) %>%
  mutate(RSCU_ervoides = CU_ervoides*n()/sum(CU_ervoides)) %>% writexl::write_xlsx(path="RSCU_values.xlsx")
#  full_join(cce, by="Cod") %>% View()

gc %>% 
  full_join(cce.df, by="Codon") %>% 
  as_tibble() %>%
  group_by(AminoAcid) %>%
  mutate(RSCU_ervoides = CU_ervoides*n()/sum(CU_ervoides)) %>% 
  filter(RSCU_ervoides == max(RSCU_ervoides)) %>% select(Codon) %>% c() -> preferred.codons

gc %>% 
  full_join(cce.df, by="Codon") %>% 
  as_tibble() %>%
  group_by(AminoAcid) %>%
  mutate(RSCU_ervoides = CU_ervoides*n()/sum(CU_ervoides)) %>% 
  filter(RSCU_ervoides >= 1.0) %>% select(Codon) %>% c() -> preferred.codons

#Figure 5A
gc %>% 
  full_join(cce.df, by="Cod") %>% 
  as_tibble() %>%
  group_by(GenCode) %>%
  mutate(RSCU_evoides = ErvoCU*n()/sum(ErvoCU)) %>% 
  full_join(cce, by="Cod") %>%
  ggplot(aes(x=Cod, y=RSCU_evoides, fill=GenCode))+
  geom_bar(stat="identity",width = 0.9)+
  facet_grid(cols = vars(GenCode), scales = "free_x", space = "free_x")+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1,vjust = 0.5, colour = "black"),
        axis.text.y = element_text( hjust=1, colour = "black"))+ 
  theme(legend.position = "none", axis.title.x = element_blank())+
  ylab("RSCU values")

gc %>% 
  full_join(cce.df, by="Codon") %>% 
  select(-Codon) %>% 
  group_by(AminoAcid) %>% 
  summarise(AAsum=sum(CU_ervoides)) %>% 
  writexl::write_xlsx(path="Encoded_amino_acid_distribution.xlsx")


melp_ervo <- MELP(ervo, subsets = list(MELP = c(names(dnaErvo) %in% ribp)))
melp_ervo

ervo <- codonTable(dnaErvo)
cce <- codonCounts(ervo)
cce <- as.data.frame(cce)
rownames(cce) <- names(dnaErvo)
cce <- as.data.frame(t(cce))
#cce$Cod <- rownames(cce)

head(cce)
colSums(cce[preferred.codons$Codon,])/colSums(cce)

melp_ervo <- data.frame(row.names =names(dnaErvo), MELP=melp_ervo, Pref.CU=colSums(cce[preferred.codons$Codon,])/colSums(cce))
melp_ervo.xy <- melp_ervo
colnames(melp_ervo.xy) <- c("y", "x")

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn(melp_ervo.xy)

# melp_ervo %>% 
#   mutate(Gene=rownames(melp_ervo)) %>% 
#   mutate(Type=ifelse(Gene %in% ribp,"Ribosomal","Other")) %>%
#   ggplot( aes(MELP, 100*Pref.CU))+
#   geom_point(aes(color=Type),size=3)+ylim(50,100)+
#   #geom_text(aes(label=ifelse(MELP>1.0 & Type == "Other", as.character(Gene), ''),hjust=0,vjust=0,color="blue"))
#   geom_text(aes(label=ifelse(MELP>1.0 , ifelse(Type == "Other", as.character(Gene), ''), '')),hjust=0,vjust=0,color="blue")+
#   theme_bw()+
#   theme(text = element_text(size=15),
#         axis.text.x = element_text(angle=0,size = 15, hjust=.5,vjust = 0.5, colour = "black"),
#         axis.text.y = element_text( hjust=1,size = 15, colour = "black"))+ 
#   scale_colour_discrete(name="Gene type")+
#   ylab("Preferred codon usage (%)")+xlab("MELP (MILC-based Expression Level Predictor)")+
#   geom_smooth(method='lm')+ 
#   geom_text(x = 1.0, y = 95, label = lm_eqn(melp_ervo.xy), parse = TRUE)

# Figure 5B
melp_ervo %>% 
  mutate(Gene=rownames(melp_ervo)) %>% 
  mutate(Type=ifelse(Gene %in% ribp,"Ribosomal","Other")) %>%
  ggplot( aes( 100*Pref.CU, MELP))+
  geom_point(aes(color=Type),size=3)+xlim(50,90)+
  #geom_text(aes(label=ifelse(MELP>1.0 & Type == "Other", as.character(Gene), ''),hjust=0,vjust=0,color="blue"))
  geom_text(aes(label=ifelse(MELP>1.0 , ifelse(Type == "Other", as.character(Gene), ''), '')),hjust=0,vjust=0,color="blue")+
  theme_bw()+
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0,size = 15, hjust=.5,vjust = 0.5, colour = "black"),
        axis.text.y = element_text( hjust=1,size = 15, colour = "black"))+ 
  scale_colour_discrete(name="Gene type")+
  xlab("Preferred codon usage (%)")+ylab("MELP (MILC-based Expression Level Predictor)")+
  geom_smooth(method='lm', se=FALSE)+ 
  geom_text(y = .8, x = 60, label = lm_eqn(melp_ervo.xy), parse = TRUE)
