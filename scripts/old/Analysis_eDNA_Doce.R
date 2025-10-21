# 0- Load libraries ----
{
  library(pheatmap)
  library(Biostrings)
  library(cowplot)
  library(stringr)
  library(DECIPHER)
  library(ggdendro)
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(writexl)
  library(future)
  library(ggh4x)
  library(ggpubr)
  library(Matrix)
  library(phyloseq)
  library(ShortRead)
  library(stringr)
  library(readxl)
  library(vegan)
  library(phyloseq)
}

# Arrumar a tabela de entrada oficial ----

FINAL_TBL_igor$`Curated ID` %>% unique()


FINAL_TBL_igor <- read_xlsx("/home/igorhan/projetos/Doce12Sdb/12Sdb/eDNA_Doce/Resultados/ASV_DC_tbl-mai24/final_asv_tbl_abd_curated.xlsx",
                       col_names = TRUE) %>% 
  tibble()


FINAL_TBL_igor <-
  FINAL_TBL_igor %>% 
  mutate("Sample new name" = case_match(`Metadata 8`,
                                        "DC1 - Ponte Nova" ~ "DC1 - Ponte Nova",
                                        "DC2 - Jequitibá" ~ "DC4 - Jequitibá",
                                        "DC3 - Matipó" ~ "DC3 - Matipó",
                                        "DC4 - Coronel Fabriciano" ~ "DC2 - Coronel Fabriciano",
                                        "DC5 - Naque SA" ~ "DC5 - Naque SA",
                                        "DC6 - Naque DC" ~ "DC8 - Naque DC",
                                        "DC7 - Baguari" ~ "DC7 - Baguari",
                                        "DC8 - Periquito" ~ "DC6 - Periquito",
                                        "DC9 - Pedra Corrida" ~ "DC9 - Pedra Corrida",
                                        )) %>% 
  mutate("Sample name2" = case_match(`Sample name`,
                                     "DC-01" ~ "DC-01",
                                     "DC-02" ~ "DC-04",
                                     "DC-03" ~ "DC-03",
                                     "DC-04" ~ "DC-02",
                                     "DC-05" ~ "DC-05",
                                     "DC-06" ~ "DC-08",
                                     "DC-07" ~ "DC-07",
                                     "DC-08" ~ "DC-06",
                                     "DC-09" ~ "DC-09"
                                     )) %>% 
  dplyr::rename("Sample name - Original" = "Sample name") %>% 
  dplyr::rename("Sample name" = "Sample name2") %>% 
  dplyr::rename("Metadata 8 - Original" = "Metadata 8") %>% 
  dplyr::rename("Metadata 8" = "Sample new name") %>%
  mutate("Metadata 9" = case_when(`Metadata 8` %in% c("DC6 - Periquito",
                                                      "DC8 - Naque DC",
                                                      "DC7 - Baguari",
                                                      "DC9 - Pedra Corrida"  ) ~  "Mainstream",
                                  TRUE ~ "Tributaries")) %>% 
  mutate("Metadata 10" = case_when(`Metadata 8` %in% c("DC5 - Naque SA",
                                                       "DC8 - Naque DC",
                                                       "DC6 - Periquito",
                                                       "DC9 - Pedra Corrida") ~  "Lentic",
                                   TRUE ~ "Lotic")) %>% 
  mutate("Metadata 11" = case_match(`Metadata 8`,
                                    "DC1 - Ponte Nova" ~ "City",
                                    "DC2 - Coronel Fabriciano" ~ "City",
                                    "DC3 - Matipó" ~ "River",
                                    "DC4 - Jequitibá" ~ "River",
                                    "DC5 - Naque SA" ~ "Affected",
                                    "DC6 - Periquito" ~ "Affected")) %>% 
  
  # renomear algumas espécies
  mutate("Curated ID" = case_when(`Curated ID` %in% c("Brycon sp. A") ~ "Brycon dulcis",
                                  `Curated ID` %in% c( "Tetragonopterus rivularis") ~ "Psalidodon cf. rivularis",
                                  `Curated ID` %in% c( "Pseudoplatystoma reticulatum") ~ "Pseudoplatystoma hybrid",
                                  TRUE ~ `Curated ID`)) %>% 
    
  mutate("Genus (BLASTn)" = case_when(`Curated ID` %in% c("Trichomycterus sp.") ~ "Trichomycterus",
                                      TRUE ~ `Genus (BLASTn)`)) %>% 
    
    mutate("Subfamily (BLASTn)" = case_when(`Curated ID` %in% c("Trichomycterus sp.") ~ "Trichomycterinae",
                                      TRUE ~ `Subfamily (BLASTn)`)) %>%
    
    mutate("Family (BLASTn)" = case_when(`Curated ID` %in% c("Trichomycterus sp.") ~ "Trichomycteridae",
                                      TRUE ~ `Family (BLASTn)`)) %>% 
    
    mutate("Suborder (BLASTn)" = case_when(`Curated ID` %in% c("Trichomycterus sp.") ~ "Loricarioidei",
                                      TRUE ~ `Suborder (BLASTn)`)) %>% 
    
    mutate("Order (BLASTn)" = case_when(`Curated ID` %in% c("Trichomycterus sp.") ~ "Siluriformes",
                                      TRUE ~ `Order (BLASTn)`)) %>% 
    
  ################# filtrar especies ################################################################################
# filter(`Final ID (BLASTn)` %in% c("Homo sapiens","Sus scrofa","Bos taurus"))
# filter(`Metadata 11` %in% c("City", "No Affected", "Affected")) %>% 

# ABD limpa por amostra
# group_by(Unique_File_name,Primer,`Expected length`,`ID status`, `Possible Metazoa`,`Read origin`
  group_by(
  # Sample,
  `Sample name`,
  `Primer expected length`,`ID status`, `Possible Metazoa`,`Read origin`,To_remove
  # ,`Possible contamination`
) %>%
  # group_by(Sample,`Expected length`,`ID status`, `Possible Metazoa`,`Read origin`,`Possible contamination`) %>%  #WWF
  mutate("Total clean sample abd."  = 0,
         "Total clean sample abd." = case_when((`ID status` %in% c("True detection") & 
                                                  `Primer expected length` %in% c("in range") & 
                                                  !To_remove %in% c("remove") &
                                                  # `Possible contamination` %in% c("True detection") &
                                                  `Possible Metazoa` == TRUE) ~ sum(`ASV absolute abundance`)
                                               # (`ID status` %in% "not IDed") ~ 0,
                                               # (`Expected length` %in% "out of range") ~ 0,
                                               # (`Possible contamination` %in% "Possible contamination") ~ 0,
                                               # (`Possible Metazoa` %in% FALSE) ~ 0
                                               # TRUE ~ 0
         )) %>% 
  ungroup() %>% 
  
  # mutate("Clean relative abd. on sample" =  case_when(`ID status` == "not IDed"  ~ 0,
  #                                                     `Expected length` == "out of range" ~ 0,
  #                                                     `ID status` == "IDed" &
  #                                                       `Expected length` == "in range" &
  #                                                       `Possible contamination` == "True detection" &
  #                                                       `Possible Metazoa` == TRUE ~ (`ASV absolute abundance`/`Total clean sample abd.`))) %>%
  mutate("Clean relative abd. on sample" =  (`ASV absolute abundance`/`Total clean sample abd.`)) %>% 
  # mutate("Unique_File_name" = File_name) %>% 
  relocate("Unique_File_name","Sample name","Primer","Curated ID","Primer expected length","ID status", "Possible Metazoa","Read origin",
           # "Possible contamination",
           "Clean relative abd. on sample","Total clean sample abd.") %>% 
  mutate("BLASTn pseudo-score" = `1_indentity` *`1_qcovhsp` /100) %>% 
  mutate("Identification" =  case_when(`BLASTn pseudo-score` >= 98 ~ `Final ID (BLASTn)`,
                                       `BLASTn pseudo-score` >= 95 & `BLASTn pseudo-score` < 98 ~ `Genus (BLASTn)`,
                                       `BLASTn pseudo-score` >= 90 & `BLASTn pseudo-score` < 95 ~ `Family (BLASTn)`,
                                       `BLASTn pseudo-score` >= 80 & `BLASTn pseudo-score` < 90 ~ `Order (BLASTn)`,
                                       `BLASTn pseudo-score` >= 60 & `BLASTn pseudo-score` < 80 ~ `Class (BLASTn)`),
         "Identification Max. taxonomy" =   case_when(`BLASTn pseudo-score` >= 98 ~ "Species",
                                                      `BLASTn pseudo-score` >= 95 & `BLASTn pseudo-score` < 98 ~ "Genus",
                                                      `BLASTn pseudo-score` >= 90 & `BLASTn pseudo-score` < 95 ~ "Family",
                                                      `BLASTn pseudo-score` >= 80 & `BLASTn pseudo-score` < 90 ~ "Order",
                                                      `BLASTn pseudo-score` >= 60 & `BLASTn pseudo-score` < 80 ~ "Class")) %>% 
  relocate("Identification","Identification Max. taxonomy","BLASTn pseudo-score") %>% 
  filter(!To_remove %in% c("remove"))

# add new metadata
metadata12_igor <- read_xlsx("/home/igorhan/projetos/Doce12Sdb/12Sdb/eDNA_Doce/data/metadados/metadata12.xlsx",
                        col_names = T)

metadata12_igor <- metadata12_igor %>% 
  rename("Metadata 12" = "Native")

FINAL_TBL_igor <- left_join(x = FINAL_TBL_igor,
                       y = metadata12_igor,
                       by = "Curated ID")

# converter a tabela final de resultados para o formato amplo, compatível com o NMDS ----  

FINAL_TBL_igor_IDs <-
  FINAL_TBL_igor %>% 
  
  # 
  # filter(`Metadata 8` %in% c("DC1 - Ponte Nova",         # REMOVER O FILTRO!!!!!!!!
  #                            "DC2 - Coronel Fabriciano",
  #                            "DC3 - Matipó",
  #                            "DC4 - Jequitibá")) %>%


  mutate("Sample name" = `Metadata 8`) %>% 
  mutate("agrupador" = `Metadata 9`) %>% 
  arrange(`Sample name`) %>% 
  ###_______________remover amostras que atrapalham a visualização_______________###
  
  filter(!To_remove %in% "remove") %>%
  
  dplyr::select(-c(      
    "Kingdom (BLASTn)",
    "Phylum (BLASTn)",
    "Subphylum (BLASTn)",
    "Class (BLASTn)",
    "Subclass (BLASTn)",
    "Order (BLASTn)",
    "Suborder (BLASTn)",
    "Family (BLASTn)",
    "Subfamily (BLASTn)",
    "Genus (BLASTn)",
    "Final ID (BLASTn)",
    "Relative abundance on sample",         
    "ASV Size (pb)",
    "Primer expected length",
    "Type",
    ###____________________      infos das IDs do BLASTn    ______________________
    "1_subject header","1_subject","1_indentity","1_qcovhsp","1_length","1_mismatches","1_gaps",
    "1_query start","1_query end","1_subject start","1_subject end","1_e-value","1_bitscore", 
    # "1_staxid",
    "2_subject header","2_subject","2_indentity","2_qcovhsp","2_length","2_mismatches","2_gaps",
    "2_query start","2_query end","2_subject start","2_subject end","2_e-value","2_bitscore", 
    # "2_staxid",
    "3_subject header","3_subject","3_indentity","3_qcovhsp","3_length","3_mismatches","3_gaps",
    "3_query start","3_query end","3_subject start","3_subject end","3_e-value","3_bitscore", 
    # "3_staxid",
    "max_tax",
    ###____________________      Sequencia da ASV    _____________________________
    "ASV (Sequence)", "ASV header",
    
    ###________________ outras infos tecnicas _________________
    "Read origin",
    "Relative abundance to all samples",
    "Sample total abundance",
    "Metadata 1",
    "Metadata 2",
    "Metadata 3",
    "Metadata 4",
    "Metadata 5",
    "Metadata 6",
    "Metadata 7",
    "Metadata 8 - Original",
    "Possible Metazoa",
    "ID status", 
    "PCR control",
    "Prop. to PCR control",
    "Prop. to Ext control",
    "Prop. to Filt control",
    ends_with("(DADA2)"),
    ends_with("(DADA2 bootstrap)")
  )) %>% 
  # colnames()
  pivot_wider(                                #pivotando as IDs de linhas pra colunas
    id_cols = c("Sample name",
                "agrupador",
                "Metadata 8",
                "Metadata 9",
                "Metadata 10",
                "Metadata 11"
                ),
    values_from ="Clean relative abd. on sample",            #utilizando abundancias Identificadas
    # values_from ="ASV rel. abd. in Sampling Unit by marker",
    values_fn = sum_uniq,
    # names_from = Identification,
    names_from = `Curated ID`,
    names_sort = TRUE,
    names_prefix = "ID_"
  ) %>% 
  relocate(c("Sample name",
             starts_with("ID_")
  )) %>%  
  mutate(dplyr::across(starts_with("ID_") , replace_na, replace = 0)) 


#### Check if you have only one row per sample name ----
FINAL_TBL_igor_IDs$`Sample name` %>% table()


FINAL_TBL_igor_IDs %>% colnames() %>% duplicated()

#### Check if your samples IDs sum to 1 ----

FINAL_TBL_igor_IDs %>% 
  select(c("Sample name",starts_with("ID_"))) %>% 
  # rowwise() %>%
  mutate("SOMA" = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% 
  relocate("SOMA")


# tem que dar 1 pra tudo
FINAL_TBL_igor_IDs %>%
  select(starts_with(match = "ID")) %>%
  rowSums() %>% plot()
# 
# # 
FINAL_TBL_igor_IDs %>%
  select(starts_with(match = "ID")) %>%
  colSums() %>% plot()


#NMDS ----

#refs
# https://rpubs.com/CPEL/NMDS

#1- prepare data for entry in vegan ----
# colnames(FINAL_TBL_igor_IDs)

all_IDs_NMDS_tbl_igor <- FINAL_TBL_igor_IDs %>% 
  mutate("Sample number" = 0) %>%
  relocate(
    `Sample number`)

#2- associate sample numbers to sample names ----
for (sample in 1:nrow(all_IDs_NMDS_tbl_igor)) {
  
  all_IDs_NMDS_tbl_igor$`Sample number`[sample] <- sample
  
}

#ordenar df usada no NMDS ----
all_IDs_NMDS_df_igor <- all_IDs_NMDS_tbl_igor %>% 
  # select(base::sort(colnames(.))) %>%
  relocate(c("Sample number",
             "Sample name",
             "agrupador",
             #           
             #           "metadata_1", 
             #           "metadata_2", 
             #           "metadata_3", 
             #           "metadata_4", 
             #           "metadata_5", 
             #           "metadata_6",
             #           "metadata_7",
                       "Metadata 8",
                       "Metadata 9",
                       "Metadata 10",
                       "Metadata 11"
  )) %>%
  as.data.frame() 

#4- name rows as Sample numbers and remove column ----
row.names(all_IDs_NMDS_df_igor) <- all_IDs_NMDS_df_igor$`Sample number`

all_IDs_NMDS_df_igor %>% dim()


colnames(all_IDs_NMDS_df_igor) %>% head()




# correct species names to avoid problems in ploting

colnames(all_IDs_NMDS_df_igor)[7:ncol(all_IDs_NMDS_df_igor)] <- colnames(all_IDs_NMDS_df_igor)[7:ncol(all_IDs_NMDS_df_igor)] %>%
  str_replace_all(pattern = " ",replacement = "_") %>% 
  str_replace_all(pattern = "\\.",replacement = "") %>% 
  str_replace_all(pattern = "\\(",replacement = "") %>% 
  str_replace_all(pattern = "\\)",replacement = "")





all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df_igor[,8:ncol(all_IDs_NMDS_df_igor)],
                                 comm = all_IDs_NMDS_df_igor[,8:ncol(all_IDs_NMDS_df_igor)],
                                 # distance = "bray"
                                 distance = "jaccard"
)

plot(all_ps_vegan_ord_meta)

dim(all_IDs_NMDS_df_igor)


# ----
#retirado daqui!!!!!!!!!!!!!!!!!!!!       https://www.rpubs.com/RGrieger/545184 

meta.envfit <- envfit(all_ps_vegan_ord_meta, 
                      all_IDs_NMDS_df_igor[,c("agrupador","Sample name","Metadata 8", "Metadata 9","Metadata 10","Metadata_11")], 
                      permutations = 999,
                      na.rm=TRUE) # this fits environmental vectors


# esse é o passo que mais demora
meta.spp.fit <- envfit(all_ps_vegan_ord_meta, all_IDs_NMDS_df_igor[,8:ncol(all_IDs_NMDS_df_igor)], permutations = 999) # this fits species vectors


site.scrs <- as.data.frame(scores_igor(all_ps_vegan_ord_meta, display = "sites")) %>% 
  mutate("Sample number" = as.double(row.names(.))) %>% 
  # left_join(y = all_IDs_NMDS_df_igor[,c("Sample",
  left_join(y = all_IDs_NMDS_df_igor[,c("Sample name",
                                   "Sample number",
                                   "agrupador",
                                   "Metadata 8",
                                   "Metadata 9",
                                   "Metadata 10",
                                   "Metadata_11"
                                   )],
            by = "Sample number") 

site.scrs

# PERMANOVA ----
adonis2(all_IDs_NMDS_df_igor[,8:ncol(all_IDs_NMDS_df_igor)] ~ site.scrs$agrupador, permutations = 10000)

# determinar centroides ----
scrs <-
  scores_igor(all_ps_vegan_ord_meta, display = "sites")

cent <-
  aggregate(scrs~`agrupador`,data = site.scrs, FUN = "mean")

#get species pvalues ----
sps_pvals <- tibble("IDs" = names(meta.spp.fit$vectors$pvals),
                    "p-value" = meta.spp.fit$vectors$pvals)


spp.scrs <- as.data.frame(scores_igor(meta.spp.fit, display = "vectors")) %>% 
  mutate("IDs" = rownames(.)) %>% 
  left_join(y = sps_pvals, by = "IDs")

sig.spp.scrs <- spp.scrs %>% 
  filter(`p-value` <=0.05)                   # selecionar para mostrar apenas sps com pval significativo


#calculate ellipses ----
NMDS <- data.frame("MDS1" = all_ps_vegan_ord_meta$points[,1], 
                   "MDS2" = all_ps_vegan_ord_meta$points[,2],
                   "agrupador"= as.factor(all_IDs_NMDS_df_igor$`agrupador`),
                   check.names = FALSE)


NMDS.mean <- aggregate(NMDS[,1:2],list(group=NMDS$`agrupador`),"mean")


# funçaõ do vegan de calcular ellipses
veganCovEllipse<-function (cov, 
                           center = c(0, 0), 
                           scale = 1, 
                           npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# dev.off()
# dev.off()
plot(all_ps_vegan_ord_meta)

ord <- ordiellipse(ord = all_ps_vegan_ord_meta, 
                   groups = all_IDs_NMDS_df_igor$`agrupador`,
                   # groups = all_IDs_NMDS_df_igor$StatEstado,
                   display = "sites",
                   kind = "ehull", conf = 0.95, label = T)


# fit <- envfit(all_ps_vegan_ord_meta~Al,varechem,perm=999,display="lc")

# vegan::ordiarrows(ord = all_ps_vegan_ord_meta,
#                   groups = 
#                     )

df_ell <- data.frame()

for(g in levels(NMDS$`agrupador`)){
  
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$`agrupador`== g,],
                                                   veganCovEllipse(cov = ord[[g]]$cov,
                                                                   center = ord[[g]]$center,
                                                                   scale = ord[[g]]$scale))),
                                "agrupador"=g))
  # StatEstado=g))
}




# fazer o grafico ----

# cores_igor ----

# devtools::install_github("karthik/wesanderson")
# devtools::install_github("gadenbuie/ggpomological")
# 
# wesanderson::wes_palette("GrandBudapest1")
# 
# scales::show_col(ggthemes::calc_pal()(12))
# scales::show_col(wesanderson::wes_palette("GrandBudapest1"))
# 
# "#004586"
# "#ff420e"
# "#ffd320"
# "#7e0021"
# "#83caff"
# "#314004"
# "#aecf00"
# "#4b1f6f"
# "#ff950e" 
# "#c5000b" 
# "#0084d1"
# "#A4A200FF"
# 
# cores_igor <- c( "Lentic" = "#ff950e",
#             "Lotic" = "#004586"
# )
# 
cores_igor <- c( "Tributaries" = "#004586",
            "Mainstream" = "#ff950e"
            )

# cores_igor <- c( "City" = "#004586",
#             "No Affected" = "#aecf00",
#             "Affected" = "#ff950e"
#             )



# Gráfico ----
manual_NMDS_plot <-
  site.scrs %>% 
  mutate(`Sample name` = case_match(`Sample name`,
                                    "DC1 - Ponte Nova" ~ "P1",
                                    "DC2 - Coronel Fabriciano" ~ "P2",
                                    "DC3 - Matipó" ~ "P3",
                                    "DC4 - Jequitibá" ~ "P4",
                                    "DC5 - Naque SA" ~ "P5",
                                    "DC6 - Periquito" ~ "P6",
                                    "DC7 - Baguari" ~ "P7",
                                    "DC8 - Naque DC" ~ "P8",
                                    "DC9 - Pedra Corrida" ~ "P9")) %>% 
  ggplot(aes(x=NMDS1, 
             y=NMDS2)) +
  #elipses ####
ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                           data = df_ell,
                           aes(x = NMDS1,
                               y = NMDS2,
                               group=`agrupador`,
                               label=`agrupador`,
                               col =`agrupador`,
                               fill =`agrupador`
                               ),
                           alpha=0.15,
                           # n = 200,
                           linetype=2,
                           # expand = unit(5, "px"),
                           expand = 0,
                           label.fontsize = 18,
                           con.cap = 0.1
                           ) +
#hulls #########
# ggforce::geom_mark_hull(aes(fill=Metadata.2),
#                         concavity = 5,
#                         expand=0,
#                         radius=0,
#                         linetype=0,
#                         alpha=0.15
#                         )+
#vetores das IDs
geom_segment(data = sig.spp.scrs, aes(x = 0,
                                      xend=NMDS1,
                                      y=0,
                                      yend=NMDS2),
             arrow = arrow(length = unit(0.1, "cm")),
             colour = "grey10",
             alpha=0.1,
             lwd=0.3) + #add vector arrows of significant species
  # #nomes das IDs
  ggrepel::geom_text_repel(data = sig.spp.scrs,
                           aes(x=NMDS1, y=NMDS2, label = IDs),
                           size=3,
                           alpha= 0.75,
                           # cex = 5,
                           direction = "both",
                           segment.size = 0.25,
                           segment.alpha=0.1,
                           max.overlaps = 100) +
  #pontos amostrais
  geom_point(aes(x=NMDS1, 
                 y=NMDS2, 
                 fill = `agrupador`,
                 col = `agrupador`,
                 # label = `Unique ID`,
                 # alpha = 0.75,
                 group = `agrupador`,
                 # shape = `agrupador`
  ), 
  stroke = 0.5,
  alpha=0.75,
  # col="#656565",
  size = 5
  # size = `metadata_8`
  )+ 
  #nomes dos pontos amostrais
  geom_text(aes(label = `Sample name`),
            # geom_text(aes(label = interaction(`Sample name`, `metadata_1`,sep = "\n")),
            hjust=0.5, 
            vjust=2.75, 
            size=5) +
  #centroides ----
geom_point(data = cent,
           aes(x=NMDS1, 
               y=NMDS2, 
               # colour = Metadata.2,
               fill = `agrupador`
           ),
           size= 1,
           colour="#222222",
           alpha=0.75,
           shape = 23
)+
  coord_fixed()+
  scale_fill_manual(values = cores_igor) +
  scale_colour_manual(values = cores_igor) +
theme_light()+ 
  # labs(colour = "Intervenção", 
  #      shape = "Local") + 
  theme(legend.position = "right", 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
  # notação do valor do stress do NMDS
  labs(
    # title  = paste0("B"),
       subtitle = paste0(
         "\n","Stress: ",format(round(all_ps_vegan_ord_meta$stress,4)), "\n", "P-value: 0.0232")) +
  theme(plot.title = element_text(size = 20)) +
  theme(plot.subtitle = element_text(size = 16)) +
  theme(legend.title = element_text(size = 12)) +
  theme(legend.text =  element_text(size = 12)) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill=guide_legend(title="Sample"),
         col=guide_legend(title="Sample"),
         shape=guide_legend(title="Sample"))

manual_NMDS_plot



#get plot area dims ---- 
plot_dims <- DeLuciatoR::get_dims(manual_NMDS_plot,
                                  maxheight = 50,
                                  maxwidth = 50,
                                  units = "cm")


#save plot ----
ggsave(file = paste0("/home/igorhan/projetos/Doce12Sdb/12Sdb/eDNA_Doce/Resultados/NMDS_manuscript","-", Sys.Date(), ".pdf"),
       plot = manual_NMDS_plot,
       device = "pdf",
       units = "cm",
       width = (plot_dims$width * 0.5),
       height = (plot_dims$height * 0.5),
       dpi = 300, limitsize = FALSE)


 # TILE PLOT ----

final_asv_tbl_igor <- FINAL_TBL_igor


# 6- Arrumando os dados para tile plot ----

final_asv_tbl_igor$`Curated ID` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
final_asv_tbl_igor$`Sample name` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()

subset(final_asv_tbl_igor, `Class (BLASTn)` == "Actinopteri")$`Order (BLASTn)` %>%
  sort() %>% unique()

# organizar a ordem das ordens
fish_ordens_igor <- c("Characiformes",
                 "Cichliformes",
                 "Cyprinodontiformes",
                 "Gymnotiformes",
                 "order of Neopterygii",
                 "Siluriformes",
                 "Synbranchiformes")

# Definir quais serao as especies e ordenar

final_asv_tbl_igor$`Sample name` %>% unique()

spp_total_igor <- final_asv_tbl_igor %>%
  filter(`Class (BLASTn)` == "Actinopteri") %>%
  filter(`Curated ID` != "") %>%
  select(`Curated ID`) %>%
  unique() %>%
  arrange(`Curated ID`) %>%
  pull()
list(spp_total_igor)

# Definir quais serao as amostras e ordenar
{
  final_asv_tbl_igor$`Metadata 8` %>% unique()
  
  samples_igor <- c("DC1 - Ponte Nova",
               "DC2 - Jequitibá",
               "DC3 - Matipó",
               "DC4 - Coronel Fabriciano",
               "DC5 - Naque SA",
               "DC6 - Naque DC",
               "DC7 - Baguari",
               "DC8 - Periquito",
               "DC9 - Pedra Corrida")
}

## 6.1- tabela so peixes ----
{
  fish_ID_tbl_igor <- final_asv_tbl_igor %>% # ids apenas os peixes a nivel de spp
    filter(`Class (BLASTn)` == "Actinopteri") %>%
    arrange(`Curated ID`) 
  # %>%
  #   filter(`ID status` != "Possible contamination",
  #          `Status curated` != "Possible contamination",
  #          `Sample name` != "DC-Br")
  
  View(fish_ID_tbl_igor)
  
}

fish_ID_tbl_igor %>% 
  select(`Metadata 8`, `Metadata 9`, `Metadata 10`) %>% View()

#7- Tile Plot das amostras ----

asv_tile_plot_igor <-
  fish_ID_tbl_igor %>% 
  mutate(`Metadata 9` = factor(`Metadata 9`, levels = c("Tributaries", "Mainstream"))) %>% 
  # mutate(`Metadata 10` = factor(`Metadata 10`, levels = c("Lotic", "Lentic"))) %>% 
  filter(!To_remove %in% "remove") %>% 
  group_by(`Metadata 8`, `Metadata 9`,
           # `Metadata 10`,
           `Order (BLASTn)`, `Curated ID`) %>% 
  summarise(RRA = sum(`Relative abundance on sample`)) %>% 
  ungroup() %>% 
  # arrange(`Clean relative abd. on sample`) %>% 
  # filter(RRA >= 0.01) %>%
  # mutate(`Curated ID` = factor(`Curated ID`, levels = rev(fish_ID_tbl_igor))) %>%
  mutate(`Order (BLASTn)` = factor(`Order (BLASTn)`)) %>%
  group_by(`Curated ID`, `Metadata 8`, `Metadata 9`,
           # `Metadata 10`,
           `Order (BLASTn)`) %>%
  
  mutate(`Metadata 8` = case_match(`Metadata 8`,
                                    "DC1 - Ponte Nova" ~ "P1",
                                    "DC2 - Coronel Fabriciano" ~ "P2",
                                    "DC3 - Matipó" ~ "P3",
                                    "DC4 - Jequitibá" ~ "P4",
                                    "DC5 - Naque SA" ~ "P5",
                                    "DC6 - Periquito" ~ "P6",
                                    "DC7 - Baguari" ~ "P7",
                                    "DC8 - Naque DC" ~ "P8",
                                    "DC9 - Pedra Corrida" ~ "P9")) %>% 
  # filter(!`Metadata 8` %in% NA) %>% 
  # select(`Metadata 8`, `Metadata 9`, `Metadata 10`) %>% View()
  
  ggplot(aes(y = `Curated ID`,
             group = `Curated ID`,
             x = `Metadata 8`,
             fill = RRA)) +
  geom_tile() +
  # geom_text(aes(label= sprintf("%0.3f", round(RRA, digits = 4))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
  #           stat = "identity",
  #           colour = "white", size = 3) +
  facet_grid(cols = vars(`Metadata 9` # Facet por ponto
  ), 
  rows = vars(`Order (BLASTn)`),
  space = "free", 
  scales = "free",
  drop = TRUE
  ) +
  scale_fill_continuous(
    trans = "log10",
    breaks = c(0.000001, 0.0087, 0.87),  # Defina os pontos de quebra desejados
    labels = c("0.01", "0.1", "1"),  # Rótulos correspondentes
    type = "viridis"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey",
                                    size = 0.2,
                                    linetype = 1),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", size = rel(1.2)),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = rel(1.2)),
    plot.title = element_text(color = "black", size = rel(1.5)),
    plot.subtitle = element_text(color = "black", size = rel(1.2)),
    strip.background = element_rect(fill = "#e4e4e4"),
    strip.text = element_text(color = "black", size = rel(1.2))) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title.x = element_text(size = 27, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.position = "top"
  ) +
  labs(x = "Sample sites",
       y = "Taxa",
       fill ='Relative \nabundance (%)',
       title = "Total taxa detected"
       # ,
       # subtitle = "Comparing Lotic x Lentic sites"
       )

# ALPHA DIVERSITY ----

# 9- calculo de indices phyloseq ----

## 9.1- função necessária para calcular abundancias ----

# Função necessária para juntar as abds da tabela em IDs
sum_uniq <- function(vec) {
  if (is.character(vec)) {
    suniq <- BiocGenerics::unique(vec)
  }
  if (is.numeric(vec)) {
    suniq <- sum(vec)
  }
  return(suniq)
}

# índices por ID ----

# nesta parte vamos recalcular as diversidades usando as identificações unicas
# a nível de "espécie", garantindo a-diversidade correta

## 9.2a- create tax table ----

#create new tax table from BLASTn identifications

blast_tax_table_id_igor <- FINAL_TBL_igor %>% 
  #selecionando colunas a se usar
  dplyr::select(c("Curated ID",
                  "Genus (BLASTn)",
                  "Subfamily (BLASTn)",
                  "Family (BLASTn)",
                  "Suborder (BLASTn)",
                  "Order (BLASTn)",
                  "Subclass (BLASTn)",
                  "Class (BLASTn)",
                  "Phylum (BLASTn)",
                  "Subphylum (BLASTn)",
                  "Kingdom (BLASTn)")) %>% 
  unique() %>% 
  as.data.frame() %>% 
  `rownames<-`(.$`Curated ID`) %>% 
  dplyr::select(-c("Curated ID")) %>% 
  as.matrix()


blast_tax_table_id_igor %>% View()

## 9.3a- creating "OTU" table ----

FINAL_TBL_igor$`Sample name` %>% unique()

# o calculo da abundancia já foi feito na tabela final, se recalcular aqui novamente, será perdido informação!!!


FINAL_ID_table_id_igor <-
  FINAL_TBL_igor %>%
  # filter(`Clean relative abd. on sample` > 0) %>% 
    dplyr::mutate("Sample name" = str_replace_all(string = `Sample name`,
                                                  pattern = "-",
                                                  replacement = "")) %>%

  #normalizar abundâncias com apenas o que restou na tabela
  # group_by(
  #   `Sample name`,
  #   `Primer expected length`,
  #   `ID status`,
  #   `Possible Metazoa`,
  #   `Read origin`
  # ) %>%
    
  # mutate("Total clean sample abd."  = 0,
  #        "Total clean sample abd." = case_when((`ID status` %in% c("True detection") & 
  #                                                 `Primer expected length` %in% c("in range") & 
  #                                                 !To_remove %in% c("remove") &
  #                                                 `Possible Metazoa` == TRUE) ~ sum(`ASV absolute abundance`)
  #        )) %>% 
  # ungroup() %>% 
  # mutate("Clean relative abd. on sample" =  (`ASV absolute abundance`/`Sample total abundance`)) %>% 
  dplyr::select(c(      
    "Sample name",      # este é nosso identificador único  de cada ponto amostral campo
    "Clean relative abd. on sample", # original
    # "ASV absolute abundance",
    # "ASVs por ID",
    # "ASV (Sequence)",
    "Curated ID"
  )) %>% 
  pivot_wider(                                #pivotando as IDs de linhas pra colunas
    id_cols = c("Curated ID"),
    values_from ="Clean relative abd. on sample",            #utilizando abundancias Identificadas
    values_fn = sum_uniq,
    names_from = "Sample name",
    # names_prefix = "DC",
    # names_from = `OTU_ID`,
    # names_from = `Genus (BLASTn)`,
    names_sort = TRUE,
    # names_prefix = "DC"
  ) |>
  mutate(across(starts_with("DC") ,~ replace_na(.,replace = 0))) |> 
  # mutate_if(is.numeric,  ~ . * 100000) |>
  # mutate_if(is.numeric, round) |>
  mutate_if(is.numeric, jaccarize) |>
  as.data.frame() |> 
  column_to_rownames(var = "Curated ID") %>% 
    dplyr::filter(rowSums(.)!= 0) %>%   ############# Remover linhas cuja soma é zero
    dplyr::select(which(!colSums(., na.rm = TRUE) %in% 0))   ############# Remover colunas cuja soma é zero

dim(FINAL_ID_table_id_igor)

FINAL_ID_table_id_igor %>% rowSums() ==0
FINAL_ID_table_id_igor %>% colSums()


## 9.4a- creating sample metadata table ----
metadata_tbl_id_igor <- FINAL_TBL_igor %>% 
  # select(c("Sample name", starts_with("Metadata"))) %>% 
  select(c("Sample name", "Metadata 8", "Metadata 9")) %>% 
  mutate("Sample name" = str_replace_all(string = `Sample name`,pattern = "-",replacement = "")) %>%
  # dplyr::group_by(`Sample name`) %>% 
    unique() %>% 
  dplyr::rename("Sample" = "Sample name",
                # "Metadata1" = "Metadata 1",
                # "Metadata2" = "Metadata 2",
                # "Metadata3" = "Metadata 3",
                # "Metadata4" = "Metadata 4",
                # "Metadata5" = "Metadata 5",
                # "Metadata6" = "Metadata 6",
                # "Metadata7" = "Metadata 7",
                "Metadata8" = "Metadata 8",
                "Metadata9" = "Metadata 9"
                # "Metadata10" = "Metadata 10",
                # "Metadata11" = "Metadata 11"
                ) %>% 
  as.data.frame() %>% 
  `rownames<-`(.$`Sample`)


## 9.5a- create phyloseq object ----
FINAL_PS_id_igor <- phyloseq::phyloseq(otu_table(FINAL_ID_table_id_igor,taxa_are_rows = T),
                                  sample_data(metadata_tbl_id_igor),
                                  tax_table(blast_tax_table_id_igor))

## 9.6a- estimate_richness(FINAL_PS) ----
div_est_phylo_id_igor <-  estimate_richness(FINAL_PS_id_igor,
                                       split = T, 
                                       measures = c("Observed", "Chao1","Obs./Chao1",
                                                    "se.chao1","ACE", "se.ACE",
                                                    "Shannon","Simpson", "InvSimpson")) %>% 
  as_tibble(rownames = "Sample") %>% #ate aq
  mutate("Obs./Chao1" = Observed/Chao1) %>% 
  # relocate("EnvMat","Status","Estado") %>% 
  tidyr::pivot_longer(cols = c("Observed", "Chao1","Obs./Chao1", "se.chao1",
                               "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson"),
                      names_to = "Indice",
                      values_to = "Values") %>% 
  mutate(Indice = factor(Indice, levels = c("Observed", "Chao1","Obs./Chao1", "se.chao1",
                                            "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson"))) %>% 
  dplyr::left_join(y = metadata_tbl_id_igor,
                   by = "Sample")

# 9.7a- conferindo se as contagens estao certas ----

FINAL_TBL_igor %>% 
  # filter(`Clean relative abd. on sample` > 0) %>% 
  select(`Sample name`,`Curated ID`, `Metadata 8`) %>% 
  unique() %>% 
  pull(`Metadata 8`) %>% 
  table()


# 10a- Plot diversidade por ponto ----

# definir cores_igor por ponto

FINAL_TBL_igor$`Metadata 8` %>% unique()

# turbo
viridis::turbo(n=9)

# cores_igor <- c("P1" = "#30123BFF",
#            "P2" = "#466BE3FF",
#            "P3" = "#28BBECFF",
#            "P4" = "#31F299FF",
#            "P5" = "#A2FC3CFF",
#            "P6" = "#EDD03AFF",
#            "P7" = "#FB8022FF",
#            "P8" = "#D23105FF",
#            "P9" = "#7A0403FF")

cores_igor <- c("P1" = "#4d7f17",
           "P2" = "#6bb120",
           "P3" = "#8ae429",
           "P4" = "#9afe2e",
           "P5" = "#aefe57",
           "P6" = "#f8ed62",
           "P7" = "#e9d700",
           "P8" = "#e9d700",
           "P9" = "#dab600")

scales::show_col(cores_igor)

## plot richness ----

# div_est_phylo_plot_id_igor_obs <-
div_est_phylo_id_igor %>% 
  filter(Indice %in% c(
    "Observed"
    # ,
    # "Chao1",
    # "Shannon",
    # "Simpson"
  )) %>% 
  mutate("Metadata8" = case_match(Metadata8,
                                  "DC1 - Ponte Nova" ~ "P1",
                                  "DC2 - Coronel Fabriciano" ~ "P2",
                                  "DC3 - Matipó" ~ "P3",
                                  "DC4 - Jequitibá" ~ "P4",
                                  "DC5 - Naque SA" ~ "P5",
                                  "DC6 - Periquito" ~ "P6",
                                  "DC7 - Baguari" ~ "P7",
                                  "DC8 - Naque DC" ~ "P8",
                                  "DC9 - Pedra Corrida" ~ "P9")) %>% 
  ggplot(aes(y=Values,
             x=Metadata8 ,
             fill = Metadata8,alpha = 0.1)) + 
  geom_col(col = "#282828",
           # alpha = 0.4,
           linewidth = 0.05) +
  # geom_jitter(aes(col = Metadata8),
  #             alpha = 0.5,
  #             height = 0,
  #             width = 0.25) +
  facet_wrap(nrow = 1,
             ~Indice,
             scales = "free")+
  # facet_grid(vars(`Metadata9`),
  #            vars(Indice),
  #            scales = "free")+
  scale_colour_manual(values = cores_igor, breaks = names(cores_igor), name = "Área de estudo")+
  scale_fill_manual(values =  cores_igor, breaks = names(cores_igor), name = "Área de estudo")+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1, size = 25),
        axis.title = element_text(size = 35),
        strip.text = element_text(size = 35),
        legend.position = "none",
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
        )+
  # labs(subtitle = "Observed species richness by sampling sites")+
  # labs(subtitle = "Diversity index by sampling sites")+
  xlab(label = "Sample sites") +
  ylab(label = "Richness") +
  guides(shape="none")

div_est_phylo_id_igor

## plot indices ----

div_est_phylo_plot_id_igor <-
  div_est_phylo_id_igor %>% 
  mutate("Metadata8" = case_match(Metadata8,
                                  "DC1 - Ponte Nova" ~ "P1",
                                  "DC2 - Coronel Fabriciano" ~ "P2",
                                  "DC3 - Matipó" ~ "P3",
                                  "DC4 - Jequitibá" ~ "P4",
                                  "DC5 - Naque SA" ~ "P5",
                                  "DC6 - Periquito" ~ "P6",
                                  "DC7 - Baguari" ~ "P7",
                                  "DC8 - Naque DC" ~ "P8",
                                  "DC9 - Pedra Corrida" ~ "P9")) %>% 
  filter(Indice %in% c("Observed",
                       "Chao1",
                       "Shannon")) %>% 
  ggplot(aes(y=Values,
             x=Metadata8 ,
             fill = Metadata8,alpha = 0.1)) +
  geom_col(col = "#282828",
           # alpha = 0.4,
           linewidth = 0.05) +
  # geom_jitter(aes(col = Metadata8),
  #             alpha = 0.5,
  #             height = 0,
  #             width = 0.25)+
  facet_wrap(nrow = 3,
             ~Indice,
             scales = "free")+
  # facet_grid(vars(`Metadata9`),
  #            vars(Indice),
  #            scales = "free")+
  scale_colour_manual(values = cores_igor, breaks = names(cores_igor), name = "Área de estudo")+
  scale_fill_manual(values =  cores_igor, breaks = names(cores_igor), name = "Área de estudo")+
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1, size = 20),
        axis.title = element_text(size = 35),
        strip.text = element_text(size = 35),
        legend.position = "none",
        panel.grid.major = element_blank())+
  # labs(subtitle = "Diversity index by sample sites")+
  xlab(label = "Sample sites") +
  ylab(label = "Value") +
  guides(shape="none")

div_est_phylo_plot_id_igor

# Tabela valores indices de diversidade ----

div_est_phylo_id_igor %>% colnames()

tbl_sup3_igor <- FINAL_TBL_igor %>% 
  select(`Curated ID`,`Final ID (BLASTn)`,
         `BLASTn pseudo-score`,
         `1_subject`, `1_indentity`, `1_qcovhsp`,
         `2_subject`, `2_indentity`, `2_qcovhsp`,
         `3_subject`, `3_indentity`, `3_qcovhsp`) %>% 
  unique() %>% 
  mutate("Final DB ID" = as.character(""))


# anosim test

anosim(x = all_IDs_NMDS_df_igor[, -c(1:7)], grouping = all_IDs_NMDS_df_igor$`Metadata 10`, permutations = 9999, distance = "bray")
anosim(x = all_IDs_NMDS_df_igor[, -c(1:7)], grouping = all_IDs_NMDS_df_igor$`Metadata 9`, permutations = 9999, distance = "bray")


# ACCUMULATION CURVE ----

######## Function to transform values into 1 and zeros
jaccarize <- function(x) {
  # Using ifelse to check each element of the vector x
  ifelse(x == 0, 0, 1)
}
###################


# transform proportions to presence/abcense 
all_IDs_NMDS_df_igor_jc <- all_IDs_NMDS_df_igor %>%
  # as_tibble()
  mutate(across(starts_with('ID_'), jaccarize))


#tirado daqui
# https://vegandevs.github.io/vegan/reference/specaccum.html

spec1 <- vegan::specaccum(comm = all_IDs_NMDS_df_igor_jc[,8:ncol(all_IDs_NMDS_df_igor_jc)]
                          # , gamma = "jack1"
                          )
spec2 <- vegan::specaccum(comm = all_IDs_NMDS_df_igor_jc[,8:ncol(all_IDs_NMDS_df_igor_jc)],
                          method = "random"
                          # , gamma = "jack1"
                          )


plot(spec1)


plot(spec1, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="lightgrey")
boxplot(spec2, col="#00CC32", add=TRUE, pch="*")


spec2 %>% str()
spec2 %>% summary()



library(summarytools)


ctable(spec2) 



tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd)
}

spec2$sites
spec2$richness
spec2$method
spec2$sd



spec2_tidy <- spec2$perm %>% 
  reshape2::melt() %>% 
  rename("Sample" = "Var1",
         "Permutation" = "Var2",
         "Value" = "value") %>% 
  as_tibble()

spec2_tidy_rch <- tibble("Richness" = spec1$richness,
                         "Sample" = spec1$sites) 


collector_plot <-
  spec2_tidy %>% 
  ggplot(aes(x = `Sample`,y = Value,group=`Sample`)) +
  # geom_point() +
  geom_boxplot(col ="#30123BFF",
               fill ="#30123BFF",notch = TRUE,
               alpha = 0.33,width = 0.5,outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.25,linetype=2) +
  geom_jitter(size=0.01,
              width = 0.25,
              # height = 0,
              col="#30123BFF") +
  geom_line(data = spec2_tidy_rch,
            inherit.aes = F,
            linewidth = 2.5,
            col ="#30123BFF",
            alpha = 0.50,
            aes(x = `Sample`, y = Richness)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 20)) +
  ylab(label = "Richness") +
  xlab(label = "Sample Sites") +
  scale_x_continuous(breaks = seq(1,12,1))
# + 
#   scale_y_continuous(limits = c(0,100), breaks = seq(20, 100, by = 20))

collector_plot

ggsave(file = paste0("/home/igorhan/projetos/Doce12Sdb/12Sdb/eDNA_Doce/Resultados/accumulation_curve","-",
                     Sys.Date(),
                     ".pdf",
                     collapse = ""),
       plot = collector_plot,
       device = "pdf",
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300)




# native tbl

native_tbl_igor <- tibble("site" = c("DC1", "DC2", "DC3", "DC4", "DC5", "DC6", "DC7", "DC8", "DC9"),
                     "native" = c(23, 30, 25, 17, 27, 22, 22, 22, 20),
                     "non-native" = c(7, 6, 5, 6, 12, 11, 12, 11, 9)
)

cores_igor <- c("native" = "#a0db8e",
           "non-native" = "#ff265b")


native_tbl_igor %>% 
  tidyr::pivot_longer(cols = c( "native", "non-native")) %>% 
  ggplot(aes(x = site,
             y = value,
             fill = name,
             col = name,
             alpha = 0.1)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_colour_manual(values = cores_igor, breaks = names(cores_igor), name = "Área de estudo")+
  scale_fill_manual(values =  cores_igor, breaks = names(cores_igor), name = "Área de estudo")+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none",
         col = "none") +
  xlab(label = "Site") +
  ylab(label = "Richness")

###############

# SIMPER test ----
tbl_natives_igor <- tibble("species" = c("ID_Aequidens_sp",
                                    "ID_Astyanax_sp",
                                    "ID_Astyanax_sp_B",
                                    "ID_Brycon_dulcis",
                                    "ID_Brycon_sp_B",
                                    "ID_Characidium_sp",
                                    "ID_Cichla_kelberi",
                                    "ID_Clarias_gariepinus",
                                    "ID_Corydoras_sp",
                                    "ID_Crenicichla_Lacustria_lacustris",
                                    "ID_Delturus_carinotus",
                                    "ID_Euryochus_thysanos",
                                    "ID_Geophagus_brasiliensis",
                                    "ID_Gymnotus_sp_A",
                                    "ID_Gymnotus_sp_B",
                                    "ID_Harttia_sp",
                                    "ID_Hemigrammus_sp",
                                    "ID_Hoplias_gr_malabaricus",
                                    "ID_Hoplias_intermedius",
                                    "ID_Hoplosternum_littorale",
                                    "ID_Hyphessobrycon_eques",
                                    "ID_Hyphessobrycon_sp_A",
                                    "ID_Hyphessobrycon_sp_B",
                                    "ID_Hypomasticus_sp_A",
                                    "ID_Hypomasticus_sp_B",
                                    "ID_Hypostomus_affinis",
                                    "ID_Hypostomus_sp",
                                    "ID_Imparfinis_sp",
                                    "ID_Knodus_moenkhausi",
                                    "ID_Knodus_sp",
                                    "ID_Lophiosilurus_alexandri",
                                    "ID_Loricariichthys_castaneus",
                                    "ID_Megaleporinus_sp",
                                    "ID_Moenkhausia_costae",
                                    "ID_Neoplecostomus_sp",
                                    "ID_Oreochromis_niloticus",
                                    "ID_Pachyurus_adspersus",
                                    "ID_Paratocinclus_sp",
                                    "ID_Pareiorhaphis_sp_A",
                                    "ID_Pareiorhaphis_sp_B",
                                    "ID_Piimelodella_sp",
                                    "ID_Pimelodus_maculatus",
                                    "ID_Poecilia_sp_B",
                                    "ID_Poeciliidae_sp_A",
                                    "ID_Prochilodus_sp",
                                    "ID_Prochilodus_vimboides",
                                    "ID_Psalidodon_cf_rivularis",
                                    "ID_Psalidodon_sp_A",
                                    "ID_Psalidodon_sp_B",
                                    "ID_Pseudauchenipterus_affinis",
                                    "ID_Pseudoplatystoma_hybrid",
                                    "ID_Pygocentrus_nattereri",
                                    "ID_Pygocentrus_sp_A",
                                    "ID_Rhamdia_quelen",
                                    "ID_Rineloricaria_sp",
                                    "ID_Salminus_brasiliensis",
                                    "ID_Serrapinnus_heterodon",
                                    "ID_Stethaprioninae",
                                    "ID_Synbranchus_sp_A",
                                    "ID_Synbranchus_sp_B",
                                    "ID_Trachelyopterus_striatulus",
                                    "ID_Trichomycterus_sp"),
                      "status" = c("non-native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "non-native",
                                   "native",
                                   "native" ,
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "non-native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "native",
                                   "native",
                                   "non-native",
                                   "non-native",
                                   "native",
                                   "native"))


COMM <- all_IDs_NMDS_df_igor


simper_tbl_igor_habitat_igor <- vegan::simper(comm = COMM[,str_detect(colnames(COMM),pattern = "^ID_")],
              group = COMM[,"Metadata 9"],
              permutations = 999)



simper_tbl_igor_river_igor <- vegan::simper(comm = COMM[,str_detect(colnames(COMM),pattern = "^ID_")],
              group = COMM[,"agrupador"],
              permutations = 999)


simper_tbl_igor_join_igor <- list(c(simper_tbl_igor_habitat_igor,
       simper_tbl_igor_river_igor))



simper_tbl_igor_join_igor[[1]] %>% names()

simper_tbl_igor <- simper_tbl_igor_join_igor

# Initialize an empty data frame to store results
simper_df_igor <- data.frame()

# Loop through each pairwise comparison and extract results
for (comparison in 1:2) {
 
  temp_df_igor <- as.data.frame(simper_tbl_igor[[1]][comparison])
  
  
  temp_df_igor <- temp_df_igor %>%  dplyr::rename_with(~ stringr::str_remove(., pattern = "^.*\\."))
  
  temp_df_igor <- temp_df_igor %>% mutate("Comparison" = simper_tbl_igor[[1]][comparison] %>% names() )
  
  # temp_df_igor$Comparison <- comparison  # Add a column for the comparison
  simper_df_igor <- dplyr::bind_rows(simper_df_igor, temp_df_igor)  # Append to the main data frame
}



simper_df_igor <- simper_df_igor %>% 
  dplyr::mutate("species" = str_replace_all(species, "[\n\r]", ""))

# Save as CSV file

# join simper_native

# 
# simp_full_igor$species %in% c("ID_Knodus_moenkhausi")
# simper_df_igor$species %in% c("ID_Knodus_moenkhausi")

simp_full_igor <- left_join(x = simper_df_igor,
                       y = tbl_natives_igor,
                       by = "species",
                       multiple = "all") 


simp_full_igor <- simp_full_igor %>% 
  rename("var1" = "ava",
         "var2" = "avb") %>%
  dplyr::select("species", "status", "var1", "var2", "Comparison") %>% 
  tidyr::pivot_longer(cols = c("var1", "var2"),
                      names_to = "Habitat",
                      values_to = "SIMPER_value") %>% 
  dplyr::mutate("Habitat" = stringr::str_remove(Habitat, pattern = "_average")) %>% 
  dplyr::mutate("Habitat" = stringr::str_replace(Habitat, pattern = "^l",replacement = "L"))

# HEATMAP ----

# create tbls ----

var_simper_igor <-
simp_full_igor %>% 
  dplyr::filter(!SIMPER_value == 0) %>%
  dplyr::filter(Comparison %in% "Tributaries_Mainstream") %>%
  dplyr::group_by(species) %>%
  slice_max(SIMPER_value, n = 1) %>%
  ungroup() %>%
  mutate("Habitat" = case_match(`Habitat`,
                                "var1" ~ "Tributaries",
                                "var2" ~ "Mainstream")) %>% 
  dplyr::select(-c(Comparison, SIMPER_value)) %>%
  
  ############## arrumar problemas com nomes  ############## 

  mutate("species" = str_remove_all(`species`,
                                    pattern = "ID_") %>% 
                     str_replace_all(pattern = "_", 
                                     replacement = " ") %>% 
                     str_replace_all(pattern = "sp", 
                                     replacement = "sp.") %>% 
                     str_replace_all(pattern = "cf",
                                     replacement = "cf.") %>% 
                     str_replace_all(pattern = "gr", 
                                     replacement = "gr.") %>% 
                     str_replace_all(pattern = "Pachyurus adsp.ersus",
                                     replacement = "Pachyurus adspersus") %>% 
                     str_replace_all(pattern = "Hemigr.ammus sp.",
                                     replacement = "Hemigrammus sp.")) %>% 
  mutate("species" = case_when(`species` %in% c("Crenicichla Lacustria lacustris") ~ "Crenicichla (Lacustria) lacustris",
                               TRUE ~ `species`)) %>% 
  BiocGenerics::unique()

species_tile_igor <- 
FINAL_ID_table_id_igor %>% 
  tibble::rownames_to_column(var = "species") %>% 
  mutate("species" = str_replace_all(`species`,
                                     pattern = "Knodus moenkhausi\r\n",
                                     replacement = "Knodus moenkhausi")) %>% 
  dplyr::rename_with(~ str_replace(., "DC0", "P"))

# plot ----

# agrupamento das especies (linhas)
simper_sps_df_igor <- var_simper_igor %>%  as.data.frame() %>% 
  select("Habitat") %>% 
  dplyr::rename("SIMPER" = "Habitat")

rownames(simper_sps_df_igor) <- var_simper_igor$species

# agrupamento dos pontos(colunas)
habitat_pontos_igor <- metadata_tbl_id_igor[2] %>% 
  dplyr::rename("Sites" = "Metadata8") %>% 
  mutate("row" = c("P4", "P3", "P2", "P5", "P8", "P7", "P6", "P9", "P1")) %>% 
  `rownames<-`(.$row) %>% 
  select(-c("row"))

# matriz de especies
species_mat_igor <- as.matrix(species_tile_igor[,2:10])

rownames(species_mat_igor) <- species_tile_igor$species

# plot

library(cowplot)

ann_colors_igor = list(
  "Sites" = c(Mainstream ="#ff950e", Tributaries =  "#004586"),
  "SIMPER" = c(Mainstream ="#f7b259", Tributaries =  "#395d80")
)

# função para deixar nomes em italico

newnames_igor <- lapply(
  rownames(species_mat_igor),
  function(x) bquote(italic(.(x))))

# Gere o heatmap

p_igor <- pheatmap(mat = species_mat_igor,
              cluster_rows = T,
              annotation_col = habitat_pontos_igor,
              labels_row = as.expression(newnames_igor),
              annotation_row = simper_sps_df_igor,
              color = colorRampPalette(c("white", "#99D17B"))(50)
              #, annotation_colors = ann_colors_igor
              )



# Salve usando cowplot e ggsave
ggsave(
  filename = "heatmap_ggsave_cowplot.pdf",
  plot = ggdraw() + draw_grob(p$gtable),
  width = 9,
  height = 10,
  dpi = 300
)
