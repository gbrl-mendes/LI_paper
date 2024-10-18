---
title: "Lagoa dos Ingleses Paper Script"
author: "Mendes, GA"
date: "05/2024"
---
  
# Loading libraries ----
{ 
  library("Biostrings")
  library("DECIPHER")
  library("phyloseq")
  library("vegan")
  library("tidyverse")
  library("dplyr")
  library("readxl")
  library("ggpubr")
  library("ggplot2")
  }

# Paths ----

{
  prjct_path <- "/home/gabriel/projetos/LI_paper"
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
  scripts_path <- paste0(prjct_path,"/scripts")
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figures")
  tbl_path <- paste0(results_path,"/tables")
  
  paths <- c(scripts_path, results_path, figs_path, tbl_path)
}

# Criando caminhos
for (dir in paths) {
  if (dir.exists(dir)) {
    print(paste("The directory", dir, "already exists!"))
  } else {
    print(paste("Making", dir, "!" ))
    dir.create(dir) 
  }
}

# ---- Post-curated IDs changes ----
## Data acquiring ---- 

t_pre_raw_results_tbl <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/curated/curated-Complete_analysis_results-2024-01-10.xlsx") %>% tibble()
t_curated_ids_tbl <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/curated/curated_lagoa_ingleses-ASVs_x_amostras-2024-01-09.xlsx") %>% tibble()
t_blast_tax <- read.csv("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/tax_blast.csv", sep = ",", check.names = FALSE)

## Adicionando a t_pre_raw_results_tbl os nomes que foram curados manualmente ----

## 1o selecionando apenas as colunas que importam

# Curated IDs
t_curated_ids_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

t_curated_ids_tbl <-
  t_curated_ids_tbl %>%
  select(c(
    "ASV header",
    "ASV (Sequence)",
    "Curated ID",
    "Final ID (BLASTn)",
    "BLASTn pseudo-score",
    "Class (BLASTn)",
    "Curated Order (BLASTn)"
    # ,
    # "Obs. Curadoria",
    # "Possible contamination",
    # "Curated Order"
  )) %>%
  unique()

View(t_curated_ids_tbl)

# Verificar se existem IDs diferentes para a mesma ASV
t_curated_ids_tbl$`ASV header`[which(t_curated_ids_tbl$`ASV header` %>% duplicated())]

# Complete analysis
t_pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

t_pre_raw_results_tbl <-
  t_pre_raw_results_tbl %>%
  select(c("Sample",
           "Unique_File_name",
           "OTU",
           "Abundance",
           "Sample total abundance",
           "Relative abundance on sample",
           "Relative abundance to all samples",
           "Obs. Curadoria",
           "Possible contamination",
           "Read origin",
           "Primer expected length",
           "ASV Size (pb)",
           "ASV header")) %>%
  mutate(Sample = str_replace(Sample, "__", "_"))

View(t_pre_raw_results_tbl)

## 2o renomeando segundo os nomes curados e reordenando

t_pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

t_raw_results_tbl <- t_pre_raw_results_tbl %>%
  left_join(t_curated_ids_tbl,
            by = c("ASV header")) %>%
  select(c("ASV header", #oriundas da t_curated_ids_tbl
           "ASV (Sequence)",
           "Curated ID",
           "Final ID (BLASTn)",
           "BLASTn pseudo-score",
           "Class (BLASTn)",
           "Curated Order (BLASTn)",
           # "Curated Order",
           "Sample", #oriundas da pre_raw_samples_tbl
           "Unique_File_name",
           "OTU",
           "Abundance",
           "Sample total abundance",
           "Relative abundance on sample",
           "Relative abundance to all samples",
           "Obs. Curadoria",
           # "Possible contamination",
           "Read origin",
           "Primer expected length",
           "ASV Size (pb)",
           "ASV header"
  ))

View(t_raw_results_tbl)

# Verificar se existem IDs diferentes para a mesma ASV
t_raw_results_tbl$`ASV header`[which(t_curated_ids_tbl$`ASV header` %>% duplicated())]

# Verificando o total de reads por amostra
raw_reads <- t_raw_results_tbl %>% group_by(Sample) %>%
  summarise(total_abd = sum(Abundance))

View(raw_reads)

# Comparar com o total que saiu do pipeline para ver se nao ha perdas:

# > raw_reads < smp_abd_ID_Final %>% group_by(Sample) %>%
#   +  summarise(total_abd = sum(Abundance))
# > print(raw_reads, n= 21)

# A tibble: 21 × 2
#   Sample             total_abd
#   <chr>                  <dbl>
# 1 L1_jan22              188658
# 2 L1_nov21               66426
# 3 L1_nov_dec_20_mi       41539
# 4 L1_out21              384154
# 5 L2_dez20              357293
# 6 L2_jan22              249320
# 7 L2_nov20              207421
# 8 L2_nov21              184214
# 9 L2_out21              530824
# 10 L3_jan22                 37
# 11 L3_nov21              49612
# 12 L3_out21             370515
# 13 L4_jan22                 32
# 14 L4_nov21             514706
# 15 L4_out21                123
# 16 STX__L1_nov21        303077
# 17 STX__L2_nov21        312437
# 18 STX__L3_nov21        278561
# 19 STX__L4_nov21          1264
# 20 br_jan_22                44
# 21 br_nov21                 37

# Verificando se ha diferencas entre a tabela de ASVs e a tabela de IDs curadas
dif_raw_pre <- t_pre_raw_results_tbl %>%
  anti_join(t_raw_results_tbl)

View(dif_raw_pre) #vazio e' bom!

# Apos ver que tudo que esta na tabela eh o que saiu do pipeline, podemos retirar 
# tudo o que estiver fora do intervalo do amplicon, identificacoes NA, brancos e
# possiveis contaminacoes

t_filt_results_tbl <-
  t_raw_results_tbl %>% 
  filter(!`Primer expected length` %in% "out of range") %>% 
  filter(!`Curated ID` == "NA") %>% 
  filter(!(`Class (BLASTn)` %in% c("Mammalia", "Aves") & `BLASTn pseudo-score` < 98)) %>% 
  filter(!(Sample %in% c("EM113_NEGPCR1",
                         "EM135c4c5_NEGPCR1",
                         "EM135c4c5_NEGPCR2",
                         "EM149_NEGPCR1",
                         "EM149_NEGPCR2",
                         "EM156_NegPCR1",
                         "S724_NEGPCR1",
                         "br_jan_22",
                         "br_nov21"))) %>% 
  # filter(!(`Possible contamination` %in% c("Possible contamination"))) %>% 
  filter(!(`Obs. Curadoria` %in% c("Possible contamination")))

View(t_filt_results_tbl)

# corrigindo o RRA

# Verificando o RRA
t_filt_results_tbl %>% 
  # t_grouped_by_ID_BLASTid %>%
  group_by(Sample) %>% 
  summarize(total_RRA = sum(`Relative abundance on sample`)) # Veja como que algumas amostras perderam muito!

# >     t_filt_results_tbl %>% 
#   +     # t_grouped_by_ID_BLASTid %>%
#   +     group_by(Sample) %>% 
#   +       summarize(total_RRA = sum(`Relative abundance on sample`))
# # A tibble: 19 × 2
# Sample           total_RRA
# <chr>                <dbl>
#   1 L1_jan22             0.307
# 2 L1_nov21             1.00 
# 3 L1_nov_dec_20_mi     0.997
# 4 L1_out21             0.999
# 5 L2_dez20             1.00 
# 6 L2_jan22             0.489
# 7 L2_nov20             1.00 
# 8 L2_nov21             1.00 
# 9 L2_out21             1.00 
# 10 L3_jan22             0.622
# 11 L3_nov21             1.00 
# 12 L3_out21             1.00 
# 13 L4_jan22             0.656
# 14 L4_nov21             1.00 
# 15 L4_out21             1.00 
# 16 STX_L1_nov21         0.949
# 17 STX_L2_nov21         0.994
# 18 STX_L3_nov21         0.910
# 19 STX_L4_nov21         0.739

t_filt_results_tbl <- t_filt_results_tbl %>%
  mutate("Curated Relative abundance to all samples" = 0,
         "Curated Relative abundance on sample" = 0,
         "Curated Sample total abundance" = 0)

abd_total <- sum(t_filt_results_tbl$Abundance)

t_filt_results_tbl <- t_filt_results_tbl %>%
  group_by(Unique_File_name,`Read origin`) %>%   
  mutate("Curated Sample total abundance" = sum(Abundance),
         "Curated Relative abundance to all samples" = round((Abundance/abd_total),digits = 10),
         "Curated Relative abundance on sample" =  round((Abundance/`Curated Sample total abundance`),digits = 10)) %>%
  relocate(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`,
           `Curated Sample total abundance`,`Curated Relative abundance to all samples`,`Curated Relative abundance on sample`) %>%
  ungroup() %>% 
  select(-c(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`)) %>% 
  rename("Sample total abundance" = "Curated Sample total abundance",
         "Relative abundance to all samples" = "Curated Relative abundance to all samples",
         "Relative abundance on sample"= "Curated Relative abundance on sample")

# Verificando se funcionou
t_filt_results_tbl %>% 
  # t_grouped_by_ID_BLASTid %>%
  group_by(Sample) %>% 
  summarize(total_RRA = sum(`Relative abundance on sample`)) # Compare com o resultado anterior!

## 3o adicionando metadados que faltaram 

t_filt_results_tbl$Sample %>% unique() %>% sort() %>% paste0(collapse = '",\n"') %>% cat()
{
  # sample
  Sample <- c("L1_nov_dec_20_mi", 
              "L2_nov20",
              "L2_dez20",
              "L1_out21",
              "L2_out21",
              "L3_out21",
              "L4_out21",
              "L1_nov21",
              "L2_nov21",
              "L3_nov21",
              "L4_nov21",
              "STX_L1_nov21",
              "STX_L2_nov21",
              "STX_L3_nov21",
              "STX_L4_nov21",
              "L1_jan22",
              "L2_jan22",
              "L3_jan22",
              "L4_jan22"
  )
  
  # expedition
  expedition <- c("Novembro e Dezembro 2020",
                  "Novembro 2020",
                  "Dezembro 2020",
                  "Outubro 2021",
                  "Outubro 2021",
                  "Outubro 2021",
                  "Outubro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Janeiro 2022",
                  "Janeiro 2022",
                  "Janeiro 2022",
                  "Janeiro 2022"
  )
  
  # year
  year <- c("2020",
            "2020",
            "2020",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2022",
            "2022",
            "2022",
            "2022"
  )
  
  # new_name
  new_name <- c("Fundação",
                "Ponte",
                "Ponte",
                "Prainha",
                "Barragem",
                "Ponte",
                "Fundação",
                "Prainha",
                "Barragem",
                "Ponte",
                "Fundação",
                "Prainha",
                "Barragem",
                "Ponte",
                "Fundação",
                "Prainha",
                "Barragem",
                "Ponte",
                "Fundação"
                )
  
  # filter
  filter <- c("MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "Sterivex",
              "Sterivex",
              "Sterivex",
              "Sterivex",
              "MCE",
              "MCE",
              "MCE",
              "MCE"
            )
  
  # run
  run <- c("run2_ago21",
           "run4_out21",
           "run4_out21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "EM156",
           "EM156",
           "EM156",
           "EM156",
           "EM156",
           "EM156",
           "EM156",
           "EM156"
           )
  
  # level
  level <- c("Full",
             "Full",
             "Full",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Low",
             "Full",
             "Full",
             "Full",
             "Full"
             )
}

# As tibble
t_metadata_tbl <- tibble(run, 
                       Sample,
                       filter,
                       new_name, 
                       expedition, 
                       year)

# Mergindo os metadados com t_filt_results_tbl

t_filt_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

t_filt_results_tbl$Sample %>% unique()

t_curated_full_tbl <- left_join(t_filt_results_tbl, t_metadata_tbl,
                              by = "Sample") %>%
  select(c("Sample",
           "run",
           "new_name",
           "expedition",
           "year",
           "filter",
           "Curated ID",
           # "Curated Order",
           # "Obs. Curadoria",
           "Final ID (BLASTn)",
           "BLASTn pseudo-score",
           "Class (BLASTn)",
           "Curated Order (BLASTn)",
           # "Order (DADA2)",
           # "Curated Max. taxonomy",
           # "Primer",
           "OTU",
           "Abundance",
           "Sample total abundance",
           "Relative abundance on sample",
           "Relative abundance to all samples",
           "Read origin",
           "Primer expected length",
           "ASV Size (pb)",
           "ASV header",
           "ASV (Sequence)"
  )) 

View(t_curated_full_tbl)

diff_cur_raw <- dplyr::anti_join(#t_curated_ids_tbl,
  t_filt_results_tbl,
  t_curated_full_tbl,
  # ver se ocorreu sem problemas o left_join
  by = "ASV header")

View(diff_cur_raw) #vazio e' bom!

# Mergindo os metadados com out_results_tbl

out_results_full <- left_join(out_results_tbl, t_metadata_tbl,
                              by = "Sample") %>%
  select(c("Sample",
           "run",
           "new_name",
           "expedition",
           "year",
           "filter",
           "Curated ID",
           "Curated Order (BLASTn)",
           # "Obs. Curadoria",
           "Final ID (BLASTn)",
           "BLASTn pseudo-score",
           "Class (BLASTn)",
           # "Order (DADA2)",
           # "Curated Max. taxonomy",
           # "Primer",
           "OTU",
           "Abundance",
           # "Sample total abundance",
           # "Relative abundance on sample",
           # "Relative abundance to all samples",
           "Read origin",
           "Primer expected length",
           "ASV Size (pb)",
           "ASV header",
           "ASV (Sequence)"
  ))

View(out_results_full)

## Dados de hits por DB ----
# Tive que rodar novamente o BLASTr por que as infos sobre a origem dos hits nao foram salvos adequadamente.
# As tabelas do pipeline geradas em 23-02 possuem esses dados completos. 

## Obter a nova tabela raw_results
{
  NW_t_pre_raw_results_tbl <- read_excel("/home/gabriel/projetos/peixes-eDNA/analises/dez_23/runs_2_4_5_EM156/results/lagoa_ingleses-Complete_analysis_results-2024-02-23.xlsx") %>% tibble()
  
  ## Adicionando o RRA correto a NW_pre_raw_results 
  {
    NW_t_pre_raw_results_tbl <- NW_t_pre_raw_results_tbl %>%
      mutate("Curated Relative abundance to all samples" = 0,
             "Curated Relative abundance on sample" = 0,
             "Curated Sample total abundance" = 0)
    
    abd_total <- sum(NW_t_pre_raw_results_tbl$Abundance)
    
    NW_t_pre_raw_results_tbl <- NW_t_pre_raw_results_tbl %>%
      group_by(Unique_File_name,`Read origin`) %>%   
      mutate("Curated Sample total abundance" = sum(Abundance),
             "Curated Relative abundance to all samples" = round((Abundance/abd_total),digits = 10),
             "Curated Relative abundance on sample" =  round((Abundance/`Curated Sample total abundance`),digits = 10)) %>%
      relocate(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`,
               `Curated Sample total abundance`,`Curated Relative abundance to all samples`,`Curated Relative abundance on sample`) %>%
      ungroup() %>% 
      select(-c(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`)) %>% 
      rename("Sample total abundance" = "Curated Sample total abundance",
             "Relative abundance to all samples" = "Curated Relative abundance to all samples",
             "Relative abundance on sample"= "Curated Relative abundance on sample")
  }
  
  ## Selecionando colunas
  NW_t_pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  
  NW_t_pre_raw_results_tbl <-
    NW_t_pre_raw_results_tbl %>%
    select(c("Sample",
             "Unique_File_name",
             "OTU",
             "Abundance",
             "Sample total abundance",
             "Relative abundance on sample",
             "Relative abundance to all samples",
             "Obs. Curadoria",
             "Possible contamination",
             "Read origin",
             "Primer expected length",
             "ASV Size (pb)",
             "ASV header")) %>%
    mutate(Sample = str_replace(Sample, "__", "_"))
  
  View(NW_t_pre_raw_results_tbl)
  
  ## Adicionando a info de contaminacao
  
  NW_contam <- t_pre_raw_results_tbl %>% 
    select(`Obs. Curadoria`,
           `Possible contamination`,
           `ASV header`)
  
  View(NW_contam)
  
  NW_t_pre_raw_results_tbl <- NW_t_pre_raw_results_tbl %>% 
    left_join(NW_contam,
              by = "ASV header")
  
  View(NW_t_pre_raw_results_tbl)
  
  # Verificar se os resultados sao os mesmos com t_pre_raw_results_tbl
  
  dif_NW_pre <- NW_t_pre_raw_results_tbl %>%
    anti_join(t_pre_raw_results_tbl)
  
  View(dif_NW_pre)
}

## Obter a nova tabela curated_IDs
NW_curated_ids_tbl <- read_excel("/home/gabriel/projetos/peixes-eDNA/analises/dez_23/runs_2_4_5_EM156/results/lagoa_ingleses-ASVs_x_amostras-2024-02-23.xlsx") %>% tibble()

## Selecionando apenas as colunas que importam para comparar

# Curated IDs

TT_curated_ids_tbl <-
  t_curated_ids_tbl %>%
  select(c(
    "ASV header",
    "ASV (Sequence)",
    # "Curated ID",
    "Final ID (BLASTn)",
    "BLASTn pseudo-score",
    "Class (BLASTn)",
    # "Curated Order (BLASTn)"
    # ,
    # "Obs. Curadoria",
    # "Possible contamination",
    # "Curated Order"
  )) %>%
  unique()

View(TT_curated_ids_tbl)

NW_curated_ids_tbl <-
  NW_curated_ids_tbl %>%
  select(c(
    "ASV header",
    "ASV (Sequence)",
    # "Curated ID",
    "Final ID (BLASTn)",
    "BLASTn pseudo-score",
    "Class (BLASTn)",
    # "Curated Order (BLASTn)"
    # ,
    # "Obs. Curadoria",
    # "Possible contamination",
    # "Curated Order"
  )) %>%
  unique()

View(NW_curated_ids_tbl)

# Verificar se os resultados sao os mesmos com t_pre_raw_results_tbl

dif_NW_pre <- NW_curated_ids_tbl %>%
  anti_join(t_curated_ids_tbl)    

View(dif_NW_pre)

# Podemos ver como a diferenca entre as tabelas e` unicamente a correcao do Final ID, mas aparentemente
# os resultados do BLASTn foram identicos. Posso entao usar a ASV_header para pegar os resultados do DB

hits_DB <- NW_curated_ids_tbl %>% 
  select(`ASV header`,
         `Final ID (BLASTn)`,
         `1_DB`,
         `2_DB`,
         `3_DB`)

View(hits_DB)

## Tabelas ----

# Criacao da lista com os possiveis nomes atribuidos as ASVs

t_curated_full_tbl %>% colnames()
t_curated_full_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

# Agrupamento da ASVs que possuem os mesmos atributos abaixo

t_grouped_by_ID_BLASTid <- t_curated_full_tbl %>%
  relocate(c("Sample",
             "expedition",
             "year",
             "Final ID (BLASTn)",
             "BLASTn pseudo-score",
             "Curated ID",
             "OTU",
             "filter",
             "Curated Order (BLASTn)",
             "Relative abundance on sample",
             "new_name",
             "Class (BLASTn)",
             # "Order (DADA2)",
             # "Curated Order",
             # "Curated Max. taxonomy",
             "Abundance",
             "ASV header",
             "ASV (Sequence)",
             "Primer expected length"
  )) %>%
  group_by(Sample, `Curated ID`,`Final ID (BLASTn)`,`BLASTn pseudo-score`, expedition, new_name, year, filter, `ASV header`, `Curated Order (BLASTn)`) %>%
  summarize("Class (BLASTn)" = unique(`Class (BLASTn)`),
            "OTU" = unique(OTU),
            # "Curated Order" = unique(`Curated Order`),
            # "Order (DADA2)" = unique(`Order (DADA2)`),
            # "Curated Max. taxonomy" = unique(`Curated Max. taxonomy`),
            "Abundance" = sum(Abundance),  # Agregando a coluna Abundance
            "ASV (Sequence)" = unique(`ASV (Sequence)`),
            "Primer expected length" = unique(`Primer expected length`)
  ) %>%
  mutate(level = case_when(`year` == 2020 ~ "Cheio", 
                           `year` == 2022 ~ "Cheio", 
                           TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o level da lagoa
  ungroup()

View(t_grouped_by_ID_BLASTid)

t_grouped_by_ID_BLASTid[which(t_grouped_by_ID_BLASTid %>% duplicated())]

# verificando o total de reads por amostra

grouped_reads <- t_grouped_by_ID_BLASTid %>% group_by(Sample) %>% 
  summarise(total_abd = sum(Abundance))

View(grouped_reads)

# comparar com o total que saiu do pipeline para ver se nao ha perdas:

# > raw_reads < smp_abd_ID_Final %>% group_by(Sample) %>% 
#   +  summarise(total_abd = sum(Abundance))
# > print(raw_reads, n= 21)

# A tibble: 21 × 2
#   Sample             total_abd
#   <chr>                  <dbl>
# 1 L1_jan22              188658
# 2 L1_nov21               66426
# 3 L1_nov_dec_20_mi       41539
# 4 L1_out21              384154
# 5 L2_dez20              357293
# 6 L2_jan22              249320
# 7 L2_nov20              207421
# 8 L2_nov21              184214
# 9 L2_out21              530824
# 10 L3_jan22                 37
# 11 L3_nov21              49612
# 12 L3_out21             370515
# 13 L4_jan22                 32
# 14 L4_nov21             514706
# 15 L4_out21                123
# 16 STX__L1_nov21        303077
# 17 STX__L2_nov21        312437
# 18 STX__L3_nov21        278561
# 19 STX__L4_nov21          1264
# 20 br_jan_22                44
# 21 br_nov21                 37

# adicionando a taxonomia do BLAST

t_blast_tax <- t_blast_tax[2:86] %>% 
  rename(`ASV (Sequence)` = "ASV")

blast_tax_less <- t_blast_tax %>%
  select(c("ASV (Sequence)",
           "Order (BLASTn)",
           "Family (BLASTn)",
           "Genus (BLASTn)"
  )) 
t_grouped_by_ID_BLASTid %>% colnames()

t_grouped_by_ID_BLASTid <- t_grouped_by_ID_BLASTid %>% 
  left_join(blast_tax_less,
            by = "ASV (Sequence)") %>% 
  mutate("Curated genus" = str_split_fixed(string = .$`Curated ID`, 
                                           pattern = " ",
                                           n = 2)[,1]) %>% 
  select(c("Sample",
           "Curated ID", 
           "Final ID (BLASTn)",
           "BLASTn pseudo-score",
           "ASV (Sequence)",
           "Class (BLASTn)",
           "Curated Order (BLASTn)",
           # "Order (BLASTn)",
           # "Curated Order",
           "Family (BLASTn)",
           "Genus (BLASTn)",
           "Curated genus",
           "expedition",
           "new_name",
           "year", 
           "filter",
           "ASV header",
           "RRA",
           "OTU", 
           "Abundance",
           "ASV (Sequence)",
           "Primer expected length",
           "level" 
  )) %>% 
  unique()

View(t_grouped_by_ID_BLASTid)

## Agrupamento da ASVs que possuem os mesmos atributos abaixo VERSAO OUT ----
{
  out_grouped_by_ID_BLASTid <- out_results_full %>%
    relocate(c("Sample",
               "expedition",
               "year",
               "Final ID (BLASTn)",
               "BLASTn pseudo-score",
               "Curated ID",
               "OTU",
               "filter",
               # "Curated Order",
               "Relative abundance on sample",
               "new_name",
               "Class (BLASTn)",
               # "Order (DADA2)",
               # "Curated Order",
               # "Curated Max. taxonomy",
               "Abundance",
               "ASV header",
               "ASV (Sequence)",
               "Primer expected length"
    )) %>%
    group_by(Sample, `Curated ID`,`Final ID (BLASTn)`,`BLASTn pseudo-score`, expedition, new_name, year, filter, `ASV header`) %>%
    summarize("RRA" = sum(`Relative abundance on sample`),
              "Class (BLASTn)" = unique(`Class (BLASTn)`),
              "OTU" = unique(OTU),
              # "Curated Order" = unique(`Curated Order`),
              # "Order (DADA2)" = unique(`Order (DADA2)`),
              # "Curated Max. taxonomy" = unique(`Curated Max. taxonomy`),
              "Abundance" = sum(Abundance),  # Agregando a coluna Abundance
              "ASV (Sequence)" = unique(`ASV (Sequence)`),
              "Primer expected length" = unique(`Primer expected length`)
    ) %>%
    mutate(level = case_when(`year` == 2020 ~ "Cheio", 
                             `year` == 2022 ~ "Cheio", 
                             TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o level da lagoa
    ungroup()
  
  View(out_grouped_by_ID_BLASTid)
  
  # adicionando a taxonomia do BLAST
  
  out_grouped_by_ID_BLASTid <- out_grouped_by_ID_BLASTid %>% 
    left_join(blast_tax_less,
              by = "ASV (Sequence)") %>% 
    mutate("Curated genus" = str_split_fixed(string = .$`Curated ID`, 
                                             pattern = " ",
                                             n = 2)[,1]) %>% 
    select(c("Sample",
             "Curated ID", 
             "Final ID (BLASTn)",
             "BLASTn pseudo-score",
             "ASV (Sequence)",
             "Class (BLASTn)",
             "Order (BLASTn)",
             # "Curated Order",
             "Family (BLASTn)",
             "Genus (BLASTn)",
             "Curated genus",
             "expedition",
             "new_name",
             "year", 
             "filter",
             "ASV header",
             "RRA",
             "OTU", 
             "Abundance",
             "ASV (Sequence)",
             "Primer expected length",
             "level" 
    )) %>% 
    unique()
  
  View(out_grouped_by_ID_BLASTid)
  
  write.csv(out_grouped_by_ID_BLASTid, paste0(tbl_path, "/out_grouped_2024.csv"))
}

## Alteracoes para o artigo 28 de maio ----

# Alteracoes baseadas nas sugestoes do Daniel, apos ver a tabela suplementar 2 da dissertacao.
# Como essa tabela é a versao wider da tabela dt_all_resume, tive que fazer as alteracoes abaixo
# na tabela t_grouped_by_ID_BLASTid para que eu possa aplicar o filtro sugerido pelo Daniel, que 
# envolve os numeros de reads apos serem agrupados como esta abaixo:

# ---- Tables ----
## Tabela suplementar 1 ----

# Tabela usada para avaliar os dados

# Tabela longer
t_dt_all_resume <- t_grouped_by_ID_BLASTid %>% 
  group_by(level, new_name, filter, year) %>%
  mutate("Abd total" = sum(Abundance)) %>%
  ungroup() %>%
  group_by(`Curated ID`) %>% 
  mutate("Abundancia total" = sum(Abundance)) %>% 
  ungroup %>% 
  group_by(`Final ID (BLASTn)`) %>% 
  mutate("Reads totais" = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(`Curated ID`,
           `Final ID (BLASTn)`,
           `new_name`, 
           # `Abundance`,
           `level`,
           `filter`,
           `year`
  ) %>%
  mutate("RRA no periodo" = (Abundance/`Abd total`)*100) %>% 
  # ungroup() %>% 
  reframe(
    "Reads totais" = `Reads totais`,
    "Abundancia total" = `Abundancia total`,
    "RRA" = sum(`RRA no periodo`),
    "ASVs" = length(unique(`ASV header`)),
    "OTUs" = length(unique(OTU)),
    "level" = level,
    "Ano" = year,
    "Filtro" = filter,
    "Ponto" = new_name,
    "Reads" = sum(Abundance),
    "Classe" = `Class (BLASTn)`,
    "Order" = `Curated Order (BLASTn)`,
    "Family" = `Family (BLASTn)`,
    "Genus" = `Curated genus`
    # ,
    # "ASV header" = `ASV header`,
    # "OTU" = OTU
  ) %>% 
  mutate(
    "RRA_formatado" = format(RRA, scientific = TRUE)
  ) %>%
  unique() 

View(dt_all_resume)

# options(scipen = 99,
#         digits = 3)


t_dt_all_resume %>%
  # dplyr::group_by(Sample) %>%    
  dplyr::group_by(Ponto, level, Filtro, Ano) %>%
  # dplyr::group_by(`Curated ID`, new_name, level, filter, year) %>%
  # summarise(total_RRA = sum(`RRA no periodo`, na.rm = TRUE))
  summarise(total_RRA = sum(RRA, na.rm = TRUE)) %>%
  ungroup()


# Teste de RRA para ver se foi calculado corretamente
dt_all_resume %>%
  filter(level == "Vazio") %>%
  filter(new_name == "Fundação") %>%
  filter(Filtro == "MCE") %>%
  filter(Ano == "2021") %>% 
  pull(RRA) %>%
  sum()  

# Tabela wider com raw data
wider_dt_all_resume <- dt_all_resume %>% 
  select(-c("RRA_formatado")) %>%
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  unite(new_name, level, Ano, Filtro, col= "ponto_level_ano_filtro") %>% 
  pivot_wider(id_cols = c("Classe","Curated ID", "Final ID (BLASTn)", "Reads totais", "Abundancia total"),
              names_from = ponto_level_ano_filtro,
              values_from =  c("RRA","ASVs","OTUs","Reads"),
              names_glue = "{ponto_level_ano_filtro}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate("Classe","Curated ID","Final ID (BLASTn)","Reads totais", "Abundancia total",
           #2020 MCE
           "Ponte_Cheio_2020_MCE_Reads",
           "Ponte_Cheio_2020_MCE_ASVs",
           "Ponte_Cheio_2020_MCE_OTUs",
           "Ponte_Cheio_2020_MCE_RRA",
           "Fundação_Cheio_2020_MCE_Reads",
           "Fundação_Cheio_2020_MCE_ASVs",
           "Fundação_Cheio_2020_MCE_OTUs",
           "Fundação_Cheio_2020_MCE_RRA",
           #2021 MCE
           "Prainha_Vazio_2021_MCE_Reads",
           "Prainha_Vazio_2021_MCE_ASVs",
           "Prainha_Vazio_2021_MCE_OTUs",
           "Prainha_Vazio_2021_MCE_RRA",
           "Barragem_Vazio_2021_MCE_Reads",
           "Barragem_Vazio_2021_MCE_ASVs",
           "Barragem_Vazio_2021_MCE_OTUs",
           "Barragem_Vazio_2021_MCE_RRA",
           "Ponte_Vazio_2021_MCE_Reads",
           "Ponte_Vazio_2021_MCE_ASVs",
           "Ponte_Vazio_2021_MCE_OTUs",
           "Ponte_Vazio_2021_MCE_RRA",
           "Fundação_Vazio_2021_MCE_Reads",
           "Fundação_Vazio_2021_MCE_ASVs",
           "Fundação_Vazio_2021_MCE_OTUs",
           "Fundação_Vazio_2021_MCE_RRA",
           #2021 Sterivex
           "Prainha_Vazio_2021_Sterivex_Reads",
           "Prainha_Vazio_2021_Sterivex_ASVs",
           "Prainha_Vazio_2021_Sterivex_OTUs",
           "Prainha_Vazio_2021_Sterivex_RRA",
           "Barragem_Vazio_2021_Sterivex_Reads",
           "Barragem_Vazio_2021_Sterivex_ASVs",
           "Barragem_Vazio_2021_Sterivex_OTUs",
           "Barragem_Vazio_2021_Sterivex_RRA",
           "Ponte_Vazio_2021_Sterivex_Reads",
           "Ponte_Vazio_2021_Sterivex_ASVs",
           "Ponte_Vazio_2021_Sterivex_OTUs",
           "Ponte_Vazio_2021_Sterivex_RRA",
           "Fundação_Vazio_2021_Sterivex_Reads",
           "Fundação_Vazio_2021_Sterivex_ASVs",
           "Fundação_Vazio_2021_Sterivex_OTUs",
           "Fundação_Vazio_2021_Sterivex_RRA",
           #2022 MCE
           "Prainha_Cheio_2022_MCE_Reads",
           "Prainha_Cheio_2022_MCE_ASVs",
           "Prainha_Cheio_2022_MCE_OTUs",
           "Prainha_Cheio_2022_MCE_RRA",
           "Barragem_Cheio_2022_MCE_Reads",
           "Barragem_Cheio_2022_MCE_ASVs",
           "Barragem_Cheio_2022_MCE_OTUs",
           "Barragem_Cheio_2022_MCE_RRA",
           "Ponte_Cheio_2022_MCE_Reads",
           "Ponte_Cheio_2022_MCE_ASVs",
           "Ponte_Cheio_2022_MCE_OTUs",
           "Ponte_Cheio_2022_MCE_RRA",
           "Fundação_Cheio_2022_MCE_Reads",
           "Fundação_Cheio_2022_MCE_ASVs",
           "Fundação_Cheio_2022_MCE_OTUs",
           "Fundação_Cheio_2022_MCE_RRA"
  ) 

wider_dt_all_resume %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

View(wider_dt_all_resume)

write.csv(wider_dt_all_resume, paste0(tbl_path, "/wider_dt_all_resume_2024.csv"))

## Tabela suplementar 2 ----

# Tabela feita com os dados filtrados, apenas peixes e reads >=100

t_grouped_filt <- t_grouped_by_ID_BLASTid %>%
  group_by(`Curated ID`,
           # `Final ID (BLASTn)`,
           `new_name`, 
           `level`,
           `filter`,
           `year`,
           Sample
  ) %>%
  reframe(
          "Sample" =  Sample,
          "Curated ID" = `Curated ID`,
          "Final ID (BLASTn)" = `Final ID (BLASTn)`,
          "BLASTn pseudo-score" = `BLASTn pseudo-score`,
          "ASV (Sequence)" = `ASV (Sequence)`,
          "Class (BLASTn)" = `Class (BLASTn)`,
          "Curated Order (BLASTn)" = `Curated Order (BLASTn)`,
          "Family (BLASTn)" = `Family (BLASTn)`,
          "Genus (BLASTn)" = `Genus (BLASTn)`,
          "Curated genus" = `Curated genus`,
          "expedition" = expedition,
          "new_name" = new_name,
          "year" = year,
          "filter" = filter,
          "ASV header" = `ASV header`,
          # "RRA" = RRA, # Nao posso usar o RRA ja que o total de read por amostra foi alterado com os filtros
          "OTU" = OTU,
          "Abundance" = Abundance,
          "Primer expected length" = `Primer expected length`,
          "level" = level,
          "ASVs" = length(unique(`ASV header`)),
          "OTUs" = length(unique(OTU)),
          "Reads" = sum(Abundance)) %>% 
  filter(`Class (BLASTn)` %in% "Actinopteri") %>% 
  filter(Reads >= 100) %>% 
  group_by(Sample) %>%
  mutate(RRA = Reads/sum(Reads)) %>% 
  ungroup() %>% 
  unique()

View(t_grouped_filt)

# Teste para ver se o RRA foi calculado corretamente

t_grouped_filt %>%
  filter(level == "Low") %>%
  filter(new_name == "Fundação") %>%
  filter(filter == "Sterivex") %>%
  filter(year == "2021") %>% 
  pull(RRA) %>%
  sum() 

# 29/08: RRA Nao foi calculado corretamente pois sairam reads com os filtros 
# Realizar a correcao.

# Tabela longer
t_dt_all_filt <- t_grouped_filt %>%
  select(-c("Final ID (BLASTn)")) %>% # Retirando a identificacao original
  group_by(level, new_name, filter, year) %>%
  mutate("Abd total" = sum(Abundance)) %>%
  ungroup() %>%
  group_by(`Curated ID`) %>% 
  mutate("Abundancia total" = sum(Abundance)) %>% 
  ungroup %>% 
  # group_by(`Final ID (BLASTn)`) %>%  # Retirando as colunas relacionadas as identificacoes originais
  # mutate("Reads totais" = sum(Abundance)) %>% # Retirando as colunas relacionadas as identificacoes originais
  # ungroup() %>% # Retirando as colunas relacionadas as identificacoes originais
  group_by(`Curated ID`,
           # `Final ID (BLASTn)`,
           `new_name`, 
           # `Abundance`,
           `level`,
           `filter`,
           `year`
  ) %>%
  mutate("RRA no periodo" = (Abundance/`Abd total`)*100) %>% 
  # ungroup() %>% 
  reframe(
    # "Reads totais" = `Reads totais`,
    "Abundancia total" = `Abundancia total`,
    "RRA" = sum(`RRA no periodo`),
    "ASVs" = length(unique(`ASV header`)),
    "OTUs" = length(unique(OTU)),
    "level" = level,
    "Ano" = year,
    "Filtro" = filter,
    # "Ponto" = new_name,
    "Reads" = sum(Abundance),
    "Classe" = `Class (BLASTn)`,
    "Order" = `Curated Order (BLASTn)`,
    "Family" = `Family (BLASTn)`,
    "Genus" = `Curated genus`
    # ,
    # "ASV header" = `ASV header`,
    # "OTU" = OTU
  ) %>% 
  mutate(
    "RRA_formatado" = format(RRA, scientific = TRUE)
  ) %>%
  unique() 

View(t_dt_all_filt)

# options(scipen = 99,
#         digits = 3)

# Teste de RRA para ver se foi calculado corretamente
t_dt_all_filt %>%
  filter(level == "Vazio") %>%
  filter(new_name == "Fundação") %>%
  filter(Filtro == "MCE") %>%
  filter(Ano == "2021") %>% 
  pull(RRA) %>%
  sum()  

# Tabela wider com os dados filtrados
t_wider_dt_filt <- t_dt_all_filt %>% 
  select(-c("RRA_formatado")) %>%
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  unite(new_name, level, Ano, Filtro, col= "ponto_level_ano_filtro") %>% 
  pivot_wider(id_cols = c("Classe",
                          "Curated ID", 
                          # "Final ID (BLASTn)",
                          # "Reads totais", 
                          "Abundancia total"),
              names_from = ponto_level_ano_filtro,
              values_from =  c("RRA","ASVs","OTUs","Reads"),
              names_glue = "{ponto_level_ano_filtro}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate("Curated ID",
           # "Final ID (BLASTn)",
           # "Reads totais", 
           "Abundancia total",
           #2020 MCE
           "Ponte_Cheio_2020_MCE_Reads", 
           "Ponte_Cheio_2020_MCE_ASVs", 
           "Ponte_Cheio_2020_MCE_OTUs",
           "Ponte_Cheio_2020_MCE_RRA",
           "Fundação_Cheio_2020_MCE_Reads",
           "Fundação_Cheio_2020_MCE_ASVs",
           "Fundação_Cheio_2020_MCE_OTUs",
           "Fundação_Cheio_2020_MCE_RRA",
           #2021 MCE
           "Prainha_Vazio_2021_MCE_Reads",
           "Prainha_Vazio_2021_MCE_ASVs",
           "Prainha_Vazio_2021_MCE_OTUs",
           "Prainha_Vazio_2021_MCE_RRA",
           "Barragem_Vazio_2021_MCE_Reads",
           "Barragem_Vazio_2021_MCE_ASVs",
           "Barragem_Vazio_2021_MCE_OTUs",
           "Barragem_Vazio_2021_MCE_RRA",
           "Ponte_Vazio_2021_MCE_Reads",
           "Ponte_Vazio_2021_MCE_ASVs",
           "Ponte_Vazio_2021_MCE_OTUs",
           "Ponte_Vazio_2021_MCE_RRA",
           "Fundação_Vazio_2021_MCE_Reads",
           "Fundação_Vazio_2021_MCE_ASVs",
           "Fundação_Vazio_2021_MCE_OTUs",
           "Fundação_Vazio_2021_MCE_RRA",
           #2021 Sterivex
           "Prainha_Vazio_2021_Sterivex_Reads",
           "Prainha_Vazio_2021_Sterivex_ASVs",
           "Prainha_Vazio_2021_Sterivex_OTUs",
           "Prainha_Vazio_2021_Sterivex_RRA",
           "Barragem_Vazio_2021_Sterivex_Reads",
           "Barragem_Vazio_2021_Sterivex_ASVs",
           "Barragem_Vazio_2021_Sterivex_OTUs",
           "Barragem_Vazio_2021_Sterivex_RRA",
           "Ponte_Vazio_2021_Sterivex_Reads",
           "Ponte_Vazio_2021_Sterivex_ASVs",
           "Ponte_Vazio_2021_Sterivex_OTUs",
           "Ponte_Vazio_2021_Sterivex_RRA",
           "Fundação_Vazio_2021_Sterivex_Reads",
           "Fundação_Vazio_2021_Sterivex_ASVs",
           "Fundação_Vazio_2021_Sterivex_OTUs",
           "Fundação_Vazio_2021_Sterivex_RRA",
           #2022 MCE
           "Prainha_Cheio_2022_MCE_Reads",
           "Prainha_Cheio_2022_MCE_ASVs",
           "Prainha_Cheio_2022_MCE_OTUs",
           "Prainha_Cheio_2022_MCE_RRA",
           "Barragem_Cheio_2022_MCE_Reads",
           "Barragem_Cheio_2022_MCE_ASVs",
           "Barragem_Cheio_2022_MCE_OTUs",
           "Barragem_Cheio_2022_MCE_RRA"
           # "Ponte_Cheio_2022_MCE_Reads", # Os dois prox pontos nao sobreviveram 
           # "Ponte_Cheio_2022_MCE_ASVs", # aos filtros por que nao tem ASVs de peixes 
           # "Ponte_Cheio_2022_MCE_OTUs",
           # "Ponte_Cheio_2022_MCE_RRA",
           # "Fundação_Cheio_2022_MCE_Reads",
           # "Fundação_Cheio_2022_MCE_ASVs",
           # "Fundação_Cheio_2022_MCE_OTUs",
           # "Fundação_Cheio_2022_MCE_RRA"
  ) 

wider_dt_filt %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

View(wider_dt_filt)

write.csv(wider_dt_filt, paste0(tbl_path, "/wider_dt_filt_2024.csv"))


## Raw informations ----

# Raw reads

# Run 2

# $ zgrep -c "^@M" ~/projetos/peixes-eDNA/raw/2020_reads/*.gz | awk -F: '{sum += $2} END {print sum}'
# 140506

# Runs 4, 5
# $ grep -c '^@' ~/projetos/peixes-eDNA/raw/combined/*.fastq.gz | awk -F: '{sum += $2} END {print sum}'
# 5528758

# Run EM156
#$ zgrep -c "^@VH" ~/projetos/peixes-eDNA/raw/ecomol/*.gz | awk -F: '{sum += $2} END {print sum}'
# 4922134

# Post-pipeline reads

# verificando o total de reads por amostra apos filtragem

filt_raw_reads <- t_filt_results_tbl %>% group_by(Sample) %>% 
  summarise(total_abd_filt = sum(Abundance))

View(t_filt_results_tbl %>% group_by(Sample) %>% 
       summarise(total_abd_filt = sum(Abundance)))

View(filt_raw_reads)
sum(filt_raw_reads$total_abd_filt)

# Taxa detected Tabela Suplementar 2
{
  ## Agora levando em consideracao apenas os dados filtrados,
  ## apenas peixes e reads >=100
  
  order_by_tax <- grouped_filt %>%
    # dt_filt_resume %>% 
    group_by(
      # `Final ID (BLASTn)`
      `Curated ID`
    ) %>%
    reframe(
      "Ordem" = `Curated Order (BLASTn)`,
      "Família" = `Family (BLASTn)`,
      "Gênero" = `Curated genus`,
      # "Maior hit do BLASTn" = `Final ID (BLASTn)`,
      "Identificação curada" = `Curated ID`,
      "Reads" = sum(Reads), # dt_all_resume && # dt_filt
      "ASVs" = sum(ASVs),
      "OTUs" = sum(OTUs),
    ) %>% 
    unique() %>% 
    select("Ordem",
           "Família",
           "Gênero",
           # "Maior hit do BLASTn",
           "Identificação curada",
           "Reads",
           "ASVs",
           "OTUs"
    )
  
  View(order_by_tax)
  
  write.csv(order_by_tax, "~/projetos/LI_paper/results/tables/order_by_tax.csv")
}

# Reads, ASVs, OTUs grouped by Order
{
  order_level <- grouped_filt %>%
    group_by(
      `Curated Order (BLASTn)`
    ) %>%
    reframe(
      # "Curated ID" = `Curated ID`,
      "Ordem" = `Curated Order (BLASTn)`,
      "Total reads" = sum(Reads), # dt_all_resume && # dt_filt
      # "Classe" = Classe,
      "ASVs" = sum(ASVs),
      "OTUs" = sum(OTUs),
      "Famílias" = length(unique(`Family (BLASTn)`)),
      "Gêneros" = length(unique(`Genus (BLASTn)`)),
      "Ids" = length(unique(`Curated ID`)),
      "Espécies" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])),
      "Nspecies" = length(unique(`Curated ID`)) - Espécies
    ) %>% 
    select(-c(`Curated Order (BLASTn)`)) %>% 
    unique() 
  
  View(order_level)
  
  write.csv(order_by_tax, "~/projetos/LI_paper/results/tables/order_level.csv")
  
}

# Exploratory tables and plots ----

# Grouped by Class 
{

analis_dt_class <-
  # t_grouped_by_ID_BLASTid %>%
  # t_filt_results_tbl %>%
  dt_all_resume %>%
  # dt_filt %>%
  group_by(
    Classe
    # `Curated ID`,
    # new_name,
    # level
  ) %>%
  reframe(
    # "Total reads" = sum(Abundance), # t_grouped_by_ID_BLASTid && # t_filt_results_tbl
    "Total reads" = sum(Reads), # dt_all_resume && # dt_filt
    "Classe" = Classe,
    "ASVs" = sum(ASVs),
    "OTUs" = sum(OTUs),
    "Ordens" = length(unique(Order)),
    "Famílias" = length(unique(Family)),
    "Gêneros" = length(unique(Genus)),
    "Espécies" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])), # IDs a level de spp
    "Nspecies" = length(unique(`Curated ID`)) - Espécies # IDs em outros niveis taxonomicos
  ) %>% 
  unique()

# View(analis_pre_eco)
View(analis_dt_class)

# Verificando ASVs que nao foram identificadas a level de spp (genus + epitetus)
{
  spp_level <- dt_all_resume %>%
    group_by(
      Classe
    ) %>% 
    reframe(
      "Espécies" = unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])
    ) %>% 
    pull(Espécies)
  
  nspp_level <- dt_all_resume %>% 
    group_by(
      Classe
    ) %>% 
    reframe(
      "Nspecies" = setdiff(unique(`Curated ID`), spp_level)
    )
  
  View(nspp_level) 
}

}

# Grouped by Order 
{

grouped_by_Order <- dt_all_resume %>%
  filter(Classe %in% "Actinopteri") %>%
  group_by(
    Order
    # Family
  ) %>%
  reframe(
    # "Curated ID" = `Curated ID`,
    # "Total reads" = sum(Reads), # dt_all_resume && # dt_filt
    # "Classe" = Classe,
    # "Ordem" = Order,
    "Famílias" = length(unique(Family)),
    "Gêneros" = length(unique(Genus)),
    "Ids" = length(unique(`Curated ID`)),
    "Espécies" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])),
    "Nspecies" = length(unique(`Curated ID`)) - `Espécies`,
    "ASVs" = sum(ASVs),
    "OTUs" = sum(OTUs)
  ) %>% 
  unique() 

View(grouped_by_Order)
}

# Grouped by Sample
{
grouped_by_Sample <- grouped_filt %>% 
  filter(expedition %in% "Novembro 2021") %>% 
  select(Sample, filter, Abundance, new_name ) %>% 
  group_by(filter, 
           Sample
  ) %>% 
  mutate("Total reads" = sum(Abundance)) %>%
  reframe(filter,
          # Sample,
          new_name,
          `Total reads`) %>% 
  unique()

View(grouped_by_Sample)
}

# Grouped by Family
{
grouped_by_Family <- grouped_filt %>% 
  unite(new_name, year, filter, sep = "_", remove = FALSE, col = "un_amostral") %>% 
  select(c(`Family (BLASTn)`, 
           `Curated Order (BLASTn)`,
           `Curated ID`, new_name, year, filter, un_amostral)) %>% 
  unique() %>%
  group_by(new_name, year, filter) %>% 
  mutate("alpha_d" = length(`Curated ID`)) %>%
  ungroup() %>% 
  group_by(`Family (BLASTn)`, new_name, year, filter) %>% 
  mutate("alpha_d_order" = length(`Curated ID`)) %>% 
  ungroup() %>% 
  mutate("alpha_d_order_%" = round(alpha_d_order/alpha_d*100, digits = 2)) %>% 
  select(-c(`Curated ID`)) %>% 
  unique() %>% 
  group_by(`Family (BLASTn)`) %>% 
  mutate("alpha_total" = sum(`alpha_d_order_%`)) %>%
  ungroup() %>% 
  group_by(`Family (BLASTn)`
           , `Curated Order (BLASTn)`
  ) %>% 
  reframe(alpha_total)

View(grouped_by_Family)
}

# Especies encontradas em apenas uma amostra
{
counts <- grouped_filt %>% 
  group_by(`Curated ID`) %>%
  # group_by(`Curated ID`, filter) %>%
  summarise(count = n_distinct(Sample),
            filter, Sample) %>% 
  unique()

unique_counts <- 
  counts %>% filter(count == 1) %>% 
  reframe(`Curated ID`)

View(unique_counts)

# Stats da Final ID

stats_final_ID <- grouped_filt %>% 
  mutate("Genus" = str_split_fixed(string = .$`Final ID (BLASTn)`, 
                                   pattern = " ",
                                   n = 2)[,1]) %>% 
  select(`Final ID (BLASTn)`, Genus)
}

# Proporcao de IDs do LGC12sDB e do NT/NCBI
{
# Usando o objeto hits_DB criado no script 012024_post_pipe.r

filt_hits_DB <- grouped_filt %>% 
  left_join(hits_DB,
            by = "ASV header") %>% 
  select(`ASV header`, `Final ID (BLASTn).x`, `Final ID (BLASTn).y`, `Curated ID`, `1_DB`, `2_DB`, `3_DB`) %>% 
  unique()

dif_hits <- hits_DB %>%
  anti_join(filt_hits_DB) 

# pre final hits_DB
group_DB <- filt_hits_DB %>% 
  select(`Final ID (BLASTn).x`,`Curated ID`, 
         `1_DB`
         # ,`2_DB`,`3_DB`
  ) %>% 
  unique()

# group_DB %>% 
#   write_xlsx(path = "/home/gabriel/projetos/lagoa_ingleses/results/tabelas/group_DB.xlsx", 
#              col_names = TRUE,
#              format_headers = TRUE)

## Upload curated_group_DB

curated_group_DB  <- read_excel("/home/gabriel/projetos/lagoa_ingleses/results/tabelas/curated_group_DB.xlsx")

# Grafico de barras 
curated_group_DB %>% 
  select(`Curated ID`,curated_1_DB) %>% 
  group_by(curated_1_DB) %>% 
  count() %>% 
  mutate(perc = n / 65 * 100) %>% 
  ggplot(aes(x = curated_1_DB, 
             y = perc,
             fill = curated_1_DB)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = seq(0, 60, 10)) +
  theme(
    panel.grid.major = element_line(color = "grey",
                                    size = 0.2,
                                    linetype = 1),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", size = rel(1.2)),
    plot.title = element_text(color = "black", size = rel(1.5)),
    plot.subtitle = element_text(color = "black", size = rel(1.2)),
  ) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.position = "none") +
  scale_fill_manual(values = c("#399283",
                               "#c0e087","#104b6d")) +
  labs(x = "Banco de dados",
       y = "Proporção %",
       title = "Origem das identificações")

# Boxplot

index_df <- grouped_filt %>%
  filter(expedition %in% c("Novembro 2021")) %>% 
  # group_by(Sample,level) %>% 
  group_by(Sample,filter) %>% 
  summarise("Species" = length(unique(`Curated ID`))) %>%
  # mutate(level = ifelse(level == "Cheio", "Full", "Empty")) %>%
  # mutate(level = factor(level, levels = c("Full", "Empty"))) %>% 
  as.data.frame()

# compare_means(Species ~ level, data=index_df, method = "t.test", paired = FALSE)
compare_means(Species ~ filter, data=index_df, method = "t.test", paired = FALSE)

boxplot_h <- grouped_filt %>% 
  filter(expedition %in% c("Novembro 2021")) %>% 
  # group_by(Sample,level) %>% 
  group_by(Sample,filter) %>% 
  reframe("Species" = length(unique(`Curated ID`))) %>%
  # mutate(level = ifelse(level == "Cheio", "Full", "Empty")) %>%
  # mutate(level = factor(level, levels = c("Full", "Empty"))) %>%
  ggplot(aes(x =
               # level,
               filter,
             color=
               # level,
               filter,
             y = 
               Species,
             alpha
  )) +
  # geom_boxplot(
  # fill = "blue",
  # color = "black",
  # size = 0.3) +
  geom_jitter(height = 0.2, width = 0.1)+ 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.major = element_line(color = "grey",
                                    size = 0.2,
                                    linetype = 1),
    # axis.ticks = element_line(color = "grey"),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", size = rel(1.2)),
    # legend.background = element_blank(),
    # legend.key = element_blank(),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = rel(1.2)),
    plot.title = element_text(color = "black", size = rel(1.5)),
    plot.subtitle = element_text(color = "black", size = rel(1.2)),
    strip.background = element_rect(fill = "#e4e4e4"),
    strip.text = element_text(color = "black", size = rel(1.2))) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 14, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
  ) +
  scale_color_manual(values = c("#ed9f0e",
                                "#3381b1")) +
  labs(x = "Filter",
       y = "Richness",
       fill ='Filter',
       # title = "Riqueza de espécies",
       # subtitle = "Filtro MCE versus Sterivex"
  ) +
  # P-value no grafico
  stat_compare_means(label.x = 0.60)

boxplot_h

# Salvar em pdf
ggsave(plot =  boxplot_h, 
       filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                        "boxplot_MCE_Sterivex", "-", Sys.Date(), ".pdf", sep = ""),
       device = "pdf",
       units = "cm",
       height = 20,
       width = 10,
       dpi = 600) 



# Proportions of identified orders 
{
  fish_ord_all <- grouped_filt %>% 
    unite(new_name, year, filter,sep = "_", remove = FALSE, col = "un_amostral") %>% 
    select(c(`Curated Order (BLASTn)`,`Curated ID`, new_name, year, filter, un_amostral)) %>% 
    unique() %>%
    group_by(new_name, year, filter) %>% 
    mutate("alpha_d" = length(`Curated ID`)) %>%
    ungroup() %>% 
    group_by(`Curated Order (BLASTn)`,new_name, year, filter) %>% 
    mutate("alpha_d_order" = length(`Curated ID`)) %>% 
    ungroup() %>% 
    mutate("alpha_d_order_%" = round(alpha_d_order/alpha_d*100, digits = 2)) %>% 
    select(-c(`Curated ID`)) %>% 
    unique() %>% 
    group_by(un_amostral) %>% 
    mutate(ano_level = case_when(
      year == "2020" ~ "2020 - Full",
      year == "2021" ~ "2021 - Low",
      year == "2022" ~ "2022 - Full"
    )) %>%
    filter(filter %in% "MCE") %>% 
    ggplot(aes(x = new_name,
               y = alpha_d_order,
               # y = length(`Curated ID`),
               # y = RRA,
               fill = `Curated Order (BLASTn)`)) +
    geom_col(position = "fill") +
    geom_text(aes( # adicionar % dentro das barras
      label = sprintf("%0.1f%%",`alpha_d_order_%`)),
      position = position_fill(vjust = 0.5),
      color = "black") +
    facet_grid(rows = vars(ano_level
                           # , 
                           # filter
    ),
    space = "free", 
    scales = "free",
    drop = TRUE) +
    coord_flip() + #virando 90o o grafico
    guides(fill = guide_legend("Orders")) + ## alterar o titulo da legenda 
    theme(
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
    ) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 19),
          strip.text = element_text(size = 15, face = "bold"),
          legend.position = "bottom",
          legend.background = element_rect(fill = "lightgray", color = NA)) +
    scale_y_continuous(labels = scales::percent) +
    # scale_colour_brewer(palette = "Oranges") +
    scale_fill_manual(values = c(
                                 # "#e41a1c",
                                 # "#377eb8",
                                 # "#4daf4a",
                                 # "#984ea3",
                                 # "#ff7f00"
                                 # "#ff7f00",
                                 # "#a64595",
                                 # "#34466d",
                                 # "#429369",
                                 # "#881c23",
                                 # "#c0710c",
                                 # "#7e80af",
                                 # "#195036",
                                 # "#c5667d",
                                 # "#cf6449",
                                 # "#53348e",
                                 "#5d9cecff",
                                 "#4fc1e9ff",
                                 "#ac92ecff",
                                 "#a0d468ff",
                                 "#ffce54ff",
                                 "#fc6e51ff",
                                 "#ed5565ff",
                                 "#48cfadff",                                  
                                 "#ec87c0ff"
                                 )) +
    labs(
         # title = "Proporções de ordens identificadas",
         # subtitle = "Todas amostras",
         x = "Sample sites",
         y = "Proportion") 
  
  fish_ord_all
  
  ggsave(plot = fish_ord_all, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2024/",
                          "alpha_ord_all_PER", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf", units = "cm", height = 20, width = 40, dpi = 600)
}

# Proportions of identified families
{
  fish_fam_all <- grouped_filt %>%
    unite(new_name, year, filter,sep = "_", remove = FALSE, col = "un_amostral") %>% 
    select(c(`Curated Order (BLASTn)`, `Family (BLASTn)`,`Curated ID`, new_name, year, filter, un_amostral)) %>% 
    unique() %>%
    group_by(new_name, year, filter) %>% 
    mutate("alpha_d" = length(`Curated ID`)) %>%
    ungroup() %>% 
    group_by(`Family (BLASTn)`,new_name, year, filter) %>% 
    mutate("alpha_d_order" = length(`Curated ID`)) %>% 
    ungroup() %>% 
    mutate("alpha_d_order_%" = round(alpha_d_order/alpha_d*100, digits = 2)) %>% 
    select(-c(`Curated ID`)) %>% 
    unique() %>% 
    group_by(un_amostral) %>% 
    mutate(ano_level = case_when(
      year == "2020" ~ "2020 - Full",
      year == "2021" ~ "2021 - Low",
      year == "2022" ~ "2022 - Full"
    )) %>%
    filter(filter %in% "MCE") %>% 
    ggplot(aes(x = new_name,
               y = alpha_d_order,
               fill = `Family (BLASTn)`)) +
    geom_col(position = "fill") +
    geom_text(aes( #adicionando % dentro das barras
      label = sprintf("%0.1f%%",`alpha_d_order_%`)),
      position = position_fill(vjust = 0.5),
      color = "black") +
    facet_grid(rows = vars(ano_level
                           # , filter
    ),
    space = "free", 
    scales = "free",
    drop = TRUE) +
    coord_flip() + #virando 90o o grafico
    guides(fill = guide_legend("Family")) + ## alterar o titulo da legenda 
    theme(
      # panel.background = element_blank(),
      # panel.grid.major = element_blank(),
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
    ) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 19),
          strip.text = element_text(size = 15, face = "bold"),
          legend.position = "bottom",
          legend.background = element_rect(fill = "lightgray", color = NA)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(
                                # "#256676",
                                #  "#41c9dc",
                                #  "#d60724",
                                #  "#a6e590",
                                #  "#69306e",
                                #  "#84ee15",
                                #  "#531ce8",
                                #  "#64903a",
                                #  "#c637bc",
                                #  "#34f199",
                                #  "#873c1a",
                                #  "#f4cacb",
                                #  "#0b5313",
                                #  "#68affc",
                                #  "#1642cd",
                                #  "#ffa8ff",
                                #  "#21a708",
                                #  "#9e73b8",
                                #  "#f3d426",
                                #  "#cc7b6f",
                                #  "#fea53b",
                                #  "#4d4815"
                                 # "#ff595e",
                                 # "#ff924c",
                                 # "#ffca3a",
                                 # "#c5ca30",
                                 # "#8ac926",
                                 # "#52a675",
                                 # "#1982c4",
                                 # "#4267ac",
                                 # "#6a4c93",
                                 # "#e2e2dfff",
                                 # "#d2d2cfff",
                                 # "#e2cfc4ff",
                                 # "#f7d9c4ff",
                                 # "#faedcbff",
                                 # "#c9e4deff",
                                 # "#c6def1ff",
                                 # "#dbcdf0ff",
                                 # "#f2c6deff",
                                 # "#f9c6c9ff",
                                 # "#c6def1ff",
                                 # "#dbcdf0ff",
                                 # "#c8ddbbff",
                                 # "#b6e2ddff",
                                 # "#e9e5afff",
                                 # "#fbdf9dff",
                                 # "#fbc99dff",
                                 # "#fbb39dff",
                                 # "#fba09dff"
                                  "#5d9cecff",
                                  "#4fc1e9ff",
                                  "#48cfadff",
                                  "#a0d468ff",
                                  "#ffce54ff",
                                  "#fc6e51ff",
                                  "#ed5565ff",
                                  "#ac92ecff",
                                  "#ec87c0ff"
                                  ),
                      # breaks = "Family (BLASTn)" +
                        # name = "Curated Order (BLASTn)"
                      ) +
    labs(
      # title = "Proporções de famílias identificadas",
      # subtitle = "Todas amostras",
      x = "Sample sites",
      y = "Proportion"
         )
  
  fish_fam_all
  
  ggsave(plot = fish_fam_all, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2024/",
                          "alpha_fam_all_PER", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf", units = "cm", height = 21, width = 40, dpi = 600)}

# Tile plots 
{
  # Defining levels 
  {
    # Defining order levels
    order_levels <- grouped_filt$`Curated Order (BLASTn)` %>% 
      sort() %>% unique()
    
    # Defining spp levels
    spp_levels <- grouped_filt$`Curated ID` %>% 
      sort() %>% unique()
    
    # Defining sample levels
    grouped_filt$Sample %>% unique()
    
    samples <- c("L1_nov_dec_20_mi",
                 "L2_dez20",
                 "L2_nov20",
                 "L1_out21",
                 "L2_out21",
                 "L3_out21",
                 "L4_out21",
                 "L1_nov21",
                 "L2_nov21",
                 "L3_nov21",
                 "L4_nov21",
                 "STX_L1_nov21",
                 "STX_L2_nov21",
                 "STX_L3_nov21",
                 "STX_L4_nov21",
                 "L1_jan22",
                 "L2_jan22",
                 "L3_jan22",
                 "L4_jan22"
    )
    
    # Defining expedition levels
    grouped_filt$expedition %>% unique()
    
    expeditions <- c("Novembro e Dezembro 2020",
                     "Novembro 2020",
                     "Dezembro 2020",
                     "Outubro 2021",
                     "Novembro 2021",
                     "Janeiro 2022")
    
  }
  
  # Tile plots 
  {
    # All years (only MCE filters) 
    tile_all <- grouped_filt %>%
      # filter(RRA >= 0.01) %>%
      mutate(`Curated ID` = factor(`Curated ID`, levels = rev(spp_levels))) %>%
      mutate(expedition = factor(expedition, levels = expeditions)) %>%
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(new_name = factor(new_name)) %>% 
      # mutate(sample = factor(sample)) %>% 
      mutate(`Curated Order (BLASTn)` = factor(`Curated Order (BLASTn)`)) %>%
      mutate(ano_level = case_when(
        year == "2020" ~ "2020 - Full",
        year == "2021" ~ "2021 - Empty",
        year == "2022" ~ "2022 - Full"
      )) %>%
      mutate(RRA = RRA*100) %>% 
      group_by(expedition, `Curated ID`, new_name, `Curated Order (BLASTn)`, level) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(expedition %in% c("Novembro 2021")) %>%
      filter(filter %in% c("MCE")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`,
                 x = new_name,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # stat = "identity",
      # colour = "black", size = 4) +
      facet_grid(cols = vars(ano_level),
                 rows = vars(`Curated Order (BLASTn)`),
                 space = "free", 
                 scales = "free",
                 drop = TRUE) +
      # scale_fill_continuous(
      #   trans = "log10",
      #   breaks = c(0.000001, 0.0087, 0.87),  # Defina os pontos de quebra desejados
      #   labels = c("0.001", "0.01", "100"),  # Rótulos correspondentes
      #   type = "viridis"
      # ) +
      scale_fill_gradientn(name = "RRA (%)",
                           # colours = c("white", "yellow", "green", "darkgreen", "blue"),
                           colours = rev(c("#30123BFF", "#4662D7FF", "darkgreen", "#72FE5EFF", "#C7EF34FF",
                                           "#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")),
                           # colours = c("white", "#FDE725FF", "#6DCD59FF", "#35B779FF", "#3E4A89FF"),
                           # colours = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                           #             "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                           #             "#B4DE2CFF", "#FDE725FF"),
                           values = scales::rescale(c(0, 0.01, 0.05, 0.25, 1, 2.5, 5, 10, 25, 50)),
                           breaks = c(0, 0.01, 0.05, 0.25, 1.00, 2.50, 5.00, 10.00, 25.00, 50.00), # Quebras exatas
                           labels = function(x) scales::number(x, accuracy = 0.01, trim = FALSE), # Precisão exata
                           limits = c(0.01, 50),
                           na.value = "#7A0403FF",
                           trans = "log10"
      ) +
      scale_colour_gradientn(name = "RRA (%)",
                             # colours = rev(c("white", "yellow", "green", "darkgreen", "blue")),
                             colours = rev(c("#30123BFF", "#4662D7FF", "darkgreen", "#72FE5EFF", "#C7EF34FF",
                                             "#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")),
                             # colours = rev(c("white", "#FDE72C5FF", "#6DCD59FF", "#35B779FF", "#3E4A89FF")),
                             # colours = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                             #             "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                             #             "#B4DE2CFF", "#FDE725FF"),
                             values = scales::rescale(c(0, 0.01, 0.05, 0.25, 1, 2.5, 5, 10, 25, 50)),
                             breaks = c(0, 0.01, 0.05, 0.25, 1.00, 2.50, 5.00, 10.00, 25.00, 50.00), # Quebras exatas
                             labels = function(x) scales::number(x, accuracy = 0.01, trim = FALSE), # Precisão exata
                             limits = c(0.01, 50),
                             na.value = "#7A0403FF",
                             trans = "log10") +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1),
        # axis.ticks = element_line(color = "grey"),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = rel(1.2)),
        # legend.background = element_blank(),
        # legend.key = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black", size = rel(1.2)),
        plot.title = element_text(color = "black", size = rel(1.5)),
        plot.subtitle = element_text(color = "black", size = rel(1.2)),
        strip.background = element_rect(fill = "#e4e4e4"),
        strip.text = element_text(color = "black", size = rel(1.2))) +
      theme(plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(size = 14, angle = 20, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 14, face = "italic"),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text = element_text(size = 14, face = "bold"),
            strip.text.y = element_text(angle=0),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
            legend.position = "right",
            # legend.key.width = unit(3, 'cm'),
            legend.key.height = unit(3, 'cm')
      ) +
      labs(fill ='RRA (%)',
           # title = "Detected species",
           # subtitle = "Use of MCE filter",
           x = "Sample sites",
           y = "Species"
      )
    
    tile_all
    
    ggsave(plot = tile_all,
           filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                            "tile_20-22_MCE", "-", Sys.Date(), ".pdf", sep = ""),
           device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)
    
    
    # MCE versus Sterivex
    tile_mce_vs_sterivex <- grouped_filt %>%
      # filter(RRA >= 0.01) %>%
      mutate(`Curated ID` = factor(`Curated ID`, levels = rev(spp_levels))) %>%
      mutate(expedition = factor(expedition, levels = expeditions)) %>%
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(new_name = factor(new_name)) %>% 
      # mutate(sample = factor(sample)) %>% 
      mutate(`Curated Order (BLASTn)` = factor(`Curated Order (BLASTn)`)) %>%
      mutate(ano_level = case_when(
        year == "2020" ~ "2020 - Cheio",
        year == "2021" ~ "2021 - Vazio",
        year == "2022" ~ "2022 - Cheio"
      )) %>%
      mutate(RRA = RRA*100) %>%
      group_by(expedition, `Curated ID`, new_name, `Curated Order (BLASTn)`, level) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(year %in% c(2021)) %>%
      filter(expedition %in% c("Novembro 2021")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`,
                 x = new_name,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # stat = "identity",
      # colour = "black", size = 4) +
      facet_grid(cols = vars(filter),
                 rows = vars(`Curated Order (BLASTn)`),
                 space = "free", 
                 scales = "free",
                 drop = TRUE) +
      # scale_fill_continuous(
      #   trans = "log10",
      #   breaks = c(0.000001, 0.0087, 0.87),  # Defina os pontos de quebra desejados
      #   labels = c("0.001", "0.01", "100"),  # Rótulos correspondentes
      #   type = "viridis"
      # ) +
      scale_fill_gradientn(name = "RRA (%)",
                           colours = rev(c("#30123BFF", "#4662D7FF", "darkgreen", "#72FE5EFF", "#C7EF34FF",
                                           "#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")),
                           # colours = c("white", "yellow", "green", "darkgreen", "purple", "red", "blue"),
                           # colours = c("white", "yellow", "green", "darkgreen", "purple","pink", "blue"),
                           # colours = c("white", "#FDE725FF", "#6DCD59FF", "#35B779FF", "#3E4A89FF"),
                           # colours = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                           #             "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                           #             "#B4DE2CFF", "#FDE725FF"),
                           values = scales::rescale(c(0, 0.01, 0.05, 0.25, 1, 2.5, 5, 10, 25, 50)),
                           breaks = c(0, 0.01, 0.05, 0.25, 1.00, 2.50, 5.00, 10.00, 25.00, 50.00), # Quebras exatas
                           labels = function(x) scales::number(x, accuracy = 0.01, trim = FALSE), # Precisão exata
                           limits = c(0.01, 50),
                           na.value = "#7A0403FF",
                           trans = "log10"
      ) +
      scale_colour_gradientn(name = "RRA (%)",
                             colours = rev(c("#30123BFF", "#4662D7FF", "darkgreen", "#72FE5EFF", "#C7EF34FF",
                                             "#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")),
                             # colours = rev(c("white", "yellow", "green", "darkgreen", "purple", "red", "blue")),
                             # colours = rev(c("white", "yellow", "green", "darkgreen", "purple","pink", "blue")),
                             # colours = rev(c("white", "#FDE72C5FF", "#6DCD59FF", "#35B779FF", "#3E4A89FF")),
                             # colours = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                             #             "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                             #             "#B4DE2CFF", "#FDE725FF"),
                             values = scales::rescale(c(0, 0.01, 0.05, 0.25, 1, 2.5, 5, 10, 25, 50)),
                             breaks = c(0, 0.01, 0.05, 0.25, 1.00, 2.50, 5.00, 10.00, 25.00, 50.00), # Quebras exatas
                             labels = function(x) scales::number(x, accuracy = 0.01, trim = FALSE), # Precisão exata
                             limits = c(0.01, 50),
                             na.value = "#7A0403FF",
                             trans = "log10") +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1),
        # axis.ticks = element_line(color = "grey"),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = rel(1.2)),
        # legend.background = element_blank(),
        # legend.key = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black", size = rel(1.2)),
        plot.title = element_text(color = "black", size = rel(1.5)),
        plot.subtitle = element_text(color = "black", size = rel(1.2)),
        strip.background = element_rect(fill = "#e4e4e4"),
        strip.text = element_text(color = "black", size = rel(1.2))) +
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14, face = "italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
        # legend.key.width = unit(3, 'cm'),
        legend.key.height = unit(3, 'cm')
      ) +
      labs(fill = "RRA%",
           # title = "Espécies detectadas",
           # subtitle = "Filtro MCE versus Sterivex",
           x = "Sample sites",
           y = "Species")
    
    tile_mce_vs_sterivex
    
    ggsave(plot = tile_mce_vs_sterivex,
           filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                            "tile_MCE_Sterivex", "-", Sys.Date(), ".pdf", sep = ""),
           device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)
    
    
    # Empty versus Full (SBG)
    tile_sbg_cheio_vs_vazio <- grouped_filt %>%
      # filter(RRA >= 0.01) %>%
      filter(filter %in% c("MCE")) %>%
      mutate(`Curated ID` = factor(`Curated ID`, levels = rev(spp_levels))) %>% 
      mutate(expedition = factor(expedition, levels = expeditions)) %>%
      mutate(Sample = factor(Sample, levels = samples)) %>%
      mutate(new_name = factor(new_name)) %>%
      # mutate(sample = factor(sample)) %>%
      mutate(`Curated Order (BLASTn)` = factor(`Curated Order (BLASTn)`)) %>%
      mutate(RRA = RRA*100) %>%
      group_by(expedition, `Curated ID`, new_name, `Curated Order (BLASTn)`, level) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(new_name %in% c("Ponte", "Fundacao")) %>%
      ggplot(aes(y = `Curated ID`,
                 group = `Curated ID`,
                 x = new_name,
                 # x = Sample,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # stat = "identity",
      # colour = "black", size = 4) +
      facet_grid(cols = vars(level),
                 labeller = labeller(level = c("Cheio" = "Full", "Vazio" = "Low")),
                 rows = vars(`Curated Order (BLASTn)`),
                 space = "free",
                 scales = "free",
                 drop = TRUE) +
      # scale_fill_continuous(
      #   trans = "log10",
      #   breaks = c(0.000001, 0.0087, 0.87),  # Defina os pontos de quebra desejados
      #   labels = c("0.001", "0.01", "100"),  # Rótulos correspondentes
      #   type = "viridis"
      # ) +
      scale_fill_gradientn(name = "RRA (%)",
                           colours = rev(c("#30123BFF", "#4662D7FF", "darkgreen", "#72FE5EFF", "#C7EF34FF",
                                           "#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")),
                           # colours = c("white", "yellow", "green", "darkgreen", "purple", "red", "blue"),
                           # colours = c("white", "yellow", "green", "darkgreen", "purple","pink", "blue"),
                           # colours = c("white", "#FDE725FF", "#6DCD59FF", "#35B779FF", "#3E4A89FF"),
                           # colours = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                           #             "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                           #             "#B4DE2CFF", "#FDE725FF"),
                           values = scales::rescale(c(0, 0.01, 0.05, 0.25, 1, 2.5, 5, 10, 25, 50)),
                           breaks = c(0, 0.01, 0.05, 0.25, 1.00, 2.50, 5.00, 10.00, 25.00, 50.00), # Quebras exatas
                           labels = function(x) scales::number(x, accuracy = 0.01, trim = FALSE), # Precisão exata
                           limits = c(0.01, 50),
                           na.value = "#7A0403FF",
                           trans = "log10"
      ) +
      scale_colour_gradientn(name = "RRA (%)",
                             colours = rev(c("#30123BFF", "#4662D7FF", "darkgreen", "#72FE5EFF", "#C7EF34FF",
                                             "#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")),
                             # colours = rev(c("white", "yellow", "green", "darkgreen", "purple", "red", "blue")),
                             # colours = rev(c("white", "yellow", "green", "darkgreen", "purple","pink", "blue")),
                             # colours = rev(c("white", "#FDE72C5FF", "#6DCD59FF", "#35B779FF", "#3E4A89FF")),
                             # colours = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                             #             "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                             #             "#B4DE2CFF", "#FDE725FF"),
                             values = scales::rescale(c(0, 0.01, 0.05, 0.25, 1, 2.5, 5, 10, 25, 50)),
                             breaks = c(0, 0.01, 0.05, 0.25, 1.00, 2.50, 5.00, 10.00, 25.00, 50.00), # Quebras exatas
                             labels = function(x) scales::number(x, accuracy = 0.01, trim = FALSE), # Precisão exata
                             limits = c(0.01, 50),
                             na.value = "#7A0403FF",
                             trans = "log10") +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1),
        # axis.ticks = element_line(color = "grey"),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = rel(1.2)),
        # legend.background = element_blank(),
        # legend.key = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black", size = rel(1.2)),
        plot.title = element_text(color = "black", size = rel(1.5)),
        plot.subtitle = element_text(color = "black", size = rel(1.2)),
        strip.background = element_rect(fill = "#e4e4e4"),
        strip.text = element_text(color = "black", size = rel(1.2))) +
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14, face = "italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
        # legend.key.width = unit(3, 'cm'),
        legend.key.height = unit(3, 'cm')
      ) +
      labs(fill = "RRA%",
           # title = "Espécies detectadas",
           # subtitle = "Filtro MCE versus Sterivex",
           x = "Sample sites",
           y = "Species")
      
    tile_sbg_cheio_vs_vazio
    
     ggsave(plot = tile_sbg_cheio_vs_vazio, 
           filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                            "tile_sbg_cheio_vs_vazio", "-", Sys.Date(), ".pdf", sep = ""),
           device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)
  }
}

# Venn diagram 
{
  # Filtering
  resume_venn <- grouped_filt %>%
    # filter(expedition %in% c("Novembro 2021")) %>%
    filter(filter %in% c("MCE")) %>%
    select(`Curated ID`, `Curated genus`, `Family (BLASTn)`, `Curated Order (BLASTn)`, filter, level) %>%
    unique()
  
  # 1a. Common spp between MCE and Sterivex
  spp_comm <- resume_venn %>% 
    group_by(`Curated ID`) %>% # Trocar para spp, genero, familia e ordem
    summarise(n_filtros = n_distinct(filter)) %>%
    filter(n_filtros >= 2) %>%
    select(`Curated ID`) # Trocar para spp, genero, familia e ordem
  View(spp_comm)
  
  # 1b. Exclusive spp MCE and Sterivex
  spp_diff <- resume_venn %>%
    group_by(`Curated ID`) %>% # Trocar para spp, genero, familia e ordem
    summarise(MCE_count = sum(as.integer(filter == "MCE")),
              STX_count = sum(as.integer(filter == "Sterivex")),
              MCE_only = MCE_count >= 1 & STX_count == 0,
              STX_only = MCE_count == 0 & STX_count >= 1) %>% 
    ungroup()
  View(spp_diff)
  
  # 2a. Common spp between Full and Low
  spp_comm_2 <- resume_venn %>% 
    group_by(`Curated ID`) %>% # Trocar para spp, genero, familia e ordem
    summarise(n_level = n_distinct(level)) %>%
    filter(n_level >= 2) %>%
    select(`Curated ID`) # Trocar para spp, genero, familia e ordem
  View(spp_comm_2)
  
  # 2b. Exclusive spp Full and Low
  spp_diff_2 <- resume_venn %>%
    group_by(`Curated ID`) %>% # Trocar para spp, genero, familia e ordem
    summarise(Low_count = sum(as.integer(level == "Vazio")),
              High_count = sum(as.integer(level == "Cheio")),
              Low_only = Low_count >= 1 & High_count == 0,
              High_only = Low_count == 0 & High_count >= 1) %>% 
    ungroup()
  View(spp_diff_2)
}
}

# Analysis ----

## NMDS ----

# NMDS (Non-Metric Multidimensional Scaling), ou Escalonamento Multidimensional Não 
# Métrico, é uma técnica de análise multivariada utilizada para visualizar a semelhança 
# ou dissimilaridade entre objetos ou observações em um conjunto de dados. O NMDS é 
# particularmente útil quando os dados não possuem uma estrutura linear clara e quando 
# a distância entre os pontos não pode ser representada de forma métrica.

# Converter planilha de identificações por ponto para o formato Amostras X IDs
{
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
  
  # Número de linhas (IDs diferentes) por amostra
  grouped_filt$Sample %>% table()
  
  grouped_filt %>% colnames()
  
  # FINAL_tbl_IDs <- t_grouped_by_ID_BLASTid %>%
  FINAL_tbl_IDs <- grouped_filt %>%
    # filter(expedition %in% c("Novembro 2021")) %>% 
    # filter(RRA >= 0.05) %>% # definicao de qual threshold de abundancia sera usado
    # filter(New_name %in% c("Ponte", "Fundacao")) %>% 
    # filter(!Sample %in% c("L2_dez20")) %>% # amostra problematica
    select(c("Sample", "Curated ID", "year", "level", "new_name", "filter", "RRA", "expedition")) %>% 
    pivot_wider(id_cols = c("Sample", "year", "level", "new_name", "filter", "expedition"),
                names_from = `Curated ID`,
                values_from = `RRA`,
                values_fn = sum_uniq,
                names_sort = TRUE,
                names_prefix = "ID_") %>% 
    relocate(c("Sample", "year", "level", "filter", "new_name", "expedition", starts_with("ID_"))) %>%  
    mutate(across(starts_with("ID_"), replace_na, replace = 0)) %>% 
    mutate("new_name" = as.factor(new_name)) %>% 
    mutate("level" = as.factor(level)) 
  
  # Verificando a soma das linhas
  FINAL_tbl_IDs %>% select(starts_with(match = "ID_")) %>% rowSums(na.rm = TRUE)
}

#  Rodando o NMDS

# Preparar dados para entrada no pacote vegan
colnames(FINAL_tbl_IDs)

FINAL_tbl_IDs$`Sample` %>% unique() %>% sort()

all_IDs_NMDS_tbl <- FINAL_tbl_IDs %>% 
  mutate("Sample number" = 0) %>% 
  relocate("Sample number" )

# Associar números de amostra aos nomes de amostra
for (sample in 1:nrow(all_IDs_NMDS_tbl)) {
  all_IDs_NMDS_tbl$`Sample number`[sample] <- sample
}

# Ordenar dataframe usado no NMDS
all_IDs_NMDS_df <- all_IDs_NMDS_tbl %>%
  # mutate(level = ifelse(level == "Cheio", "Full", "Empty")) %>%
  as.data.frame()

# Nomear linhas como números de amostra
row.names(all_IDs_NMDS_df) <- all_IDs_NMDS_df$`Sample number`
all_IDs_NMDS_df <- all_IDs_NMDS_df 

# Corrigir nomes das espécies para evitar problemas na plotagem
colnames(all_IDs_NMDS_df)
colnames(all_IDs_NMDS_df)[8:ncol(all_IDs_NMDS_df)] <- colnames(all_IDs_NMDS_df)[8:ncol(all_IDs_NMDS_df)] %>%
  str_replace_all(pattern = " ", replacement = "_") %>% 
  str_replace_all(pattern = "\\.", replacement = "") %>% 
  str_replace_all(pattern = "\\(", replacement = "") %>% 
  str_replace_all(pattern = "\\)", replacement = "")

# Executando o NMDS
# Esta e a funcao que faz o NMDS. Para ela fornece-se apenas
# as colunas relativas as especies nas amostras
all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df[,8:ncol(all_IDs_NMDS_df)],
                                 comm = all_IDs_NMDS_df[8:ncol(all_IDs_NMDS_df)],
                                 # distance = "bray", # Leva em conta o RRA
                                 distance = "jaccard", # Considera apenas presença/ausenciaa
                                 halfchange = TRUE 
)

dim(all_IDs_NMDS_df)

# Fazer o fit das variáveis ambientais

# meta.envfit <- envfit(all_ps_vegan_ord_meta, all_IDs_NMDS_df[, c("level", "Sample")],
meta.envfit <- envfit(all_ps_vegan_ord_meta, all_IDs_NMDS_df[, c("filter", "Sample")],
                      permutations = 999,
                      na.rm=TRUE) # this fits environmental vectors

# Espécies 

# Fazer o fit das espécies para identificar a significância delas na explicação dos agrupamentos. 
# Esse é o passo que mais demora quando com muitas amostras e espécies.
meta.spp.fit <- envfit(all_ps_vegan_ord_meta, 
                       all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)], 
                       permutations = 999) # this fits species vectors

# Obter valores de p para as espécies
sps_pvals <- tibble("IDs" = names(meta.spp.fit$vectors$pvals),
                    "p-value" = meta.spp.fit$vectors$pvals)

spp.scrs <- as.data.frame(scores(meta.spp.fit, display = "vectors")) %>%
  mutate("IDs" = rownames(.)) %>%
  left_join(y = sps_pvals, by = "IDs")

{
  spp.scrs$IDs <- gsub("ID_", "", spp.scrs$IDs)
  spp.scrs$IDs <- gsub("_", " ", spp.scrs$IDs)
  } 

# Selecionar espécies significativas
sig.spp.scrs <- spp.scrs 
# %>%
#   filter(`p-value` <=
#            0.05) # definir p-value aqui!

# Pontos amostrais

# Definir os valores de NMDS1 e NMDS2, e os metadados de cada amostra/ponto amostral
site.scrs <- as.data.frame(scores(all_ps_vegan_ord_meta, display = "sites")) %>%
  mutate("Sample number" = as.double(row.names(.))) %>%
  left_join(y = all_IDs_NMDS_df[, c("Sample number",
                                    "Sample",
                                    "level",
                                    "new_name",
                                    "filter"
  )], 
  by = "Sample number")

# Determinar centroides
scrs <- scores(all_ps_vegan_ord_meta, display = "sites")

cent <- aggregate(scrs ~ filter, data = site.scrs, FUN = "mean")

# Calcular elipses
NMDS <- data.frame("MDS1" = all_ps_vegan_ord_meta$points[, 1],
                   "MDS2" = all_ps_vegan_ord_meta$points[, 2],
                   "filter" = as.factor(all_IDs_NMDS_df$filter), check.names = FALSE)

NMDS.mean <- aggregate(NMDS[, 1:2], list(group = NMDS$filter), "mean")


# Elipses
{
  plot(all_ps_vegan_ord_meta)
  
  # Sobrepor as elipses
  ord <- ordiellipse(ord = all_ps_vegan_ord_meta, 
                     groups = all_IDs_NMDS_df$filter,
                     display = "sites",
                     kind = "ehull", conf = 0.95, label = T)
  
  
  # Funcao do vegan de calcular elipses
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ell <- data.frame()
  
  for(g in levels(NMDS$filter)){
    print(g)
    df_ell <- 
      rbind(df_ell, 
            cbind(as.data.frame(with(NMDS[NMDS$filter==g,],
                                     veganCovEllipse(
                                       ord[[g]]$cov,
                                       ord[[g]]$center,
                                       ord[[g]]$scale))),
                  filter=g))}
}

# Configurar plot NMDS

# Definir paleta de cores
{
  # paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")
  # paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[1:3]
  # paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[c(1:3, 2,7,9,10)]
  # 
  # my_cols <- c("BAM" = "#0F6B99FF",
  #              "PARAC" = "#6B990FFF",
  #              "SAM" = "#99540FFF",
  #              "SFC" = "#A3CC51FF",
  #              "SFI" = "#7EC3E5FF",
  #              "SFM" = "#E5B17EFF") 
  
  cores <- c("#f1c232",
             "#3381b1")
  
}

# Plot
{
  NMDS_MCE_Sterivex <-
    ggplot(data = site.scrs,
           aes(x=NMDS1, 
               y=NMDS2)) +
    
    # Elipses 
    ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                               data = df_ell,
                               aes(x = NMDS1,
                                   y = NMDS2,
                                   group = filter,
                                   label = filter,
                                   col = filter,
                                   fill = filter
                               ),
                               alpha=0.10,
                               # n = 200,
                               linetype=2,
                               expand = 0,
                               label.fontsize = 14,
                               con.cap = 0.1
    ) +
    # Niveis amostrais
    geom_point(aes(x=NMDS1,
                   y=NMDS2,
                   fill = filter,
                   # label = `Sampling unit`,  #descomentar se quiser exibir os nomes dos `Sampling sites`s
                   col=filter,
                   group = filter,
                   shape = new_name),
               stroke = 0.5,
               alpha = 0.75,
               size = 3) +
    
    # Nomes dos `Sampling sites`s amostrais
    ggrepel::geom_text_repel(aes(label = new_name),  #descomentar bloco se quiser exibir os nomes dos pontos
                             # hjust=0.5,
                             # vjust=2.75,
                             size=4.5,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha=0.1,
                             # min.segment.length = 1.5,
                             force = 3,
                             max.overlaps = 100,
                             fontface = "bold") +
    
    # Vetores das IDs
    geom_segment(data = sig.spp.scrs, aes(x = 0,
                                          xend = NMDS1,
                                          y = 0,
                                          yend = NMDS2),
                 arrow = arrow(length = unit(0.1, "cm")),
                 colour = "grey30",
                 alpha = 0.5,
                 lwd = 0.3) + #add vector arrows of significant species
    
    # Nomes das IDs
    ggrepel::geom_text_repel(data = sig.spp.scrs,
                             aes(x=NMDS1, y=NMDS2, label = IDs),
                             size = 3.5,
                             alpha = 0.75,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha = 0.1,
                             max.overlaps = 100) +
    
    # Centroides
    geom_point(data = cent,
               aes(x = NMDS1,
                   y = NMDS2, # colour = filter,
                   fill = filter
               ),
               size = 8,
               colour = "#222222",
               alpha = 0.75,
               shape  = 13
    ) + 
    coord_fixed(expand = c(0.5))+
    theme_light() +
    theme(
      panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      panel.grid.major = element_line(color = "grey",
                                      size = 0.2,
                                      linetype = 1),
      # axis.ticks = element_line(color = "grey"),
      axis.ticks = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = rel(1.2)),
      # legend.background = element_blank(),
      # legend.key = element_blank(),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black", size = rel(1.2)),
      plot.title = element_text(color = "black", size = rel(1.5)),
      plot.subtitle = element_text(color = "black", size = rel(1.2))) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          # legend.position = "bottom"
    ) +
    labs(x = "NMDS1",
         y = "NMDS2",
         title = "Composição da ictiofauna",
         # subtitle = "Filtro MCE versus Sterivex") +
         subtitle = "Stress: ",format(round(all_ps_vegan_ord_meta$stress))) +
    scale_fill_manual(values = cores) +
    scale_colour_manual(values = cores) 
  # + scale_fill_manual(values = viridis::turbo(n = 5)) # cores das formas das samples/pontos
  
  
  NMDS_MCE_Sterivex
  
  
  # Salvar em pdf 
  ggsave(plot = NMDS_MCE_Sterivex, 
         filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                          "NDMS_MCE_STX_filt", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 25,
         width = 25,
         dpi = 600)
}

## PCoA ----

# Dados
## Usando a tabela all_IDs_NMDS_df gerada para o NDMDS

## Executando o PCoA para MCE versus Sterivex

# Filtrando os dados

mce_stx_IDs_NMDS_df <- all_IDs_NMDS_df %>% 
  filter(expedition %in% c("Novembro 2021")) %>% 
  select(where(~ any(. != 0))) 

colnames(mce_stx_IDs_NMDS_df[,8:ncol(mce_stx_IDs_NMDS_df)])

#Criando a Matriz de distancia
pcOa_dist <- vegan::vegdist(x = mce_stx_IDs_NMDS_df[,8:ncol(mce_stx_IDs_NMDS_df)],
                            method = "jaccard",
                            binary = TRUE)

pcOa <- cmdscale(pcOa_dist, eig = TRUE)

ordiplot(pcOa, display = 'sites', type = 'text')

# Espécies 

# Fazer o fit das espécies para identificar a significância delas na explicação dos agrupamentos. 
# Esse é o passo que mais demora quando com muitas amostras e espécies.
meta.spp.fit <- envfit(pcOa, 
                       mce_stx_IDs_NMDS_df[,8:ncol(mce_stx_IDs_NMDS_df)], 
                       permutations = 999) # this fits species vectors

# Obter valores de p para as espécies
sps_pvals <- tibble("IDs" = names(meta.spp.fit$vectors$pvals),
                    "p-value" = meta.spp.fit$vectors$pvals)

spp.scrs <- as.data.frame(scores(meta.spp.fit, display = "vectors")) %>%
  mutate("IDs" = rownames(.)) %>%
  left_join(y = sps_pvals, by = "IDs")

{
  spp.scrs$IDs <- gsub("ID_", "", spp.scrs$IDs)
  spp.scrs$IDs <- gsub("_", " ", spp.scrs$IDs)
  } 

# Selecionar espécies significativas
sig.spp.scrs <- spp.scrs 
# %>%
#   filter(`p-value` <=
#            0.05) # definir p-value aqui!

# Pontos amostrais

#Definir os valores de PCoA1 e PCoA2, e os metadados de cada amostra/ponto amostral
site.scrs <- as.data.frame(scores(pcOa, display = "sites"))

colnames(site.scrs) <- c("PCoA1", "PCoA2")

site.scrs <- site.scrs %>% 
  mutate("Sample number" = as.double(row.names(.))) %>%
  left_join(y = mce_stx_IDs_NMDS_df[, c("Sample number",
                                    "Sample",
                                    "level",
                                    "new_name",
                                    "filter",
                                    "expedition"
  )],
  by = "Sample number")

# Determinar centroides
scrs <- scores(pcOa, display = "sites")

cent <- aggregate(scrs ~ filter, data = site.scrs, FUN = "mean")

# Calcular elipses
PCoA <- data.frame("PCoA1" = pcOa$points[, 1],
                   "PCoA2" = pcOa$points[, 2],
                   "filter" = as.factor(mce_stx_IDs_NMDS_df$filter), check.names = FALSE)

PCoA.mean <- aggregate(PCoA[, 1:2], list(group = PCoA$filter), "mean")

# Elipses

# plot(all_ps_vegan_ord_meta)
ordiplot(pcOa, display = 'sites', type = 'text')

# Sobrepor as elipses
ord <- ordiellipse(ord = pcOa, 
                   groups = mce_stx_IDs_NMDS_df$filter,
                   display = "sites",
                   kind = "ehull", conf = 0.95, label = T)


# Funcao do vegan de calcular elipses
veganCovEllipse <-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()

for(g in levels(PCoA$filter)){
  print(g)
  df_ell <- 
    rbind(df_ell, 
          cbind(as.data.frame(with(PCoA[PCoA$filter==g,],
                                   veganCovEllipse(
                                     ord[[g]]$cov,
                                     ord[[g]]$center,
                                     ord[[g]]$scale))),
                filter=g))}

# Configurar plot PCoA

# Definir paleta de cores

cores <- c("#f1c232",
           "#3381b1")

# Plot
{
  PCoA_MCE_Sterivex <-
    ggplot(data = site.scrs,
           aes(x=PCoA1, 
               y=PCoA2)) +
    
    # Elipses 
    ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                               data = df_ell,
                               aes(x = Dim1,
                                   y = Dim2,
                                   group = filter,
                                   label = filter,
                                   col = filter,
                                   fill = filter
                               ),
                               alpha=0.10,
                               # n = 200,
                               linetype=2,
                               expand = 0,
                               label.fontsize = 14,
                               con.cap = 0.1
    ) +
    # Niveis amostrais
    geom_point(aes(x=PCoA1,
                   y=PCoA2,
                   fill = filter,
                   # label = `Sampling unit`,  #descomentar se quiser exibir os nomes dos `Sampling sites`s
                   col=filter,
                   group = filter,
                   shape = new_name),
               stroke = 0.5,
               alpha = 0.75,
               size = 3) +
    
  # Nomes dos Sampling sites amostrais
    ggrepel::geom_text_repel(aes(label = new_name),  #descomentar bloco se quiser exibir os nomes dos pontos
                             # hjust=0.5,
                             # vjust=2.75,
                             size=4.5,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha=0.1,
                             # min.segment.length = 1.5,
                             force = 3,
                             max.overlaps = 100,
                             fontface = "bold") +

  # Vetores das IDs
  geom_segment(data = sig.spp.scrs, aes(x = 0,
                                        xend = Dim1,
                                        y = 0,
                                        yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "grey30",
               alpha = 0.5,
               lwd = 0.3) + #add vector arrows of significant species

  # Nomes das IDs
  ggrepel::geom_text_repel(data = sig.spp.scrs,
                           aes(x=Dim1, y=Dim2, label = IDs),
                           size = 3.5,
                           alpha = 0.75,
                           direction = "both",
                           segment.size = 0.25,
                           segment.alpha = 0.1,
                           max.overlaps = 25) +
  
  # Centroides
  geom_point(data = cent,
             aes(x = Dim1,
                 y = Dim2, # colour = filter,
                 fill = filter
             ),
             size = 8,
             colour = "#222222",
             alpha = 0.75,
             shape  = 13
  ) + 
    coord_fixed(expand = c(0.5))+
    theme_light() +
    theme(
      panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      panel.grid.major = element_line(color = "grey",
                                      size = 0.2,
                                      linetype = 1),
      # axis.ticks = element_line(color = "grey"),
      axis.ticks = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = rel(1.2)),
      # legend.background = element_blank(),
      # legend.key = element_blank(),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black", size = rel(1.2)),
      plot.title = element_text(color = "black", size = rel(1.5)),
      plot.subtitle = element_text(color = "black", size = rel(1.2))) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          # legend.position = "bottom"
    ) +
    labs(x = "PCoA1 53.34%",
         y = "PCoA2 18,77%",
         # title = "Composição da ictiofauna",
         # subtitle = "Filtro MCE versus Sterivex") +
         # subtitle = "Stress: ",format(round(all_ps_vegan_ord_meta$stress))
    ) +
    scale_fill_manual(values = cores) +
    scale_colour_manual(values = cores) 
  # + scale_fill_manual(values = viridis::turbo(n = 5)) # cores das formas das samples/pontos
  
  PCoA_MCE_Sterivex
}

# Salvar em pdf 
ggsave(plot = PCoA_MCE_Sterivex, 
       filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                        "PCoA_MCE_Sterivex_clean",
                        # "PCoA_MCE_Sterivex_full",
                        "-", Sys.Date(), ".pdf", sep = ""),
       device = "pdf",
       units = "cm",
       height = 25,
       width = 25,
       dpi = 600)

## Executando o PCoA para Full versus Low

# Filtrando os dados

mce_IDs_NMDS_df <- all_IDs_NMDS_df %>% 
  filter(filter %in% c("MCE")) %>% 
  select(where(~ any(. != 0))) 

colnames(mce_IDs_NMDS_df[,8:ncol(mce_IDs_NMDS_df)])

#Criando a Matriz de distancia
pcOa_dist <- vegan::vegdist(x = mce_IDs_NMDS_df[,8:ncol(mce_IDs_NMDS_df)],
                            method = "jaccard",
                            binary = TRUE)

pcOa <- cmdscale(pcOa_dist, eig = TRUE)

ordiplot(pcOa, display = 'sites', type = 'text')

# Espécies 

# Fazer o fit das espécies para identificar a significância delas na explicação dos agrupamentos. 
# Esse é o passo que mais demora quando com muitas amostras e espécies.
meta.spp.fit <- envfit(pcOa, 
                       mce_IDs_NMDS_df[,8:ncol(mce_IDs_NMDS_df)], 
                       permutations = 999) # this fits species vectors

# Obter valores de p para as espécies
sps_pvals <- tibble("IDs" = names(meta.spp.fit$vectors$pvals),
                    "p-value" = meta.spp.fit$vectors$pvals)

spp.scrs <- as.data.frame(scores(meta.spp.fit, display = "vectors")) %>%
  mutate("IDs" = rownames(.)) %>%
  left_join(y = sps_pvals, by = "IDs")

{
  spp.scrs$IDs <- gsub("ID_", "", spp.scrs$IDs)
  spp.scrs$IDs <- gsub("_", " ", spp.scrs$IDs)
  } 

# Selecionar espécies significativas
sig.spp.scrs <- spp.scrs 
# %>%
#   filter(`p-value` <=
#            0.05) # definir p-value aqui!

# Pontos amostrais

#Definir os valores de PCoA1 e PCoA2, e os metadados de cada amostra/ponto amostral
site.scrs <- as.data.frame(scores(pcOa, display = "sites"))

colnames(site.scrs) <- c("PCoA1", "PCoA2")

site.scrs <- site.scrs %>% 
  mutate("Sample number" = as.double(row.names(.))) %>%
  left_join(y = mce_IDs_NMDS_df[, c("Sample number",
                                    "Sample",
                                    "level",
                                    "new_name",
                                    # "filter",
                                    "expedition"
  )],
  by = "Sample number")

# Determinar centroides
scrs <- scores(pcOa, display = "sites")

cent <- aggregate(scrs ~ level, data = site.scrs, FUN = "mean")

# Calcular elipses
PCoA <- data.frame("PCoA1" = pcOa$points[, 1],
                   "PCoA2" = pcOa$points[, 2],
                   # "filter" = as.factor(mce_IDs_NMDS_df$filter), check.names = FALSE)
                   "level" = as.factor(mce_IDs_NMDS_df$level), check.names = FALSE)

PCoA.mean <- aggregate(PCoA[, 1:2], list(group = PCoA$level), "mean")

# Elipses

# plot(all_ps_vegan_ord_meta)
ordiplot(pcOa, display = 'sites', type = 'text')

# Sobrepor as elipses
ord <- ordiellipse(ord = pcOa, 
                   groups = mce_IDs_NMDS_df$level,
                   display = "sites",
                   kind = "ehull", conf = 0.95, label = T)


# Funcao do vegan de calcular elipses
veganCovEllipse <-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()

for(g in levels(PCoA$level)){
  print(g)
  df_ell <- 
    rbind(df_ell, 
          cbind(as.data.frame(with(PCoA[PCoA$level==g,],
                                   veganCovEllipse(
                                     ord[[g]]$cov,
                                     ord[[g]]$center,
                                     ord[[g]]$scale))),
                level=g))}

# Configurar plot PCoA

# Definir paleta de cores

cores <- c("#3381b1",
           "#f1c232")

# Plot
{
  PCoA_MCE <-
    ggplot(data = site.scrs,
           aes(x=PCoA1, 
               y=PCoA2)) +
    
    # Elipses 
    ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                               data = df_ell,
                               aes(x = Dim1,
                                   y = Dim2,
                                   group = level,
                                   label = level,
                                   col = level,
                                   fill = level
                               ),
                               alpha=0.10,
                               # n = 200,
                               linetype=2,
                               expand = 0,
                               label.fontsize = 14,
                               con.cap = 0.1
    ) +
    # Niveis amostrais
    geom_point(aes(x=PCoA1,
                   y=PCoA2,
                   fill = level,
                   # label = `Sampling unit`,  #descomentar se quiser exibir os nomes dos `Sampling sites`s
                   col=level,
                   group = level,
                   shape = new_name),
               stroke = 0.5,
               alpha = 0.75,
               size = 3) +
    
    # Nomes dos Sampling sites amostrais
    ggrepel::geom_text_repel(aes(label = new_name),  #descomentar bloco se quiser exibir os nomes dos pontos
                             # hjust=0.5,
                             # vjust=2.75,
                             size=4.5,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha=0.1,
                             # min.segment.length = 1.5,
                             force = 3,
                             max.overlaps = 100,
                             fontface = "bold") +
    
    # Vetores das IDs
    geom_segment(data = sig.spp.scrs, aes(x = 0,
                                          xend = Dim1,
                                          y = 0,
                                          yend = Dim2),
                 arrow = arrow(length = unit(0.1, "cm")),
                 colour = "grey30",
                 alpha = 0.5,
                 lwd = 0.3) + #add vector arrows of significant species
    
    # Nomes das IDs
    ggrepel::geom_text_repel(data = sig.spp.scrs,
                             aes(x=Dim1, y=Dim2, label = IDs),
                             size = 3.5,
                             alpha = 0.75,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha = 0.1,
                             max.overlaps = 25) +
    
    # Centroides
    geom_point(data = cent,
               aes(x = Dim1,
                   y = Dim2, # colour = filter,
                   fill = level
               ),
               size = 8,
               colour = "#222222",
               alpha = 0.75,
               shape  = 13
    ) + 
    coord_fixed(expand = c(0.5))+
    theme_light() +
    theme(
      panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      panel.grid.major = element_line(color = "grey",
                                      size = 0.2,
                                      linetype = 1),
      # axis.ticks = element_line(color = "grey"),
      axis.ticks = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = rel(1.2)),
      # legend.background = element_blank(),
      # legend.key = element_blank(),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black", size = rel(1.2)),
      plot.title = element_text(color = "black", size = rel(1.5)),
      plot.subtitle = element_text(color = "black", size = rel(1.2))) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          # legend.position = "bottom"
    ) +
    labs(x = "PCoA1 53.34%",
         y = "PCoA2 18,77%",
         # title = "Composição da ictiofauna",
         # subtitle = "Filtro MCE versus Sterivex") +
         # subtitle = "Stress: ",format(round(all_ps_vegan_ord_meta$stress))
    ) +
    scale_fill_manual(values = cores) +
    scale_colour_manual(values = cores) 
  # + scale_fill_manual(values = viridis::turbo(n = 5)) # cores das formas das samples/pontos
  
  PCoA_MCE
}

# Salvar em pdf 
ggsave(plot = PCoA_MCE, 
       filename = paste("/home/gabriel/projetos/LI_paper/results/figures/",
                        "PCoA_MCE_Sterivex_clean",
                        # "PCoA_MCE_Sterivex_full",
                        "-", Sys.Date(), ".pdf", sep = ""),
       device = "pdf",
       units = "cm",
       height = 25,
       width = 25,
       dpi = 600)

## Permanova ----

## Primeiro e necessario gerar uma matriz de distancias para plotar as amostras
## e em seguida testar a assumpcao de que os dados dispersos formam grupos bem definidos.

# Distance matrix

pcOa_mtx <-  as.data.frame(as.matrix(pcOa_dist)) # gerado anteriormente para realizar o plot PCoA

#Assumptions

Dispersion <- betadisper(pcOa_dist, group = all_IDs_NMDS_df$filter, type = "centroid") # testando a dispercao

plot(Dispersion)

anova(Dispersion)

#Test

perma_test <- adonis2(pcOa_mtx ~ as.factor(all_IDs_NMDS_df$filter), data = pcOa_mtx,
                      permutations=9999)
perma_test

## Correlation ----

# Selecionando as informacoes necessarias

View(grouped_filt) # Tabela principal
colnames(grouped_filt)

resume_tbl <- grouped_filt %>% 
  # select(Sample, `Curated ID`, filter, RRA) %>% 
  filter(expedition %in% c("Novembro 2021")) %>% 
  group_by(Sample, filter) %>% 
  mutate("Sample total abundance" = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(`Curated ID`, `filter`) %>% 
  reframe("ID rel. abd. on sample" =  sum(Abundance)/mean(`Sample total abundance`)) %>% 
  ungroup()

View(resume_tbl)

# Pivotar a tabela para realizar a correlacao

IDs_samples_abd <- clean_res_tbl_abd %>% 
  select(-c(`Unique sample name`)) %>% 
  pivot_wider(id_cols = c("Sample","curated final ID"), names_from = "Primer",values_from = "ID rel. abd. on sample")  %>% 
  mutate_if(is.numeric, replace_na, replace = 0)

pivot_IDs_abd <- resume_tbl %>%
  pivot_wider(id_cols = c("Curated ID"), 
              names_from = "filter", 
              values_from = "ID rel. abd. on sample") %>% 
  mutate_if(is.numeric, replace_na, replace = 0)

# Plotar as correlacoes

# MCE Vs. Sterivex ----
options(scipen=999)
corr_MCE_STX <- pivot_IDs_abd %>% 
  ggplot(aes(x = MCE,
             y = Sterivex
  )) +
  geom_point(aes()) +
  coord_fixed(xlim= c(-0.1,1.1),
              ylim= c(-0.1,1.1)) +
  geom_smooth(method = "lm",se = T) +
  scale_x_log10(oob = scales::squish_infinite) +
  scale_y_log10(oob = scales::squish_infinite) +
  stat_cor(method="pearson") 

corr_MCE_STX

# Salvar em pdf 
ggsave(plot = corr_MCE_STX, 
       filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2024/",
                        "corr_MCE_STX", "-", Sys.Date(), ".pdf", sep = ""),
       device = "pdf",
       units = "cm",
       height = 25,
       width = 15,
       dpi = 600)

### Species accumulation curves ----
  
grouped_filt %>% View()
  
  wider_grouped <- grouped_filt %>%
    select(c(Sample, 
             new_name,
             level,
             `Curated ID`, 
             Reads)) %>% 
    group_by(Sample, `Curated ID`, Reads) %>% 
    unique %>% 
    pivot_wider(names_from = `Curated ID`,
                values_from = Reads,
                names_sort = TRUE) 
  
  wider_grouped <- wider_grouped %>% as.data.frame()
    
  rownames(wider_grouped) <- wider_grouped$Sample
  
  wider_grouped <- wider_grouped %>% 
    select(-c(Sample))
  
  wider_grouped <- wider_grouped %>% 
    mutate(across(everything(), ~replace_na(.x, 0)))
  
  ## Function to transform values into 1 and zeros
  jaccarize <- function(x) {
    # Using ifelse to check each element of the vector x
    ifelse(x == 0, 0, 1)
  }
  
  wider_grouped_jac <- wider_grouped %>% 
    mutate(across(c(3:ncol(wider_grouped)), jaccarize))
  
   grouped_accum_1 <- specaccum(wider_grouped_jac[3:ncol(wider_grouped_jac)]
                              ,
                              method = "rarefaction"
                              # method = "random"
                              )
   
   grouped_accum_2 <- specaccum(wider_grouped_jac[3:ncol(wider_grouped_jac)]
                                ,
                                # method = "rarefaction"
                                method = "random"
   )
  
  plot(grouped_accum)
  
  
    all_IDs_NMDS_df %>% colnames() %>% duplicated()
  all_IDs_NMDS_df$ID_Apareiodon_sp2 %>% jspecarize()
  
  
  # transform proportions to presence/abcense 
  all_IDs_NMDS_df_jc <- all_IDs_NMDS_df %>%
    # as_tibble()
    mutate(across(starts_with('ID_'), jaccarize))
  
  
  #tirado daqui
  # https://vegandevs.github.io/vegan/reference/specaccum.html
  
  spec1 <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[,4:ncol(all_IDs_NMDS_df_jc)])
  spec2 <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[,4:ncol(all_IDs_NMDS_df_jc)],
                            method = "random")
  
  
  plot(spec1)
  
  
  plot(grouped_accum_1, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="lightgrey")
  boxplot(grouped_accum_2, col="#00CC32", add=TRUE, pch="*")
  
  
  grouped_accum_2 %>% str()
  grouped_accum_2 %>% summary()
  
  
  
  library(summarytools)
  
  
  ctable(grouped_accum_2) 
  
  
  
  tidy_specaccum <- function(x) {
    data.frame(
      site = x$sites,
      richness = x$richness,
      sd = x$sd)
  }
  
  grouped_accum_2$sites
  grouped_accum_2$richness
  grouped_accum_2$method
  grouped_accum_2$sd
  
  
  
  grouped_accum_2_tidy <- grouped_accum_2$perm %>% 
    reshape2::melt() %>% 
    rename("Sample" = "Var1",
           "Permutation" = "Var2",
           "Value" = "value") %>% 
    as_tibble()
  
  grouped_accum_2_tidy_rch <- tibble("Richness" = grouped_accum_1$richness,
                           "Sample" = grouped_accum_1$sites)
  
  
  collector_plot <- grouped_accum_2_tidy %>% ggplot(aes(x = Sample,y = Value,group=Sample))+
    # geom_point() +
    geom_boxplot(col ="#005602",
                 fill ="#005602",notch = TRUE,
                 alpha = 0.33,width = 0.5,outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.25,linetype=2) +
    geom_jitter(size=0.01,
                width = 0.25,
                # height = 0,
                col="#005602") +
    geom_line(data = grouped_accum_2_tidy_rch,
              inherit.aes = F,
              linewidth = 2.5,
              col ="#005602",
              alpha = 0.50,
              aes(x = Sample, y = Richness)) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", size = rel(0.25))) +
    ylab(label = "Richness") +
    xlab(label = "Sites/Samples") +
    scale_x_continuous(breaks = seq(1,16,1))
  
  
  # https://r-charts.com/ggplot2/grid/#google_vignette
  
  
  
  ggsave(file = paste0(figs_path,"/",
                       prjct_rad,"-",
                       "Coletor-",Sys.Date(),
                       # "-num_taxa_per_class-id80",
                       ".pdf",
                       collapse = ""),
         plot = collector_plot,
         device = "pdf",
         width = 16,
         height = 12,
         units = "cm",
         dpi = 300)
  
  collector_plot
  
  # spec1$freq
  # spec2$perm
  
  
  #creating a dataframe for ggplot2
  data <- data.frame("Sites"=grouped_accum_1$sites, "Richness1"=grouped_accum_1$richness,"SD1"=grouped_accum_1$sd,
                     "Richness2"=grouped_accum_1$richness, "SD2"=grouped_accum_1$sd
  )
  
  
  
  
  
  ggplot() +
    geom_point(data=data, aes(x=Sites, y=Richness1)) +
    geom_line(data=data, aes(x=Sites, y=Richness2)) +
    geom_ribbon(data=data ,aes(x=Sites, ymin=(Richness1-2*SD1),ymax=(Richness2+2*SD2)),alpha=0.2)
  
  
  
  specpool(x, pool, smallsample = TRUE)
  estimateR(x, ...)
  specpool2vect(X, index = c("jack1","jack2", "chao", "boot","Species"))
  poolaccum(x, permutations = 100, minsize = 3)
  estaccumR(x, permutations = 100, parallel = getOption("mc.cores"))
  ## S3 method for class 'poolaccum'
  summary(object, display, alpha = 0.05, ...)
  ## S3 method for class 'poolaccum'
  plot(x, alpha = 0.05, type = c("l","g"), ...)
  
  

### Collector curve per sample  ----
  
  all_IDs_NMDS_df %>% colnames() %>% duplicated()
  all_IDs_NMDS_df$ID_Apareiodon_sp2 %>% jspecarize()
  
  library(reshape2)
  library(ggrepel)
  
  
  #selecione 
  
  # P1
  P1_spec <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[c(1:2),4:ncol(all_IDs_NMDS_df_jc)])
  # P2
  P2_spec <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[c(3:4),4:ncol(all_IDs_NMDS_df_jc)])
  # P3
  P3_spec <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[c(5:6),4:ncol(all_IDs_NMDS_df_jc)])
  # P4
  P4_spec <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[c(7:8),4:ncol(all_IDs_NMDS_df_jc)])
  # P5
  P5_spec <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[c(9:10),4:ncol(all_IDs_NMDS_df_jc)])
  # P6
  P6_spec <- vegan::specaccum(comm = all_IDs_NMDS_df_jc[c(11:12),4:ncol(all_IDs_NMDS_df_jc)])
  
  
  P1_spec$sd
  P2_spec$sd
  P3_spec$sd
  P4_spec$sd
  P5_spec$sd
  P6_spec$sd
  
  
  P1_tbl <- P1_spec1$richness %>% melt(value.name = "Richness") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P1")
  P2_tbl <- P2_spec1$richness %>% melt(value.name = "Richness") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P2")
  P3_tbl <- P3_spec1$richness %>% melt(value.name = "Richness") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P3")
  P4_tbl <- P4_spec1$richness %>% melt(value.name = "Richness") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P4")
  P5_tbl <- P5_spec1$richness %>% melt(value.name = "Richness") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P5")
  P6_tbl <- P6_spec1$richness %>% melt(value.name = "Richness") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P6")
  
  
  spec_per_site <- bind_rows(P1_tbl, P2_tbl, P3_tbl, P4_tbl, P5_tbl, P6_tbl) %>% 
    mutate(Replicates = as.numeric(Replicates))
  
  
  P1_tbl_sd <- P1_spec1$sd %>% melt(value.name = "SD") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P1")
  P2_tbl_sd <- P2_spec1$sd %>% melt(value.name = "SD") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P2")
  P3_tbl_sd <- P3_spec1$sd %>% melt(value.name = "SD") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P3")
  P4_tbl_sd <- P4_spec1$sd %>% melt(value.name = "SD") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P4")
  P5_tbl_sd <- P5_spec1$sd %>% melt(value.name = "SD") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P5")
  P6_tbl_sd <- P6_spec1$sd %>% melt(value.name = "SD") %>% as_tibble(rownames = "Replicates") %>% mutate("Site" = "P6")
  
  
  spec_per_site_sd <- bind_rows(P1_tbl_sd, P2_tbl_sd, P3_tbl_sd, P4_tbl_sd, P5_tbl_sd, P6_tbl_sd) %>% 
    mutate(Replicates = as.numeric(Replicates))
  
  
  
  
  
  
  spec_per_site_complete <- spec_per_site %>% left_join(y = spec_per_site_sd,by = c("Site","Replicates"))
  
  # spec_per_site_complete$Replicates[spec_per_site_complete$Site == "P1" & spec_per_site_complete$Replicates == 1] <- 0.7
  # spec_per_site_complete$Replicates[spec_per_site_complete$Site == "P2" & spec_per_site_complete$Replicates == 1] <- 0.8
  # spec_per_site_complete$Replicates[spec_per_site_complete$Site == "P3" & spec_per_site_complete$Replicates == 1] <- 0.9
  # spec_per_site_complete$Replicates[spec_per_site_complete$Site == "P4" & spec_per_site_complete$Replicates == 1] <- 1
  # spec_per_site_complete$Replicates[spec_per_site_complete$Site == "P5" & spec_per_site_complete$Replicates == 1] <- 1.1
  # spec_per_site_complete$Replicates[spec_per_site_complete$Site == "P6" & spec_per_site_complete$Replicates == 1] <- 1.2
  
  
  
  labels_tbl <- spec_per_site_complete %>% filter(Replicates == 2)
  
  
  
  collector_plot_per_sample <- spec_per_site_complete %>%
    ggplot(aes(
      x = Replicates,
      y = Richness,
      group = Site,
      col = Site,
      fill = Site,alpha = 0.75)) +
    geom_line(linewidth = 2) +
    geom_label_repel(
      data = labels_tbl,
      inherit.aes = F,
      aes(
        x = Replicates,
        y = Richness,
        label = Site,
        col = Site,
        alpha = 0.75),
      nudge_x = 0.05) +
    geom_errorbar(aes(ymin=Richness-SD, 
                      ymax=Richness+SD), 
                  width= 0.05,
                  linetype = 2,
                  linewidth = 1,
                  alpha = 0.5) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", size = rel(0.25))) +
    ylab(label = "Richness") +
    xlab(label = "Sites/Samples") +
    scale_x_continuous(breaks = seq(1,2,1)) + 
    scale_colour_manual(values = cores) + 
    guides(col = "none",
           alpha = "none",
           lable = "none")
  
  
  
  # https://r-charts.com/ggplot2/grid/#google_vignette
  
  
  collector_plot_per_sample
  
  ggsave(file = paste0(figs_path,"/",
                       prjct_rad,"-",
                       "Coletor-per_sample",Sys.Date(),
                       # "-num_taxa_per_class-id80",
                       ".pdf",
                       collapse = ""),
         plot = collector_plot_per_sample,
         device = "pdf",
         width = 16,
         height = 12,
         units = "cm",
         dpi = 300)
  
  
  
  
  



