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
  library("readxl")}

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

pre_raw_results_tbl <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/curated/curated-Complete_analysis_results-2024-01-10.xlsx") %>% tibble()
curated_ids_tbl <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/curated/curated_lagoa_ingleses-ASVs_x_amostras-2024-01-09.xlsx") %>% tibble()
blast_tax <- read.csv("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/tax_blast.csv", sep = ",", check.names = FALSE)

## Adicionando o RRA correto a pre_raw_results ----
{
  pre_raw_results_tbl <- pre_raw_results_tbl %>%
    mutate("Curated Relative abundance to all samples" = 0,
           "Curated Relative abundance on sample" = 0,
           "Curated Sample total abundance" = 0)
  
  abd_total <- sum(pre_raw_results_tbl$Abundance)
  
  pre_raw_results_tbl <- pre_raw_results_tbl %>%
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

## Adicionando a pre_raw_results_tbl os nomes que foram curados manualmente ----

## 1o selecionando apenas as colunas que importam

# Curated IDs
curated_ids_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

curated_ids_tbl <-
  curated_ids_tbl %>%
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

View(curated_ids_tbl)

# Verificar se existem IDs diferentes para a mesma ASV
curated_ids_tbl$`ASV header`[which(curated_ids_tbl$`ASV header` %>% duplicated())]

# Complete analysis
pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

pre_raw_results_tbl <-
  pre_raw_results_tbl %>%
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

View(pre_raw_results_tbl)

## 2o renomeando segundo os nomes curados e reordenando

pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

raw_results_tbl <- pre_raw_results_tbl %>%
  left_join(curated_ids_tbl,
            by = c("ASV header")) %>%
  select(c("ASV header", #oriundas da curated_ids_tbl
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

View(raw_results_tbl)

# Verificar se existem IDs diferentes para a mesma ASV
raw_results_tbl$`ASV header`[which(curated_ids_tbl$`ASV header` %>% duplicated())]

# Verificando o total de reads por amostra
raw_reads <- raw_results_tbl %>% group_by(Sample) %>%
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
dif_raw_pre <- pre_raw_results_tbl %>%
  anti_join(raw_results_tbl)

View(dif_raw_pre) #vazio e' bom!

# Apos ver que tudo que esta na tabela eh o que saiu do pipeline, podemos retirar 
# tudo o que estiver fora do intervalo do amplicon, identificacoes NA, brancos e
# possiveis contaminacoes

filt_results_tbl <-
  raw_results_tbl %>% 
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

View(filt_results_tbl)

# ASVs que foram retiradas
out_results_tbl <-
  setdiff(raw_results_tbl,filt_results_tbl)

# verificando o total de reads por amostra apos filtragem

filt_raw_reads <- filt_results_tbl %>% group_by(Sample) %>% 
  summarise(total_abd_filt = sum(Abundance))

View(filt_raw_reads)

stats_filt_reads <- raw_reads %>% 
  left_join(filt_raw_reads, 
            by = c("Sample")) %>% 
  mutate("Proportion left" = total_abd_filt / total_abd * 100)

View(stats_filt_reads)

# corrigindo o RRA

# Verificando o RRA
filt_results_tbl %>% 
  # grouped_by_ID_BLASTid %>%
  group_by(Sample) %>% 
  summarize(total_RRA = sum(`Relative abundance on sample`)) # Veja como que algumas amostras perderam muito!

# >     filt_results_tbl %>% 
#   +     # grouped_by_ID_BLASTid %>%
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

filt_results_tbl <- filt_results_tbl %>%
  mutate("Curated Relative abundance to all samples" = 0,
         "Curated Relative abundance on sample" = 0,
         "Curated Sample total abundance" = 0)

abd_total <- sum(filt_results_tbl$Abundance)

filt_results_tbl <- filt_results_tbl %>%
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
filt_results_tbl %>% 
  # grouped_by_ID_BLASTid %>%
  group_by(Sample) %>% 
  summarize(total_RRA = sum(`Relative abundance on sample`)) # Compare com o resultado anterior!

## 3o adicionando metadados que faltaram 

filt_results_tbl$Sample %>% unique() %>% sort() %>% paste0(collapse = '",\n"') %>% cat()
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
}

# As tibble
metadata_tbl <- tibble(run, 
                       Sample,
                       filter,
                       new_name, 
                       expedition, 
                       year)

# Mergindo os metadados com raw_results_tbl

filt_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

filt_results_tbl$Sample %>% unique()

curated_full_tbl <- left_join(filt_results_tbl, metadata_tbl,
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

View(curated_full_tbl)

diff_cur_raw <- dplyr::anti_join(#curated_ids_tbl,
  filt_results_tbl,
  curated_full_tbl,
  # ver se ocorreu sem problemas o left_join
  by = "ASV header")

View(diff_cur_raw) #vazio e' bom!

# Mergindo os metadados com out_results_tbl

out_results_full <- left_join(out_results_tbl, metadata_tbl,
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
           "Sample total abundance",
           "Relative abundance on sample",
           "Relative abundance to all samples",
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
  NW_pre_raw_results_tbl <- read_excel("/home/gabriel/projetos/peixes-eDNA/analises/dez_23/runs_2_4_5_EM156/results/lagoa_ingleses-Complete_analysis_results-2024-02-23.xlsx") %>% tibble()
  
  ## Adicionando o RRA correto a NW_pre_raw_results 
  {
    NW_pre_raw_results_tbl <- NW_pre_raw_results_tbl %>%
      mutate("Curated Relative abundance to all samples" = 0,
             "Curated Relative abundance on sample" = 0,
             "Curated Sample total abundance" = 0)
    
    abd_total <- sum(NW_pre_raw_results_tbl$Abundance)
    
    NW_pre_raw_results_tbl <- NW_pre_raw_results_tbl %>%
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
  NW_pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  
  NW_pre_raw_results_tbl <-
    NW_pre_raw_results_tbl %>%
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
  
  View(NW_pre_raw_results_tbl)
  
  ## Adicionando a info de contaminacao
  
  NW_contam <- pre_raw_results_tbl %>% 
    select(`Obs. Curadoria`,
           `Possible contamination`,
           `ASV header`)
  
  View(NW_contam)
  
  NW_pre_raw_results_tbl <- NW_pre_raw_results_tbl %>% 
    left_join(NW_contam,
              by = "ASV header")
  
  View(NW_pre_raw_results_tbl)
  
  # Verificar se os resultados sao os mesmos com pre_raw_results_tbl
  
  dif_NW_pre <- NW_pre_raw_results_tbl %>%
    anti_join(pre_raw_results_tbl)
  
  View(dif_NW_pre)
}

## Obter a nova tabela curated_IDs
NW_curated_ids_tbl <- read_excel("/home/gabriel/projetos/peixes-eDNA/analises/dez_23/runs_2_4_5_EM156/results/lagoa_ingleses-ASVs_x_amostras-2024-02-23.xlsx") %>% tibble()

## Selecionando apenas as colunas que importam para comparar

# Curated IDs

TT_curated_ids_tbl <-
  curated_ids_tbl %>%
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

# Verificar se os resultados sao os mesmos com pre_raw_results_tbl

dif_NW_pre <- NW_curated_ids_tbl %>%
  anti_join(curated_ids_tbl)    

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

curated_full_tbl %>% colnames()
curated_full_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

# Agrupamento da ASVs que possuem os mesmos atributos abaixo

grouped_by_ID_BLASTid <- curated_full_tbl %>%
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
  mutate(Nivel = case_when(`year` == 2020 ~ "Cheio", 
                           `year` == 2022 ~ "Cheio", 
                           TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o nivel da lagoa
  ungroup()

View(grouped_by_ID_BLASTid)

grouped_by_ID_BLASTid[which(grouped_by_ID_BLASTid %>% duplicated())]

# verificando o total de reads por amostra

grouped_reads <- grouped_by_ID_BLASTid %>% group_by(Sample) %>% 
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

blast_tax <- blast_tax[2:86] %>% 
  rename(`ASV (Sequence)` = "ASV")

blast_tax_less <- blast_tax %>%
  select(c("ASV (Sequence)",
           "Order (BLASTn)",
           "Family (BLASTn)",
           "Genus (BLASTn)"
  )) 
grouped_by_ID_BLASTid %>% colnames()

grouped_by_ID_BLASTid <- grouped_by_ID_BLASTid %>% 
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
           "Nivel" 
  )) %>% 
  unique()

View(grouped_by_ID_BLASTid)

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
    mutate(Nivel = case_when(`year` == 2020 ~ "Cheio", 
                             `year` == 2022 ~ "Cheio", 
                             TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o nivel da lagoa
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
             "Nivel" 
    )) %>% 
    unique()
  
  View(out_grouped_by_ID_BLASTid)
  
  write.csv(out_grouped_by_ID_BLASTid, paste0(tbl_path, "/out_grouped_2024.csv"))
}

## Alteracoes para o artigo 28 de maio ----

# Alteracoes baseadas nas sugestoes do Daniel, apos ver a tabela suplementar 2 da dissertacao.
# Como essa tabela é a versao wider da tabela dt_all_resume, tive que fazer as alteracoes abaixo
# na tabela grouped_by_ID_BLASTid para que eu possa aplicar o filtro sugerido pelo Daniel, que 
# envolve os numeros de reads apos serem agrupados como esta abaixo:

grouped_filt <- grouped_by_ID_BLASTid %>%
  group_by(`Curated ID`,
           `Final ID (BLASTn)`,
           `new_name`, 
           `Nivel`,
           `filter`,
           `year`
  ) %>%
  reframe("RRA" = RRA,
          "Sample" =  Sample,
          "ASVs" = length(unique(`ASV header`)),
          "OTUs" = length(unique(OTU)),
          "Nivel" = Nivel,
          "Ano" = year,
          "Filtro" = filter,
          "Reads" = sum(Abundance),
          "Classe" = `Class (BLASTn)`,
          "Order" = `Curated Order (BLASTn)`,
          "Family" = `Family (BLASTn)`,
          "Genus" = `Curated genus`,
          "expedition" = expedition
  ) %>%
  unique() %>%
  filter(Reads >= 100)

View(grouped_filt)

# ---- Tables ----
## Tabela suplementar 1 ----

# Tabela usada para avaliar os dados

# Tabela longer
dt_all_resume <- grouped_by_ID_BLASTid %>% 
  group_by(Nivel, new_name, filter, year) %>%
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
           `Nivel`,
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
    "Nivel" = Nivel,
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

View(dt_all_resume)

# options(scipen = 99,
#         digits = 3)

# Teste de RRA para ver se foi calculado corretamente
dt_all_resume %>%
  filter(Nivel == "Vazio") %>%
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
  unite(new_name, Nivel, Ano, Filtro, col= "ponto_nivel_ano_filtro") %>% 
  pivot_wider(id_cols = c("Classe","Curated ID", "Final ID (BLASTn)", "Reads totais", "Abundancia total"),
              names_from = ponto_nivel_ano_filtro,
              values_from =  c("RRA","ASVs","OTUs","Reads"),
              names_glue = "{ponto_nivel_ano_filtro}_{.value}") %>%
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

# Tabela longer
dt_all_filt <- grouped_by_ID_BLASTid %>%
  filter(`Class (BLASTn)` %in% "Actinopteri") %>% # Filtro de apenas peixes
  select(-c("Final ID (BLASTn)")) %>% # Retirando a identificacao original
  group_by(Nivel, new_name, filter, year) %>%
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
           `Nivel`,
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
    "Nivel" = Nivel,
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
  filter(Reads >= 100) %>% # Filtro de abundancia >= 100
  unique() 

View(dt_all_filt)

# options(scipen = 99,
#         digits = 3)

# Teste de RRA para ver se foi calculado corretamente
dt_all_filt %>%
  filter(Nivel == "Vazio") %>%
  filter(new_name == "Fundação") %>%
  filter(Filtro == "MCE") %>%
  filter(Ano == "2021") %>% 
  pull(RRA) %>%
  sum()  

# Tabela wider com os dados filtrados
wider_dt_filt <- dt_all_filt %>% 
  select(-c("RRA_formatado")) %>%
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  unite(new_name, Nivel, Ano, Filtro, col= "ponto_nivel_ano_filtro") %>% 
  pivot_wider(id_cols = c("Classe",
                          "Curated ID", 
                          # "Final ID (BLASTn)",
                          # "Reads totais", 
                          "Abundancia total"),
              names_from = ponto_nivel_ano_filtro,
              values_from =  c("RRA","ASVs","OTUs","Reads"),
              names_glue = "{ponto_nivel_ano_filtro}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate("Classe",
           "Curated ID",
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

## Tabela definitiva 1 ---- 

# desconsiderar o codigo abaixo. Ainda nao sei como ela deve ser
# ela ficou redundante com a construcao da Tabela Suplementar 2

# Apenas e' a Tabela suplementar 1 apenas com peixes

# Tabela longer
dt_filt_resume <- grouped_by_ID_BLASTid %>%
  filter(`Class (BLASTn)` %in% "Actinopteri") %>%
  group_by(Nivel, new_name, filter, year) %>%
  mutate("Abd total" = sum(Abundance)) %>%
  ungroup() %>%
  group_by(`Curated ID`) %>% 
  mutate("Abundancia total" = sum(Abundance)) %>% 
  ungroup %>% 
  group_by(`Final ID (BLASTn)`) %>% 
  mutate("Reads totais" = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(`Curated ID`,
           `new_name`, 
           `Nivel`,
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
    "Nivel" = Nivel,
    "Ano" = year,
    "Filtro" = filter,
    # "Ponto" = new_name,
    "Reads" = sum(Abundance),
    "Classe" = `Class (BLASTn)`,
    "Order" = `Curated Order (BLASTn)`,
    "Family" = `Family (BLASTn)`,
    "Genus" = `Curated genus`,
  ) %>% 
  mutate(
    "RRA_formatado" = format(RRA, scientific = TRUE)
  ) %>%
  unique()

View(dt_filt_resume)

# Tabela wider com os dados filtrados
wider_dt_filt <- dt_filt_resume %>% 
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  # unite(new_name, Nivel,col= "ponto_nivel") %>% 
  # pivot_wider(id_cols = c("Curated ID"),
  #             names_from = ponto_nivel,
  #             values_from =  c("RRA","ASVs","OTUs","Abundance"),
  #             names_glue = "{ponto_nivel}_{.value}") %>%
  # select(sort(colnames(.))) %>% 
  # relocate("Curated ID") 
  unite(new_name, Nivel, Ano, Filtro, col= "ponto_nivel_ano_filtro") %>% 
  pivot_wider(id_cols = c("Classe","Curated ID"
                          # , "Reads totais"
                          , "Abundancia total"
  ),
  names_from = ponto_nivel_ano_filtro,
  values_from =  c("RRA","ASVs","OTUs","Reads"),
  names_glue = "{ponto_nivel_ano_filtro}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate("Classe","Curated ID",
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
           # ,
           # "Ponte_Cheio_2022_MCE_Reads",
           # "Ponte_Cheio_2022_MCE_ASVs",
           # "Ponte_Cheio_2022_MCE_OTUs",
           # "Ponte_Cheio_2022_MCE_RRA",
           # "Fundação_Cheio_2022_MCE_Reads",
           # "Fundação_Cheio_2022_MCE_ASVs",
           # "Fundação_Cheio_2022_MCE_OTUs",
           # "Fundação_Cheio_2022_MCE_RRA"
  ) 

View(wider_dt_filt)

wider_dt_filt$`Curated ID` %>% unique() %>% sort()

fish_ID_tbl$`Curated ID` %>% unique() %>% sort()

write.csv(wider_dt_filt, paste0(tbl_path, "/wider_dt_filt_2024.csv"))

## Raw statistics ----

# Raw reads

# Run 2
# ~/projetos/peixes-eDNA/raw/2020_reads$ zgrep -c "^@M" *.gz

# Runs 4, 5
# ~/projetos/peixes-eDNA/raw/combined$ grep -c '^@' *.fastq.gz	

# Run EM156
#~/projetos/peixes-eDNA/raw/ecomol$ zgrep -c "^@VH" *.gz	

# Post-pipeline reads

# verificando o total de reads por amostra apos filtragem

filt_raw_reads <- filt_results_tbl %>% group_by(Sample) %>% 
  summarise(total_abd_filt = sum(Abundance))

View(filt_raw_reads)
sum(filt_raw_reads$total_abd_filt)


# Analise das amostras ----

## Esta tabela sumariza, para cada classe, em cada ano, a quantidade de reads, 
## o numero de ASVs, de OTUs, de ids a nivel de spp e ids em outros niveis taxonomicos

# Grouped by Class
analis_dt_class <-
  # grouped_by_ID_BLASTid %>%
  # filt_results_tbl %>%
  dt_all_resume %>%
  # dt_filt %>%
  group_by(
    Classe
    # `Curated ID`,
    # new_name,
    # Nivel
  ) %>%
  reframe(
    # "Total reads" = sum(Abundance), # grouped_by_ID_BLASTid && # filt_results_tbl
    "Total reads" = sum(Reads), # dt_all_resume && # dt_filt
    "Classe" = Classe,
    "ASVs" = sum(ASVs),
    "OTUs" = sum(OTUs),
    "Ordens" = length(unique(Order)),
    "Famílias" = length(unique(Family)),
    "Gêneros" = length(unique(Genus)),
    "Espécies" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])), # IDs a nivel de spp
    "Nspecies" = length(unique(`Curated ID`)) - Espécies # IDs em outros niveis taxonomicos
  ) %>% 
  unique()

# View(analis_pre_eco)
View(analis_dt_class)

# Verificando ASVs que nao foram identificadas a nivel de spp (genus + epitetus)
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

### Taxa detected Tabela Suplementar 2 ----

order_by_tax <- dt_all_resume %>%
  # dt_filt_resume %>% 
  filter(Classe %in% "Actinopteri") %>% 
  group_by(
    # `Final ID (BLASTn)`
    `Curated ID`
  ) %>%
  reframe(
    "Ordem" = Order,
    "Família" = Family,
    "Gênero" = Genus,
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

# write.csv(order_by_tax, "~/projetos/lagoa_ingleses/results/tabelas/order_by_tax.csv")

### Grouped by Order ----

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

### Grouped by Sample ----

  grouped_by_Sample <- fish_ID_tbl %>% 
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

### Grouped by Family ----

  grouped_by_Family <- fish_ID_tbl %>% 
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

### Especies encontradas em apenas uma amostra ----

counts <- fish_ID_tbl %>% 
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

stats_final_ID <- fish_ID_tbl %>% 
  mutate("Genus" = str_split_fixed(string = .$`Final ID (BLASTn)`, 
                                   pattern = " ",
                                   n = 2)[,1]) %>% 
  select(`Final ID (BLASTn)`, Genus)

### Proporcao de IDs do LGC12sDB e do NT/NCBI ----

  # Usando o objeto hits_DB criado no script 012024_post_pipe.r
  
  filt_hits_DB <- fish_ID_tbl %>% 
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

  ### Curva do coletor
  
  ou, pras curvas do coletor, pra tudo ou por grupo de amostras
  
  
  #Collector curve ----
  
  
  {r eval=FALSE, echo=TRUE}
  
  all_IDs_NMDS_df %>% colnames() %>% duplicated()
  all_IDs_NMDS_df$ID_Apareiodon_sp2 %>% jspecarize()
  
  
  
  ######## Function to transform values into 1 and zeros
  jaccarize <- function(x) {
    # Using ifelse to check each element of the vector x
    ifelse(x == 0, 0, 1)
  }
  ###################
  
  
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
  
  
  collector_plot <- spec2_tidy %>% ggplot(aes(x = Sample,y = Value,group=Sample))+
    # geom_point() +
    geom_boxplot(col ="#005602",
                 fill ="#005602",notch = TRUE,
                 alpha = 0.33,width = 0.5,outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.25,linetype=2) +
    geom_jitter(size=0.01,
                width = 0.25,
                # height = 0,
                col="#005602") +
    geom_line(data = spec2_tidy_rch,
              inherit.aes = F,
              linewidth = 2.5,
              col ="#005602",
              alpha = 0.50,
              aes(x = Sample, y = Richness)) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", size = rel(0.25))) +
    ylab(label = "Richness") +
    xlab(label = "Sites/Samples") +
    scale_x_continuous(breaks = seq(1,12,1))
  
  
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
  data <- data.frame("Sites"=spec1$sites, "Richness1"=spec1$richness,"SD1"=spec1$sd,
                     "Richness2"=spec1$richness, "SD2"=spec1$sd,
  )
  
  
  
  
  
  ggplot() +
    geom_point(data=data, aes(x=Sites, y=Richness)) +
    geom_line(data=data, aes(x=Sites, y=Richness)) +
    geom_ribbon(data=data ,aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)
  
  
  
  specpool(x, pool, smallsample = TRUE)
  estimateR(x, ...)
  specpool2vect(X, index = c("jack1","jack2", "chao", "boot","Species"))
  poolaccum(x, permutations = 100, minsize = 3)
  estaccumR(x, permutations = 100, parallel = getOption("mc.cores"))
  ## S3 method for class 'poolaccum'
  summary(object, display, alpha = 0.05, ...)
  ## S3 method for class 'poolaccum'
  plot(x, alpha = 0.05, type = c("l","g"), ...)
  
  
  
  
  
  
  
  
  
  #----
  
  
  #Collector curve per sample  ----
  
  
  {r eval=FALSE, echo=TRUE}
  
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
  
  
  
  
  #----
  



