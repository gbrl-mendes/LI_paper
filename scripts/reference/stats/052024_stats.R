"Estatísticas"
"Mendes, GA; Hilário, OH"
"12/2023"

# Pacotes ----
{
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(base)
  library(writexl)
}

# Obtencao dos dados ----

# Usando os objetos criados pelo codigo 012024_post_pipe.r:
# 1. grouped_by_ID_BLASTid


# Tabela suplementar 1 ----

## Tabela usada para avaliar os dados

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

write.csv(wider_dt_all_resume, "~/projetos/lagoa_ingleses/results/tabelas/wider_dt_all_resume_2024.csv")


# Tabela suplementar 3 ----

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

write.csv(wider_dt_all_resume, "~/projetos/lagoa_ingleses/results/tabelas/wider_dt_all_resume_2024.csv")

# Tabela definitiva 1 ----

## Apenas e' a Tabela suplementar 1 apenas com peixes

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

write.csv(wider_dt_filt, "~/projetos/lagoa_ingleses/results/tabelas/wider_dt_filt_2024.csv")


# Raw statistics ----

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
    "Espécies" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])),
    "Nspecies" = length(unique(`Curated ID`)) - Espécies
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
      "Nspecies" = setdiff(unique(`Curated ID`), spp)
    )
  
  View(nspp_level) 
}

# Grouped by Order
{
  order_level <- 
    dt_all_resume %>%
    filter(Classe %in% "Actinopteri") %>%
    group_by(
      # Order
    ) %>%
    reframe(
      # "Curated ID" = `Curated ID`,
      "Total reads" = sum(Reads), # dt_all_resume && # dt_filt
      # "Classe" = Classe,
      "ASVs" = sum(ASVs),
      "OTUs" = sum(OTUs),
      "Famílias" = length(unique(Family)),
      "Gêneros" = length(unique(Genus)),
      "Ids" = length(unique(`Curated ID`)),
      "Espécies" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])),
      "Nspecies" = length(unique(`Curated ID`)) - Espécies
    ) %>% 
    unique() 
  
  View(order_level)
}

# Taxa detected Tabela Suplementar 2

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

write.csv(order_by_tax, "~/projetos/lagoa_ingleses/results/tabelas/order_by_tax.csv")

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
    }

# Grouped by Family
{
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
}

# Especies encontradas em apenas uma amostra

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

# Proporcao de IDs do LGC12sDB e do NT/NCBI

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

group_DB %>% 
  write_xlsx(path = "/home/gabriel/projetos/lagoa_ingleses/results/tabelas/group_DB.xlsx", 
                    col_names = TRUE,
                    format_headers = TRUE)

## Upload curated_group_DB

curated_group_DB  <- read_excel("/home/gabriel/projetos/lagoa_ingleses/results/tabelas/curated_group_DB.xlsx")

## Grafico de barras
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

##################

# Grouped by Curated ID
options(scipen = 999)

grouped_sampled_filt_tbl %>% 
  mutate("Total reads" = sum(Reads)) %>% 
  group_by(`Curated ID`) %>% 
  reframe("Total species reads" = sum(Reads),
            "Relative abundance for species" = `Total species reads` / `Total reads` * 100,
          "Order" = Order) %>% 
  relocate(Order, `Curated ID`) %>% 
  unique() %>% 
  View()

# Grouped by Order

grouped_sampled_filt_tbl %>% 
  mutate("Total reads" = sum(Reads)) %>% 
  group_by(Order) %>% 
  reframe("Total order reads" = sum(Reads),
          "Relative abundance for order" = `Total order reads` / `Total reads` * 100,
          "Total ASVs" = sum(ASVs),
          "Total Species" = n_distinct(`Curated ID`)) %>% 
  unique() %>% 
  View()
