grouped2_filt_tbl <- filt_tbl %>% 
  # creating total abundance by period
  group_by(Sample) %>%
  mutate("Abd total" = sum(Abundance)) %>%
  ungroup() %>%
  # creating total abundance by Curated ID
  group_by(`Curated ID`) %>% 
  mutate("Abundancia total" = sum(Abundance)) %>% 
  ungroup() %>% 
  # calculating RRA by period
  group_by(`Curated ID`, Sample) %>%
  mutate("RRA no periodo" = (Abundance/`Abd total`)*100) %>% 
  reframe("Curated ID" = `Curated ID`,
          "Order" = `Order (BLASTn)`,
          "Family" = `Family (BLASTn)`,
          "Genus" = `Genus (BLASTn)`,
          "Site" = new_name,
          "Level" = level,
          "Year" = year,
          "Filter" = filter,
          "Reads" = sum(Abundance),
          "Total abundance" = `Abundancia total`, 
          "RRA" = sum(`RRA no periodo`),
          "ASVs" = length(unique(`ASV header`)),
          "OTUs" = length(unique(OTU))) %>% 
  mutate("RRA_formatado" = format(RRA, scientific = TRUE)) %>%
  unique()

View(grouped2_filt_tbl)

# Teste para ver se o RRA foi calculado corretamente

grouped2_filt_tbl %>%
  # group_by(Sample) %>%    
  group_by(Sample) %>%
  # group_by(`Curated ID`, new_name, level, filter, year) %>%
  # summarise(total_RRA = sum(`RRA no periodo`, na.rm = TRUE))
  summarise(total_RRA = sum(RRA, na.rm = TRUE)) %>%
  ungroup()

# Tabela wider com os dados filtrados
wider2_filt_tbl <- grouped2_filt_tbl %>% 
  select(-c("RRA_formatado")) %>%
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  unite(Site, Level, Year, Filter, expedicol= "site_level_year_filter") %>% 
  pivot_wider(id_cols = c("Curated ID", "Total abundance"),
              names_from = site_level_year_filter,
              values_from =  c("RRA","ASVs","OTUs","Reads"),
              names_glue = "{site_level_year_filter}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate("Curated ID",
           "Total abundance",
           #2020 MCE
           "Ponte_Full_2020_MCE_Reads", 
           "Ponte_Full_2020_MCE_ASVs", 
           "Ponte_Full_2020_MCE_OTUs",
           "Ponte_Full_2020_MCE_RRA",
           "Fundação_Full_2020_MCE_Reads",
           "Fundação_Full_2020_MCE_ASVs",
           "Fundação_Full_2020_MCE_OTUs",
           "Fundação_Full_2020_MCE_RRA",
           #2021 MCE
           "Prainha_Low_2021_MCE_Reads",
           "Prainha_Low_2021_MCE_ASVs",
           "Prainha_Low_2021_MCE_OTUs",
           "Prainha_Low_2021_MCE_RRA",
           "Barragem_Low_2021_MCE_Reads",
           "Barragem_Low_2021_MCE_ASVs",
           "Barragem_Low_2021_MCE_OTUs",
           "Barragem_Low_2021_MCE_RRA",
           "Ponte_Low_2021_MCE_Reads",
           "Ponte_Low_2021_MCE_ASVs",
           "Ponte_Low_2021_MCE_OTUs",
           "Ponte_Low_2021_MCE_RRA",
           "Fundação_Low_2021_MCE_Reads",
           "Fundação_Low_2021_MCE_ASVs",
           "Fundação_Low_2021_MCE_OTUs",
           "Fundação_Low_2021_MCE_RRA",
           #2021 Sterivex
           "Prainha_Low_2021_Sterivex_Reads",
           "Prainha_Low_2021_Sterivex_ASVs",
           "Prainha_Low_2021_Sterivex_OTUs",
           "Prainha_Low_2021_Sterivex_RRA",
           "Barragem_Low_2021_Sterivex_Reads",
           "Barragem_Low_2021_Sterivex_ASVs",
           "Barragem_Low_2021_Sterivex_OTUs",
           "Barragem_Low_2021_Sterivex_RRA",
           "Ponte_Low_2021_Sterivex_Reads",
           "Ponte_Low_2021_Sterivex_ASVs",
           "Ponte_Low_2021_Sterivex_OTUs",
           "Ponte_Low_2021_Sterivex_RRA",
           "Fundação_Low_2021_Sterivex_Reads",
           "Fundação_Low_2021_Sterivex_ASVs",
           "Fundação_Low_2021_Sterivex_OTUs",
           "Fundação_Low_2021_Sterivex_RRA",
           #2022 MCE
           "Prainha_Full_2022_MCE_Reads",
           "Prainha_Full_2022_MCE_ASVs",
           "Prainha_Full_2022_MCE_OTUs",
           "Prainha_Full_2022_MCE_RRA",
           "Barragem_Full_2022_MCE_Reads",
           "Barragem_Full_2022_MCE_ASVs",
           "Barragem_Full_2022_MCE_OTUs",
           "Barragem_Full_2022_MCE_RRA"
           # "Ponte_Full_2022_MCE_Reads", # Os dois prox pontos nao sobreviveram 
           # "Ponte_Full_2022_MCE_ASVs", # aos filtros por que nao tem ASVs de peixes 
           # "Ponte_Full_2022_MCE_OTUs",
           # "Ponte_Full_2022_MCE_RRA",
           # "Fundação_Full_2022_MCE_Reads",
           # "Fundação_Full_2022_MCE_ASVs",
           # "Fundação_Full_2022_MCE_OTUs",
           # "Fundação_Full_2022_MCE_RRA"
  ) 

wider_filt_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

View(wider_filt_tbl)