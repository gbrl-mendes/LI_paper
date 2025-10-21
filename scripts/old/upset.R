library(UpSetR)

upset_data <- resume_venn %>%
  select(`Curated ID`, Site, Level) %>% 
  view()

upset_matrix <- upset_data %>%
  mutate(Site_Level = paste(Site, Level, sep = "_"), Presence = 1) %>%
  select(`Curated ID`, Site_Level, Presence) %>%  # Apenas as colunas necessárias
  pivot_wider(names_from = Site_Level, values_from = Presence, values_fill = list(Presence = 0)) %>%
  as.data.frame()

upset(upset_matrix, sets = colnames(upset_matrix)[-1], order.by = "freq", keep.order = TRUE)
