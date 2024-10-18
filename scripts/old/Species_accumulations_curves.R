## Este codigo e´ originado do script LI_paper.Rmd.
## Retirei este trecho do script original pois ficou
## redundante após a criacao do plot nested com todas
## as curvas em uma so´ figura. Porem este codigo ainda
## pode ser util, portanto salvei em um arquivo separado.
## G. Mendes 01/10/24

# Todas amostras MCE Individuals

## ASVs
{
  # Data re-formatting using ASVs
  wider_grouped_sampled_filt_ASV <- grouped_sampled_filt_tbl %>%
    select(c(Sample,
             Site,
             Level,
             `Curated ID`,
             ASVs,
             Expedition,
             Filter,
             Level)) %>%
    group_by(Sample, `Curated ID`, ASVs) %>%
    unique %>%
    pivot_wider(names_from = `Curated ID`,
                values_from = ASVs,
                names_sort = TRUE) %>% as.data.frame()
  
  rownames(wider_grouped_sampled_filt_ASV) <- wider_grouped_sampled_filt_ASV$Sample
  
  wider_grouped_sampled_filt_ASV <- wider_grouped_sampled_filt_ASV %>% 
    select(-c(Sample))
  
  wider_grouped_sampled_filt_ASV <- wider_grouped_sampled_filt_ASV %>% 
    mutate(across(everything(), ~replace_na(.x, 0)))
  
  wider_grouped_sampled_filt_ASV <- wider_grouped_sampled_filt_ASV %>% 
    filter(Filter %in% "MCE")
  
  grouped_accum_1 <- specaccum(comm = wider_grouped_sampled_filt_ASV[,5:ncol(wider_grouped_sampled_filt_ASV)]
                               ,
                               method = "rarefaction"
                               # method = "random"
  )
  plot(grouped_accum_1$individuals, grouped_accum_1$richness)
  
  # Getting into a dataframe
  grouped_accum_1_df <- tibble("Sites" = c(0, grouped_accum_1$sites),
                               "Individuals" = c(0, grouped_accum_1$individuals),
                               "Richness" = c(0, grouped_accum_1$richness),
                               "sd" = c(0, grouped_accum_1$sd))
  
  collector_plot_MCE_ASV <-
    grouped_accum_1_df %>% 
    ggplot(aes(x = Individuals,
               y = Richness)) +
    geom_line(linewidth = 2.5,
              col ="#4CAF50",
              alpha = 0.50) +
    geom_point(size = 2, color = "#4CAF50") +
    # geom_errorbar(aes(ymin = Richness - sd, ymax = Richness + sd), width = 0.2) +
    geom_ribbon(aes(ymin = Richness - sd, ymax = Richness + sd), fill = "#005602", alpha = 0.2) +
    # theme_minimal() +
    labs(x = "ASVs", y = "Riqueza de espécies") +
    ylab(label = "Richness") +
    xlab(label = "ASVs") +
    ggtitle("Collector's curve MCE samples (2020 - 2022)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  collector_plot_MCE_ASV
  
  
  ggplot(spec_accum_df, aes(x = Sites, y = Richness)) +
    geom_line(color = "blue") +              # Curva de acumulação
    geom_point(size = 2, color = "blue") +   # Pontos da curva
    geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), width = 0.2) +  # Barra de erro
    
    ggtitle("Curva de Acumulação de Espécies (Collector's Curve) com Desvio Padrão") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file = paste0(figs_path,"/", "collector_plot_MCE_all", Sys.Date(), ".pdf", collapse = ""),
         plot = collector_plot_MCE,
         device = "pdf",
         width = 16,
         height = 12,
         units = "cm",
         dpi = 300)
}

## Reads
{
  # Data re-formatting using Reads
  wider_grouped_sampled_filt_Reads <- grouped_sampled_filt_tbl %>%
    select(c(Sample,
             Site,
             Level,
             `Curated ID`,
             Reads,
             Expedition,
             Filter,
             Level)) %>%
    group_by(Sample, `Curated ID`, Reads) %>%
    unique %>%
    pivot_wider(names_from = `Curated ID`,
                values_from = Reads,
                names_sort = TRUE) %>% as.data.frame()
  
  rownames(wider_grouped_sampled_filt_Reads) <- wider_grouped_sampled_filt_Reads$Sample
  
  wider_grouped_sampled_filt_Reads <- wider_grouped_sampled_filt_Reads %>% 
    select(-c(Sample))
  
  wider_grouped_sampled_filt_Reads <- wider_grouped_sampled_filt_Reads %>% 
    mutate(across(everything(), ~replace_na(.x, 0)))
  
  wider_grouped_sampled_filt_Reads <- wider_grouped_sampled_filt_Reads %>% 
    filter(Filter %in% "MCE")
  
  grouped_accum_1 <- specaccum(comm = wider_grouped_sampled_filt_Reads[,5:ncol(wider_grouped_sampled_filt_Reads)]
                               ,
                               method = "rarefaction"
                               # method = "random"
  )
  plot(grouped_accum_1$individuals, grouped_accum_1$richness)
  
  # Getting into a dataframe
  grouped_accum_1_df <- tibble("Sites" = c(0, grouped_accum_1$sites),
                               "Individuals" = c(0, grouped_accum_1$individuals),
                               "Richness" = c(0, grouped_accum_1$richness),
                               "sd" = c(0, grouped_accum_1$sd))
  
  collector_plot_MCE_Reads <-
    grouped_accum_1_df %>% 
    ggplot(aes(x = Individuals,
               y = Richness)) +
    geom_line(linewidth = 2.5,
              col ="#4CAF50",
              alpha = 0.50) +
    geom_point(size = 2, color = "#4CAF50") +
    # geom_errorbar(aes(ymin = Richness - sd, ymax = Richness + sd), width = 0.2) +
    geom_ribbon(aes(ymin = Richness - sd, ymax = Richness + sd), fill = "#005602", alpha = 0.2) +
    # theme_minimal() +
    labs(x = "Reads", y = "Riqueza de espécies") +
    ylab(label = "Richness") +
    xlab(label = "Reads") +
    ggtitle("Collector's curve MCE samples (2020 - 2022)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  collector_plot_MCE_Reads
  
  
  ggplot(spec_accum_df, aes(x = Sites, y = Richness)) +
    geom_line(color = "blue") +              # Curva de acumulação
    geom_point(size = 2, color = "blue") +   # Pontos da curva
    geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), width = 0.2) +  # Barra de erro
    
    ggtitle("Curva de Acumulação de Espécies (Collector's Curve) com Desvio Padrão") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file = paste0(figs_path,"/", "collector_plot_MCE_all", Sys.Date(), ".pdf", collapse = ""),
         plot = collector_plot_MCE,
         device = "pdf",
         width = 16,
         height = 12,
         units = "cm",
         dpi = 300)
}