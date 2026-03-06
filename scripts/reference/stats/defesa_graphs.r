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

# Graficos Stats

classe <- c("Actinopteri", "Aves", "Mammalia")
reads_totais <- c(3553061, 351, 171429) 
reads_relativas <- c(0.9538, 0.00009, 0.04602)
ASVs <- c(563, 9, 131)
OTUs <- c(359, 9, 97)
Ordens <- c(10, 4, 7)
Familias <- c(22, 5, 10)
Generos <- c(40, 6, 11)
Especies <- c(44, 6, 11)

process <- as_tibble(data.frame(classe,
                                reads_totais,
                                reads_relativas,
                                ASVs,
                                OTUs,
                                Ordens,
                                Familias,
                                Generos,
                                Especies))
process %>% 
ggplot(aes(x = classe,
           y = reads_totais,
           # fill = `Family (BLASTn)`
           )) +
  # geom_bar(position = "fill") +
  geom_bar(stat = "identity") +
  geom_text(aes( #adicionando % dentro das barras
    label = sprintf("%1f%%",`reads_relativas`)),
    # position = position_fill(vjust = 1),
    color = "black") +
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
  # scale_fill_manual(values = c("#256676",
  #                              "#41c9dc",
  #                              "#d60724",
  #                              "#a6e590",
  #                              "#69306e",
  #                              "#84ee15",
  #                              "#531ce8",
  #                              "#64903a",
  #                              "#c637bc",
  #                              "#34f199",
  #                              "#873c1a",
  #                              "#f4cacb",
  #                              "#0b5313",
  #                              "#68affc",
  #                              "#1642cd",
  #                              "#ffa8ff",
  #                              "#21a708",
  #                              "#9e73b8",
  #                              "#f3d426",
  #                              "#cc7b6f",
  #                              "#fea53b",
  #                              "#4d4815"),
                    breaks = "Family (BLASTn)" +
                      name = "Curated Order (BLASTn)") +
  labs(x = "Local",
       y = "Proporção",
       title = "Proporções de famílias identificadas",
       subtitle = "Todas amostras")

fish_fam_all                                