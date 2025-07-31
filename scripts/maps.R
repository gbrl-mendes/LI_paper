# Ingleses Lake map ----

library(ggmap)

sample_sites <- tibble("Site" = c("SP01",
                                  "SP02",
                                  "SP03",
                                  "SP04"
                                  ),
                       "Lat" = c(-20.163846, 
                                 -20.174130,
                                 -20.179039, 
                                 -20.178366
                                 ),
                       "Long" = c(-43.955173, 
                                  -43.950338, 
                                  -43.959647, 
                                  -43.943100
                                  ))

map_limits <- tibble("Points" = c("P1", "P2", "P3", "P4"),
                       "Lat" = c(-20.16097, -20.15964, -20.19662, -20.19613),
                       "Long" = c( -43.98486, -43.93147, -43.98657, -43.92576))

# Calcule as coordenadas mínimas e máximas
lat_min <- min(map_limits$Lat)
lat_max <- max(map_limits$Lat)
long_min <- min(map_limits$Long)
long_max <- max(map_limits$Long)

# Registro Stadia Maps
register_stadiamaps(key = "4de1f9bb-f927-4729-8da2-b8c98b43485e")

# Defina a região do mapa com uma margem
bbox <- c(
  left = long_min - 0.00025,   # Longitude mínima
  bottom = lat_min - 0.00025, # Latitude mínima
  right = long_max + 0.00025,  # Longitude máxima
  top = lat_max + 0.00025     # Latitude máxima
)

bbox <- c(
  left = long_min,   # Longitude mínima
  bottom = lat_min, # Latitude mínima
  right = long_max,  # Longitude máxima
  top = lat_max     # Latitude máxima
)

# Obtenha o mapa de fundo (Stadiamap)
mapa <- get_stadiamap(bbox, zoom = 17, maptype = "outdoors")
ggmap(mapa)

# Network graph ----

library(tidygraph)
library(ggraph)

# edge table
edge_tbl <- resume_venn %>% 
  select(`Curated ID`, Site, Level) %>% 
  distinct() %>% 
  print()

# graph
graph <- edge_tbl[,1:2] %>% 
  as_tbl_graph(directed = FALSE)
print(graph)

# plotting graph
graph %>% ggraph(layout = "auto") +
  geom_edge_link(color = "gray", alpha = 0.7) +
  geom_node_point(aes(color = ifelse(name %in% unique(edge_tbl$`Curated ID`), "Species", "Sample site")), size = 5) +
  geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 3) +
  theme_void() +
  labs(color = "Type") +
  coord_fixed()

# Plotting map + network graph ----
edges <- graph %>%
  activate(edges) %>%
  as_tibble()

# 2. Combinar as arestas com as coordenadas dos nós de origem e destino
edges_with_coords <- edges %>%
  left_join(layout_g %>% as_tibble(), by = c("from" = ".ggraph.orig_index")) %>%  # Coordenadas de origem
  left_join(layout_g %>% as_tibble(), by = c("to" = ".ggraph.orig_index"))        # Coordenadas de destino
