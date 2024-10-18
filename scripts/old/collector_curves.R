library(iNEXT)
library(ecodados)
library(ggplot2)
library(vegan)
library(nlme)
library(dplyr)
library(piecewiseSEM)

# Interpolation method ----

## Rarefaction based on number of individuals 

## Definindo os dados STX vs MCE
wider2_grouped_STX_MCE <- grouped_sampled_filt_tbl %>%
  select(`Curated ID`, 
         Filter, 
         Reads) %>%
  group_by(`Curated ID`, 
           Filter) %>%
  summarize(Reads = sum(Reads)) %>%
  pivot_wider(names_from = Filter, 
              values_from = Reads, 
              names_sort = TRUE) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  column_to_rownames("Curated ID")

## Analysis with Interpolation method
int_curve_STX_MCE <- iNEXT(wider_grouped_STX_MCE, q = 0,
                   datatype = "abundance", endpoint = 8000)

ggiNEXT(int_curve_STX_MCE, type = 1) +
  geom_vline(xintercept = 166, lty = 2) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  # scale_colour_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  # scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  labs(x = "Número de indivíduos", y = " Riqueza de espécies")

## Rarefaction based on number of samples 

# Function to transform values into 1 and zeros
jaccarize <- function(x) {
  # Using ifelse to check each element of the vector x
  ifelse(x == 0, 0, 1)
  }

wider2_grouped_STX_MCE_jac <- wider2_grouped_STX_MCE %>%
  mutate(across(everything(), jaccarize))


data_inext <- as.incfreq(t(wider2_grouped_STX_MCE_jac)) 

resultados_incidencia <- iNEXT(data_inext, q = 0, datatype = "incidence_freq", 
                               endpoint = 28)

## Gráfico
ggiNEXT(resultados_incidencia, type = 1) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  # scale_colour_manual(values = "darkorange") +
  # scale_fill_manual(values = "darkorange") +
  labs(x = "Número de amostras", y = " Riqueza de espécies")
