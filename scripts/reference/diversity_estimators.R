library(iNEXT)
library(ggplot2)
library(vegan)
library(nlme)
library(dplyr)
library(ggpubr)

# Diversity estimators ----

## Warning! 
## This diversity estimators only uses data of species presence/absence! 

## All samples ----
{
  ### Data obtainment ----
  
  # All samples
  wider_grouped_sampled_filt[5:ncol(wider_grouped_sampled_filt)]
  
  ### Chao2 method ----
  
  # Estimation
  est_chao2 <- poolaccum(wider_grouped_sampled_filt[5:ncol(wider_grouped_sampled_filt)], permutations = 100)
  summary(est_chao2, display = "chao")
  
  est_chao2 <- poolaccum(wider_grouped_sampled_filt[5:ncol(wider_grouped_sampled_filt)], 
                         minsize = 0, permutations = 100)
  summary(est_chao2, display = "chao")
  
  # Data reshaping for ggplot2
  resultados_chao2 <- summary(est_chao2, display = c("S", "chao"))
  res_chao2 <- cbind(resultados_chao2$chao[, 1:4], resultados_chao2$S[, 2:4])
  res_chao2 <- as.data.frame(res_chao2)
  colnames(res_chao2) <- c("Samples", "Chao2", "C_inferior", "C_superior",
                           "Richness", "R_inferior", "R_superior")
  
  ### Jackknife1 method ----
  
  # Estimation
  est_jack1 <- poolaccum(wider_grouped_sampled_filt[5:ncol(wider_grouped_sampled_filt)], permutations = 100)
  summary(est_jack1, display = "jack1")
  
  # Data reshaping for ggplot2
  resultados_jack1 <- summary(est_jack1, display = c("S", "jack1"))
  res_jack1 <- cbind(resultados_jack1$jack1[, 1:4], resultados_jack1$S[, 2:4])
  res_jack1 <- as.data.frame(res_jack1)
  colnames(res_jack1) <- c("Samples", "JACK1", "JACK1_inferior", "JACK1_superior", 
                           "Richness", "R_inferior", "R_superior")

  
  ### Jackknife 2 method ----
  
  # Estimation
  est_jack2 <- poolaccum(wider_grouped_sampled_filt[5:ncol(wider_grouped_sampled_filt)], permutations = 100)
  summary(est_jack2, display = "jack2")
  
  # Data reshaping for ggplot2
  resultados_jack2 <- summary(est_jack2, display = c("S", "jack2"))
  res_jack2 <- cbind(resultados_jack2$jack2[, 1:4], resultados_jack2$S[, 2:4])
  res_jack2 <- as.data.frame(res_jack2)
  colnames(res_jack2) <- c("Samples", "jack2", "jack2_inferior", "jack2_superior", 
                           "Richness", "R_inferior", "R_superior")
  
  ### Merging the estimators in one plot ----
  
  res_long <- bind_rows(
    res_jack2 %>% select(Samples, Richness, Inferior = R_inferior, Superior = R_superior) %>% mutate(Method = "Obs"), 
    res_jack2 %>% select(Samples, Richness = jack2, Inferior = jack2_inferior, Superior = jack2_superior) %>% mutate(Method = "Jack2"),
    res_jack1 %>% select(Samples, Richness = JACK1, Inferior = JACK1_inferior, Superior = JACK1_superior) %>% mutate(Method = "Jack1"),
    res_chao2 %>% select(Samples, Richness = Chao2, Inferior = C_inferior, Superior = C_superior) %>% mutate(Method = "Chao2"))
  
  # Individually plotting diversity estimators
  
  jack1_plot <- 
    res_long %>% 
    filter(Method %in% c("Obs", "Jack1")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "Chao2"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "Chao2" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "Chao2")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  jack2_plot <-
    res_long %>% 
    filter(Method %in% c("Obs", "Jack2")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack2" = "#FF7F0E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "Chao2"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "Chao2" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "Chao2")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  chao2_plot <-
    res_long %>% 
    filter(Method %in% c("Obs", "Chao2")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Chao2" = "#2CA02C")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "Chao2"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "Chao2" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "Chao2")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  # Merging the plots
  estimators_plot_all <-
    ggarrange(jack1_plot, jack2_plot, chao2_plot,
              ncol = 3, nrow = 1,
              common.legend = TRUE, legend = "right")
  
  }

## All MCE samples ----
{
  ### Data obtainment ----
  
  # All MCE samples
  wider_grouped_sampled_filt_mce
  
  ### chao2_mce method ----
  
  # Estimation
  
  est_chao2_mce <- poolaccum(wider_grouped_sampled_filt_mce, 
                         minsize = 0, permutations = 100)
  summary(est_chao2_mce, display = "chao")
  
  # Data reshaping for ggplot2
  resultados_chao2_mce <- summary(est_chao2_mce, display = c("S", "chao"))
  res_chao2_mce <- cbind(resultados_chao2_mce$chao[, 1:4], resultados_chao2_mce$S[, 2:4])
  res_chao2_mce <- as.data.frame(res_chao2_mce)
  colnames(res_chao2_mce) <- c("Samples", "chao2_mce", "C_inferior", "C_superior",
                           "Richness", "R_inferior", "R_superior")
  
  ### Jackknife1 method ----
  
  # Estimation
  est_jack1_mce <- poolaccum(wider_grouped_sampled_filt_mce, minsize = 0, permutations = 100)
  summary(est_jack1_mce, display = "jack1")
  
  # Data reshaping for ggplot2
  resultados_jack1_mce <- summary(est_jack1_mce, display = c("S", "jack1"))
  res_jack1_mce <- cbind(resultados_jack1_mce$jack1[, 1:4], resultados_jack1_mce$S[, 2:4])
  res_jack1_mce <- as.data.frame(res_jack1_mce)
  colnames(res_jack1_mce) <- c("Samples", "JACK1", "JACK1_inferior", "JACK1_superior", 
                           "Richness", "R_inferior", "R_superior")
  
  
  ### Jackknife 2 method ----
  
  # Estimation
  est_jack2_mce <- poolaccum(wider_grouped_sampled_filt_mce, minsize = 3, permutations = 100)
  summary(est_jack2_mce, display = "jack2")
  
  # Data reshaping for ggplot2
  resultados_jack2_mce <- summary(est_jack2_mce, display = c("S", "jack2"))
  res_jack2_mce <- cbind(resultados_jack2_mce$jack2[, 1:4], resultados_jack2_mce$S[, 2:4])
  res_jack2_mce <- as.data.frame(res_jack2_mce)
  colnames(res_jack2_mce) <- c("Samples", "jack2", "jack2_inferior", "jack2_superior", 
                           "Richness", "R_inferior", "R_superior")
  
  ### Merging the estimators in one plot ----
  
  res_long_mce <- bind_rows(
    res_jack2_mce %>% select(Samples, Richness, Inferior = R_inferior, Superior = R_superior) %>% mutate(Method = "Obs"), 
    res_jack2_mce %>% select(Samples, Richness = jack2, Inferior = jack2_inferior, Superior = jack2_superior) %>% mutate(Method = "Jack2"),
    res_jack1_mce %>% select(Samples, Richness = JACK1, Inferior = JACK1_inferior, Superior = JACK1_superior) %>% mutate(Method = "Jack1"),
    res_chao2_mce %>% select(Samples, Richness = chao2_mce, Inferior = C_inferior, Superior = C_superior) %>% mutate(Method = "chao2_mce"))
  
  # Individually plotting diversity estimators
  
  jack1_plot_mce <- 
    res_long_mce %>% 
    filter(Method %in% c("Obs", "Jack1")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_mce"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_mce" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_mce")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  jack2_plot_mce <-
    res_long_mce %>% 
    filter(Method %in% c("Obs", "Jack2")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack2" = "#FF7F0E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_mce"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_mce" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_mce")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  chao2_plot_mce <-
    res_long_mce %>% 
    filter(Method %in% c("Obs", "chao2_mce")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "chao2_mce" = "#2CA02C")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_mce"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_mce" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_mce")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  # Merging the plots
  estimators_plot_mce <-
    ggarrange(jack1_plot_mce, jack2_plot_mce, chao2_plot_mce,
              ncol = 3, nrow = 1,
              common.legend = TRUE, legend = "right")
  
}

## MCE 2021 samples ----
{
  ### Data obtainment ----
  
  # All MCE samples
  wider_grouped_sampled_filt_mce21
  
  ### chao2 method ----
  
  # Estimation
  
  est_chao2_mce21 <- poolaccum(wider_grouped_sampled_filt_mce21, minsize = 0, permutations = 100)
  summary(est_chao2_mce21, display = "chao")
  
  # Data reshaping for ggplot2
  resultados_chao2_mce21 <- summary(est_chao2_mce21, display = c("S", "chao"))
  res_chao2_mce21 <- cbind(resultados_chao2_mce21$chao[, 1:4], resultados_chao2_mce21$S[, 2:4])
  res_chao2_mce21 <- as.data.frame(res_chao2_mce21)
  colnames(res_chao2_mce21) <- c("Samples", "chao2_mce21", "C_inferior", "C_superior",
                               "Richness", "R_inferior", "R_superior")
  
  ### Jackknife1 method ----
  
  # Estimation
  est_jack1_mce21 <- poolaccum(wider_grouped_sampled_filt_mce21, minsize = 0, permutations = 100)
  summary(est_jack1_mce21, display = "jack1")
  
  # Data reshaping for ggplot2
  resultados_jack1_mce21 <- summary(est_jack1_mce21, display = c("S", "jack1"))
  res_jack1_mce21 <- cbind(resultados_jack1_mce21$jack1[, 1:4], resultados_jack1_mce21$S[, 2:4])
  res_jack1_mce21 <- as.data.frame(res_jack1_mce21)
  colnames(res_jack1_mce21) <- c("Samples", "JACK1", "JACK1_inferior", "JACK1_superior", 
                               "Richness", "R_inferior", "R_superior")
  
  
  ### Jackknife 2 method ----
  
  # Estimation
  est_jack2_mce21 <- poolaccum(wider_grouped_sampled_filt_mce21, minsize = 0, permutations = 100)
  summary(est_jack2_mce21, display = "jack2")
  
  # Data reshaping for ggplot2
  resultados_jack2_mce21 <- summary(est_jack2_mce21, display = c("S", "jack2"))
  res_jack2_mce21 <- cbind(resultados_jack2_mce21$jack2[, 1:4], resultados_jack2_mce21$S[, 2:4])
  res_jack2_mce21 <- as.data.frame(res_jack2_mce21)
  colnames(res_jack2_mce21) <- c("Samples", "jack2", "jack2_inferior", "jack2_superior", 
                               "Richness", "R_inferior", "R_superior")
  
  ### Merging the estimators in one plot ----
  
  res_long_mce21 <- bind_rows(
    res_jack2_mce21 %>% select(Samples, Richness, Inferior = R_inferior, Superior = R_superior) %>% mutate(Method = "Obs"), 
    res_jack2_mce21 %>% select(Samples, Richness = jack2, Inferior = jack2_inferior, Superior = jack2_superior) %>% mutate(Method = "Jack2"),
    res_jack1_mce21 %>% select(Samples, Richness = JACK1, Inferior = JACK1_inferior, Superior = JACK1_superior) %>% mutate(Method = "Jack1"),
    res_chao2_mce21 %>% select(Samples, Richness = chao2_mce21, Inferior = C_inferior, Superior = C_superior) %>% mutate(Method = "chao2_mce21"))
  
  # Individually plotting diversity estimators
  
  jack1_plot_mce21 <- 
    res_long_mce21 %>% 
    filter(Method %in% c("Obs", "Jack1")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_mce21"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_mce21" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_mce21")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  jack2_plot_mce21 <-
    res_long_mce21 %>% 
    filter(Method %in% c("Obs", "Jack2")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack2" = "#FF7F0E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_mce21"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_mce21" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_mce21")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  chao2_plot_mce21 <-
    res_long_mce21 %>% 
    filter(Method %in% c("Obs", "chao2_mce21")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "chao2_mce21" = "#2CA02C")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_mce21"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_mce21" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_mce21")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  # Merging the plots
  estimators_plot_mce21 <-
    ggarrange(jack1_plot_mce21, jack2_plot_mce21, chao2_plot_mce21,
              ncol = 3, nrow = 1,
              common.legend = TRUE, legend = "right")
  
  estimators_plot_mce21
}

## Sterivex samples ----
{
  ### Data obtainment ----
  
  # All Sterivex samples
  wider_grouped_sampled_filt_stx
  
  ### chao2 method ----
  
  # Estimation
  
  est_chao2_stx <- poolaccum(wider_grouped_sampled_filt_stx, minsize = 0, permutations = 100)
  summary(est_chao2_stx, display = "chao")
  
  # Data reshaping for ggplot2
  resultados_chao2_stx <- summary(est_chao2_stx, display = c("S", "chao"))
  res_chao2_stx <- cbind(resultados_chao2_stx$chao[, 1:4], resultados_chao2_stx$S[, 2:4])
  res_chao2_stx <- as.data.frame(res_chao2_stx)
  colnames(res_chao2_stx) <- c("Samples", "chao2_stx", "C_inferior", "C_superior",
                                 "Richness", "R_inferior", "R_superior")
  
  ### Jackknife1 method ----
  
  # Estimation
  est_jack1_stx <- poolaccum(wider_grouped_sampled_filt_stx, minsize = 0, permutations = 100)
  summary(est_jack1_stx, display = "jack1")
  
  # Data reshaping for ggplot2
  resultados_jack1_stx <- summary(est_jack1_stx, display = c("S", "jack1"))
  res_jack1_stx <- cbind(resultados_jack1_stx$jack1[, 1:4], resultados_jack1_stx$S[, 2:4])
  res_jack1_stx <- as.data.frame(res_jack1_stx)
  colnames(res_jack1_stx) <- c("Samples", "JACK1", "JACK1_inferior", "JACK1_superior", 
                                 "Richness", "R_inferior", "R_superior")
  
  ### Jackknife 2 method ----
  
  # Estimation
  est_jack2_stx <- poolaccum(wider_grouped_sampled_filt_stx, minsize = 3, permutations = 100)
  summary(est_jack2_stx, display = "jack2")
  
  # Data reshaping for ggplot2
  resultados_jack2_stx <- summary(est_jack2_stx, display = c("S", "jack2"))
  res_jack2_stx <- cbind(resultados_jack2_stx$jack2[, 1:4], resultados_jack2_stx$S[, 2:4])
  res_jack2_stx <- as.data.frame(res_jack2_stx)
  colnames(res_jack2_stx) <- c("Samples", "jack2", "jack2_inferior", "jack2_superior", 
                                 "Richness", "R_inferior", "R_superior")
  
  ### Merging the estimators in one plot ----
  
  res_long_stx <- bind_rows(
    res_jack1_stx %>% select(Samples, Richness, Inferior = R_inferior, Superior = R_superior) %>% mutate(Method = "Obs"), 
    res_jack2_stx %>% select(Samples, Richness = jack2, Inferior = jack2_inferior, Superior = jack2_superior) %>% mutate(Method = "Jack2"),
    res_jack1_stx %>% select(Samples, Richness = JACK1, Inferior = JACK1_inferior, Superior = JACK1_superior) %>% mutate(Method = "Jack1"),
    res_chao2_stx %>% select(Samples, Richness = chao2_stx, Inferior = C_inferior, Superior = C_superior) %>% mutate(Method = "chao2_stx"))
  
  # Individually plotting diversity estimators
  
  jack1_plot_stx <- 
    res_long_stx %>% 
    filter(Method %in% c("Obs", "Jack1")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_stx"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_stx" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_stx")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  jack2_plot_stx <-
    res_long_stx %>% 
    filter(Method %in% c("Obs", "Jack2")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "Jack2" = "#FF7F0E")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_stx"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_stx" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_stx")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  chao2_plot_stx <-
    res_long_stx %>% 
    filter(Method %in% c("Obs", "chao2_stx")) %>% 
    ggplot(aes(y = Richness,
               x = Samples,
               color = Method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3,
               show.legend = TRUE) +
    # scale_color_manual(values = c("Obs" = "#3366CC", "chao2_stx" = "#2CA02C")) +
    scale_color_manual(limits = c("Obs", "Jack1", "Jack2", "chao2_stx"),
                       values = c("Obs" = "#3366CC", "Jack1" = "#B82E2E", "Jack2" = "#FF7F0E", "chao2_stx" = "#2CA02C"),
                       labels = c("Observed", "Jackknife1", "Jackknife2", "chao2_stx")) +
    geom_line() +
    geom_errorbar(aes(ymin = Inferior, ymax = Superior), position = position_dodge(width = 0.5), width = 0.2)
  
  # Merging the plots
  estimators_plot_stx <-
    ggarrange(jack1_plot_stx, jack2_plot_stx, chao2_plot_stx,
              ncol = 3, nrow = 1,
              common.legend = TRUE, legend = "right")
  
  estimators_plot_stx
}

# Extrapolation  ----

## All samples ----
{
  ### Data obtainment ----
  wider_data_all_jac <- as.data.frame(lapply(wider_grouped_sampled_filt[5:ncol(wider_grouped_sampled_filt)], jaccarize)) 
  
  ### Extrapolation ----
  dados_inext <- as.incfreq(t(wider_data_all_jac)) 
  incid_all <- iNEXT(dados_inext, q = 0, datatype = "incidence_freq", 
                                        endpoint = 30)
  
  ### Plot ----
  extra_plot_mce21 <-
    ggiNEXT(incid_all, type = 1) +
    scale_linetype_discrete(labels = c("Interpolated", "Extrapolated")) +
    scale_colour_manual(values = "darkorange") +
    scale_fill_manual(values = "darkorange") +
    labs(x = "Sample number", y = " Species richness")
  extra_plot_mce21
}

## MCE21 samples ----
{
  ### Data obtainment ----
  wider_data_mce21_jac <-  as.data.frame(lapply(wider_grouped_sampled_filt_mce21, jaccarize))
  
  ### Extrapolation ----
  dados_inext_mce21 <- as.incfreq(t(wider_data_mce21_jac)) 
  incid_mce21 <- iNEXT(dados_inext_mce21, q = 0, datatype = "incidence_freq",
                                              endpoint = 10)
  
  ### Plot ----
  extra_plot_all <-
    ggiNEXT(incid_mce21, type = 1) +
    scale_linetype_discrete(labels = c("Interpolated", "Extrapolated")) +
    scale_colour_manual(values = "darkorange") +
    scale_fill_manual(values = "darkorange") +
    labs(x = "Sample number", y = " Species richness")
  extra_plot_all
}


## Sterivex samples ----
{
  ### Data obtainment ----
  wider_data_stx_jac <-  as.data.frame(lapply(wider_grouped_sampled_filt_stx, jaccarize))
  
  ### Extrapolation ----
  dados_inext_stx <- as.incfreq(t(wider_data_stx_jac)) 
  incid_stx <- iNEXT(dados_inext_stx, q = 0, datatype = "incidence_freq",
                       endpoint = 36)
  
  ### Plot ----
  extra_plot_stx <-
    ggiNEXT(incid_stx, type = 1) +
    scale_linetype_discrete(labels = c("Interpolated", "Extrapolated")) +
    scale_colour_manual(values = "darkorange") +
    scale_fill_manual(values = "darkorange") +
    labs(x = "Sample number", y = " Species richness")
  extra_plot_stx
}


## Combined MCE21 and STX samples ----
{
  ### Data obtainment ----
  wider_data_stx_mce_jac <-  as.data.frame(lapply(wider_grouped_sampled_filt_stx_mce, jaccarize))
  
  ### Extrapolation ----
  dados_inext_stx_mce <- as.incfreq(t(wider_data_stx_mce_jac)) 
  incid_stx_mce <- iNEXT(dados_inext_stx_mce, q = 0, datatype = "incidence_freq",
                     endpoint = 30)
  
  ### Plot ----
  extra_plot_stx_mce <-
    ggiNEXT(incid_stx_mce, type = 1) +
    scale_linetype_discrete(labels = c("Interpolated", "Extrapolated")) +
    scale_colour_manual(values = "darkorange") +
    scale_fill_manual(values = "darkorange") +
    labs(x = "Sample number", y = " Species richness")
  extra_plot_stx_mce
}


combined_plot <- geom_line(data = extra_plot_all$data,
                           aes(x = x, y = y, color = extra_plot_all$data$Assemblage)) +
  geom_line(data = extra_plot_stx_mce$data,
            aes(x = x, y = y, color = extra_plot_stx_mce$data$Assemblage))


