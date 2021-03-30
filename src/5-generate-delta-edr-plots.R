####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: jackknifing
# 5-generate-delta-edr-plots.R
# Created March 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(ggplot2)
library(GGally)
library(ggpubr)
theme_set(theme_pubclean())

species <- list.files(path = "data/distance")
sapply(paste0("plots/mod", seq(1,5), "/"), dir.create)

for (sp in species)
{
  data_dir <- paste0("data/distance/", sp, "/")
  
  ####### Model 1 Plot ###############################
  load(paste0(data_dir, "tau/tau_1.rda"))
  tau_1 <- tau_1[!duplicated(tau_1$Project), ]
  n_row <- nrow(tau_1)
  
  p1 <- ggplot() +
    geom_point(data = tau_1[2:n_row, ],
               aes(x = Project, y = tau)) +
    geom_errorbar(data = tau_1[2:n_row, ],
                  aes(x = Project,
                      ymin = tau_2.5,
                      ymax = tau_97.5,
                      width = 0.2)) +
    geom_hline(data = tau_1[1, ],
               aes(yintercept = tau)) +
    geom_hline(data = tau_1[1, ],
               aes(yintercept = tau_2.5),
               linetype = "dashed") +
    geom_hline(data = tau_1[1, ],
               aes(yintercept = tau_97.5),
               linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("Species: ", sp, ", Model 1")) +
    NULL
  
  png(filename = paste0("plots/mod1/", sp, ".png"),
      width = 15,
      height = 10,
      units = "in",
      res = 300)
  print(p1)
  dev.off()
  
}