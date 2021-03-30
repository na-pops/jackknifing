####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: jackknifing
# 4-calculate-tau.R
# Created March 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(MASS)
library(foreach)
library(doParallel)

####### Read Data #################################

#dis <- rm_non_sp(read.csv("../results/coefficients/distance.csv"))
family <- read.csv("../utilities/NACC_list_species.csv")[, c("common_name",
                                                             "family")]
ibp <- read.csv("../utilities/IBP-Alpha-Codes20.csv")[, c("SPEC",
                                                          "COMMONNAME")]
#load("../results/var-covar/dis_vcv_list.rda")

####### Set Constants #############################

#n_cores <- max(dis$model)

####### Simulate tau ##############################

forest_coverage <- c(0,1)
roadside <- c(rep(1, length(forest_coverage)),
              rep(0, length(forest_coverage)))

species <- list.files(path = "data/distance")

cluster <- makeCluster(12, type = "PSOCK")
registerDoParallel(cluster)

foreach (sp = species, .packages = 'MASS') %dopar%
{
  dir <- paste0("data/distance/", sp, "/")
  dir.create(paste0(dir, "tau"))
  dis <- read.csv(paste0(dir, "coef/distance.csv"))
  load(paste0(dir, "coef/project_list.rda"))
  load(paste0(dir, "coef/dis_vcv_list.rda"))
  
  for (i in 1:5)
  {
    dis_reduced <- dis[which(dis$model == i), ]
    sim_data <- data.frame(Project = rep(projects,
                                         each = length(roadside)),
                           Intercept = rep(1,
                                           times = length(projects) *
                                             length(roadside)),
                           Roadside = rep(roadside, length(projects)),
                           Forest = rep(forest_coverage, 2 * length(projects)))
    sim_data$Interaction <- sim_data$Forest * sim_data$Roadside
    
    tau <- NULL
    tau_low <- NULL
    tau_high <- NULL
    
    for (p in projects)
    {
      design <- sim_data[which(sim_data$Project == p),
                         c("Intercept", "Roadside", "Forest", "Interaction")]
      
      coefficients <- as.numeric(dis_reduced[which(dis_reduced$Project == p),
                                             c("intercept", "road",
                                               "forest", "roadforest")])
      zeros_indices <- which(is.na(coefficients)) - 1
      if (length(zeros_indices) > 0)
      {
        coefficients <- coefficients[-which(is.na(coefficients))]    
      }
      
      vcv <- dis_vcv_list[[i]][[p]]
      
      # Simulate a bunch of possible coefficients
      # Try-catch surrounding for the odd non positive definite var-covar
      sim_coef <-tryCatch(
        {
          rbind(coefficients, MASS::mvrnorm(10^4, coefficients, vcv))
          #sim_coef <- rbind(coefficients, MASS::mvrnorm(10^4, coefficients, vcv))      
        },
        error = function(e)
        {
          return(NA)
        }
      )
      
      if (is.na(sim_coef))
      {
        # In the case of an error, just output the calculated tau and NA
        # for the upper and lower
        if (length(zeros_indices) > 0)
        {
          coef_zeros <- c(coefficients, rep(0, length(zeros_indices)))
          id <- c(seq_along(coefficients), zeros_indices + 0.5)
          coef_zeros <- coef_zeros[order(id)]    
        }else
        {
          coef_zeros <- coefficients
        }
        
        tau_pred <- exp(as.matrix(design) %*% (coef_zeros))
        tau <- c(tau, as.numeric(tau_pred[,1]))
        tau_low <- c(tau_low, rep(NA, nrow(design)))
        tau_high <- c(tau_high, rep(NA, nrow(design)))    
      }else
      {
        # Add columns of zeros back in to where NA coefficients were previously
        # See https://stackoverflow.com/a/1495204/5665609 for explanation
        if (length(zeros_indices) > 0)
        {
          coef_zeros <- cbind(sim_coef, matrix(0,
                                               ncol = length(zeros_indices),
                                               nrow = nrow(sim_coef)))
          id <- c(seq_along(sim_coef[1,]), zeros_indices + 0.5)
          coef_zeros <- coef_zeros[,order(id)]    
        }else
        {
          coef_zeros <- sim_coef
        }
        
        tau_pred <- exp(as.matrix(design) %*% t(coef_zeros))
        tau <- c(tau, as.numeric(tau_pred[,1]))
        
        # Calculate quantiles
        tau_pred <- tau_pred[,-1]
        tau_low <- c(tau_low,
                     as.numeric(apply(tau_pred,
                                      1,
                                      quantile,
                                      probs = c(0.025),
                                      na.rm = TRUE)))
        tau_high <- c(tau_high,
                      as.numeric(apply(tau_pred,
                                       1,
                                       quantile,
                                       probs = c(0.975),
                                       na.rm = TRUE)))       
      }
      
    }
    sim_data$tau <- tau
    sim_data$tau_2.5 <- tau_low
    sim_data$tau_97.5 <- tau_high
    
    # ####### Simulate q ################################  
    # radius_values <- seq(50, 400, by = 50)
    # radius <- rep(radius_values, times = nrow(sim_data))
    # 
    # sim_data <- sim_data[rep(seq_len(nrow(sim_data)),
    #                          each = length(radius_values)), ]
    # sim_data$Radius <- radius
    # 
    # sim_data$q <- ifelse(sim_data$Radius == "Inf",
    #                      1,
    #                      ((sim_data$tau ^ 2) / (sim_data$Radius ^ 2)) *
    #                        (1 - exp(-(sim_data$Radius ^ 2) /
    #                                   sim_data$tau ^ 2)))
    # 
    # sim_data$q_2.5 <- ifelse(sim_data$Radius == "Inf",
    #                          1,
    #                          ((sim_data$tau_2.5 ^ 2) / (sim_data$Radius ^ 2)) *
    #                            (1 - exp(-(sim_data$Radius ^ 2) /
    #                                       sim_data$tau_2.5 ^ 2)))
    # 
    # sim_data$q_97.5 <- ifelse(sim_data$Radius == "Inf",
    #                           1,
    #                           ((sim_data$tau_97.5 ^ 2) / (sim_data$Radius ^ 2)) *
    #                             (1 - exp(-(sim_data$Radius ^ 2) /
    #                                        sim_data$tau_97.5 ^ 2)))
    
    ####### Save Results ############################ 
    var_name <- paste0("tau_", i)
    assign(var_name, sim_data)
    #obj <-  eval(parse(text = paste0("tau_", i)))
    save(list = var_name,
         file = paste0(dir,
                       "tau/",
                       var_name,
                       ".rda"))
  }
}

stopCluster(cluster)