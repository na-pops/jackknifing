####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: jackknifing
# 3-extract-coefficients.R
# Created March 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(detect)
library(doParallel)
library(foreach)

####### Set Constants #############################

n_dis_models <- 5

####### Distance Model Selection ##################

species <- list.files(path = "data/distance")

sp_coef_list <- vector(mode = "list", length = length(species))
names(sp_coef_list) <- species

foreach (sp = species, .packages = 'detect') %dopar%
{
  dir <- paste0("data/distance/", sp, "/")
  projects <- substr(list.files(path = paste0(dir)),
                     start = 1,
                     stop = nchar(list.files(path = paste0(dir))) - 4)
  
  dist_coef <- data.frame(Project = rep(projects, each = n_dis_models),
                          n = NA, 
                          model = rep(seq(1,n_dis_models), times = length(projects)), 
                          intercept = NA,
                          road = NA,
                          forest = NA,
                          roadforest = NA)
  
  dis_vcv_list <- vector(mode = "list", length = n_dis_models)
  proj_list <- vector(mode = "list", length = length(projects))
  names(proj_list) <- projects
  for (m in 1:n_dis_models)
  {
    dis_vcv_list[[m]] <- proj_list
  }
  
  for (s in projects)
  {
    load(file = paste0(dir, s, ".rda"))
    dist_coef[which(dist_coef$Project == s), "n"] <- nrow(distance_list[[1]]$Y)
    
    for (m in 1:n_dis_models)
    {
      coef <- coef(distance_list[[m]])
      if (m == 1) {
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"intercept"] = coef[1]
      }
      if (m == 2) {
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"intercept"] = coef[1]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"road"] = coef[2]
      }
      if (m == 3) {
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"intercept"] = coef[1]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"forest"] = coef[2]
      }
      if (m == 4) {
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"intercept"] = coef[1]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"road"] = coef[2]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"forest"] = coef[3]
      }
      if (m == 5) {
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"intercept"] = coef[1]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"road"] = coef[2]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"forest"] = coef[3]
        dist_coef[which(dist_coef$Project == s & dist_coef$model == m),"roadforest"] = coef[4]
      }    
      
      dis_vcv_list[[m]][[s]] <- distance_list[[m]]$vcov
    }
  }  
  
  dir.create(path = paste0(dir, "coef"))
  write.table(dist_coef, 
              file = paste0(dir, "coef/distance.csv"),
              sep = ",",
              row.names = FALSE)
  
  save(dis_vcv_list, file = paste0(dir, "coef/dis_vcv_list.rda"))
  
  save(projects, file = paste0(dir, "coef/project_list.rda"))
}




