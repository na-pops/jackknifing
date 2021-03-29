####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: jackknifing
# 1-combine-subset-data.R
# Created March 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(reshape2)

source("../utilities/get-data.R")

####### Set Constants #############################

sp <- c("PAWR", "WIWR", "WTSP", "BARS", "NRWS", "HOLA",
        "ACFL", "BAWW", "BTNW", "RBNU", "WOTH")

####### Read Data #################################

time <- read.csv("../covariates/survey/time_lookup.csv")
dist <- read.csv("../covariates/survey/distance_lookup.csv")

# Get project names
project_list <- read.table("../utilities/proj-list")
n_proj <- nrow(project_list)

# Create combined count data frame
project_counts <- vector('list', n_proj)
project_samples <- vector('list', n_proj)
for (i in 1:n_proj)
{
  p <- project_list$V1[i]
  # Get counts
  data_dir <- paste0("../project-",
                     p,
                     "/output/",
                     p,
                     "_counts.rda")
  project_counts[[i]] <- data.frame(get_data(data_dir))
  
  # Get samples
  data_dir <- paste0("../project-",
                     p,
                     "/output/",
                     p,
                     "_samples.rda")
  project_samples[[i]] <- data.frame(get_data(data_dir))
}
project_counts <- do.call(rbind, project_counts)
project_samples <- do.call(rbind, project_samples)

# Create combined landcover covariate df and temporal covariate df
landcover_covariates <- vector('list', n_proj)
temporal_covariates <- vector('list', n_proj)
for (i in 1:n_proj)
{
  p <- project_list$V1[i]
  
  # Landcover covariates
  data_dir <- paste0("../covariates/landcover/project-",
                     p,
                     ".csv")
  landcover_covariates[[i]] <- read.csv(data_dir)
  
  # Temporal covariates
  data_dir <- paste0("../covariates/temporal/project-",
                     p,
                     ".csv")
  temporal_covariates[[i]] <- read.csv(data_dir)
}
landcover_covariates <- do.call(rbind, landcover_covariates)
temporal_covariates <- do.call(rbind, temporal_covariates)

####### Wrangle Data ##############################

project_counts_reduced <- project_counts[which(project_counts$Species %in% sp), ]
project_samples_reduced <- project_samples[which(project_samples$Sample_ID %in% 
                                                   project_counts_reduced$Sample_ID), ]

# Create time count matrix
time_only <- project_counts_reduced[, c(1:3, 7,8 )]
time_only <- time_only[-which(time_only$Time_Method == "ZZ"), ]
time_only <- time_only[-which(is.na(time_only$Time_Level)), ]
time_count_matrix <- dcast(time_only,
                           Sample_ID + Species + Time_Method ~ as.numeric(Time_Level),
                           value.var = "Abundance",
                           fun.aggregate = sum)
time_count_matrix <- merge(x = time_count_matrix,
                           y = project_samples_reduced[, c("Sample_ID", "Project")],
                           by = "Sample_ID")


# Create distance count matrix
dist_only <- project_counts_reduced[, c(1:5)]
dist_only <- dist_only[-which(dist_only$Distance_Method == "ZZ"), ]
if (any(is.na(dist_only$Distance_Level)))
{
  dist_only <- dist_only[-which(is.na(dist_only$Distance_Level)), ]  
}

dist_count_matrix <- dcast(dist_only,
                           Sample_ID + Species + Distance_Method ~ as.numeric(Distance_Level),
                           value.var = "Abundance",
                           fun.aggregate = sum)
dist_count_matrix <- merge(x = dist_count_matrix,
                           y = project_samples_reduced[, c("Sample_ID", "Project")],
                           by = "Sample_ID")

# Create time removal design matrix
time <- time[-c(1), ]
time_design <- dcast(time, Method + Max_Duration ~ Level, value.var = "End_Duration")

# Create distance design matrix
dist <- dist[-c(1), ]
dist_design <- dcast(dist, Method + Max_Distance ~ Level, value.var = "End_Distance")

####### Output Data ###############################
save(project_counts_reduced, file = "data/counts.rda")
save(project_samples_reduced, file = "data/samples.rda")
save(time_count_matrix, file = "data/time_count_matrix.rda")
save(dist_count_matrix, file = "data/dist_count_matrix.rda")
save(landcover_covariates, file = "data/landcover_covariates.rda")
save(temporal_covariates, file = "data/temporal_covariates.rda")
save(time_design, file = "data/time_design.rda")
save(dist_design, file = "data/dist_design.rda")