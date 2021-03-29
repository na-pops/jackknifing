####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: jackknifing
# 2-jackknife-distance.R
# Created March 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(detect)
library(plyr)
library(doParallel)
library(foreach)

####### Read Data #################################

load(file = "data/dist_count_matrix.rda")
load(file = "data/landcover_covariates.rda")
load(file = "data/dist_design.rda")
na_sp_list <- read.csv("../utilities/IBP-Alpha-Codes20.csv")

####### Data Wrangling ############################

# Drop most of the landcover covariates for now
covariates_reduced <- landcover_covariates[, c("Sample_ID", "roaddist", "ForestOnly_5x5")]
covariates_reduced$ForestOnly_5x5 <- covariates_reduced$ForestOnly_5x5 / 25

# Drop invalid road distances, and create roadside variable
covariates_reduced <- covariates_reduced[which(covariates_reduced$roaddist > 0), ]
covariates_reduced$roadside <- ifelse(covariates_reduced$roaddist < 30, 1, 0)

# Drop method I
dist_count_matrix <- dist_count_matrix[-which(dist_count_matrix$Distance_Method == "I"), ]
dist_design <- dist_design[-which(dist_design$Method == "I"), ]

max_bands <- ncol(dist_design) - 2
count_names <- c("Sample_ID", "Species", "Distance_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands),
                 "Project")
names(dist_count_matrix) <- count_names

design_names <- c("Distance_Method", "Max_Distance",
                  paste0(rep("Interval", times = max_bands), 1:max_bands))
names(dist_design) <- design_names

# Join data
count_design <- plyr::join(dist_count_matrix, dist_design,
                           by = "Distance_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Distance_Method", "Project", count_names[5:length(count_names)-1])]
design <- count_design[,c("Sample_ID", "Species", "Distance_Method", "Project", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[5:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Add in landcover covariates to counts matrix
covars <- plyr::join(counts, covariates_reduced, by = "Sample_ID", type = "left", match = "all")

# Remove rows with NA covariates
to_remove <- which(is.na(covars$roadside))
counts <- counts[-c(to_remove), ]
design <- design[-c(to_remove), ]
covars <- covars[-c(to_remove), ]

to_remove <- which(is.na(covars$ForestOnly_5x5))
counts <- counts[-c(to_remove), ]
design <- design[-c(to_remove), ]
covars <- covars[-c(to_remove), ]

########### Full Model ############################

# Build matrices by species
species_all <- sort(as.character(unique(counts$Species)))
landbirds <- na_sp_list[which(na_sp_list$LANDBIRD == "TRUE"), "SPEC"]
species <- species_all[which(species_all %in% landbirds)]

# Create species directories
sapply(species, function (x) dir.create(paste0("data/distance/", x)))

for (s in species)
{
  if (nrow(counts[counts$Species == s, ]) >= 50)
  {
    assign(paste0("Y_",s), as.matrix(counts[counts$Species==s, col_names]))
    assign(paste0("D_",s), as.matrix(design[design$Species==s, col_names]))
    assign(paste0("C_",s), subset(covars,
                                  Species==s,
                                  select = c(19:20)))
  }else
  {
    species <- species[!(species %in% s)]
  }
}

input_list <- vector(mode="list", length=length(species))
names(input_list) <- species
for (s in species)
{
  input_list[[s]] <- list(Y=get(paste0("Y_",s)), D=get(paste0("D_",s)), C=get(paste0("C_",s)))
}

cluster <- makeCluster(12, type = "PSOCK")
registerDoParallel(cluster)

foreach(sp = names(input_list), .packages = 'detect') %dopar%
  {
    x <- input_list[[sp]]
    m1 = cmulti(x$Y | x$D ~ 1, type="dis")
    m2 = cmulti(x$Y | x$D ~ x$C$roadside, type="dis")
    m3 = cmulti(x$Y | x$D ~ x$C$ForestOnly_5x5, type="dis")
    m4 = cmulti(x$Y | x$D ~ x$C$roadside + x$C$ForestOnly_5x5, type="dis")
    m5 = cmulti(x$Y | x$D ~ x$C$roadside * x$C$ForestOnly_5x5, type="dis")
    distance_list <- list(m1, m2, m3, m4, m5)
    save(distance_list, file = paste0("data/distance/", sp, "/0full.rda"))    
  }
stopCluster(cluster)

########### Jackknifing #############################
project_list <- unique(counts$Project)

cluster <- makeCluster(24, type = "PSOCK")
registerDoParallel(cluster)

foreach(p = project_list, .packages = 'detect') %dopar%
{
  temp <- counts[which(counts$Project == p), ]
  species_red <- sort(as.character(unique(temp$Species)))
  counts_red <- counts[-which(counts$Project == p), ]
  design_red <- design[-which(design$Project == p), ]
  covars_red <- covars[-which(covars$Project == p), ]
  
  for (s in species_red)
  {
    assign(paste0("Y_",s), as.matrix(counts_red[counts_red$Species==s, col_names]))
    assign(paste0("D_",s), as.matrix(design_red[design_red$Species==s, col_names]))
    assign(paste0("C_",s), subset(covars_red,
                                  Species==s,
                                  select = c(19:20)))
  }
  
  input_list <- vector(mode="list", length=length(species_red))
  names(input_list) <- species_red
  for (s in species_red)
  {
    input_list[[s]] <- list(Y=get(paste0("Y_",s)), D=get(paste0("D_",s)), C=get(paste0("C_",s)))
  }
  
  for (sp in species_red)
  {
    x <- input_list[[sp]]
    m1 = cmulti(x$Y | x$D ~ 1, type="dis")
    m2 = cmulti(x$Y | x$D ~ x$C$roadside, type="dis")
    m3 = cmulti(x$Y | x$D ~ x$C$ForestOnly_5x5, type="dis")
    m4 = cmulti(x$Y | x$D ~ x$C$roadside + x$C$ForestOnly_5x5, type="dis")
    m5 = cmulti(x$Y | x$D ~ x$C$roadside * x$C$ForestOnly_5x5, type="dis")
    distance_list <- list(m1, m2, m3, m4, m5)
    save(distance_list, file = paste0("data/distance/", sp, "/", p, ".rda"))  
  }
}

stopCluster(cluster)