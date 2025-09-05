
library(tidyverse)
library(betaC)
library(vegan)

all.files <- list.files(path = "gamma_dataset",
                        pattern = "*.txt",
                        full.names = TRUE,
                        recursive = TRUE)
filenames.non.integer <- list.files(path = "gamma_dataset/non-integer",
                                    pattern = "*.txt")

for (i in 1:length(all.files)) assign(
  x = basename(all.files[[i]]), 
  value = read.table(file = all.files[[i]],
                     header = TRUE,
                     sep = " ",
                     fill = TRUE,
                     row.names = NULL))
all.files <- basename(all.files)

for (i in 1:length(filenames.non.integer)) {
  print(filenames.non.integer[[i]])
  stu1 <- noquote(filenames.non.integer[[i]])
  stu2 <- get(stu1, 1)
  dat_abund.0 <- as.data.frame(stu2[ , -(1:14)])
  a <- ceiling(dat_abund.0)
  b <- cbind(stu2[ , (1:14)], a)
  assign(filenames.non.integer[i], b)
}

# Checking n of each study
n <- matrix(0, 67, 1)
for (x in 1:length(all.files)) {
  stu1 <- noquote(all.files[x])
  stu2 <- get(stu1, 1)
  n[x,] <- nrow(stu2)
}
rownames(n) <- all.files
n

#Beta_spie function ENS method
beta_SPIE <- function(x) {
  x = as.matrix(x)
  total = colSums(x)
  gamma_value = as.numeric(betaC::calc_PIE(total, ENS = T))
  alpha_value = mean(apply(x, 1, betaC::calc_PIE, ENS = T))
  beta = gamma_value / alpha_value
  return(beta)
}

# #########count spatial block number in each dataset########
# n.block <- list()
# for (m in 1:length(all.files)) {
#   print(all.files[[m]])
#   stu1 <- noquote(all.files[[m]])
#   Dataset <- get(stu1, 1)
#   Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
#     sf::st_set_crs(4326)
#   
#   block <- Dataset |>
#     group_by(Spatial_block) |> 
#     summarise(Count = n()) |> 
#     sf::st_drop_geometry()
#   
#   Dataset_id <- distinct(Dataset,Dataset_id)
#   n.block[[m]] <- cbind.data.frame(Dataset_id,block)
# }

##########calculate the minmun N for Gamma_Sn and C_target for Beta_C########

SC <- list()

for (m in 1:length(all.files)) {
  print(all.files[[m]])
  stu1 <- noquote(all.files[[m]])
  Dataset <- get(stu1, 1)

  Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
  stu4 <- cbind.data.frame(Block, Dataset)
  spl.block <- split(stu4, Block) #spliting in datasets according to block
  blocos <- summary(as.factor(Block))
  n.blocks <- length(blocos)
  # 
  # save.metrics <- matrix(ncol = 2)
  # colnames(save.metrics) <- c("C.coms", "N")
  # 
  NS <- matrix(ncol = 1)
  colnames(NS) <- c("NS")
  
  C.coms <- matrix(ncol = 1)
  colnames(C.coms) <- c("C.coms")
  
  # start the loop between blocks
  for (t in 1:n.blocks) {
    Dataset <- spl.block[[t]] 
    
    Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
      sf::st_set_crs(4326)
    
    lu_area <- dplyr::bind_cols(
      Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
      area = Dataset |>
        dplyr::group_by(Land_use) |>
        dplyr::summarise() |>
        sf::st_convex_hull() |>
        sf::st_area() |>
        as.numeric()
    )
    #plots number counting
    sites_number <- Dataset |> 
      group_by(Land_use) |> 
      summarise(length(Plot_id)) |> 
      sf::st_drop_geometry()
    
    # check if min area site also have the least sites
    sites_min <- sites_number[which.min(sites_number$`length(Plot_id)`), "Land_use"] |> dplyr::pull()
    threshold_sites <- min(sites_number$`length(Plot_id)`)
    
    threshold_lu <- lu_area[which.min(lu_area$area), "Land_use"] |> dplyr::pull() 
    threshold_area <- min(lu_area$area)
    # only do next calculation if min area land use same with also min sites land use
    if (threshold_sites > 2 & sites_min == threshold_lu) {
      
      #calculate the metrics in min_area_land_use
      threshold_data <- Dataset |>
        dplyr::filter(Land_use == threshold_lu) |> 
        sf::st_drop_geometry()
      meta_info_min <- threshold_data |>
        slice(1) |>
        select(Dataset_id, Spatial_block,
               Temporal_block, Plot_id, Replicate_id,
               Taxa, Continent, Land_use) |> 
        sf::st_drop_geometry()
      
      comm1 <- threshold_data |> 
        dplyr::select(-(1:13))
      soma1 <- (colSums(comm1))
      tsoma1 <- t(soma1)
      N1 <- rowSums(tsoma1)
      NS <- rbind(na.omit(NS),N1)
      C.comm1 <-C_target(comm1,factor = 2)
      C.coms <- rbind(na.omit(C.coms),C.comm1)
      
      #prepare for the next step
      n_sites_ref_lu <- Dataset |>
        dplyr::filter(Land_use == threshold_lu) |>
        nrow()
      
      result <- Dataset |>
        dplyr::select(Plot_id, Land_use) |> 
        sf::st_drop_geometry()
      
      #calculate the metric in the other land-use categories
      add_meta <- list()
      
      for (i in lu_area |>
           dplyr::filter(area > threshold_area) |>
           dplyr::select(Land_use) |>
           dplyr::pull()) {
        temp_data <- Dataset |>
          dplyr::filter(Land_use == i)
        meta_info <- temp_data |>
          slice(1) |>
          select(Dataset_id,Block,Spatial_block,
                 Temporal_block,Plot_id, Replicate_id,
                 Taxa,Continent,Land_use) |> 
          sf::st_drop_geometry()
        
        n_sites_i_lu <- nrow(temp_data)
        
        # All combinations of sites down to the number of sites the smallest lu has
        if (choose(length(temp_data$Plot_id), n_sites_ref_lu) < 200) {
          combs <- combn(temp_data$Plot_id,
                         n_sites_ref_lu,
                         simplify = FALSE)
        }
        
        if (choose(length(temp_data$Plot_id),  n_sites_ref_lu) >= 200) {
          combs <- list()
          set.seed(19)
          for (q in 1:200) {
            combs[[q]] <- sort(sample(temp_data$Plot_id, n_sites_ref_lu))
          }
          combs <- unique(combs)
        }
        
        # randomly sample 100 combinations
        if (length(combs) > 100) {    
          set.seed(19)
          combs <- sample(combs, size = 100)
        }
        
        # calculating the areas
        areas <- combs |>
          purrr::map_dbl(.f = function(comb) {
            temp_data |>
              dplyr::filter(Plot_id %in% comb) |> 
              sf::st_union() |>
              sf::st_convex_hull() |> 
              sf::st_cast("POLYGON") |>
              sf::st_area()
          })
        
        # Drop the geometry column
        temp_data <- temp_data |> 
          sf::st_drop_geometry()
        # purrr::map_dbl(
        #   .x = combs[areas |> dplyr::between(left = threshold * 0.75, right = threshold * 1.25)],
        #   .f = function(comb) {
        #     mean(temp_data[comb, 28L:109L])
        #   })
        combs <- combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                               right = threshold_area * 1.5)]
        

        for (comb in combs) {
          comm2 <- temp_data |> 
            dplyr::filter(Plot_id %in% comb) |> 
            dplyr::select(-(1:13),)
          soma2 <- (colSums(comm2))
          tsoma2 <- t(soma2)
          N2 <- rowSums(tsoma2)
          NS <- rbind(na.omit(NS),N2)
          C.comm2 <-C_target(comm2,factor = 2)
          C.coms <- rbind(na.omit(C.coms), C.comm2)
        }
      }
    }
  }
  
  NS.min <- t(as.data.frame(min(NS)))
  C.coms.min <- t(as.data.frame(min(C.coms)))
  Dataset_id <- distinct(Dataset,Dataset_id)
  SC[[m]] <- cbind.data.frame(Dataset_id,NS.min,C.coms.min)
}

SC <- data.frame(matrix(unlist(SC), nrow=67, byrow=TRUE),stringsAsFactors=FALSE)

colnames(SC) <-c("Dataset_id","N.min","Target.c")

# create list for saving data
save_all_data <- list()

######## start the loop at study level to calculate all metrics#########
for (m in 1:length(all.files)) {
  print(all.files[[m]])
  stu1 <- noquote(all.files[[m]])
  Dataset <- get(stu1, 1)
  dataset_output <- list()
  add_meta_min <- list()
  Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
  stu4 <- cbind.data.frame(Block, Dataset)
  spl.block <- split(stu4, Block) #spliting in datasets according to block
  blocos <- summary(as.factor(Block))
  n.blocks <- length(blocos)
  C.min <- SC[m, "Target.c"]
  N.min <-  SC[m, "N.min"]
  # start the loop between blocks
  for (t in 1:n.blocks) {
    Dataset <- spl.block[[t]] 
    
    Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
      sf::st_set_crs(4326)
    
    lu_area <- dplyr::bind_cols(
      Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
      area = Dataset |>
        dplyr::group_by(Land_use) |>
        dplyr::summarise() |>
        sf::st_convex_hull() |>
        sf::st_area() |>
        as.numeric()
    )
    #plots number counting
    sites_number <- Dataset |> 
      group_by(Land_use) |> 
      summarise(length(Plot_id)) |> 
      sf::st_drop_geometry()
    
    # check if min area site also have the least sites
    sites_min <- sites_number[which.min(sites_number$`length(Plot_id)`), "Land_use"] |> dplyr::pull()
    threshold_sites <- min(sites_number$`length(Plot_id)`)
    
    threshold_lu <- lu_area[which.min(lu_area$area), "Land_use"] |> dplyr::pull() 
    threshold_area <- min(lu_area$area)
    # only do next calculation if min area land use same with also min sites land use
    if (threshold_sites > 2 & sites_min == threshold_lu) {
      
      #calculate the metrics in min_area_land_use
      threshold_data <- Dataset |>
        dplyr::filter(Land_use == threshold_lu) |> 
        sf::st_drop_geometry()
      meta_info_min <- threshold_data |>
        slice(1) |>
        select(Dataset_id,Block, Spatial_block,
               Temporal_block, Plot_id, Replicate_id,
               Taxa, Continent, Land_use) |> 
        sf::st_drop_geometry()
      
      comm <- threshold_data |> 
        dplyr::select(-(1:13))
      soma <- (colSums(comm))
      tsoma <- t(soma)
      gamma_N <- rowSums(tsoma)
      gamma_S <- sum(tsoma > 0) # sum only TRUE (which means occurence)
      gamma_Spie <- diversity(tsoma, index = "invsimpson")
      if (is.na(N.min)) {
        gamma_Sn <- NA
      } else {
        gamma_Sn <- as.vector(rarefy(tsoma,N.min))
      }
      alpha_S <-  rowSums(comm > 0) 
      beta_S <- gamma_S/mean(alpha_S) 
      alpha_Spie <- diversity(comm, index = "invsimpson")
      beta_Spie <- gamma_Spie/mean(alpha_Spie) 
      alpha_N <- rowSums(comm)
      min_N <- min(c(alpha_N)) # minimum N in the pair of streams
      alpha_Sn <- rarefy(comm, sample = min_N)
      gamma_Sn.b <- rarefy(tsoma, sample = min_N) 
      beta_Sn <- gamma_Sn.b/mean(alpha_Sn)
      #beta.c 
      if (is.na(C.min)) {
        bet.c <- "NA"
      } else if (C.min < 0.5) { #define at least 0.5 target_C
        # C.min <- 0.5
        betac.calcu <-beta_C(comm,C.min,extrapolation = T,interrupt = T)
        bet.c <- c(betac.calcu)#get one vector
      } else {
        betac.calcu <-beta_C(comm,C.min,extrapolation = T,interrupt = T)
        bet.c <- c(betac.calcu)#get one vector
      }
      #beta.spie Tory function
      beta_Spie.T<-beta_SPIE(comm)
      #save data at this min area land-use
      metrics_min <- matrix(ncol = 9)
      metrics_min <- t(c(gamma_N, gamma_S, gamma_Spie,gamma_Sn, beta_S, 
                         beta_Sn, beta_Spie,beta_Spie.T,bet.c))
      colnames(metrics_min) <- c("gamma_N", "gamma_S", "gamma_Spie","gamma_Sn",
                                 "beta_S","beta_Sn","beta_Spie","beta_Spie.T","beta_C")
      
      add_meta_min[[t]] <- cbind(meta_info_min,metrics_min)
      
      #prepare for the next step
      n_sites_ref_lu <- Dataset |>
        dplyr::filter(Land_use == threshold_lu) |>
        nrow()
      
      result <- Dataset |>
        dplyr::select(Plot_id, Land_use) |> 
        sf::st_drop_geometry()
      
      #calculate the metric in the other land-use categories
      add_meta <- list()
      
      for (i in lu_area |>
           dplyr::filter(area > threshold_area) |>
           dplyr::select(Land_use) |>
           dplyr::pull()) {
        
        temp_data <- Dataset |>
          dplyr::filter(Land_use == i)
        meta_info <- temp_data |>
          slice(1) |>
          select(Dataset_id,Block,Spatial_block,
                 Temporal_block,Plot_id, Replicate_id,
                 Taxa,Continent,Land_use) |> 
          sf::st_drop_geometry()
        
        n_sites_i_lu <- nrow(temp_data)
        
        # All combinations of sites down to the number of sites the smallest lu has
        if (choose(length(temp_data$Plot_id), n_sites_ref_lu) < 200) {
          combs <- combn(temp_data$Plot_id,
                         n_sites_ref_lu,
                         simplify = FALSE)
        }
        
        if (choose(length(temp_data$Plot_id),  n_sites_ref_lu) >= 200) {
          combs <- list()
          set.seed(19)
          for (q in 1:200) {
            combs[[q]] <- sort(sample(temp_data$Plot_id, n_sites_ref_lu))
          }
          combs <- unique(combs)
        }
        
        # randomly sample 100 combinations
        if (length(combs) > 100) {    
          set.seed(19)
          combs <- sample(combs, size = 100)
        }
        
        # calculating the areas
        areas <- combs |>
          purrr::map_dbl(.f = function(comb) {
            temp_data |>
              dplyr::filter(Plot_id %in% comb) |> 
              sf::st_union() |>
              sf::st_convex_hull() |> 
              sf::st_cast("POLYGON") |>
              sf::st_area()
          })
        
        # Drop the geometry column
        temp_data <- temp_data |> 
          sf::st_drop_geometry()
        # purrr::map_dbl(
        #   .x = combs[areas |> dplyr::between(left = threshold * 0.75, right = threshold * 1.25)],
        #   .f = function(comb) {
        #     mean(temp_data[comb, 28L:109L])
        #   })
        combs <- combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                               right = threshold_area * 1.5)]
        
        save.metrics <- matrix(ncol = 9)
        colnames(save.metrics) <- c("gamma_N", "gamma_S", "gamma_Spie","gamma_Sn",
                                    "beta_S","beta_Sn","beta_Spie","beta_Spie.T","beta_C")
        for (comb in combs) {
          comm <- temp_data |> 
            dplyr::filter(Plot_id %in% comb) |> 
            dplyr::select(-(1:13),)
          # S = vegan::specnumber(comm)
          soma <- (colSums(comm))
          tsoma <- t(soma)
          gamma_N <- rowSums(tsoma)
          gamma_S <- sum(tsoma > 0) # sum only TRUE (which means occurence)
          gamma_Spie <- diversity(tsoma, index = "invsimpson")
          if (is.na(N.min)) {
            gamma_Sn <- "NA"
          } else {
            gamma_Sn <- as.vector(rarefy(tsoma,N.min))
          }
          alpha_S <-  rowSums(comm > 0) 
          beta_S <- gamma_S/mean(alpha_S) 
          alpha_Spie <- diversity(comm, index = "invsimpson")
          beta_Spie <- gamma_Spie/mean(alpha_Spie) 
          alpha_N <- rowSums(comm)
          min_N <- min(c(alpha_N)) # minimum N in the pair of streams
          alpha_Sn <- rarefy(comm, sample = min_N)
          gamma_Sn.b <- rarefy(tsoma, sample = min_N) 
          beta_Sn <- gamma_Sn.b/mean(alpha_Sn)
          #beta.c 
          if (is.na(C.min)) {
            bet.c <- "NA"
          } else if (C.min < 0.5) { #define at least 0.5 target_C
            # C.min <- 0.5
            betac.calcu <-beta_C(comm,C.min,extrapolation = T,interrupt = T)
            bet.c <- c(betac.calcu)#get one vector
          } else {
            betac.calcu <-beta_C(comm,C.min,extrapolation = T,interrupt = T)
            bet.c <- c(betac.calcu)#get one vector
          }
          #beta.spie Tory function
          beta_Spie.T<-beta_SPIE(comm)
          
          metrics <- c(gamma_N, gamma_S, gamma_Spie,gamma_Sn, beta_S, 
                       beta_Sn, beta_Spie,beta_Spie.T,bet.c)
          save.metrics <- rbind(na.omit(save.metrics),metrics)
        }
        metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
        add_meta[[i]] <- cbind(meta_info,metrics_ave)
      }
      dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta_min,add_meta)))
    }
  }
  save_all_data[[m]] <- remove_rownames(bind_rows(c(dataset_output)))
}
metrics_loop <- remove_rownames(bind_rows(c(save_all_data)))

#######For dataset_9###########

#NV has the least area but not the minimum sites number.
#can choose the least number of sites from every land-use and then cut each
#land-use into minimum land-use area.
Dataset <- data.frame(Dataset_9.txt,stringsAsFactors = FALSE)
Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <- split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)

NS <- matrix(ncol = 1)
colnames(NS) <- c("NS")

C.coms <- matrix(ncol = 1)
colnames(C.coms) <- c("C.coms")

# start the loop between blocks for N.min and C.min
for (t in 1:n.blocks) {
  #plots number counting
  
  Dataset <- spl.block[[t]] 
  
  Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(4326)
  
  sites_number <- Dataset |> 
    group_by(Land_use) |> 
    summarise(length(Plot_id)) |> 
    sf::st_drop_geometry()
  
  threshold_sites <- min(sites_number$`length(Plot_id)`)
  
  lu_area <- dplyr::bind_cols(
    Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
    area = Dataset |>
      dplyr::group_by(Land_use) |>
      dplyr::summarise() |>
      sf::st_convex_hull() |>
      sf::st_area() |>
      as.numeric()
  )
  
  threshold_area <- min(lu_area$area)
  
  threshold_sites <- (threshold_sites)-1 #minus one, so the land use with least
  #sites number can be more flexiable

  for (i in unique(Dataset$Land_use)) {
    temp_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |>
      sf::st_drop_geometry()
   
    # All combinations of sites down to the number of sites the smallest lu has
    if (choose(length(temp_data$Plot_id), threshold_sites) < 200) {
      combs <- combn(temp_data$Plot_id,
                     threshold_sites,
                     simplify = FALSE)
    }
    
    if (choose(length(temp_data$Plot_id),  threshold_sites) >= 200) {
      combs <- list()
      set.seed(19)
      for (q in 1:200) {
        combs[[q]] <- sort(sample(temp_data$Plot_id,threshold_sites))
      }
      combs <- unique(combs)
    }
    
    # randomly sample 100 combinations
    if (length(combs) > 100) {
      set.seed(19)
      combs <- sample(combs, size = 100)
    }
    
    # calculating the areas
    areas <- combs |>
      purrr::map_dbl(
        .f = function(comb) {
          temp_data |>
            dplyr::filter(Plot_id %in% comb) |>
            sf::st_union() |>
            sf::st_convex_hull() |>
            sf::st_cast("POLYGON") |>
            sf::st_area()
        }
      )
    
    # Drop the geometry column
    temp_data <- temp_data |>
      sf::st_drop_geometry()

    combs <-
      combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                    right = threshold_area * 1.5)]
    
    for (comb in combs) {
      comm2 <- temp_data |>
        dplyr::filter(Plot_id %in% comb) |>
        dplyr::select(-(1:13), )
      comm2 <- comm2[rowSums(comm2) >= 5,] 
      soma2 <- (colSums(comm2))
      tsoma2 <- t(soma2)
      N2 <- rowSums(tsoma2)
      NS <- rbind(na.omit(NS), N2)
      C.comm2 <- C_target(comm2, factor = 2)
      C.coms <- rbind(na.omit(C.coms), C.comm2)
    }
  }
}

N.min <- t(as.data.frame(min(NS)))
C.min <- t(as.data.frame(min(C.coms)))
#give C.min value 0.5 since it is below 0.5
# C.min <- 0.5

#calculate the metrics
dataset_output <- list()
for (t in 1:n.blocks) {
  #plots number countint
  Dataset <- spl.block[[t]] 
  Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(4326)
  sites_number <- Dataset |> 
    group_by(Land_use) |> 
    summarise(length(Plot_id)) |> 
    sf::st_drop_geometry()
  threshold_sites <- min(sites_number$`length(Plot_id)`)
  lu_area <- dplyr::bind_cols(
    Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
    area = Dataset |>
      dplyr::group_by(Land_use) |>
      dplyr::summarise() |>
      sf::st_convex_hull() |>
      sf::st_area() |>
      as.numeric()
  )
  threshold_area <- min(lu_area$area)
  threshold_sites <- (threshold_sites)-1 
  #minus one, so the land use with least sites number can be more flexiable
  add_meta <- list()
  for (i in unique(Dataset$Land_use)) {
    temp_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Block,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |>
      sf::st_drop_geometry()
    
    # All combinations of sites down to the number of sites the smallest lu has
    if (choose(length(temp_data$Plot_id), threshold_sites) < 200) {
      combs <- combn(temp_data$Plot_id,
                     threshold_sites,
                     simplify = FALSE)
    }
    
    if (choose(length(temp_data$Plot_id),  threshold_sites) >= 200) {
      combs <- list()
      set.seed(19)
      for (q in 1:200) {
        combs[[q]] <- sort(sample(temp_data$Plot_id,threshold_sites))
      }
      combs <- unique(combs)
    }
    
    # randomly sample 100 combinations
    if (length(combs) > 100) {
      set.seed(19)
      combs <- sample(combs, size = 100)
    }
    
    # calculating the areas
    areas <- combs |>
      purrr::map_dbl(
        .f = function(comb) {
          temp_data |>
            dplyr::filter(Plot_id %in% comb) |>
            sf::st_union() |>
            sf::st_convex_hull() |>
            sf::st_cast("POLYGON") |>
            sf::st_area()
        }
      )
    # Drop the geometry column
    temp_data <- temp_data |>
      sf::st_drop_geometry()
    combs <-
      combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                    right = threshold_area * 1.5)]
    for (comb in combs) {
      comm <- temp_data |> 
        dplyr::filter(Plot_id %in% comb) |> 
        dplyr::select(-(1:13),) 
      comm <- comm[rowSums(comm) >= 5,] #if rowsums < 5, can not calculate betaC
      soma <- (colSums(comm))
      tsoma <- t(soma)
      gamma_N <- rowSums(tsoma)
      gamma_S <- sum(tsoma > 0) # sum only TRUE (which means occurence)
      gamma_Spie <- diversity(tsoma, index = "invsimpson")
      N.min <- as.numeric(N.min)
      gamma_Sn <- as.vector(rarefy(tsoma,N.min))
      alpha_S <-  rowSums(comm > 0) 
      beta_S <- gamma_S/mean(alpha_S) 
      alpha_Spie <- diversity(comm, index = "invsimpson")
      beta_Spie <- gamma_Spie/mean(alpha_Spie) 
      alpha_N <- rowSums(comm)
      min_N <- min(c(alpha_N)) # minimum N in the pair of streams
      alpha_Sn <- rarefy(comm, sample = min_N)
      gamma_Sn.b <- rarefy(tsoma, sample = min_N) 
      beta_Sn <- gamma_Sn.b/mean(alpha_Sn)
      #beta.c 
      betac.calcu <-beta_C(comm,C.min,extrapolation = T,interrupt = T)
      bet.c <- c(betac.calcu)#get one vector
      #beta.spie Tory function
      beta_Spie.T<-beta_SPIE(comm)
      metrics <- c(gamma_N, gamma_S, gamma_Spie,gamma_Sn, beta_S, 
                   beta_Sn, beta_Spie,beta_Spie.T,bet.c)
      save.metrics <- rbind(na.omit(save.metrics),metrics)
    }
    metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
    add_meta[[i]] <- cbind(meta_info,metrics_ave)
  }
  dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta)))
}
metrics_9 <- remove_rownames(bind_rows(c(dataset_output)))

#######For dataset_17###########
#first three urban sites as one gamma group

Dataset <- data.frame(Dataset_17.txt,stringsAsFactors = FALSE)
Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <- split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)

NS <- matrix(ncol = 1)
colnames(NS) <- c("NS")

C.coms <- matrix(ncol = 1)
colnames(C.coms) <- c("C.coms")

# start the loop between blocks for N.min and C.min
for (t in 1:n.blocks) {
  #plots number counting
  Dataset <- spl.block[[t]]
  for (i in seq(1, nrow(Dataset), by = 3)) {
    temp_data <- Dataset[i+0:2,] |>
    sf::st_drop_geometry()
    comm2 <- temp_data |>
      dplyr::select(-(1:15),)
    soma2 <- (colSums(comm2))
    tsoma2 <- t(soma2)
    N2 <- rowSums(tsoma2)
    NS <- rbind(na.omit(NS), N2)
    C.comm2 <- C_target(comm2, factor = 2)
    C.coms <- rbind(na.omit(C.coms), C.comm2)
  }
}
N.min <- as.numeric(t(as.data.frame(min(NS))))
C.min <- as.numeric(t(as.data.frame(min(C.coms))))

#calculate the metrics
dataset_output <- list()
for (t in 1:n.blocks) {
  #plots number counting
  Dataset <- spl.block[[t]]
  add_meta <- list()
  for (i in Dataset$Land_use) {
    save.metrics <- matrix(ncol = 9)
    colnames(save.metrics) <- c("gamma_N", "gamma_S", "gamma_Spie","gamma_Sn",
                                "beta_S","beta_Sn","beta_Spie","beta_Spie.T","beta_C")
    lu_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Block,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |> sf::st_drop_geometry()
    for (q in seq(1, nrow(lu_data), by = 3)) {
      temp_data <- lu_data[q+0:2,] |>
        sf::st_drop_geometry()
      comm <- temp_data |>
        dplyr::select(-(1:15),)
      soma <- (colSums(comm))
      tsoma <- t(soma)
      gamma_N <- rowSums(tsoma)
      gamma_S <- sum(tsoma > 0) # sum only TRUE (which means occurence)
      gamma_Spie <- diversity(tsoma, index = "invsimpson")
      N.min <- as.numeric(N.min)
      gamma_Sn <- as.vector(rarefy(tsoma,N.min))
      alpha_S <-  rowSums(comm > 0) 
      beta_S <- gamma_S/mean(alpha_S) 
      alpha_Spie <- diversity(comm, index = "invsimpson")
      beta_Spie <- gamma_Spie/mean(alpha_Spie) 
      alpha_N <- rowSums(comm)
      min_N <- min(c(alpha_N)) # minimum N in the pair of streams
      alpha_Sn <- rarefy(comm, sample = min_N)
      gamma_Sn.b <- rarefy(tsoma, sample = min_N) 
      beta_Sn <- gamma_Sn.b/mean(alpha_Sn)
      #beta.c 
      betac.calcu <-beta_C(comm,C.min,extrapolation = T,interrupt = T)
      bet.c <- c(betac.calcu)#get one vector
      #beta.spie Tory function
      beta_Spie.T<-beta_SPIE(comm)
      metrics <- c(gamma_N, gamma_S, gamma_Spie,gamma_Sn, beta_S, 
                   beta_Sn, beta_Spie,beta_Spie.T,bet.c)
      save.metrics <- rbind(na.omit(save.metrics),metrics)
    }
    metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
    add_meta[[i]] <- cbind(meta_info,metrics_ave)
  }
  dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta)))
}
metrics_17<-remove_rownames(bind_rows(c(dataset_output)))




#######For dataset_33########
Dataset <- Dataset_33.txt

block <- Dataset |>
  group_by(Spatial_block,Land_use) |>
  summarise(Count = n()) |>
  sf::st_drop_geometry()

sites_number <- Dataset |>
  group_by(Spatial_block, Land_use) |>
  summarise(length(Plot_id)) |>
  sf::st_drop_geometry()
Block <-
  paste(Dataset$Spatial_block, Dataset$Temporal_block, sep = "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <-
  split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)
#for block 1: can choose 
Dataset <- spl.block[[1]] |> 
  slice(3:14)

Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
  sf::st_set_crs(4326)

threshold_sites <- 4 #see if we can cut 5 sites with similar spatial extent
#tested 5 doesn't work, so we try 4 now. And we found some could work
area <- list()
#calculate and choose the least area
for (i in Dataset$Land_use) {
  temp_data <- Dataset |>
    dplyr::filter(Land_use == i)
  
  # All combinations of sites down to the number of sites the smallest lu has
  combs <- combn(temp_data$Plot_id,
                   threshold_sites,
                   simplify = FALSE)
  # calculating the areas
  areas <- combs |>
    purrr::map_dbl(
      .f = function(comb) {
        temp_data |>
          dplyr::filter(Plot_id %in% comb) |>
          sf::st_union() |>
          sf::st_convex_hull() |>
          sf::st_cast("POLYGON") |>
          sf::st_area()
      }
    )
  area[[i]] <- areas
}

threshold_area <- min(area[[1]])
#the min in agriculture can compare with some combination in urban
#Roel: This is manual selection. Because we choose whatever is best for us.

NS <- matrix(ncol = 1)
colnames(NS) <- c("NS")
C.coms <- matrix(ncol = 1)
colnames(C.coms) <- c("C.coms")
for (i in Dataset$Land_use) {
  temp_data <- Dataset |>
    dplyr::filter(Land_use == i)
  # All combinations of sites down to the number of sites the smallest lu has
    combs <- combn(temp_data$Plot_id,
                   threshold_sites,
                   simplify = FALSE)
  # calculating the areas
  areas <- combs |>
    purrr::map_dbl(
      .f = function(comb) {
        temp_data |>
          dplyr::filter(Plot_id %in% comb) |>
          sf::st_union() |>
          sf::st_convex_hull() |>
          sf::st_cast("POLYGON") |>
          sf::st_area()
      }
    )
  # Drop the geometry column
  temp_data <- temp_data |>
    sf::st_drop_geometry()
  combs <-
    combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                  right = threshold_area * 1.5)]
  for (comb in combs) {
    comm2 <- temp_data |>
      dplyr::filter(Plot_id %in% comb) |>
      dplyr::select(-(1:13),)
    soma2 <- (colSums(comm2))
    tsoma2 <- t(soma2)
    N2 <- rowSums(tsoma2)
    NS <- rbind(na.omit(NS), N2)
    C.comm2 <- C_target(comm2, factor = 2)
    C.coms <- rbind(na.omit(C.coms), C.comm2)
  }
}
N.min <- as.numeric(t(as.data.frame(min(NS))))
C.min <- as.numeric(t(as.data.frame(min(C.coms))))

add_meta <- list()
for (i in Dataset$Land_use) {
  save.metrics <- matrix(ncol = 9)
  colnames(save.metrics) <- c("gamma_N", "gamma_S", "gamma_Spie","gamma_Sn",
                              "beta_S","beta_Sn","beta_Spie","beta_Spie.T","beta_C")
   temp_data <- Dataset |>
    dplyr::filter(Land_use == i)
  meta_info <- temp_data |>
    slice(1) |>
    select(
      Dataset_id,
      Block,
      Spatial_block,
      Temporal_block,
      Plot_id, Replicate_id,
      Taxa,
      Continent,
      Land_use
    ) |>
    sf::st_drop_geometry()
  
  # All combinations of sites down to the number of sites the smallest lu has
    combs <- combn(temp_data$Plot_id,
                   threshold_sites,
                   simplify = FALSE)
  # calculating the areas
  areas <- combs |>
    purrr::map_dbl(
      .f = function(comb) {
        temp_data |>
          dplyr::filter(Plot_id %in% comb) |>
          sf::st_union() |>
          sf::st_convex_hull() |>
          sf::st_cast("POLYGON") |>
          sf::st_area()
      }
    )
  # Drop the geometry column
  temp_data <- temp_data |>
    sf::st_drop_geometry()
  combs <-
    combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                  right = threshold_area * 1.5)]
  for (comb in combs) {
    comm <- temp_data |>
      dplyr::filter(Plot_id %in% comb) |>
      dplyr::select(-(1:13), )
    # S = vegan::specnumber(comm)
    # comm <- comm[rowSums(comm) >= 5,] #if rowsums < 5, can not calculate betaC
    soma <- (colSums(comm))
    tsoma <- t(soma)
    gamma_N <- rowSums(tsoma)
    gamma_S <-
      sum(tsoma > 0) # sum only TRUE (which means occurence)
    gamma_Spie <- diversity(tsoma, index = "invsimpson")
    N.min <- as.numeric(N.min)
    gamma_Sn <- as.vector(rarefy(tsoma, N.min))
    alpha_S <-  rowSums(comm > 0)
    beta_S <- gamma_S / mean(alpha_S)
    alpha_Spie <- diversity(comm, index = "invsimpson")
    beta_Spie <- gamma_Spie / mean(alpha_Spie)
    alpha_N <- rowSums(comm)
    min_N <- min(c(alpha_N)) # minimum N in the pair of streams
    alpha_Sn <- rarefy(comm, sample = min_N)
    gamma_Sn.b <- rarefy(tsoma, sample = min_N)
    beta_Sn <- gamma_Sn.b / mean(alpha_Sn)
    #beta.c
    betac.calcu <-
      beta_C(comm,
             C.min,
             extrapolation = T,
             interrupt = T)
    bet.c <- c(betac.calcu)#get one vector
    #beta.spie Tory function
    beta_Spie.T <- beta_SPIE(comm)
    metrics <- c(
      gamma_N,
      gamma_S,
      gamma_Spie,
      gamma_Sn,
      beta_S,
      beta_Sn,
      beta_Spie,
      beta_Spie.T,
      bet.c
    )
    save.metrics <- rbind(na.omit(save.metrics), metrics)
  }
  metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
  add_meta[[i]] <- cbind(meta_info, metrics_ave)
}
metrics_33 <- remove_rownames(bind_rows(c(add_meta)))



#######For dataset_36,37,39,40,49,172,175,176,178########

subdata <- all.files[c(10,36,37,39,40,46,60,61,63)]

SC <- list()
for (m in 1:length(subdata)) {
  print(subdata[[m]])
  stu1 <- noquote(subdata[[m]])
  Dataset <- get(stu1, 1)
  Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
  stu4 <- cbind.data.frame(Block, Dataset)
  spl.block <- split(stu4, Block) #spliting in datasets according to block
  blocos <- summary(as.factor(Block))
  n.blocks <- length(blocos)
  NS <- matrix(ncol = 1)
  colnames(NS) <- c("NS")
  C.coms <- matrix(ncol = 1)
  colnames(C.coms) <- c("C.coms")
  # start the loop between blocks for N.min and C.min
  for (t in 1:n.blocks) {
    #plots number counting
    Dataset <- spl.block[[t]] 
    Dataset <- sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
      sf::st_set_crs(4326)
    sites_number <- Dataset |> 
      group_by(Land_use) |> 
      summarise(length(Plot_id)) |> 
      sf::st_drop_geometry()
    threshold_sites <- min(sites_number$`length(Plot_id)`)
    lu_area <- dplyr::bind_cols(
      Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
      area = Dataset |>
        dplyr::group_by(Land_use) |>
        dplyr::summarise() |>
        sf::st_convex_hull() |>
        sf::st_area() |>
        as.numeric()
    )
    threshold_area <- min(lu_area$area)
    threshold_sites <- (threshold_sites)-1 #minus one, so the land use with least
    #sites number can be more flexiable
    for (i in unique(Dataset$Land_use)) {
      temp_data <- Dataset |>
        dplyr::filter(Land_use == i)
      meta_info <- temp_data |>
        slice(1) |>
        select(
          Dataset_id,
          Spatial_block,
          Temporal_block,
          Plot_id, Replicate_id,
          Taxa,
          Continent,
          Land_use
        ) |>
        sf::st_drop_geometry()
        combs <- list()
        set.seed(19)
        for (q in 1:200) {
          combs[[q]] <- sort(sample(temp_data$Plot_id,threshold_sites))
        }
        combs <- unique(combs)
      # randomly sample 100 combinations
        if (length(combs) > 100) {    
          set.seed(19)
          combs <- sample(combs, size = 100)
        }
      # calculating the areas
      areas <- combs |>
        purrr::map_dbl(
          .f = function(comb) {
            temp_data |>
              dplyr::filter(Plot_id %in% comb) |>
              sf::st_union() |>
              sf::st_convex_hull() |>
              sf::st_cast("POLYGON") |>
              sf::st_area()
          }
        )
      # Drop the geometry column
      temp_data <- temp_data |>
        sf::st_drop_geometry()
      #filter the area that meet our criteria
      combs <-
        combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                      right = threshold_area * 1.5)]
      for (comb in combs) {
        comm2 <- temp_data |>
          dplyr::filter(Plot_id %in% comb) |>
          dplyr::select(-(1:13), )
        soma2 <- (colSums(comm2))
        tsoma2 <- t(soma2)
        N2 <- rowSums(tsoma2)
        NS <- rbind(na.omit(NS), N2)
        C.comm2 <- C_target(comm2, factor = 2)
        C.coms <- rbind(na.omit(C.coms), C.comm2)
      }
    }
  }
  NS.min <- t(as.data.frame(min(NS)))
  C.coms.min <- t(as.data.frame(min(C.coms)))
  Dataset_id <- distinct(Dataset,Dataset_id)
  SC[[m]] <- cbind.data.frame(Dataset_id,NS.min,C.coms.min)
}

SC <- data.frame(matrix(unlist(SC), nrow=9, byrow=TRUE),stringsAsFactors=FALSE)

colnames(SC) <-c("Dataset_id","N.min","Target.c")

#calculate the metrics
save_all_data <- list()
for (m in 1:length(subdata)) {
  print(subdata[[m]])
  stu1 <- noquote(subdata[[m]])
  Dataset <- get(stu1, 1)
  Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
  stu4 <- cbind.data.frame(Block, Dataset)
  spl.block <- split(stu4, Block) #spliting in datasets according to block
  blocos <- summary(as.factor(Block))
  n.blocks <- length(blocos)
  C.min <- SC[m, "Target.c"]
  N.min <-  SC[m, "N.min"]
  dataset_output <- list()
  for (t in 1:n.blocks) {
    Dataset <- spl.block[[t]]
    Dataset <-
      sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
      sf::st_set_crs(4326)
    sites_number <- Dataset |>
      group_by(Land_use) |>
      summarise(length(Plot_id)) |>
      sf::st_drop_geometry()
    threshold_sites <- min(sites_number$`length(Plot_id)`)
    lu_area <- dplyr::bind_cols(
      Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
      area = Dataset |>
        dplyr::group_by(Land_use) |>
        dplyr::summarise() |>
        sf::st_convex_hull() |>
        sf::st_area() |>
        as.numeric()
    )
    threshold_area <- min(lu_area$area)
    threshold_sites <- (threshold_sites) - 1
    #minus one, so the land use with least sites number can be more flexiable
    add_meta <- list()
    for (i in unique(Dataset$Land_use)) {
      save.metrics <- matrix(ncol = 9)
      colnames(save.metrics) <- c("gamma_N", "gamma_S", "gamma_Spie","gamma_Sn",
                                  "beta_S","beta_Sn","beta_Spie","beta_Spie.T","beta_C")
      temp_data <- Dataset |>
        dplyr::filter(Land_use == i)
      meta_info <- temp_data |>
        slice(1) |>
        select(
          Dataset_id,
          Block,
          Spatial_block,
          Temporal_block,
          Plot_id, Replicate_id,
          Taxa,
          Continent,
          Land_use
        ) |>
        sf::st_drop_geometry()
      # All combinations of sites down to the number of sites the smallest lu has
        combs <- list()
        set.seed(19)
        for (q in 1:200) {
          combs[[q]] <- sort(sample(temp_data$Plot_id, threshold_sites))
        }
        combs <- unique(combs)

      # randomly sample 100 combinations
        if (length(combs) > 100) {    
          set.seed(19)
          combs <- sample(combs, size = 100)
        }
      
      # calculating the areas
      areas <- combs |>
        purrr::map_dbl(
          .f = function(comb) {
            temp_data |>
              dplyr::filter(Plot_id %in% comb) |>
              sf::st_union() |>
              sf::st_convex_hull() |>
              sf::st_cast("POLYGON") |>
              sf::st_area()
          }
        )
      # Drop the geometry column
      temp_data <- temp_data |>
        sf::st_drop_geometry()
      combs <-
        combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                      right = threshold_area * 1.5)]
      for (comb in combs) {
        comm <- temp_data |>
          dplyr::filter(Plot_id %in% comb) |>
          dplyr::select(-(1:13), )
        # S = vegan::specnumber(comm)
        # comm <- comm[rowSums(comm) >= 5,] #if rowsums < 5, can not calculate betaC
        soma <- (colSums(comm))
        tsoma <- t(soma)
        gamma_N <- rowSums(tsoma)
        gamma_S <-
          sum(tsoma > 0) # sum only TRUE (which means occurence)
        gamma_Spie <- diversity(tsoma, index = "invsimpson")
        N.min <- as.numeric(N.min)
        gamma_Sn <- as.vector(rarefy(tsoma, N.min))
        alpha_S <-  rowSums(comm > 0)
        beta_S <- gamma_S / mean(alpha_S)
        alpha_Spie <- diversity(comm, index = "invsimpson")
        beta_Spie <- gamma_Spie / mean(alpha_Spie)
        alpha_N <- rowSums(comm)
        min_N <- min(c(alpha_N)) # minimum N in the pair of streams
        alpha_Sn <- rarefy(comm, sample = min_N)
        gamma_Sn.b <- rarefy(tsoma, sample = min_N)
        beta_Sn <- gamma_Sn.b / mean(alpha_Sn)
        #beta.c
        betac.calcu <-
          beta_C(comm,
                 C.min,
                 extrapolation = T,
                 interrupt = T)
        bet.c <- c(betac.calcu)#get one vector
        #beta.spie Tory function
        beta_Spie.T <- beta_SPIE(comm)
        metrics <- c(
          gamma_N,
          gamma_S,
          gamma_Spie,
          gamma_Sn,
          beta_S,
          beta_Sn,
          beta_Spie,
          beta_Spie.T,
          bet.c
        )
        save.metrics <- rbind(na.omit(save.metrics), metrics)
      }
      metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
      add_meta[[i]] <- cbind(meta_info, metrics_ave)
    }
    dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta)))
  }
  save_all_data[[m]] <- remove_rownames(bind_rows(c(dataset_output)))
}
metrics_group <- remove_rownames(bind_rows(c(save_all_data)))

#######For dataset_38##########
Dataset <- Dataset_38.txt
Block <-
  paste(Dataset$Spatial_block, Dataset$Temporal_block, sep = "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <-
  split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)
NS <- matrix(ncol = 1)
colnames(NS) <- c("NS")
C.coms <- matrix(ncol = 1)
colnames(C.coms) <- c("C.coms")
# start the loop between blocks for N.min and C.min
for (t in 1:n.blocks) {
  #plots number counting
  Dataset <- spl.block[[t]]
  Dataset <-
    sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(4326)
  sites_number <- Dataset |>
    group_by(Land_use) |>
    summarise(length(Plot_id)) |>
    sf::st_drop_geometry()
  threshold_sites <- min(sites_number$`length(Plot_id)`)
  lu_area <- dplyr::bind_cols(
    Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
    area = Dataset |>
      dplyr::group_by(Land_use) |>
      dplyr::summarise() |>
      sf::st_convex_hull() |>
      sf::st_area() |>
      as.numeric()
  )
  threshold_area <- min(lu_area$area)
  threshold_sites <-
    (threshold_sites) - 15 #minus one, so the land use with least
  #sites number can be more flexiable
  for (i in unique(Dataset$Land_use)) {
    temp_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |>
      sf::st_drop_geometry()
    combs <- list()
    set.seed(19)
    for (q in 1:200) {
      combs[[q]] <- sort(sample(temp_data$Plot_id, threshold_sites))
    }
    combs <- unique(combs)
    # randomly sample 100 combinations
    if (length(combs) > 100) {
      set.seed(19)
      combs <- sample(combs, size = 100)
    }
    # calculating the areas
    areas <- combs |>
      purrr::map_dbl(
        .f = function(comb) {
          temp_data |>
            dplyr::filter(Plot_id %in% comb) |>
            sf::st_union() |>
            sf::st_convex_hull() |>
            sf::st_cast("POLYGON") |>
            sf::st_area()
        }
      )
    # Drop the geometry column
    temp_data <- temp_data |>
      sf::st_drop_geometry()
    #filter the area that meet our criteria
    combs <-
      combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                    right = threshold_area * 1.5)]
    for (comb in combs) {
      comm2 <- temp_data |>
        dplyr::filter(Plot_id %in% comb) |>
        dplyr::select(-(1:13),)
      soma2 <- (colSums(comm2))
      tsoma2 <- t(soma2)
      N2 <- rowSums(tsoma2)
      NS <- rbind(na.omit(NS), N2)
      C.comm2 <- C_target(comm2, factor = 2)
      C.coms <- rbind(na.omit(C.coms), C.comm2)
    }
  }
}
N.min <- t(as.data.frame(min(NS)))
C.min <- t(as.data.frame(min(C.coms)))

#metrics
Dataset <- Dataset_38.txt
Block <-
  paste(Dataset$Spatial_block, Dataset$Temporal_block, sep = "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <-
  split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)

dataset_output <- list()
for (t in 1:n.blocks) {
  Dataset <- spl.block[[t]]
  Dataset <-
    sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(4326)
  sites_number <- Dataset |>
    group_by(Land_use) |>
    summarise(length(Plot_id)) |>
    sf::st_drop_geometry()
  threshold_sites <- min(sites_number$`length(Plot_id)`)
  lu_area <- dplyr::bind_cols(
    Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
    area = Dataset |>
      dplyr::group_by(Land_use) |>
      dplyr::summarise() |>
      sf::st_convex_hull() |>
      sf::st_area() |>
      as.numeric()
  )
  threshold_area <- min(lu_area$area)
  threshold_sites <- (threshold_sites) - 15
  #minus one, so the land use with least sites number can be more flexiable
  add_meta <- list()
  for (i in unique(Dataset$Land_use)) {
    save.metrics <- matrix(ncol = 9)
    colnames(save.metrics) <-
      c(
        "gamma_N",
        "gamma_S",
        "gamma_Spie",
        "gamma_Sn",
        "beta_S",
        "beta_Sn",
        "beta_Spie",
        "beta_Spie.T",
        "beta_C"
      )
    temp_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Block,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |>
      sf::st_drop_geometry()
    # All combinations of sites down to the number of sites the smallest lu has
    combs <- list()
    set.seed(19)
    for (q in 1:200) {
      combs[[q]] <- sort(sample(temp_data$Plot_id, threshold_sites))
    }
    combs <- unique(combs)
    
    # randomly sample 100 combinations
    if (length(combs) > 100) {
      set.seed(19)
      combs <- sample(combs, size = 100)
    }
    
    # calculating the areas
    areas <- combs |>
      purrr::map_dbl(
        .f = function(comb) {
          temp_data |>
            dplyr::filter(Plot_id %in% comb) |>
            sf::st_union() |>
            sf::st_convex_hull() |>
            sf::st_cast("POLYGON") |>
            sf::st_area()
        }
      )
    # Drop the geometry column
    temp_data <- temp_data |>
      sf::st_drop_geometry()
    combs <-
      combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                    right = threshold_area * 1.5)]
    for (comb in combs) {
      comm <- temp_data |>
        dplyr::filter(Plot_id %in% comb) |>
        dplyr::select(-(1:13),)
      # S = vegan::specnumber(comm)
      # comm <- comm[rowSums(comm) >= 5,] #if rowsums < 5, can not calculate betaC
      soma <- (colSums(comm))
      tsoma <- t(soma)
      gamma_N <- rowSums(tsoma)
      gamma_S <-
        sum(tsoma > 0) # sum only TRUE (which means occurence)
      gamma_Spie <- diversity(tsoma, index = "invsimpson")
      N.min <- as.numeric(N.min)
      gamma_Sn <- as.vector(rarefy(tsoma, N.min))
      alpha_S <-  rowSums(comm > 0)
      beta_S <- gamma_S / mean(alpha_S)
      alpha_Spie <- diversity(comm, index = "invsimpson")
      beta_Spie <- gamma_Spie / mean(alpha_Spie)
      alpha_N <- rowSums(comm)
      min_N <- min(c(alpha_N)) # minimum N in the pair of streams
      alpha_Sn <- rarefy(comm, sample = min_N)
      gamma_Sn.b <- rarefy(tsoma, sample = min_N)
      beta_Sn <- gamma_Sn.b / mean(alpha_Sn)
      #beta.c
      betac.calcu <-
        beta_C(comm,
               C.min,
               extrapolation = T,
               interrupt = T)
      bet.c <- c(betac.calcu)#get one vector
      #beta.spie Tory function
      beta_Spie.T <- beta_SPIE(comm)
      metrics <- c(
        gamma_N,
        gamma_S,
        gamma_Spie,
        gamma_Sn,
        beta_S,
        beta_Sn,
        beta_Spie,
        beta_Spie.T,
        bet.c
      )
      save.metrics <- rbind(na.omit(save.metrics), metrics)
    }
    metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
    add_meta[[i]] <- cbind(meta_info, metrics_ave)
  }
  dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta)))
}
metrics_38 <- remove_rownames(bind_rows(c(dataset_output)))



#######For dataset_213##########
Dataset <- Dataset_213.txt
Block <-
  paste(Dataset$Spatial_block, Dataset$Temporal_block, sep = "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <-
  split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)
NS <- matrix(ncol = 1)
colnames(NS) <- c("NS")
C.coms <- matrix(ncol = 1)
colnames(C.coms) <- c("C.coms")
# start the loop between blocks for N.min and C.min
for (t in 1:n.blocks) {
  #plots number counting
  Dataset <- spl.block[[t]]
  Dataset <-
    sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(4326)
  sites_number <- Dataset |>
    group_by(Land_use) |>
    summarise(length(Plot_id)) |>
    sf::st_drop_geometry()
  threshold_sites <- min(sites_number$`length(Plot_id)`)
  lu_area <- dplyr::bind_cols(
    Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
    area = Dataset |>
      dplyr::group_by(Land_use) |>
      dplyr::summarise() |>
      sf::st_convex_hull() |>
      sf::st_area() |>
      as.numeric()
  )
  threshold_area <- min(lu_area$area)
  threshold_sites <-
    (threshold_sites) - 3 #minus one, so the land use with least
  #sites number can be more flexiable
  for (i in unique(Dataset$Land_use)) {
    temp_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |>
      sf::st_drop_geometry()
    combs <- list()
    set.seed(19)
    for (q in 1:200) {
      combs[[q]] <- sort(sample(temp_data$Plot_id, threshold_sites))
    }
    combs <- unique(combs)
    # randomly sample 100 combinations
    if (length(combs) > 100) {
      set.seed(19)
      combs <- sample(combs, size = 100)
    }
    # calculating the areas
    areas <- combs |>
      purrr::map_dbl(
        .f = function(comb) {
          temp_data |>
            dplyr::filter(Plot_id %in% comb) |>
            sf::st_union() |>
            sf::st_convex_hull() |>
            sf::st_cast("POLYGON") |>
            sf::st_area()
        }
      )
    # Drop the geometry column
    temp_data <- temp_data |>
      sf::st_drop_geometry()
    #filter the area that meet our criteria
    combs <-
      combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                    right = threshold_area * 1.5)]
    for (comb in combs) {
      comm2 <- temp_data |>
        dplyr::filter(Plot_id %in% comb) |>
        dplyr::select(-(1:13),)
      soma2 <- (colSums(comm2))
      tsoma2 <- t(soma2)
      N2 <- rowSums(tsoma2)
      NS <- rbind(na.omit(NS), N2)
      C.comm2 <- C_target(comm2, factor = 2)
      C.coms <- rbind(na.omit(C.coms), C.comm2)
    }
  }
}
N.min <- t(as.data.frame(min(NS)))
C.min <- t(as.data.frame(min(C.coms)))

#metrics
Dataset <- Dataset_213.txt
Block <-
  paste(Dataset$Spatial_block, Dataset$Temporal_block, sep = "_")
stu4 <- cbind.data.frame(Block, Dataset)
spl.block <-
  split(stu4, Block) #spliting in datasets according to block
blocos <- summary(as.factor(Block))
n.blocks <- length(blocos)

dataset_output <- list()
for (t in 1:n.blocks) {
  Dataset <- spl.block[[t]]
  Dataset <-
    sf::st_as_sf(Dataset, coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(4326)
  sites_number <- Dataset |>
    group_by(Land_use) |>
    summarise(length(Plot_id)) |>
    sf::st_drop_geometry()
  threshold_sites <- min(sites_number$`length(Plot_id)`)
  lu_area <- dplyr::bind_cols(
    Dataset |> dplyr::group_by(Land_use) |> dplyr::group_keys(),
    area = Dataset |>
      dplyr::group_by(Land_use) |>
      dplyr::summarise() |>
      sf::st_convex_hull() |>
      sf::st_area() |>
      as.numeric()
  )
  threshold_area <- min(lu_area$area)
  threshold_sites <- (threshold_sites) - 3
  #minus one, so the land use with least sites number can be more flexiable
  add_meta <- list()
  for (i in unique(Dataset$Land_use)) {
    save.metrics <- matrix(ncol = 9)
    colnames(save.metrics) <-
      c(
        "gamma_N",
        "gamma_S",
        "gamma_Spie",
        "gamma_Sn",
        "beta_S",
        "beta_Sn",
        "beta_Spie",
        "beta_Spie.T",
        "beta_C"
      )
    temp_data <- Dataset |>
      dplyr::filter(Land_use == i)
    meta_info <- temp_data |>
      slice(1) |>
      select(
        Dataset_id,
        Block,
        Spatial_block,
        Temporal_block,
        Plot_id, Replicate_id,
        Taxa,
        Continent,
        Land_use
      ) |>
      sf::st_drop_geometry()
    # All combinations of sites down to the number of sites the smallest lu has
    combs <- list()
    set.seed(19)
    for (q in 1:200) {
      combs[[q]] <- sort(sample(temp_data$Plot_id, threshold_sites))
    }
    combs <- unique(combs)
    
    # randomly sample 100 combinations
    if (length(combs) > 100) {
      set.seed(19)
      combs <- sample(combs, size = 100)
    }
    
    # calculating the areas
    areas <- combs |>
      purrr::map_dbl(
        .f = function(comb) {
          temp_data |>
            dplyr::filter(Plot_id %in% comb) |>
            sf::st_union() |>
            sf::st_convex_hull() |>
            sf::st_cast("POLYGON") |>
            sf::st_area()
        }
      )
    # Drop the geometry column
    temp_data <- temp_data |>
      sf::st_drop_geometry()
    combs <-
      combs[areas |> dplyr::between(left = threshold_area * 0.5,
                                    right = threshold_area * 1.5)]
    for (comb in combs) {
      comm <- temp_data |>
        dplyr::filter(Plot_id %in% comb) |>
        dplyr::select(-(1:13),)
      # S = vegan::specnumber(comm)
      # comm <- comm[rowSums(comm) >= 5,] #if rowsums < 5, can not calculate betaC
      soma <- (colSums(comm))
      tsoma <- t(soma)
      gamma_N <- rowSums(tsoma)
      gamma_S <-
        sum(tsoma > 0) # sum only TRUE (which means occurence)
      gamma_Spie <- diversity(tsoma, index = "invsimpson")
      N.min <- as.numeric(N.min)
      gamma_Sn <- as.vector(rarefy(tsoma, N.min))
      alpha_S <-  rowSums(comm > 0)
      beta_S <- gamma_S / mean(alpha_S)
      alpha_Spie <- diversity(comm, index = "invsimpson")
      beta_Spie <- gamma_Spie / mean(alpha_Spie)
      alpha_N <- rowSums(comm)
      min_N <- min(c(alpha_N)) # minimum N in the pair of streams
      alpha_Sn <- rarefy(comm, sample = min_N)
      gamma_Sn.b <- rarefy(tsoma, sample = min_N)
      beta_Sn <- gamma_Sn.b / mean(alpha_Sn)
      #beta.c
      betac.calcu <-
        beta_C(comm,
               C.min,
               extrapolation = T,
               interrupt = T)
      bet.c <- c(betac.calcu)#get one vector
      #beta.spie Tory function
      beta_Spie.T <- beta_SPIE(comm)
      metrics <- c(
        gamma_N,
        gamma_S,
        gamma_Spie,
        gamma_Sn,
        beta_S,
        beta_Sn,
        beta_Spie,
        beta_Spie.T,
        bet.c
      )
      save.metrics <- rbind(na.omit(save.metrics), metrics)
    }
    metrics_ave <- t(as.data.frame(colMeans(save.metrics)))
    add_meta[[i]] <- cbind(meta_info, metrics_ave)
  }
  dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta)))
}
metrics_213 <- remove_rownames(bind_rows(c(dataset_output)))


######Manually corrected sites##########
all.files.manul <- list.files(path = "manually-gamma-txt",
                        pattern = "*.txt",
                        full.names = TRUE,
                        recursive = TRUE)
filenames.non.integer <- list.files(path = "manually-gamma-txt/non-integer",
                                    pattern = "*.txt")

for (i in 1:length(all.files.manul)) assign(
  x = basename(all.files.manul[[i]]), 
  value = read.table(file = all.files.manul[[i]],
                     header = TRUE,
                     sep = " ",
                     fill = TRUE,
                     row.names = NULL))
all.files.manul <- basename(all.files.manul)

for (i in 1:length(filenames.non.integer)) {
  print(filenames.non.integer[[i]])
  stu1 <- noquote(filenames.non.integer[[i]])
  stu2 <- get(stu1, 1)
  dat_abund.0 <- as.data.frame(stu2[ , -(1:14)])
  a <- ceiling(dat_abund.0)
  b <- cbind(stu2[ , (1:14)], a)
  assign(filenames.non.integer[i], b)
}

# Checking n of each study
n <- matrix(0, 27, 1)
for (x in 1:length(all.files.manul)) {
  stu1 <- noquote(all.files.manul[x])
  stu2 <- get(stu1, 1)
  n[x,] <- nrow(stu2)
}
rownames(n) <- all.files.manul
n

SC <- list()
for (m in 1:length(all.files.manul)) {
  print(all.files.manul[[m]])
  stu1 <- noquote(all.files.manul[[m]])
  Dataset <- get(stu1, 1)
  Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
  stu4 <- cbind.data.frame(Block, Dataset)
  spl.block <- split(stu4, Block) #spliting in datasets according to block
  blocos <- summary(as.factor(Block))
  n.blocks <- length(blocos)
  NS <- matrix(ncol = 1)
  colnames(NS) <- c("NS")
  C.coms <- matrix(ncol = 1)
  colnames(C.coms) <- c("C.coms")
  # start the loop between blocks for N.min and C.min
  for (t in 1:n.blocks) {
    #plots number counting
    Dataset <- spl.block[[t]] 
    for (i in unique(Dataset$Land_use)) {
      temp_data <- Dataset |>
        dplyr::filter(Land_use == i)
      meta_info <- temp_data |>
        slice(1) |>
        select(
          Dataset_id,
          Spatial_block,
          Temporal_block,
          Plot_id, Replicate_id,
          Taxa,
          Continent,
          Land_use
        ) 
        comm2 <- temp_data |> 
          dplyr::select(-(1:15), )
        soma2 <- (colSums(comm2))
        tsoma2 <- t(soma2)
        N2 <- rowSums(tsoma2)
        NS <- rbind(na.omit(NS), N2)
        C.comm2 <- C_target(comm2, factor = 2)
        C.coms <- rbind(na.omit(C.coms), C.comm2)
    }
  }
  NS.min <- t(as.data.frame(min(NS)))
  C.coms.min <- t(as.data.frame(min(C.coms)))
  Dataset_id <- distinct(Dataset,Dataset_id)
  SC[[m]] <- cbind.data.frame(Dataset_id,NS.min,C.coms.min)
}

SC <- data.frame(matrix(unlist(SC), nrow=27, byrow=TRUE),stringsAsFactors=FALSE)

colnames(SC) <-c("Dataset_id","N.min","Target.c")

#calculate the metrics
save_all_data <- list()
for (m in 1:length(all.files.manul)) {
  print(all.files.manul[[m]])
  stu1 <- noquote(all.files.manul[[m]])
  Dataset <- get(stu1, 1)
  Block <- paste(Dataset$Spatial_block, Dataset$Temporal_block, sep= "_")
  stu4 <- cbind.data.frame(Block, Dataset)
  spl.block <- split(stu4, Block) #spliting in datasets according to block
  blocos <- summary(as.factor(Block))
  n.blocks <- length(blocos)
  C.min <- SC[m, "Target.c"]
  N.min <-  SC[m, "N.min"]
  dataset_output <- list()
  for (t in 1:n.blocks) {
    Dataset <- spl.block[[t]]
    #minus one, so the land use with least sites number can be more flexiable
    add_meta <- list()
    for (i in unique(Dataset$Land_use)) {
      save.metrics <- matrix(ncol = 9)
      colnames(save.metrics) <- c("gamma_N", "gamma_S", "gamma_Spie","gamma_Sn",
                                  "beta_S","beta_Sn","beta_Spie","beta_Spie.T",
                                  "beta_C")
      temp_data <- Dataset |>
        dplyr::filter(Land_use == i)
      meta_info <- temp_data |>
        slice(1) |>
        select(
          Dataset_id,
          Block,
          Spatial_block,
          Temporal_block,
          Plot_id, Replicate_id,
          Taxa,
          Continent,
          Land_use
        ) 
      # All combinations of sites down to the number of sites the smallest lu has
        comm <- temp_data |>
          dplyr::select(-(1:15), )
        # S = vegan::specnumber(comm)
        # comm <- comm[rowSums(comm) >= 5,] #if rowsums < 5, can not calculate betaC
        soma <- (colSums(comm))
        tsoma <- t(soma)
        gamma_N <- rowSums(tsoma)
        gamma_S <-
          sum(tsoma > 0) # sum only TRUE (which means occurence)
        gamma_Spie <- diversity(tsoma, index = "invsimpson")
        N.min <- as.numeric(N.min)
        gamma_Sn <- as.vector(rarefy(tsoma, N.min))
        alpha_S <-  rowSums(comm > 0)
        beta_S <- gamma_S / mean(alpha_S)
        alpha_Spie <- diversity(comm, index = "invsimpson")
        beta_Spie <- gamma_Spie / mean(alpha_Spie)
        alpha_N <- rowSums(comm)
        min_N <- min(c(alpha_N)) # minimum N in the pair of streams
        alpha_Sn <- rarefy(comm, sample = min_N)
        gamma_Sn.b <- rarefy(tsoma, sample = min_N)
        beta_Sn <- gamma_Sn.b / mean(alpha_Sn)
        #beta.c
        betac.calcu <-
          beta_C(comm,
                 C.min,
                 extrapolation = T,
                 interrupt = T)
        bet.c <- c(betac.calcu)#get one vector
        #beta.spie Tory function
        beta_Spie.T <- beta_SPIE(comm)
        metrics <- c(
          gamma_N,
          gamma_S,
          gamma_Spie,
          gamma_Sn,
          beta_S,
          beta_Sn,
          beta_Spie,
          beta_Spie.T,
          bet.c
        )
      save.metrics <- rbind(na.omit(save.metrics), metrics)
      add_meta[[i]] <- cbind(meta_info,save.metrics)
    }
    dataset_output[[t]] <- remove_rownames(bind_rows(c(add_meta)))
  }
  save_all_data[[m]] <- remove_rownames(bind_rows(c(dataset_output)))
}
metrics_manual <- remove_rownames(bind_rows(c(save_all_data)))


##########Combined all metrics results#########
gamma_bet_metrics <- bind_rows(metrics_loop,metrics_group,metrics_manual,
                                metrics_9,metrics_33,metrics_38,metrics_213)

gamma_bet_metrics <- gamma_bet_metrics |> 
  drop_na(gamma_N)

gamma_bet_metrics |> 
  group_by(Dataset_id) |>
  summarise(Count = n()) |> 
  filter(Count==1)
x<- gamma_bet_metrics |> 
  group_by(Dataset_id,Land_use) |>
  summarise(Count = n())
x |> group_by(Dataset_id) |>
  summarise(Count = n()) |> 
  filter(Count==1)

gamma_bet_metrics <- gamma_bet_metrics |> 
  filter(!Dataset_id %in% c("17","23","31","41","86","167","177"))
gamma_bet_metrics <-  bind_rows(gamma_bet_metrics,metrics_17)
write.table(gamma_bet_metrics,"gamma_bet_metrics.txt")

