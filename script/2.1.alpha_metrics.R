#******* Code to obtain MoB metrics at alpha scale using spatial extent procedure ********#

#------------------ 1. Data handling ------------------------#
library(tidyverse)
library(brms)
library(devtools)
library(vegan)
library(SpadeR)
library(plyr)

setwd("~/share/groups/synthesis/Minghua/Across_scale_analysis/data")

all.files <- list.files(path = "alpha_txt",
                        pattern = "*.txt",
                        full.names = TRUE,
                        recursive = TRUE)
filenames.non.integer <- list.files(path = "txt/non-integer",
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
n <- matrix(0,207,1)
for (x in 1:length(all.files)){
  stu1 <-noquote(all.files[x])
  stu2 <- get(stu1,1)
  n[x,] <- nrow(stu2)}
rownames(n) <- all.files; n


# for (i in 1:length(all.files.alpha)){
#   stu1 <-noquote(all.files.alpha[i])
#   stu2 <- get(stu1,1)
#   if (length (unique(stu2$Sample_effort))!=1 )
#   {print(all.files.alpha[i])}
#   }

##------------- 2. Calculating MoB statistics at alpha-scale ---------------#
div_list <- list()
all.files.alpha<-all.files
for (i in 1:length(all.files.alpha)){
  print(all.files.alpha[i])
  stu1 <-noquote(all.files.alpha[i])
  stu2 <- get(stu1,1)
  colnames(stu2)[13] <- "Land_use"
  stu2 <- stu2 %>% filter(Sample_id !="NA")
  # name.i <- unlist(strsplit(stu1, split='.', fixed=TRUE))[1]
  # Study <- rep(name.i, times=nrow(stu2))
  # taxa.i <- metadata.alpha[i, "Taxa"]
  # Taxa <- rep(taxa.i, times=nrow(stu2))
  # stu3 <- cbind(Study,Taxa, stu2)
  dat_abund.0 <- as.data.frame(stu2 [ ,-(1:14)])
  dat_abund <- dat_abund.0 [rowSums(dat_abund.0)>0, ]
  stu4 <- stu2 [rowSums(dat_abund.0)>0, ]
  Block <- paste(stu4$Spatial_block, stu4$Temporal_block, sep= "_")
  stu5 <- cbind.data.frame(Block, stu4)
  if (length (unique(stu5$Sample_effort))==1 ) {
    N <- rowSums(dat_abund)
    S <- rowSums(dat_abund > 0)
    Spie <- diversity(dat_abund, index = "invsimpson")
    minimo <- min(N)
    Sn <- rarefy(dat_abund, sample= minimo) 
    indices <- cbind(N, S, Spie, Sn)
    info <- stu5 [ ,c("Taxa", "Land_use", "Latitude","Longitude",
                      "Dataset_id","Site_id", "Plot_id", "Sample_id","Replicate_id", 
                      "Block", "Spatial_block", "Temporal_block","Continent")]
    div_indi <- cbind(info, indices)
    div_list[[all.files.alpha[i]]] <- div_indi
  } else {
    stu5 <- stu5 %>% dplyr::mutate(row_id = row_number(),.before = "Block")
    div_indi <- list()
    for (m in unique(stu5$row_id)) {
      minsameff.m <- min(stu5$Sample_effort,na.rm=T)
      row_id.m<-stu5 %>%
          filter(row_id==m)
      indiave <- data.frame()
        if(!is.na(row_id.m$Sample_effort)){
            taxondata <- as.data.frame(row_id.m[ ,-(1:16)])#check the number
            taxondata <- pivot_longer(taxondata, cols = tidyselect::everything(), 
                                      names_to = "Taxon", values_to = "Abundance")
            pooled<- rep(taxondata$Taxon,taxondata$Abundance)
            npooled<-length(pooled)
            sameff<-row_id.m$Sample_effort
            incidescom<-data.frame()
            for (n in 1:100) {
              resample<-replicate(1, sample(pooled,size=ceiling((minsameff.m/sameff)*npooled),replace=TRUE,prob = NULL),)
              resample.0<-pivot_wider(as.data.frame(table(resample)),names_from = resample, values_from = Freq)
              resample.abun<-resample.0 [rowSums(resample.0)>0, ]
              N<-rowSums(resample.abun) ##it is equal to size
              S<-rowSums(resample.abun>0)
              Spie<-diversity(resample.abun, index = "invsimpson") 
              minimo <- min(N)
              Sn <- rarefy(resample.abun, sample= minimo)
              indices <- cbind(N, S, Spie, Sn)
              incidescom<-rbind(incidescom,indices)
            }
            indiave<-rbind(indiave,colMeans(incidescom))
        }
        colnames(indiave) <- c('N','S','Spie','Sn')
        info <- row_id.m [ ,c("Taxa", "Land_use", "Latitude","Longitude",
                          "Dataset_id","Site_id", "Plot_id", "Sample_id","Replicate_id", 
                          "Block", "Spatial_block", "Temporal_block","Continent")]
        div_indi[[m]] <- cbind(info, indices)
        
    }  
    div_list[[all.files.alpha[i]]] <- remove_rownames(bind_rows(c(div_indi)))
  }
}

div_list

all.stu <- ldply(div_list, data.frame) #extracting from list to data.frame
mob_data_alpha <- all.stu 
dim(mob_data_alpha)

setwd("~/share/groups/synthesis/Minghua/Across_scale_analysis/results/temp")

write.table(mob_data_alpha,"alpha_metrics.txt")

