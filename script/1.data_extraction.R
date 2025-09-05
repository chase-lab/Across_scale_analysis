setwd("C:/Users/fx44ojik.USR/Dropbox/iDiv project/About dataset/Data analysis/20231130_analysis/database")

plot_data<-readxl::read_xlsx(path = "Land-use-and-freshwater-database.xlsx", sheet = 2)

observations<-readxl::read_xlsx(path = "Land-use-and-freshwater-database.xlsx", sheet = 3)

landuse<-readxl::read_xlsx(path = "land_use_category.xlsx", sheet = 1)

landuse<-landuse |> 
  dplyr::select(Dataset_id,Site_id,Plot_id,Land_use)

lu_included <- filter(landuse, Land_use %in% c("Natural vegetation","Forestry", "Agriculture",
                                               "Urban", "Mining"))

plot_data <- dplyr::select(plot_data,Plot_id,Site_id,Dataset_id,Stream_order,Latitude,Longitude,Continent)

dat_plots<-left_join(lu_included,plot_data,by=c("Plot_id","Site_id","Dataset_id"))

observations <- dplyr::select(observations,Sample_id,Plot_id,Dataset_id,Spatial_block,Temporal_block,
                              Replicate_id,Start_year,Sample_effort,Sample_area,Taxon,Abundance,
                              Abundance_converted,Abundance_Type)

all_dat<-left_join(dat_plots,observations,by=c("Plot_id","Dataset_id"))

str(all_dat)

all_dat$Land_use <-gsub('Natural vegetation','Natural_vegetation',all_dat$Land_use)

all_dat$Stream_order<-gsub('/','NA',all_dat$Stream_order)

all_dat <- all_dat |> 
  mutate(Abundance = Abundance) |> 
  mutate(Abundance = ifelse(Abundance_converted !="NA",Abundance_converted, Abundance))

setwd("txt")

all_dat$Abundance <- as.numeric(all_dat$Abundance)
all_dat$Dataset_id <- as.character(all_dat$Dataset_id)

all_dat  <- all_dat |> 
  mutate(integercheck = Abundance %% 1 == 0)

noninteger <- all_dat |>  filter(integercheck == FALSE ) |> 
  distinct(Dataset_id)

yesinteger <-  all_dat |> 
  distinct(Dataset_id) |> 
  filter(!Dataset_id %in% noninteger$Dataset_id)

for (n in unique(yesinteger$Dataset_id)) {
  Dataset.n<-all_dat |>
    filter(Dataset_id==n)
  Dataset.n<-Dataset.n |>
    as_tibble() |>
    group_by(Sample_id,Plot_id, Site_id,Dataset_id,Replicate_id, Spatial_block,Temporal_block,Sample_effort,Latitude,Longitude,Taxa,Stream_order, Land_use,Continent,Taxon) |>
    dplyr::summarise(Abundance = sum(Abundance)) |> 
    ungroup() |>
    tidyr::pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)
  form<-sprintf('Dataset_%s.txt', n)
  write.table(Dataset.n,file=form,sep=" ",row.names = FALSE,quote =FALSE,fileEncoding = "UTF-8" )
}

setwd("non-integer")

for (n in unique(noninteger$Dataset_id)) {
  Dataset.n<-all_dat |>
    filter(Dataset_id==n)
  Dataset.n<-Dataset.n |>
    as_tibble() |>
    group_by(Sample_id,Plot_id, Site_id,Dataset_id,Replicate_id, Spatial_block,Temporal_block,Sample_effort,Latitude,Longitude,Taxa,Stream_order, Land_use,Continent,Taxon) |>
    dplyr::summarise(Abundance = sum(Abundance)) |> 
    ungroup() |>
    tidyr::pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)
  form<-sprintf('Dataset_%s.txt', n)
  write.table(Dataset.n,file=form,sep=" ",row.names = FALSE,quote =FALSE,fileEncoding = "UTF-8" )
}
