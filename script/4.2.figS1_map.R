
world <- map_data('world')

Metadata <- read_csv("data/Metadata.csv")

#choose alpha dataset_id, code at  3.1alpha_model.R

metadata_map <- Metadata |> 
  filter(Dataset_id %in% alpha_dataset_id)

str(metadata_map)

taxa_color= c('Macroinvertebrate' = '#B93102',
                 'Fish' = '#274659',
                 'Zooplankton' = '#CAAE10',
                 'Algae' = '#fb9a99',
                 'Amphibian' = '#F2790F')

ggplot() +
  geom_polygon(data=world, 
               aes(long, lat, group = group), colour=NA, fill='#f0f0f0', size=0) +
  geom_point(data = metadata_map,
             aes(x = Longitude, y = Latitude, colour = Taxa),
             size = 3,alpha = 0.5
  ) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(name = 'Longitude', breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(name = 'Latitude', breaks = c(0, -23.5, 23.5, -60, 60)) +
  scale_size_area(name = 'Number of locations',
                  trans = 'log10') +# name = 'Duration (years)', breaks = c(10,40,160)
  scale_colour_manual(name = 'Biological group', values = taxa_color) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = 'black', size = 0.1), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top',
        plot.margin = unit(c(0,0,0,0), units = 'mm'),
        legend.margin = margin(),
        legend.box.spacing = unit(c(0,0,0,0), units = 'mm'),
        legend.text = element_text(size = 10, face = 'plain'),
        legend.title = element_text(size = 11, face = 'bold')) +
  guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5),
         shape = guide_legend(title.position = 'top', title.hjust = 0.5, size = 3),
         size = guide_legend(title.position = 'top', title.hjust = 0.5))

ggsave('figures/FigS1_map.pdf',
       width = 184, height = 150, units = 'mm')
