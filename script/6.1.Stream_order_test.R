library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)

# column for stream order #ORD_STRA 
#for alpha S: S_alpha_alt_3
#for alpha Spie: Spie_alpha_alt_2
getwd()

read.table("~/Models-alpha_beta_gamma/mob_data_alpha_so.txt", 
           quote="\"", comment.char="")

load("Spie_alpha_alt_2.Rdata")

load("S_alpha_alt_3.Rdata")

## alpha S#####
resid<-residuals(S_alpha_alt_3,
                 type = 'pearson',
                 method = 'predict') %>% 
  dplyr::as_tibble() %>% bind_cols(S_alpha_alt_3$data) 

fitted <- fitted(S_alpha_alt_3, re_formula = NA)
predict <- predict(S_alpha_alt_3)

resid$fitted <- fitted[,'Estimate']
resid$predict <- predict[,'Estimate']

plot(resid$Estimate ~ resid$fitted,ylab = 'Pearson residual')

plot(resid$Estimate ~ resid$Land_use,ylab = 'Pearson residual',xlab="Land_use")

#draw study level estimates
Forest <- S_alpha_alt_3 %>%
  spread_draws(
    b_Land_useForest,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111) 

Forest <- Forest %>%
  mutate(Forest=b_Land_useForest+r_Dataset_id) %>%
  filter(Land_use=="Land_useForest")

Agriculture <- S_alpha_alt_3 %>% 
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111) 

Agriculture <- Agriculture %>% 
  mutate(Agriculture=b_Land_useAgriculture+r_Dataset_id) %>% 
  filter(Land_use=="Land_useAgriculture")

Urban <- S_alpha_alt_3 %>% 
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111)
Urban <- Urban %>% 
  mutate(Urban=b_Land_useUrban+r_Dataset_id) %>% 
  filter(Land_use=="Land_useUrban")

Forestry <- S_alpha_alt_3 %>% 
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111)
Forestry <- Forestry %>% 
  mutate(Forestry=b_Land_useForestry+r_Dataset_id) %>% 
  filter(Land_use=="Land_useForestry")

# Mining <- S_alpha_alt_3 %>% 
#   spread_draws(
#     b_Land_useMining,
#     r_Dataset_id[Dataset_id,Land_use],
#     ndraws = 1000,seed = 111)
# Mining <- Mining %>% 
#   mutate(Mining=b_Land_useMining+r_Dataset_id) %>% 
#   filter(Land_use=="Land_useMining")

Forest$AN<-(Agriculture$Agriculture)-
  (Forest$Forest)
Forest$UN<-(Urban$Urban)-
  (Forest$Forest)
Forest$FN<-(Forestry$Forestry)-
  (Forest$Forest)
# Forest$MN<-(Mining$Mining)-
#   (Forest$Forest)
# Forest$UA<-(Urban$Urban)-
#   (Agriculture$Agriculture)

Comparisons <- Forest %>% 
  pivot_longer(cols = c("AN","UN","FN"),
               names_to = "Comparison",
               values_to = "Ratio")

Comparisons$Comparison<- factor(Comparisons$Comparison,
                                levels = c("AN","UN","FN"),
                                labels = c('Agriculture/Natural vegetation',
                                           'Urban/Natural vegetation',
                                           'Forestry/Natural vegetation'
                                ))
mob_data_alpha_so <- mob_data_alpha_so %>% 
  select(c("Dataset_id","ORD_STRA"))

Comparisons <- left_join(Comparisons,mob_data_alpha_so,by='Dataset_id')

ComparisonsA0 <- Comparisons

#stream order

Comparisons$ORD_STRA <- as.character(Comparisons$ORD_STRA)


ggplot() +
  geom_density_ridges_gradient(data = Comparisons,
                               aes(x = Ratio,
                                   y = ORD_STRA))

alpha_stream_order_S <- ggplot() +
  geom_density_ridges_gradient(data = Comparisons,
                               aes(x = Ratio,
                                   y = ORD_STRA
                                   # ,
                                   # fill = factor(after_stat(x) > 0)
                                   # ,
                                   # fill = stat(quantile)
                               ),
                               # quantiles = c(0.025, 0.25, 0.75, 0.975),
                               # calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0)+
  xlim(-1,1)+
  geom_vline(data = Comparisons,
             aes(xintercept = mean(Ratio)),
             size = 0.5,
             alpha = 0.5,
             lty = 2) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(data = Comparisons,
             aes(x = Ratio,
                 y = ORD_STRA),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18, colour = 'black')+
  labs(y = 'S',
       x = 'Log Ratio'
  ) +
  geom_text(
    data = Comparisons %>%
      filter(Ratio < 0) %>%
      group_by(ORD_STRA) %>%
      summarise(Count = n()) %>%
      mutate(percentages = Count / 1000) %>%
      ungroup() %>%
      distinct(ORD_STRA, percentages, .keep_all = T),
    aes(
      x = -0.85,
      y = ORD_STRA,
      label = paste(percentages)
    ),
    size = 3,
    nudge_y = 0.5,
    parse = T
  )+ theme(legend.position="none")+
  # + guides(fill=guide_legend(title=" "))+
  scale_fill_manual(values = c("black","gray"),labels =c("Below 0","Above 0"))+
  scale_y_discrete(labels = scales::wrap_format(9))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 7))

###### alpha Spie#####

Forest <- Spie_alpha_alt_2 %>%
  spread_draws(
    b_Land_useForest,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111) 

Forest <- Forest %>%
  mutate(Forest=b_Land_useForest+r_Dataset_id) %>%
  filter(Land_use=="Land_useForest")

Agriculture <- Spie_alpha_alt_2 %>% 
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111) %>% 
  print(n=100)
Agriculture <- Agriculture %>% 
  mutate(Agriculture=b_Land_useAgriculture+r_Dataset_id) %>% 
  filter(Land_use=="Land_useAgriculture")

Urban <- Spie_alpha_alt_2 %>% 
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111)
Urban <- Urban %>% 
  mutate(Urban=b_Land_useUrban+r_Dataset_id) %>% 
  filter(Land_use=="Land_useUrban")

Forestry <- Spie_alpha_alt_2 %>% 
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id,Land_use],ndraws = 1000,seed = 111)
Forestry <- Forestry %>% 
  mutate(Forestry=b_Land_useForestry+r_Dataset_id) %>% 
  filter(Land_use=="Land_useForestry")

Forest$AN<-(Agriculture$Agriculture)-
  (Forest$Forest)
Forest$UN<-(Urban$Urban)-
  (Forest$Forest)
Forest$FN<-(Forestry$Forestry)-
  (Forest$Forest)


Comparisons <- Forest %>% 
  pivot_longer(cols = c("AN","UN","FN"),
               names_to = "Comparison",
               values_to = "Ratio")

Comparisons$Comparison<- factor(Comparisons$Comparison,
                                levels = c("AN","UN","FN"),
                                labels = c('Agriculture/Forest',
                                           'Urban/Forest',
                                           'Forestry/Forest'
                                ))

Comparisons <- left_join(Comparisons, mob_data_alpha_so,by='Dataset_id')

ComparisonsA2 <- Comparisons

###plot S and Spie in different stream order

Comparisons0 <- ComparisonsA0 %>% 
  rename(RatioA0 = Ratio) 

Comparisons0 <- Comparisons0 %>% 
  ungroup() %>% 
  select(c("ORD_STRA","RatioA0"))

Comparisons2 <- ComparisonsA2 %>% 
  rename(RatioA2 = Ratio)

Comparisons2 <- Comparisons2 %>% 
  ungroup() %>% 
  select(c("ORD_STRA","RatioA2"))

Comparisons02 <-  cbind(Comparisons0,Comparisons2)

Comparisons02 <- Comparisons02 %>% select(-3,)

Comparisons02 <- Comparisons02 %>%
  pivot_longer(cols = starts_with("ratio"), 
               names_to = "ratio_type", 
               values_to = "Ratio")

Comparisons02 %>% 
  group_by(ORD_STRA) %>% 
  summarise(n=n())

 

Comparisons02$ORD_STRA <- as.factor(Comparisons02$ORD_STRA)

n_total <- Comparisons02 %>%  
  group_by(ORD_STRA,ratio_type) %>% 
  summarise(n_total=n()) %>% 
  ungroup() 

n_count <- Comparisons02 %>%
  filter(Ratio < 0) %>%
  group_by(ORD_STRA,ratio_type) %>%
  summarise(n_count = n())

n_total_count <- left_join(n_total,n_count,by=c("ORD_STRA","ratio_type"))

#plot
alpha_stream_order <- ggplot() +
  geom_density_ridges_gradient(data = Comparisons02,
                               aes(x = Ratio,
                                   # y = Comparison,
                                   y = ORD_STRA,
                                   fill = ratio_type
                                   # y = fct_reorder(ORD_STRA, Ratio, .fun = mean),
                                   # fill = factor(after_stat(x) > 0)
                               ),scale = 0.9, alpha = 0.1,
                               linetype = 0)+
  geom_point(data = Comparisons02,
             aes(x = Ratio,
                 y = ORD_STRA,
                 color = ratio_type),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             shape = 21,  fill = "white",
             stroke = 2,
             size = 4.5)+
  xlim(-1.0,1.0)+
  geom_text(
    data = n_total_count  %>%
      mutate(d = n_count / n_total) %>%
      mutate(percentages = sprintf('%.2f',d))   %>%
      distinct(ORD_STRA, ratio_type, percentages, .keep_all = T),
    aes(x=-0.85,
        # y = ORD_STRA,
        y = as.numeric(ORD_STRA) + ifelse(ratio_type == "RatioA0", 0.55, 0.25),
        color = ratio_type,
        label = percentages),
    size = 4.5, parse = T,show.legend = FALSE
  )+
  geom_vline(xintercept = 0, size = 0.5,linetype = 2)+
  labs(y = 'Basin size',
       x = 'Effect size') +
  # guides(fill=guide_legend(title=" "))+
  scale_fill_manual(values = c("RatioA0"=scales::alpha("#F2790F",0.7),
                               "RatioA2"=scales::alpha( "#CAAE10",0.7)),
                    name = " "
                    ,
                    # labels = c("RatioA0" = " ", 
                    #            "RatioA2" = " ")
                    labels = c("RatioA0" = "S",
                               "RatioA2" = "Spie")
  )+
  scale_color_manual(values = c("RatioA0"=scales::alpha("#F2790F",0.9),
                                "RatioA2"=scales::alpha( "#CAAE10",0.9)),
                     name = " "
                     ,
                     # labels = c("RatioA0" = " ", 
                     #            "RatioA2" = " ")
                     labels = c("RatioA0" = "S",
                                "RatioA2" = "Spie")
  )+
  scale_y_discrete(labels = scales::wrap_format(9))+
  theme_minimal() +
  guides(color = guide_legend(byrow = FALSE))+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.key.spacing.x = unit(1, "cm"))+
  # theme(legend.position = "none")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5, color = "gray"),
        panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        legend.text = element_text(size = 16))

ggsave("alpha_stream_order.pdf",alpha_stream_order)
