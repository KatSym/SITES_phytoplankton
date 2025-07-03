library(tidyverse)
library(vegan)

# load data
load("./data/sites_FC_phyto.RData")


div.e <- erk.dat %>% 
  # drop the lake
  filter(Treatment != "ERK",
         # remove classes that don't represent one taxonomic group or are zooplankton
         !label %in% c("Amoebida_761", "Cil_other", "Ciliate_0000000", "Multispecies", "OTHER", 
                       "ROTIFERA", "FILAMENTS", "Chlcolony_0000000", "Heliozoa_0000000",
                       "ProCol000_1001946")) %>% 
  mutate(Treatment = droplevels(Treatment)) %>%
  group_by(ExpDay, Treatment, mesocosm, Name, label) %>%
  summarise(count = n(),
            mean_abd = mean(AbdDiameter),
            mean_len = mean(Length),
            mean_biovol = mean(biovolMS),
            mean_sa = mean(surfacearea),
            mean_prob = mean(ProbabilityScore),
            volume_run = mean(volume_run),
            concfactor = mean(concfactor)
  ) %>% 
  mutate(abund = 1000*(count/volume_run)/concfactor, # in individuals/L
         biovolume = mean_biovol * abund) %>% # um3 / L 

select(ExpDay, Treatment, mesocosm, label, abund, biovolume)  

biov.e <- div.e %>%   
  group_by(ExpDay, Treatment, mesocosm, label) %>% 
  summarise(biov = mean(biovolume)) %>% 
  pivot_wider(id_cols = c("ExpDay", "Treatment", "mesocosm"), 
              names_from = label,
              values_from = biov
  ) %>% 
  replace(is.na(.), 0)

abun.e <- div.e %>%   
  group_by(ExpDay, Treatment, mesocosm, label) %>% 
  summarise(abun = mean(abund)) %>% 
  pivot_wider(id_cols = c("ExpDay", "Treatment", "mesocosm"), 
              names_from = label,
              values_from = abun
  ) %>% 
  replace(is.na(.), 0)

He.biov <- diversity(biov.e[4:34]) 
He.abun<- diversity(abun.e[4:34]) 

div.erk <- div.e %>% 
  group_by(Name, ExpDay, Treatment, mesocosm, label) %>% 
  summarise(sp.rich = n()) %>% 
  group_by(ExpDay, Treatment, mesocosm) %>% 
  summarise(sp.rich = sum(sp.rich)) %>% 
  ungroup() %>% 
  mutate(Hbv = He.biov,
         Hab = He.abun,
         Jab = Hab/log(sp.rich),
         Jbv = Hbv/log(sp.rich)) %>% 
  
  pivot_longer(cols = c("Hbv", "Hab", "Jab", "Jbv", "sp.rich"), 
              names_to = "index",
              values_to = "vals"
  )

div.erk %>% 
  ggplot(., aes(x = ExpDay, y = vals,  colour = Treatment))+
  geom_point(size=0.8, alpha=0.6, 
             position = position_jitterdodge()) +
  stat_summary(
    aes(group = Treatment, colour = Treatment),
    geom = "line",
    size = .7,
    fun.y = "mean",
  )+
  facet_wrap(~index, scales = "free_y")+
  # scale_fill_manual(values = trt.cols)+
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme")) +
  theme_bw() +
  labs(
    # title = expression(.$LabelTrue),
    y = "Abundance /L",
    x = "Date")

