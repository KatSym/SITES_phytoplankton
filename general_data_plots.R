library(tidyverse)
library(ggpubr)

# load data
load("./data/sites_FC_phyto.RData")


erken <- erk.dat %>% 
  # drop the lake
  filter(Treatment != "ERK",
         ProbabilityScore > 0.5,
         # remove classes that don't represent one taxonomic group or are zooplankton
         !label %in% c("Amoebida_761", "Cil_other", "Ciliate_0000000", "Multispecies", "OTHER", 
                       "ROTIFERA", 
                       # "FILAMENTS", 
                       # "Chlcolony_0000000", 
                       "Heliozoa_0000000", "ProCol000_1001946")) %>% 
  mutate(Treatment = droplevels(Treatment)) %>%
  group_by(ExpDay, Treatment, mesocosm, label) %>%
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
  group_by(ExpDay, Treatment, mesocosm) %>%
  summarise(abundance = sum(abund),
            biovol = sum(biovolume)) %>% 
  pivot_longer(cols = c("abundance", "biovol"), 
               names_to = "metric",
               values_to = "vals")





erk <- erken %>% 
  ggplot(., aes(x = ExpDay, y = vals,  colour = Treatment))+
  geom_point(size=0.8, alpha=0.6, 
             position = position_jitterdodge()) +
  stat_summary(
    aes(group = Treatment, colour = Treatment),
    geom = "line",
    size = .7,
    fun.y = "mean",
  )+
  facet_wrap(~metric, scales = "free_y")+
  # scale_fill_manual(values = trt.cols)+
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme")) +
  theme_bw() +
  labs(
    title = "Erken",
    y = "Abundance /L",
    x = "Experimental day")


# bolmen -----


bolmen <- bol.dat %>% 
  # drop the lake
  filter(Treatment != "BOL",
         ProbabilityScore > 0.5,
         # remove classes that don't represent one taxonomic group or are zooplankton
         !label %in% c("Amoebida_761", "Ciliate_0000000", "Multispecies", "OTHER", 
                       "ROTIFERA", 
                       # "FILAMENTS", 
                       # "Chlcolony_0000000", 
                       "Heliozoa_0000000", "ProCol000_1001946", "OliVor000_7351754")) %>% 
  mutate(Treatment = droplevels(Treatment)) %>%
  group_by(ExpDay, Treatment, mesocosm, label) %>%
  summarise(count = n(),
            mean_abd = mean(AbdDiameter),
            mean_len = mean(Length),
            mean_biovol = mean(biovolume),
            mean_sa = mean(surfacearea),
            mean_prob = mean(ProbabilityScore),
            volume_run = mean(volume_run),
            concfactor = mean(concfactor)
  ) %>% 
    mutate(abund = 1000*(count/volume_run)/concfactor, # in individuals/L
         biovolume = mean_biovol * abund) %>% # um3 / L 
  group_by(ExpDay, Treatment, mesocosm) %>%
  summarise(abundance = sum(abund),
            biovol = sum(biovolume)) %>% 
  pivot_longer(cols = c("abundance", "biovol"), 
               names_to = "metric",
               values_to = "vals")





bol <- bolmen %>% 
  ggplot(., aes(x = ExpDay, y = vals,  colour = Treatment))+
  geom_point(size=0.8, alpha=0.6, 
             position = position_jitterdodge()) +
  stat_summary(
    aes(group = Treatment, colour = Treatment),
    geom = "line",
    size = .7,
    fun.y = "mean",
  )+
  ggh4x::facet_grid2(~metric, scales = "free_y", independent = "y") +  
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme")) +
  theme_bw() +
  labs(
    title = "Bolmen",
    y = "Abundance /L",
    x = "Experimental day")

b1 <- bol + annotate("rect",
               xmin = 19, xmax = 21, 
               ymin = -Inf, ymax = Inf,  fill = "red", alpha=.1)

e1 <- erk + annotate("rect",
               xmin = 35, xmax = 37, 
               ymin = -Inf, ymax = Inf,  fill = "red", alpha=.1)

ggarrange(e1 + rremove("xlab"),
          b1, nrow = 2)



## Functional groups

# from the models.R script
edat.f <- edat %>% 
  mutate(biovolume = gr_biovol/vol.offset,
         abundance = count/vol.offset, 
         .keep = "unused") %>% 
  pivot_longer(cols = c("abundance", "biovolume"), 
               names_to = "metric",
               values_to = "vals") %>% 
  ungroup() %>% 
  mutate(fun_grp = as.factor(fun_grp))

erk.f <- edat.f %>% 
  ggplot(., aes(x = ExpDay, y = vals,  colour = Treatment))+
  geom_point(size=0.8, alpha=0.6, 
             position = position_jitterdodge()) +
  stat_summary(
    aes(group = Treatment, colour = Treatment),
    geom = "line",
    size = .7,
    fun.y = "mean",
  )+
  ggh4x::facet_grid2(metric~fun_grp, scales = "free_y", independent = "y",
             # space = "free", 
             # axes = "all_y"
             )+
  # scale_fill_manual(values = trt.cols)+
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme")) +
  theme_bw() +
  labs(
    title = "Erken",
    y = "um3 / L                                                 ind./L",
    x = "Experimental day")



bolmen.f <- bol.dat %>% 
  # drop the lake
  filter(Treatment != "BOL",
         ProbabilityScore > 0.5,
         # remove classes that don't represent one taxonomic group or are zooplankton
         !label %in% c("Amoebida_761", "Ciliate_0000000", "Multispecies", "OTHER", 
                       "ROTIFERA", 
                       # "FILAMENTS", 
                       # "Chlcolony_0000000", 
                       "Heliozoa_0000000", "ProCol000_1001946", "OliVor000_7351754")) %>% 
  mutate(vol.offset = volume_run*concfactor, .keep = "unused") %>% 
  group_by(ExpDay, Treatment, mesocosm, Name, label) %>%
  summarise(count = n(),
            mean_abd = mean(AbdDiameter),
            mean_len = mean(Length),
            mean_biovol = mean(biovolume),
            mean_sa = mean(surfacearea),
            mean_prob = mean(ProbabilityScore),
            vol.offset = mean(vol.offset)
  ) %>% 
  mutate(Treatment = droplevels(Treatment),
         class = substr(label, 1, 3),
         fun_grp = case_when(class == "Bac" ~ "VI",
                             label == "ChlAst000_2639495" ~ "I",
                             label == "Chlcolony_0000000" ~ "I",
                             label == "ChlDes000_2652440" ~ "IV",
                             label == "ChlEud000_0000000" ~ "VII",
                             label == "ChlPed000_2641883" ~ "IV",
                             label == "ChlSel000_0000000" ~ "IV",
                             label == "ChrDin000_3194946" ~ "II",
                             label == "ChrMal000_3195066" ~ "II",
                             label == "ChrSyn000_0000000" ~ "II",
                             label == "ConSta000_2647648" ~ "IV",
                             label == "CyaAphflo_7690242" ~ "III",
                             label == "CyaDollem_7901634" ~ "III",
                             label == "CyaMic000_3217749" ~ "VII",
                             label == "CyaPla000_3218374" ~ "III",
                             label == "CyaWornae_3217010" ~ "VII",
                             label == "DinCerhir_7598904" ~ "V",
                             label == "FILAMENTS" ~ "III",
                             label == "RapGonsem_3194431" ~ "IV",
                             label == "ZygClo000_2646356" ~ "IV",
                             label == "ZygCloacu_2646408" ~ "IV",
                             label == "ZygCos000_2648709" ~ "IV"),
         .keep = "unused") %>% 
  drop_na() %>% 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) %>% 
  summarise(count = sum(count),
            biovol = sum(mean_biovol),
            vol.offset = mean(vol.offset)) %>% 
   mutate(biovolume = count*biovol/vol.offset,
         abundance = count/vol.offset, 
         .keep = "unused") %>% 
  pivot_longer(cols = c("abundance", "biovolume"), 
               names_to = "metric",
               values_to = "vals") %>% 
  ungroup() %>% 
  mutate(fun_grp = as.factor(fun_grp))


bol.f <- bolmen.f %>% 
  ggplot(., aes(x = ExpDay, y = vals,  colour = Treatment))+
  geom_point(size=0.8, alpha=0.6, 
             position = position_jitterdodge()) +
  stat_summary(
    aes(group = Treatment, colour = Treatment),
    geom = "line",
    size = .7,
    fun.y = "mean",
  )+
  ggh4x::facet_grid2(metric~fun_grp, scales = "free_y", independent = "y",
                     # space = "free", 
                     # axes = "all_y"
  )+
  # scale_fill_manual(values = trt.cols)+
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme")) +
  theme_bw() +
  labs(
    title = "Bolmen",
    y = "um3 / L                                                 ind./L",
    x = "Experimental day")
