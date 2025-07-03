library(tidyverse)
library(brms)
library(tidybayes)

# load data
load("./data/sites_FC_phyto.RData")
functional.erk <- read.csv("./data/Erken_functional.csv", header = T, sep = ",")

edat <- erk.dat %>% 
  # drop the lake
  filter(Treatment != "ERK") %>% 
  mutate(Treatment = droplevels(Treatment),
         vol.offset = volume_run*concfactor,
         .keep = "unused") %>%
  group_by(ExpDay, Treatment, mesocosm, Name, label) %>%
  summarise(count = n(),
            mean_abd = mean(AbdDiameter),
            mean_len = mean(Length),
            mean_biovol = mean(biovolMS),
            mean_sa = mean(surfacearea),
            mean_prob = mean(ProbabilityScore),
            vol.offset = mean(vol.offset)
            ) %>% 
  
  left_join(., functional.erk, by = join_by(label == Code)) %>%
  mutate(class = substr(label, 1, 3),
         fun_grp = case_when(label == "FILAMENTS" ~ "III",
                             label == "CyaPla000_3218374" ~ "III",
                             class == "Bac" ~ "VI",
                             .default = KRUK_MBFG)) %>% 
  select(-taxon, -KRUK_MBFG, -Name) %>% 
  
  drop_na() %>% 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) %>% 
  summarise(count = sum(count),
            vol.offset = mean(vol.offset)) 



 

ptm <- proc.time()
m1.e <- brm(bf(count ~ Treatment * ExpDay +
                 (ExpDay|mesocosm) +
                 offset(vol.offset)),
            family = "negbinomial",
            chains = 4,
            iter = 4000,
            cores = 4,
            # control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = edat,
            file = "models/250702_Erk-count-offset.m",
            file_refit = "on_change"
            ) 
proc.time() - ptm
summary(m1.e)

plot(conditional_effects(m1.e, effects = "ExpDay:Treatment"), points = TRUE)

ptm <- proc.time()
m2.e <- brm(bf(count ~ Treatment * ExpDay + fun_grp * Treatment +
                 (ExpDay|mesocosm) +
                 offset(vol.offset)),
            family = "negbinomial",
            chains = 4,
            iter = 4000,
            cores = 4,
            # control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = edat,
            file = "models/250702_Erk_FunctGrp-count-offset.m",
            file_refit = "on_change"
) 
proc.time() - ptm

summary(m2.e)
plot(conditional_effects(m2.e, offset = T), ask = FALSE, points = TRUE)



trt.cols <- c(`C`= "#000000", #black - C
              `D`= "#0a8754", #g - D
              `I`= "#4472ca", #blu - I
              `E`= "#e84855") 

edat.pred <- edat %>% 
  group_by(mesocosm, ExpDay, Treatment) %>%
  add_epred_draws(m2.e,
                  re_formula = NA,
  ) 

edat.pred %>% 
  ggplot(
    ., 
         aes(x = ExpDay,
                y = count,
                colour = Treatment)
  ) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = c(0.95)
                  # position = position_dodge(.5),
                  # linewidth = 2
  )+
  geom_point(
    edat,
         mapping = aes(x = ExpDay,
                 y = count,
                 colour = Treatment),
             alpha = .35
             # position = position_jitterdodge(dodge.width = .5),
  ) +
  facet_wrap(~ fun_grp
             , scales = "free_y"
  )+
  scale_color_manual(values = trt.cols)
  # scale_fill_manual(values = trt.cols)+
  # labs(title = "Heterotroph", 
  #      y = "<span style='font-size: 15pt'>Abundance </span>
  #        <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
  #      x = NULL)+
  # theme(axis.title.y = ggtext::element_markdown(),
  #       axis.title.x = element_text(size = 15),
  #       plot.title = element_text(size=16, face="italic")
  # ) 
