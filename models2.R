library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)
library(see)

load("./data/sites_FC_phyto.RData")


functional.erk <- read.csv("./data/Erken_functional.csv", header = T, sep = ",")
species_traits <- read.csv("./data/species_traits.csv", header = T, sep = ",")

# taxon level data
edat_tax <- erk.dat %>% 
  filter(# drop the lake
         Treatment != "ERK",
         # remove labels with low confidence
         ProbabilityScore > 0.5
  ) %>% 
  mutate(Treatment = factor(Treatment, levels = c("C","D","I","E")),
         vol.offset = volume_run*concfactor,
         label = ifelse(label == "ConSta000_2647648", "ZygSta000_2647648", label),
         .keep = "unused") %>%
  group_by(ExpDay, Treatment, mesocosm, Name, label) %>%
  summarise(count = n(),
            # average the two flowcam runs
            abd = mean(AbdDiameter),
            len = mean(Length),
            cellvol = mean(biovolMS), # cubic micrometers
            sa = mean(surfacearea),
            prob = mean(ProbabilityScore),
            vol.offset = mean(vol.offset) # mL
            ) %>% 
  ungroup() %>% 
  left_join(., functional.erk, by = join_by(label == Code)) %>%
  mutate(class = substr(label, 1, 3),
         fun_grp = case_when(label == "FILAMENTS" ~ "III",
                             label == "CyaPla000_3218374" ~ "III",
                             class == "Bac" ~ "VI",
                             .default = KRUK_MBFG)) %>% 
  select(-label, -KRUK_MBFG, -Name) %>% 
  drop_na() %>% 
  mutate(taxonvol = count*cellvol) 


# community level data
edat_tot = edat_tax %>% 
  group_by(ExpDay, Treatment, mesocosm) %>% 
  summarise(count = sum(count),
            biovol = sum(taxonvol),
            vol.offset = mean(vol.offset)) %>% 
  mutate(dens = count/vol.offset,
         voldens = biovol/vol.offset # cubic micro per mL
         ) %>% 
  ungroup()


# group level data
edat_fgroup = edat_tax %>% 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) %>% 
  summarise(count = sum(count),
            biovol = sum(taxonvol),
            vol.offset = mean(vol.offset)) %>% 
  mutate(dens = count/vol.offset,
         voldens = biovol/vol.offset,
         ExpDay = as.factor(ExpDay)) %>% 
  ungroup()


### total biovolume #####

mtot = brm(
  bf(
    log(voldens) ~ s(ExpDay, by = Treatment, k = 5) +
      (ExpDay|mesocosm)
  ),
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.99),
  seed = 543,
  backend = "cmdstanr", 
  data = edat_tot
) 
pp_check(mtot, ndraws = 100)
summary(mtot)
plot(conditional_effects(mtot, effects = "ExpDay:Treatment",
                         re_formula = NULL), points = T)


unique(edat_tot[,1:3])  %>% 
  add_epred_draws(mtot, re_formula = NULL) %>%
  ggplot(aes(x = ExpDay, y = log(voldens), fill = Treatment)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  geom_point(data = edat_tot, 
             color = "black", alpha = .5, size = 1,
             position = position_jitter(.9)) +
  scale_fill_manual(values = c("#00000030",
                               "#0a875430",
                               "#4472ca30",
                               "#e8485530")) +
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab(expression(log(mu*m^3/mL))) +
  xlab("Experimental Day") +
  theme_light() +
  #theme(legend.position = "none") +
  facet_wrap(vars(Treatment))

library(emmeans)
library(modelbased)


contrasts = estimate_contrasts(
  mtot,
  contrast = "Treatment",
  by = "ExpDay",
  method = "trt.vs.ctrl",
  length = 50,
  backend = "emmeans"
)
# Add Contrast column by concatenating
contrasts$Contrast = paste(contrasts$Level1, "-", contrasts$Level2)

ggplot(contrasts, aes(x = ExpDay, y = Difference, )) +
  # Add line and CI band
  geom_line(aes(color = Contrast)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = Contrast), alpha = 0.2) +
  # Add line at 0, indicating no difference
  geom_hline(yintercept = 0, linetype = "dashed") +
  # Colors
  theme_light()
