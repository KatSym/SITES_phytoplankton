library(tidyverse)
library(brms)
library(emmeans)
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
         voldens = biovol/vol.offset) %>% 
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
  seed = 123,
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


mtot %>% emmeans("Treatment", by = "ExpDay",
                 at = list(ExpDay = c(0,4,12,20,28,36))
                 #at = list(ExpDay = seq_range(edat_tot$ExpDay, n = 31))
                 ) %>% 
  contrast(method = "trt.vs.ctrl") %>%
  gather_emmeans_draws() %>% 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  stat_lineribbon(aes(y = .value), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  scale_fill_manual(values = c(#"#00000030",
                               "#0a875430",
                               "#4472ca30",
                               "#e8485530")) +
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab("Difference") +
  xlab("Experimental Day") +
  theme_light() +
  facet_wrap(vars(contrast))


#### group-level biovolume #####################################################

mgroup = brm(
  bf(
    log(voldens) ~ s(ExpDay, by = interaction(fun_grp, 
                                              Treatment), 
                     k = 5) +
      (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
    sigma ~ fun_grp
  ),
  family = gaussian(),
  prior = prior(normal(0,5), class = "b")+
          prior(exponential(2), class = "sd"),
  init = 0,
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 12),
  seed = 543,
  backend = "cmdstanr", 
  data = edat_fgroup
) 
saveRDS(mgroup, "mgroup6000.RDS")
mgroup = readRDS("mgroup6000.RDS")
pp_check(mgroup, ndraws = 100)
summary(mgroup)

unique(edat_fgroup[,1:4])  %>% 
  add_epred_draws(mgroup, re_formula = NULL) %>%
  ggplot(aes(x = ExpDay, y = log(voldens), fill = Treatment)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  geom_point(data = edat_fgroup, 
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
  facet_grid(vars(Treatment), 
             vars(fun_grp))


default_prior(bf(
  log(voldens) ~ s(ExpDay, by = interaction(fun_grp, Treatment), k = 5) +
    (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
  sigma ~ fun_grp
),
family = gaussian(),
data = edat_fgroup)

######### compositional changes ################################################

edat_comp = edat_fgroup %>% 
  
  select(ExpDay, Treatment, mesocosm, fun_grp, voldens) %>% 
  pivot_wider(names_from = "fun_grp",
              values_from = "voldens") %>% 
  mutate(tot_biov = rowSums(across(c(I,II,III,IV,V,VI,VII)), na.rm = T)) %>% 
  mutate(across(c(I, II, III, IV, V, VI, VII),
                ~ . / tot_biov)) %>% 
  mutate(V = replace_na(V, 1e-06)) %>% 
  mutate(tot_biov = rowSums(across(c(I,II,III,IV,V,VI,VII)), na.rm = T)) %>% 
  mutate(across(c(I, II, III, IV, V, VI, VII),
                ~ . / tot_biov))

edat_comp$Y = with(edat_comp, cbind(I,II,III,IV,V,VI,VII))

mcomp = brm(
  bf(
    Y ~ s(ExpDay, by = Treatment, k = 5) +
       (ExpDay | mesocosm)
  ),
  family = dirichlet(),
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.95),
  seed = 543,
  backend = "cmdstanr", 
  data = edat_comp
) 
pp_check(mcomp, ndraws = 100)
summary(mcomp)
conditional_effects(mcomp, effects = "ExpDay",
                    re_formula = NA,
                    conditions = data.frame(Treatment = c("C","D","I","E")),
                    categorical = T)

mcomp = brm(
  bf(
    Y ~ s(ExpDay, by = Treatment, k = 5) +
      (ExpDay | mesocosm)
  ),
  family = logistic_normal(),
  prior = prior("lkj(3)", "lncor"),
  chains = 4,
  iter = 12000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.99),
  seed = 543,
  backend = "cmdstanr", 
  data = edat_comp
) 
pp_check(mgroup, ndraws = 100)
summary(mgroup)
