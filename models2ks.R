library(tidyverse)
library(brms)
library(tidybayes)
library(emmeans)
library(modelr)
# library(see)

load("./data/sites_FC_phyto.RData")

# ----- Erken ------
functional.erk <- read.csv("./data/Erken_functional.csv", header = T, sep = ",")

# taxon level data
edat_tax <- erk.dat  |>  
  filter(# drop the lake
         Treatment != "ERK",
         # remove labels with low confidence
         ProbabilityScore > 0.5
  ) |> 
  mutate(Treatment = factor(Treatment, levels = c("C","D","I","E")),
         vol.offset = volume_run*concfactor,
         label = case_when(label == "ConSta000_2647648" ~ "ZygSta000_2647648",
                           label == "Heliozoa_0000000" ~ "ChlGolrad_2641100",
                           .default = label),
         .keep = "unused") |>
  group_by(ExpDay, Treatment, mesocosm, Name, label) |>
  summarise(# get the count of each taxon
            count = n(),
            # get the average properties of ech taxon per run
            abd = mean(AbdDiameter),
            len = mean(Length),
            cellvol = mean(biovolMS), # cubic micrometers
            sa = mean(surfacearea),
            prob = mean(ProbabilityScore),
            vol.offset = mean(vol.offset) # mL
            ) |> 
  group_by(ExpDay, Treatment, mesocosm, label) |>
  summarise(count = mean(count),
            # average the two flowcam runs
            abd = mean(abd),
            len = mean(len),
            cellvol = mean(cellvol), # cubic micrometers
            sa = mean(sa),
            prob = mean(prob),
            vol.offset = mean(vol.offset) # mL
            ) |> 
  ungroup() |> 
  left_join(functional.erk, by = join_by(label == Code)) |>
  mutate(class = substr(label, 1, 3),
         fun_grp = case_when(label == "FILAMENTS" ~ "III",
                             label == "ChlGolrad_2641100" ~ "I",
                             label == "CyaPla000_3218374" ~ "III",
                             class == "Bac" ~ "VI",
                             .default = KRUK_MBFG)) |> 
  select(-taxon, -KRUK_MBFG) |> 
  drop_na() |> 
  mutate(taxonvol = count*cellvol) 


# community level data
edat_tot = edat_tax |> 
  group_by(ExpDay, Treatment, mesocosm) |> 
  summarise(count = sum(count),
            biovol = sum(taxonvol),
            vol.offset = mean(vol.offset)) |> 
  mutate(dens = count/vol.offset,
         biovoldens = biovol/vol.offset # cubic micro per mL
         ) |> 
  ungroup()


# functional group level data
edat_fgroup = edat_tax |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(count = sum(count),
            biovol = sum(taxonvol),
            vol.offset = mean(vol.offset)) |> 
  mutate(dens = count/vol.offset,
         biovoldens = biovol/vol.offset,
         # ExpDay = as.factor(ExpDay),
         fun_grp = as.factor(fun_grp)) |> 
  ungroup()

## ---- models ----
### total biovolume ####

mtot = brm(
  bf(
    log10(biovoldens) ~ s(ExpDay, by = Treatment, k = 5) +
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


mtot |> emmeans("Treatment", by = "ExpDay",
                at = list(ExpDay = c(0,4,12,20,28,36))) |> 
  contrast(method = "trt.vs.ctrl") |>
  gather_emmeans_draws() |> 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
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

### functional group biovolume ####

mgroup = brm(
  bf(
    log10(biovoldens) ~ fun_grp*Treatment + 
      s(ExpDay, by = interaction(fun_grp, Treatment), k = 5) +
      (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
    sigma ~ fun_grp
  ),
  family = gaussian(),
  # prior = prior(normal(0,5), class = "b")+
  #   prior(exponential(2), class = "sd"),
  init = 0,
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.99
                 # max_treedepth = 12
                 ),
  seed = 543,
  backend = "cmdstanr",
  data = edat_fgroup
) # 79 min

saveRDS(mgroup, "models/Erk_mfgroup_fcdat.rds")
mgroup.e <- readRDS("models/Erk_mfgroup_fcdat.rds")
# mgroup.e <- readRDS("mgroup6000_fcdat.RDS")

pp_check(mgroup.e, ndraws = 100)
summary(mgroup.e) 

# default_prior(bf(
#   log10(biovoldens) ~ s(ExpDay, by = interaction(fun_grp, Treatment), k = 5) +
#     (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
#   sigma ~ fun_grp
#                 ),
#   family = gaussian(),
#   data = edat_fgroup)


#### plots ----

unique(edat_fgroup[,1:4])  |> 
  add_epred_draws(mgroup.e, re_formula = NULL) |>
  ggplot(aes(x = ExpDay, y = log10(biovoldens), fill = Treatment)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  geom_point(data = edat_fgroup, 
             color = "black", alpha = .5, size = 1,
             position = position_jitter(.9)) +
  scale_fill_manual(values = c("#00000030",
                               "#45766130",
                               "#5e8fcd30",
                               "#e79f4f30"))+
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab(expression(log(mu*m^3/mL))) +
  xlab("Experimental Day") +
  theme_light() +
  #theme(legend.position = "none") +
  facet_grid(vars(Treatment), 
             vars(fun_grp))


efgroup.pred <- edat_fgroup |> 
  data_grid(Treatment = unique(edat_fgroup$Treatment),
            ExpDay = seq_range(edat_fgroup$ExpDay, n = 25),
            fun_grp = unique(edat_fgroup$fun_grp),
            mesocosm = unique(edat_fgroup$mesocosm)) |> 
  add_epred_draws(mgroup.e, re_formula = NULL) 


(plot1 <- edat_fgroup %>%
    ggplot(aes(x = ExpDay,
               y = log10(biovoldens),
               colour = Treatment)
    ) +
    geom_point(
      # edat_fgroup, mapping = aes(x = ExpDay,
      #                              y = log10(biovoldens),
      #                              colour = Treatment),
               alpha = .3,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    stat_lineribbon(efgroup.pred, mapping  = aes(y = (.epred),
                        fill = Treatment),
                    point_interval = "mean_qi",
                    .width = c(0.95),
                    # alpha = .5
                    # position = position_dodge(.5),
                    linewidth = .5
    )+
    facet_wrap(~ fun_grp
               , scales = "free_y"
    )+
    scale_color_manual(values = c("#000000",
                                  "#457661",
                                  "#5e8fcd",
                                  "#e79f4f")) +
    scale_fill_manual(values = c("#00000030",
                                 "#45766130",
                                 "#5e8fcd30",
                                 "#e79f4f30"))+
    labs(y = expression(log(mu*m^3/mL)),
         x = "Experimental day")+
    theme_bw())

contrasts.e <- efgroup.pred |> 
  ungroup() |> 
  mutate(ExpDay %in% unique(edat_comp$ExpDay)) |> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 
  as.data.frame() |>
  mutate(ExpDay = as.factor(ExpDay))


mgroup.e |> emmeans("Treatment", by = c("ExpDay", "fun_grp"),
                at = list(ExpDay = c(0,4,12,20,28,36))) |> 
  contrast(method = "trt.vs.ctrl") |>
  gather_emmeans_draws() |> 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "darkgrey") +
  stat_lineribbon(aes(y = .value, group = fun_grp), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  scale_fill_manual(values = c(
    # "#00000030",
                               "#45766130",
                               "#5e8fcd30",
                               "#e79f4f30"), guide = F) +
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab("Difference") +
  xlab("Experimental Day") +
  theme_bw() +
  ggh4x::facet_grid2(fun_grp ~ contrast) 


contrasts.e |> 
  # filter(fun_grp == "I") |> 
  rowwise() |> 
  mutate(sign = ifelse(between(0, .lower, .upper), "no", "yes")) |>  
  ungroup() |> 
  ggplot(aes(x = .epred, y = Treatment, 
             colour = sign,
             # alpha = sign
  )) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper),
                  linewidth  = .4,
                  # fatten = 3,
                  # width = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  ggh4x::facet_grid2(fun_grp ~ ExpDay, 
                     scales = "free", 
                     independent = "x",
                     # space = "free", 
                     # axes = "all_y"
  )+
  scale_x_continuous(breaks = 0)+
  scale_color_manual(values = c("grey50", "darkred"), guide = F)+
  # scale_alpha_manual(values = c(.3, 1))+
  ylab("contrast") +
  xlab("estimate") +
  theme_bw()


### compositional changes ###########################

edat_comp <- edat_fgroup |> 
  select(ExpDay, Treatment, mesocosm, fun_grp, biovoldens) |> 
  pivot_wider(names_from = "fun_grp",
              values_from = "biovoldens") |> 
  # select(-V) |> 
  mutate(tot_biov = rowSums(across(c(I,II,III,IV,V,VI,VII)), na.rm = T),
         across(c(I, II, III, IV,V,VI, VII),
                ~ . / tot_biov),
         V = replace_na(V, 1e-06),
         tot_biov = rowSums(across(c(I,II,III,IV,V,VI,VII)), na.rm = T),
         across(c(I, II, III, IV,V,VI, VII),
                ~ . / tot_biov)
         )

# make a 'list' column with all percentages
edat_comp$Y = with(edat_comp, cbind(I,II,III,IV,V,VI,VII))

# mcomp = brm(
#   bf(
#     Y ~ Treatment+ s(ExpDay, by = Treatment, k = 5) +
#       (ExpDay | mesocosm)
#   ),
#   family = dirichlet(),
#   chains = 4,
#   iter = 6000,
#   warmup = 3000,
#   cores = 4,
#   control = list(adapt_delta = 0.95),
#   seed = 543,
#   backend = "cmdstanr",
#   data = edat_comp
# ) # 35 min
pp_check(mcomp.e, ndraws = 100)
summary(mcomp.e)

 saveRDS(mcomp, "models/Erk_funct-comp_all1.rds") # 36 min

mcomp.e <- readRDS("models/Erk_funct-comp_all1.rds")

conditional_effects(mcomp.e, effects = "ExpDay",
                    re_formula = NULL,
                    conditions = data.frame(Treatment = c("C","D","I","E")),
                    categorical = T, points = T) 


edat_comp_plt <- edat_comp |> 
  pivot_longer(c(I, II, III, IV, V, VI, VII), 
               names_to = "fun_grp", 
               values_to = "perc") |> 
  select(-Y) |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(perc = sum(perc)) |> 
  ungroup()

comp.pred <- edat_comp_plt |>  
  data_grid(Treatment = unique(edat_comp$Treatment),
            ExpDay = seq_range(edat_comp$ExpDay, n = 25),
            mesocosm = unique(edat_comp$mesocosm)) |> 
  add_epred_draws(mcomp.e, re_formula = NULL) |> 
  rename(fun_grp = .category) |> 
  ungroup()


(plot <- comp.pred %>%
    ggplot(aes(x = ExpDay,
               y = .epred,
               colour = fun_grp)
    ) +
    stat_lineribbon(aes(y = (.epred),
                        fill = fun_grp),
                    .width = c(0.95),
                    alpha = .5
                    # position = position_dodge(.5),
                    # linewidth = 2
    )+
    geom_point(edat_comp_plt, mapping = aes(x = ExpDay,
                                   y = perc,
                                   colour = fun_grp),
               alpha = .35,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    facet_wrap(~ Treatment
               , scales = "free_y"
    )+
    # scale_color_manual(values = trt.cols) +
    # scale_fill_manual(values = fil.cols) +
labs(y = "Biovolume %",
     x = "Experimental day")+
    theme_bw())

comp_contrasts.e <- comp.pred |> 
  filter(ExpDay %in% unique(edat_comp$ExpDay)) |> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 
  as.data.frame() |> 
  mutate(ExpDay = as.factor(ExpDay)) 

comp_contrasts.e |> 
  # filter(fun_grp == "I") |> 
  rowwise() |> 
  mutate(sign = ifelse(between(0, .lower, .upper), "no", "yes")) |>  
  ungroup() |> 
  ggplot(aes(x = .epred, y = Treatment, 
             colour = sign,
             # alpha = sign
  )) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper),
                  linewidth  = .4,
                  # fatten = 3,
                  # width = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  ggh4x::facet_grid2(fun_grp ~ ExpDay, 
                     scales = "free", 
                     independent = "x",
                     # space = "free", 
                     # axes = "all_y"
  )+
  scale_x_continuous(breaks = 0)+
  scale_color_manual(values = c("grey50", "darkred"), guide = F)+
  # scale_alpha_manual(values = c(.3, 1))+
  ylab("contrast") +
  xlab("estimate") +
  theme_bw()




# ---- Bolmen -----

functional.bol <- read.csv("./data/Bolmen_functional.csv", header = T, sep = ",")


# taxon level data
bdat_tax <- bol.dat  |>  
  filter(# drop the lake
    Treatment != "BOL",
    # remove labels with low confidence
    ProbabilityScore > 0.5
  ) |> 
  mutate(Treatment = factor(Treatment, levels = c("C","D","I","E")),
         vol.offset = volume_run*concfactor,
         .keep = "unused") |>
  group_by(ExpDay, Treatment, mesocosm, Name, label) |>
  summarise(count = n(),
            # average the two flowcam runs
            abd = mean(AbdDiameter),
            len = mean(Length),
            cellvol = mean(biovolume), # cubic micrometers
            sa = mean(surfacearea),
            prob = mean(ProbabilityScore),
            vol.offset = mean(vol.offset) # mL
  ) |> 
  ungroup() |> 
  left_join(functional.bol, by = "label") |>
  select(-Name) |> 
  drop_na() |> 
  mutate(taxonvol = count*cellvol) 


# community level data
bdat_tot = bdat_tax |> 
  group_by(ExpDay, Treatment, mesocosm) |> 
  summarise(count = sum(count),
            biovol = sum(taxonvol),
            vol.offset = mean(vol.offset)) |> 
  mutate(dens = count/vol.offset,
         biovoldens = biovol/vol.offset # cubic micro per mL
  ) |> 
  ungroup()


# functional group level data
bdat_fgroup = bdat_tax |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(count = sum(count),
            biovol = sum(taxonvol),
            vol.offset = mean(vol.offset)) |> 
  mutate(dens = count/vol.offset,
         biovoldens = biovol/vol.offset,
         fun_grp = as.factor(fun_grp)) |> 
  ungroup()

## ---- models ----
### total biovolume ####

mtot = brm(
  bf(
    log10(biovoldens) ~ s(ExpDay, by = Treatment, k = 5) +
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
  data = bdat_tot
) 
pp_check(mtot, ndraws = 100)
summary(mtot) 


mtot |> emmeans("Treatment", by = "ExpDay",
                # at = list(ExpDay = c(0,4,12,20,28,36))
                at = list(ExpDay = seq_range(edat_tot$ExpDay, n = 31))
) |> 
  contrast(method = "trt.vs.ctrl") |>
  gather_emmeans_draws() |> 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
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



### functional group biovolume ####


# mgroup = brm(
#   bf(
#     log10(biovoldens) ~ fun_grp*Treatment + 
#       s(ExpDay, by = interaction(fun_grp, Treatment), k = 5) +
#       (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
#     sigma ~ fun_grp
#   ),
#   family = gaussian(),
#   # prior = prior(normal(0,4), class = "b")+
#   #   prior(exponential(2), class = "sd"),
#   init = 0,
#   chains = 4,
#   iter = 6000,
#   warmup = 3000,
#   cores = 4,
#   control = list(adapt_delta = 0.99
#                  # max_treedepth = 12
#                  ),
#   seed = 543,
#   backend = "cmdstanr",
#   data = bdat_fgroup
# ) # 34 min

saveRDS(mgroup, "models/Bol_mgroup_fcdat.rds")
mgroup.b = readRDS("models/Bol_mgroup_fcdat.rds")
pp_check(mgroup.b, ndraws = 100)
summary(mgroup.b) 


# default_prior(bf(
#   log10(biovoldens) ~ s(ExpDay, by = interaction(fun_grp, Treatment), k = 5) +
#     (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
#   sigma ~ fun_grp
# ),
# family = gaussian(),
# data = bdat_fgroup)

#### plots -----

unique(bdat_fgroup[,1:4])  |> 
  add_epred_draws(mgroup.b, re_formula = NULL) |>
  ggplot(aes(x = ExpDay, y = log10(biovoldens), fill = Treatment)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  geom_point(data = bdat_fgroup, 
             color = "black", alpha = .5, size = 1,
             position = position_jitter(.9)) +
  scale_fill_manual(values = c("#00000030",
                               "#45766130",
                               "#5e8fcd30",
                               "#e79f4f30")) +
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab(expression(log(mu*m^3/mL))) +
  xlab("Experimental Day") +
  theme_light() +
  #theme(legend.position = "none") +
  facet_grid(vars(Treatment), 
             vars(fun_grp))

bfgroup.pred <- bdat_fgroup |> 
  data_grid(Treatment = unique(edat_fgroup$Treatment),
            ExpDay = seq_range(edat_fgroup$ExpDay, n = 37),
            fun_grp = unique(edat_fgroup$fun_grp),
            mesocosm = unique(edat_fgroup$mesocosm)) |> 
  add_epred_draws(mgroup.b, re_formula = NULL) 


(plot2 <- bdat_fgroup %>%
    ggplot(aes(x = ExpDay,
               y = log10(biovoldens),
               colour = Treatment)
    ) +
    geom_point(
      # edat_fgroup, mapping = aes(x = ExpDay,
      #                              y = log10(biovoldens),
      #                              colour = Treatment),
      alpha = .3,
      position = position_jitterdodge(dodge.width = .5),
    ) +
    stat_lineribbon(bfgroup.pred, mapping  = aes(y = (.epred),
                                                 fill = Treatment),
                    point_interval = "mean_qi",
                    .width = c(0.95),
                    # alpha = .5
                    # position = position_dodge(.5),
                    linewidth = .5
    )+
    facet_wrap(~ fun_grp
               , scales = "free_y"
    )+
    scale_color_manual(values = c("#000000",
                                  "#457661",
                                  "#5e8fcd",
                                  "#e79f4f")) +
    scale_fill_manual(values = c("#00000030",
      "#45766130",
      "#5e8fcd30",
      "#e79f4f30"))+
    labs(y = expression(log(mu*m^3/mL)),
         x = "Experimental day")+
    theme_bw())


mgroup.b |> emmeans("Treatment", by = c("ExpDay", "fun_grp"),
                    at = list(ExpDay = c(0,4,12,20,28,36))) |> 
  contrast(method = "trt.vs.ctrl") |>
  gather_emmeans_draws() |> 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "darkgrey") +
  stat_lineribbon(aes(y = .value, group = fun_grp), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  scale_fill_manual(values = c(#"#00000030",
    "#45766130",
    "#5e8fcd30",
    "#e79f4f30"), guide = F) +
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab("Difference") +
  xlab("Experimental Day") +
  theme_bw() +
  ggh4x::facet_grid2(fun_grp ~ contrast) 


contrasts <- bfgroup.pred |> 
  filter(ExpDay %in% unique(bdat_fgroup$ExpDay)) |> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 
  as.data.frame() |> 
  mutate(ExpDay = as.factor(ExpDay)) 
  
contrasts |> 
  # filter(fun_grp == "I") |> 
  rowwise() |> 
  mutate(sign = ifelse(between(0, .lower, .upper), "no", "yes")) |>  
  ungroup() |> 
  ggplot(aes(x = .epred, y = Treatment, 
             colour = sign,
             # alpha = sign
             )) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper),
                  linewidth  = .4,
                  # fatten = 3,
                  # width = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  ggh4x::facet_grid2(fun_grp ~ ExpDay, 
                     scales = "free", 
                     independent = "x",
                     # space = "free", 
                     # axes = "all_y"
  )+
  scale_x_continuous(breaks = 0)+
  scale_color_manual(values = c("grey50", "darkred"), guide = F)+
  # scale_alpha_manual(values = c(.3, 1))+
  ylab("contrast") +
  xlab("estimate") +
  theme_bw()

### compositional changes ###########################

bdat_comp <- bdat_fgroup |> 
  select(ExpDay, Treatment, mesocosm, fun_grp, biovoldens) |> 
  pivot_wider(names_from = "fun_grp",
              values_from = "biovoldens") |> 
  # select(-V) |> 
  mutate(tot_biov = rowSums(across(c(I,II,III,IV,
                                     V,
                                     VI,VII)), na.rm = T),
         across(c(I, II, III, IV, 
                  V,
                  VI, VII),
                ~ . / tot_biov),
         I = replace_na(I, 1e-06),
         II = replace_na(II, 1e-06),
         V = replace_na(V, 1e-06),
         tot_biov = rowSums(across(c(I,II,III,IV,
                                     V,
                                     VI,VII)), na.rm = T),
         across(c(I, II, III, IV, 
                  V,
                  VI, VII),
                ~ . / tot_biov)
  )

# make a 'list' column with all percentages
bdat_comp$Y = with(bdat_comp, cbind(I,II,III,IV,
                                    V,
                                    VI,VII))

mcomp = brm(
  bf(
    Y ~ Treatment + s(ExpDay, by = Treatment, k = 5) +
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
  data = bdat_comp
) # 35 min
pp_check(mcomp, ndraws = 100)
summary(mcomp.b)

saveRDS(mcomp, "models/Bol_funct-comp_all1.rds") # 26 min


mcomp.b <- readRDS("models/Bol_funct-comp_all1.rds")

conditional_effects(mcomp.b, effects = "ExpDay",
                    re_formula = NA,
                    conditions = data.frame(Treatment = c("C","D","I","E")),
                    categorical = T, points = T)


bdat_comp_plt <- bdat_comp |> 
  pivot_longer(c(I, II, III, IV, V, VI, VII), 
               names_to = "fun_grp", 
               values_to = "perc") |> 
  select(-Y) |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(perc = sum(perc)) |> 
  ungroup()

compb.pred <- bdat_comp_plt |>  
  data_grid(Treatment = unique(bdat_comp$Treatment),
            ExpDay = seq_range(bdat_comp$ExpDay, n = 25),
            mesocosm = unique(bdat_comp$mesocosm)) |> 
  add_epred_draws(mcomp.b, re_formula = NULL) 


(plot <- compb.pred %>%
    ggplot(aes(x = ExpDay,
               y = .epred,
               colour = .category)
    ) +
    stat_lineribbon(aes(y = (.epred),
                        fill = .category),
                    .width = c(0.95),
                    point_interval = "mean_qi",
                    alpha = .5
                    # position = position_dodge(.5),
                    # linewidth = 2
    )+
    geom_point(bdat_comp_plt, mapping = aes(x = ExpDay,
                                   y = perc,
                                   colour = fun_grp),
               alpha = .35,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    facet_wrap(~ Treatment
               , scales = "free_y"
    )+
    # scale_color_manual(values = trt.cols) +
    # scale_fill_manual(values = fil.cols) +
    labs(y = "Biovolume %",
         x = "Experimental day")+
    theme_bw())


comp_contrasts.b <- compb.pred |> 
  filter(ExpDay %in% unique(bdat_comp$ExpDay)) |> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 
  as.data.frame() |> 
  mutate(ExpDay = as.factor(ExpDay)) 

comp_contrasts.b |> 
  # filter(fun_grp == "I") |> 
  rowwise() |> 
  mutate(sign = ifelse(between(0, .lower, .upper), "no", "yes")) |>  
  ungroup() |> 
  ggplot(aes(x = .epred, y = Treatment, 
             colour = sign,
             # alpha = sign
  )) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper),
                  linewidth  = .4,
                  # fatten = 3,
                  # width = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  ggh4x::facet_grid2(fun_grp ~ ExpDay, 
                     scales = "free", 
                     independent = "x",
                     # space = "free", 
                     # axes = "all_y"
  )+
  scale_x_continuous(breaks = 0)+
  scale_color_manual(values = c("grey50", "darkred"), guide = F)+
  # scale_alpha_manual(values = c(.3, 1))+
  ylab("contrast") +
  xlab("estimate") +
  theme_bw()
