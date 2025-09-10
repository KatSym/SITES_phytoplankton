library(tidyverse)
library(brms)
library(tidybayes)
library(emmeans)
library(modelr)
# library(see)

load("./data/sites_FC_phyto.RData")

# for plots
trt.cols <- c(`C`= "#000000", #black - C
              `D`= "#0a8754", #g - D
              `I`= "#4472ca", #blu - I
              `E`= "#e84855") 

fil.cols <- c(`C`= "#00000030", #black - C
              `D`= "#0a875430", #g - D
              `I`= "#4472ca30", #blu - I
              `E`= "#e8485530") 


# ----- Erken ------
functional.erk <- read.csv("./data/Erken_functional.csv", header = T, sep = ",") |> 
  mutate(KRUK_MBFG = case_when(Code == "CyaDollem_7901634" ~ "VIII",
                               Code == "CyaDol000_5428755" ~ "VIII",
                               Code == "CyaGloech_0000000" ~ "VIII",
                               .default = KRUK_MBFG))

large_phyto <- read.csv("data/Erken_large_phyto.csv") |> 
  pivot_longer(cols = 4:7, names_sep = "_", names_to = c("taxon", "variable"), values_to = "val") |> 
  pivot_wider(id_cols = c(day, Treatment, mesocosm, taxon), names_from = variable, values_from = val) |> 
  mutate(label = ifelse(taxon == "gloe", "CyaGloech_0000000", "BacFracro_3192403"),
         fun_grp = ifelse(taxon == "gloe", "VIII", "VI"), .keep = "unused") |> 
  rename(ExpDay = day,
         dens = abund,
         biovoldens = biovol)

small_phyto <- read.csv("data/WRONG_small_phyto.csv") |> 
  rename(dens = abund,
         biovoldens = biovolume) %>% 
  select(-c(Species, Genus, size_fract, vol))

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
            # get the average properties of each taxon per run
            abd = mean(AbdDiameter),
            len = mean(Length),
            cellvol = mean(biovolMS), # cubic micrometers
            sa = mean(surfacearea),
            prob = mean(ProbabilityScore),
            vol.offset = mean(vol.offset) # mL
            ) |> 
  group_by(ExpDay, Treatment, mesocosm, label) |>
            # average the two flowcam runs
  summarise(count = mean(count),
            cellvol = mean(cellvol), # cubic micrometers
            vol.offset = mean(vol.offset) # mL
            ) |> 
  mutate(taxonvol = count*cellvol,
         dens = count/vol.offset,
         biovoldens = taxonvol/vol.offset) %>% 
  ungroup() |> 
  left_join(functional.erk, by = join_by(label == Code)) |>
  mutate(class = substr(label, 1, 3),
         fun_grp = case_when(label == "FILAMENTS" ~ "III",
                             label == "ChlGolrad_2641100" ~ "I",
                             label == "CyaPla000_3218374" ~ "III",
                             label == "ZygSta000_2647648" ~ "IV",
                             label == "CyaGlo000_filaments" ~ "VIII",
                             class == "Bac" ~ "VI",
                             .default = KRUK_MBFG)) |> 
  # select(-taxon, -KRUK_MBFG) |> 
  select(-taxon, -KRUK_MBFG, -vol.offset, -count, -taxonvol, -cellvol, -class) |> 
  drop_na()
  

edat_tax <- edat_tax %>% 
  rbind(large_phyto) %>% 
  rbind(small_phyto) %>% 
  group_by(ExpDay, Treatment, mesocosm, fun_grp, label) |>
  summarise(dens = sum(dens),
            biovoldens = sum(biovoldens)) %>%  # cubic micrometers
  ungroup()


# community level data
edat_tot = edat_tax |> 
  group_by(ExpDay, Treatment, mesocosm) |> 
  summarise(dens = sum(dens),
            biovoldens = sum(biovoldens)) |> 
  ungroup()


# functional group level data
edat_fgroup = edat_tax |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(dens = sum(dens),
            biovoldens = sum(biovoldens)) |> 
  mutate(fun_grp = as.factor(fun_grp)) |> 
  ungroup() %>% 
  mutate(Treatment = factor(Treatment, levels = c("C", "D", "I", "E")),
         fun_grp = factor(fun_grp, 
                          levels = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII")),
         # replace 0s with NAs - model doesn't run otherwise
         # across(biovoldens,  ~replace(.x, .x == 0, NA))
         )
  
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
saveRDS(mtot, "models/Erk_total_biovol.rds")
pp_check(mtot, ndraws = 100)


mtot |> emmeans("Treatment", by = "ExpDay",
                at = list(ExpDay = c(0,4,12,20,28,36))) |> 
  contrast(method = "trt.vs.ctrl") |>
  gather_emmeans_draws() |> 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
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

etot.pred <- edat_tot |> 
  data_grid(Treatment = unique(edat_fgroup$Treatment),
            ExpDay = seq_range(edat_fgroup$ExpDay, n = 20),
            fun_grp = unique(edat_fgroup$fun_grp),
            mesocosm = unique(edat_fgroup$mesocosm)) |> 
  add_epred_draws(mtot, re_formula = NULL) %>% 
  ungroup()

(p1 <- edat_tot %>%
    ggplot(aes(x = ExpDay,
               y = log10(biovoldens),
               colour = Treatment)
    ) +
    geom_point(
      alpha = .3,
      position = position_jitterdodge(dodge.width = .5),
    ) +
    stat_lineribbon(etot.pred, mapping  = aes(y = (.epred),
                                                 fill = Treatment),
                    point_interval = "mean_qi",
                    .width = c(0.95),
                    # alpha = .5
                    # position = position_dodge(.5),
                    linewidth = .7
    )+
    # scale_color_manual(values = c("#000000",
    #                               "#457661",
    #                               "#5e8fcd",
    #                               "#e79f4f")) +
    # scale_fill_manual(values = c("#00000030",
    #                              "#45766130",
    #                              "#5e8fcd30",
    #                              "#e79f4f30"))+
    scale_color_manual(values = trt.cols) +
    scale_fill_manual(values = fil.cols)+
    labs(y = expression(log(mu*m^3/mL)),
         x = "Experimental day")+
    theme_bw() + 
    theme(strip.text = element_text(face="bold")))

### functional group biovolume ####

mgroup.e = brm(
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
  # opencl = opencl(c(0, 0)),
  iter = 6000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.99
                 # max_treedepth = 12
                 ),
  seed = 543,
  backend = "cmdstanr",
  data = edat_fgroup) # 79 min
beepr::beep(1)

saveRDS(mgroup.e, "models/Erk_mfgroup_fcdat_wrong-small.rds")
mgroup.e <- readRDS("models/Erk_mfgroup_fcdat_wrong-small.rds")
# mgroup.e <- readRDS("models/Erk_mfgroup_fcdat.rds")

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
  add_epred_draws(mgroup.e, re_formula = NULL) %>% 
  ungroup() %>% 
  mutate(Treatment = factor(Treatment, levels = c("C", "D", "I", "E")),
         fun_grp - factor(fun_grp, 
                          levels = c("I", "II", "III", "IV", "V", "VI", "VII", "Gloe"))) |>
  ggplot(aes(x = ExpDay, y = log10(biovoldens), fill = Treatment)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  geom_point(data = edat_fgroup, 
             color = "black", alpha = .6, size = 1,
             position = position_jitter(.9)) +
  scale_fill_manual(values = c("#00000050",
                               "#45766150",
                               "#5e8fcd50",
                               "#e79f4f50"
                               ))+
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab(expression(log(mu*m^3/mL))) +
  xlab("Experimental Day") +
  theme_light() +
  #theme(legend.position = "none") +
  facet_grid(vars(Treatment), 
             vars(fun_grp))


efgroup.pred <- edat_fgroup |> 
  data_grid(Treatment = unique(edat_fgroup$Treatment),
            ExpDay = seq_range(edat_fgroup$ExpDay, n = 20),
            fun_grp = unique(edat_fgroup$fun_grp),
            mesocosm = unique(edat_fgroup$mesocosm)) |> 
  add_epred_draws(mgroup.e, re_formula = NULL) %>% 
  ungroup()


(plot1 <- edat_fgroup %>%
    ggplot(aes(x = ExpDay,
               y = log10(biovoldens),
               colour = Treatment)
    ) +
    geom_point(
               alpha = .3,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    stat_lineribbon(efgroup.pred, mapping  = aes(y = (.epred),
                        fill = Treatment),
                    point_interval = "mean_qi",
                    .width = c(0.95),
                    # alpha = .5
                    # position = position_dodge(.5),
                    linewidth = .7
    )+
    facet_wrap(~ fun_grp
               , scales = "free_y",
               nrow = 4, ncol = 2
    )+
    # scale_color_manual(values = c("#000000",
    #                               "#457661",
    #                               "#5e8fcd",
    #                               "#e79f4f")) +
    # scale_fill_manual(values = c("#00000030",
    #                              "#45766130",
    #                              "#5e8fcd30",
    #                              "#e79f4f30"))+
    scale_color_manual(values = trt.cols) +
    scale_fill_manual(values = fil.cols)+
    labs(y = expression(log(mu*m^3/mL)),
         x = "Experimental day")+
    theme_bw() + 
    theme(strip.text = element_text(face="bold")))




# treatment - control differences
mgroup.e |> emmeans("Treatment", by = c("ExpDay", "fun_grp"),
                    at = list(ExpDay = c(0,4,12,20,28,36))) |> 
  contrast(method = "trt.vs.ctrl") |>
  gather_emmeans_draws() |> 
  ggplot(aes(x = ExpDay, y = .value, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey20") +
  stat_lineribbon(aes(y = .value, group = fun_grp), .width = c(0.95), #alpha = .5,
                  point_interval = "mean_qi",
                  linewidth = .25) + 
  scale_fill_manual(values = c(
                                # "#00000030",
                                "#0a875440", 
                                "#4472ca40", 
                                "#e8485540"
                              ), guide = "none") +
  scale_x_continuous(breaks = c(0,4,12,20,28,36)) +
  ylab("Difference") +
  xlab("Experimental Day") +
  theme_bw() +
  ggh4x::facet_grid2(fun_grp ~ contrast, 
                     scales = "free_y", 
                     independent = "y") 

# all contrasts
contrasts.e <- edat_fgroup |> 
  data_grid(Treatment = unique(edat_fgroup$Treatment),
            ExpDay = unique(edat_fgroup$ExpDay),
            fun_grp = unique(edat_fgroup$fun_grp),
            mesocosm = unique(edat_fgroup$mesocosm)) |> 
  add_epred_draws(mgroup.e, re_formula = NULL) %>% 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 
  as.data.frame() |>
  mutate(ExpDay = as.factor(ExpDay))

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
  # relocate("Gloe", .after = "VII") %>% 
  # select(-V) |> 
  mutate(
    # Gloe = ifelse(Gloe == 0, NA, Gloe),
    #      Treatment = factor(Treatment, levels = c("C", "D", "I", "E")),
    tot_biov = rowSums(across(c(I,II,III,IV,V,VI,VII, VIII)), na.rm = T),
    across(c(I, II, III, IV,V,VI, VII, VIII),
           ~ . / tot_biov)) %>% 
    mutate(V = ifelse(V == 0, 1e-06, V),
    # Gloe = replace_na(Gloe, 1e-06),
    V = replace_na(V, 1e-06)
    ) %>% 
  mutate(tot_biov = rowSums(across(c(I,II,III,IV,V,VI,VII, VIII)), na.rm = T),
         across(c(I, II, III, IV,V,VI, VII, VIII),
                ~ . / tot_biov) )
              

# make a 'list' column with all percentages
edat_comp$Y = with(edat_comp, cbind(I,II,III,IV,V,VI,VII, Gloe))

mcomp = brm(
  bf(
    Y ~ Treatment+ s(ExpDay, by = Treatment, k = 5) +
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
) # 35 min
summary(mcomp)

 saveRDS(mcomp, "models/Erk_funct-comp_wrong-small.rds") # 36 min
 mcomp.e <- readRDS("models/Erk_funct-comp_wrong-small.rds")
mcomp.e <- readRDS("models/Erk_funct-comp_all1.rds")

conditional_effects(mcomp.e, effects = "ExpDay",
                    re_formula = NULL,
                    conditions = data.frame(Treatment = c("C","D","I","E")),
                    categorical = T, points = T) 


edat_comp_plt <- edat_comp |> 
  pivot_longer(c(I, II, III, IV, V, VI, VII, Gloe), 
               names_to = "fun_grp", 
               values_to = "perc") |> 
  select(-Y) |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(perc = sum(perc)) |> 
  ungroup()

comp.pred <- edat_comp_plt |>  
  data_grid(Treatment = unique(edat_comp$Treatment),
            ExpDay = seq_range(edat_comp$ExpDay, n = 20),
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
    facet_wrap(~ Treatment,
               scales = "free_y",
               nrow = 1,
               labeller = labeller(Treatment = c("C" = "Control", 
                                                 "D" = "Daily", 
                                                 "I" = "Intermittent", 
                                                 "E" = "Extreme"))
    )+
    ggokabeito::scale_fill_okabe_ito(name = "Functional groups") +
    ggokabeito::scale_color_okabe_ito(guide = F) +
    # scale_color_manual(values = trt.cols) +
    # scale_fill_manual(values = fil.cols) +
    labs(y = "Biovolume %",
         x = "Experimental day")+
    theme_bw() +
    theme(strip.text = element_text(face="bold"))
)

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


bol_small <- read.csv("data/bolmen_counts.csv") |> 
  select(Species, Day, Mesocosm, Treatment, cell_abundanceperml, totalbvperml) |> 
  rename(ExpDay = Day,
         mesocosm = Mesocosm, 
         dens = cell_abundanceperml,
         biovoldens = totalbvperml) |> 
  mutate(fun_grp = case_when(Species == "Aphanocapsa" ~ "VII",
                           Species == "Aphanothece" ~ "VII",
                           Species == "Asteroccoccus" ~ "I",
                           Species == "Chlamydomonas" ~ "I",
                           Species == "Chroomonas" ~ "I",
                           Species == "Chroomonas acuta" ~ "I",
                           Species == "Chrysochromulina parva" ~ "I",
                           Species == "Closterium" ~ "IV",
                           Species == "Coelastrum" ~ "IV",
                           Species == "Coenochloris sp " ~ "I",#
                           Species == "Cosmarium" ~ "IV",
                           Species == "Cosmarium bioculatum" ~ "IV",
                           Species == "Cosmarium crenatum" ~ "IV",
                           Species == "Cosmarium meneghinii" ~ "IV",
                           Species == "Cosmarium regnelli" ~ "IV",
                           Species == "Crucigenia tetrapedia" ~ "IV",
                           Species == "Cyclotella" ~ "VI",
                           Species == "Desmodesmus" ~ "IV",
                           Species == "Dictyosphaerium pulchellum" ~ "I",
                           Species == "Cryptomonas" ~ "V",
                           Species == "Cyanodictyon reticulatum" ~ "I",#
                           Species == "Elakatothrix" ~ "IV", 
                           Species == "Euastrum binale" ~ "IV",
                           Species == "Eudorina" ~ "VII",
                           Species == "Eunotia" ~ "VI",
                           Species == "Frustulia" ~ "VI",
                           Species == "Golenkinia" ~ "IV",
                           Species == "Gomphonema" ~ "VI",
                           Species == "Gonium pectorale" ~ "I",# ot IV
                           Species == "Kirchneriella" ~ "IV",
                           Species == "Kirchneriella lunaris" ~ "IV",
                           Species == "Mallomonas" ~ "II",
                           Species == "Merismopedia" ~ "I",
                           Species == "Micractinium" ~ "I",#
                           Species == "Micractinium pusillum" ~ "I",
                           Species == "Monoraphidium" ~ "IV",
                           Species == "Monoraphidium contortum" ~ "IV",
                           Species == "Navicula" ~ "VI",
                           Species == "Oocystis" ~ "IV",
                           Species == "Paulschulzia" ~ "I",#
                           Species == "Pediastrum tetras" ~ "IV",
                           Species == "Pennate diatom" ~ "VI",
                           Species == "Pseudanabaena" ~ "III",
                           Species == "Rhodomonas" ~ "I",
                           Species == "Scenedesmus" ~ "IV",#
                           Species == "Selenastrum " ~ "I",
                           Species == "Selenastrum bibraianum" ~ "I",#
                           Species == "Small dinoflagellate" ~ "V",#
                           Species == "Synedra" ~ "VI",
                           Species == "Synura" ~ "II",
                           Species == "Trachelomonas" ~ "V"#
                           ),
         label = case_when(Species == "Aphanocapsa" ~ "CyaAph000_0000000",#
                           Species == "Aphanothece" ~ "CyaAph000_3216476",#
                           Species == "Asteroccoccus" ~ "ChlAst000_2639495",#
                           Species == "Chlamydomonas" ~ "ChlChl000_5271044",
                           Species == "Chroomonas" ~ "CryChr000_3202572",
                           Species == "Chroomonas acuta" ~ "CryChracu_0000000",#
                           Species == "Chrysochromulina parva" ~ "PryChrpar_3202214",# 
                           Species == "Closterium" ~ "ZygClo000_2646356",
                           Species == "Coelastrum" ~ "ChlCoe000_7749802",
                           Species == "Coenochloris sp " ~ "ChlCoe000_2641720",#
                           Species == "Cosmarium" ~ "ZygCos000_2648709",
                           Species == "Cosmarium bioculatum" ~ "ZygCosbio_5274000",
                           Species == "Cosmarium crenatum" ~ "ZygCoscre_0000000",#
                           Species == "Cosmarium meneghinii" ~ "ConCosmen_0000000",# different order?
                           Species == "Cosmarium regnelli" ~ "ConCosreg_5274551",#
                           Species == "Crucigenia tetrapedia" ~ "TreCrutet_8342650",#
                           Species == "Cyclotella" ~ "BacCyc000_3193095",
                           Species == "Desmodesmus" ~ "ChlDes000_2652440",
                           Species == "Dictyosphaerium pulchellum" ~ "TreDicpul_2642006",#
                           Species == "Cryptomonas" ~ "CryCry000_3202413",
                           Species == "Cyanodictyon reticulatum" ~ "CyaCyaret_0000000",#
                           Species == "Elakatothrix" ~ "KleEla000_0000000", #
                           Species == "Euastrum binale" ~ "ConEuabin_2648463",#
                           Species == "Eudorina" ~ "ChlEud000_0000000",
                           Species == "Eunotia" ~ "BacEun000_0000000",#
                           Species == "Frustulia" ~ "BacFru000_0000000",#
                           Species == "Golenkinia" ~ "ChlGol000_2641093",#
                           Species == "Gomphonema" ~ "BacGom000_7592153",
                           Species == "Gonium pectorale" ~ "ChlGonpec_2639286",#
                           Species == "Kirchneriella" ~ "ChlKir000_2641223",
                           Species == "Kirchneriella lunaris" ~ "ChlKirlun_5271638",
                           Species == "Mallomonas" ~ "ChrMal000_3195066",
                           Species == "Merismopedia" ~ "CyaMer000_0000000",
                           Species == "Micractinium" ~ "TreMic000_0000000",#
                           Species == "Micractinium pusillum" ~ "TreMicpus_2641082",
                           Species == "Monoraphidium" ~ "ChlMon000_0000000",
                           Species == "Monoraphidium contortum" ~ "ChlMoncon_2641582",
                           Species == "Navicula" ~ "BacNav000_2634825",
                           Species == "Oocystis" ~ "TreOoc000_2641454",
                           Species == "Paulschulzia" ~ "ChlPau000_2639505",
                           Species == "Pediastrum tetras" ~ "ChlPedtet_0000000",#
                           Species == "Pennate diatom" ~ "Bac000000_0000000",
                           Species == "Pseudanabaena" ~ "CyaPse000_3218042",#
                           Species == "Rhodomonas" ~ "CryRho000_3202546",
                           Species == "Scenedesmus" ~ "ChlSyn000_2640664",#
                           Species == "Selenastrum " ~ "ChlSel000_0000000",
                           Species == "Selenastrum bibraianum" ~ "ChlSelbib_2641190",#
                           Species == "Small dinoflagellate" ~ "Din000000_0000000",#
                           Species == "Synedra" ~ "BacSyn000_3192703",
                           Species == "Synura" ~ "ChrSyn000_0000000",
                           Species == "Trachelomonas" ~ "EugTra000_3208889"#
         ),
         .keep = "unused"
         )


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
  summarise(
    # get taxon count
    count = n(),
    # get the average properties of each taxon per run
    abd = mean(AbdDiameter),
    len = mean(Length),
    cellvol = mean(biovolume), # cubic micrometers
    sa = mean(surfacearea),
    prob = mean(ProbabilityScore),
    vol.offset = mean(vol.offset) # mL
  ) |> 
  group_by(ExpDay, Treatment, mesocosm, label) |>
  # average the two flowcam runs
  summarise(count = mean(count),
            cellvol = mean(cellvol), # cubic micrometers
            vol.offset = mean(vol.offset) # mL
  ) |> 
  mutate(taxonvol = count*cellvol,
         dens = count/vol.offset,
         biovoldens = taxonvol/vol.offset,
         .keep = "unused") %>% 
  ungroup() |> 
  left_join(functional.bol, by = "label") |>
  drop_na() |> 
  select(-taxonvol, -class)


bdat_tax <- bdat_tax %>% 
  rbind(bol_small) %>% 
  group_by(ExpDay, Treatment, mesocosm, fun_grp, label) |>
  summarise(dens = sum(dens),
            biovoldens = sum(biovoldens)) %>%  # cubic micrometers
  ungroup()



# community level data
bdat_tot = bdat_tax |> 
  group_by(ExpDay, Treatment, mesocosm) |> 
  summarise(dens = sum(dens),
         biovoldens = sum(biovoldens, na.rm = T) # cubic micro per mL
  ) |> 
  ungroup()


# functional group level data
bdat_fgroup = bdat_tax |> 
  group_by(ExpDay, Treatment, mesocosm, fun_grp) |> 
  summarise(dens = sum(dens),
            biovoldens = sum(biovoldens, na.rm = T) # cubic micro per mL
  ) |>
  ungroup() |> 
  mutate(Treatment = factor(Treatment, levels = c("C", "D", "I", "E")),
         fun_grp = factor(fun_grp, 
                          levels = c("I", "II", "III", "IV", "V", "VI", "VII"))
  )

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


mgroup = brm(
  bf(
    log10(biovoldens) ~ fun_grp*Treatment +
      s(ExpDay, by = interaction(fun_grp, Treatment), k = 5) +
      (ExpDay + fun_grp + ExpDay:fun_grp | mesocosm),
    sigma ~ fun_grp
  ),
  family = gaussian(),
  # prior = prior(normal(0,4), class = "b")+
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
  data = bdat_fgroup
) # 34 min

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
  data_grid(Treatment = unique(bdat_fgroup$Treatment),
            ExpDay = seq_range(bdat_fgroup$ExpDay, n = 20),
            fun_grp = unique(bdat_fgroup$fun_grp),
            mesocosm = unique(bdat_fgroup$mesocosm)) |> 
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
                    linewidth = .7
    )+
    facet_wrap(~ fun_grp
               , scales = "free_y",
               ncol = 2
    )+
    # scale_color_manual(values = c("#000000",
    #                               "#457661",
    #                               "#5e8fcd",
    #                               "#e79f4f")) +
    # scale_fill_manual(values = c("#00000030",
    #   "#45766130",
    #   "#5e8fcd30",
    #   "#e79f4f30"))+
    scale_color_manual(values = trt.cols) +
    scale_fill_manual(values = fil.cols)+
    labs(y = expression(log(mu*m^3/mL)),
         x = "Experimental day")+
    theme_bw() + 
    theme(strip.text = element_text(face="bold")))


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
         VII = replace_na(VII, 1e-06),
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
summary(mcomp.b)

saveRDS(mcomp, "models/Bol_funct-comp_all.rds") # 26 min


mcomp.b <- readRDS("models/Bol_funct-comp_all.rds")

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
            ExpDay = seq_range(bdat_comp$ExpDay, n = 20),
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
               , scales = "free_y",
               nrow = 1,
               labeller = labeller(Treatment = c("C" = "Control", 
                                                 "D" = "Daily", 
                                                 "I" = "Intermittent", 
                                                 "E" = "Extreme"))
    )+
    ggokabeito::scale_fill_okabe_ito() + 
    ggokabeito::scale_color_okabe_ito() +
    labs(y = "Biovolume %",
         x = "Experimental day")+
    theme_bw()+
    theme(strip.text = element_text(face="bold")))


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
