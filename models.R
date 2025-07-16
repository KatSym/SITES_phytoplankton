library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)

# load data
# --- Erken ---
load("./data/sites_FC_phyto.RData")
functional.erk <- read.csv("./data/Erken_functional.csv", header = T, sep = ",")

edat <- erk.dat %>% 
  # drop the lake
  filter(Treatment != "ERK",
         ProbabilityScore > 0.5
         ) %>% 
  mutate(Treatment = droplevels(Treatment),
         vol.offset = volume_run*concfactor,
         label = ifelse(label == "ConSta000_2647648", "ZygSta000_2647648", label),
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
            biovol = sum(mean_biovol),
            vol.offset = mean(vol.offset)) %>% 
  mutate(gr_biovol = count*biovol/vol.offset) |> 
  ungroup()


edat.tr <- erk.dat %>% 
  # drop the lake
  filter(Treatment != "ERK",
         ProbabilityScore > 0.5
  ) %>% 
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
         label = ifelse(label == "ConSta000_2647648", "ZygSta000_2647648", label),
         fun_grp = case_when(label == "FILAMENTS" ~ "III",
                             label == "CyaPla000_3218374" ~ "III",
                             class == "Bac" ~ "VI",
                             .default = KRUK_MBFG)) %>% 
  select(-taxon, -KRUK_MBFG, -Name) %>% 
  
  drop_na()  %>% 
  full_join(., sp_trt, by = join_by(label == Code))
 

ptm <- proc.time()
m1.e <- brm(bf(count ~ Treatment * ExpDay +
                 (ExpDay|mesocosm) +
                 offset(log(vol.offset))),
            family = poisson(),
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
pp_check(m1.e, ndraws = 100)
plot(conditional_effects(m1.e, effects = "ExpDay:Treatment"), points = TRUE)

ptm <- proc.time()
m2.e <- brm(bf(count ~ Treatment + ExpDay + fun_grp +
                 Treatment : ExpDay +
                 ExpDay : fun_grp +
                 fun_grp : Treatment +
                 (ExpDay|mesocosm) +
                 offset(log(vol.offset))),
            family = negbinomial(),
            chains = 4,
            iter = 4000,
            cores = 4,
            control = list(adapt_delta = 0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = edat,
            file = "models/Erk_FunctGrp-count-offset_20250708.m",
            file_refit = "on_change"
) 
proc.time() - ptm
pp_check(m2.e, ndraws = 100) + scale_x_log10()
summary(m2.e)
plot(conditional_effects(m2.e, offset = T), ask = FALSE, points = TRUE)
plot(conditional_effects(m2.e, effects = "ExpDay:Treatment"), points = TRUE)


trt.cols <- c(`C`= "#000000", #black - C
              `D`= "#0a8754", #g - D
              `I`= "#4472ca", #blu - I
              `E`= "#e84855") 

fil.cols <- c(`C`= "#00000030", #black - C
              `D`= "#0a875430", #g - D
              `I`= "#4472ca30", #blu - I
              `E`= "#e8485530") 
  


m2e.pred <- edat |> 
  data_grid(Treatment = unique(edat$Treatment),
            fun_grp = unique(edat$fun_grp),
            ExpDay = unique(edat$ExpDay),
            vol.offset = 1) |> 
  add_epred_draws(m2.e, re_formula = NA) 


m2e.pred|> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 

  ggplot(aes(x = .epred, y = Treatment)) +
  stat_pointinterval(point_interval = "mean_hdi",
                     .width =  0.95, linewidth  = .4, 
                     # interval_colour = ExpDay, 
                     # point_color = ExpDay
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  facet_wrap(~ fun_grp, scales = "free")+
  ylab("contrast") +
  xlab("difference") +
  theme_minimal()

(plot <- m2e.pred %>%
  ggplot(aes(x = ExpDay,
             y = .epred,
             colour = Treatment)
  ) +
  stat_lineribbon(aes(y = (.epred),
                      fill = Treatment),
                  .width = c(0.95)
                  # alpha = .3
                  # position = position_dodge(.5),
                  # linewidth = 2
  )+
  geom_point(edat, mapping = aes(x = ExpDay,
                                 y = count/vol.offset,
                                 colour = Treatment),
             alpha = .35,
             position = position_jitterdodge(dodge.width = .5),
  ) +
  facet_wrap(~ fun_grp
             , scales = "free_y"
  )+
  scale_color_manual(values = trt.cols) +
  scale_fill_manual(values = fil.cols) )
  labs(title = "Heterotroph",
       y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
       x = NULL)
#   # theme(axis.title.y = ggtext::element_markdown(),
#   #       axis.title.x = element_text(size = 15),
#   #       plot.title = element_text(size=16, face="italic")
#   # ) 

draws <- as_draws_df(m2.e) %>% 
  select(-contains("r_mesocosm")) %>% 
  pivot_longer(cols = contains("b_"), names_to = "interaction", values_to = "intrcepts") %>% 
  ggplot(aes(x = intrcepts, y = interaction)) +
  stat_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dotted")
  


ptm <- proc.time()
m3.e <- brm(bf(gr_biovol ~ Treatment + ExpDay + fun_grp +
                 Treatment : ExpDay +
                 ExpDay : fun_grp +
                 fun_grp : Treatment +
                 (ExpDay|mesocosm)),
            family = lognormal(),
            chains = 4,
            iter = 4000,
            cores = 4,
            control = list(adapt_delta = 0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = edat,
            file = "models/Erk_FunctGrp-biovol-offset_20250709.m",
            file_refit = "on_change"
) 
proc.time() - ptm
pp_check(m3.e, ndraws = 100) + scale_x_log10()
summary(m3.e)


m3e.pred <- edat |> data_grid(Treatment = unique(edat$Treatment),
                              fun_grp = unique(edat$fun_grp),
                              ExpDay = unique(edat$ExpDay)
) |> 
  add_epred_draws(m3.e, re_formula = NA)  


m3e.pred|> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  mean_qi() |> 
  # filter(ExpDay == 12) |> 
  ggplot(aes(x = .epred, y = Treatment)) +
  stat_pointinterval(point_interval = "mean_hdi",
                     .width =  0.95, linewidth  = .4, 
                     # interval_colour = ExpDay, 
                     # point_color = ExpDay
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  facet_wrap(~ fun_grp, scales = "free")+
  ylab("contrast") +
  xlab("difference") +
  theme_minimal()
  

(plot <- m3e.pred %>%
    ggplot(aes(x = ExpDay,
               y = .epred,
               colour = Treatment)
    ) +
    stat_lineribbon(aes(y = (.epred),
                        fill = Treatment),
                    .width = c(0.95)
                    # alpha = .3
                    # position = position_dodge(.5),
                    # linewidth = 2
    )+
    geom_point(edat, mapping = aes(x = ExpDay,
                                   y = gr_biovol,
                                   colour = Treatment),
               alpha = .35,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    facet_wrap(~ fun_grp
               , scales = "free_y"
    )+
    scale_color_manual(values = trt.cols) +
    scale_fill_manual(values = fil.cols) )
labs(title = "Heterotroph",
     y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
     x = NULL)




# --- Bolmen ----
ptm <- proc.time()
m1.b <- brm(bf(count ~ Treatment + ExpDay + fun_grp +
                 Treatment : ExpDay +
                 ExpDay : fun_grp +
                 fun_grp : Treatment +
                 (ExpDay|mesocosm) +
                 offset(log(vol.offset))),
            family = negbinomial(),
            chains = 4,
            iter = 4000,
            cores = 4,
            control = list(adapt_delta = 0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = bolmen.f,
            file = "models/Bol_FunctGrp-count-offset_20250716.m",
            file_refit = "on_change"
) 
proc.time() - ptm

pp_check(m1.b, ndraws = 100) + scale_x_log10()
summary(m1.b)


m1b.pred <- bolmen.f |> 
  data_grid(Treatment = unique(edat$Treatment),
            fun_grp = unique(edat$fun_grp),
            ExpDay = unique(edat$ExpDay),
            vol.offset = 1
  ) |> 
  add_epred_draws(m1.b, re_formula = NA) 


m1b.pred |> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  # filter(ExpDay == 12) |> 
  mean_qi() |> 
  ggplot(aes(x = .epred, y = Treatment)) +
  stat_pointinterval(point_interval = "mean_hdi",
                     .width =  0.95, linewidth  = .4, 
                     # interval_colour = ExpDay, 
                     # point_color = ExpDay
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  facet_wrap(~ fun_grp, scales = "free")+
  ylab("contrast") +
  xlab("difference") +
  theme_minimal()

(plot <- m1b.pred %>%
    ggplot(aes(x = ExpDay,
               y = .epred,
               colour = Treatment)
    ) +
    stat_lineribbon(aes(y = (.epred),
                        fill = Treatment),
                    .width = c(0.95)
                    # alpha = .3
                    # position = position_dodge(.5),
                    # linewidth = 2
    )+
    geom_point(bolmen.f, mapping = aes(x = ExpDay,
                                   y = count/vol.offset,
                                   colour = Treatment),
               alpha = .35,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    facet_wrap(~ fun_grp
               , scales = "free_y"
    )+
    scale_color_manual(values = trt.cols) +
    scale_fill_manual(values = fil.cols) )
labs(title = "Heterotroph",
     y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
     x = NULL)


ptm <- proc.time()
m2.b <- brm(bf(gr_biovol ~ Treatment + ExpDay + fun_grp +
                 Treatment : ExpDay +
                 ExpDay : fun_grp +
                 fun_grp : Treatment +
                 (ExpDay|mesocosm) ),
            family = lognormal(),
            chains = 4,
            iter = 4000,
            cores = 4,
            control = list(adapt_delta = 0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = bolmen.f,
            file = "models/Bol_FunctGrp-biov-offset_20250716.m",
            # file_refit = "on_change"
) 
proc.time() - ptm

pp_check(m2.b, ndraws = 100) + scale_x_log10()
summary(m2.b)

m2b.pred <- bolmen.f |> 
  data_grid(Treatment = unique(edat$Treatment),
            fun_grp = unique(edat$fun_grp),
            ExpDay = unique(edat$ExpDay)
  ) |> 
  add_epred_draws(m2.b, re_formula = NA) 


m2b.pred |> 
  compare_levels(.epred, by = Treatment,
                 comparison = pairwise) |> 
  # filter(ExpDay == 12) |> 
  mean_qi() |> 
  ggplot(aes(x = .epred, y = Treatment)) +
  stat_pointinterval(point_interval = "mean_hdi",
                     .width =  0.95, linewidth  = .4, 
                     # interval_colour = ExpDay, 
                     # point_color = ExpDay
  ) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  facet_wrap(~ fun_grp, scales = "free")+
  ylab("contrast") +
  xlab("difference") +
  theme_minimal()


(plot <- m2b.pred %>%
    ggplot(aes(x = ExpDay,
               y = .epred,
               colour = Treatment)
    ) +
    stat_lineribbon(aes(y = (.epred),
                        fill = Treatment),
                    .width = c(0.95)
                    # alpha = .3
                    # position = position_dodge(.5),
                    # linewidth = 2
    )+
    geom_point(bolmen.f, mapping = aes(x = ExpDay,
                                       y = gr_biovol,
                                       colour = Treatment),
               alpha = .35,
               position = position_jitterdodge(dodge.width = .5),
    ) +
    facet_wrap(~ fun_grp
               , scales = "free_y"
    )+
    scale_color_manual(values = trt.cols) +
    scale_fill_manual(values = fil.cols) )
labs(title = "Heterotroph",
     y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
     x = NULL)

# traits models ----

ptm <- proc.time()
m2tr.e <- brm(bf(count ~ Treatment + ExpDay + Paff + Iaff + mu +
                 Treatment : ExpDay +
                   Paff : Treatment +
                   Iaff : Treatment +
                   mu : Treatment +
                 (ExpDay|mesocosm) +
                 offset(log(vol.offset))),
            family = negbinomial(),
            chains = 4,
            iter = 4000,
            cores = 4,
            control = list(adapt_delta = 0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = edat.tr,
            file = "models/Erk_traits-count-offset_20250710.m",
            file_refit = "on_change"
) 
proc.time() - ptm
pp_check(m2tr.e, ndraws = 100) + scale_x_log10()
summary(m2tr.e)
