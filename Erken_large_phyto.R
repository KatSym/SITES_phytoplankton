library(readxl)

# gloeotrichia counts
gloe <- read_xlsx("data/gloeotrichia_colonies_etc.xlsx",
                  sheet = 1,
                  range = "A1:G107") |>
  rename(day = `day of experiment`, 
         Treatment = treatment,
         sample_vol_ml = `sample volume (ml)`) |> 
  select(day, mesocosm, Treatment, sample_vol_ml, gloeotrichia) |> 
  filter(Treatment != "LE") |> 
  mutate(mesocosm = as.numeric(mesocosm),
         day = as.numeric(day))


# fragilaria counts
frag <- read_xlsx("data/gloeotrichia_colonies_etc.xlsx",
                  sheet = 2,
                  range = "A2:J91") |> 
  filter(treatment != "LE") |> 
  rename(day = `day of experiment`, 
         Treatment = treatment,
         sample_vol_ml = `sample volume (ml)`,
         diatoms_chamber = `diatoms (in chamber)`,
         diatoms_snake = `diatoms (in snake)`) |> 
  select(day, mesocosm, Treatment, sample_vol_ml, diatoms_chamber, diatoms_snake) |> 
  mutate(day = as.numeric(day),
         mesocosm = as.numeric(mesocosm))

# take counts for day 36 from measurement data
frag_d36 <- read.csv("data/fragilaria_measurements.csv", header = T) |> 
  group_by(day, Treatment, mesocosm) |> 
  summarise(diatoms_chamber = n()) |> 
  ungroup() |> 
  filter(Treatment != "ERK") |> 
  filter(day == 36)
# correct values from the notes
frag_d36$diatoms_chamber[frag_d36$day == 36 & frag_d36$mesocosm == 8] <- 69
frag_d36$diatoms_chamber[frag_d36$day == 36 & frag_d36$mesocosm == 10] <- 845
frag_d36$diatoms_chamber[frag_d36$day == 36 & frag_d36$mesocosm == 15] <- 102
frag_d36$diatoms_chamber[frag_d36$day == 36 & frag_d36$mesocosm == 16] <- 132

# combine the two
frag_count <- full_join(frag, frag_d36, by = c("day", "Treatment", "mesocosm")) |> 
  mutate(diatoms_chamber = coalesce(diatoms_chamber.y, diatoms_chamber.x), 
         # get count from both methods
         frag_count = coalesce(diatoms_chamber, diatoms_snake),
         .keep = "unused"
         ) |> 
  select(-diatoms_chamber)

frag_meas <- read.csv("data/fragilaria_measurements.csv", header = T) |> 
  select(-Key, -File) |> 
  # need to give a reference for that
  mutate(height = 2.5,
         frag_vol = width*length*height) |> 
  group_by(day, Treatment, mesocosm) |> 
  # get mean volume of fragilaria colonies per mesocosm
  summarise(frag_vol = mean(frag_vol)) |> 
  ungroup() |> 
  add_row(day = 36, Treatment = "I", mesocosm = 3, frag_vol = 0) |> 
  add_row(day = 28, Treatment = "I", mesocosm = 14, frag_vol = NA)

gloe_diam <- read.csv("data/gloeotrichia_diameter_Karlsson.csv", header = T) |> 
  select(-contains("X")) |> 
  mutate(gloe_vol = 4/3*pi*(Diameter_um3/2)^3,
         day = case_when(Date == "01/07/2001" ~ 0,
                         Date == "09/07/2001" ~ 4,
                         Date == "17/07/2001" ~ 12,
                         Date == "24/07/2001" ~ 20,
                         Date == "01/08/2001" ~ 28,
                         Date == "16/08/2001" ~ 36,
                         Date == "01/07/2000" ~ 0,
                         Date == "09/07/2000" ~ 4,
                         Date == "17/07/2000" ~ 12,
                         Date == "24/07/2000" ~ 20,
                         Date == "01/08/2000" ~ 28,
                         Date == "08/08/2000" ~ 36)) |> 
  group_by(day) |> 
  summarise(gloe_vol = mean(gloe_vol)) |> 
  ungroup() |> 
  summarise_all(mean)
  


erk_large_phyto <- gloe |> 
  left_join(frag_count, by = c("day", "Treatment", "mesocosm")) |> 
  # left_join(volv, by = c("day", "Treatment", "mesocosm", "sample_vol_ml")) |> 
  left_join(frag_meas, by = c("day", "Treatment", "mesocosm")) |> 
  # left_join(gloe_diam, by = "day") |> 
  mutate(gloe_vol = gloe_diam$gloe_vol) |> 
  rename(gloe_count = gloeotrichia,
         sample_vol_ml = sample_vol_ml.x) |> 
  # calculate abundance (colonies per mL) and biovolume (um3 per mL)
  mutate(frag_abund = frag_count/(sample_vol_ml),
         frag_biovol = frag_abund * frag_vol,
         gloe_abund = gloe_count/(sample_vol_ml),
         gloe_biovol = gloe_abund * gloe_vol,
         .keep = "unused") |> 
  select(-sample_vol_ml.y)



write.csv(erk_large_phyto, "data/Erken_large_phyto.csv", row.names = F, quote = F)























frag_LW <- read.csv("data/fragilaria_measurements.csv", header = T) |> 
  select(-Key, -File, -mesocosm) |> 
  filter(Treatment != "ERK") |> 
  group_by(day, Treatment) |> 
  summarise_all(mean) |> 
  ungroup()
