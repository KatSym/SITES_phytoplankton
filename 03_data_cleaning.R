library(tidyverse)

#  ERKEN -----
## load required data ----
erken <- readRDS("./data/Erken_phyto_all_20250626.rds") 

metadata.erk <- read.csv("./data/Metadata_Erken.csv", 
                     header=T, sep=",") %>% 
  # filter(ID %in% unique(data$Name)) %>% 
  select(ID, vol_FC..ml., concvol..ml., concfactor) %>% 
  rename(volume_run = vol_FC..ml.,
         concentr_vol = concvol..ml.)


# transform data into something usable 
erk.dat <- erken %>% 
  # remove files that should not be used (from the notes book)
  filter(!Name %in% c("E10_20220804_4x_1000in125_02a",
                      "E11_20220812_4x_500in125_03.csv",
                      "E11_20220812_4x_500in125_03a.csv", 
                      "E11_20220812_4x_500in125_03b.csv")) %>% 
  mutate(AbdDiameter = as.numeric(AbdDiameter),
         Length = as.numeric(Length),
         biovolMS = as.numeric(biovolMS),
         surfacearea = as.numeric(surfacearea),
         ProbabilityScore = as.numeric(ProbabilityScore),
         label = ifelse(LabelTrue %in% c(NA , "", "EVERYTHING"), LabelPredicted, LabelTrue),
         Name = case_when(Name == "E014_20220804_4x_500in125_01" ~ "E14_20220804_4x_500in125_01",
                          Name == "E014_20220804_4x_500in125_02" ~ "E14_20220804_4x_500in125_02",
                          .default = Name),
         .keep = "unused"
         ) %>%
  # filter out unwanted classes
  filter(!label %in% c(     "",
                            "detritus", 
                            "CyaGlo000_filaments",
                            "CyaGloech_0000000",
                            "Mesozoo_0000000"
                          )) %>% 
  # remove '_updated' from sample name (don't know how it's there)
  mutate(Name = str_remove(Name, "_updated")) %>% 

  separate(Name, into = c("mes", "day", "magn", "conc", "replicate"), 
           sep = "_", remove = F) %>% 
  select(-Date, -magn, -conc) %>% 
  mutate(Date = as.Date(day, "%Y%m%d"),
    mesocosm = as.numeric(gsub("[A-Za-z]", "", mes)),
    Treatment = case_when(mesocosm == 1 | mesocosm == 7 | mesocosm == 10 | mesocosm == 16 ~ "C",
                          mesocosm == 2 | mesocosm == 8 | mesocosm == 11 | mesocosm == 13 ~ "D",
                          mesocosm == 3 | mesocosm == 5 | mesocosm == 12 | mesocosm == 14 ~ "I",
                          mesocosm == 4 | mesocosm == 6 | mesocosm == 9 | mesocosm == 15 ~ "E",
                          .default = "ERK"),
    ExpDay = case_when(Date == "2022-07-07" ~ 0,
                       Date == "2022-07-11" ~ 4,
                       Date == "2022-07-15" ~ 8,
                       Date == "2022-07-19" ~ 12,
                       Date == "2022-07-23" ~ 16,
                       Date == "2022-07-27" ~ 20,
                       Date == "2022-07-31" ~ 24,
                       Date == "2022-08-04" ~ 28,
                       Date == "2022-08-08" ~ 32,
                       Date == "2022-08-12" ~ 36),
    Treatment = factor(Treatment, levels = c("C", "D", "I", "E", "ERK")),
    ProbabilityScore = ifelse(is.na(ProbabilityScore) == TRUE, 1, ProbabilityScore), # or maybe make it 1
    # correct wrong class
    label = ifelse(label == "RapGonsem_3194431", "Cil_other", label)
          ) %>% 
  # add metadata to calculate abundance
  left_join(., metadata.erk, by = join_by(Name == ID)) %>% 
  select(-c(mes, day, concentr_vol))


# BOLMEN -----
## load data ----
bolmen <- readRDS("./data/Bolmen_FCphyto_20250630-biovol.rds") 

metadata.bol <- read.csv("./data/Metadata_Bolmen.csv", 
                         header=T, sep=",") %>% 
  # filter(ID %in% unique(data$Name)) %>% 
  select(ID, vol_FC..ml., concvol..ml., concfactor) %>% 
  rename(volume_run = vol_FC..ml.,
         concentr_vol = concvol..ml.)
  



bol.dat <- bolmen %>% 
  mutate(biovolume = as.numeric(biovolume),
         surfacearea = as.numeric(surfacearea),
         ProbabilityScore = as.numeric(ProbabilityScore),
         label = ifelse(LabelTrue %in% c(NA , "", "EVERYTHING"), LabelPredicted, LabelTrue),
         .keep = "unused") %>% 
  # filter out unwanted classes
  filter(!label %in% c("detritus", 
                            "CyaGloech_0000000"
  )) %>% 
  separate(Name, into = c("mes", "day", "magn", "conc", "replicate"), 
           sep = "_", remove = F) %>% 
  select( -magn, -conc) %>% 
  mutate(Date = as.Date(day, "%Y%m%d"),
         mesocosm = as.numeric(gsub("[A-Za-z]", "", mes)),
         Treatment = case_when(mesocosm == 1 | mesocosm == 7 | mesocosm == 10 | mesocosm == 16 ~ "C",
                               mesocosm == 2 | mesocosm == 8 | mesocosm == 11 | mesocosm == 13 ~ "D",
                               mesocosm == 3 | mesocosm == 5 | mesocosm == 12 | mesocosm == 14 ~ "I",
                               mesocosm == 4 | mesocosm == 6 | mesocosm == 9 | mesocosm == 15 ~ "E",
                               .default = "BOL"),
         ExpDay = case_when(Date == "2022-07-07" ~ 0,
                            Date == "2022-07-11" ~ 4,
                            Date == "2022-07-15" ~ 8,
                            Date == "2022-07-19" ~ 12,
                            Date == "2022-07-23" ~ 16,
                            Date == "2022-07-27" ~ 20,
                            Date == "2022-07-31" ~ 24,
                            Date == "2022-08-04" ~ 28,
                            Date == "2022-08-08" ~ 32,
                            Date == "2022-08-12" ~ 36),
         Treatment = factor(Treatment, levels = c("C", "D", "I", "E", "BOL")),
         ProbabilityScore = ifelse(is.na(ProbabilityScore) == TRUE, 1, ProbabilityScore), # or maybe make it 1
           ) %>% 
  # add metadata to calculate abundance
  left_join(., metadata.bol, by = join_by(Name == ID)) %>% 
  select(-c(mes, day, concentr_vol))

save(erk.dat, bol.dat, file = "./data/sites_FC_phyto.RData")
