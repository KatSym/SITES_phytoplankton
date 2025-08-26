library(tidyverse)

comcomp <- read.csv("data/comm_comp_GAL.csv", header = T,
                    skip = 5) %>% 
  # some "data cleaning" to match my other data structures
  mutate(samplDay = as.numeric(substring(Sample_date, 3, 3)),
         ExpDay = case_when(samplDay == 0 ~ 0,
                            samplDay == 3 ~ 12,
                            samplDay == 5 ~ 20,
                            samplDay == 9 ~ 36), .keep = "unused") %>% 
  separate(Sample_ID, c("mesocosm", "blah"), "_") %>% 
  mutate(Treatment = substring(blah, 1, 1)) %>% 
  # abundance calculation - cells/mL
  mutate(abund = ((Counts*Utt_Area)/(No_LF*fov))/sed_volume) %>% # in cells/mL
  select(ExpDay, Treatment, mesocosm, Species, abund) %>% 
  filter(Species %in% c("Ankistrodesmus", "Ankyra", "Chlamydomonas", "Chlorella", "Chlorophyceae", "Chroomonas", 
                        "Circle", "Coelastrum", "Comasiella", "Cosmarium", "Cryptomonas", 
                        "Elakatothrix", "Kirchneriella", "Desmodesmus", "Gloeocapsa", "Koliella", "Heliozoa",
                        "Monoraphidium", "Ochromonas", "Oocystis", "Phaester", "Merismopedia", 
                        "Quadrigula", "Rhodomonas sp", "Paulschulzsia", "Quadrigula", "Snowella", "Stauridium", 
                        "Stephanodiscus", "Tetraedron")) %>% 
  mutate(size_fract = "small",
         # add label column to match with the flowcam data
         label = case_when(Species == "Ankistrodesmus" ~ "ChlAnk000_9516573",
                           Species == "Ankyra" ~ "ChlAnkjud_2642105",
                           Species == "Chlamydomonas" ~ "ChlChl000_5271004",
                           Species == "Chlorella" ~ "ChlEud000_0000000",
                           Species == "Heliozoa" ~ "ChlGolrad_2641100",
                           Species == "Chlorophyceae" ~ "ChlSph000_2639562",
                           Species == "Chroomonas" ~ "CryChr000_3202572",
                           Species == "Circle" ~ "BacCyc000_3193095",
                           Species == "Coelastrum" ~ "ChlCoe000_7749802",
                           Species == "Comasiella" ~ "ChlScearc_2640846",
                           Species == "Cosmarium" ~ "ZygCos000_2648709",
                           Species == "Cryptomonas" ~ "CryCry000_3202413",
                           Species == "Elakatothrix" ~ "KleElagel_2641001",
                           Species == "Kirchneriella" ~ "ChlKir000_2641223",
                           Species == "Desmodesmus" ~ "ChlDes000_2652440",
                           Species == "Gloeocapsa" ~ "CyaGlo000_3217637",
                           Species == "Koliella" ~ "TreKollon_2638730",
                           Species == "Monoraphidium" ~ "ChlMon000_0000000",
                           Species == "Ochromonas" ~ "CryCry000_3202413",
                           Species == "Oocystis" ~ "TreOoc000_2641454",
                           Species == "Phaester" ~ "PryChr000_3202213",
                           Species == "Merismopedia" ~ "CyaMer000_0000000",
                           Species == "Rhodomonas sp" ~ "CryRho000_3202546",
                           Species == "Paulschulzsia" ~ "ChlSph000_2639562",
                           Species == "Quadrigula" ~ "ChlQua000_2641259",
                           Species == "Snowella" ~ "CyaSno000_3217087",
                           Species == "Stauridium" ~ "ChlSta000_4339881",
                           Species == "Stephanodiscus" ~ "BacSte000_3193250",
                           Species == "Tetraedron" ~ "ChlTet000_2672378"),
         
         # correct species to match with other datasets
         Species = case_when(Species == "Ankyra" ~ "Ankyra judayi",
                             Species == "Chlorella" ~ "Eudorina",
                             Species == "Heliozoa" ~ "Golenkinia radiata",
                             Species == "Chlorophyceae" ~ "Sphaerocystis",
                             Species == "Circle" ~ "Cyclotella",
                             Species == "Comasiella" ~ "Scenedesmus arcuatus",
                             Species == "Elakatothrix" ~ "Elakatothrix gelatinosa",
                             Species == "Ochromonas" ~ "Cryptomonas",
                             Species == "Phaester" ~ "Chrysochromulina",
                             Species == "Rhodomonas sp" ~ "Rhodomonas",
                             Species == "Paulschulzsia" ~ "Sphaerocystis",
                             .default = Species),
         
         fun_grp = case_when(Species == "Ankistrodesmus" ~ "IV",
                             Species == "Ankyra judayi" ~ "IV",
                             Species == "Chlamydomonas" ~ "I",
                             Species == "Eudorina" ~ "VII",
                             Species == "Golenkinia radiata" ~ "IV",
                             Species == "Sphaerocystis" ~ "I",
                             Species == "Chroomonas" ~ "I",
                             Species == "Cyclotella" ~ "VI",
                             Species == "Coelastrum" ~ "IV",
                             Species == "Scenedesmus arcuatus" ~ "IV",
                             Species == "Cosmarium" ~ "IV",
                             Species == "Cryptomonas" ~ "V",
                             Species == "Elakatothrix gelatinosa" ~ "IV",
                             Species == "Kirchneriella" ~ "IV",
                             Species == "Desmodesmus" ~ "IV",
                             Species == "Gloeocapsa" ~ "VII",
                             Species == "Koliella" ~ "IV",
                             Species == "Monoraphidium" ~ "IV",
                             Species == "Ochromonas" ~ "V",
                             Species == "Oocystis" ~ "IV",
                             Species == "Chrysochromulina" ~ "I",
                             Species == "Merismopedia" ~ "I",
                             Species == "Rhodomonas" ~ "I",
                             Species == "Paulschulzsia" ~ "I",
                             Species == "Quadrigula" ~ "IV",
                             Species == "Snowella" ~ "I",
                             Species == "Stauridium" ~ "IV",
                             Species == "Stephanodiscus" ~ "VI",
                             Species == "Tetraedron" ~ "IV")) |> 
  separate(Species, " ", into = c("Genus", "sp"), remove = F) |> 
  select(-sp)

# THESE ARE CELL VOLUMES
vols <- read.csv("data/Erken_small_phyto_vol.csv") |> 
  filter(!Minstorlek %in% c(15, 20, 25)) |> 
  mutate(biovol_um = Biovolym..mm3.l.*1e9,
         # in um^3
         cell_vol = biovol_um/density..celler.l.) |> 
  group_by(taxon) |> 
  # mean of 2 sampling in beginning of July for 3 years
  summarise(cellvol = mean(cell_vol))


nomp_vols <- read.csv("data/nomp_subset.csv", header = T, check.names = F) |> 
  select(-`Calculated_volume_\xb5m3/counting_unit`) |> 
  rename(vol = `Calculated_volume_\xb5m3 (with formula) - NOT IMPORTED, NOT handled by ICES`) |> 
  group_by(Genus) |>
  summarise_all(mean) |> 
  select(Genus, vol) |> 
  ungroup()
  # readxl::read_xlsx("data/bvol_nomp_version_2024.xlsx", 
  #                              sheet = 1, range = "A1:AM3847") |> 
  # discard(~all(is.na(.))) |> 
  # filter(Genus %in% comcomp$Genus,
  #        Species %in% comcomp$Species) |> 
  # select(-contains(c("AphiaID", "Length", "Height", "Carbon", "Diameter", "Width", "Synonyms" )),
  #        -c(List, `HELCOM area`, Division, Class, Order, Author, `WORMS Rank`, Trophy)) |> 
  # filter(!(Genus %in% c("Merismopedia", "Snowella", "Elakatothrix") & Unit == "cell")) |> 
  # group_by(Species) |> 
  # summarise_all(mean)

small_phyto <- comcomp |> 
  left_join(nomp_vols, by = "Genus") |> 
  mutate(biovolume = abund*vol)

write.csv(small_phyto, "data/WRONG_small_phyto.csv", row.names = F, quote = F)
