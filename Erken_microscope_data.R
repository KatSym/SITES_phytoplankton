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
         Species = case_when(
           # Species == "Ankyra" ~ "Ankyra judayi",
                             Species == "Chlorella" ~ "Eudorina",
                             Species == "Heliozoa" ~ "Golenkinia",
                             Species == "Chlorophyceae" ~ "Sphaerocystis",
                             Species == "Circle" ~ "Centrales",
                             Species == "Comasiella" ~ "Scenedesmus",
                             # Species == "Elakatothrix" ~ "Elakatothrix gelatinosa",
                             Species == "Ochromonas" ~ "Cryptomonas",
                             Species == "Phaester" ~ "Chrysochromulina",
                             Species == "Rhodomonas sp" ~ "Plagioselmis",
                             Species == "Paulschulzsia" ~ "Sphaerocystis",
                             Species == "Stephanodiscus" ~ "Centrales",
                             .default = Species), 
         
         # Species = ifelse(Species == "Rhodomonas sp", "Rhodomonas", Species),
         fun_grp = case_when(Species == "Ankistrodesmus" ~ "IV",
                             Species == "Ankyra" ~ "IV",
                             Species == "Chlamydomonas" ~ "I",
                             Species == "Eudorina" ~ "VII",
                             Species == "Golenkinia" ~ "IV",
                             Species == "Sphaerocystis" ~ "I",
                             Species == "Chroomonas" ~ "I",
                             # Species == "Cyclotella" ~ "VI",
                             Species == "Coelastrum" ~ "IV",
                             Species == "Scenedesmus" ~ "IV",
                             Species == "Cosmarium" ~ "IV",
                             Species == "Cryptomonas" ~ "V",
                             Species == "Elakatothrix" ~ "IV",
                             Species == "Kirchneriella" ~ "IV",
                             Species == "Desmodesmus" ~ "IV",
                             Species == "Gloeocapsa" ~ "VII",
                             Species == "Koliella" ~ "IV",
                             Species == "Monoraphidium" ~ "IV",
                             Species == "Ochromonas" ~ "V",
                             Species == "Oocystis" ~ "IV",
                             Species == "Chrysochromulina" ~ "I",
                             Species == "Merismopedia" ~ "I",
                             Species == "Plagioselmis" ~ "I",
                             Species == "Sphaerocystis" ~ "I",
                             Species == "Quadrigula" ~ "IV",
                             Species == "Snowella" ~ "I",
                             Species == "Stauridium" ~ "IV",
                             Species == "Centrales" ~ "VI",
                             Species == "Tetraedron" ~ "IV")) |> 
  separate(Species, " ", into = c("Genus", "sp"), remove = F) |> 
  select(-sp)

# # THESE ARE CELL VOLUMES
# vols <- read.csv("data/Erken_small_phyto_vol.csv") |> 
#   filter(!Minstorlek %in% c(15, 20, 25)) |> 
#   mutate(biovol_um = Biovolym..mm3.l.*1e9,
#          # in um^3
#          cell_vol = biovol_um/density..celler.l.) |> 
#   group_by(taxon) |> 
#   # mean of 2 sampling in beginning of July for 3 years
#   summarise(cellvol = mean(cell_vol))


# get average number of cells per colony
nomp <- read.csv("data/nomp_subset.csv", header = T, check.names = F) |> 
  select(-`Calculated_volume_\xb5m3/counting_unit`) |> 
  rename(nomp_vol = `Calculated_volume_\xb5m3 (with formula) - NOT IMPORTED, NOT handled by ICES`,
         cells_per_unit = `No_of_cells/counting_unit`) |> 
  group_by(Genus) |>
  summarise_all(mean) |> 
  mutate(cells_per_colony = round(cells_per_unit), .keep = "unused") |> 
  # select(Genus, vol) |> 
  ungroup() |> 
  select(Genus, nomp_vol, cells_per_colony) 




erk20 <- readxl::read_xlsx("data/Phytoplankton 1960-20xx SITES.xlsx", sheet = 1)[,1:14] |> 
  separate(ScientificName, " ", into = c("spec", "sp","var"), remove = F) |> 
  select(-c(var, sp, `Species Flag (SFLAG)`, Details, Phylum, Class, Order), -contains("Depth")) |> 
  filter(spec %in% comcomp$Species 
         # & spec == "Plagioselmis"
         ) |> 
  mutate(Date = as.character(Date),
         month = parse_number(substr(Date, 6,7))) |>
  # filter(month == 7:8) |>
  mutate(cell_vol =(`Biovol.(µm³/l)`/`Density (cells/l)`), .keep = "unused") |> 
  group_by(spec) |> 
  summarise(cell_vol = mean(cell_vol)) |> 
  # add cells per colony
  full_join(nomp, by = join_by("spec" == "Genus")) |> 
  mutate(cells_per_colony = ifelse(is.na(cells_per_colony), 1, cells_per_colony))



#
erk_vol_extra <- readxl::read_xlsx("data/Phytoplankton 1960-20xx SITES.xlsx", sheet = 1)[,1:14] |>
  separate(ScientificName, " ", into = c("spec", "sp","var"), remove = F) |>
  # filter(spec == "Chlorophyceae")
  filter(spec == "Chlorophyceae" & Details == "colony in gell") |> 
  mutate(cell_vol = 1e9*(`Biovol.(µm³/l)`/`Density (cells/l)`), .keep = "unused") |> 
  select(spec, cell_vol) |> 
  summarise_all(mean)
 




small_phyto <- comcomp |> 
  left_join(erk20, by = join_by("Genus" == "spec")) |> 
  mutate(
    # make average cells per colony an even number for cyanos
    cells_per_colony = case_when(cells_per_colony == 21 ~ 22,
                                 cells_per_colony == 57 ~ 58,
                                 # Species == "Ankistrodesmus" ~ 4,
                                  Species == "Kirchneriella" ~ mean(4:32),
                                  Species == "Quadrigula" ~ mean(2:8),
                                  Species == "Gloeocapsa" ~ 1,
                                 .default = cells_per_colony),
    # or width
    diameter = case_when(Species == "Golenkinia" ~ mean(6:26), # Ellis & Machlis 1968 https://bsapubs.onlinelibrary.wiley.com/doi/abs/10.1002/j.1537-2197.1968.tb07416.x
                         Species == "Ankistrodesmus" ~ 2.5, # AlgeaBase
                         Species == "Chroomonas" ~ mean(6:10), # Hoef-Emden 2018 https://www.sciencedirect.com/science/article/abs/pii/S1434461018300294
                         Species == "Gloeocapsa" ~ mean(2.5:14), # https://www.sciencedirect.com/science/article/abs/pii/B9780127415505500040
                         Species == "Kirchneriella" ~ mean(3:6), #AglaeBase or ?  http://protist.i.hosei.ac.jp/pdb/images/chlorophyta/kirchneriella/lunaris/lunaris8.html
                         Species == "Quadrigula" ~ mean(1:8), # AlgeaBase
                         Species == "Tetraedron" ~ mean(2:4), # Nordic 
                         ), 
    
    d2 = case_when(Species == "Chroomonas" ~ mean(8:12),
    ),
    length = case_when(Species == "Ankistrodesmus" ~ mean(15:100),
                       Species == "Chroomonas" ~ mean(12:18),
                       Species == "Kirchneriella" ~ mean(8:12),
                       Species == "Quadrigula" ~ mean(7:45),
                       Species == "Tetraedron" ~ mean(8:10)
                       ),
    
    # volumes seem too big
    calc_vol = case_when(Species == "Golenkinia" ~ (1/6)*pi*(diameter^3), 
                         Species == "Ankistrodesmus" ~ (2/15)*pi*(diameter^2)*length, # spindle
                         Species == "Chroomonas" ~ (1/6)*pi*length*diameter*d2, # flattened ellipsoid
                         Species == "Gloeocapsa" ~ (1/6)*pi*(diameter^3),
                         Species == "Kirchneriella" ~ (1/6)*pi*(diameter^2)*length, # rotational ellipsoid
                         Species == "Quadrigula" ~ (2/15)*pi*(diameter^2)*length,
                         Species == "Tetraedron" ~ diameter*length,
                         ),
    cell_vol = ifelse(Species == "Sphaerocystis", erk_vol_extra$cell_vol, cell_vol),
    indiv_vol = coalesce(cell_vol, calc_vol),
    biovolume = abund*indiv_vol*cells_per_colony,
    # .keep = "unused"
    )  
  select(-c(10:17))
  

write.csv(small_phyto, "data/WRONG_small_phyto2.csv", row.names = F, quote = F)

