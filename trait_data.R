library(readxl)
library(tidyverse)

tr1 <- read_xlsx("./data/traits.xlsx",
                 sheet = 1,  
                 range = "A2:D95",
                 col_names = T) %>% 
  group_by(Species) %>% 
  summarise(max_drowth = mean(`Max. growth`))

tr2 <- read_xlsx("./data/traits.xlsx",
             sheet = 2,  
             range = "A1:B120",
             col_names = T) %>% 
  group_by(Species) %>% 
  summarise(predation = mean(`Susceptibility to predation (-)`))

tr3 <- read_xlsx("./data/traits.xlsx",
             sheet = 3,  
             range = "A2:B29",
             col_names = T) %>% 
  group_by(Species) %>% 
  summarise(Paff = mean(`Phosphate affinity`))

tr4 <- read_xlsx("./data/traits.xlsx",
             sheet = 4,  
             range = "A1:C26",
             col_names = T)

tr5 <- read_xlsx("./data/traits.xlsx",
                 sheet = 5,  
                 range = "A1:F42",
                 col_names = T) %>% 
  select(-α, -Group, -Size)

trt <- tr1 %>% full_join(., tr2, by = "Species") %>% 
  full_join(., tr3, by = "Species") %>% 
  full_join(., tr4, by = "Species") %>% 
  full_join(., tr5, by = "Species") %>% 
  mutate(Psaff = as.numeric(Psaff),
         growth_max = coalesce(max_drowth, µmax.x, µmax.y)) %>% 
  separate(Species, into = c("genus", "sp"), sep = " ", remove = F) %>% 
  select(-sp) %>% 
  filter(!Species %in% c("Dinobryon sertularia", "Euglena gracilis", "Gomphonema truncatum", 
                         "Mallomonas cratis", "Monoraphidium griffithii", "Mougeotia thylespora",
                         "Nitzschia actinastroides", "Phormidium mucicola", "Rhodomonas lens",
                         "Scenedesmus abundans", "Scenedesmus acutus", "Scenedesmus crassus",                        
                         "Scenedesmus dimorphus", "Scenedesmus protuberans", "Synedra radians",
                         "Synedra rumpens"),
         genus != c(NA, "Stephanodiscus")) %>% 
  group_by(genus) %>% 
  summarise_all(mean)

trt_pr <- read.csv("./data/traits_prager.csv", header = T) %>%   
  separate(species, into = c("genus", "sp"), sep = " ", remove = F) %>% 
  select(-sp)

sp_cd <- read.csv("./data/species-codes.csv", header =T) %>% 
  separate(Species, into = c("genus", "sp"), sep = " ", remove = F) %>% 
  select(-sp) %>% 
  mutate(genus = case_when(Species == "Dolichospermum sp." ~ "Dolichospermum2",
                           .default = genus)) %>% 
  filter(!genus %in% c("Amoebida", "Ciliate", "Coleps", "Conochilus", "Euchlanis", "Heliozoa", "Kellicottia",
                          "Keratella", "Mesozooplankton", "Nauplii", "Polyarthra", "Trichocerca"))

trt2 <- trt_pr %>% 
  filter(genus %in% sp_cd$genus,
         !species %in% c("Fragilaria sp.", "Fragilaria famelica", "Fragilaria radians", 
                         "Fragilaria rumpens", "Fragilaria tenera")) %>% 
  mutate(genus = case_when(species == "Dolichospermum circinale" |
                             species == "Dolichospermum flosaquae" |
                             species == "Dolichospermum planctonicum" ~ "Dolichospermum2",
                           .default = genus)) %>% 
  group_by(genus) %>% 
  summarise_all(mean)



sp_trt <- sp_cd %>% 
  full_join(., trt2, by = "genus") %>% 
  mutate(genus = case_when(Species == "Closterium aciculare" ~ "Closterium2",
                             Species == "Closterium sp." ~ "Closterium4",
                           .default = genus)) %>% 
  group_by(genus) %>% 
  summarise(Paff = mean(a_p),
            Iaff = mean(a_i),
            mu = mean(mu_p),
            Code = ifelse(length(Code)>1, Code[[1]], Code)) %>% 
  select(-genus)

# write.csv(sp_trt, "./data/species_traits.csv",row.names = F, quote = F)







  
