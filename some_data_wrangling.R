library(tidyverse)
library(plotly)
library(gapminder)

data <- readRDS("./data/play_data.rds") 

metadata <- read.csv("./data/Metadata_Erken.csv", 
                     header=T, sep=";") %>% 
  filter(ID %in% data$Name) %>% 
  select(ID, vol_FC..ml., finalfactor) %>% 
  rename(volume_run = vol_FC..ml.)



dat <- data %>% 
  # filter out unwanted classes
  filter(!LabelTrue %in% c(
                            "detritus", 
                            "CyaGlo000_filaments",
                            "CyaGloech_0000000",
                            "Mesozoo_0000000"
                          )) %>% 
  # remove '_updated' from sample name (don't know how it got there)
  mutate(Name = str_remove(Name, "_updated")) %>% 
  group_by(Date, Name, LabelTrue) %>% 
  summarise(count = n(),
            mean_abd = mean(AbdDiameter)) %>% 
  separate(Name, into = c("mesocosm", "day", "magn", "conc", "replicate"), 
           sep = "_", remove = F) %>% 
  select(-Date, -magn, -conc) %>% 
  mutate(Date = as.Date(day, "%Y%m%d")) %>% 
  left_join(., metadata, by =join_by(Name == ID)) %>% 
  mutate(
    
    # calculate abundance
    abund = 1000*count*finalfactor, # in individuals/L

    # calculate biovolume
    # vol.calc = (4/3)*pi*(mean_abd/2)^3,
    # biovol = (sum(vol.calc)/volume_run)* finalfactor,  # in um^3/mL
    Mes_ID = as.numeric(gsub("\\D", "", mesocosm)),
    Treatment = case_when(Mes_ID == 1 | Mes_ID == 7 | Mes_ID == 10 | Mes_ID == 16 ~ "C",
                          Mes_ID == 2 | Mes_ID == 8 | Mes_ID == 11 | Mes_ID == 13 ~ "D",
                          Mes_ID == 3 | Mes_ID == 5 | Mes_ID == 12 | Mes_ID == 14 ~ "I",
                          Mes_ID == 4 | Mes_ID == 6 | Mes_ID == 9 | Mes_ID == 15 ~ "E",
                          .default = "ERK"),
    
    .keep = "unused"
          )
# %>% 
#   
#   group_by(Date, Treatment, Mes_ID, LabelTrue) %>% 
#   summarise(mean(abund))


trt.cols <- c(`C`= "#000000", #black - C
              `D`= "#0a8754", #g - D
              `I`= "#4472ca", #blu - I
              `E`= "#e84855", #r - E
              `ERK` = "#ffb703") 



dat %>% 
  filter(LabelTrue == "ZygCloaci_2646560") %>% 
  ggplot(., aes(x = as.factor(Date), y = abund,  fill = Treatment))+
  geom_jitter(aes(color = Treatment), size=0.4, alpha=0.8) +
  geom_boxplot() +
  scale_fill_manual(values = trt.cols)


