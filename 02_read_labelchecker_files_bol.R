library(tools)
library(tidyverse)

file.names.training  <- list.files(path ="C:/Users/symiakaki/Desktop/fc/bolmen_data_with_biovol", 
                          pattern = "LabelChecker", 
                          full.names = T, recursive=T) 
# file.names.all  <- list.files(path ="data/LC", 
#                           pattern = "LabelChecker", 
#                           full.names = T, recursive=T) 

file.names.all  <- list.files(path ="C:/Users/symiakaki/Desktop/fc/bolmen_LCs_20250625-all/", 
                              pattern = "LabelChecker", 
                              full.names = T, recursive=T) 


# start file loop
result.table.all  <- NULL


# file = file.names.all[13]

for(i in 1:length(file.names.all)) {
  
  # read in .csv file
  file = file.names.all[i]
  print(paste0("working on ", basename(file)))
  

  # subset data
  data <- read_csv(file, 
                   # col_types = cols(.default = "c"),
                    guess_max = Inf, #head=T, sep="," 
                                      col_types = cols(PreprocessingTrue = col_character(),
                                                       LabelPredicted = col_character(),
                                                       ProbabilityScore = col_double(),
                                                       LabelTrue = col_character())
                   ) %>% 
    filter(is.na(PreprocessingTrue) | PreprocessingTrue == "object") %>% 
    select(Name, Date, Uuid, AbdDiameter, Length, any_of(c("BiovolumeHSosik", "BiovolumeMS")), any_of(c("SurfaceAreaHSosik", "SurfaceAreaMS")),
           LabelPredicted, ProbabilityScore, LabelTrue) %>% 
    rename(biovolume = any_of(c("BiovolumeHSosik", "BiovolumeMS")),
           surfacearea = any_of(c("SurfaceAreaHSosik", "SurfaceAreaMS")))
# 
#   data <- read_csv(file, col_names = T, 
#                    col_types = cols(PreprocessingTrue = col_character(),
#                                     LabelTrue = col_character(),
#                                     LabelPredicted = col_character(),
#                                     ProbabilityScore = col_double()),
#                    guess_max = 10000)
#     
    
  result.table.all <- rbind(result.table.all, data)
}
unique(result.table.all$Name)


result.table.tr  <- NULL
for(i in 1:length(file.names.training)) {
  
  # read in .csv file
  file = file.names.training[i]
  print(paste0("working on ", basename(file)))
  
  
  # subset data
  data <- read_csv(file, 
                   # col_types = cols(.default = "c"),
                   # guess_max = Inf, #head=T, sep="," 
                   # col_types = cols(PreprocessingTrue = col_character(),
                   #                  LabelTrue = col_character())
  ) %>% 
    filter(is.na(PreprocessingTrue) | PreprocessingTrue == "object") %>% 
    select(Name, Uuid, AbdDiameter, Length, LabelTrue)
  # 
  #   data <- read_csv(file, col_names = T, 
  #                    col_types = cols(PreprocessingTrue = col_character(),
  #                                     LabelTrue = col_character(),
  #                                     LabelPredicted = col_character(),
  #                                     ProbabilityScore = col_double()),
  #                    guess_max = 10000)
  #     
  
  result.table.tr <- rbind(result.table.tr, data)
}

unique(result.table.tr$Name)

result.table.tr <- result.table.tr %>% 
  mutate(Name = str_remove(Name, "_updated"))



#sort(unique(result.table.all$LabelTrue))

result.table.all$AbdDiameter <- as.numeric(result.table.all$AbdDiameter)
result.table.all$Length <- as.numeric(result.table.all$Length)

all.data <- result.table.all %>% 
  left_join(.,result.table.tr, by = c("Name", "Uuid", "AbdDiameter", "Length")) %>% 
  mutate(LabelTrue=coalesce(LabelTrue.x,LabelTrue.y),
         .keep = "unused") %>% 
  select(-Date) 
  

sort(unique(all.data$LabelTrue))
length(unique(all.data$Name))

saveRDS(all.data, file = "./data/Bolmen_FCphyto_20250630-biovol.rds")




result.table.all %>% filter(LabelTrue == "CHLORO_Coelastrum") %>% 
  pull(LabelPredicted)
  
  
result.table.all$Name[result.table.all$LabelTrue== "prob_gloeo_strand"]
