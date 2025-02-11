library(tools)
library(tidyverse)

file.names  <- list.files(path ="C:/Users/symiakaki/Desktop/fc/corrected_training_erken", 
                          pattern = "clean", 
                          full.names = T, recursive=T) 

# start file loop
result.table  <- NULL

# file = file.names[3]

for(i in 1:length(file.names)) {
  
  # read in .csv file
  file = file.names[i]
  print(paste0("working on ", basename(file)))
  

  # subset data
  data <- read.csv(file, head=T, sep=",")  %>% 
    filter(PreprocessingTrue %in% c("object", ""),
           !LabelTrue %in% c(""
                             # , "detritus", 
                             # "CyaGlo000_filaments", 
                             # "CyaGloech_0000000"
                             )
           ) %>% 
    select(Name, Date, AbdDiameter, LabelTrue) 

  result.table <- rbind(result.table, data)
}

sort(unique(result.table$LabelTrue))

saveRDS(result.table, file = "./data/play_data.rds")






