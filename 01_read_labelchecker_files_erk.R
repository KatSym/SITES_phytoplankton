library(tools)
library(tidyverse)

file.names.training  <- list.files(path ="C:/Users/symiakaki/Desktop/fc/corrected_training_erken", 
                          pattern = "clean", 
                          full.names = T, recursive=T) 
# file.names.all  <- list.files(path ="data/LC", 
#                           pattern = "LabelChecker", 
#                           full.names = T, recursive=T) 

file.names.all  <- list.files(path ="C:/Users/symiakaki/Desktop/fc/erken_data_with_biovol", 
                              pattern = "LabelChecker", 
                              full.names = T, recursive=T) 

biovol.names <-  list.files(path ="C:/Users/symiakaki/Desktop/fc/biovol_LC_Erken/", 
                                               pattern = "LabelChecker", 
                                               full.names = T, recursive=T)

# start file loop
result.table.all  <- NULL
for(i in 1:length(file.names.all)) {
  
  # read in .csv file
  file = file.names.all[i]
  print(paste0("working on ", basename(file)))
  
  
  # subset data
  data <- read_csv(file, col_types = cols(.default = "c"),
                   guess_max = Inf#head=T, sep="," 
                   #colClasses = c(PreprocessingTrue = "character",
                   #                 LabelTrue = "character",
                   #                 LabelPredicted = "character",
                   #                 ProbabilityScore = "numeric"
                   #               ),
                   # stringsAsFactors = FALSE, check.names = FALSE
  ) %>% 
    filter(is.na(PreprocessingTrue) | PreprocessingTrue == "object") %>% 
    select(Name, Date, Uuid, AbdDiameter, Length, BiovolumeHSosik, SurfaceAreaHSosik,
           LabelPredicted, ProbabilityScore, LabelTrue)
  
  
  result.table.all <- rbind(result.table.all, data)
}


# training data
result.table.tr  <- NULL
for(i in 1:length(file.names.training)) {
  
  # read in .csv file
  file = file.names.training[i]
  print(paste0("working on ", basename(file)))
  

  # subset data
  data <- read_csv(file, col_types = cols(.default = "c"),
                   guess_max = Inf#head=T, sep="," 
                   #colClasses = c(PreprocessingTrue = "character",
                   #                 LabelTrue = "character",
                   #                 LabelPredicted = "character",
                   #                 ProbabilityScore = "numeric"
                   #               ),
                   # stringsAsFactors = FALSE, check.names = FALSE
                   ) %>% 
    filter(is.na(PreprocessingTrue) | PreprocessingTrue == "object") %>% 
    select(Name, Date, Uuid, AbdDiameter, Length, LabelPredicted, LabelTrue, ProbabilityScore)
     
    
  result.table.tr <- rbind(result.table.tr, data)
}


#sort(unique(result.table.all$LabelTrue))



all.data <- result.table.all %>% 
  left_join(result.table.tr, by = c("Name", "Uuid", "Date", "AbdDiameter", "Length")) %>% 
  mutate(LabelTrue=coalesce(LabelTrue.x,LabelTrue.y),
         LabelPredicted = LabelPredicted.x,
         .keep = "unused") %>% 
  select(-LabelPredicted.y,  -ProbabilityScore.y) %>% 
  rename(biovolMS = BiovolumeHSosik,
         surfacearea = SurfaceAreaHSosik,
         ProbabilityScore = ProbabilityScore.x)
  
 

all.data$LabelTrue[all.data$LabelTrue == "DET"] <- "detritus"
all.data$LabelTrue[all.data$LabelTrue == "Detritus"] <- "detritus"
all.data$LabelTrue[all.data$LabelTrue == "Thin_stripy_FIL"] <- "FILAMENTS"
all.data$LabelTrue[all.data$LabelTrue == "Thick_FIL"] <- "FILAMENTS"
all.data$LabelTrue[all.data$LabelTrue == "Amoebida_ 761"] <- "Amoebida_761"
all.data$LabelTrue[all.data$LabelTrue == "EurKercoc_1002590"] <- "ROTIFERA"
all.data$LabelTrue[all.data$LabelTrue == "EurKel000_1002568"] <- "ROTIFERA"
all.data$LabelTrue[all.data$LabelTrue == "EurTrisim_1002086"] <- "ROTIFERA"
all.data$LabelTrue[all.data$LabelTrue == "EurPol000_1003101"] <- "ROTIFERA"
all.data$LabelTrue[all.data$LabelTrue == "EVERYTHING"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "duplicate"] <- ""
# all.data$LabelTrue[all.data$LabelTrue == "Multispecies"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "else"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "ROTIFERS"] <- "ROTIFERA"
all.data$LabelTrue[all.data$LabelTrue == "CopNaup_0000000"] <- "Mesozoo_0000000"
all.data$LabelTrue[all.data$LabelTrue == "ZOO"] <- "Mesozoo_0000000"

all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Acanthoceras"] <- "BacAca000_0000000"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Closterium"] <- "ZygCloaci_2646560"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Eudorina"] <- "ChlEud000_0000000"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Closterium-like"] <- "ZygCloaci_2646560"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Pediastrum"] <- "ChlPed000_2641883"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Apiocystis"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Coelastrum"] <- "ChlCoe000_0000000"
all.data$LabelTrue[all.data$LabelTrue == "CYANO_Dolichospermum"] <- "CyaDol000_5428755"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Asterococcus"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Staurastrum"] <- "ConSta000_2647648"
all.data$LabelTrue[all.data$LabelTrue == "CHLORO_Dictyosphaerium"] <- "TreDic000_2641996"

all.data$LabelTrue[all.data$LabelTrue == "CYANO_Microcystis"] <- "CyaMic000_3217749"
all.data$LabelTrue[all.data$LabelTrue == "CYANO_Woronichinia"] <- "CyaWor000_3216993"

all.data$LabelTrue[all.data$LabelTrue == "DIATOM_Fragillaria"] <- "BacFracro_3192403"
all.data$LabelTrue[all.data$LabelTrue == "DIATOM_Stephanodiscus"] <- "BacSte000_3193250"
all.data$LabelTrue[all.data$LabelTrue == "DIATOM_Acanthoceras"] <- "BacAca000_0000000"
all.data$LabelTrue[all.data$LabelTrue == "DIATOM_centric_chain"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "DIATOM_Nitszchia"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "DIATOM_Asterionella_formosa"] <- "BacAstfor_3192686"

all.data$LabelTrue[all.data$LabelTrue == "prob_gloeo_strand"] <- ""
all.data$LabelTrue[all.data$LabelTrue == "DINO_Ceratium_hirundinella"] <- "DinCerhir_7598904"
all.data$LabelTrue[all.data$LabelTrue == "CILIATE"] <- "Ciliate_0000000"




sort(unique(all.data$LabelTrue))
length(unique(all.data$Name))

saveRDS(all.data, file = "./data/Erken_phyto_all_20250626.rds")


