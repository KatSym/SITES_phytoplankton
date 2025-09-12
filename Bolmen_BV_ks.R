# load in packages
library("tidyverse")
library("readxl")

# read in microscope counts from Bolmen 2022
counts <- read.csv("C:/Users/symiakaki/Downloads/phyto.csv") |> 
  mutate(cell_count = Count*Cells, .keep = "unused") |> 
  group_by(Day, Treatment, Mesocosm, Species) |> 
  summarise(count = sum(cell_count),
            fov = n())

# correcting spellings
counts$Species = str_replace(counts$Species, c('Asterococcus'), 'Asteroccoccus')
counts$Species = str_replace(counts$Species, c('Kirchinella'), 'Kirchneriella')
counts$Species = str_replace(counts$Species, c('Pseudanabaeba'), 'Pseudanabaena')
counts$Species = str_replace(counts$Species, c('Coenochloris sp '), 'Coenochloris')
counts$Species = str_replace(counts$Species, c('Selenastrum '), 'Selenastrum')
counts$Species = str_replace(counts$Species, c('Selenastrumbibraianum'), 'Selenastrum bibraianum')


# read in measurement data from Bolmen 2022
data <- read.csv("C:/Users/symiakaki/Downloads/Measurements_bolmen.csv")

# calculating the biovolume for each phytoplankton. Shape formulas and values for hidden dimensions are detailed in:
# European Committee for Standardization (2015) Water quality- Guidance on the estimation of phytoplankton biovolume, European standard BS EN 16695:2015
# https://standards.iteh.ai/catalog/standards/cen/bcc87031-164e-45b9-933a-7db83d4658f4/en-16695-2015?srsltid=AfmBOopH8HXHtfXrao8lTR37Kd-artS7z1FsuUXaTZRmM7aTf1omdF6x

# Biovolume of Aphanocapsa shape= Sphere
Aphanocapsa <- data[ which(data$Species=='Aphanocapsa'), ]
Aphanocapsa$diameter <- (Aphanocapsa$x + Aphanocapsa$y) /2 # 2 measurements recorded for diameter so took the mean
Aphanocapsa$BV <- (1/6)*(pi)*(Aphanocapsa$diameter^3)
Aphanocapsa$diameter <- NULL

# Biovolume of Aphanothece shape= Spheroid
Aphanothece <- data[ which(data$Species=='Aphanothece'), ]
Aphanothece$BV <- (1/6)*(pi)*(Aphanothece$x^2)*(Aphanothece$y)

# Biovolume of Asteroccocus shape= Spheroid
Asteroccoccus <- data[ which(data$Species=='Asteroccoccus'), ]
Asteroccoccus$BV <- (1/6)*(pi)*(Asteroccoccus$x^2)*(Asteroccoccus$y)

# Biovolume of Chlamydomonas shape= Spheroid
Chlamydomonas <- data[ which(data$Species=='Chlamydomonas'), ]
Chlamydomonas$BV <- (1/6)*(pi)*(Chlamydomonas$x^2)*(Chlamydomonas$y)

# Biovolume of Chroomonas shape= ellipsoid v= 1/6*pi*d*d2*h [d2= 0.8*d1]
Chroomonas <- data[ which(data$Species=='Chroomonas'), ]
Chroomonas$BV <- (1/6)*(pi)*(Chroomonas$x)*(Chroomonas$x*0.8)*(Chroomonas$y)

# Biovolume of Chroomonas acuta shape= ellipsoid v= 1/6*pi*d*d2*h [d2= 0.8*d1]
Chroomonas_ac <- data[ which(data$Species=='Chroomonas acuta'), ]
Chroomonas_ac$BV <- (1/6)*(pi)*(Chroomonas_ac$x)*(Chroomonas_ac$x*0.8)*(Chroomonas_ac$y)

# Biovolume of Chrysochromulina parva shape= ellipsoid v= 1/6*pi*d*d2*h [d2= 0.6*d1]
Chrysochromulina_parva <- data[ which(data$Species=='Chrysochromulina parva'), ]
Chrysochromulina_parva$BV <- (1/6)*(pi)*(Chrysochromulina_parva$x)*(Chrysochromulina_parva$x*0.6)*(Chrysochromulina_parva$y)

# Biovolume of Closterium shape= spindle
Closterium <- data[ which(data$Species=='Closterium'), ]
Closterium$BV <- (2/15)*(pi)*(Closterium$x^2)*(Closterium$y)

# Biovolume of Cosmarium bioculatum shape= ellipsoid v= 1/6*pi*d*d2*h [d2= 0.6*d1]
Cosmarium_bio <- data[ which(data$Species=='Cosmarium bioculatum'), ]
Cosmarium_bio$BV <- (1/6)*(pi)*(Cosmarium_bio$x)*(Cosmarium_bio$x*0.6)*(Cosmarium_bio$y)

# Biovolume of Cosmarium regnelli shape= ellipsoid v= 1/6*pi*d*d2*h [d2= 0.6*d1]
Cosmarium_reg <- data[ which(data$Species=='Cosmarium regnelli'), ]
Cosmarium_reg$BV <- (1/6)*(pi)*(Cosmarium_reg$x)*(Cosmarium_reg$x*0.6)*(Cosmarium_reg$y)

# Biovolume of Crucigenia tetrapedia shape= tetrahedron 
Crucigenia_tet <- data[ which(data$Species=='Crucigenia tetrapedia'), ]
Crucigenia_tet$BV <- (1/12)*sqrt(3)*(Crucigenia_tet$x^2)*(Crucigenia_tet$y)

# Biovolume of Cryptomonas shape= Ellipsoid v= 1/6*pi*d*d2*h [d2= 0.8*d1]
Cryptomonas <- data[ which(data$Species=='Cryptomonas'), ]
Cryptomonas$BV <- (1/6)*(pi)*(Cryptomonas$x)*(Cryptomonas$x*0.8)*(Cryptomonas$y)

# Biovolume of 	Cyanodictyon reticulatum shape= Spheroid 
Cyanodictyon_reticulatum <- data[ which(data$Species=='Cyanodictyon reticulatum'), ]
Cyanodictyon_reticulatum$BV <- (1/6)*(pi)*(Cyanodictyon_reticulatum$x^2)*(Cyanodictyon_reticulatum$y)

# Biovolume of 	Cyclotella shape= Cylinder h=0.61*d
Cyclotella <- data[ which(data$Species=='Cyclotella'), ]
Cyclotella$BV <- (1/4)*(pi)*(Cyclotella$x^2)*(Cyclotella$x*0.61)

# Biovolume of 	desmodesmus shape= Spheroid
Desmodesmus <- data[ which(data$Species=='Desmodesmus'), ]
Desmodesmus$BV <- (1/6)*(pi)*(Desmodesmus$x^2)*(Desmodesmus$y)

# Biovolume of 	Dictyospharium shape= Spheroid
Dictyospharium <- data[ which(data$Species=='Dictyosphaerium pulchellum'), ]
Dictyospharium$BV <- (1/6)*(pi)*(Dictyospharium$x^2)*(Dictyospharium$y)

# Biovolume of 	Elakatothrix shape= Spindle
Elakatothrix <- data[ which(data$Species=='Elakatothrix'), ]
Elakatothrix$BV <- (2/15)*(pi)*(Elakatothrix$x^2)*(Elakatothrix$y)

# Biovolume of 	Eunotia shape= (Elliptic cylinder * 1.03) r2=0.5*h h=2*r2
Eunotia <- data[ which(data$Species=='Eunotia'), ]
Eunotia$BV <- ((pi)*(Eunotia$x/2)*(Eunotia$y)*(Eunotia$y*0.5))*1.03

# Biovolume of Golenkenia shape= Sphere
Golenkinia <- data[ which(data$Species=='Golenkinia'), ]
Golenkinia$diameter <- (Golenkinia$x + Golenkinia$y) /2# 2 measurements recorded for diameter so took the mean
Golenkinia$BV <- (1/6)*(pi)*(Golenkinia$diameter^3)
Golenkinia$diameter <- NULL
# Biovolume of Gonium shape= Spheroid
Gonium <- data[ which(data$Species=='Gonium pectorale'), ]
Gonium$BV <- (1/6)*(pi)*(Gonium$x^2)*(Gonium$y)

# Biovolume of 	Kirchneriella shape= Spindle
Kirchneriella <- data[ which(data$Species=='Kirchneriella'), ]
Kirchneriella$BV <- (2/15)*(pi)*(Kirchneriella$x^2)*(Kirchneriella$y)

# Biovolume of Mallomonas shape= ellipsoid d2=0.8*d1
Mallomonas <- data[ which(data$Species=='Mallomonas'), ]
Mallomonas$BV <- (1/6)*(pi)*(Mallomonas$x)*(Mallomonas$x*0.8)*(Mallomonas$y)

# Biovolume of Merismopedia shape= sphere
Merismopedia <- data[ which(data$Species=='Merismopedia'), ]
Merismopedia$diameter <- (Merismopedia$x + Merismopedia$y) /2# 2 measurements recorded for diameter so took the mean
Merismopedia$BV <- (1/6)*(pi)*(Merismopedia$diameter^3)
Merismopedia$diameter <- NULL
# Biovolume of Monoraphidium shape= spindle
Monoraphidium <- data[ which(data$Species=='Monoraphidium contortum'), ]
Monoraphidium$BV <- (2/15)*(pi)*(Monoraphidium$x^2)*(Monoraphidium$y)

# Biovolume of 	Navicula shape= (Elliptic cylinder) h=0.85*d2
Navicula <- data[ which(data$Species=='Navicula'), ]
Navicula$BV <- ((pi)*(Navicula$x/2)*(Navicula$y/2)*(Navicula$y*0.85))

# Biovolume of Oocystis shape= Spheroid
Oocystis <- data[ which(data$Species=='Oocystis'), ]
Oocystis$BV <- (1/6)*(pi)*(Oocystis$x^2)*(Oocystis$y)

# Biovolume of pediastrum shape= cuboid*0.67 c=1*b
Pediastrum <- data[ which(data$Species=='Pediastrum tetras'), ]
Pediastrum$BV <- ((Pediastrum$x)*(Pediastrum$y)*(Pediastrum$x))*0.67

# Biovolume of 	Pseudanabaena shape= Cylinder 
Pseudanabaena <- data[ which(data$Species=='Pseudanabaena'), ]
Pseudanabaena$BV <- (1/4)*(pi)*(Pseudanabaena$x^2)*(Pseudanabaena$y)

# Biovolume of Rhodomonas shape= ellipsoid d2=0.9*d1
Rhodomonas <- data[ which(data$Species=='Rhodomonas'), ]
Rhodomonas$BV <- (1/6)*(pi)*(Rhodomonas$x)*(Rhodomonas$x*0.9)*(Rhodomonas$y)

# Biovolume of Selenastrum bibraianum shape= spindle
Selenastrum_bib <- data[ which(data$Species=='Selenastrum bibraianum'), ]
Selenastrum_bib$BV <- (2/15)*(pi)*(Selenastrum_bib$x^2)*(Selenastrum_bib$y)

# Biovolume of Small dinoflag shape=  cone with a half sphere h=1.2*d
Dinoflag_sml <- data[ which(data$Species=='Small dinoflagellate'), ]
Dinoflag_sml$BV <- (1/12)*(pi)*(Dinoflag_sml$x^2)*((Dinoflag_sml$y)+1/2*Dinoflag_sml$x)

# Biovolume of 	Synedra shape= (Elliptic cylinder) d2= 1*h h=1*d2
Synedra <- data[ which(data$Species=='Synedra'), ]
Synedra$BV <- (1/4*(pi)*(Synedra$x)*(Synedra$y)*(Synedra$y))

# Biovolume of Small dinoflag shape=  cone with a half sphere h=1.2*d
Synura <- data[ which(data$Species=='Synura'), ]
Synura$BV <- (1/12)*(pi)*(Synura$x^2)*((Synura$y)+1/2*Synura$x)

# Biovolume of Trachelomonas shape= Sphere
Trachelomonas <- data[ which(data$Species=='Trachelomonas'), ]
Trachelomonas$diameter <- (Trachelomonas$x + Trachelomonas$y) /2# 2 measurements recorded for diameter so took the mean
Trachelomonas$BV <- (1/6)*(pi)*(Trachelomonas$diameter^3)
Trachelomonas$diameter <- NULL


## Biovolumes for taxa without measurements - from databases/literature

# Biovolume of Coelastrum shape= Sphere
Coelastrum <- data.frame(Species = "Coelastrum")
Coelastrum$diameter <- mean(2:30) # AgleaBase
Coelastrum$BV <- (1/6)*(pi)*(Coelastrum$diameter^3)

# Biovolume of Coenochloris shape= Sphere
Coenochloris<- data.frame(Species = "Coenochloris")
Coenochloris$diameter <- mean(5:10) # https://oak.go.kr/central/journallist/journaldetail.do?article_seq=16531
Coenochloris$BV <- (1/6)*(pi)*(Coenochloris$diameter^3)

# Biovolume of Cosmarium shape= ellipsoid d2 = ?
Cosmarium<- data.frame(Species = "Cosmarium")
Cosmarium$len <- mean(17:80) # https://doaj.org/article/af196a8c39284085851a087a5efc9aa2
Cosmarium$wid <- mean(14:57)
Cosmarium$BV <- (1/6)*(pi)*Cosmarium$len*Cosmarium$wid*(Cosmarium$len/2) # d2 arbitrary

# Biovolume of Cosmarium crenatum shape= ellipsoid 
# no info, biovolume same as Cosmarium
Cosmarium_cre<- data.frame(Species = "Cosmarium crenatum")
Cosmarium_cre$BV <- Cosmarium$BV


# Biovolume of Cosmarium meneghinii shape= ellipsoid d2 = ?
Cosmarium_men<- data.frame(Species = "Cosmarium meneghinii")
Cosmarium_men$len <- mean(14:21) # https://www.outerhebridesalgae.uk/desmids/desmid-species.php?id=728
Cosmarium_men$wid <- mean(11:18)
Cosmarium_men$d2 <- mean(7:10)
Cosmarium_men$BV <- (1/6)*(pi)*Cosmarium_men$len*Cosmarium_men$wid*Cosmarium_men$d2 

# Biovolume of Euastrum binale shape= ellipsoid d2 = ?
Euastrum_bin <- data.frame(Species = "Euastrum binale")
Euastrum_bin$len <- mean(15:30) # http://protist.i.hosei.ac.jp/pdb/images/Chlorophyta/Euastrum/binale/index.html
Euastrum_bin$wid <- mean(12.5:21)
Euastrum_bin$BV <- (1/6)*(pi)*Euastrum_bin$len*Euastrum_bin$wid*(Euastrum_bin$len/2) # d2 arbitrary

# Biovolume of Eudorina shape= Sphere
Eudorina<- data.frame(Species = "Eudorina")
Eudorina$diameter <- mean(8:10) # http://protist.i.hosei.ac.jp/pdb/images/chlorophyta/eudorina/elegans/sp_5.html#:~:text=A%20dividing%20cell.%20elegans%20Ehrenberg:%20A%20colony,(Illustrations%20of%20The%20Japanese%20Fresh%2Dwater%20Algae%2C%201977).
Eudorina$BV <- (1/6)*(pi)*(Eudorina$diameter^3)

# Biovolume of Frustulia shape= Elliptic cylinder
Frustulia<- data.frame(Species = "Frustulia")
Frustulia$len <- mean(26:195) # https://diatoms.org/genera/frustulia
Frustulia$wid <- mean(6:31)
Frustulia$BV <- (pi)*(Frustulia$len/2)*(Frustulia$wid/2)*(Frustulia$wid*0.85) # same as Navicula

# Biovolume of Gomphonema shape= Elliptic cylinder
Gomphonema<- data.frame(Species = "Gomphonema")
Gomphonema$len <- mean(10:118) # https://diatoms.org/genera/gomphonema
Gomphonema$wid <- mean(3:18)
Gomphonema$BV <- (pi)*(Gomphonema$len/2)*(Gomphonema$wid/2)*(Gomphonema$wid*0.85) # same as Navicula

# Biovolume of Micractinium pusillum shape= Sphere
Micractinium_pus<- data.frame(Species = "Micractinium pusillum")
Micractinium_pus$diameter <- mean(3:13) # https://nordicmicroalgae.org/taxon/micractinium-pusillum/
Micractinium_pus$BV <- (1/6)*(pi)*(Micractinium_pus$diameter^3)

# Biovolume of Micractinium shape= Sphere
Micractinium<- data.frame(Species = "Micractinium")
Micractinium$diameter <- mean(3:13) # same as M. pusillum
Micractinium$BV <- (1/6)*(pi)*(Micractinium$diameter^3)

# Biovolume of Paulschulzia  shape= Sphere
Paulschulzia<- data.frame(Species = "Paulschulzia")
Paulschulzia$diameter <- mean(5:13) # https://www.algaebase.org/search/genus/detail/?genus_id=Te751094cce081366
Paulschulzia$BV <- (1/6)*(pi)*(Paulschulzia$diameter^3)

# Biovolume of Scenedesmus shape= Spheroid
Scenedesmus <- data.frame(Species = "Scenedesmus")
Scenedesmus$len <- mean(5:15) # https://www.sciencedirect.com/science/article/abs/pii/B978012741550550012X
Scenedesmus$wid <- mean(2:10) # https://www.sciencedirect.com/science/article/pii/S0960852421003941
Scenedesmus$BV <- (1/6)*(pi)*(Scenedesmus$len^2)*(Scenedesmus$wid)

# Biovolume of Selenastrum shape= Spindle
Selenastrum <- data.frame(Species = "Selenastrum")
Selenastrum$len <- mean(7:42) # https://www.algaebase.org/search/genus/detail/?genus_id=43450
Selenastrum$wid <- mean(1.5:8) 
Selenastrum$BV <- (2/15)*(pi)*(Selenastrum$len^2)*(Selenastrum$wid)

# putting all the data back together

df_list <- list(Aphanocapsa, Aphanothece, Asteroccoccus, Chlamydomonas, Chroomonas, Chroomonas_ac, 
                Chrysochromulina_parva, Closterium, Cosmarium_bio, Cosmarium_reg, Crucigenia_tet,
                Cryptomonas, Cyanodictyon_reticulatum, Cyclotella, Desmodesmus, Dictyospharium,
                Elakatothrix, Eunotia, Golenkinia, Gonium, Kirchneriella, Mallomonas, Merismopedia,
                Monoraphidium, Navicula, Oocystis, Pediastrum, Pseudanabaena, Rhodomonas, Selenastrum_bib,
                Dinoflag_sml, Synedra, Synura, Trachelomonas)


df <- df_list %>% reduce(full_join, by=c('Species', 'Meas_number', 'x', 'y', 'Day', 'Treatment', 'Shape', 'BV'))

exrtabv.list <- list(Coelastrum, Coenochloris, Cosmarium, Cosmarium_cre, Cosmarium_men,
                  Euastrum_bin, Eudorina, Frustulia, Gomphonema, Micractinium, 
                  Micractinium_pus, Paulschulzia, Scenedesmus, Selenastrum)
extra_bv <- exrtabv.list |> reduce(full_join, by = c("Species", "BV")) |> 
  select(Species, BV)


# calculating mean biovolume for each species per day and treatment 
Mean_BV <- df %>%
  group_by(Day, Treatment, Species) %>%
  summarise(Biovolume = mean(BV)) |> 
  ungroup()

# merging counts and biovolumes
Counts_bv <- merge(counts, Mean_BV, by = c("Day", "Treatment", "Species"), all.x = TRUE) |> 
  full_join(extra_bv, by = "Species") |> 
  mutate(Biovolume = coalesce(Biovolume, BV), 
         .keep = "unused")

# checking for species with no biovolumes
unique(Counts_bv$Species[is.na(Counts_bv$Biovolume)])

# for the day/treatment/species combination that we don't have measurements 
# take mean taxon biovolume
Counts_bv <- Counts_bv |> 
  mutate(Biovolume = case_when(Species == "Kirchneriella lunaris" ~ mean(Kirchneriella$BV),
                               Species == "Kirchneriella" ~ mean(Kirchneriella$BV),
                               Species == "Monoraphidium" ~ mean(Monoraphidium$BV),
                               Species == "Aphanothece" ~ mean(Aphanothece$BV),
                               Species == "Chroomonas" ~ mean(Chroomonas$BV),
                               Species == "Crucigenia tetrapedia" ~ mean(Crucigenia_tet$BV),
                               Species == "Desmodesmus" ~ mean(Desmodesmus$BV),
                               Species == "Eunotia" ~ mean(Eunotia$BV),
                               Species == "Mallomonas" ~ mean(Mallomonas$BV),
                               Species == "Merismopedia" ~ mean(Merismopedia$BV),
                               Species == "Navicula" ~ mean(Navicula$BV),
                               Species == "Pseudanabaena" ~ mean(Pseudanabaena$BV),
                               Species == "Synura" ~ mean(Synura$BV),
                               Species == "Asteroccoccus" ~ mean(Asteroccoccus$BV),
                               Species == "Chrysochromulina parva" ~ mean(Chrysochromulina_parva$BV),
                               Species == "Oocystis" ~ mean(Oocystis$BV),
                               Species == "Trachelomonas" ~ mean(Trachelomonas$BV),
                               Species == "Chlamydomonas" ~ mean(Chlamydomonas$BV),
                               Species == "Cryptomonas" ~ mean(Cryptomonas$BV),
                               Species == "Cyanodictyon reticulatum" ~ mean(Cyanodictyon_reticulatum$BV),
                               Species == "Cyclotella" ~ mean(Cyclotella$BV),
                               Species == "Rhodomonas" ~ mean(Rhodomonas$BV),
                               Species == "Synedra" ~ mean(Synedra$BV),
                               Species == "Golenkinia" ~ mean(Golenkinia$BV),
                               Species == "Gonium pectorale" ~ mean(Gonium$BV),
                               Species == "Pediastrum tetras" ~ mean(Pediastrum$BV),
                               Species == "Aphanocapsa" ~ mean(Aphanocapsa$BV),
                               Species == "Elakatothrix" ~ mean(Elakatothrix$BV),
                               Species == "Aphanocapsa" ~ mean(Aphanocapsa$BV),
                               Species == "Dictyosphaerium pulchellum" ~ mean(Dictyospharium$BV),
                               Species == "Closterium" ~ mean(Closterium$BV),
                               Species == "Cosmarium bioculatum" ~ mean(Cosmarium_bio$BV),
                               Species == "Monoraphidium contortum" ~ mean(Monoraphidium$BV),
                               Species == "Selenastrum bibraianum" ~ mean(Selenastrum_bib$BV),
                               Species == "Small dinoflagellate" ~ mean(Dinoflag_sml$BV),
                               .default = Biovolume))


# calculate abundance and biovolume per ml
small_bol <- Counts_bv |> 
  mutate(fov_area = 0.625,
         chamber_area = 452.39,
         sample_vol = 25, # mL
         
         cells_in_chamber = (count*chamber_area)/(fov*fov_area),
         abundance = cells_in_chamber/25,
         biovolume = Biovolume*abundance) |> 
  select(Day, Treatment, Mesocosm, Species, abundance, biovolume)
  
write.csv(small_bol, "data/Bolmen_microscope_counts.csv", row.names = F, quote = F)
