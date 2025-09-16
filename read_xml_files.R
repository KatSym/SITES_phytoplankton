library(xml2)
library(tidyverse)

read_spreadsheetml <- function(file_path) {
  xml <- read_xml(file_path)
  
  sheets <- xml_find_all(xml, ".//ss:Worksheet", 
                         ns = c(ss = "urn:schemas-microsoft-com:office:spreadsheet"))
  
  if (length(sheets) < 2) {
    warning(paste("File", file_path, "has less than 2 sheets. Skipping."))
    return(NULL)
  }
  
  extract_sheet <- function(sheet) {
    rows <- xml_find_all(sheet, ".//ss:Row", 
                         ns = c(ss = "urn:schemas-microsoft-com:office:spreadsheet"))
    data <- lapply(rows, function(row) {
      cells <- xml_find_all(row, ".//ss:Cell/ss:Data", 
                            ns = c(ss = "urn:schemas-microsoft-com:office:spreadsheet"))
      sapply(cells, xml_text)
    })
    df <- do.call(rbind, lapply(data, function(x) {
      length(x) <- 2
      return(x)
    }))
    as.data.frame(df, stringsAsFactors = FALSE)
  }
  
  # Read both sheets
  df1_raw <- extract_sheet(sheets[[1]])
  df2_raw <- extract_sheet(sheets[[2]])
  
  label1 <- trimws(df1_raw[1, 2])
  label2 <- trimws(df2_raw[1, 2])
  
  # Normalize labels
  label1 <- if (grepl("Fragilaria W", label1, ignore.case = TRUE)) "Fragilaria W" else
    if (grepl("Frag(ilaria)? L", label1, ignore.case = TRUE)) "Fragilaria L" else label1
  
  label2 <- if (grepl("Fragilaria W", label2, ignore.case = TRUE)) "Fragilaria W" else
    if (grepl("Frag(ilaria)? L", label2, ignore.case = TRUE)) "Fragilaria L" else label2
  
  # Clean data: drop first 16 rows, rename and convert to numeric
  df1 <- df1_raw[-c(1:14), ]
  df2 <- df2_raw[-c(1:14), ]
  
  names(df1) <- c("Key", "Value")
  names(df2) <- c("Key", "Value")
  
  df1$Value <- as.numeric(df1$Value)
  df2$Value <- as.numeric(df2$Value)
  
  # Prepare a full data frame with consistent column names
  merged <- full_join(df1, df2, by = "Key", suffix = c("_1", "_2"))
  
  # Assign values based on header labels
  result <- merged %>%
    transmute(
      Key,
      `Fragilaria_W` = case_when(
        label1 == "Fragilaria W" ~ Value_1,
        label2 == "Fragilaria W" ~ Value_2,
        TRUE ~ NA_real_
      ),
      `Frag_L` = case_when(
        label1 == "Fragilaria L" ~ Value_1,
        label2 == "Fragilaria L" ~ Value_2,
        TRUE ~ NA_real_
      ),
      File = tools::file_path_sans_ext(basename(file_path))
    )
  
  return(result)
}


a <- read_spreadsheetml("./data/diatoms_lw/D20/2024_1_30__11_48_39_D20_E01_Fragilaria.xml")

xml_files <- list.files(path = "data/diatoms_lw/", pattern = "\\.xml$", 
                        full.names = TRUE, recursive = T)

all_data <- lapply(xml_files, read_spreadsheetml)

combined_data <- do.call(rbind, all_data) |> 
  separate(File, into = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"), sep = "_", remove = F) |> 
  select(-c("s1", "s2", "s3", "s4", "s5", "s6", "s7")) |> 
  separate(s8, into = c("day", "mes", "st1", "st2"), sep = "-") |> 
  select(-c("st1", "st2")) |> 
  # create more mess to fix the mess i.e.
  # have consistent mesocosm and day naming
  mutate(day = ifelse(day =="", s9, day),
         day = ifelse(File == "2024_10_22__8_40_18___32_LE1_Fragilaria_Fragilaria W", "32", day),
         mes = ifelse(is.na(mes), s10, mes),
         mes = case_when(mes == "Frag L" ~ s9,
                              mes == "Fragilaria W" ~ s9, 
                              mes == "Fragilaria" ~ s9,
                              File == "2024_10_22__8_40_18___32_LE1_Fragilaria_Fragilaria W" ~ "LE",
                              .default = mes)) |> 
  select(-s9, -s10) |> 
  mutate(mesocosm = case_when(mes == "D2" ~ "8",
                              mes == "D3" ~ "11",
                              mes == "E3" ~ "9",
                              mes == "I4" ~ "14",
                              mes == "LE" ~ "ERK",
                              mes == "LE1" ~ "ERK",
                              .default = mes),
         mesocosm = parse_number(mesocosm),
         day = parse_number(day),
         Treatment = case_when(mesocosm == 1 | mesocosm == 7 | mesocosm == 10 | mesocosm == 16 ~ "C",
                               mesocosm == 2 | mesocosm == 8 | mesocosm == 11 | mesocosm == 13 ~ "D",
                               mesocosm == 3 | mesocosm == 5 | mesocosm == 12 | mesocosm == 14 ~ "I",
                               mesocosm == 4 | mesocosm == 6 | mesocosm == 9 | mesocosm == 15 ~ "E",
                               .default = "ERK"),
         .keep = "unused") |> 
  # fix data issues from the notes file
  mutate(Frag_L = case_when(File == "2024_1_30__10_27_33_D00-E01_Fragilaria" ~ 357.1,
                            .default = Frag_L)) |> 
  rename(width = Fragilaria_W,
         length = Frag_L)

write.csv(combined_data, "data/fragilaria_measurements.csv", row.names = F, quote = F)
