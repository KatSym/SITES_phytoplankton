library(xml2)
library(dplyr)
library(tools)

read_spreadsheetml <- function(file_path) {
  xml <- read_xml(file_path)
  
  # Get all worksheets
  sheets <- xml_find_all(xml, ".//ss:Worksheet", 
                         ns = c(ss = "urn:schemas-microsoft-com:office:spreadsheet"))
  
  if (length(sheets) < 2) {
    warning(paste("File", file_path, "has less than 2 sheets. Skipping."))
    return(NULL)
  }
  
  # Function to extract data from a worksheet
  extract_sheet <- function(sheet) {
    rows <- xml_find_all(sheet, ".//ss:Row", 
                         ns = c(ss = "urn:schemas-microsoft-com:office:spreadsheet"))
    
    data <- lapply(rows, function(row) {
      cells <- xml_find_all(row, ".//ss:Cell/ss:Data", 
                            ns = c(ss = "urn:schemas-microsoft-com:office:spreadsheet"))
      sapply(cells, xml_text)
    })
    
    df <- do.call(rbind, lapply(data, function(x) {
      length(x) <- 2  # Ensure two columns (pad with NA if needed)
      return(x)
    }))
    
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    names(df) <- c("Key", "Value")
    return(df)
  }
  
  df1 <- extract_sheet(sheets[[1]])[-c(1:16),]
  df2 <- extract_sheet(sheets[[2]])[-c(1:16),]
  
  # Merge the sheets by Key
  df_merged <- df1 %>%
    rename(Value1 = Value) %>%
    inner_join(df2 %>% rename(Value2 = Value), by = "Key") %>%
    mutate(Value1 = as.numeric(Value1),
           Value2 = as.numeric(Value2),
           File = file_path_sans_ext(basename(file_path)))  # Add filename column
  
  return(df_merged)
}


xml_files <- list.files(path = "data/diatoms_lw/", pattern = "\\.xml$", 
                        full.names = TRUE, recursive = T)

all_data <- lapply(xml_files, read_spreadsheetml)

combined_data <- do.call(rbind, all_data) |> 
  mutate(length = pmax(Value1, Value2),
         width = pmin(Value1, Value2),
         .keep = "unused"
         )



a <- read_spreadsheetml("./data/diatoms_lw/D04/2024_1_31__11_46_24_D04-E01-big-chamber_Frag L.xml")
