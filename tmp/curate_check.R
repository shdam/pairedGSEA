library(dplyr)
library(readxl)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)

# Load list of files into dataframe
files2df <- function(filename, name = "filename"){
  df <- filename %>% 
    list.files(full.names = TRUE) %>% 
    as.tibble() %>% 
    filter(!str_detect(value, "_delete")) %>% 
    mutate(geo = value %>% str_extract("\\d+\\_GSE\\d+")) %>% 
    separate(geo, into = c("id", "geo"), sep = "_", convert = TRUE)
  colnames(df) <- c(name, "id", "geo")
  return(df)
  
}
df1 <- files2df("~/Dropbox/05_0_individual_study_excel_org", name = "shd")
df2 <- files2df("~/Dropbox/05_1_individual_study_excel_KVS", name = "kvs")

# Join the two curated lists
df <- df1 %>% 
  full_join(df2, by = c("id", "geo")) %>% 
  filter(!(id > 125))
df$compare <- NA

# Comparison of each file
items <- df %>% filter(!is.na(shd) & !is.na(kvs)) %>% pull(id)
for(item in items){
  print(item)
  shd <- df %>% filter(id == item) %>% 
    pull(shd) %>% 
    read_xlsx() %>% 
    filter(!is.na(`comparison (baseline_v_condition)`)) %>% 
    pull(`comparison (baseline_v_condition)`) %>% 
    sort()
  kvs <- df %>% filter(id == item) %>% 
    pull(kvs) %>% 
    read_xlsx() %>% 
    filter(!is.na(`comparison (baseline_v_condition)`)) %>% 
    pull(`comparison (baseline_v_condition)`) %>% 
    sort()
  
  if(length(shd) == length(kvs)){
    if(all(shd == kvs)){
      df[df$id == item,]$compare <- TRUE
    } else{
      df[df$id == item,]$compare <- FALSE
    }
  } else{
    df[df$id == item,]$compare <- FALSE
  }
    
}

dfFalse <- filter(df, !compare)

# Function that extracts differences
getComparison <- function(df, com_id){
  print(paste("File id:", com_id))
  shd <- df %>% filter(id == com_id) %>% 
    pull(shd) %>% 
    read_xlsx() %>% 
    filter(!is.na(`comparison (baseline_v_condition)`)) %>% 
    select(`comparison (baseline_v_condition)`, `comparison_title (empty_if_not_okay)`) %>% 
    arrange(`comparison (baseline_v_condition)`)
  kvs <- df %>% filter(id == com_id) %>% 
    pull(kvs) %>% 
    read_xlsx() %>% 
    filter(!is.na(`comparison (baseline_v_condition)`)) %>% 
    select(`comparison (baseline_v_condition)`, `comparison_title (empty_if_not_okay)`) %>% 
    arrange(`comparison (baseline_v_condition)`)
  return(list("shd"=shd, "kvs"=kvs))
}

# Manual examination of differences
getComparison(df, dfFalse$id[1])

# Deleted files
dfNA <- filter(df, is.na(compare))
dfNA$id
