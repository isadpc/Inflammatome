# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(biomaRt)
library(httr)
library(jsonlite)

# Inflammation marker (gold set) -----------------------------------------------

#PMIDs extracted on 28/04/2023 using EDirect in UNIX console
#Text-mining API (Jensen lab) accessed on 01/05/2023

# Load PMIDs from 28/04/2023
PMIDs <- read_tsv(paste0("data/", "PMIDs_markers_inflammation.txt"), col_names = c("PMIDs"))
PMIDs_string <- paste(as.character(PMIDs$PMIDs), collapse = " ")

call <- "https://api.jensenlab.org/Textmining" #previously "https://api11.jensenlab.org/Textmining"
body <- list(documents = PMIDs_string, format = "json", limit = "100", type2 = "9606")

# Formulate request

API_request <- POST(call, body = body, verbose())
http_status(API_request)
http_type(API_request)

# Extract results

API_output <- content(API_request, "text", encoding = "ISO-8859-1")
API_output_json <- jsonlite::fromJSON(API_output)

API_output_DF <- base::as.data.frame(API_output_json) %>% 
  dplyr::select(-TRUE.)

API_output_DF_score <- API_output_DF %>% 
  dplyr::select(ends_with("score")) %>% 
  pivot_longer(everything(), names_to = "ENSP.ID", values_to = "Textmining.score") %>% 
  separate_wider_delim("ENSP.ID", delim = ".", names = c("ENSP.ID", "trash")) %>% 
  dplyr::select(-"trash")

API_output_DF_name <- API_output_DF %>% 
  dplyr::select(ends_with("name")) %>% 
  pivot_longer(everything(), names_to = "ENSP.ID", values_to = "Protein.name") %>% 
  separate_wider_delim("ENSP.ID", delim = ".", names = c("ENSP.ID", "trash")) %>% 
  dplyr::select(-"trash")

API_output_DF_final <- merge(API_output_DF_name, API_output_DF_score, by = "ENSP.ID") %>% 
  arrange(desc(Textmining.score)) %>% 
  remove_rownames()

# Load API output from 01/05/2023
API_output_DF_final <- read_tsv(paste("data/", "Textmining_Top100.tsv", sep = ""))

# Load identifier list from STRING (STRING V11.5, from 03/05/2023)
STRING_ENSP_to_ENSG <- read_tsv(paste0("data/", "human.aliases.filtered.txt"), col_names = c("SPLIT", "ENSG.ID", "data.type")) %>% 
  separate("SPLIT", c("tax", "ENSP.ID")) %>% 
  dplyr::select(-c(tax, data.type)) 

# Create annotated inflammation marker list
Inflamm_Top100_marker <- API_output_DF_final %>% 
  left_join(STRING_ENSP_to_ENSG, by = c("ENSP.ID")) %>% 
  dplyr::rename("Gene.name" = "Protein.name")

# Load annotation list from 03/05/2023
Inflamm_Top100_marker <- read_tsv(paste0("data/", "Inflammation_markers_Top100.tsv"))
