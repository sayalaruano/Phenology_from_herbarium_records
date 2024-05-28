# Import libraries
library(taxize)
library(tidyverse)

# Load data
df <- read.csv("List_all_species_phenology_proj_UDLA2022_with_synonyms.csv")

# Create an empty df
df_authnames <- data.frame(Scientific_name_auth=character())

# Create an index for the original names
index = 0

# Retrieve scientific names with authors

for(i in df["Scientific_names"][1:2908,]) {
  # Obtain the record inf from tropicos
  temp <- tp_search(i)
  
  # Extract scientific name with authors
  name_with_auth <- temp$scientificnamewithauthors[1]
  
  # Condition to validate that exists a scientific name with authors
  if(is_empty(name_with_auth) == TRUE){
    
    index <- index+1
    print(i)
    print(index)
    next
    
  }else{
    
    # Append scientific name with authors to the complete list
    df_authnames <- rbind(df_authnames, name_with_auth)
    index <- index+1
  }
  
}

colnames(df_authnames) <- c("Scientific_name_auth")

# Export data
write_csv(df_authnames, "List_50_species_withauth_and_synonyms_phenology_proj_UDLA2022_def.csv")

sum(duplicated(df_authnames))
