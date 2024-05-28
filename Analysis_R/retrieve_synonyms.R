# Import libraries
library(taxize)
library(tidyverse)

# Load data
df <- read.csv("../List_all_species_phenology_proj_UDLA2022_filtered_def.csv")

# Create a column for the number of synonyms 
df[, "Number_of_synonyms"] <- NA

# Create a dataframe for original names and synonyms
df_synonyms <- df["Scientific_names"]

# Create an index for the original names
index = 1

for(i in df["Scientific_names"][1:444,]) {
   # Obtain the synonyms with taxize
   temp <- synonyms(i, db = "tropicos", 1)
   
   # Convert list into a dataframe
   synonyms <- synonyms_df(temp)
   
   # Condition to validate that exists synonyms or not
   if(is_empty(synonyms) == TRUE){
     
     df[index, "Number_of_synonyms"] <- 0
     
     index <- index+1
     
     next
     
   }else{
     # Extract only the scientific name of synonyms
     names_synonyms <- synonyms["scientificname"]
     
     colnames(names_synonyms) <- c("Scientific_names")
     
     # Add number of synonyms 
     df[index, "Number_of_synonyms"] <- nrow(names_synonyms)
     
     # Append synonyms to the complete list
     df_synonyms <- rbind(df_synonyms, names_synonyms)
     
     index <- index+1
   }
  
}

# Export data
write_csv(df, "List_all_species_phenology_proj_UDLA2022_n_synonyms.csv")

write_csv(df_synonyms, "List_all_species_phenology_proj_UDLA2022_with_synonyms.csv")





