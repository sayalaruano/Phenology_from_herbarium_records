# Import libraries
library(taxize)
library(tidyverse)

# Load data
df <- read.csv("../List_species/List_all_species_phenology_proj_UDLA2022_filtered_def.csv")

# Create a column for the number of synonyms 
df[, "Number_of_synonyms"] <- NA

# Create a dataframe for original names and synonyms
df_synonyms <- df["Scientific_names"]

df_synonyms[, "Original_names"] <- df["Scientific_names"]

# Create an index for the original names
index = 1

for(i in df["Scientific_names"][1:438,]) {
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
     
     # Create a df of the size of syn with the orig name
     names_synonyms["Original_names"] <- rep(i, each=nrow(names_synonyms))
     
     # Change column names 
     colnames(names_synonyms) <- c("Scientific_names", "Original_names")
     
     # Add number of synonyms 
     df[index, "Number_of_synonyms"] <- nrow(names_synonyms)
     
     # Append synonyms to the complete list
     df_synonyms <- rbind(df_synonyms, names_synonyms)
     
     index <- index+1
   }
  
}

# Export data
write_csv(df, "List_all_species_phenology_proj_UDLA2022_number_synonyms.csv")

write_csv(df_synonyms, "List_all_species_phenology_proj_UDLA2022_with_synonyms.csv")





