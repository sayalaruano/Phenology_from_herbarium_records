README

INITIAL LIST

The initial list comprised 636 unique species, stored in the "List_all_species_phenology_proj_UDLA2022.csv" file. 

FIRST FILTERS

I filtered this initial list to delete 17 records without species information, 1 null record, and 37 records with the genus "Indeterminada"; 
obtaining a list of 581 species, available in the "List_all_species_phenology_proj_UDLA2022_filtered.csv" file.

ADDITIONAL FILTERS

First, I removed the aff and cf words from the names.

affinis (aff) - 9
conferatur (cf) - 73

Then, I deleted the duplicates after rmeoving the mentioned words

duplicated records: 35

After removing duplicates: 546

I created a list with the 47 records that have cf and aff in their names, which is available in the "List_47records_cff_aff.csv" file.

Also, I added to the list records of species of plants with information about subspecies and variants.  

subspecies (subsp) - 4
variant (var) - 5

After adding species with subsp and var: 554

This list is available in the "List_all_species_phenology_proj_UDLA2022_filtered_aff_cf_subsp_var.csv" file.

FINAL FILTER

I deleted all the records identified at genus level, which looks like "Banara 1inECPI_VERD_02", and 4 extra duplicates. 

Final list size: 444

The final list is available in the "List_all_species_phenology_proj_UDLA2022_filtered_def.csv" file.

SCRIPTS FOR DATA PROCESSING

All the code used to apply the filters of this part is available in the "EDA_UDLA2022_Andeanforestproj.ipynb" file.
