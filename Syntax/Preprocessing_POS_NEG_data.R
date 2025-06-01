chooseCRANmirror(ind=1) #selecting cran mirror option
packages = c("BiocManager", "tidyverse",'stringr','KODAMA', "pheatmap")
for(lib in packages){
if(!lib %in% installed.packages()){
if(lib %in% available.packages()[,1]){
install.packages(lib)} else{ BiocManager::install(lib)}}}

suppressMessages(sapply(packages, require, character = TRUE))	

################################################################################
# Positive ion mode MS data
################################################################################

#---- Data Loading ----

df_pos <- read.csv("Data/POS_Height_1_2024_09_24_05_19_14.csv", stringsAsFactors = FALSE)


#---- Run Quality Check ----
df_pos_run <- df_pos %>%
  filter(Metabolite.name != "Unknown",
         !grepl("no MS2:", Metabolite.name),
         MS.MS.matched != "FALSE")

colnames(df_pos_run)


# Creating a new dataframe with the specified columns

df_pos_run_selected <- df_pos_run %>%
  select(Average.Rt.min., Adduct.type, Metabolite.name, 
         contains('dda'))
		 
# Replace 0 with NA in numeric/intensity columns (excluding metadata)
df_pos_run_selected[ , 4:ncol(df_pos_run_selected)] <- lapply(
  df_pos_run_selected[ , 4:ncol(df_pos_run_selected)], 
  function(x) ifelse(x == 0, NA, x)
)		 
	 

# Verify the selected column names

colnames(df_pos_run_selected)

# Rename samples using lookup table
sample_names <- read.csv("Data/Sample_Names.csv", header = FALSE, stringsAsFactors = FALSE)
names(sample_names) <- c("Original", "Renamed")
rename_lookup <- setNames(sample_names$Renamed, sample_names$Original)
df_pos_run_selected <- df_pos_run_selected %>% rename_with(~ rename_lookup[gsub("_sb.cn.pos_dda1", "", .x)], contains("sample"))


# Explicitly rename blank and pool columns, if needed

df_pos_run_selected <- df_pos_run_selected %>%
  rename(
    Blank_1 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda1,
    Blank_2 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda2,
    Blank_3 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda3,
    Blank_4 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda4,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda2
  )

write.csv(df_pos_run_selected, "df_pos_run_selected.csv" ,row.names=FALSE) 

colnames(df_pos_run_selected)

replicate_cols <- colnames(df_pos_run_selected)[!(colnames(df_pos_run_selected) %in% c("Average.Rt.min.", "Adduct.type", "Metabolite.name"))]


# Reshape to long format
df_long <- df_pos_run_selected %>%
  select(Average.Rt.min., Adduct.type, Metabolite.name, all_of(replicate_cols)) %>%
  pivot_longer(cols = all_of(replicate_cols), names_to = "Sample", values_to = "Area") %>%
  mutate(
    Group = str_remove(Sample, "_\\d+$"),  # Extract base group name
    IsMissing = ifelse(is.na(Area) | Area == 0, 1, 0)  # Mark missing or zero values
  )

# Count missing values per metabolite per group
missing_counts <- df_long %>%
  group_by(Metabolite.name, Average.Rt.min., Adduct.type, Group) %>%
  summarise(MissingReplicates = sum(IsMissing), .groups = "drop")

# Step 1: Identify metabolites where all groups have <2 missing replicates
metabolites_to_keep <- missing_counts %>%
  group_by(Metabolite.name, Average.Rt.min., Adduct.type) %>%
  summarise(MaxMissing = max(MissingReplicates), .groups = "drop") %>%
  filter(MaxMissing < 2)

# Step 2: Filter the original dataset
df_pos_filtered_strict <- df_pos_run_selected %>%
  semi_join(metabolites_to_keep, by = c("Metabolite.name", "Average.Rt.min.", "Adduct.type"))
  

# Step 1: Calculate average intensity across all blank samples for each metabolite
df_pos_with_blank_avg <- df_pos_filtered_strict %>%
  rowwise() %>%
  mutate(Average_Blank = mean(c_across(starts_with("Blank_")), na.rm = TRUE)) %>%
  ungroup()

# Step 2: Identify sample columns to correct (excluding metadata, blanks, pools)
sample_columns <- setdiff(names(df_pos_with_blank_avg), 
                          c("Metabolite.name", "Adduct.type", "Average.Rt.min.",
                            grep("^Blank_", names(df_pos_with_blank_avg), value = TRUE),
                            grep("^Pool_", names(df_pos_with_blank_avg), value = TRUE),
                            "Average_Blank"))

# Step 3: Subtract average blank and set negative values to 0
df_pos_bg_corrected <- df_pos_with_blank_avg %>%
  mutate(across(all_of(sample_columns), 
                ~ ifelse((. - Average_Blank) < 0, 0, . - Average_Blank))) %>%
  select(-starts_with("Blank_"), -starts_with("Pool_"), -Average_Blank)
  
  
  
  
df_pos_filtered_80 <- df_pos_bg_corrected %>%
  mutate(Zero_Percentage = rowMeans(select(., -Metabolite.name, -Adduct.type, -Average.Rt.min.) == 0, na.rm = TRUE)) %>%
  filter(Zero_Percentage <= 0.20) %>%
  select(-Zero_Percentage)

# Step 2: Prepare data for heatmap
pos_heatmap_data <- df_pos_filtered_80 %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.)

# Step 3: Convert to binary matrix: 0 stays 0, all others become 1
pos_binary_heatmap_data <- pos_heatmap_data %>%
  mutate(across(everything(), ~ ifelse(. == 0, 0, 1))) %>%
  as.matrix()

# Step 4: Define color palette and breaks
color_palette <- c("black", "lightblue")  # Black = 0, Light blue = non-zero
breaks <- c(-0.5, 0.5, 1.5)  # For binary heatmap

# Step 5: Plot the heatmap
pheatmap(
  pos_binary_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,
  breaks = breaks,
  main = "Two-Color Heatmap of Zero Patterns in df_filtered_80",
  show_rownames = FALSE,
  show_colnames = TRUE
)




# Perform TIC normalization using the method = "sum"
pos_norm_TIC <- normalization(df_pos_filtered_80 %>% select(-Metabolite.name, -Adduct.type, -Average.Rt.min.), method = "sum")$newXtrain

# Combine normalized data with metadata
df_pos_TIC_normalized <- cbind(
  df_pos_filtered_80 %>% select(Metabolite.name, Adduct.type , Average.Rt.min.),
  as.data.frame(pos_norm_TIC)
)


df_pos_normalized <- df_pos_TIC_normalized

write.csv( df_pos_normalized, "df_pos_normalized.csv")

# Check for any negative values across numeric columns
any_negative <- any(df_pos_normalized  %>% 
                      select(where(is.numeric)) %>% 
                      unlist() < 0, na.rm = TRUE)

# Print result
if (any_negative) {
  print("There are negative values in the dataframe.")
} else {
  print("No negative values found in the dataframe.")
}


################################################################################
# Negative ion mode MS data
################################################################################

#---- Data Loading ----


df_neg <- read.csv("Data/NEG_Height_1_2024_09_25_03_48_22.csv", stringsAsFactors = FALSE)


#---- Run Quality Check ----
df_neg_run <- df_neg %>%
  filter(Metabolite.name != "Unknown",
         !grepl("no MS2:", Metabolite.name),
         MS.MS.matched != "FALSE")

colnames(df_neg_run)


# Creating a new dataframe with the specified columns

df_neg_run_selected <- df_neg_run %>%
  select(Average.Rt.min., Adduct.type, Metabolite.name, 
         contains('dda'))
		 
# Replace 0 with NA in numeric/intensity columns (excluding metadata)
df_neg_run_selected[ , 4:ncol(df_neg_run_selected)] <- lapply(
  df_neg_run_selected[ , 4:ncol(df_neg_run_selected)], 
  function(x) ifelse(x == 0, NA, x)
)		 
	 

# Verify the selected column names

colnames(df_neg_run_selected)

# Rename samples using lookup table
sample_names <- read.csv("Data/Sample_Names.csv", header = FALSE, stringsAsFactors = FALSE)
names(sample_names) <- c("Original", "Renamed")
rename_lookup <- setNames(sample_names$Renamed, sample_names$Original)
df_neg_run_selected <- df_neg_run_selected %>% rename_with(~ rename_lookup[gsub("_sb.cn.neg_dda1", "", .x)], contains("sample"))


# Explicitly rename blank and pool columns, if needed

df_neg_run_selected <- df_neg_run_selected %>%
  rename(
    Blank_1 = X240830_1137_hydrolysates_blank_sb.cn.neg_dda1,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.neg_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.neg_dda2
  )

write.csv(df_neg_run_selected, "df_neg_run_selected.csv" ,row.names=FALSE) 

colnames(df_neg_run_selected)

replicate_cols_neg <- colnames(df_neg_run_selected)[
  !(colnames(df_neg_run_selected) %in% c("Average.Rt.min.", "Adduct.type", "Metabolite.name"))
]
# Reshape to long format
df_long <- df_neg_run_selected %>%
  select(Average.Rt.min., Adduct.type, Metabolite.name, all_of(replicate_cols_neg)) %>%
  pivot_longer(cols = all_of(replicate_cols_neg), names_to = "Sample", values_to = "Area") %>%
  mutate(
    Group = str_remove(Sample, "_\\d+$"),  # Extract base group name
    IsMissing = ifelse(is.na(Area) | Area == 0, 1, 0)  # Mark missing or zero values
  )

# Count missing values per metabolite per group
missing_counts <- df_long %>%
  group_by(Metabolite.name, Average.Rt.min., Adduct.type, Group) %>%
  summarise(MissingReplicates = sum(IsMissing), .groups = "drop")

# Step 1: Identify metabolites where all groups have <2 missing replicates
metabolites_to_keep <- missing_counts %>%
  group_by(Metabolite.name, Average.Rt.min., Adduct.type) %>%
  summarise(MaxMissing = max(MissingReplicates), .groups = "drop") %>%
  filter(MaxMissing < 2)

# Step 2: Filter the original dataset
df_neg_filtered_strict <- df_neg_run_selected %>%
  semi_join(metabolites_to_keep, by = c("Metabolite.name", "Average.Rt.min.", "Adduct.type"))
  
  
  
# Step 1: Calculate average intensity across all blank samples for each metabolite
df_neg_with_blank_avg <- df_neg_filtered_strict %>%
  rowwise() %>%
  mutate(Average_Blank = mean(c_across(starts_with("Blank_")), na.rm = TRUE)) %>%
  ungroup()

# Step 2: Identify sample columns to correct (excluding metadata, blanks, pools)
sample_columns <- setdiff(names(df_neg_with_blank_avg), 
                          c("Metabolite.name", "Adduct.type", "Average.Rt.min.",
                            grep("^Blank_", names(df_neg_with_blank_avg), value = TRUE),
                            grep("^Pool_", names(df_neg_with_blank_avg), value = TRUE),
                            "Average_Blank"))

# Step 3: Subtract average blank and set negative values to 0
df_neg_bg_corrected <- df_neg_with_blank_avg %>%
  mutate(across(all_of(sample_columns), 
                ~ ifelse((. - Average_Blank) < 0, 0, . - Average_Blank))) %>%
  select(-starts_with("Blank_"), -starts_with("Pool_"), -Average_Blank)
  
  
  
  
df_neg_filtered_80 <- df_neg_bg_corrected %>%
  mutate(Zero_Percentage = rowMeans(select(., -Metabolite.name, -Adduct.type, -Average.Rt.min.) == 0, na.rm = TRUE)) %>%
  filter(Zero_Percentage <= 0.20) %>%
  select(-Zero_Percentage)

# Step 2: Prepare data for heatmap
neg_heatmap_data <- df_neg_filtered_80 %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.)

# Step 3: Convert to binary matrix: 0 stays 0, all others become 1
neg_binary_heatmap_data <- neg_heatmap_data %>%
  mutate(across(everything(), ~ ifelse(. == 0, 0, 1))) %>%
  as.matrix()

# Step 4: Define color palette and breaks
color_palette <- c("black", "lightblue")  # Black = 0, Light blue = non-zero
breaks <- c(-0.5, 0.5, 1.5)  # For binary heatmap

# Step 5: Plot the heatmap
pheatmap(
  neg_binary_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,
  breaks = breaks,
  main = "Two-Color Heatmap of Zero Patterns in df_filtered_80",
  show_rownames = FALSE,
  show_colnames = TRUE
)




# Perform TIC normalization using the method = "sum"
neg_norm_TIC <- normalization(df_neg_filtered_80 %>% select(-Metabolite.name, -Adduct.type, -Average.Rt.min.), method = "sum")$newXtrain

# Combine normalized data with metadata
df_neg_TIC_normalized <- cbind(
  df_neg_filtered_80 %>% select(Metabolite.name, Adduct.type , Average.Rt.min.),
  as.data.frame(neg_norm_TIC)
)


df_neg_normalized <- df_neg_TIC_normalized

write.csv( df_neg_normalized, "df_neg_normalized.csv")

# Check for any negative values across numeric columns
any_negative <- any(df_neg_normalized  %>% 
                      select(where(is.numeric)) %>% 
                      unlist() < 0, na.rm = TRUE)

# Print result
if (any_negative) {
  print("There are negative values in the dataframe.")
} else {
  print("No negative values found in the dataframe.")
}




###############################################################################

# Combine the two data frames
df_combined <- rbind(df_neg_normalized, df_pos_normalized)

# Check the combined data frame
print(dim(df_combined))  # Check dimensions (rows and columns)
head(df_combined)        # Preview the first few rows

write.csv(df_combined, "df_combined_POS_NEG.csv")

###############################################################################
