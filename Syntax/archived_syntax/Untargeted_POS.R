################################################################################
#Libraries Loading ----

library(tidyverse)
library(pheatmap)
# library(vsn)
library(patchwork)
library(naniar)
# library(simputation)
# library(missForest)
library(KODAMA)
library(viridis)
library(ggrepel)


################################################################################
#Data Loading ----

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
         starts_with("X240830_1137_hydrolysates_blank"),
         starts_with("X240830_1137_hydrolysates_pool"),
         starts_with("sample") & ends_with("_sb.cn.pos_dda1"))

# Verify the selected column names

colnames(df_pos_run_selected)


#Rename Samples
# Create a lookup table for renaming columns

sample_names <- read.csv("Data/Sample_Names.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(sample_names) <- c("Original", "Renamed")  # Assign meaningful column names

# Create a lookup table for renaming

rename_lookup <- setNames(sample_names$Renamed, sample_names$Original)


# Rename columns in df_run_selected using the lookup table

df_pos_run_selected <- df_pos_run_selected %>%
  rename_with(~ rename_lookup[gsub("_sb.cn.pos_dda1", "", .x)], .cols = contains("sample"))

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

# Trim any leading or trailing spaces in column names

colnames(df_pos_run_selected) <- trimws(colnames(df_pos_run_selected))

# Verify renamed column names

colnames(df_pos_run_selected)
write.csv(df_pos_run_selected, "df_pos_run_selected.csv" ) 



# Converting selected columns to long format

df_pos_run_selected_long <- df_pos_run_selected %>%
  pivot_longer(
    cols = -c(Average.Rt.min., Adduct.type, Metabolite.name), # Exclude metadata columns
    names_to = "Samples",
    values_to = "Area"
  )


### Raw Data distribution ----

# Create the box plot

ggplot(df_pos_run_selected_long, aes(x = Samples, y = Area)) +
  geom_boxplot() +
  geom_jitter(width=0.15)+
  theme(axis.text.x = element_text(angle = 90))+  
  labs(x = "Samples", y = "Area", title = "Box Plot of Area by Sample")

# Create the box plot with log-transformed data

ggplot(df_pos_run_selected_long, aes(x = Samples, y = log(Area))) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +  # Added alpha for better visualization of points
  theme(axis.text.x = element_text(angle = 90))+  
  labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log-Transformed Area by Sample")


### PCA ----
# 1. Prepare data for PCA
# Select only the columns representing samples and exclude metadata columns

df_pos_pca_data <- df_pos_run_selected %>%
  select(-Average.Rt.min., -Adduct.type, -Metabolite.name) %>%
  t() %>%
  as.data.frame()

# Standardize the data (mean=0, variance=1) for PCA

df_pos_pca_data_scaled <- scale(df_pos_pca_data)

# 2. Perform PCA

pca_pos_result <- prcomp(df_pos_pca_data_scaled, center = TRUE, scale. = TRUE)

# 3. Extract sample names and group by ignoring replicate numbers

sample_names <- rownames(df_pos_pca_data)
group_names <- str_extract(sample_names, "^[^_]+")  # Extract the part before "_"

# 4. Create a data frame with PCA results and group labels

pca_pos_data <- data.frame(Sample = sample_names, 
                       Group = group_names,
                       PC1 = pca_pos_result$x[,1], 
                       PC2 = pca_pos_result$x[,2])

# 5. Plot the PCA results with colors based on sample groups

ggplot(pca_pos_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of Samples by Group", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "right")



################################################################################
### Missing values ----

#### Missing vlues % ----

# Calculate missing value percentage for each sample
pos_missing_percentage <- df_pos_run_selected %>%
  select(-Average.Rt.min., -Adduct.type, -Metabolite.name) %>%  # Exclude metadata columns
  summarise(across(everything(), ~ mean(is.na(.) | . == 0) * 100)) %>%  # Calculate percentage of missing values
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Missing_Percentage")  # Convert to long format

# Plot missing value percentage as a bar chart

ggplot(pos_missing_percentage, aes(x = Sample, y = Missing_Percentage)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Missing Value Percentage by Sample",
       x = "Sample Names",
       y = "Missing Value Percentage (%)")


################################################################################
# Step 1: Identify sample columns and separate group names from replicate numbers

# Define the sample columns (exclude metadata columns like "Average.Rt.min." and "Metabolite.name")

pos_sample_mv_columns <- names(df_pos_run_selected)[!(names(df_pos_run_selected) %in% c("Average.Rt.min.", "Adduct.type", "Metabolite.name", "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2"))]

# Convert the data to long format to make it easier to work with groups and replicates

df_pos_mv_handle_long <- df_pos_run_selected %>%
  pivot_longer(cols = all_of(pos_sample_mv_columns),
               names_to = "Sample",
               values_to = "Area") %>%
  mutate(Group = sub("_\\d+$", "", Sample),       # Extract the base sample group name (everything before the last underscore and digit)
         Replicate = as.numeric(sub(".*_", "", Sample))  # Extract replicate number (the last number after the underscore)
  )

# Step 2: Count the non-missing values in each group for each metabolite

# Group by metabolite and group, and count non-missing values in each group

df_pos_group_counts <- df_pos_mv_handle_long %>%
  group_by(Metabolite.name, Group) %>%
  summarise(NonMissingCount = sum(!is.na(Area)), .groups = "drop")

# Step 3: Filter metabolites that meet the threshold (at least 2 non-missing values in any group)

# Find metabolites with at least 2 replicates detected in any group

pos_metabolites_to_keep <- df_pos_group_counts %>%
  filter(NonMissingCount >= 2) %>%
  pull(Metabolite.name) %>%
  unique()  # Unique list of metabolites that meet the criteria

# Step 4: Filter the original df_run_selected based on the metabolites_to_keep list

df_pos_filtered <- df_pos_run_selected %>%
  filter(Metabolite.name %in% pos_metabolites_to_keep)


# Step 5: Count the number of identified metabolites (non-missing values) per sample
# Select sample columns and count non-NA values per column

pos_metabolite_counts <- df_pos_filtered %>%
  select(-Average.Rt.min., -Adduct.type, -Metabolite.name) %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Count")

# Step 6: Plot a bar plot for the number of identified metabolites per sample

ggplot(pos_metabolite_counts, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Number of Identified Metabolites per Sample", x = "Sample", y = "Count of Identified Metabolites")



# missing value visualization

# Step 1: Create a binary matrix where 1 represents missing values (0 in this case) and 0 represents present values

pos_missing_matrix <- df_pos_filtered %>%
  select(-Average.Rt.min., -Adduct.type, -Metabolite.name) %>%
  mutate(across(everything(), ~ ifelse(. == 0, 1, 0))) %>%  # Mark 0 as missing (1), others as present (0)
  as.matrix()

# Double-check the binary matrix for any missing values represented as 1
table(pos_missing_matrix)


# Define color and breaks explicitly for binary data
color_palette <- c("white", "black")
breaks <- c(-0.5, 0.5, 1.5)

# Plot the heatmap
pheatmap(pos_missing_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = color_palette,
         breaks = breaks,
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "Heatmap of Missing Values Across Metabolites and Samples")


# Missing value identification using library(naniar) 

# Convert 0 values to NA in the dataset temporarily for missing data analysis

df_pos_temp <- df_pos_filtered %>%
  mutate(across(everything(), ~ ifelse(. == 0, NA, .)))

vis_miss(df_pos_temp, cluster = TRUE)  # Visualize clustered missing values

# Upset plot for missing data combinations
gg_miss_upset(df_pos_temp)

# 
# # Count the occurrences of each type of "missing" value in df_run_selected
# pos_missing_summary <- df_pos_run_selected %>%
#   summarise(
#     Total_Zero = sum(across(where(is.numeric), ~ sum(. == 0, na.rm = TRUE))),
#     Total_NA = sum(is.na(.)),
#     Total_NaN = sum(across(where(is.numeric), ~ sum(is.nan(.))))
#   )
# 
# # Display the summary of missing values
# print(pos_missing_summary)
# 
# # Plot missingness percentages by groups
# ggplot(pos_missing_summary, aes(x = Sample, y = Missing_Percentage, fill = Type)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Missingness by Group and Sample",
#        x = "Sample Names",
#        y = "Missing Value Percentage (%)")


################################################################################
# DDA; these are random observation; not everything was quantified
################################################################################

# If blanks have high intensities for certain metabolites, subtracting them first avoids propagating those effects through normalization.
# When using scaling techniques (e.g., Pareto scaling or mean centering) for multivariate analysis (e.g., PCA), it is usually better to subtract blanks first.



# Step 1: Calculate Average Blank Values

# Calculate the average of blank intensities
df_pos_with_blank_avg <- df_pos_filtered %>%
  rowwise() %>%
  mutate(Average_Blank = mean(c_across(starts_with("Blank_")), na.rm = TRUE))


# Step 2: Subtract Blank Values from Samples
# Define the sample columns (exclude metadata and blank columns)
sample_columns <- setdiff(names(df_pos_with_blank_avg), 
                          c("Metabolite.name", "Adduct.type", "Average.Rt.min.", "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2", "Average_Blank"))

## Subtract blanks and replace negatives with NA or 0
## Negative values can be replaced with NA to treat them as missing data for downstream analyses or with zero if they signify an absence of signal rather than an artifact.

# # Subtract blanks and replace negatives with NA
# df_blank_corrected <- df_with_blank_avg %>%
#   mutate(across(all_of(sample_columns), ~ ifelse((. - Average_Blank) < 0, NA, . - Average_Blank))) %>%
#   select(-starts_with("Blank_"), -starts_with("Pool_"), -Average_Blank)  # Drop Blank, Pool, and Average_Blank columns

# Subtract blanks and replace negatives with zero
df_pos_blank_corrected <- df_pos_with_blank_avg %>%
  mutate(across(all_of(sample_columns), ~ ifelse((. - Average_Blank) < 0, 0, . - Average_Blank))) %>%
  select(-starts_with("Blank_"), -starts_with("Pool_"), -Average_Blank)  # Drop Blank, Pool, and Average_Blank columns

################################################################################
# Check for any negative values across numeric columns
any_negative <- any(df_pos_blank_corrected %>% 
                      select(where(is.numeric)) %>% 
                      unlist() < 0, na.rm = TRUE)

# Print result
if (any_negative) {
  print("There are negative values in the dataframe.")
} else {
  print("No negative values found in the dataframe.")
}
################################################################################

# Identify metabolites with >20% zero values in samples
pos_metabolites_to_keep <- df_pos_blank_corrected %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.) %>%
  summarise(across(everything(), ~ mean(. == 0, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(),
               names_to = "Sample",
               values_to = "Zero_Percentage") %>%
  filter(Zero_Percentage <= 0.20) %>%
  pull(Sample)

# Filter the dataset to keep only those metabolites
df_pos_filtered_80 <- df_pos_blank_corrected %>%
  select(Metabolite.name,
         Adduct.type,
         Average.Rt.min.,
         all_of(pos_metabolites_to_keep))


################################################################################
## Visualize missing values
# Prepare the dataset for heatmap
# Remove metadata columns (Metabolite.name, Average.Rt.min.)
pos_heatmap_data <- df_pos_filtered_80 %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.)

# Create a binary-like matrix for coloring: 0 stays as 0, all other values become 1
pos_binary_heatmap_data <- pos_heatmap_data %>%
  mutate(across(everything(), ~ ifelse(. == 0, 0, 1))) %>%
  as.matrix()

# Define color palette and breaks for binary data
color_palette <- c("black", "lightblue") # Black for 0, light blue for all other values
breaks <- c(-0.5, 0.5, 1.5) # Ensure two distinct bins: 0 and non-zero values

# Plot the heatmap
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



################################################################################

# Follow the recent protocol: Pakkir Shah, Abzer K., et al. "Statistical analysis of feature-based molecular networking results from non-targeted metabolomics data." Nature protocols (2024): 1-71. https://www.nature.com/articles/s41596-024-01046-3

# Many feature extraction software programs, such as MZmine 3, often generate tables with missing values denoted as ‘NA’, ‘NaN’ or 0. This means that for several m/z and RT traces in a given sample, there may not be a peak detected and therefore no value is available. However, many statistical approaches, such as PCA, require numerical values for each observation. Hence, these features with missing values need to be removed or imputed. In this section, we handle the zero values in our blank-removed feature quantification table. https://www.nature.com/articles/s41596-024-01046-3

# To handle log transformation when dataset contains zeros, we can use one of the following approaches to avoid undefined values (since log ⁡ ( 0 ) log(0) is undefined):


# Step 1: Identify the second smallest non-zero value
# Flatten all numeric data, filter out zeros, and get the second smallest value
pos_non_zero_values <- df_pos_filtered_80 %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  .[. > 0] %>%           # Keep only non-zero values
  sort()

pos_second_min_value <- unique(pos_non_zero_values)[2]  # The second smallest non-zero value

# Print the second smallest non-zero value
print(paste("The second smallest non-zero value is:", pos_second_min_value))

# Step 2: Replace zeros with random values between 1 and the second smallest non-zero value
set.seed(141222)  # Set seed for reproducibility

df_pos_zero_imputed <- df_pos_filtered_80 %>%
  mutate(across(where(is.numeric), ~ ifelse(. == 0,                                             runif(1, min = 1, max = pos_second_min_value), .)))

# Step 3: Verify the replacement of zeros
# Count zeros before and after imputation
pos_zeros_before <- sum(df_pos_filtered_80 == 0, na.rm = TRUE)
pos_zeros_after <- sum(df_pos_zero_imputed == 0, na.rm = TRUE)

print(paste("Number of zeros before imputation:", pos_zeros_before))
print(paste("Number of zeros after imputation:", pos_zeros_after))


################################################################################
# Normalization

# Perform TIC normalization using the method = "sum"
pos_norm_TIC <- normalization(df_pos_zero_imputed %>% select(-Metabolite.name, -Adduct.type, -Average.Rt.min.), method = "sum")$newXtrain

# Combine normalized data with metadata
df_pos_TIC_normalized <- cbind(
  df_pos_zero_imputed %>% select(Metabolite.name, Adduct.type , Average.Rt.min.),
  as.data.frame(pos_norm_TIC)
)

# Verify the sums of each sample
pos_TIC_after <- colSums(df_pos_TIC_normalized %>% select(-Metabolite.name, -Adduct.type, -Average.Rt.min.))
print("TIC values after normalization (should be ~1):")
print(pos_TIC_after)

# Check for NAs
print(paste("Number of NA values in Normalized data:", sum(is.na(pos_norm_TIC) == TRUE)))


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
# visualization
# Prepare data for boxplot / Log2
# Prepare data for plotting (before imputation)
df_pos_long_before <- df_pos_filtered_80 %>%
  pivot_longer(cols = -c(Metabolite.name, Adduct.type, Average.Rt.min.), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Stage = "Before Imputation", Log2_Intensity = log2(Intensity + 1))  # Add log2 transformation


# Prepare data for plotting (after imputation)
df_pos_long_after_impute <- df_pos_zero_imputed %>%
  pivot_longer(cols = -c(Metabolite.name, Adduct.type, Average.Rt.min.), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Stage = "After Imputation", Log2_Intensity = log2(Intensity + 1))  # Add log2 transformation

# Prepare data for plotting (after normalization)
df_pos_long_after_norm <- df_pos_normalized %>%
  pivot_longer(cols = -c(Metabolite.name, Adduct.type, Average.Rt.min.), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Stage = "After Normalization", Log2_Intensity = log2(Intensity + 1))  # Add log2 transformation


# Boxplot for before imputation
plot_before <- ggplot(df_pos_long_before, aes(x = Sample, y = Log2_Intensity)) +
  geom_boxplot(outlier.alpha = 0.5, fill = "lightblue") +
  labs(title = "Before Imputation", x = "Samples", y = "Log2(Intensity)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels
  )

# Boxplot for after imputation
plot_after_impute <- ggplot(df_pos_long_after_impute, aes(x = Sample, y = Log2_Intensity)) +
  geom_boxplot(outlier.alpha = 0.5, fill = "lightcoral") +
  labs(title = "After Imputation", x = "Samples", y = "Log2(Intensity)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels
  )

# Boxplot for after normalization
plot_after_norm <- ggplot(df_pos_long_after_norm, aes(x = Sample, y = Log2_Intensity)) +
  geom_boxplot(outlier.alpha = 0.5, fill = "lightgreen") +
  labs(title = "After Normalization", x = "Samples", y = "Log2(Intensity)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels
  )

# Combine the three plots using patchwork
combined_plot <- plot_before +plot_after_impute + plot_after_norm +
  plot_layout(ncol = 1) & theme(legend.position = "none")

# Display the combined plot
print(combined_plot)

################################################################################
## Let's combined pos and neg

## Then follow downstream

################################################################################
## Scaling

# Select numeric data (exclude metadata)
df_pos_scaled <- df_pos_normalized %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.) %>% 
  scale() %>%
  as.data.frame()

# Combine scaled data with metadata
df_pos_scaled <- cbind(
  df_pos_normalized %>% select(Metabolite.name, Adduct.type, Average.Rt.min.),
  df_pos_scaled
)

# Check scaling results
summary(df_pos_scaled %>% select(-Metabolite.name, -Adduct.type, -Average.Rt.min.))

################################################################################
## PCA

# Prepare data for PCA (exclude metadata)
pca_data <- df_pos_scaled %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.) %>%
  t() %>%
  as.data.frame()

# Assign rownames to the PCA data
rownames(pca_data) <- colnames(df_pos_normalized %>% select(-Metabolite.name,-Adduct.type, -Average.Rt.min.))

# Perform PCA on scaled data
pca_result <- prcomp(pca_data, center = FALSE, scale. = FALSE)

# Create a data frame for plotting PCA results
pca_plot_data <- as.data.frame(pca_result$x[, 1:2])  # Extract PC1 and PC2
pca_plot_data$Group <- sapply(rownames(pca_data), function(x) strsplit(x, "_")[[1]][1])
pca_plot_data$Sample <- rownames(pca_data)

# Predefined solid colors for 9 groups
solid_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"
)

# Dynamically map colors to groups
group_colors <- setNames(solid_colors, unique(pca_plot_data$Group))

# Step 1: Calculate variance explained by each PC
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(explained_variance[1] * 100, 2)  # PC1 percentage
pc2_var <- round(explained_variance[2] * 100, 2)  # PC2 percentage

# Step 2: Create the PCA plot with PC values in axis labels

pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.5) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PCA Plot with Percent Variance Explained",
    x = paste0("Principal Component 1 (", pc1_var, "%)"),
    y = paste0("Principal Component 2 (", pc2_var, "%)"),
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

# Print the PCA plot
print(pca_plot)

# Save the PCA plot as PDF and JPEG
ggsave("PCA_Plot_Custom_Colors.pdf", plot = pca_plot, device = "pdf", width = 10, height = 8, dpi = 300)
ggsave("PCA_Plot_Custom_Colors.jpeg", plot = pca_plot, device = "jpeg", width = 10, height = 8, dpi = 300)


################################################################################
## Heat Map

# Extract numeric data for heatmap
heatmap_data <- df_scaled %>% select(-Metabolite.name, -Average.Rt.min.) %>% as.matrix()

color_palette <- colorRampPalette(c( "#A6CEE3" , "#33A02C", "#E31A1C"))(50)

pheatmap(
  heatmap_data,
  cluster_rows = TRUE,    # Perform hierarchical clustering on rows (metabolites).
  cluster_cols = TRUE,    # Perform hierarchical clustering on columns (samples).
  color = color_palette,  # Apply the defined color palette.
  main = "Heatmap of Scaled Data",  # Title of the heatmap.
  show_rownames = FALSE,  # Optional: Hide row names (metabolites) to improve visualization.
  show_colnames = TRUE    # Optional: Show column names (sample names).
)



################################################################################
# Summarize multiple copies of metabolites coming from different libraries

# Step 1: Create a new column 'Metabolite.name.mod' by extracting text before the first ';'
df_scaled_v2 <- df_scaled %>%
  mutate(Metabolite.name.mod = str_extract(Metabolite.name, "^[^;]+"))

# Step 2: Group by 'Metabolite.name.mod' and calculate the average for all sample columns
df_scaled_v2_averaged <- df_scaled_v2 %>%
  group_by(Metabolite.name.mod) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()


## Heat Map

# Extract numeric data for heatmap
heatmap_data_v2 <- df_scaled_v2_averaged %>% select(-Metabolite.name.mod, -Average.Rt.min.) %>% as.matrix()

color_palette <- colorRampPalette(c( "#A6CEE3" , "#33A02C", "#E31A1C"))(50)

pheatmap(
  heatmap_data_v2,
  cluster_rows = TRUE,    # Perform hierarchical clustering on rows (metabolites).
  cluster_cols = TRUE,    # Perform hierarchical clustering on columns (samples).
  color = color_palette,  # Apply the defined color palette.
  main = "Heatmap of Scaled Data",  # Title of the heatmap.
  show_rownames = FALSE,  # Optional: Hide row names (metabolites) to improve visualization.
  show_colnames = TRUE    # Optional: Show column names (sample names).
)



################################################################################
# Find Top metabolites

# Summarize multiple copies of metabolites coming from different libraries

# Step 1: Create a new column 'Metabolite.name.mod' by extracting text before the first ';'
df_normalized_v2 <- df_normalized %>%
  mutate(Metabolite.name.mod = str_extract(Metabolite.name, "^[^;]+"))

# Step 2: Group by 'Metabolite.name.mod' and calculate the average for all sample columns
df_normalized_v2_averaged <- df_normalized_v2 %>%
  group_by(Metabolite.name.mod) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()


# Calculate the Average Intensity Across All Samples

# Exclude metadata columns
metabolite_data <- df_normalized_v2_averaged %>% 
  select(-Metabolite.name.mod, -Average.Rt.min.)


df_normalized_v2_averaged <- df_normalized_v2_averaged %>% 
  mutate(Median_Intensity = apply(metabolite_data, 1, median, na.rm = TRUE))

# Sort and select the top metabolites by median intensity
highest_metabolites <- df_normalized_v2_averaged %>%
  arrange(desc(Median_Intensity)) %>%
  select(Metabolite.name.mod, Median_Intensity)

# Print the top metabolites by median intensity
print(highest_metabolites)


write.csv(df_normalized_v2_averaged, "df_normalized_v2_averaged.csv")


################################################################################

#  find the top metabolites for each group

# Step 1: Extract group names from sample columns
sample_columns <- colnames(df_normalized_v2_averaged)[3:(ncol(df_normalized_v2_averaged) - 1)] # Exclude metadata columns
group_names <- sapply(sample_columns, function(x) strsplit(x, "_")[[1]][1])

# Step 2: Rename columns with group information
colnames(df_normalized_v2_averaged)[3:(ncol(df_normalized_v2_averaged) - 1)] <- group_names

# Step 3: Reshape the dataframe to long format
df_normalized_v2_averaged_long <- df_normalized_v2_averaged %>%
  pivot_longer(
    cols = -c(Metabolite.name.mod, Average.Rt.min., Median_Intensity), 
    names_to = "Group", 
    values_to = "Intensity"
  )

# Step 4: Group by Metabolite.name.mod and Group, calculate mean intensity
group_summary <- df_normalized_v2_averaged_long %>%
  group_by(Group, Metabolite.name.mod) %>%
  summarise(Mean_Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Step 5: For each group, find the top 25 metabolites
top_25_metabolites <- group_summary %>%
  group_by(Group) %>%
  slice_max(order_by = Mean_Intensity, n = 50) %>%
  arrange(Group, desc(Mean_Intensity))

# Step 6: View results
print(top_25_metabolites)

# Optional: Save the results to a CSV file
write.csv(top_25_metabolites, "Top_25_Metabolites_Per_Group.csv", row.names = FALSE)


################################################################################

# Create a presence/absence matrix for metabolites across groups, dropping Mean_Intensity
presence_matrix <- top_25_metabolites %>%
  select(Group, Metabolite.name.mod) %>%  # Only keep the relevant columns
  mutate(Presence = 1) %>%  # Add a column for presence
  pivot_wider(
    names_from = Group,  # Groups as columns
    values_from = Presence,  # Use presence/absence
    values_fill = 0  # Fill missing values with 0
  )

# Check the resulting matrix
print(presence_matrix)


# Convert to a matrix for visualization
presence_matrix_matrix <- presence_matrix %>%
  select(-Metabolite.name.mod) %>%
  as.matrix()
rownames(presence_matrix_matrix) <- presence_matrix$Metabolite.name.mod



library(ComplexUpset)

# Prepare data for UpSet plot
presence_matrix_for_upset <- presence_matrix %>%
  pivot_longer(-Metabolite.name.mod, names_to = "Group", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  mutate(Value = TRUE) %>%
  select(-Presence) %>%
  pivot_wider(names_from = Group, values_from = Value, values_fill = FALSE)

# Plot the UpSet
upset(
  presence_matrix_for_upset,
  colnames(presence_matrix_for_upset)[-1],
  name = "Metabolite Groups",
  base_annotations = list(
    'Intersection size' = intersection_size(counts = TRUE)
  )
)


# Visualize the presence/absence matrix using a heatmap
pheatmap(
  presence_matrix_matrix,
  cluster_rows = TRUE,  # Cluster metabolites
  cluster_cols = TRUE,  # Cluster groups
  color = colorRampPalette(c("white", "darkred"))(50),
  main = "Presence of Top 25 Metabolites Across Groups",
  show_rownames = FALSE,  # Hide row names
  show_colnames = TRUE    # Show column (group) names
)


################################################################################







