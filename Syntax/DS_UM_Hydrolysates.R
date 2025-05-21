library(tidyverse)   
library(ggplot2)
library(RColorBrewer)  
library(ggrepel)    
library(pheatmap)    
library(ComplexUpset)


## For upstream analysis syntax refer to "Untargted_POS_NEG.R" file

df_combined <- read.csv("Data/df_combined_POS_NEG.csv") 
Met_Category <- read.csv("Data/HMDB_metabolites.csv")

# Remove the 'X' column
df_combined <- df_combined %>% select(-X)

str(df_combined)
colnames(df_combined)

str(Met_Category)
colnames(Met_Category)

# Clean 'Metabolite.name' by removing extra text after the first semicolon
df_combined <- df_combined %>%
  mutate(Metabolite.name = str_trim(str_remove(Metabolite.name, ";.*")))

df_combined <- df_combined %>%
  mutate(Metabolite.name = case_when(
    Metabolite.name == "Inosine-5-monophosphate" ~ "Inosine-5'-monophosphate",
    Metabolite.name == "2'-Deoxyadenosine-5'-monophosphate" ~ "2'-Deoxyadenosine 5'-monophosphate",
    Metabolite.name == "(-)-Riboflavin" ~ "Riboflavin",
    Metabolite.name == "D-(+)-Pantothenic Acid" ~ "Pantothenic acid",
    Metabolite.name == "2'-Deoxyadenosine-5'-Monophosphate" ~ "2'-Deoxyadenosine 5'-Monophosphate",
    Metabolite.name == "Aconitic Acid (Not Validated, Isomer Of 271)" ~ "Aconitic Acid",
    Metabolite.name == "Adenosine-3-Monophosphate" ~ "Adenosine 3'-Monophosphate",
    Metabolite.name == "Azelaic Acid (Not Validated)" ~ "Azelaic Acid",
    Metabolite.name == "Citric Acid (Not Validated, Isomer Of 227)" ~ "Citric Acid",
    Metabolite.name == "Flavin??Adenine??Dinucleotide" ~ "Flavin Adenine Dinucleotide",
    Metabolite.name == "Inosine-5-Monophosphate" ~ "Inosine-5'-Monophosphate",
    Metabolite.name == "Isoleucine (Not Validated)" ~ "Isoleucine",
    Metabolite.name == "Xanthosine (Not Validated)" ~ "Xanthosine",
    TRUE ~ Metabolite.name
  ))


# Remove unwanted columns
df_filtered <- df_combined %>%
  select(-Adduct.type, -Average.Rt.min.)

# Average all numeric columns by Metabolite.name
df_avg <- df_filtered %>%
  group_by(Metabolite.name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# # Join averaged data with metabolite categories
# df_avg <- df_avg %>%
#   left_join(Met_Category, by = c("Metabolite.name" = "Metabolite"))
# 
# unmatched_names <- df_avg %>%
#   filter(is.na(Category)) %>%
#   distinct(Metabolite.name) %>%
#   pull(Metabolite.name)





join_with_category <- function(df_avg, met_category_df,
                               avg_col = "Metabolite.name", cat_col = "Metabolite") {
  df_avg <- df_avg %>%
    mutate(join_key = tolower(str_trim(!!sym(avg_col))))  # create join key from metabolite name
  
  met_category_df <- met_category_df %>%
    mutate(join_key = tolower(str_trim(!!sym(cat_col))))  # create join key from reference table
  
  df_joined <- df_avg %>%
    left_join(met_category_df %>% select(join_key, Category), by = "join_key") %>%
    select(-join_key)
  
  return(df_joined)
}



df_avg <- join_with_category(df_avg, Met_Category)


df_avg %>%
  filter(is.na(Category)) %>%
  distinct(Metabolite.name)


str(df_avg)
### Important

# This file has normalized data.
# For PCA and Heatmap scale the data as below

###############################################################################
## Scaling

# Matrix of intensity values only
df_scaled <- df_avg %>%
  select(-Metabolite.name, -Category) %>%
  scale() %>%
  as.data.frame()

# Add back metabolite names and category
df_scaled <- df_scaled %>%
  mutate(Metabolite.name = df_avg$Metabolite.name,
         Category = df_avg$Category) %>%
  relocate(Metabolite.name, .before = 1)



################################################################################
## PCA

# Prepare data for PCA (exclude metadata)
df_avg_scaled_matrix <- df_avg %>%
  select(-Metabolite.name, -Category) %>%
  t() %>%
  as.data.frame()

# Run PCA
pca_result <- prcomp(df_avg_scaled_matrix, center = FALSE, scale. = FALSE)

# Extract first two PCs for plotting
df_pca_plot <- as.data.frame(pca_result$x[, 1:2])
df_pca_plot$Group <- sapply(rownames(df_pca_plot), function(x) strsplit(x, "_")[[1]][1])
df_pca_plot$Sample <- rownames(df_pca_plot)

# Define colors for up to 9 groups (adjust if more)
pca_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"
)
group_color_map <- setNames(pca_colors, unique(df_pca_plot$Group))

# Calculate variance explained
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(explained_var[1] * 100, 2)
pc2_var <- round(explained_var[2] * 100, 2)

# Create PCA plot
pca_plot <- ggplot(df_pca_plot, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.5) +
  scale_color_manual(values = group_color_map) +
  labs(
    title = "PCA Plot of Scaled Metabolite Abundance",
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

# Display PCA plot
print(pca_plot)

# Save PCA plot as PDF and JPEG
ggsave("PCA_Plot_Scaled_Metabolites.pdf", plot = pca_plot, device = "pdf", width = 10, height = 8, dpi = 300)
ggsave("PCA_Plot_Scaled_Metabolites.jpeg", plot = pca_plot, device = "jpeg", width = 10, height = 8, dpi = 300)



## PCA and Loading plots


# Extract PCA loadings
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$Metabolite.name <- df_avg$Metabolite.name
loadings$Category <- df_avg$Category

# Variance explained
pc1_var <- round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 2)
pc2_var <- round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 2)

# Plot
loading_plot <- ggplot(loadings, aes(x = PC1, y = PC2, color = Category)) +
  geom_point(alpha = 0.85, size = 3) +
  labs(
    title = "PCA Loading Plot by Metabolite Category",
    x = paste0("PC1 Loadings (", pc1_var, "%)"),
    y = paste0("PC2 Loadings (", pc2_var, "%)")
  ) +
  scale_color_brewer(palette = "Set3") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

print(loading_plot)

# Save
ggsave("PCA_Loading_Plot_By_Category.pdf", plot = loading_plot, device = "pdf", width = 10, height = 8, dpi = 300)
ggsave("PCA_Loading_Plot_By_Category.jpeg", plot = loading_plot, device = "jpeg", width = 10, height = 8, dpi = 300)


################################################################################
## Heat Map


# Prepare the intensity matrix
# Rows = metabolites, Columns = samples
m.int <- df_avg %>%
  select(-Metabolite.name, -Category) %>%
  as.matrix()

# Assign rownames as metabolite names
rownames(m.int) <- df_avg$Metabolite.name

# Create row annotation (metabolite category)
row_annotation <- data.frame(Category = df_avg$Category)
rownames(row_annotation) <- rownames(m.int)

# Define category colors
# Automatically generate colors for each unique category
unique_categories <- unique(df_avg$Category)
col_cat2 <- list(Category = setNames(brewer.pal(n = length(unique_categories), name = "Set3"), unique_categories))

# Plot the heatmap
p.heat <- pheatmap(
  m.int,
  scale = "row",
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "#E31A1C"))(50),
  annotation_row = row_annotation,
  angle_col = 315,
  annotation_colors = col_cat2
)

# Save to PDF
pdf("metaboliteheat.pdf", height = 15, width = 8)
print(p.heat)
dev.off()




################################################################################
# Find Top metabolites
# Reshape df_avg to long format
df_long <- df_avg %>%
  pivot_longer(cols = -c(Metabolite.name, Category),
               names_to = "Sample", values_to = "Intensity")

# Extract sample group (everything before last "_")
df_long <- df_long %>%
  mutate(Group = str_extract(Sample, ".*(?=_[0-9]+$)"))

# Average intensity per metabolite per group
df_group_avg <- df_long %>%
  group_by(Group, Metabolite.name, Category) %>%
  summarise(Mean_Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# For each group, get top 50 metabolites
top50_metabolites <- df_group_avg %>%
  group_by(Group) %>%
  slice_max(order_by = Mean_Intensity, n = 50) %>%
  arrange(Group, desc(Mean_Intensity)) %>%
  ungroup()

# View result
print(top50_metabolites)



################################################################################

# Create binary presence matrix
upset_data <- top50_metabolites %>%
  select(Group, Metabolite.name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Group, values_from = value, values_fill = 0)

# Store metabolite names separately (optional)
metabolite_names <- upset_data$Metabolite.name

# Drop the name column to retain just the presence/absence matrix
upset_matrix <- upset_data %>%
  select(-Metabolite.name)

ComplexUpset::upset(
  data = upset_data,
  intersect = colnames(upset_matrix),
  name = "Top 50 Metabolites",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(vjust = -0.5, size = 3)
    )
  ),
  width_ratio = 0.2
)




