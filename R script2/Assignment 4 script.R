#################################################################################
# Geographic Patterns of Speciation in Oncorhynchus (Pacific Salmon & Trout)
# Testing the Allopatric Speciation Hypothesis
#
# Testing whether sister species in Oncorhynchus occupy separate
# geographic regions, consistent with allopatric speciation via vicariance
#
# Author: Mahnoor Nizamani
#
# Nov 28th 2025
#
# GitHub: https://github.com/mahnoor-nizz/BINF6210-Assignment-4
# ###############################################################################




# Load required packages -------------------------------------------------------


library(ape)           # phylogenetic analysis
library(phangorn)      # model testing and ML trees
library(DECIPHER)      # alignment with codon models
library(Biostrings)    # sequence manipulation
library(rentrez)       # NCBI data access
library(rgbif)         # GBIF occurrence data
library(sf)            # spatial data handling
library(rnaturalearth) # base maps
library(ggplot2)       # visualization
library(ggtree)        # phylogenetic tree plotting
library(phytools)      # phylogenetic comparative methods
library(dplyr)         # data manipulation
library(tidyr)         # data reshaping
library(viridis)       # color palettes


##### Loading data #####


# Downloaded COI sequences from NCBI (dec 1, 2025)
# Searched for Oncorhynchus COI sequences using this line of code

search_results <- entrez_search(
  db = "nucleotide",
  term = "Oncorhynchus[Organism] AND (COI[Gene] OR cytochrome oxidase I[Gene]) AND 400:800[SLEN]",
  retmax = 500)


# Used cat to see the number of sequences there were

# cat("Found", search_results$count, "sequences\n")
# Found 1202 sequences

# Since there were a lot I fetched sequences in batches in a for loop and I only downloaded 500 sequences.


fetch_sequences <- function(ids, batch_size = 100) {
  all_seqs <- list()

  for (i in seq(1, length(ids), by = batch_size)) {
    batch_ids <- ids[i:min(i + batch_size - 1, length(ids))]
    batch_fetch <- entrez_fetch(
      db = "nucleotide",
      id = batch_ids,
      rettype = "fasta")

    all_seqs[[length(all_seqs) + 1]] <- batch_fetch
    Sys.sleep(0.5)
  } # to not overwhelm NCBI servers
  return(paste(all_seqs, collapse = ""))
}

# Downloaded sequences
onco_seqs_raw <- fetch_sequences(search_results$ids)

# Wrote them to file
writeLines(onco_seqs_raw, "Data/oncorhynchus_coi_raw.fasta")




# Now sequences can be read into R without having to download!
onco_seqs <- readDNAStringSet("Data/oncorhynchus_coi_raw.fasta")

cat("loaded in", length(onco_seqs), "sequences")
# Should be 500 sequences



# Extract species names from sequence headers & FASTA headers
extract_species <- function(seq_names) {
  species_names <- sapply(seq_names, function(x) {
    match <- regmatches(x, regexpr("Oncorhynchus [a-z]+", x, ignore.case = TRUE)) # ensured species & genus both listed
    if (length(match) > 0) {
      return(match)
    } else {
      return(NA)
    }
  })
  return(species_names)
}


species_list <- extract_species(names(onco_seqs))


# Exclude "Oncorhynchus sp" as its not a defined species
keep_idx <- !grepl("^Oncorhynchus sp$", species_list, ignore.case = TRUE)

onco_seqs <- onco_seqs[keep_idx]
species_list <- species_list[keep_idx]



names(onco_seqs) <- paste0("seq", 1:length(onco_seqs), "|", species_list) # removed seq(n) + | that was present before the name.

##### Quality control and filtering #####

# Check sequence lengths
seq_lengths <- width(onco_seqs)
summary(seq_lengths)

# Plotted sequence length distribution and saved figure as pdf for use in storyboard
pdf("Figures/01_raw_sequence_lengths.pdf", width = 8, height = 6)
hist(seq_lengths, breaks = 50, col = "steelblue",
  main = "Distribution of Sequence Lengths (Raw Data)",
  xlab = "Sequence Length (bp)", ylab = "Frequency")

abline(v = median(seq_lengths), col = "red", lwd = 2, lty = 2)

dev.off()



# Filtered sequences by length, keeping sequences within a reasonable range
length_threshold_min <- 500
length_threshold_max <- 700
onco_seqs_filt <- onco_seqs[seq_lengths >= length_threshold_min &
  seq_lengths <= length_threshold_max]

cat(length(onco_seqs_filt), "sequences kept after filtering for length")
# 491 sequences kept after filtering for length


# Remove sequences with >5% N's
n_content <- letterFrequency(onco_seqs_filt, "N") / width(onco_seqs_filt)
onco_seqs_filt <- onco_seqs_filt[n_content < 0.05]


cat(length(onco_seqs_filt), "sequences kept after filtering for N's")
# 491 sequences kept after filtering for N's


# Removed duplicate sequences and selected representative sequences, keeping one sequence per species
species_filt <- extract_species(names(onco_seqs_filt))
species_table <- table(species_filt)
print(species_table)

# Function to select representative sequences per species
select_representatives <- function(seqs, species_names, max_per_species = 3) {
  unique_species <- unique(species_names[!is.na(species_names)])
  selected_seqs <- DNAStringSet()

  for (sp in unique_species) {
    sp_seqs <- seqs[which(species_names == sp)]
    # Take up to max_per_species sequences
    n_select <- min(length(sp_seqs), max_per_species)
    selected_seqs <- c(selected_seqs, sp_seqs[1:n_select])
  }
  return(selected_seqs)
}

onco_seqs_final <- select_representatives(onco_seqs_filt, species_filt, max_per_species = 1)
cat("Final dataset:", length(onco_seqs_final), "sequences and",
  length(unique(extract_species(names(onco_seqs_final)))), "species")

# Write filtered sequences to data file so I don't have to run this all again if i need to clear my environment later
writeXStringSet(onco_seqs_final, "Data/oncorhynchus_coi_filtered.fasta")


##### Sequence Alignment #####

# Aligned sequences using DECIPHER



# Read sequences back
seqs_to_align <- readDNAStringSet("data/oncorhynchus_coi_filtered.fasta")

# Aligned
aligned_seqs <- AlignSeqs(seqs_to_align, iterations = 2, refinements = 2, processors = 2)

# Wrote alignment to data file
writeXStringSet(aligned_seqs, "Data/oncorhynchus_coi_aligned.fasta")

# alignment quality check
cat("Alignment length:", width(aligned_seqs)[1], "bp")
# Alignment length: 658 bp


##### Phylogenetic Analysis #####


# Converted to phydat format for phangorn
align_phyDat <- as.phyDat(as.matrix(aligned_seqs))



# Model selection

# Tested different nucleotide substitution models, tarting with a NJ tree for model testing
dist_matrix <- dist.ml(align_phyDat)
nj_tree <- NJ(dist_matrix)

model_test <- modelTest(align_phyDat, tree = nj_tree,
  model = c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "TPM2",
    "TPM3", "TIM1", "TIM2", "TIM3", "TVMe", "TVM", "SYM", "GTR"))


# Saved model test results so I didn't have to run that again as it takes a while.
write.csv(model_test, "Results/model_test_results.csv", row.names = FALSE)

# Selected best model by AIC (lowest score is the best)
best_model <- model_test$Model[which.min(model_test$AIC)]

cat("Best model by AIC is", best_model)
# Best model by AIC is TIM2+I




# Maximum likelihood tree construction

# Started by fitting the best model
fit <- pml(nj_tree, data = align_phyDat, model = best_model)

# Optimized parameters using phangorn optim.pml
fit_optimized <- optim.pml(fit, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = FALSE, rearrangement = "NNI")

# where optNni = TRUE, optimizes topology
# optBf = TRUE, optimizes base frequencies
# and optQ = TRUE,  optimizes rate matrix
# (this is just for my own reference later in case I forget) also, I'm not sure if this step was necessary as the changes weren't very drastic but oh well.


# Save the ML tree
write.tree(fit_optimized$tree, "results/oncorhynchus_ML_tree.newick")

# Bootstrap analysis (this takes a while to run would be faster if I used 100 rather than 1000)
bootstrap_trees <- bootstrap.pml(fit_optimized, bs = 1000, optNni = TRUE, multicore = TRUE, mc.cores = 2)

# Calculated bootstrap support vals
ml_tree_boot <- plotBS(fit_optimized$tree, bootstrap_trees, type = "phylogram", bs.col = "red", bs.adj = c(1, -0.2))

# Save bootstrapped tree
write.tree(ml_tree_boot, "results/oncorhynchus_ML_tree_bootstrap.newick")

# Plot tree with bootstrap support
pdf("figures/02_phylogenetic_tree_bootstrap.pdf", width = 10, height = 12)
plot(ml_tree_boot, cex = 0.7, main = "Oncorhynchus ML Phylogeny (Bootstrap Support)")
nodelabels(ml_tree_boot$node.label,
  adj = c(1, -0.2),
  frame = "none",
  cex = 0.6,
  col = "red")
add.scale.bar()
dev.off()





##### Geographic Data Acquisition #####


# Downloaded occurrence data from GBIF

# Extracted species list from the tree
tree_species <- unique(extract_species(ml_tree_boot$tip.label))
tree_species <- tree_species[!is.na(tree_species)]

cat("Species in phylogeny:", length(tree_species))

# Set function to download GBIF data for each species
download_gbif_occurrences <- function(species_name) {
  occ_data <- occ_search(scientificName = species_name, limit = 1000,
    hasCoordinate = TRUE, hasGeospatialIssue = FALSE)

  if (is.null(occ_data$data)) {
    return(data.frame())
  }

  df <- occ_data$data %>% # Extract columns
    select(species, decimalLongitude, decimalLatitude,
      coordinateUncertaintyInMeters, year) %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
  return(df)
}


# Download data for all 8 species
all_occurrences <- lapply(tree_species, function(sp) {
  cat("downloading data for", sp, "...") # (not sure where I saw this but i thought it was a nice addition when you have to sit and wait for data to load in)
  Sys.sleep(1) # to not overwhelm GBIF servers
  download_gbif_occurrences(sp)
})

# Combined a into single dataframe
occurrences_df <- bind_rows(all_occurrences)


# Saved data
write.csv(occurrences_df, "Data/oncorhynchus_occurrences_raw.csv", row.names = FALSE)


cat(nrow(occurrences_df), "occurrence records downloaded")
# 7730 occurrence records downloaded




# Cleaned and filtered occurrence data

# Removed duplicate coordinates per species, filtered out uncertain coordinates
occurrences_clean <- occurrences_df %>%
  distinct(species, decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 10000) %>%
  # I focused on the North Pacific region
  filter(decimalLongitude >= -180, decimalLongitude <= -100,
    decimalLatitude >= 30, decimalLatitude <= 70)

# Converted that to a spatial object
occurrences_sf <- st_as_sf(occurrences_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326)

# And saved the cleaned data
write.csv(occurrences_clean, "Data/oncorhynchus_occurrences_clean.csv", row.names = FALSE)
st_write(occurrences_sf, "Data/oncorhynchus_occurrences.shp", append = FALSE)

# Printed a summary by species to view my data
occ_summary <- occurrences_clean %>%
  group_by(species) %>%
  summarise(n_records = n(), lon_range = paste(round(min(decimalLongitude), 1),
    round(max(decimalLongitude), 1), sep = " to "),
  lat_range = paste(round(min(decimalLatitude), 1),
    round(max(decimalLatitude), 1), sep = " to "))

print(occ_summary)




##### Sister Species Identification #####

# Identified sister species pairs from phylogeny

# Function to find sister pairs
find_sister_pairs <- function(tree) {
  sister_pairs <- list()
  pair_count <- 0
  # internal nodes
  internal_nodes <- unique(tree$edge[, 1])

  for (node in internal_nodes) {
    # descendants of node
    descendants <- tree$edge[tree$edge[, 1] == node, 2]
    # Check if the node has exactly 2 tip descendants (sister pair)
    if (length(descendants) == 2 && all(descendants <= length(tree$tip.label))) {
      pair_count <- pair_count + 1
      sister_pairs[[pair_count]] <- list(node = node,
        species1 = tree$tip.label[descendants[1]],
        species2 = tree$tip.label[descendants[2]],
        species1_clean = extract_species(tree$tip.label[descendants[1]]),
        species2_clean = extract_species(tree$tip.label[descendants[2]]))
    }
  }
  return(sister_pairs)
}


# Applied the function
sister_pairs <- find_sister_pairs(ml_tree_boot)

cat("Identified", length(sister_pairs), "sister pairs")
# Identified 3 sister pairs


# Created dataframe for sister pairs
sister_df <- data.frame(pair_id = 1:length(sister_pairs),
  species1 = sapply(sister_pairs, function(x) x$species1_clean),
  species2 = sapply(sister_pairs, function(x) x$species2_clean))

print(sister_df)
write.csv(sister_df, "Results/sister_pairs.csv", row.names = FALSE)




##### Geographic Range Overlap Analysis #####


# Calculated range overlap between sister species

# Function to calculate range overlap using convex hulls
# convex hull draws the smallest polygon that contains all occurrence points for a species.

calculate_range_overlap <- function(sp1_coords, sp2_coords) {
  # Need at least 3 points
  if (nrow(sp1_coords) < 3 || nrow(sp2_coords) < 3) {
    return(list(overlap_area = NA, overlap_prop = NA,
      sp1_area = NA, sp2_area = NA))
  }

  # Created convex hulls
  sp1_hull <- st_convex_hull(st_union(sp1_coords))
  sp2_hull <- st_convex_hull(st_union(sp2_coords))

  # Calculated areas (km^2, using equal area projection)
  sp1_hull_proj <- st_transform(sp1_hull, crs = "+proj=laea +lat_0=45 +lon_0=-100")
  sp2_hull_proj <- st_transform(sp2_hull, crs = "+proj=laea +lat_0=45 +lon_0=-100")

  sp1_area <- as.numeric(st_area(sp1_hull_proj)) / 1e6  # converts to km^2
  sp2_area <- as.numeric(st_area(sp2_hull_proj)) / 1e6

  # Calculated intersection
  intersection <- st_intersection(sp1_hull_proj, sp2_hull_proj)

  if (length(intersection) == 0 || st_is_empty(intersection)) {
    overlap_area <- 0
  } else {
    overlap_area <- as.numeric(st_area(intersection)) / 1e6
  }

  # Calculated proportion of overlap
  overlap_prop <- overlap_area / min(sp1_area, sp2_area)

  return(list(
    overlap_area = overlap_area,
    overlap_prop = overlap_prop,
    sp1_area = sp1_area,
    sp2_area = sp2_area
  ))
}

# Calculated overlap for all sister pairs
sister_overlap_results <- list()

for (i in 1:nrow(sister_df)) {
  sp1 <- sister_df$species1[i]
  sp2 <- sister_df$species2[i]


  # Occurrence data
  sp1_data <- occurrences_sf %>% filter(species == sp1)
  sp2_data <- occurrences_sf %>% filter(species == sp2)

  if (nrow(sp1_data) > 0 && nrow(sp2_data) > 0) {
    overlap <- calculate_range_overlap(sp1_data, sp2_data)

    sister_overlap_results[[i]] <- data.frame(pair_id = i, species1 = sp1,
      species2 = sp2, n_occ_sp1 = nrow(sp1_data),
      n_occ_sp2 = nrow(sp2_data),
      range_area_sp1_km2 = overlap$sp1_area,
      range_area_sp2_km2 = overlap$sp2_area,
      overlap_area_km2 = overlap$overlap_area,
      overlap_proportion = overlap$overlap_prop)
  }
}


# Combined the results
sister_overlap_df <- bind_rows(sister_overlap_results)
write.csv(sister_overlap_df, "Results/sister_pair_overlap_analysis.csv", row.names = FALSE)

# Summary statistics
print(summary(sister_overlap_df$overlap_proportion))



# Compared to non-sister pairs

# Calculated overlap for all pairwise combinations

all_species <- unique(occurrences_clean$species)
non_sister_comparisons <- list()
count <- 0

for (i in 1:(length(all_species) - 1)) {
  for (j in (i + 1):length(all_species)) {
    sp1 <- all_species[i]
    sp2 <- all_species[j]

    # checked to see if its is not a sister pair
    is_sister <- any((sister_df$species1 == sp1 & sister_df$species2 == sp2) |
      (sister_df$species1 == sp2 & sister_df$species2 == sp1))

    if (!is_sister) {
      count <- count + 1
      sp1_data <- occurrences_sf %>% filter(species == sp1)
      sp2_data <- occurrences_sf %>% filter(species == sp2)

      if (nrow(sp1_data) >= 3 && nrow(sp2_data) >= 3) {
        overlap <- calculate_range_overlap(sp1_data, sp2_data)

        non_sister_comparisons[[count]] <- data.frame(
          species1 = sp1,
          species2 = sp2,
          overlap_proportion = overlap$overlap_prop
        )
      }
    }
  }
}

non_sister_df <- bind_rows(non_sister_comparisons)



# Statistical test: Are sister species less overlapping than expected? (wilcox test)

if (nrow(sister_overlap_df) > 0 && nrow(non_sister_df) > 0) {
  wilcox_test <- wilcox.test(sister_overlap_df$overlap_proportion,
    non_sister_df$overlap_proportion,
    alternative = "less")


  cat("Wilcoxon test: Sister vs Non-sister overlap")
  print(wilcox_test)

  # Saved the test results
  test_results <- data.frame(
    comparison = "Sister vs Non-sister",
    test = "Wilcoxon rank-sum",
    sister_median = median(sister_overlap_df$overlap_proportion, na.rm = TRUE),
    non_sister_median = median(non_sister_df$overlap_proportion, na.rm = TRUE),
    p_value = wilcox_test$p.value,
    interpretation = ifelse(wilcox_test$p.value < 0.05,
      "Sister species have significantly lower overlap (allopatry supported)",
      "No significant difference in overlap")
  )

  write.csv(test_results, "results/statistical_test_results.csv", row.names = FALSE)
}





##### Figure Creation/Visualizations #####

# Map of all species distributions

# Started with a world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Focused it on the North Pacific (relevant for what I want to analyze)
pacific_extent <- c(xmin = -180, xmax = -100, ymin = 30, ymax = 70)

pdf("Figures/03_species_distributions_map.pdf", width = 12, height = 10)
ggplot() +
  geom_sf(data = world, fill = "gray90", color = "gray50") +
  geom_sf(data = occurrences_sf, aes(color = species), size = 1.5, alpha = 0.6) +
  coord_sf(xlim = c(pacific_extent["xmin"], pacific_extent["xmax"]),
    ylim = c(pacific_extent["ymin"], pacific_extent["ymax"])) +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  labs(title = "Geographic Distribution of Oncorhynchus Species",
    x = "Longitude", y = "Latitude",
    color = "Species") +
  theme(legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"))
dev.off()




# Sister pair overlap comparison plots

# Prepared data for plotting
plot_data <- rbind(
  data.frame(group = "Sister pairs",
    overlap = sister_overlap_df$overlap_proportion),
  data.frame(group = "Non-sister pairs",
    overlap = non_sister_df$overlap_proportion)
)

pdf("Figures/04_sister_vs_nonsister_overlap.pdf", width = 8, height = 6)
ggplot(plot_data, aes(x = group, y = overlap, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Sister pairs" = "#E64B35",
    "Non-sister pairs" = "#4DBBD5")) +
  labs(title = "Geographic Range Overlap: Sister vs Non-sister Species",
    subtitle = "Lower overlap in sister pairs supports allopatric speciation",
    x = "", y = "Proportion of Range Overlap",
    caption = paste0("Wilcoxon test p-value: ",
      round(wilcox_test$p.value, 4))) +
  theme_minimal() +
  theme(legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11))
dev.off()





# Phylogeny with geographic information

# Created a color scheme for species based on longitude centroid
species_centroids <- occurrences_clean %>%
  group_by(species) %>%
  summarise(
    centroid_lon = mean(decimalLongitude),
    centroid_lat = mean(decimalLatitude)
  )

# Matched to points at tree tips
tip_data <- ml_tree_boot$tip.label %>%
  data.frame(label = ., stringsAsFactors = FALSE) %>%
  mutate(species = extract_species(label)) %>%
  left_join(species_centroids, by = "species")

# %<+% attaches a data frame to the tree
tree_plot <- ggtree(ml_tree_boot, layout = "rectangular") %<+% tip_data

pdf("figures/05_phylogeny_with_geography.pdf", width = 10, height = 12)
tree_plot +
  geom_tiplab(aes(label = species), size = 3, offset = 0.001) +
  geom_tippoint(aes(color = centroid_lon), size = 3) +
  geom_text2(aes(subset = !isTip, label = label),
    data = fortify(ml_tree_boot),
    size = 2.5, hjust = 1.2, vjust = -0.3) +
  scale_color_viridis_c(option = "plasma", name = "Longitude\nCentroid") +
  theme_tree2() +
  labs(
    title = "Oncorhynchus Phylogeny with Geographic Distribution",
    subtitle = "Tip colors show longitudinal centroid of species ranges"
  )

dev.off()





# Individual sister pair maps

# Created maps for each sister pair
for (i in 1:nrow(sister_df)) {
  sp1 <- sister_df$species1[i]
  sp2 <- sister_df$species2[i]


  sp1_data <- occurrences_sf %>% filter(species == sp1)
  sp2_data <- occurrences_sf %>% filter(species == sp2)

  if (nrow(sp1_data) > 0 && nrow(sp2_data) > 0) {
    pair_data <- rbind(
      data.frame(species = sp1, geometry = st_geometry(sp1_data)),
      data.frame(species = sp2, geometry = st_geometry(sp2_data))
    )
    pair_sf <- st_as_sf(pair_data)

    pdf(paste0("Figures/06_sister_pair_", i, "_",
      gsub(" ", "_", sp1), "_vs_", gsub(" ", "_", sp2), ".pdf"),
    width = 10, height = 8)

    p <- ggplot() +
      geom_sf(data = world, fill = "gray90", color = "gray50") +
      geom_sf(data = pair_sf, aes(color = species), size = 2, alpha = 0.7) +
      coord_sf(xlim = c(pacific_extent["xmin"], pacific_extent["xmax"]),
        ylim = c(pacific_extent["ymin"], pacific_extent["ymax"])) +
      scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
      theme_minimal() +
      labs(title = paste("Sister Pair", i, ":", sp1, "vs", sp2),
        subtitle = paste("Overlap proportion:",
          round(sister_overlap_df$overlap_proportion[i], 3)),
        x = "Longitude", y = "Latitude") +
      theme(plot.title = element_text(size = 12, face = "bold"))

    print(p)
    dev.off()
  }
}




# Reproduceability checkpoints to ensure we have the right data

cat("Number of sequences downloaded:", search_results$count)
cat("Number of sequences after filtering:", length(onco_seqs_final))
cat("Number of species in phylogeny:", length(unique(tree_species)))
cat("Number of occurrence records:", nrow(occurrences_clean))
cat("Number of sister pairs identified:", nrow(sister_df))
