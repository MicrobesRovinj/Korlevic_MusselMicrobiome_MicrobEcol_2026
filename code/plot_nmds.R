#################################################################################################################
# plot_nmds.R
#
# Script to generate NMDS plots.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared
#               data/raw/metadata.tsv
#               code/functions/format_labels.R
#               code/functions/custom_avg_bray.R
#               data/raw/colour_environment.R
#               data/raw/shape_site.R
#               data/raw/theme.R
#               data/raw/colour_month.R
#               data/raw/shape_location.R
# Produces: results/numerical/permanova.Rdata
#           results/figures/nmds.jpg
#
#################################################################################################################

# Load OTU/sample data
shared <- read_tsv(file = "data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared")

# Load metadata
metadata <- read_tsv("data/raw/metadata.tsv")

# Customise metadata using the custom function
metadata <- format_labels(x = metadata)

# Further customise metadata
metadata <- metadata %>%
  # Combine sampling month and year into one column
  mutate(month_year = paste0(month, " ", year), .after = year)

# Join metadata with OTU/sample data
shared_metadata <- shared %>%
  # Rename the column "Group" to "id"
  rename(id = Group) %>%
  # Convert the sample ID column to character type
  mutate(id = as.character(x = id)) %>%
  # Join metadata with OTU/sample data
  right_join(x = metadata, y = ., by = c("id" = "id"))

# Calculate the minimum number of sequences per sample
min_seq <- shared_metadata %>%
  # Select columns that contain OTU abundances
  select(starts_with(match = "Otu")) %>%
  # Sum the sequences in each sample
  rowSums() %>%
  # Find the minimum value
  min()

#################################################################################################################
# PERMANOVA
#################################################################################################################

# Environments

# Calculate the mean Bray–Curtis dissimilarity from 100 subsamples of the
# OTU/sample data (rarefaction) using the custom function
environments <- custom_avg_bray(shared_metadata = shared_metadata,
                                filter_environment = NULL,
                                sample = min_seq, subsamples = 100)
# Set seed
set.seed(19800101)
# Calculate PERMANOVA
(permanova_environments <- adonis2(
  formula = environments$avg_bray ~ environment,
  data = environments$shared_metadata, permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_environments <- betadisper(
  d = environments$avg_bray,
  group = environments$shared_metadata$environment,
  add = "lingoes"))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_environments <- permutest(betadisper_environments,
                                     permutations = 999))
# Set seed
set.seed(19800101)
# Calculate pairwise PERMANOVA
(pairwise_permanova_environments <- pairwise.adonis(
  x = environments$avg_bray,
  factors = environments$shared_metadata$environment,
  p.adjust.m = "bonferroni", perm = 999))

# Seawater

# Calculate the mean Bray–Curtis dissimilarity from 100 subsamples of the
# OTU/sample data (rarefaction) using the custom function
seawater <- custom_avg_bray(shared_metadata = shared_metadata,
                            filter_environment = "Seawater",
                            sample = min_seq, subsamples = 100)
# Set seed
set.seed(19800101)
# Calculate PERMANOVA (sites)
(permanova_seawater_site <- adonis2(formula = seawater$avg_bray ~ site,
                                    data = seawater$shared_metadata,
                                    permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_seawater_site <- betadisper(d = seawater$avg_bray,
                                        group = seawater$shared_metadata$site))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_seawater_site <- permutest(betadisper_seawater_site,
                                     permutations = 999))

# Set seed
set.seed(19800101)
# Calculate PERMANOVA (locations)
(permanova_seawater_location <- adonis2(formula = seawater$avg_bray ~ location,
                                       data = seawater$shared_metadata,
                                       permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_seawater_location <- betadisper(
  d = seawater$avg_bray, group = seawater$shared_metadata$location))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_seawater_location <- permutest(betadisper_seawater_location,
                                          permutations = 999))

# Sediment

# Calculate the mean Bray–Curtis dissimilarity from 100 subsamples of the
# OTU/sample data (rarefaction) using the custom function
sediment <- custom_avg_bray(shared_metadata = shared_metadata,
                            filter_environment = "Sediment",
                            sample = min_seq, subsamples = 100)
# Set seed
set.seed(19800101)
# Calculate PERMANOVA (sites)
(permanova_sediment_site <- adonis2(formula = sediment$avg_bray ~ site,
                                    data = sediment$shared_metadata,
                                    permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_sediment_site <- betadisper(
  d = sediment$avg_bray, group = sediment$shared_metadata$site))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_sediment_site <- permutest(betadisper_sediment_site,
                                      permutations = 999))
# Set seed
set.seed(19800101)
# Calculate PERMANOVA (locations)
(permanova_sediment_location <- adonis2(formula = sediment$avg_bray ~ location,
                                        data = sediment$shared_metadata,
                                        permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_sediment_location <- betadisper(
  d = sediment$avg_bray, group = sediment$shared_metadata$location))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_sediment_location <- permutest(betadisper_sediment_location,
                                          permutations = 999))

# Gills

# Calculate the mean Bray–Curtis dissimilarity from 100 subsamples of the
# OTU/sample data (rarefaction) using the custom function
gills <- custom_avg_bray(shared_metadata = shared_metadata,
                         filter_environment = "Gills",
                         sample = min_seq, subsamples = 100)
# Set seed
set.seed(19800101)
# Calculate PERMANOVA (locations)
(permanova_gills_location <- adonis2(formula = gills$avg_bray ~ location,
                                     data = gills$shared_metadata,
                                     permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_gills_location <- betadisper(
  d = gills$avg_bray, group = gills$shared_metadata$location, add = "lingoes"))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_gills_location <- permutest(betadisper_gills_location,
                                       permutations = 999))
# Set seed
set.seed(19800101)
# Calculate PERMANOVA (sampling months)
(permanova_gills_month_year <- adonis2(formula = gills$avg_bray ~ month_year,
                                       data = gills$shared_metadata,
                                       permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_gills_month_year <- betadisper(
  d = gills$avg_bray,
  group = gills$shared_metadata$month_year,
  add = "lingoes"))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_gills_month_year <- permutest(betadisper_gills_month_year,
                                         permutations = 999))
# Set seed
set.seed(19800101)
# Calculate pairwise PERMANOVA (sampling months)
(pairwise_permanova_gills_month_year <- pairwise.adonis(
  x = gills$avg_bray, factors = gills$shared_metadata$month_year,
  p.adjust.m = "bonferroni", perm = 999))

# Gills (excluding samples from January 2021)

# Calculate the mean Bray–Curtis dissimilarity from 100 subsamples of the
# OTU/sample data (rarefaction) using the custom function
gills_no_jan_2021 <- shared_metadata %>%
  # Filter out gill samples from January 2021
  filter(!(environment == "Gills" & month_year == "January 2021")) %>%
  # Calculate the mean Bray–Curtis dissimilarity from 100 subsamples
  custom_avg_bray(filter_environment = "Gills",
                  sample = min_seq, subsamples = 100)
# Set seed
set.seed(19800101)
# Calculate PERMANOVA (sampling months)
(permanova_gills_no_jan_2021_month_year <- adonis2(
  formula = gills_no_jan_2021$avg_bray ~ month_year,
  data = gills_no_jan_2021$shared_metadata,
  permutations = 999))
# Calculate multivariate homogeneity of group dispersions
(betadisper_gills_no_jan_2021_month_year <- betadisper(
  d = gills_no_jan_2021$avg_bray,
  group = gills_no_jan_2021$shared_metadata$month_year,
  add = "lingoes"))
# Set seed
set.seed(19800101)
# Perform a permutation test for multivariate homogeneity of group
# dispersions
(permutest_gills_no_jan_2021_month_year <- permutest(
  betadisper_gills_no_jan_2021_month_year,
  permutations = 999))
# Set seed
set.seed(19800101)
# Calculate pairwise PERMANOVA (sampling months)
(pairwise_permanova_gills_no_jan_2021_month_year <- pairwise.adonis(
  x = gills_no_jan_2021$avg_bray,
  factors = gills_no_jan_2021$shared_metadata$month_year,
  p.adjust.m = "bonferroni", perm = 999))

# Combine the PERMANOVA results
permanova <- list(
  environments = list(
    permanova = permanova_environments,
    betadisper = betadisper_environments,
    betadisper_permutest = permutest_environments,
    pairwise_permanova = pairwise_permanova_environments),
  seawater = list(
    site = list(permanova = permanova_seawater_site,
                betadisper = betadisper_seawater_site,
                betadisper_permutest = permutest_seawater_site),
    location = list(permanova = permanova_seawater_location,
                    betadisper = betadisper_seawater_location,
                    betadisper_permutest = permutest_seawater_location)),
  sediment = list(
    site = list(permanova = permanova_sediment_site,
                betadisper = betadisper_sediment_site,
                betadisper_permutest = permutest_sediment_site),
    location = list(permanova = permanova_sediment_location,
                    betadisper = betadisper_sediment_location,
                    betadisper_permutest = permutest_sediment_location)),
  gills = list(
    location = list(permanova = permanova_gills_location,
                    betadisper = betadisper_gills_location,
                    betadisper_permutest = permutest_gills_location),
    month_year = list(
      permanova = permanova_gills_month_year,
      betadisper = betadisper_gills_month_year,
      betadisper_permutest = permutest_gills_month_year,
      pairwise_permanova = pairwise_permanova_gills_month_year),
    month_year_no_jan_2021 = list(
      permanova = permanova_gills_no_jan_2021_month_year,
      betadisper = betadisper_gills_no_jan_2021_month_year,
      betadisper_permutest = permutest_gills_no_jan_2021_month_year,
      pairwise_permanova = pairwise_permanova_gills_no_jan_2021_month_year)))

# Save the combined PERMANOVA results
save(permanova, file = "results/numerical/permanova.Rdata")

#################################################################################################################
# NMDS (environments and sites)
#################################################################################################################

# Set seed
set.seed(19800101)
# Calculate NMDS
nmds <- metaMDS(comm = environments$avg_bray, k = 2)

# Extract point coordinates
coordinates <- nmds %>%
  # Get site scores
  scores(display = "sites") %>%
  # Convert the matrix to a tibble
  as_tibble(rownames = "id")

# Join metadata with point coordinates
metadata_coordinates <- inner_join(x = metadata, y = coordinates,
                                   by = c("id" = "id"))

# Load plot customisation data
source(file = "data/raw/colour_environment.R")
source(file = "data/raw/shape_site.R")

# Generate plot
p_environments_sites <- metadata_coordinates %>%
  # Initialise a ggplot object and define aesthetic mappings
  ggplot(mapping = aes(x = NMDS1, y = NMDS2, fill = environment)) +
  # Add points, define aesthetic mappings, and specify the appearance
  # of the points
  geom_point(mapping = aes(shape = site), size = 5, stroke = 0.5) +
  # Add ellipses, define aesthetic mappings, and specify the appearance
  # of the ellipses
  stat_ellipse(mapping = aes(colour = environment), geom = "polygon",
               level = 0.95, alpha = 0.25, show.legend = FALSE) +
  # Customise the continuous scale of the x-axis
  scale_x_continuous(breaks = c(-0.4, 0, 0.4)) +
  # Customise the continuous scale of the y-axis
  scale_y_continuous(breaks = c(-0.4, 0, 0.4)) +
  # Specify fill colours for the points
  scale_fill_manual(name = "Environment",
                    guide = guide_legend(
                      override.aes = list(shape = 22, size = 7)),
                    values = colour_environment,
                    breaks = names(x = colour_environment)) +
  # Specify the shape of the points
  scale_shape_manual(name = "Site",
                     values = shape_site,
                     breaks = names(x = shape_site)) +
  # Specify line colours
  scale_colour_manual(name = "Environment",
                      values = colour_environment,
                      breaks = names(x = colour_environment)) +
  # Define axis titles and plot title
  labs(title = paste0("Environments and Sites (Stress = ",
                      format(x = round(x = nmds$stress, digits = 2),
                             nsmall = 2), ")"),
       x = "NMDS 1",
       y = "NMDS 2") +
  # Add the PERMANOVA results
  annotate(geom = "text", x = I(x = 0.42), y = I(x = 0.20),
           label = "underline(bold('PERMANOVA'))",
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  annotate(geom = "text", x = I(x = 0.42), y = I(x = 0.15),
           label = paste0(
             "bold('Environments:')~bolditalic('R'^'2')~bold('=')~",
             "bold('", round(x = permanova_environments$R2[1], digits = 2),
             ",')~bolditalic('p')~bold('",
             format_p_values(x = permanova_environments$`Pr(>F)`[1],
                             space = TRUE), "*')"),
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  annotate(geom = "text", x = I(x = 0.42), y = I(x = 0.10),
           label = paste0(
             "bold('Sites (Seawater):')~bolditalic('R'^'2')~bold('=')~",
             "bold('", round(x = permanova_seawater_site$R2[1], digits = 2),
             ",')~bolditalic('p')~bold('= ",
             format_p_values(x = permanova_seawater_site$`Pr(>F)`[1],
                             space = TRUE), "')"),
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  annotate(geom = "text", x = I(x = 0.42), y = I(x = 0.05),
           label = paste0(
             "bold('Sites (Sediment):')~bolditalic('R'^'2')~bold('=')~",
             "bold('", round(x = permanova_sediment_site$R2[1], digits = 2),
             ",')~bolditalic('p')~bold('= ",
             format_p_values(x = permanova_sediment_site$`Pr(>F)`[1],
                             space = TRUE), "')"),
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  # Set the plot aspect ratio and define axis limits
  coord_fixed(ratio = 1, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  # Use a general custom theme
  theme +
  # Add additional plot customisations
  theme(panel.border = element_rect(fill = NA))

#################################################################################################################
# NMDS (gill microbiome)
#################################################################################################################

# Set seed
set.seed(19800101)
# Calculate NMDS
nmds <- metaMDS(comm = gills$avg_bray, k = 2, trymax = 200)

# Extract the point coordinates
coordinates <- nmds %>%
  # Get the site scores
  scores(display = "sites") %>%
  # Convert the matrix to a tibble
  as_tibble(rownames = "id")

# Join metadata with point coordinates
metadata_coordinates <- inner_join(x = metadata, y = coordinates,
                                   by = c("id" = "id"))

# Load the plot customisation data
source(file = "data/raw/colour_month.R")
source(file = "data/raw/shape_location.R")

# Generate plot
p_gills <- metadata_coordinates %>%
  # Initialise a ggplot object and define the aesthetic mappings
  ggplot(mapping = aes(x = NMDS1, y = NMDS2,
                       fill = month_year)) +
  # Add points, define the aesthetic mappings, and specify the appearance
  # of the points
  geom_point(mapping = aes(shape = location), size = 5, stroke = 0.5) +
  # Add ellipses, define the aesthetic mappings, and specify the appearance
  # of the ellipses
  stat_ellipse(mapping = aes(colour = month_year), geom = "polygon",
               level = 0.95, alpha = 0.25, show.legend = FALSE) +
  # Customise the continuous scale of the x-axis
  scale_x_continuous(breaks = c(-0.4, 0, 0.4, 0.8)) +
  # Customise the continuous scale of the y-axis
  scale_y_continuous(breaks = c(-0.8, -0.4, 0, 0.4)) +
  # Specify fill colours for the points
  scale_fill_manual(name = "Month",
                    guide = guide_legend(
                      override.aes = list(shape = 22, size = 7), order = 1),
                    values = colour_month,
                    breaks = names(x = colour_month)) +
  # Specify the shape of the points
  scale_shape_manual(name = "Location",
                     values = shape_location,
                     breaks = names(x = shape_location)) +
  # Specify line colours
  scale_colour_manual(name = "Month",
                      values = colour_month,
                      breaks = names(x = colour_month)) +
  # Define axis titles and plot title
  labs(title = paste0("Gill Microbiome (Stress = ",
                      format(x = round(x = nmds$stress, digits = 2),
                             nsmall = 2), ")"),
       x = "NMDS 1",
       y = "NMDS 2") +
  # Add the PERMANOVA results
  annotate(geom = "text", x = I(x = 0.52), y = I(x = 0.15),
           label = "underline(bold('PERMANOVA'))",
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  annotate(geom = "text", x = I(x = 0.52), y = I(x = 0.10),
           label = paste0(
             "bold('Locations:')~bolditalic('R'^'2')~bold('=')~",
             "bold('", round(x = permanova_gills_location$R2[1],
                             digits = 2),
             ",')~bolditalic('p')~bold('",
             format_p_values(
               x = permanova_gills_location$`Pr(>F)`[1],
               space = TRUE), "')"),
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  annotate(geom = "text", x = I(x = 0.52), y = I(x = 0.05),
           label = paste0(
             "bold('Months:')~bolditalic('R'^'2')~bold('=')~",
             "bold('", round(x = permanova_gills_no_jan_2021_month_year$R2[1],
                             digits = 2),
             ",')~bolditalic('p')~bold('",
             format_p_values(
               x = permanova_gills_no_jan_2021_month_year$`Pr(>F)`[1],
               space = TRUE), "')"),
           family = "Times",  size = 14 / .pt, hjust = 0, parse = TRUE) +
  # Set the plot aspect ratio and define axis limits
  coord_fixed(ratio = 1, xlim = c(-0.5, 0.8), ylim = c(-0.8, 0.5)) +
  # Use a general custom theme
  theme +
  # Add additional plot customisations
  theme(panel.border = element_rect(fill = NA))

# Combine plots
p <- plot_grid(p_environments_sites, NULL, p_gills, nrow = 1,
               align = "hv", rel_widths = c(1, -0.04, 1))

# Save
ggsave(filename = "results/figures/nmds.jpg", plot = p,
       width = 297 * 1.35, height = 210 * 0.8, units = "mm")

