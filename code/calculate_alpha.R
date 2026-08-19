#################################################################################################################
# calculate_alpha.R
#
# Script to calculate richness estimators and diversity indices.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared
#               data/raw/metadata.tsv
#               code/functions/format_labels.R
#               code/functions/custom_alpha_parameters.R
# Produces: results/numerical/alpha.Rdata
#
#################################################################################################################

# Load OTU/sample data
shared <- read_tsv(file = "data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared")

# Load metadata
metadata <- read_tsv("data/raw/metadata.tsv")

# Customise metadata using the custom function
metadata <- format_labels(x = metadata)

# Subsample the OTU/sample data 100 times (rarefaction) using the custom
# function and, after each subsampling, calculate alpha diversity parameters
# (richness estimators and diversity indices)
alpha_subsampling <- map(.x = c(1 : 100),
    .f = ~ custom_alpha_parameters(shared = shared, subsample = .x)) %>%
  # Combine list elements into a tibble
  list_rbind()

# Calculate the mean value for each alpha diversity parameter from the
# results of all subsamples
alpha_mean_sd <- alpha_subsampling %>%
  # Group by sample ID
  group_by(id) %>%
  # Calculate the mean and standard deviation
  summarise(S.obs_mean = mean(x = S.obs),
            S.obs_sd = sd(x = S.obs),
            S.chao1_mean = mean(x = S.chao1),
            S.chao1_sd = sd(x = S.chao1),
            S.ACE_mean = mean(x = S.ACE),
            S.ACE_sd = sd(x = S.ACE),
            eshannon_mean = mean(x = eshannon),
            eshannon_sd = sd(x = eshannon),
            invsimpson_mean = mean(x = invsimpson),
            invsimpson_sd = sd(x = invsimpson))

# Join metadata with the calculated richness estimators and diversity indices
alpha <- inner_join(x = metadata, y = alpha_mean_sd, by = c("id" = "id"))

# Pivot the data longer
alpha <- alpha %>%
  pivot_longer(cols = matches(match = "_(mean|sd)$"),
               names_to = c("parameter", "mean_sd"),
               names_pattern = "^(.+)_(mean|sd)$",
               values_to = "value")

# Pivot the data wider to place mean and standard deviation values in
# separate columns
alpha <- alpha %>%
  pivot_wider(names_from = mean_sd, values_from = value)

# Sort factor levels of the variable parameter
alpha <- alpha %>%
  mutate(parameter = factor(x = parameter, levels = c("S.obs",
                                                      "S.chao1",
                                                      "S.ACE",
                                                      "eshannon",
                                                      "invsimpson")))

# Save
save(alpha, file = "results/numerical/alpha.Rdata")

