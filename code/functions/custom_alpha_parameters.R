#################################################################################################################
# custom_alpha_parameters.R
#
# Function to customise the calculation of alpha diversity parameters (observed number of OTUs, Chao1, ACE,
# exponential of the Shannon diversity index, and Inverse Simpson diversity index).
#
#################################################################################################################

custom_alpha_parameters <- function(shared = NULL, subsample = NULL) {
  
  # Ensure that the required arguments are provided
  if(is.null(shared) || is.null(subsample)) {
    stop("'shared' and 'subsample' must be provided.")
  }
  
  # Prepare the OTU/sample data to generate a randomly rarefied version of the
  # community data
  shared <- shared %>%
    # Rename the column "Group" to "id"
    rename(id = Group) %>%
    # Keep the column containing sample IDs and remove abundance columns
    # that contain only 0
    select(id, starts_with(match = "Otu") & where(fn = ~ any(.x != 0))) %>%
    # Add sample IDs to row names (input for the vegan package)
    column_to_rownames(var = "id")
  
  # Create a randomly rarefied version of the community data
  rarefied <- shared %>%
    # Apply vegan's rrarefy() function to generate a randomly rarefied
    # version of the community data
    rrarefy(., sample = min(rowSums(x = .))) %>%
    # Convert the output to a tibble
    as_tibble(.name_repair = "check_unique", rownames = NA) %>%
    # Add sample IDs from row names as a column
    rownames_to_column(var = "id") %>%
    # Keep the column containing sample IDs and remove abundance columns
    # that contain only 0
    select(id, starts_with(match = "Otu") & where(fn = ~ any(.x != 0)))
  
  # Add sample IDs to row names (input for the vegan package)
  rarefied <- rarefied %>%
    column_to_rownames(var = "id")
  
  # Calculate the observed number of OTUs and the estimators Chao1 and ACE
  estimators <- rarefied %>%
    # Calculate the estimators
    estimateR() %>%
    # Transpose the calculated data
    t() %>%
    # Convert the matrix to a tibble
    as_tibble(.name_repair = "check_unique", rownames = NA) %>%
    # Add sample IDs from the row names as a column
    rownames_to_column("id")
  
  # Calculate the diversity indices (Shannon entropy and inverse Simpson)
  shannon <- rarefied %>%
    # Calculate the Shannon entropy
    diversity(index = "shannon") %>%
    # Convert the vector to a tibble
    enframe(name = "id", value = "shannon")
  invsimpson <- rarefied %>%
    # Calculate the inverse Simpson diversity index
    diversity(index = "invsimpson") %>%
    # Convert the vector to a tibble
    enframe(name = "id", value = "invsimpson")
  
  # Transform the Shannon entropy to the effective number of OTUs
  # (http://www.loujost.com/Statistics%20and%20Physics/Diversity%20and%20Similarity/EffectiveNumberOfSpecies.htm)
  eshannon <- mutate(shannon, shannon = exp(shannon)) %>%
    # Rename the column "shannon" to "eshannon"
    rename(eshannon = shannon)
  
  # Join the richness estimators and the diversity indices
  alpha <- inner_join(x = estimators, y = eshannon, by = c("id" = "id")) %>%
    inner_join(y = invsimpson, by = c("id" = "id"))
  
  # Add subsample number
  alpha <- alpha %>%
    mutate(subsample = subsample, .before = "id")
  
  # Return the output
  alpha
  
}

