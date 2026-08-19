#################################################################################################################
# custom_avg_bray.R
#
# Function to customise the calculation of average Bray–Curtis dissimilarity coefficients.
#
#################################################################################################################

custom_avg_bray <- function(shared_metadata = NULL, filter_environment = NULL,
                            sample = NULL, subsamples = NULL) {
  
  # Ensure that the required arguments are provided
  if(is.null(shared_metadata) || is.null(sample) || is.null(subsamples)) {
    stop("'shared_metadata', 'sample', and 'subsamples' must be provided.")
  }
  
  # Filter samples from the environment specified by the argument
  # "filter_environment"
  if(!is.null(filter_environment)) {
    shared_metadata <- shared_metadata %>%
      # Filter samples from the specified environment
      filter(environment == {{ filter_environment }}) %>%
      # Remove columns that contain only NAs
      select(where(fn = ~ all(!is.na(.x)))) %>%
      # Remove OTU abundance columns that contain only 0
      select(!starts_with(match = "Otu"), starts_with(match = "Otu") &
               where(fn = ~ any(.x != 0)))
  }
  
  # Calculate the mean Bray–Curtis dissimilarity from subsamples of the
  # OTU/sample data (rarefaction)
  avg_bray <- shared_metadata %>%
    # Select the column that contains sample IDs and the columns that contain
    # OTU abundances
    select(id, starts_with(match = "Otu")) %>%
    # Add sample IDs to row names (input for the vegan package)
    column_to_rownames(var = "id") %>%
    # Calculate the mean Bray–Curtis dissimilarity from the number of
    # subsamples specified by the argument "subsamples"; random subsampling
    # is performed at the subsampling depth defined by the argument "sample"
    avgdist(sample = sample, istfun = vegdist, meanfun = mean,
            iterations = subsamples, dmethod = "bray")
  
  # Return the output
  list(shared_metadata = shared_metadata,
       avg_bray = avg_bray)
  
}

