#################################################################################################################
# custom_num_to_text.R
#
# Convert integers to text if they are less than 10.
#
#################################################################################################################

custom_num_to_text <- function(x = NULL) {
  
  # Ensure that the required argument is provided
  if (is.null(x)) {
    stop("Input x cannot be NULL.")
  }
  
  # Convert the integer to text
  result <- case_when(x == 1 ~ "one",
                      x == 2 ~ "two",
                      x == 3 ~ "three",
                      x == 4 ~ "four",
                      x == 5 ~ "five",
                      x == 6 ~ "six",
                      x == 7 ~ "seven",
                      x == 8 ~ "eight",
                      x == 9 ~ "nine",
                      TRUE ~ as.character(x = x))
  
  # Return the output
  return(result)
  
}

