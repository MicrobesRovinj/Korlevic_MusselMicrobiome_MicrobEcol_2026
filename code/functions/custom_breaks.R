#################################################################################################################
# custom_breaks.R
#
# Function to set custom axis breaks.
#
#################################################################################################################

custom_breaks <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (m > 10     && m <= 20)    seq(0, 20, by = 5) else
  if (m > 20     && m <= 40)    seq(0, 40, by = 10) else
  if (m > 40     && m <= 60)    seq(0, 60, by = 10) else
  if (m > 60     && m <= 80)    seq(0, 80, by = 20) else
  if (m > 80     && m <= 100)   seq(0, 100, by = 25) else
  if (m > 100    && m <= 200)   seq(0, 200, by = 50) else
  if (m > 200    && m <= 600)   seq(0, 600, by = 100) else
  if (m > 600    && m <= 800)   seq(0, 800, by = 200) else
  if (m > 800    && m <= 1000)  seq(0, 1000, by = 200) else
  if (m > 1000   && m <= 1200)  seq(0, 1200, by = 300) else
  if (m > 1200   && m <= 2000)  seq(0, 2000, by = 500) else
  if (m > 2000   && m <= 2500)  seq(0, 2500, by = 500) else
  if (m > 2500   && m <= 5000)  seq(0, 5000, by = 1000) else
  if (m > 5000   && m <= 8000)  seq(0, 8000, by = 2000) else
  if (m > 8000   && m <= 10000) seq(0, 10000, by = 2000) else
  if (m > 10000  && m <= 15000) seq(0, 15000, by = 3000) else
  if (m > 15000  && m <= 20000) seq(0, 20000, by = 4000) else
  if (m > 20000  && m <= 25000) seq(0, 25000, by = 5000) else
  if (m > 25000  && m <= 35000) seq(0, 35000, by = 7000) else
  if (m > 35000  && m <= 40000) seq(0, 40000, by = 8000) else
  if (m > 40000  && m <= 75000) seq(0, 75000, by = 15000) else
  if (m > 75000)                seq(0, 100000, by = 20000) else
    numeric(length = 0)
}

