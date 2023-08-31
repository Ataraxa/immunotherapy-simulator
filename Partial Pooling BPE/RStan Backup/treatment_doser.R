treatment_doser <- function(t, t_array, delay, last) {
  for (t_in in 1:length(t_array)) {
    if (((t_in + delay) < t) && (t < (t_in + delay + last))) {
      return(TRUE)
    }
  }
  return(FALSE)
}


