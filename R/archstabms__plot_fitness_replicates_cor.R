
#' archstabms__plot_fitness_replicates_cor
#'
#' Plot fitness replicate correlations.
#'
#' @param input_dt list of folder paths to fitness data (required)
#' @param output_file output file (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__plot_fitness_replicates_cor <- function(
  input_dt, 
  output_file,
  colour_scheme
  ){
  
  fitness_dt <- input_dt

  #Plot replicate fitness correlation
  num_rep = 3
  archstabms__ggpairs_binhex(
    input_dt = fitness_dt[STOP==F & STOP_readthrough==F,.SD,.SDcols = paste0('fitness', c(1:num_rep), "_uncorr")], 
    output_file = output_file,
    xlab = "Fitness",
    ylab = "Fitness",
    width = 3,
    height = 3,
    label_size = 2,
    bins = 20,
    plot_colours = c("lightgrey", "black"))
  
}
