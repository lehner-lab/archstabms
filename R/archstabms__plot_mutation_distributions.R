
#' archstabms__plot_mutation_distributions
#'
#' Plot mutation distributions.
#'
#' @param input_dt data.table with fitness estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__plot_mutation_distributions <- function(
  input_dt,
  report_outpath,
  dataset_name,
  colour_scheme
  ){

  #Plot data.table
  plot_dt <- copy(input_dt)

  #Maximum theoretical mutation order
  max_order <- sum(sapply(apply(as.data.table(t(plot_dt[, strsplit(aa_seq, "")])), 2, table), length)>1)

  #Median mutation order
  median_order <- plot_dt[,median(Nham_aa)]

  #Density plots
  plot_dt[, mut_order := Nham_aa]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Nham_aa)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = max_order, linetype = 2) +
    ggplot2::geom_vline(xintercept = median_order) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Number of AA substitutions") +
    ggplot2::ylab("Count")
  ggplot2::ggsave(file.path(report_outpath, paste0("mutation_density_Nham_aa_", dataset_name, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
}
