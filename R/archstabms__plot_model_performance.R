
#' archstabms__plot_model_performance
#'
#' Plot model performance.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param trait trait name and suffix for output files (required)
#' @param highlight_colour colour for highlights (default:red)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__plot_model_performance <- function(
  input_dt,
  report_outpath,
  trait,
  highlight_colour = "red"
  ){

  #Plot observed versus predicted fitness
  plot_dt <- input_dt[mut_order>0 & mut_order<15,.(mut_order, observed_fitness = fitness, predicted_fitness = fold_1)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),mut_order], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_wrap(~mut_order, scales = "free", nrow = 3) + 
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_", trait, "_mutorderfacet_model1.pdf")), d, width = 15, height = 9, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order
  plot_dt <- input_dt[,.(observed_fitness = fitness, predicted_fitness = fold_1, sigma)]
  #Maximum explainable variance
  max_r2 <- 1 - plot_dt[sigma<plot_dt[,quantile(sigma, 0.75)],sum(sigma^2)/(var(observed_fitness)*.N)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2/max_r2, 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_", trait, "_model1.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(observed_fitness = fitness, predicted_fitness = fitness_pred, sigma)]
  #Maximum explainable variance
  max_r2 <- 1 - plot_dt[sigma<plot_dt[,quantile(sigma, 0.75)],sum(sigma^2)/(var(observed_fitness)*.N)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2/max_r2, 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_", trait, "_allmodels_test.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(observed_fitness = fitness, predicted_fitness = fitness_pred, sigma)]
  #Maximum explainable variance
  max_r2 <- 1 - plot_dt[sigma<plot_dt[,quantile(sigma, 0.75)],sum(sigma^2)/(var(observed_fitness)*.N)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2/max_r2, 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_", trait, "_allmodels_test_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(observed_fitness = fitness, predicted_fitness = fitness_pred, sigma)]
  #Maximum explainable variance
  max_r2 <- 1 - plot_dt[sigma<plot_dt[,quantile(sigma, 0.75)],sum(sigma^2)/(var(observed_fitness)*.N)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2/max_r2, 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_", trait, "_allmodels_test_highres.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(observed_fitness = fitness, predicted_fitness = fitness_pred, sigma)]
  #Maximum explainable variance
  max_r2 <- 1 - plot_dt[sigma<plot_dt[,quantile(sigma, 0.75)],sum(sigma^2)/(var(observed_fitness)*.N)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2/max_r2, 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_", trait, "_allmodels_test_highres_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(predicted_fitness = fitness_pred, residual_fitness = fitness-fitness_pred)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(predicted_fitness, residual_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::geom_smooth(color = highlight_colour, formula = y ~ s(x, bs = "cs"), method = "gam") +
    ggplot2::xlab("Predicted fitness") +
    ggplot2::ylab("Residual fitness") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  # d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_residuals_predicted_scatter_", trait, "_allmodels_test.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

}
