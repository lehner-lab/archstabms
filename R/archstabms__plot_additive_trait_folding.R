
#' archstabms__plot_additive_trait_folding
#'
#' Plot folding additive trait.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param RT constant (default:0.001987*(273+24))
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__plot_additive_trait_folding <- function(
  mochi_outpath,
  input_dt,
  RT = 0.001987*(273+24),
  report_outpath,
  dataset_name,
  colour_scheme
  ){

  #Plot data.table
  plot_dt <- copy(input_dt)

  #Return if no global epistasis
  if(!'fold_1_additive_trait0' %in% names(plot_dt)){
    return()
  }

  #Folding dG
  plot_dt[, f_dg_pred := fold_1_additive_trait0*RT]
  plot_dt[, f_ddg_pred := f_dg_pred-plot_dt[mut_order==0,f_dg_pred]]

  #Model data.table (for geom_line)
  folding_energy_range <- plot_dt[mut_order>0,range(f_dg_pred)]
  folding_energy_grid <- seq(folding_energy_range[1], folding_energy_range[2], (folding_energy_range[2]-folding_energy_range[1])/500)
  #Predicted fitness
  pred_fitness_list <- archstabms__predict_fitness(
    mochi_outpath = mochi_outpath,
    dataset_name = dataset_name,
    folding_energy = folding_energy_grid,
    binding_energy = NA,
    RT = RT)
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    f_dg_pred = rep(folding_energy_grid, 2),
    f_ddg_pred = rep(folding_energy_grid, 2)-plot_dt[mut_order==0,f_dg_pred][1],
    fitness = pred_fitness_list[["fitness_folding"]],
    mut_order = rep(c(1, 2), each = length(folding_energy_grid)))
  
  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(f_dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab(expression("Folding "*Delta*"G (inferred)")) +
    ggplot2::ylab("Abundance fitness\n(observed)") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,f_dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_folding_binhex_", dataset_name, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(f_dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab(expression("Folding "*Delta*"G (inferred)")) +
    ggplot2::ylab("Abundance fitness\n(observed)") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,f_dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_folding_binhex_", dataset_name, "_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(f_dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab(expression("Folding "*Delta*"G (inferred)")) +
    ggplot2::ylab("Abundance fitness\n(observed)") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,f_dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_folding_binhex_", dataset_name, "_highres.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(f_dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab(expression("Folding "*Delta*"G (inferred)")) +
    ggplot2::ylab("Abundance fitness\n(observed)") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,f_dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_folding_binhex_", dataset_name, "_highres_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  # #Binhex no facet - xlim [-3, 6]
  # plot_xlim <- c(-3, 6)
  # d <- ggplot2::ggplot(plot_dt[mut_order>0 & f_dg_pred>plot_xlim[1] & f_dg_pred<plot_xlim[2]],ggplot2::aes(f_dg_pred, fitness)) +
  #   ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
  #   ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
  #   ggplot2::xlab(expression(Delta*"G Folding (inferred)")) +
  #   ggplot2::ylab("Fitness (Abundance)") +
  #   ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  #   ggplot2::geom_vline(xintercept = 0, linetype = 2) +
  #   ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,f_dg_pred][1], color = "grey", linetype = 2) +
  #   ggplot2::geom_line(data = pred_fitness_dt[mut_order==1 & f_dg_pred>plot_xlim[1] & f_dg_pred<plot_xlim[2]], color = colour_scheme[["shade 0"]][[1]]) +
  #   ggplot2::coord_cartesian(xlim = plot_xlim) +
  #   ggplot2::theme_classic()
  # ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_folding_binhex_xlim_", dataset_name, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  # #Binhex no facet - xlim [-2, 7]
  # plot_xlim <- c(-2, 7)
  # d <- ggplot2::ggplot(plot_dt[mut_order>0 & f_ddg_pred>plot_xlim[1] & f_ddg_pred<plot_xlim[2]],ggplot2::aes(f_ddg_pred, fitness)) +
  #   ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
  #   ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
  #   ggplot2::xlab(expression(Delta*Delta*"G Folding (inferred)")) +
  #   ggplot2::ylab("Fitness (Abundance)") +
  #   ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  #   ggplot2::geom_vline(xintercept = 0, linetype = 2) +
  #   ggplot2::geom_line(data = pred_fitness_dt[mut_order==1 & f_ddg_pred>plot_xlim[1] & f_ddg_pred<plot_xlim[2]], color = colour_scheme[["shade 0"]][[1]]) +
  #   ggplot2::coord_cartesian(xlim = plot_xlim) +
  #   ggplot2::theme_classic()
  # ggplot2::ggsave(file.path(report_outpath, paste0("ddG_observed_folding_binhex_xlim_", dataset_name, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  # #Scatter facet highlighting confidenct dGs
  # d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(f_dg_pred, fitness, colour = f_ddg_pred_conf)) +
  #   ggplot2::geom_point(data = plot_dt[f_ddg_pred_conf==F & mut_order>0], size = 1, alpha = 1/4) +
  #   ggplot2::geom_point(data = plot_dt[f_ddg_pred_conf==T & mut_order>0], size = 1, alpha = 1/4) +
  #   ggplot2::xlab(expression(Delta*"G Folding (inferred)")) +
  #   ggplot2::ylab("Fitness (Abundance)") +
  #   ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  #   ggplot2::geom_vline(xintercept = 0, linetype = 2) +
  #   ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,f_dg_pred][1], color = "grey", linetype = 2) +
  #   ggplot2::geom_line(data = pred_fitness_dt, color = colour_scheme[["shade 0"]][[1]]) +
  #   ggplot2::labs(color = "Confident\nfree energy") +
  #   ggplot2::facet_wrap(~mut_order, nrow = 1) +
  #   ggplot2::theme_classic()
  # if(!is.null(colour_scheme)){
  #   d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  # }
  # ggplot2::ggsave(file.path(report_outpath, "dG_observed_folding_scatter_conf_facet.pdf"), d, width = 8, height = 4, useDingbats=FALSE)

  # #Scatter facet highlighting confidenct ddGs
  # d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(f_ddg_pred, fitness, colour = f_ddg_pred_conf)) +
  #   ggplot2::stat_binhex(bins = 20, size = 0) +
  #   ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
  #   ggplot2::xlab(expression(Delta*Delta*"G Folding (inferred)")) +
  #   ggplot2::ylab("Fitness (Abundance)") +
  #   ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  #   ggplot2::geom_vline(xintercept = 0, linetype = 2) +
  #   ggplot2::geom_line(data = pred_fitness_dt, color = colour_scheme[["shade 0"]][[1]]) +
  #   ggplot2::labs(color = "Confident\nfree energy") +
  #   ggplot2::facet_wrap(f_ddg_pred_conf~mut_order) +
  #   ggplot2::theme_classic()
  # if(!is.null(colour_scheme)){
  #   d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  # }
  # ggplot2::ggsave(file.path(report_outpath, "ddG_observed_folding_binhex_conf_facet.pdf"), d, width = 8, height = 4, useDingbats=FALSE)

}
