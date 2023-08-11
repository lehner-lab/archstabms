
#' archstabms__plot_ddPCAvalidation_scatter
#'
#' Plot correlation with validation data.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param ddPCA_inpath path to ddPCA free energies (required)
#' @param trait trait name to plot: folding or binding (required)
#' @param report_outpath output path for scatterplots (required)
#' @param highlight_colour colour for highlights (default:red)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__plot_ddPCAvalidation_scatter <- function(
  input_dt,
  ddPCA_inpath,
  trait,
  report_outpath,
  highlight_colour = "red",
  RT = 0.001987*(273+24)
  ){

  #Check validation data exists
  if(is.na(ddPCA_inpath)){return()}

  #Load data
  coef_dt_ddpca <- fread(ddPCA_inpath)

  #Merge
  plot_dt <- merge(
    input_dt[trait_name == trait & id != "WT"], 
    coef_dt_ddpca[,.(id_ref = id, ddg_pred = .SD[[1]]*RT),,.SDcols = paste0(tolower(trait), "_coefficient")], 
    id = "id_ref")

  #Plot folding ddG validation scatter
  plot_dt[, ddg_pred_comb := .SD[[1]],,.SDcols = "mean_kcal/mol"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ddg_pred, ddg_pred_comb)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_abline(lty = 2) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(ddg_pred, ddg_pred_comb), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1), color = "black") +
    # ggrepel::geom_text_repel(ggplot2::aes(label = id_ref), show.legend = T, 
    #   max.overlaps = Inf) +
    ggplot2::theme_classic() +
    ggplot2::xlab(bquote(.(trait) ~ Delta*Delta*"G (kcal/mol) - ddPCA")) +
    ggplot2::ylab(bquote(.(trait) ~ Delta*Delta*"G (kcal/mol)"))
  ggplot2::ggsave(file.path(report_outpath, paste0("ddG_ddpca_scatter_order1_", trait, ".pdf")), d, width = 3, height = 3, useDingbats=FALSE)

}
