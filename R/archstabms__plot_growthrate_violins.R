
#' archstabms__plot_growthrate_violins
#'
#' Plot growth rate violins.
#'
#' @param input_dt data.table with fitness estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param trait trait name and suffix for output files (required)
#' @param suffix suffix for output files (default:"")
#' @param threshold growth rate threshold (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__plot_growthrate_violins <- function(
  input_dt,
  report_outpath,
  dataset_name,
  trait,
  suffix="",
  threshold,
  colour_scheme
  ){

  #Abundance growthrate violins with hamming distance
  plot_dt <- copy(input_dt)
  plot_dt[, mut_order := as.factor(Nham_aa)]
  #Maximum mut_order in violin plot
  text_dt <- plot_dt[is.na(WT),.(label = round(sum(growthrate>threshold)/.N*100, 0), growthrate = 0.17),mut_order]
  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=plot_dt[, max(Nham_aa)]))
  # d <- ggplot2::ggplot(plot_dt[is.na(WT) & Nham_aa<29],ggplot2::aes(mut_order, growthrate, fill = mut_order, color = mut_order)) +
  d <- ggplot2::ggplot(plot_dt[is.na(WT)],ggplot2::aes(mut_order, growthrate, fill = mut_order, color = mut_order)) +
    ggplot2::geom_violin() +
    ggplot2::scale_colour_manual(values = cc) +
    ggplot2::scale_fill_manual(values = cc) +
    ggplot2::geom_hline(yintercept = plot_dt[WT==T,growthrate]) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Growth rate")
  ggplot2::ggsave(file.path(report_outpath, paste0("growthrate_violins_", dataset_name, "_", trait, "_Nham_aa", suffix, ".pdf")), d, width = 8, height = 3, useDingbats=FALSE)
  ggplot2::ggsave(file.path(report_outpath, paste0("growthrate_violins_", dataset_name, "_", trait, "_Nham_aa_big", suffix, ".pdf")), d, width = 10, height = 7, useDingbats=FALSE)
  ggplot2::ggsave(file.path(report_outpath, paste0("growthrate_violins_", dataset_name, "_", trait, "_Nham_aa_small", suffix, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Abundance growthrate histograms vs hamming distance
  plot_dt <- input_dt[is.na(WT),.(count = .N), Nham_aa]
  plot_dt[, mut_order := as.factor(Nham_aa)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, count)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes()) +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Frequency") +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("growthrate_violins_", dataset_name, "_", trait, "_Nham_aa", suffix, "_hist.pdf")), d, width = 8, height = 2, useDingbats=FALSE)
  ggplot2::ggsave(file.path(report_outpath, paste0("growthrate_violins_", dataset_name, "_", trait, "_Nham_aa_big", suffix, "_hist.pdf")), d, width = 10, height = 2, useDingbats=FALSE)
  ggplot2::ggsave(file.path(report_outpath, paste0("growthrate_violins_", dataset_name, "_", trait, "_Nham_aa_small", suffix, "_hist.pdf")), d, width = 4, height = 2, useDingbats=FALSE)

}
