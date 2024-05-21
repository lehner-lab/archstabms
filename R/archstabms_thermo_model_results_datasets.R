
#' archstabms_thermo_model_results_datasets
#'
#' Evaluate thermo model results of multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param literature_free_energies_list list of paths to literature free energies (default:list())
#' @param ddPCA_free_energies list of paths to ddPCA free energies (default:list())
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
archstabms_thermo_model_results_datasets <- function(
  dataset_names,
  literature_free_energies_list = list(),
  ddPCA_free_energies_list = list(),
  base_dir,
  output_dir,
  stagenum,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  for(i in dataset_names){
    #order 1
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o1")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 1 -linear
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, paste0("001_archstabms_linear_model_", i, "o1")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1-linear"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 1 -ddPCA single background linear
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, paste0("001_archstabms_ddPCA_singlebackground_linear_model_", i, "o1")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1-ddPCAsinglebackgroundlinear"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 1 -ddPCA linear
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, paste0("001_archstabms_ddPCA_linear_model_", i, "o1")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1-ddPCAlinear"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 1 -ddPCA biophysical
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, paste0("001_archstabms_ddPCA_biophysical_model_", i, "o1")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1-ddPCA"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 2
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o2")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o2"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 2 -linear
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, paste0("001_archstabms_linear_model_", i, "o2")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o2-linear"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 3
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o3")),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o3"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
    #order 3 -linear
    archstabms_thermo_model_results(
      mochi_outpath = file.path(base_dir, paste0("001_archstabms_linear_model_", i, "o3")),
      ddPCA_free_energies = ddPCA_free_energies_list[[i]],
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o3-linear"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)

    # #order 1 - ensemble
    # archstabms_thermo_model_results(
    #   mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o1-ensemble")),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme,
    #   execute = execute)
    # #order 1 -linear - ensemble
    # archstabms_thermo_model_results(
    #   mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o1-linear-ensemble")),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o1-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme,
    #   execute = execute)
    # #order 2 - ensemble
    # archstabms_thermo_model_results(
    #   mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o2-ensemble")),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o2-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme,
    #   execute = execute)
    # #order 2 -linear - ensemble
    # archstabms_thermo_model_results(
    #   mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o2-linear-ensemble")),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o2-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme,
    #   execute = execute)
    # #order 3 - ensemble
    # archstabms_thermo_model_results(
    #   mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o3-ensemble")),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o3-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme,
    #   execute = execute)
    # #order 3 -linear - ensemble
    # archstabms_thermo_model_results(
    #   mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o3-linear-ensemble")),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o3-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme,
    #   execute = execute)
  }

  ### 2nd order model comparisons
  ###########################

  outpath <- archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_datasets"), stagenum=stagenum, base_dir=output_dir)
  
  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Load data
  i <- "CM6"
  #CM6 2nd order biophysical model
  coef_dt <- fread(file.path(archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o2"), stagenum=stagenum, base_dir=output_dir), "model_coefficients.txt"))
  coef_dt[, model := "biophysical"]
  #CM6 2nd order linear model
  coef_dt_linear <- fread(file.path(archstabms__format_dir(dir_suffix=paste0("_archstabms_thermo_model_results_", i, "o2-linear"), stagenum=stagenum, base_dir=output_dir), "model_coefficients.txt"))
  coef_dt_linear[, model := "linear"]
  #Merge
  coef_dt <- rbind(
    coef_dt[,.(id, coef_order, mean, std, model)],
    coef_dt_linear[,.(id, coef_order, mean, std, model)])
  #Significance
  coef_dt[, nlog10p_value := -log10(pnorm(abs(mean), 0, std, lower.tail = F))]
  # theoretical <- -log10(rank(p.values)/length(p.values))

  #Cumulative distribution of P-values
  plot_dt <- coef_dt[coef_order==2]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(nlog10p_value, color = model)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::theme_classic() +
    ggplot2::xlab("-log10(P-value)") +
    ggplot2::ylab("ECDF")
  d <- d + ggplot2::scale_colour_manual(values = c(unlist(colour_scheme[["shade 0"]][[1]]), 'black'))
  ggplot2::ggsave(file.path(outpath, paste0("coupling_ecdf_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)

}
