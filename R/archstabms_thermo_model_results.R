
#' archstabms_thermo_model_results
#'
#' Evaluate thermo model results.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param temperature temperature in degrees celcuis (default:30)
#' @param literature_free_energies path to literature free energies (default:NA)
#' @param ddPCA_free_energies path to ddPCA free energies (default:NA)
#' @param position_offset residue position offset (default:0)
#' @param position_dict residue position dictionary (default:{})
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_thermo_model_results <- function(
  mochi_outpath,
  temperature = 30,
  literature_free_energies = NA,
  ddPCA_free_energies = NA,
  position_offset = 0,
  position_dict = list(),
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed or mochi_outpath doesn't exist
  if(!execute | !dir.exists(mochi_outpath)){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: archstabms_thermo_model_results for", domain_name), "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Load model results
  model_results <- archstabms__get_model_results(
    input_folder = mochi_outpath, 
    # input_dt = fitness_dt, 
    RT = RT)
  pred_dt <- model_results[['pred']]
  coef_dt <- model_results[['coef']]

  #Translate positions
  if(length(position_dict)!=0){
    for(i in rev(names(position_dict))){
      #Translate Pos_ref
      for(j in c(paste0('(^)', i, '($)'), paste0('(^)', i, '(_)'), paste0('(_)', i, '(_)'), paste0('(_)', i, '($)'))){
        coef_dt[, Pos_ref := gsub(j, paste0("\\1", position_dict[[i]], "\\2"), Pos_ref)]
      }
      #Translate id_ref
      coef_dt[, id_ref := gsub(paste0('([A-Za-z])', i, '([A-Za-z])'), paste0('\\1', position_dict[[i]], '\\2'), id_ref)]
    }
  }

  #Call confident ddGs
  coef_dt <- archstabms__define_confident_free_energies(
    input_dt = coef_dt, 
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]])

  #Plot model performance
  dataset_names <- names(pred_dt)[grepl("Abundance_|Binding_", names(pred_dt))]
  for(i in dataset_names){
    archstabms__plot_model_performance(
      input_dt = pred_dt[get(i)==1,], 
      report_outpath = outpath, 
      trait = sapply(strsplit(i, "_"), '[', 1),
      highlight_colour = colour_scheme[["shade 0"]][[1]])
  }

  #Plot 2D global epistasis (folding energy vs. folding fitness) - all phenotypes
  dataset_names <- names(pred_dt)[grepl("Abundance_", names(pred_dt))]
  for(i in dataset_names){
    archstabms__plot_additive_trait_folding(
      mochi_outpath = mochi_outpath,
      input_dt = pred_dt[get(i)==1,], 
      RT = RT,
      report_outpath = outpath,
      dataset_name = gsub("Abundance_", "", i),
      colour_scheme = colour_scheme)
  }

  #Plot 3D global epistasis (folding+binding energies vs. binding fitness) 
  dataset_names <- names(pred_dt)[grepl("Binding_", names(pred_dt))]
  for(i in dataset_names){
    archstabms__plot_additive_trait_binding(
      mochi_outpath = mochi_outpath,
      input_dt = pred_dt[get(i)==1,], 
      RT = RT,
      report_outpath = outpath,
      dataset_name = gsub("Binding_", "", i),
      colour_scheme = colour_scheme)
  }

  #Plot correlation with ddPCA data - folding ddGs
  archstabms__plot_ddPCAvalidation_scatter(
    input_dt = coef_dt, 
    ddPCA_inpath = ddPCA_free_energies, 
    trait = "Folding",
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]],
    RT = RT)

  #Plot correlation with ddPCA data - binding ddGs
  if("Binding" %in% coef_dt[,trait_name]){
    archstabms__plot_ddPCAvalidation_scatter(
      input_dt = coef_dt, 
      ddPCA_inpath = ddPCA_free_energies, 
      trait = "Binding",
      report_outpath = outpath, 
      highlight_colour = colour_scheme[["shade 0"]][[1]],
      RT = RT)
  }

  #Save dGs and ddGs
  write.table(pred_dt, 
    file = file.path(outpath, "model_predictions.txt"), 
    quote = F, sep = "\t", row.names = F)
  write.table(coef_dt, 
    file = file.path(outpath, "model_coefficients.txt"), 
    quote = F, sep = "\t", row.names = F)
}