
#' archstabms__predict_fitness
#'
#' Predict fitness from folding and binding energies for 3-state model.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param dataset_name dataset name  (required)
#' @param folding_energy folding energies for 3-state model (required)
#' @param binding_energy binding energies for 3-state model (required for binding fitness prediction)
#' @param RT constant (default:0.001987*(273+24))
#' @param phenotype_names Phenotype names (default:c("Abundance", "Binding"))
#'
#' @return A list with folding and binding fitness predictions
#' @export
#' @import data.table
archstabms__predict_fitness <- function(
  mochi_outpath,
  dataset_name,
  folding_energy,
  binding_energy,
  RT = 0.001987*(273+24),
  phenotype_names = c("Abundance", "Binding")
  ){

  #Linear transformation parameters
  linears_file_abundance <- file.path(mochi_outpath, "task_1", 'weights', paste0("linears_weights_", phenotype_names[1], "_", dataset_name, ".txt"))
  linears_params_abundance <- as.list(fread(linears_file_abundance)[1,2:3])

  #Predicted fitness
  pred_list <- list()

  #Folding fitness
  if(!is.null(folding_energy)){
    pred_list[["fraction_folded"]] <- archstabms__fraction_folded(folding_energy, RT)
    pred_list[["fitness_folding"]] <- pred_list[["fraction_folded"]] * linears_params_abundance[["kernel"]] + linears_params_abundance[["bias"]]
  }

  #Binding fitness
  if(!is.null(folding_energy) & !is.null(binding_energy) & length(folding_energy) == length(binding_energy)){
    #Linear transformation parameters
    linears_file_binding <- file.path(mochi_outpath, "task_1", 'weights', paste0("linears_weights_", phenotype_names[2], "_", dataset_name, ".txt"))
    linears_params_binding <- as.list(fread(linears_file_binding)[1,2:3])
    #Binding fitness
    pred_list[["fraction_bound"]] <- archstabms__fraction_bound(folding_energy, binding_energy, RT)
    pred_list[["fitness_binding"]] <- pred_list[["fraction_bound"]] * linears_params_binding[["kernel"]] + linears_params_binding[["bias"]]
  }

  #Return
  return(pred_list)
}