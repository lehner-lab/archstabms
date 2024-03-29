
#' archstabms__get_model_results
#'
#' Get model results and coefficients.
#'
#' @param input_folder path to input folder with model results (required)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return data.table with model results
#' @export
#' @import data.table
archstabms__get_model_results <- function(
  input_folder, 
  RT = 0.001987*(273+24)
  ){

  #Final model task folder
  task_folder <- ""
  for(i in list.files(input_folder, pattern = "task")){
    if(file.exists(file.path(input_folder, i, 'predictions', "predicted_phenotypes_all.txt"))){
      task_folder <- i
      break
    }    
  }

  #Predicted phenotype
  pred_dt <- fread(file.path(input_folder, task_folder, 'predictions', "predicted_phenotypes_all.txt"))
  #Mutation order
  pred_dt[, mut_order := Nham_aa]
  #Predicted phenotype for held out test data
  for(i in 1:10){pred_dt[Fold==i, fitness_pred := .SD[[1]],,.SDcols = paste0('fold_', i)]}

  #Model coefficients - folding
  coef_dt <- fread(file.path(input_folder, task_folder, 'weights', "weights_Folding.txt"), colClasses = list('character' = c('Pos', 'Pos_ref')))
  #Model coefficients - binding
  if(file.exists(file.path(input_folder, task_folder, 'weights', "weights_Binding.txt"))){
    coef_dt <- rbind(
      coef_dt,
      fread(file.path(input_folder, task_folder, 'weights', "weights_Binding.txt"), colClasses = list('character' = c('Pos', 'Pos_ref'))))
  }

  #Initial model coefficients
  if(task_folder!="task_1"){
    #Model coefficients - folding
    temp_dt <- fread(file.path(input_folder, "task_1", 'weights', "weights_Folding.txt"), colClasses = list('character' = c('Pos', 'Pos_ref')))
    #Model coefficients - binding
    if(file.exists(file.path(input_folder, "task_1", 'weights', "weights_Binding.txt"))){
      temp_dt <- rbind(
        temp_dt,
        fread(file.path(input_folder, "task_1", 'weights', "weights_Binding.txt"), colClasses = list('character' = c('Pos', 'Pos_ref'))))
    }
    temp_dt <- temp_dt[!Pos_ref %in% coef_dt[,Pos_ref],.SD,,.SDcols = c("id", "id_ref", "Pos", "Pos_ref", "trait_name")]
    coef_dt <- rbind(coef_dt, temp_dt, fill = T)
  }

  #Coefficient order
  coef_dt[, coef_order := sapply(strsplit(Pos, "_"), length)]

  return(list(
    'pred' = pred_dt, 
    'coef' = coef_dt))

  # #Mutated residues
  # input_dt[, mut1 := sapply(strsplit(id, ","), '[', 1)]
  # input_dt[, mut2 := sapply(strsplit(id, ","), '[', 2)]
  # input_dt[, mut_pos1 := as.integer(substr(mut1, 2, nchar(mut1)-1))]
  # input_dt[, mut_pos2 := as.integer(substr(mut2, 2, nchar(mut2)-1))]

  # #Load predicted fitness
  # pred_list <- list()
  # for(i in c(0, 1:n_bootstraps)){
  #   prefix <- ""
  #   if(i!=0){prefix <- file.path("bootstrap", paste0(bootstrap_prefix, i))}
  #   temp_dt <- fread(file.path(input_folder, prefix, "predicted_fitness_0.txt"))
  #   temp_dt[, id := input_dt[,id]]
  #   temp_dt[, dataset_binding := input_dt[,dataset_binding]]
  #   temp_dt[, f_dg_pred := additive_trait_folding*RT]
  #   temp_dt[, f_ddg_pred := f_dg_pred-temp_dt[id=="-0-",f_dg_pred][1]]
  #   temp_dt[, b_dg_pred := additive_trait_binding*RT]
  #   temp_dt[, b_ddg_pred := b_dg_pred-temp_dt[id=="-0-",b_dg_pred][1]]
  #   pred_list[[as.character(i)]] <- temp_dt
  # }
  # #Final predictions
  # pred_dt <- pred_list[[1]]

  # #Add bootstrap standard deviations
  # boot_dt <- do.call("cbind", pred_list[2:(n_bootstraps+1)])
  # pred_dt[, predicted_fitness_sd := boot_dt[,apply(.SD, 1, sd),,.SDcols = names(boot_dt)[grep("predicted_fitness", names(boot_dt))]]]
  # pred_dt[, f_ddg_pred_sd := boot_dt[,apply(.SD, 1, sd),,.SDcols = names(boot_dt)[grep("f_ddg_pred", names(boot_dt))]]]
  # pred_dt[, b_ddg_pred_sd := boot_dt[,apply(.SD, 1, sd),,.SDcols = names(boot_dt)[grep("b_ddg_pred", names(boot_dt))]]]
  # pred_dt[, mut_order := nchar(gsub("0", "", seq))]

  # #Remove variants with missing coefficients
  # pred_dt <- pred_dt[!(b_ddg_pred==0 & b_ddg_pred_sd==0 & id!="-0-")]

  # #Load all model weights (dGs)
  # dg_list <- list()
  # for(i in c(0, 1:n_bootstraps)){
  #   prefix <- ""
  #   if(i!=0){prefix <- file.path("bootstrap", paste0(bootstrap_prefix, i))}
  #   temp_dt <- fread(file.path(input_folder, prefix, "predicted_fitness_0.txt"))[input_dt[,id]=="-0-"]
  #   WT_coef_folding <- temp_dt[,additive_trait_folding][1]
  #   WT_coef_binding <- temp_dt[,additive_trait_binding][1]
  #   temp_dt <- fread(file.path(input_folder, prefix, "model_weights_0.txt"))
  #   temp_dt[, f_dg_pred := (folding_coefficient + WT_coef_folding)*RT]
  #   temp_dt[, f_ddg_pred := folding_coefficient*RT]
  #   temp_dt[, b_dg_pred := (binding_coefficient + WT_coef_binding)*RT]
  #   temp_dt[, b_ddg_pred := binding_coefficient*RT]
  #   dg_list[[as.character(i)]] <- temp_dt
  # }
  # #Final coefficients
  # coef_dt <- dg_list[[1]]

  # #Add bootstrap standard deviations
  # boot_dt <- do.call("cbind", dg_list[2:(n_bootstraps+1)])
  # coef_dt[, f_ddg_pred_sd := apply(boot_dt[,.SD,,.SDcols = grep("f_ddg_pred", names(boot_dt))], 1, sd)]
  # coef_dt[, b_ddg_pred_sd := apply(boot_dt[,.SD,,.SDcols = grep("b_ddg_pred", names(boot_dt))], 1, sd)]

  # #Coefficient ids
  # if(!"id" %in% names(coef_dt)){
  #   coef_names <- unlist(fread(file.path(input_folder, "", "feature_matrix_colnames.txt"), header = F))
  #   coef_dt[, name := gsub("pos", "", coef_names)]
  #   all_mut_pos <- unique(unlist(as.data.frame(input_dt[id!="-0-",.(mut_pos1, mut_pos2)])))
  #   all_mut_pos <- all_mut_pos[order(all_mut_pos)]
  #   all_mut_pos <- rev(all_mut_pos[!is.na(all_mut_pos)])
  #   coef_dt[name!="(Intercept)", Pos := all_mut_pos[as.integer(sapply(lapply(strsplit(name, "\\."), unlist), "[", 1))]]
  #   coef_dt[name!="(Intercept)", Mut := sapply(lapply(strsplit(name, "\\."), unlist), "[", 2)]
  #   coef_dt[name!="(Intercept)", WT_AA := sapply(as.list(Pos), function(x){substr(input_dt[id=="-0-", aa_seq][1], x, x)})]
  #   coef_dt[name!="(Intercept)", id := paste0(WT_AA, Pos, Mut)]
  #   coef_dt[name=="(Intercept)", id := "-0-"]
  #   coef_dt[, mut_order := 1]
  #   coef_dt[id == "-0-", mut_order := 0]
  # }else{
  #   coef_dt[id == "WT", id := "-0-"]
  #   coef_dt[, mut_order := 1]
  #   coef_dt[id == "-0-", mut_order := 0]
  # }
 
  # #Append to predicted fitness data.table
  # pred_dt <- rbind(pred_dt, coef_dt[!id %in% pred_dt[,id],.(f_dg_pred, f_ddg_pred, f_ddg_pred_sd, b_dg_pred, b_ddg_pred, b_ddg_pred_sd, id, mut_order)], fill = T)
  # return(pred_dt)
}
