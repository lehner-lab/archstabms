
#' archstabms_ddPCA_biophysical_model
#'
#' Biophysical model from shallow ddPCA data.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param temperature temperature in degrees celcuis (default:30)
#' @param ddpca_outpath path to ddPCA data (required)
#' @param position_offset residue position offset (required)
#' @param load_mochi_data python (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_ddPCA_biophysical_model <- function(
  mochi_outpath,
  temperature = 30,
  ddpca_outpath,
  position_offset,
  load_mochi_data,
  outpath
  ){

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: archstabms_ddPCA_biophysical_model for", domain_name), "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Load model data and coerce to data.table
  mochi_dict <- load_mochi_data(file.path(mochi_outpath, "task_1", "saved_models", "data.pbz2"))
  mm <- data.table(py_to_r(mochi_dict[['Xohi']]))
  mm[, fitness := data.table(py_to_r(mochi_dict[['fitness']])[,'fitness'])]

  #Variant table
  vtable <- data.table(py_to_r(mochi_dict[['vtable']]))
  vtable[, paste0(mochi_dict[['phenotype_names']]) := 1]
  vtable_names <- names(vtable)
  vtable <- vtable[,.SD,,.SDcols = (!names(vtable) %in% c("WT", "nt_seq"))]
  vtable[, WT := NA]
  vtable[vtable[,Nham_aa]==0, WT := T]
  vtable[, nt_seq := NA]
  vtable <- vtable[,.SD,,.SDcols = vtable_names]

  #CV groups
  cvgroups <- data.table(py_to_r(mochi_dict[['cvgroups']]))
  cvgroups[vtable[,WT]==T, fold := 0]
  fold <- unlist(cvgroups[,fold])
  cvgroups <- cvgroups[,.SD,,.SDcols = names(cvgroups)!="fold"]
  cvgroups[, fold := fold]

  #ddPCA linear model - fitness
  #Load energies
  energies <- fread(file.path(ddpca_outpath, "model_weights_0.txt"))[!is.na(folding_coefficient),.(id, folding_coefficient)]
  energies[, WT_AA := substr(id, 1, 1)]
  energies[id!="WT", Mut := substr(id, nchar(id), nchar(id))]
  energies[id!="WT", Pos := as.integer(substr(id, 2, nchar(id)-1))-position_offset]
  energies[id!="WT", id := paste0(WT_AA, Pos, Mut)]

  mm_ddpca <- mm[apply(mm[,.SD,,.SDcols = !grepl('fitness', names(mm))], 1, sum)==1,]
  mm_ddpca_list <- list()
  for(i in names(mm)[!grepl("fitness", names(mm))]){
    if(i=='WT'){
      mm_ddpca_list[[i]] <- copy(mm_ddpca)
      mm_ddpca_list[[i]][, fitness := energies[id==i,folding_coefficient]]
      mm_ddpca_list[[i]][, paste0(i) := 1]
    }else if(energies[id==i,.N]!=0){
      mm_ddpca_list[[i]] <- copy(mm_ddpca)
      mm_ddpca_list[[i]][, fitness := energies[id==i,folding_coefficient]+energies[id=='WT',folding_coefficient]]
      mm_ddpca_list[[i]][, paste0(i) := 1]
    }
  }
  mm_ddpca <- rbindlist(mm_ddpca_list)
  #Fit model
  ddpca_lm <- lm(fitness ~.+0,data = mm_ddpca)

  #Model results
  coef_dt <- data.table()
  for(i in 1:10){
    #Predict
    vtable[, paste0('fold_', i) := archstabms__predict_fitness_ddPCA(mochi_outpath = ddpca_outpath, folding_energy = predict(ddpca_lm, newdata = mm), binding_energy = NULL, RT = 1)[["fitness_folding"]]]
    #Coefficients
    if(length(coef_dt)==0){
      coef_dt <- data.table(
        "id" = names(ddpca_lm$coefficients),
        "fold" = unlist(ddpca_lm$coefficients))
      names(coef_dt)[2] <- paste0("fold_", i)
    }else{
      coef_dt[, paste0("fold_", i) := unlist(ddpca_lm$coefficients)]
    }
  }

  #Add remaining columns to vtable
  vtable[, mean := rowMeans(.SD),,.SDcols = paste0("fold_", 1:10)]
  vtable[, std := apply(.SD, 1, sd),,.SDcols = paste0("fold_", 1:10)]
  vtable[, ci95 := std*1.96*2]
  vtable[, Fold := cvgroups[, fold]]
  vtable[WT==T, Fold := NA]

  #Add remaining columns to coef_dt
  coef_dt[, id_ref := id]
  coef_dt[, Pos := gsub("[a-zA-Z]", "", id)]
  coef_dt[, Pos_ref := Pos]
  coef_dt[, n := 10]
  coef_dt[, mean := rowMeans(.SD),,.SDcols = paste0("fold_", 1:10)]
  coef_dt[, std := apply(.SD, 1, sd),,.SDcols = paste0("fold_", 1:10)]
  coef_dt[, ci95 := std*1.96*2]
  coef_dt[, trait_name := 'Folding']
  coef_dt[, paste0("mean_kcal/mol") := mean*RT]
  coef_dt[, paste0("std_kcal/mol") := std*RT]
  coef_dt[, paste0("ci95_kcal/mol") := ci95*RT]

  #Position dictionary
  all_positions <- coef_dt[, unlist(lapply(strsplit(Pos, "_"), as.integer))]
  all_positions <- unique(all_positions[order(all_positions)])
  position_dict <- as.list(all_positions+position_offset)
  names(position_dict) <- all_positions

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

  #Reorder
  coef_dt <- coef_dt[,.SD,,.SDcols = c(
    "id", "id_ref", "Pos", "Pos_ref", paste0("fold_", 1:10),
    "n", "mean", "std", "ci95", "trait_name", "mean_kcal/mol", "std_kcal/mol", "ci95_kcal/mol")]

  #Save
  dir.create(file.path(outpath, "task_1", "predictions"), recursive = T)
  dir.create(file.path(outpath, "task_1", "weights"), recursive = T)
  write.table(vtable, file = file.path(outpath, "task_1", "predictions", "predicted_phenotypes_all.txt"), sep = "\t", row.names = F, quote = F, na = "")
  write.table(coef_dt, file = file.path(outpath, "task_1", "weights", "weights_Folding.txt"), sep = "\t", row.names = F, quote = F, na = "")
}
