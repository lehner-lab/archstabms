
#' archstabms_fit_linear_model
#'
#' Fit linear model to combinatorial data.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param position_offset residue position offset (required)
#' @param load_mochi_data python (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_fit_linear_model <- function(
  mochi_outpath,
  position_offset,
  load_mochi_data,
  outpath
  ){

  #Return if mochi_outpath doesn't exist
  if(!dir.exists(mochi_outpath)){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: archstabms_fit_linear_model for", domain_name), "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Load model data and coerce to data.table
  mochi_dict <- load_mochi_data(file.path(mochi_outpath, "task_1", "saved_models", "data.pbz2"))
  mm <- data.table(py_to_r(mochi_dict[['Xohi']]))
  mm[, fitness := data.table(py_to_r(mochi_dict[['fitness']])[,'fitness'])]

  #Number of genotypes per double
  if(sum(grepl("_", names(mm)))!=0){
    ng_doubles <- mm[,median(apply(.SD, 2, sum)),,.SDcols = grepl("_", names(mm))]
    print(paste0("Number of genotypes per double (", basename(mochi_outpath), "): ", ng_doubles))
  }

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

  #Fit linear models using 10-fold cross-validation
  coef_dt <- data.table()
  for(i in 1:10){
    #Fit model
    cv_lm <- lm(fitness ~.+0,data = mm[cvgroups[,fold]!=i,])
    #Predict
    vtable[, paste0('fold_', i) := predict(cv_lm, newdata = mm)]
    #Coefficients
    if(length(coef_dt)==0){
      coef_dt <- data.table(
        "id" = names(cv_lm$coefficients),
        "fold" = unlist(cv_lm$coefficients))
      names(coef_dt)[2] <- paste0("fold_", i)
    }else{
      coef_dt[, paste0("fold_", i) := unlist(cv_lm$coefficients)]
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
  coef_dt[, paste0("mean_kcal/mol") := mean]
  coef_dt[, paste0("std_kcal/mol") := std]
  coef_dt[, paste0("ci95_kcal/mol") := ci95]

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
