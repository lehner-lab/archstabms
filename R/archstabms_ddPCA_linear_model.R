
#' archstabms_ddPCA_linear_model
#'
#' Linear model from from shallow ddPCA data.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param ddpca_outpath path to ddPCA data (required)
#' @param position_offset residue position offset (required)
#' @param load_mochi_data python (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_ddPCA_linear_model <- function(
  mochi_outpath,
  ddpca_outpath,
  position_offset,
  load_mochi_data,
  outpath
  ){

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: archstabms_ddPCA_linear_model for", domain_name), "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

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
  #Load variants
  load(ddpca_outpath)
  all_variants <- all_variants[STOP==F & STOP_readthrough==F]
  wt_seq <- all_variants[WT==T,unlist(strsplit(aa_seq, ""))]
  #Feature matrix
  feat_dt <- as.data.table(t(all_variants[, strsplit(aa_seq, "")]))
  #WT reference = 0
  for(i in 1:length(wt_seq)){
    feat_dt[get(paste0("V",i))==wt_seq[i], paste0("V",i) := "0"]
  }
  names(feat_dt) <- paste0(wt_seq, gsub("V", "", names(feat_dt)))
  feat_dt <- feat_dt[,.SD,,.SDcols = names(feat_dt)[apply(feat_dt!="0", 2, sum)!=0]]
  feat_dt[, fitness := all_variants[,fitness]]
  #Model matrix
  mm_ddpca <- model.matrix(fitness~., feat_dt)
  mm_ddpca <- data.table(mm_ddpca, fitness = all_variants[,fitness])
  names(mm_ddpca)[1] <- 'WT'
  #Fit model
  ddpca_lm <- lm(fitness ~.+0,data = mm_ddpca)

  #Format mm to match mm_ddpca
  names(mm)[2:(dim(mm)[2]-1)] <- paste0(
    substr(names(mm)[2:(dim(mm)[2]-1)], 1, 1),
    as.integer(substr(names(mm)[2:(dim(mm)[2]-1)], 2, nchar(names(mm)[2:(dim(mm)[2]-1)])-1))+position_offset,
    substr(names(mm)[2:(dim(mm)[2]-1)], nchar(names(mm)[2:(dim(mm)[2]-1)]), nchar(names(mm)[2:(dim(mm)[2]-1)])))
  mm_orig <- copy(mm)
  for(i in names(mm_ddpca)[!names(mm_ddpca) %in% names(mm)]){
    mm[, paste0(i) := 0]
  }

  #Model results
  coef_dt <- data.table()
  mm_predict <- predict(ddpca_lm, newdata = mm[,.SD,,.SDcols = names(mm_ddpca)])
  for(i in 1:10){
    #Predict
    vtable[, paste0('fold_', i) := mm_predict]
    #Coefficients
    if(length(coef_dt)==0){
      coef_dt <- data.table(
        "id" = names(ddpca_lm$coefficients)[names(ddpca_lm$coefficients) %in% names(mm_orig)],
        "fold" = unlist(ddpca_lm$coefficients)[names(ddpca_lm$coefficients) %in% names(mm_orig)])
      names(coef_dt)[2] <- paste0("fold_", i)
    }else{
      coef_dt[, paste0("fold_", i) := unlist(ddpca_lm$coefficients)[names(ddpca_lm$coefficients) %in% names(mm_orig)]]
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

  #Translate positions
  coef_dt[id!='WT', Pos := as.character(as.integer(Pos_ref)-position_offset)]
  coef_dt[id!='WT', id := gsub(Pos_ref, Pos, id),id]

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
