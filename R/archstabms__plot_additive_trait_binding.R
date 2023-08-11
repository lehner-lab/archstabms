
#' archstabms__plot_additive_trait_binding
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
archstabms__plot_additive_trait_binding <- function(
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
  if(!'fold_1_additive_trait1' %in% names(plot_dt)){
    return()
  }

  #Folding dG
  plot_dt[, f_dg_pred := fold_1_additive_trait0*RT]
  plot_dt[, f_ddg_pred := f_dg_pred-plot_dt[mut_order==0,f_dg_pred]]

  #Binding dG
  plot_dt[, b_dg_pred := fold_1_additive_trait1*RT]
  plot_dt[, b_ddg_pred := b_dg_pred-plot_dt[mut_order==0,b_dg_pred]]

  #WT dGs
  f_dg_WT <- plot_dt[mut_order==0,f_dg_pred][1]
  b_dg_WT <- plot_dt[mut_order==0,b_dg_pred][1]

  #Number of grid points
  num_grid <- 15

  ### All data
  ###########################

  #Model data.table (for geom_line)
  folding_energy_range <- plot_dt[mut_order>0,range(f_dg_pred)]
  binding_energy_range <- plot_dt[mut_order>0,range(b_dg_pred)]
  folding_energy_grid <- seq(folding_energy_range[1], folding_energy_range[2], (folding_energy_range[2]-folding_energy_range[1])/num_grid)
  binding_energy_grid <- seq(binding_energy_range[1], binding_energy_range[2], (binding_energy_range[2]-binding_energy_range[1])/num_grid)
  
  energy_grid_dt <- as.data.table(expand.grid(folding_energy_grid = folding_energy_grid, binding_energy_grid = binding_energy_grid))

  #Predicted fitness
  pred_fitness_list <- archstabms__predict_fitness(
    mochi_outpath = mochi_outpath,
    dataset_name = dataset_name,
    folding_energy = energy_grid_dt[,folding_energy_grid],
    binding_energy = energy_grid_dt[,binding_energy_grid],
    RT = RT)
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    f_dg_pred = energy_grid_dt[,folding_energy_grid],
    b_dg_pred = energy_grid_dt[,binding_energy_grid],
    fitness = pred_fitness_list[["fitness_binding"]])

  Cairo::CairoPDF(file = file.path(report_outpath, "dG_observed_binding_scatter.pdf"))
  plot3D::persp3D(
    x = folding_energy_grid, 
    y = binding_energy_grid, 
    z = matrix(data=unlist(pred_fitness_dt[,fitness]), nrow=length(folding_energy_grid), ncol=length(folding_energy_grid)), 
    r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
    xlab = "dG Folding",
    ylab = "dG Binding",
    zlab = "Fitness (Binding)")
  plot3D::scatter3D(
    x = plot_dt[,f_dg_pred], 
    y = plot_dt[,b_dg_pred], 
    z = plot_dt[,fitness], 
    add = T, col = "black", alpha = 0.2, cex = 0.2)
  dev.off()

  Cairo::CairoPDF(file = file.path(report_outpath, "ddG_observed_binding_scatter.pdf"))
  plot3D::persp3D(
    x = folding_energy_grid-f_dg_WT, 
    y = binding_energy_grid-b_dg_WT, 
    z = matrix(data=pred_fitness_dt[,fitness], nrow=length(folding_energy_grid), ncol=length(folding_energy_grid)), 
    r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
    xlab = "ddG Folding",
    ylab = "ddG Binding",
    zlab = "Fitness (Binding)")
  plot3D::scatter3D(
    x = plot_dt[,f_dg_pred]-f_dg_WT, 
    y = plot_dt[,b_dg_pred]-b_dg_WT, 
    z = plot_dt[,fitness], 
    add = T, col = "black", alpha = 0.2, cex = 0.2)
  dev.off()

}
