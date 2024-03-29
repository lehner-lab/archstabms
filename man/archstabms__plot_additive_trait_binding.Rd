% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archstabms__plot_additive_trait_binding.R
\name{archstabms__plot_additive_trait_binding}
\alias{archstabms__plot_additive_trait_binding}
\title{archstabms__plot_additive_trait_binding}
\usage{
archstabms__plot_additive_trait_binding(
  mochi_outpath,
  input_dt,
  RT = 0.001987 * (273 + 24),
  report_outpath,
  dataset_name,
  colour_scheme
)
}
\arguments{
\item{mochi_outpath}{path to MoCHI thermo model fit results (required)}

\item{input_dt}{data.table with model free energy estimates (required)}

\item{RT}{constant (default:0.001987*(273+24))}

\item{report_outpath}{output path for scatterplots (required)}

\item{dataset_name}{dataset name and suffix for output files (required)}

\item{colour_scheme}{colour scheme file (default:ggplot colours)}
}
\value{
Nothing
}
\description{
Plot folding additive trait.
}
