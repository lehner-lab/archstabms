% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archstabms__plot_growthrate_violins.R
\name{archstabms__plot_growthrate_violins}
\alias{archstabms__plot_growthrate_violins}
\title{archstabms__plot_growthrate_violins}
\usage{
archstabms__plot_growthrate_violins(
  input_dt,
  report_outpath,
  dataset_name,
  trait,
  suffix = "",
  threshold,
  colour_scheme
)
}
\arguments{
\item{input_dt}{data.table with fitness estimates (required)}

\item{report_outpath}{output path for scatterplots (required)}

\item{dataset_name}{dataset name and suffix for output files (required)}

\item{trait}{trait name and suffix for output files (required)}

\item{suffix}{suffix for output files (default:"")}

\item{threshold}{growth rate threshold (required)}

\item{colour_scheme}{colour scheme file (required)}
}
\value{
Nothing
}
\description{
Plot growth rate violins.
}
