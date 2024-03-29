% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archstabms__secondary_structure_from_PDB.R
\name{archstabms__secondary_structure_from_PDB}
\alias{archstabms__secondary_structure_from_PDB}
\title{archstabms__secondary_structure_from_PDB}
\usage{
archstabms__secondary_structure_from_PDB(input_file, chain = "A")
}
\arguments{
\item{input_file}{path to PDB file (required)}

\item{chain}{chain id (default:A)}
}
\value{
data.table with secondary structure
}
\description{
Get secondary structure.
}
