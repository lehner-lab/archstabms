Welcome to the GitHub repository for the following publication: [The genetic architecture of protein stability (Faure AJ et al., 2023)]()

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Table Of Contents

* **1. [Required Software](#required-software)**
* **2. [Required Data](#required-data)**
* **3. [Installation Instructions](#installation-instructions)**
* **4. [Usage](#usage)**

# Required Software

To run the archstabms pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/)** (Biostrings, Cairo, bio3d, data.table, GGally, ggplot2, ggrepel, hexbin, plot3D, ppcor, reshape2, reticulate, ROCR, scales)

# Required Data

Fitness scores, inferred free energy changes and required miscellaneous files should be downloaded from **[here]()** and unzipped in your project directory (see 'base_dir' option) i.e. where output files should be written.

# Installation Instructions

Make sure you have git and conda installed and then run (expected install time <5min):

```
# Install dependencies manually (preferably in a fresh conda environment)
conda install -c conda-forge bioconductor-biostrings cairo r-base>4.0.0 r-bio3d r-cairo r-data.table r-devtools r-ggally r-ggplot2 r-ggrepel r-hexbin r-plot3d r-ppcor r-reshape2 r-reticulate r-rocr r-roxygen2 r-scales

# Open an R session and install the archstabms R package
devtools::install_github("lehner-lab/archstabms")
```

# Usage

The top-level function **archstabms()** is the recommended entry point to the pipeline and by default reproduces the figures and results from the computational analyses described in the following publication: [The genetic architecture of protein stability (Faure AJ et al., 2023)](). See [Required Data](#required-data) for instructions on how to obtain all required data and miscellaneous files before running the pipeline. Expected run time <20min.

```
library(archstabms)
archstabms(base_dir = "MY_PROJECT_DIRECTORY")
```
