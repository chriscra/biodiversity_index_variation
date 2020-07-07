### Libraries used for Zambia Biodiversity project
pkgs <- c("raster", "rgdal", "sf", "sp", "gdalUtils", "rasterVis", "data.table",
          "tictoc", "fasterize", "viridis", "lwgeom", "magrittr",
          "parallel", "reshape2", "cowplot", "tidyverse",
          "knitr", "kableExtra", "spatstat", "rmarkdown","bookdown", "pryr") # main package names

pkgs_aux <- c("citr", "rticles", "mapview", "gdata", "GISTools",
          "rgeos", "RColorBrewer", "stringr",
          "maps", "rangeBuilder", "rnaturalearth", "spatialEco", "smoothr",
          "ggnewscale", "plotly", "scales", "patchwork", "biomod2") # aux package names
#install.packages(pkgs)
#install.packages(pkgs_aux)
#update.packages()

inst <- lapply(pkgs, library, character.only = TRUE) # load them
inst_aux <- lapply(pkgs_aux, library, character.only = TRUE) # load them

lyndon_pkgs <- c("dtraster", "lmisc", "agroEcoTradeoff")
lyndon_inst <- lapply(lyndon_pkgs, library, character.only = TRUE) # load them
#install.packages(lyndon_pkgs)

# generate .bib document from packages I've loaded and used,
# to then add to Mendeley to make it into the library.bib doc
# that is source file for the notebook's citations.
# getwd()
write_bib(
  c(lyndon_pkgs, pkgs, "base"),
  file = "Zambia_packages_2020.bib")

# note that some of the packages will say something like this:

##There are binary versions available but the source versions are later:
##  binary source needs_compilation
##dplyr  0.8.2  0.8.3              TRUE
##haven  2.1.0  2.1.1              TRUE
##sf     0.7-4  0.7-5              TRUE
##Do you want to install from sources the packages which need compilation? (Yes/no/cancel)

# usually this isn't a good idea to do.

# In order to get the packages to install and compile, I downloaded the latest version of R, and at their advice, downloaded clang, gfortran, and XQuartz.

# to install Lyndon's packages, and the agroEcoTradeoff model.
# check out https://github.com/PrincetonUniversity/agroEcoTradeoff
# run the following three lines in the terminal, after setting the directory to
# wherever you want it to save the installer and the agroEcoTradeoff folder.



# you then install the agroEcoTradeoff package with the following lines:
# (just remove the comment #s)
########install.packages("devtools", repos = 1:2)
#install.packages("devtools")
#library(devtools)
#install_github("ldemaz/dtraster")
#install_github("PrincetonUniversity/lmisc")
#install_github("PrincetonUniversity/agroEcoTradeoff@devel")

