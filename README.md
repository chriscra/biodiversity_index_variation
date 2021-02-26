# Consequences of under-explored variation in biodiversity indices used for land-use prioritization

[![DOI](https://zenodo.org/badge/277938120.svg)](https://zenodo.org/badge/latestdoi/277938120)

This repository houses scripts with data and analyses for:
> Crawford C.L.\*, Estes L.D., Searchinger T.D., and Wilcove, D.S. 2021. Consequences of under-explored variation in biodiversity indices used for land-use prioritization. *Ecological Applications.*

\*@chriscra, ccrawford@princeton.edu, Robertson Hall, Princeton University, Princeton, NJ

We explore how variation in the design of biodiversity indices affects the outcome of land-use prioritization. We use the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) land-use prioritization model (Estes and Spiegel 2016), a trade-off model designed to identify areas for agricultural expansion that meet a given production target at the least environmental cost, and apply the  model to a case study in Zambia. The [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) model allows users to minimize four constraints -- 1) biodiversity loss, 2) total agricultural area (maximizing yields), 3) carbon loss, and 4) transportation costs. Our analysis focuses on how biodiversity loss is modeled: specifically, we assess agreement between the least biodiverse areas in Zambia as identified by biodiversity indices that vary in their construction. We explore results for a wide range of criteria and methods that biologists and land-use planners have used, including: published composite indices, vertebrate taxonomic groups, metrics of species richness, methods for combining layers, and spatial resolutions.

This is a living repository. See the public release and Zenodo archive for the code used at time of publication. 

## Components of this Repository

This repository has one primary directory housing the scripts used in our analysis. Because the data sources are all either publicly available or available upon request from the appropriate organizations, this repository does not include the raw data itself.

The five primary `.Rmd` scripts are numbered based on the order in which parts of the analysis were conducted. Note that these `.Rmd` scripts are intended to be run interactively chunk by chunk (sometimes line by line), not knit together all at once. They also include code for data visualization and exploratory data analysis; not all code is strictly required. These files include:

- **1_Start.Rmd** serves as a useful working starting place, loading required packages, custom functions, and useful GIS components. This script was originally meant for development, and is not strictly necessary for reproducing the analysis.
- **2_BD_input_prep.Rmd** outlines the production of various biodiversity inputs explored in this project directly from data sources.
- **3_Model_Runs.Rmd** takes these biodiversity inputs and runs them through the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) model, applying seven primary model weighting specifications (distributing weight between biodiversity and other constraints like yield, carbon loss, and travel cost). This also includes model runs for a range of robustness checks.
- **4_Analyses.Rmd** analyzes the model outputs, most notably extracting and merging model conversion recommendations, calculating weighted Jaccard Similarity values, and producing a range of summary statistics. This script also contains code to produce a wide range of figures, including those in our main text and a larger number for the Appendix S1.
- **5_SI.Rmd** is a stand-alone document to produce our supplementary material, "Appendix S1."

Scripts starting with "cc_" serve as utility scripts, housing the packages (**cc_libraries.R**), custom functions (**cc_functions.R**), and file path names (**cc_pathnames.R**) used in the analysis.

## Downloading agroEcoTradeoff

Note that in order for these scripts to run, one must download the `agroEcoTradeoff` package, which must be downloaded directly. See: https://github.com/PrincetonUniversity/agroEcoTradeoff. To do accomplish this in a Unix shell environment (e.g. bash), make sure gdal and wget are installed. One can use brew to install gdal, and then install wget. ([See here for an example of how to install wget for Mac OS](https://stackoverflow.com/questions/33886917/how-to-install-wget-in-macos).) 

After that, run the following three lines of code from the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) GitHub in your shell environment, which will then install the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) package:

```
wget https://github.com/PrincetonUniversity/agroEcoTradeoff/raw/master/installer.sh
chmod +x installer.sh
./installer.sh
```

The first line uses wget to download the installer.sh script from the  [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) github.
The second line changes the scripts permissions to allow it to be executed. "chmod" is a command for modifying the file's permissions. "+x" sets execute permissions. The script name is placed last, so that the shell environment knows which file to modify.
The third line runs the script itself.

In order for the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) model to work correctly, the working directory of your R project *must* be set to your "agroEcoTradeoff/" directory. If you're working in a git repository, your repository can have a different name from the folder within which your R project lives, but the folder on your machine must be named "agroEcoTradeoff/".  You can create a scripts folder within your "agroEcoTradeoff/" directory, and go from there. The data that the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) model pulls from is housed in a folder within external/data/folder_name. Users provide the "folder_name" as the "input key" that points the model towards the folder you want when you run `tradeoff_mod(input_key = "folder_name")` function.
