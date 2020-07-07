# Authors
Christopher L. Crawford^[*1], Lyndon D. Estes^[2], Timothy D. Searchinger^[1], and David S. Wilcove^[1, 3]

*Corresponding Author, @chriscra, ccrawford@princeton.edu, Robertson Hall, Princeton University, Princeton, NJ  
1. Princeton School of Public and International Affairs, Princeton University, Princeton, NJ  
2. Graduate School of Geography, Clark University, Worcester, MA  
3. Department of Ecology & Evolutionary Biology, Princeton University, Princeton, NJ  

# Description

This repository houses scripts with data and analyses for our paper "Consequences of under-explored variation in biodiversity indices used for land-use prioritization." We explore how different methods for designing biodiversity indices affect the outcomes of a land-use prioritization model designed to identify areas for agricultural expansion so as to meet a given production target while minimizing the environmental effects of that expansion, applying the [`agroEcoTradeoff`](https://github.com/PrincetonUniversity/agroEcoTradeoff) model to a case study in Zambia. The `agroEcoTradeoff` model allows users to minimize four constraints -- 1) biodiversity loss, 2) total agricultural area (maximizing yields), 3) carbon loss, and 4) transportation costs. Our primary analysis focuses on how biodiversity loss is modelled: we assess congruence between the least biodiverse areas in Zambia as identified using different ways of representing biodiversity.

## Downloading agroEcoTradeoff

Note that in order for these scripts to run, one must download the `agroEcoTradeoff` package, which must be downloaded directly. See: https://github.com/PrincetonUniversity/agroEcoTradeoff. In order to do this in terminal, first you have to install wget (and before that, gdal): https://stackoverflow.com/questions/33886917/how-to-install-wget-in-macos. You can use brew to do this, and then install wget. After that, run the following three lines of code from Lyndon's github, which installed the `agroEcoTradeoff` package:

wget https://github.com/PrincetonUniversity/agroEcoTradeoff/raw/master/installer.sh  
chmod +x installer.sh  
./installer.sh  

The first one uses wget to download the installer.sh script from Lyndon's github. (wget does a similar thing to curl, and I'm not sure why one is preferred over another.)
The second line changes the scripts permissions to allow it to be executed. "chmod" is a command for modifying the file's permissions. "+x" sets execute permissions. The script name is placed last, so that terminal knows which file to modify.
The third line runs the script itself.

In order for the `agroEcoTradeoff` model to work correctly, the working directory of your R project *must* be set to your agroEcoTradeoff/ directory. If you're working in a git repository, you can have it be a different name from the folder within which your R project lives.  You can create a scripts folder within your agroEcoTradeoff/ directory, and go from there. The data that the agroEcoTradeoff model pulls from is housed in a folder within external/data/folder_name. You use the name of this folder as the "input key" that points the model towards the folder you want when you run `tradeoff_mod(input_key = "folder_name")` function.

## Analyses

## Data

