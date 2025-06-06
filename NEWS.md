# MMVmalaria 1.4.0

MMVmalaria 1.4.0 updates to depend on IQRtools 2025.05

# Documentation
* Suggestion to use CRAN snapshot https://packagemanager.posit.co/cran/2025-04-22 added to README

# Dependencies
* IQRtools (>= 2025.05 from >= 2024.09 )
* MMVbase (>=1.2.1 from ==1.2.0)
* Imports added: survival 

# Build.sh 
Build.sh script added to build the R package

# .gitignore
Updated to be more functional and in line with MMVbase, but allowing .Rdata files 

# Functions
Confusion Matrix functions used to evaluate true/false positives/negatives added to MMVmalaria
* ConfusionMatrix() added to 09-KeyParameters.R 
* SummarizeConfusionMatrix() added to 09-KeyParameters.R

# Vignettes
PredictHumanDoses.Rmd
* IQRtools::IQRmodel(modelFile) explicitly called in vignette (previously IQRmodel(modelFile))