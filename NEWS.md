# MMVmalaria 1.3.0

MMVmalaria 1.3.0 updates to depend on IQRtools 2025.05

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
All
* Extensive documentation pass to ensure all arguments are well captured by Roxygen2 standards

Confusion Matrix functions used to evaluate true/false positives/negatives added to MMVmalaria
* ConfusionMatrix() added to 09-KeyParameters.R 
* SummarizeConfusionMatrix() added to 09-KeyParameters.R
MMVbase:: prefixed to all functions exported from MMVbase that were previously in MMVmalaria 
predictDose_Generic
* Adjusted to import and declare IQRtools:: for IQRmodel() 

03-DataPreparation.R 
* import_SCIDpkpdDataCombo_TAD : Argument PDrangeDoseTimeRange removed as unused
* load_AfricanPediatricMalariaPopulation2to5 REMOVED 

19-IQRaux.R
* verifyArg() Removed, now in IQRtools
* createSimIQRmodelSpec() Removed, now in IQRtools
* split_eventTable() Removed, now in IQRtools
* getDefaultEventTable() Removed, now in IQRtools
* addTminToEventTable() Removed, now in IQRtools 
* sim_IQRmodel_consistent() Removed, now in IQRtools as sim_IQRmodel_strict() 

98-OldFunctions.R 
* REMOVED

# Vignettes
PredictHumanDoses.Rmd
* IQRtools::IQRmodel(modelFile) explicitly called in vignette (previously IQRmodel(modelFile))
* Calls to MalariaPopulation.Rdata adjusted 
* Definition of WT0 from WEIGHT_kg added 

# Tests
test-predictHumanDose.R
* Calls to MalariaPopulation.Rdata adjusted 

# Examples
findMinMMV_ex
* Adjusted findMinMMV to MMVbase::find_MinMMV to reflect function being moved to MMVbase and renamed 