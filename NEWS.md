# MMVmalaria 1.5.0

# Pull request template
* Added

# Tests
* Unit test added test-modelLibraries to test that models within defined model libraries (inst/modelLibrary and modelLibrary_tobs0) run without error. 

# Model library
* model_PKPD_3cpt_ClearanceModel and model_PKPD_3cpt_abs0_ClearanceModel added to model library. These are the simple EMAX clearance model (for SCID modelling) with the Tobs0 specification. 
* Pass made on Tobs0 models to remove text chunks describing history of PLerr and Tobs0 models
* Baseline section of Tobs0 models corrected where they incorrectly referred to the PLerr methodology
* Clearance models description of "zombie state" updated to better reflect the underlying biology
* Formatting adjustment made to all library models to improve readability and standardize model notes section
* Alpha renamed to Alpha1 in all appropriate models (GDPI models) to avoid conflict with IQRtools reserving "alpha" string. 
* _tobs0 removed from names of model file .txt in modelLibrary_tobs0 as the library name makes this clear without needing _tobs0 in the file name.
