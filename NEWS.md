# MMVmalaria 1.6.0

# New functions
* `getExcursionsAboveMIC()` added: detects each contiguous excursion (time interval) during which concentration stays continuously above the MIC, one row per excursion with its own onset/end time - needed to detect and visualize gaps (troughs dropping back below MIC between doses) for multi-dose regimens, where `getTimeAboveMIC()`'s summed total can hide a real interruption in coverage.
* `evaluateDoseCriterion_ContinuousTimeAboveMIC()` added: an alternative to `evaluateDoseCriterion_TimeAboveMIC()` for `predictDose_Generic()`, for multi-dose regimens where continuous coverage (no gap) matters, not just the summed total time above MIC. Targets the length of the *last* uninterrupted excursion above MIC (the one anchored at the profile's global peak and its decline after the final dose) - not the first, which would converge every target duration onto the same "merge every early gap" dose regardless of how long a coverage window was actually requested.

# Tests
* `tests/testthat/test-keyParameters.R` added, covering `getExcursionsAboveMIC()` and `evaluateDoseCriterion_ContinuousTimeAboveMIC()`, including an asymmetric-excursion-length case that specifically distinguishes first-vs-last excursion targeting.



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
