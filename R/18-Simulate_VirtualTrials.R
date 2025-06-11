#' Simulate virtual clinical trials for a set of experimental setting and dose combinations
#'
#' @md
#'
#' @details For shortness, the word 'scenario' is used to indicate a combination of an experimental setting and a dose.
#'  Having two types of covariates, namely, individual and experimental, allows to use the same
#'  individuals for different experimental conditions, and, therefore, to see the difference between
#'  experimental conditions with more clarity.
#'
#' @param projectPathList List containing the path of all models to use for sampling. Each element of the list should point to an 'IQRnlmeProject', 'IQRsysProject' and/or a 'GPF' file (see [IQRtools::GPF]).
#' @param modelFile Path to the PKPD model equation file to be used for the simulations.
#' @param outputFolder a character string denoting a path to a directory where a .rds file containing the resulting list summarizing each simulated trial is to be saved.
#' @param regressorExtra a character vector containing the names of additional regressors to be used in the simulations. Default: NULL.
#' @param Doses Numeric vector of the doses to tes,t or a dosing table containing columns "DoseID", "TIME", "ADM",
#' "AMT", "TINF", such that multiple dosing regimen correspond to different DoseID values. Regimens with time-varying dose
#' are supported unless, the dose is a covariate of the model. Note that the TINF column is neglected and overriden by
#' 1e-4. If there is a need to specify infusion time for a given INPUT, this can be done by specifying the Tk0 argument.
#' @param T_FirstDose Numeric defining time of first dose. Might be diferent than initial time of simulation (Default: 0).
#' @param nbrDoses Numeric Vector. If length is 1, all doses in 'Doses' will have the same number of doses, otherwise should be of length 'Doses'.
#' @param timeInterval Numeric defining the interval between doses (Note: Sssumed to be the same for all doses)
#' @param DOSEcovariate a character vector containing the names of the covariates to be used as dose covariate. Default: `c("INPUT1"="DOSELEVEL")`.
#' @param N_Trial
#' @param N_SubjPerTrial
#' @param IndCovSample Dataframe or Function to generate individual covariate (e.g. WT0, BMI, SEX, etc). Default: NULL.
#' @param ExpCovSample Dataframe of experimental covariates to account for (e.g. Formulation, Fasted/Fed, etc). Default: NULL.
#' @param FLAG_SAMPLE
#' @param replace
#' @param setting a named list with elements among the following:
#' * FLAGreturnSimPKPD A logical (default: FALSE) indicating if simPKPD should be included in the resulting list.
#' Note that setting this argument to TRUE is memory intensive and can cause memory crashes.
#' * FLAGsaveSimPKPDToFile a logical (default FALSE) indicating if the simPKPD argument should be saved in
#' `outputFolder` in the result .rds file for the corresponding trial. Setting this to TRUE while outputFolder is NULL
#' will raise in an error.
#' * FLAGreuseStoredTrialSummaryFiles a logical (default FALSE) indicating if previously stored .rds trial summary files in `outputFolder` should
#' be reused instead of rerunning each trial simulation.
#' * fun_SummarizeTrial a character string denoting the name of a summarizing function to be
#' called on the simulation results for an experimental setting, a dose, and a virtual trial population. Default:
#' `"MMVmalaria::summarizeTrial_Basic"`. This function should accept at least a data.frame
#' and a character string argument, namely a data.frame containing the PKPD model simulation result and individual model parameters as
#' returned from the function [IQRtools::sim_IQRmodel()]) and a filename containing the path to the trial .rds file on the disk.
#' The .rds file can be loaded into a list and used to get the individual covariate values for each subject in the trial. Example code:
#' `trial <- readRDS(trialFileName); simPKPDcovs <- merge(simPKPD, trial$covariates, by = intersect(names(simPKPD), names(trial$covariates)))`.
#' Any additional arguments of the function can be specified in the list argument
#' `args_SummarizeTrial`. Note that these arguments will have the same values over all trials. The function is expected to return
#' a named list with at least three data.frame members: summaryPKPD.ByTrial, summaryClinEnd.ByIndiv, summaryClinEnd.ByTrial. In addition,
#' the function can store result objects into .rds files for post-processing.
#' See [summarizeTrial_Basic()] for an example.
#' * args_SummarizeTrial a named list or NULL (default) containing additional arguments
#' to be passed to `fun_SummarizeTrial`. See documentation of [summarizeTrial_Basic()]
#' for possibilities.
#' * .paropts a named list specifying arguments passed to [foreach::foreach()]. By default only this member is set to
#' `.paropts = list(.inorder = TRUE,
#'                       .packages = c ("IQRtools","MMVmalaria"))`.
#' @param Tk0 a named character vector with names INPUT1, INPUT2, etc, containing the names of corresponding Tk0 parameters.
#' Note that the value of the Tk0 parameter for an input will overwrite the value of TINF in the dosing object.
#' Default `c("INPUT1"="Tk0")`. INPUTs in this argument are assumed to have 0-order absorption, i.e. Tk0 is
#' the time for the dose to be fully absorbed at a constant rate. This is suitable to simulate infusion at constant rate.
#' The Tk0 parameters should be either present in the GPF file or available as regressors in the regressorExtra argument.
#' @param Fpediatric
#' @param Nparallel
#' @param simMode To chose which mode to use for simulations between "Clinical" (New set of individuals sampled for each scenario - Default) and
#' "Compare" (Same set of individuals used between each scenario adjusted by `ExpCovSample` and Dose-regimen as defined in `Doses`.).
#' @param FLAGrunSimulation Run simulations after PK and PD parameters sampling
#'
#' @return a named list with the following members:
#'
#' * summaryPKPD.ByTrial a data.frame resulting from rbinding all summaryPKPD.ByTrial data.frames returned by
#' `setting$fun_SummarizeTrial`.
#' * summaryClinEnd.ByIndiv a data.frame resulting from rbinding all summaryPKPD.ByIndiv data.frames returned by
#' `setting$fun_SummarizeTrial`.
#' * summaryClinEnd.ByTrial a data.frame resulting from rbinding all summaryClinEnd.ByTrial data.frames returned by
#' `setting$fun_SummarizeTrial`.
#' * Covariates a data.frame
#' * PKPDparameters a data.frame
#' * eventTable a data.frame
#' * eventData a data.frame
#' * regression
#' * abs0inputs
#' * abs0Tk0param
#' * setting
#'
#'
#' @export
#' @family Simulations
#' @author Mohammed H. Cherkaoui, \email{cherkaouim@@mmv.org}, Venelin Mitov \email{venelin.mitov@@intiquan.com}
simulate_VirtualTrials <- function(
    projectPathList,
    modelFile      ,
    outputFolder,
    regressorExtra                                = NULL,
    Doses,
    T_FirstDose                                   = 0,
    nbrDoses                                      = 1,
    timeInterval                                  = 24,
    DOSEcovariate                                 = c("INPUT1"="DOSELEVEL"),  # Indicate name for the dose covariate
    N_Trial                                       = 100,
    N_SubjPerTrial                                = 10,
    IndCovSample                                  = NULL,
    ExpCovSample                                  = NULL,
    FLAG_SAMPLE                                   = 1,
    replace                                       = TRUE,   # Setting for the function sample().
    setting                                       = list(.paropts = list(.inorder = TRUE,
                                                                         .packages = c ("IQRtools","MMVmalaria"))),
    Tk0                                           = c("INPUT1"="Tk0"),
    Fpediatric                                    = NULL,
    Nparallel                                     = 1,
    simMode                                       = c("Clinical","Compare")[1],
    FLAGrunSimulation                             = TRUE) {

  # Simulations of virtual trials:
  # Two types of covariates are accounted for:
  #   - IndCovSample:   Covariates that corresponds to individual covariates such as WT0 or AGE
  #   - ExpCovSample:   Covariates that corresponds to experimental conditions.


  #---------------------------------------------------#
  # STEP 0: Run some basic check ----
  #---------------------------------------------------#

  #-------------------------------------------#
  ## 0.a: Validate setting  ----
  #-------------------------------------------#

  setting <- validateSetting_SimulateVirtualTrials(setting)


  #-------------------------------------------#
  ## 0.b: Create outputFolder ----
  #-------------------------------------------#
  if(is.null(outputFolder) || identical(outputFolder, "") || !is.character(outputFolder)) {
    stop("simulate_VirtualTrials: outputFolder should be a valid directory path.")
  } else if(!dir.exists(outputFolder)) {
    dir.create(file.path(outputFolder))
  } else if(dir.exists(outputFolder) && !setting$FLAGreuseStoredTrialSummaryFiles) {
    rmdirStatus <- try(unlink(outputFolder, recursive = TRUE, force = TRUE), silent = TRUE)
    if(inherits(rmdirStatus, "try-error") || rmdirStatus == 1) {
      stop("simulate_VirtualTrials: FLAGreuseStoredTrialSummaryFiles is set to FALSE, but failing to delete existing outputFolder ", outputFolder)
    } else {
      dir.create(file.path(outputFolder))
    }
  }

  #-------------------------------------------#
  ## 0.c: Save setting to output folder ----
  #-------------------------------------------#

  saveRDS(setting, file = file.path(outputFolder, "setting.rds"))

  #---------------------------------------------------#
  # STEP 1: Preparations ----
  #---------------------------------------------------#

  # Message:
  start_time <- Sys.time()
  cat("\n\n----------Running Simulate_VirtualTrials()----------\n")
  cat("\tSTEP 1: Checking all variables\n")


  # Check DOSEcovariate:
  DOSEcovariate.Names   <- NULL
  DOSEcovariate.Numeric <- NULL
  if(!is.null(DOSEcovariate)){
    # Check if character:
    if(!is.character(DOSEcovariate)){
      DOSEcovariate <- c("INPUT1"="DOSELEVEL")
      DOSEcovariate.Names   <- "INPUT1"
      DOSEcovariate.Numeric <- 1
      warning("Variable 'DOSEcovariate' is not a character vector: It was set to c('INPUT1'='DOSELEVEL') - Adjust if needed")
    }

    # Check if it has names:
    if(!is.null(names(DOSEcovariate))){
      # As it has names, it needs to have valid names:
      DOSEcovariate.Names <- names(DOSEcovariate)
      DOSEcovariate.Numeric <- as.numeric(gsub("INPUT","",DOSEcovariate.Names))

      # Check that DOSEcovariate.Numeric is numeric
      if(any(is.na(DOSEcovariate.Numeric))){
        stop("'DOSEcovariate' does not have valid names: It should be defined as DOSEcovariate <- c(INPUT1='DOSELEVELx1', INPUT2='DOSELEVELx2)")
      }

    }else{
      names(DOSEcovariate)  <- paste0("INPUT", 1:length(DOSEcovariate))
      DOSEcovariate.Names   <- names(DOSEcovariate)
      DOSEcovariate.Numeric <- 1:length(DOSEcovariate)
    }
  }

  # Check Tk0:
  if (!is.null(Tk0)){
    # Check if character:
    if(!is.character(Tk0)){
      Tk0 <- NULL
      warning("Variable 'Tk0' is not a character vector: It was set to NULL - Adjust if needed")
    }

    # Check if it has names:
    if(!is.null(names(Tk0))){
      # As it has names, it needs to have valid names:
      Tk0.Names <- names(Tk0)
      Tk0.Numeric <- as.numeric(gsub("INPUT","",Tk0.Names))

      # Check that Tk0.Numeric is numeric
      if(any(is.na(Tk0.Numeric))){
        stop("'Tk0' does not have valid names: It should be defined as Tk0 <- c(INPUT1='Tk0x1', INPUT2='Tk0x2)")
      }
      rm(Tk0.Names, Tk0.Numeric)

    }else{
      names(Tk0) <- paste0("INPUT", 1:length(Tk0))
    }
  }

  # If 'Doses' is a numeric vector:
  if (is.numeric(Doses)){
    # Number of dose groups to simulate:
    nbrDoseRegimens <- length(Doses)

    # Check if Doses & nbrDoses have same length:
    #   Only if length(nbrDoses)>1
    if (length(nbrDoses)>1 & length(nbrDoses)!=nbrDoseRegimens){
      stop("'nbreDoses' should be of length 1 or same as 'Doses'.")
    } else if(length(nbrDoses)==1){
      nbrDoses <- rep(nbrDoses,nbrDoseRegimens)
    }

    # Check if Doses & timeInterval have same length:
    #   Only if length(timeInterval)>1
    if (length(timeInterval)>1 & length(timeInterval)!=nbrDoseRegimens){
      stop("'timeInterval' should be of length 1 or same as 'Doses'.")
    } else if(length(timeInterval)==1){
      timeInterval <- rep(timeInterval,nbrDoseRegimens)
    }

    # if 'Doses' is a dosing object (i.e. a data.frame)
  } else if (is.data.frame(Doses)){
    if (!all(c("DoseID", "TIME", "ADM", "AMT", "TINF") %in% names(Doses))){
      stop("If 'Doses' is a data frame, it should contain the columns DoseID, TIME, ADM, AMT, TINF: Please adjust 'Doses'.")
    }else{
      # Check that TINF column is not specified to anything different than 1e-4, since such values are neglected. Tk0 argument
      # should be used instead.
      if(length(unique(Doses$TINF)) != 1 || unique(Doses$TINF) != 1e-4) {
        warning("The TINF column in the Doses argument will be overriden by 1e-4. Consider specifying the Tk0 argument if you need to simulate infusion or 0-order absorption time for some of the INPUTs.")
      }
      # number of dosing regimen
      nbrDoseRegimens <- length(unique(Doses$DoseID))
    }

    # Convert to data.table:
    data.table::setDT(Doses, key = c("DoseID","ADM"))

    # If USUBJID is present in Doses, the number of individual should be equal to N_SubjPerTrial
    if("USUBJID" %in% names(Doses) || "IndivID" %in% names(Doses)){
      # Change name of USUBJID to match the IndivID of each simulated tiral
      data.table::setnames(Doses,"USUBJID","IndivID",skip_absent = TRUE)

      # Check if length is consitant
      if(length(unique(Doses$IndivID))!=N_SubjPerTrial){
        stop("'USUBJID' is present in 'Doses', but the number of individual is not equal to 'N_SubjPerTrial': Please adjust 'Doses'")
      }
      # else if(any(is.na(match(unique(Doses$IndivID), 1:N_SubjPerTrial)))){
      #   stop("'USUBJID' in 'Doses' does not match the simulated USUBJID: Please adjust 'Doses', 'USUBJID' shoudl be formed as rep(1:N_Trial, each = N_SubjPerTrial)*10^(nchar(N_SubjPerTrial))+rep(1:N_SubjPerTrial, times = N_Trial)")
      # }
    }
  }else{
    stop("'Doses' should be a numeric vector or a data frame: Not coded for dosing object yet.")
  }

  # check replace:
  if(nrow(IndCovSample)<N_SubjPerTrial && !replace){
    replace <- TRUE
    warning("'replace' was set to TRUE and nrow(IndCovSample)<N_SubjPerTrial")
  }

  # Number of simulation ID:
  N_Subj   <- N_Trial*N_SubjPerTrial
  maxID    <- N_Subj*nbrDoseRegimens
  if(!is.null(ExpCovSample)){
    maxID <- maxID*nrow(ExpCovSample)
  }
  maxIDrun <- ceiling(maxID/N_SubjPerTrial)

  # Check if 'FLAG_SAMPLE' is null:
  if(is.null(FLAG_SAMPLE)){
    stop("'Please define 'FLAG_SAMPLE': Should be 0, 1, 2, 3, 4 or 5")
  }

  # Check which mode to be used for the simulations:
  if(simMode=="Compare"){
    FLAGclinicalSimMode <- FALSE
    ScenIDsimMode       <- NULL
  }else{
    FLAGclinicalSimMode <- TRUE
    ScenIDsimMode       <- "ScenID"
  }

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #---------------------------------------------------#
  # STEP 2: Handling of models ----
  #---------------------------------------------------#

  # Message:
  start_time <- Sys.time()
  cat("\tSTEP 2: Checking model(s)\n")

  # Load equation of the PKPD model:
  modelIQR <- IQRmodel(modelFile)
  IQRtools::export_IQRmodel(modelIQR, filename = file.path(outputFolder, "model.txt"))

  # Combined all model into one GPF model:
  if(IQRtools::is_GPF(projectPathList)){
    fullModelGPF <- projectPathList
  }else{
    modelGPFlist <- lapply(projectPathList, function(x){
      IQRtools::load_GPF(x)
    })
    fullModelGPF <- IQRtools:::combineGPF(modelGPFlist)
  }
  IQRtools::export_GPF(fullModelGPF, filename = file.path(outputFolder, "GPF.xlsx"))

  # Regressor Names:
  idx_mpar         <- which(fullModelGPF$estimates$TYPE=="MODEL PARAMETER")
  regressorNameGPF <- fullModelGPF$estimates$PARAMETER[idx_mpar]

  # Add Extra Regressor:
  regressorName <- c(regressorNameGPF, regressorExtra)
  idx <- duplicated(regressorName)
  if (any(idx)){
    message(MMVbase::collapseMMV(x = regressorName[idx]), " regressor are already present in the model(s) defined in 'projectPathList': Remove them from 'regressorExtra'.")
    regressorName <- unique(regressorName)
  }

  # Keep only the parameters present in the model:
  if (!all(regressorName %in% names(modelIQR$parameters))){
    message("The parameters ", MMVbase::collapseMMV(x = setdiff(regressorName, names(modelIQR$parameters))), " are present in the model(s) defined in 'projectPathList' but not present in 'modelFile': They were removed.")
    regressorName <- intersect(regressorName, names(modelIQR$parameters))
  }

  # Get Ind. Covariate as dataFrame to check later:
  if (!is.null(IndCovSample)){
    if (is.function(IndCovSample)){
      IndCovSampleTest <- IndCovSample(2)
    } else{
      IndCovSampleTest <- IndCovSample[1:2, , drop = FALSE]
    }
    #   If USUBJID is present in IndCovSample, it is renamed to not mix with the USUBJID generated in simulate_VirtualTrials
    data.table::setnames(IndCovSampleTest,c("USUBJID"),c("USUBJID.IndCov"),skip_absent=TRUE)

  }else{
    IndCovSampleTest <- NULL
  }
  IndCovName  <- names(IndCovSampleTest)
  IndCovClass <- sapply(names(IndCovSampleTest), function(x){class(IndCovSampleTest[[x]])})

  # Check that Fpediatric is defined in  IndCovSample or ExpCovSample:
  if (!is.null(Fpediatric)){
    if (!xor(Fpediatric %in% names(IndCovSampleTest), Fpediatric %in% names(ExpCovSample))){
      stop("Fpediatric=", Fpediatric, " is not available in 'IndCovSample' nor 'ExpCovSample', or is in both: It should be defined in 'IndCovSample' or 'ExpCovSample' exclusively.")
    }
  }

  # Check PK & PD model covariates:
  covPKPDmodel <- IQRtools:::getNamesCovariates(fullModelGPF)

  # Covariate FLAGs:
  if (length(covPKPDmodel) == 0){
    FLAGcovPKPD <- FALSE
    FLAGcovDOSE <- FALSE
  }else if (length(covPKPDmodel) == 1 & covPKPDmodel[1] == "") {
    FLAGcovPKPD <- FALSE
    FLAGcovDOSE <- FALSE
  } else {
    FLAGcovPKPD <- TRUE
    if (!all(setdiff(covPKPDmodel,"") %in% c(names(IndCovSampleTest), names(ExpCovSample), DOSEcovariate))){
      # case a covariate not supplied by user
      stop(paste0("Please give values to sample from for the following covariates:\n", paste0(covPKPDmodel, collapse = ", ")))
    }

    # Make sure that DOSEcovariate is not in IndCovSample or ExpCovSample:
    if (any(DOSEcovariate %in% names(IndCovSampleTest))){
      stop(MMVbase::collapseMMV(DOSEcovariate[DOSEcovariate %in% names(IndCovSampleTest)]),
           " in 'DOSEcovariate' should not be included in 'IndCovSample'.")
    }
    if (any(DOSEcovariate %in% names(ExpCovSample))){
      stop(MMVbase::collapseMMV(DOSEcovariate[DOSEcovariate %in% names(ExpCovSample)]),
           " in 'DOSEcovariate' should not be included in 'ExpCovSample'.")
    }

    # Check if DOSEcovariate is a covariate of the model:
    if (DOSEcovariate %in% covPKPDmodel){
      FLAGcovDOSE <- TRUE
    } else{
      FLAGcovDOSE <- FALSE
    }
  }

  # If dose is used as covariate check that dose=0 is not simulated,
  # but set to 1e-12:
  if (FLAGcovDOSE) {
    if (is.numeric(Doses)){
      idx0 <- (Doses==0)
      if (any(idx0)) {
        message("Dose is a covariate in the model. Dose of 0 cannot be simulated as continuous covariates are log-transfomred.\nSimulation with a dose of 1e-12 value instead.")
        Doses[idx0] <- 1e-12
      }
    }else if (is.data.frame(Doses)){
      idx0 <- (Doses$AMT==0)
      if (any(idx0)) {
        message("Dose is a covariate in the model. Dose of 0 cannot be simulated as continuous covariates are log-transfomred.\nSimulation with a dose of 1e-12 value instead.")
        Doses$AMT[idx0] <- 1e-12
      }
    }
  }

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #---------------------------------------------------#
  # STEP 3: Create Dosing Object ----
  #---------------------------------------------------#

  # Message:
  start_time <- Sys.time()
  cat("\tSTEP 3: Create Dosing Object\n")

  # ~~~~~~~~~~~~~~~#
  # Create Dosing Object:
  #   If 'Doses' is a numeric vector
  colIndivID <- NULL
  if(is.numeric(Doses)){
    dosing <- data.table(DoseID = 1:nbrDoseRegimens)
    dosing <- dosing[,{
      DoseID_k <- DoseID[1]
      list(TIME = seq(T_FirstDose, by=timeInterval[[DoseID_k]], length.out=nbrDoses[[DoseID_k]]),
           ADM  = 1,
           AMT  = Doses[[DoseID_k]],
           TINF = 0)},
      by = "DoseID"]

    #   If 'Doses' is a data frame
  }else if(is.data.frame(Doses)){
    if("IndivID" %in% names(Doses)){
      colIndivID <- "IndivID"
    }
    dosing <- data.table:::subset.data.table(Doses, select = c("DoseID", colIndivID, "TIME", "ADM", "AMT", "TINF"))
  }
  #   Dose Group for output
  Doses.Groupe <- dosing[ADM==1,list(Dose=AMT[1], nbrDoses=length(AMT)), by = c("DoseID", colIndivID)]
  #   Set dosing key
  setkey(dosing,"DoseID", "TIME")

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #---------------------------------------------------#
  # STEP 4: Parameter Sampling ----
  #---------------------------------------------------#
  #-------------------------------------------#
  ## 4.a: Manage Experimental Covariate  ----
  #-------------------------------------------#
  # NOTE: Dose and ExpCovSample are covariate that we want to control
  #       within a trial to differentiate cohorts.
  # At the end of 4.a, we have a data.frame called ExpCovSampleExpanded with one row per combination of a
  # tupplet of experimental covariates and a dose-regimen, i.e. one row per scenario.
  # TODO Implement a check that fails if dose is a covariate and the amount is varying during a dose regimen.

  # Message:
  start_time <- Sys.time()
  cat("\tSTEP 4: Covariate and Parameter Sampling\n")
  cat("\t   4.a: Manage Experimental Covariate if 'ExpCovSample' or Dose is a covariate\n")

  # Add dose as a covariate:
  #   Get Dose as a data.table
  if(is.numeric(Doses) & length(DOSEcovariate)<=1){
    Doses.Cov <- data.table(DoseID = 1:nbrDoseRegimens,
                            AMT    = Doses,
                            key    = c("DoseID"))
  }else if(is.numeric(Doses) & length(DOSEcovariate)>1){
    stop("'DOSEcovariate' has 2 or more element, meaning there are multiple site of administration. In this case, 'Doses' can not be a numerical vector and should be a dosing object.")
  }else if(is.data.frame(Doses)){
    # Get the dose-level of each dosing regimen as a covariate, assuming that each regimen has a constant dose.
    # (not implemented time-varying doses).

    # Keep on the first dose:
    #   This means that if the doses changes with time it will not be reflect in the covariate
    Doses.Cov <- Doses[,{
      idx_Tmin <- which.min(TIME)
      list(AMT = AMT[idx_Tmin])
    },keyby=c("DoseID","ADM")]
    # Keep on the ADM of interest
    Doses.Cov <- Doses.Cov[ADM %in% DOSEcovariate.Numeric,c("DoseID","ADM","AMT")]
    Doses.Cov <- data.table::dcast(Doses.Cov,DoseID~ADM, value.var="AMT")
    if(ncol(Doses.Cov)-1<length(DOSEcovariate)){
      idx_Missing <- !(DOSEcovariate.Numeric %in% names(Doses.Cov))
      Doses.Cov[,as.character(DOSEcovariate.Numeric[idx_Missing])] <- 0
      data.table::setcolorder(Doses.Cov, c("DoseID", as.character(DOSEcovariate.Numeric)))
    }
  }
  names(Doses.Cov) <- c("DoseID", DOSEcovariate)
  #   Convert ExpCovSample to  data.table and add ExpID
  if (!is.null(ExpCovSample)){
    ExpCovSample <- data.table::setDT(ExpCovSample)
    ExpCovSample[,
                 ExpID := (1:.N)]
    setkey(ExpCovSample, "ExpID")
  }else{
    ExpCovSample <- data.table(ExpID = 1,
                               key   = c("ExpID"))
  }
  ExpCovSampleExpanded <- ExpCovSample[ , c(.SD, Doses.Cov), by = "ExpID"]

  # Get number of scenario:
  N_Scenario <- nrow(unique(ExpCovSampleExpanded[,c("ExpID","DoseID")]))

  # Add ScenID:
  ExpCovSampleExpanded[,
                       ScenID := max(DoseID)*(ExpID-1) + DoseID]
  #   Get vector as it might not be 1 to N_Scenario
  ScenID_Vec <- unique(ExpCovSampleExpanded$ScenID)

  # ExpID/USUBJID:
  #   If USUBJID is present in ExpCovSample, it is renamed to not mix with the USUBJID generated in simulate_VirtualTrials
  data.table::setnames(ExpCovSampleExpanded,c("USUBJID"),c("USUBJID.ExpCov"),skip_absent=TRUE)

  # Set key:
  setkey(ExpCovSampleExpanded, "ScenID", "ExpID", "DoseID")

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #-------------------------------------------#
  ## 4.b: Sample Individual Covariate for each trial  ----
  #-------------------------------------------#

  # NOTE: IndivCovSample are the covariate that we have not control within a scenario

  # At the end of this step we have a data.frame IndCovSampled, with one row per individual in a scenario of a trial
  # In Clinical mode the number of rows is (N_Trial*N_SubjPerTrial*N_Scenario)
  # In compare mode the number of rows is (N_Trial*N_SubjPerTrial).
  #
  # Message:
  start_time <- Sys.time()
  cat("\t   4.b: Sample individual covariates\n")

  # Assign memory for Ind. Covariate:
  #   If Clinical simulation mode
  if(FLAGclinicalSimMode){
    IndCovariate <- data.table::data.table(ScenID  = rep(ScenID_Vec      , each  = N_Trial*N_SubjPerTrial),
                                           TrialID = rep(rep(1:N_Trial   , each  = N_SubjPerTrial),times=N_Scenario),
                                           IndivID = rep(1:N_SubjPerTrial, times = N_Scenario*N_Trial))
    IndCovariate <- IndCovariate[,USUBJID := ScenID*10^(nchar(N_Trial)+nchar(N_SubjPerTrial))+TrialID*10^(nchar(N_SubjPerTrial))+IndivID]
    data.table::setkey(IndCovariate, "ScenID", "TrialID", "IndivID", "USUBJID")

    #   If Compare simulation mode
  }else{
    IndCovariate <- data.table::data.table(TrialID = rep(1:N_Trial       , each  = N_SubjPerTrial),
                                           IndivID = rep(1:N_SubjPerTrial, times = N_Trial))
    IndCovariate <- IndCovariate[,USUBJID := TrialID*10^(nchar(N_SubjPerTrial))+IndivID]
    data.table::setkey(IndCovariate, "TrialID", "IndivID", "USUBJID")
  }
  #   Add columns and adjust Classes:
  if (!is.null(IndCovName)){
    IndCovariate <- IndCovariate[,(IndCovName):=lapply(IndCovName, function(col_k){
      if (IndCovClass[[col_k]]=="numeric"){
        out <- 0
      }else if(IndCovClass[[col_k]]=="integer"){
        out <- 0L
      }else{
        out <- "a"
      }
      out
    })]

    # Adjust Class Type Based on Covariates:
    for (IndCovName_k in IndCovName){
      IndCovariate[,(IndCovName_k) := eval(parse(text=paste0("as.",IndCovClass[[IndCovName_k]],"(", IndCovName_k,")")))]
    }
  }

  # Sample around Trial:
  # TODO Consider disabling IndCovSample being a function
  # TODO Consider rewriting without pre-reservation of memory.
  if (!is.null(IndCovName)){
    IndCovariate[,(IndCovName):={
      # Sample:
      if (is.data.frame(IndCovSample)){
        if(nrow(IndCovSample)==N_SubjPerTrial){
          IndCovSampled <- IndCovSample
        }else{
          IndCovSampled <- as.data.frame(IndCovSample[sample(nrow(IndCovSample), N_SubjPerTrial, replace=replace), ])
        }
        names(IndCovSampled) <- IndCovName
      } else if (is.function(IndCovSample)){
        IndCovSampled <- IndCovSample(N_SubjPerTrial)
      } else if (is.null(IndCovSample)){
        IndCovSampled <- NULL
      } else{
        stop("'IndCovSample' should be NULL, a data frame or a function.")
      }

      # f USUBJID is present in IndCovSampled, it is renamed to not mix with the USUBJID generated in simulate_VirtualTrials
      data.table::setnames(IndCovSampled,c("USUBJID"),c("USUBJID.IndCov"),skip_absent=TRUE)

      # data.table output:
      lapply(IndCovName, function(col_k){IndCovSampled[[col_k]]})
    },
    by=c(ScenIDsimMode,"TrialID")]
  }

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #-------------------------------------------#
  ## 4.c: Concatenate Exp. and Ind. covariates ----
  #-------------------------------------------#
  # At the end of this step we have a data.frame Covariates with
  # Number of experiments * Number of dosing regimens * nrow(IndivCovSampled).
  # Message:
  start_time <- Sys.time()
  cat("\t   4.c: Concatenate Exp. and Ind. covariates\n")

  # Generate Covariate for output:
  # TODO Logically the same merge should be used for both modes.
  if(FLAGclinicalSimMode){
    Covariates <- merge(ExpCovSampleExpanded,IndCovariate, by = "ScenID")
  }else{
    Covariates <- ExpCovSampleExpanded[ , c(.SD, IndCovariate), by = "ScenID"]
  }
  setkey(Covariates,"ScenID", "ExpID", "DoseID", "TrialID", "IndivID", "USUBJID")

  # This ID column will be used for joining with simPKPD tables
  Covariates[, ID:=1:.N]

  setcolorder(Covariates, c("ID",
                            "ScenID",
                            "ExpID",
                            "DoseID",
                            DOSEcovariate,
                            "TrialID",
                            "IndivID",
                            "USUBJID",
                            setdiff(c(names(Covariates)),
                                    c("ID",
                                      "ScenID",
                                      "DoseID",
                                      DOSEcovariate,
                                      "ExpID",
                                      "TrialID",
                                      "IndivID",
                                      "USUBJID"))))

  # Split covariates into a list (needed for step 6).
  Covariates_List <- split(Covariates, by = c("ScenID", "ExpID", "DoseID", "TrialID"))

  # Save covariate data chunks for each trial to trial files
  for(covdata in Covariates_List) {
    saveToTrialFile(outputFolder, covdata, listMemberName = "covariates", setting$FLAGreuseStoredTrialSummaryFiles)
  }

  # save some memory, because Covariates_List is no more needed
  rm(Covariates_List)

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")

  #-------------------------------------------#
  ## 4.d: Sample PKPD parameters  ----
  #-------------------------------------------#

  # Message:
  start_time <- Sys.time()
  cat("\t   4.d: Sample PKPD parameters\n")


  if(FLAGclinicalSimMode){
    #-------------------------------------------#
    ### 4.d.1: Clinical Mode  ----

    # Assign memory for PK/PD parameters:
    PKPDparameters <- data.table::data.table(ScenID  = rep(ScenID_Vec      , each  = N_Trial*N_SubjPerTrial),
                                             TrialID = rep(rep(1:N_Trial   , each  = N_SubjPerTrial),times=N_Scenario),
                                             IndivID = rep(1:N_SubjPerTrial, times = N_Scenario*N_Trial))
    PKPDparameters <- PKPDparameters[,USUBJID := ScenID*10^(nchar(N_Trial)+nchar(N_SubjPerTrial))+TrialID*10^(nchar(N_SubjPerTrial))+IndivID]
    PKPDparameters <- PKPDparameters[,(regressorName):=lapply(regressorName, function(x){0})]
    data.table::setkey(PKPDparameters, "ScenID", "TrialID", "IndivID", "USUBJID")

    # Adjust Class Type Based on Covariates:
    for (IndCovName_k in intersect(IndCovName,names(PKPDparameters))){
      PKPDparameters[,(IndCovName_k) := eval(parse(text=paste0("as.",IndCovClass[[IndCovName_k]],"(", IndCovName_k,")")))]
    }

    # Sample PKPD parameters:
    #   Dummy sampling to have the warning messages only once
    if (FLAG_SAMPLE %in% c(0,1)){
      PKPDdummy <- IQRtools::sample_GPF(fullModelGPF,
                                        Nsamples    = 2,
                                        FLAG_SAMPLE = FLAG_SAMPLE,
                                        covariates  = Covariates[c(1,2),],
                                        # covariates  = Covariates[TrialID==Covariates$TrialID[1],],
                                        FLAGid      = FALSE,
                                        verbose     = FALSE)$indParamValues
    }else if (FLAG_SAMPLE %in% c(2,3,4)){
      PKPDdummy <- IQRtools::sample_GPF(fullModelGPF,
                                        Nsamples    = 2,
                                        FLAG_SAMPLE = FLAG_SAMPLE,
                                        covariates  = Covariates[c(1,2),],
                                        # covariates  = Covariates[TrialID==Covariates$TrialID[1],],
                                        FLAGid      = FALSE,
                                        verbose     = FALSE)$popParamValues
    }else if (FLAG_SAMPLE %in% c(5)){
      PKPDdummy <- IQRtools::sample_GPF(fullModelGPF,
                                        Nsamples    = 2,
                                        FLAG_SAMPLE = FLAG_SAMPLE,
                                        covariates  = Covariates[c(1,2),],
                                        # covariates  = Covariates[TrialID==Covariates$TrialID[1],],
                                        FLAGid      = FALSE,
                                        verbose     = FALSE)$typicalIndParamValues

    }else{
      stop("Not a valid value for 'FLAG_SAMPLE': It should be 0, 1, 2, 3, 4 or 5")
    }
    #   Adjust Class Type Based on Sampled Parameters:
    for (PKPDpara_k in names(PKPDdummy)){
      PKPDparameters[,(PKPDpara_k) := eval(parse(text=paste0("as.",class(PKPDdummy[[PKPDpara_k]]),"(", PKPDpara_k,")")))]
    }
    #   Sample for real now
    PKPDparameters[,(regressorNameGPF):={

      # Get Trial ID:
      TrialID_k <- TrialID[1]

      # Sample:
      base::suppressWarnings(
        base::suppressMessages(
          if (FLAG_SAMPLE %in% c(0,1)){
            PKPDparameters_k <- IQRtools::sample_GPF(fullModelGPF,
                                                     Nsamples    = N_Scenario*N_SubjPerTrial,
                                                     FLAG_SAMPLE = FLAG_SAMPLE,
                                                     covariates  = Covariates[TrialID==TrialID_k,],
                                                     FLAGid      = FALSE,
                                                     verbose     = FALSE)$indParamValues
          }else if (FLAG_SAMPLE %in% c(2,3,4)){
            PKPDparameters_k <- IQRtools::sample_GPF(fullModelGPF,
                                                     Nsamples    = N_Scenario*N_SubjPerTrial,
                                                     FLAG_SAMPLE = FLAG_SAMPLE,
                                                     covariates  = Covariates[TrialID==TrialID_k,],
                                                     FLAGid      = FALSE,
                                                     verbose     = FALSE)$popParamValues
          }else if (FLAG_SAMPLE %in% c(5)){
            # TODO Remove
            PKPDparameters_k <- IQRtools::sample_GPF(fullModelGPF,
                                                     Nsamples    = N_Scenario*N_SubjPerTrial,
                                                     FLAG_SAMPLE = FLAG_SAMPLE,
                                                     covariates  = Covariates[TrialID==TrialID_k,],
                                                     FLAGid      = FALSE,
                                                     verbose     = FALSE)$typicalIndParamValues

          }else{
            stop("Not a valid value for 'FLAG_SAMPLE': It should be 0, 1, 2, 3, 4 or 5")
          }
        )
      )

      # data.table output:
      lapply(regressorNameGPF,function(col_k){PKPDparameters_k[[col_k]]})
    },
    by=c("TrialID")]

    # Add ExpID and DoseID:
    PKPDparameters <- PKPDparameters[ , c(.SD, Covariates[,c("ExpID","DoseID")])]

    # Free some memory:
    invisible(gc(verbose=FALSE))


  }else{
    #-------------------------------------------#
    ### 4.d.2: Compare Mode  ----
    #-------------------------------------------#

    # Assign memory for PK/PD parameters:
    PKPDparameters <- data.table::data.table(TrialID = rep(1:N_Trial       , each  = N_SubjPerTrial),
                                             IndivID = rep(1:N_SubjPerTrial, times = N_Trial))
    PKPDparameters <- PKPDparameters[,USUBJID := TrialID*10^(nchar(N_SubjPerTrial))+IndivID]
    PKPDparameters <- PKPDparameters[,(regressorName):=lapply(regressorName, function(x){0})]
    data.table::setkey(PKPDparameters, "TrialID", "IndivID", "USUBJID")

    # Adjust Class Type Based on Covariates:
    for (IndCovName_k in intersect(IndCovName,names(PKPDparameters))){
      PKPDparameters[,(IndCovName_k) := eval(parse(text=paste0("as.",IndCovClass[[IndCovName_k]],"(", IndCovName_k,")")))]
    }

    # Sample PKPD parameters:
    #   Dummy sampling to have the warning messages only onc
    # TODO: Make sure to use Covariates instead of IndCovariate. Covariates includes also the experimental and regimen covariates.
    # This includes removing step 4.e.
    if (FLAG_SAMPLE %in% c(0,1)){
      PKPDdummy <- IQRtools::sample_GPF(fullModelGPF,
                                        Nsamples    = 2,
                                        FLAG_SAMPLE = FLAG_SAMPLE,
                                        covariates  = IndCovariate[c(1,2),],
                                        # covariates  = IndCovariate[TrialID==IndCovariate$TrialID[1],],
                                        FLAGid      = FALSE,
                                        verbose     = FALSE)$indParamValues
    }else if (FLAG_SAMPLE %in% c(2,3,4)){
      PKPDdummy <- IQRtools::sample_GPF(fullModelGPF,
                                        Nsamples    = 2,
                                        FLAG_SAMPLE = FLAG_SAMPLE,
                                        covariates  = IndCovariate[c(1,2),],
                                        # covariates  = IndCovariate[TrialID==IndCovariate$TrialID[1],],
                                        FLAGid      = FALSE,
                                        verbose     = FALSE)$popParamValues
    }else if (FLAG_SAMPLE %in% c(5)){
      PKPDdummy <- IQRtools::sample_GPF(fullModelGPF,
                                        Nsamples    = 2,
                                        FLAG_SAMPLE = FLAG_SAMPLE,
                                        covariates  = IndCovariate[c(1,2),],
                                        # covariates  = IndCovariate[TrialID==IndCovariate$TrialID[1],],
                                        FLAGid      = FALSE,
                                        verbose     = FALSE)$typicalIndParamValues

    }else{
      stop("Not a valid value for 'FLAG_SAMPLE': It should be 0, 1, 2, 3, 4 or 5")
    }
    #   Adjust Class Type Based on Sampled Parameters:
    for (PKPDpara_k in names(PKPDdummy)){
      PKPDparameters[,(PKPDpara_k) := eval(parse(text=paste0("as.",class(PKPDdummy[[PKPDpara_k]]),"(", PKPDpara_k,")")))]
    }
    #   Sample for real now
    PKPDparameters[,(regressorNameGPF):={

      # Get Trial ID:
      TrialID_k <- TrialID[1]

      # Sample:
      base::suppressWarnings(
        base::suppressMessages(
          if (FLAG_SAMPLE %in% c(0,1)){
            PKPDparameters_k <- IQRtools::sample_GPF(fullModelGPF,
                                                     Nsamples    = N_SubjPerTrial,
                                                     FLAG_SAMPLE = FLAG_SAMPLE,
                                                     covariates  = IndCovariate[TrialID==TrialID_k,],
                                                     FLAGid      = FALSE,
                                                     verbose     = FALSE)$indParamValues
          }else if (FLAG_SAMPLE %in% c(2,3,4)){
            PKPDparameters_k <- IQRtools::sample_GPF(fullModelGPF,
                                                     Nsamples    = N_SubjPerTrial,
                                                     FLAG_SAMPLE = FLAG_SAMPLE,
                                                     covariates  = IndCovariate[TrialID==TrialID_k,],
                                                     FLAGid      = FALSE,
                                                     verbose     = FALSE)$popParamValues
          }else if (FLAG_SAMPLE %in% c(5)){
            PKPDparameters_k <- IQRtools::sample_GPF(fullModelGPF,
                                                     Nsamples    = N_SubjPerTrial,
                                                     FLAG_SAMPLE = FLAG_SAMPLE,
                                                     covariates  = IndCovariate[TrialID==TrialID_k,],
                                                     FLAGid      = FALSE,
                                                     verbose     = FALSE)$typicalIndParamValues

          }else{
            stop("Not a valid value for 'FLAG_SAMPLE': It should be 0, 1, 2, 3, 4 or 5")
          }
        )
      )

      # data.table output:
      lapply(regressorNameGPF,function(col_k){PKPDparameters_k[[col_k]]})
    },
    by=c("TrialID")]

    # Free some memory:
    invisible(gc(verbose=FALSE))
  }

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #-------------------------------------------#
  ### 4.e: (Compare mode only) Duplicate PKPDparameters by Experimental Setting  ----
  #-------------------------------------------#
  # TODO: This step should disappear if using Covariates everywhere in 4.d.

  # Message:
  start_time <- Sys.time()
  cat("\t   4.e: Duplicate 'PKPDparameters' by Experimental Setting\n")


  if(FLAGclinicalSimMode){
    cat("\t\t\tThis step was skipped as the simulations are done within the clinical mode\n")
  }else{
    # Get the pop PKPD parameters for each scenario
    PKPDparameters <- ExpCovSampleExpanded[,
                                           c(.SD, PKPDparameters),
                                           .SDcols = c("ExpID", "DoseID"),
                                           by = "ScenID"]

    # pop PK parameters for the reference scenario:
    PKPDparameters.POP0        <- fullModelGPF$estimates$VALUE[idx_mpar]
    names(PKPDparameters.POP0) <- fullModelGPF$estimates$PARAMETER[idx_mpar]
    PKPDparameters.POP0        <- as.data.frame(t(PKPDparameters.POP0))

    # pop PKPD parameters for each scenario:
    base::suppressWarnings(
      base::suppressMessages(
        PKPDparameters.POP <- data.table::setDT(IQRtools::sample_GPF(fullModelGPF,
                                                                     covariates  = ExpCovSampleExpanded,
                                                                     FLAG_SAMPLE = 5,
                                                                     verbose     = FALSE)$typicalIndParamValues )
      )
    )

    # Add IDs:
    PKPDparameters.POP <- cbind(ExpCovSampleExpanded[,.SD,.SDcols = c("ScenID", "ExpID", "DoseID")],
                                PKPDparameters.POP)

    # Update Indiv PKPD Parameters:
    #   Transform to normal scale
    PKPDparameters      <- IQRtools:::transformParamToNormal(fullModelGPF, PKPDparameters     , columnNames = regressorNameGPF)
    PKPDparameters.POP  <- IQRtools:::transformParamToNormal(fullModelGPF, PKPDparameters.POP , columnNames = regressorNameGPF)
    PKPDparameters.POP0 <- IQRtools:::transformParamToNormal(fullModelGPF, PKPDparameters.POP0, columnNames = regressorNameGPF)
    #   Correct
    PKPDparameters[,(regressorNameGPF):={
      ScenID_k <- ScenID[1]
      lapply(regressorNameGPF, function(reg_k){
        get(reg_k) + PKPDparameters.POP[ScenID==ScenID_k,get(reg_k)] - PKPDparameters.POP0[[reg_k]]
      })},
      by = "ScenID"]
    #   Transform back to original scale
    PKPDparameters <- IQRtools:::untransformParamFromNormal(fullModelGPF, PKPDparameters, columnNames = regressorNameGPF)
    PKPDparameters <- setDT(PKPDparameters)

    # Remove Variables:
    rm(PKPDparameters.POP, PKPDparameters.POP0)
  }

  # Set Keys:
  setkey(PKPDparameters,"ScenID", "ExpID", "DoseID", "TrialID","IndivID","USUBJID")
  setcolorder(PKPDparameters, c("ScenID",
                                "ExpID",
                                "DoseID",
                                "TrialID",
                                "IndivID",
                                "USUBJID",
                                regressorName))

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #-------------------------------------------#
  ## 4.f: Adjust PKPDparameters  ----
  #-------------------------------------------#

  # Message:
  start_time <- Sys.time()
  cat("\t   4.f: Adjust 'PKPDparameters'\n")

  #-------------------------------------------#
  ### 4.f.1: Adjust Parameters with 'IndCovariate'  ----
  # TODO this case is already covered by 4.f.2
  #   Only if needed
  regressorName_Ind <- intersect(regressorName, names(IndCovariate))
  if (length(regressorName_Ind)>0){
    PKPDparameters[,(regressorName_Ind):={
      lapply(regressorName_Ind,function(col_k){
        Covariates[,get(col_k)]
      })}]
  }

  #-------------------------------------------#
  ### 4.f.2: Adjust Parameters with 'ExpCovSample' ----
  #   Only if needed
  regressorName_Exp <- intersect(regressorName, names(Covariates))
  if (length(regressorName_Exp)>0){
    # TODO check that nrow(Covariates) == nrow(PKPDparameters)
    PKPDparameters[,(regressorName_Exp):={
      lapply(regressorName_Exp,function(col_k){
        Covariates[,get(col_k)]
      })}]
  }

  #-------------------------------------------#
  ### 4.f.3: Add Fpediatric if necessary ----
  #   Only if needed
  if (!is.null(Fpediatric)){
    if (Fpediatric %in% names(Covariates)){
      # Adjust Fpediatric first:
      setnames(Covariates, Fpediatric, "Fpediatric")

      # Add Fpediatric to PKPDparameters:
      PKPDparameters[,Fpediatric:={
        .(Covariates[,Fpediatric])}]
    }
  }

  # Reorder the columns PKPDparameters:
  PKPDparameters[, ID:=(1:.N)]
  setcolorder(PKPDparameters, c("ID",
                                "ScenID",
                                "ExpID",
                                "DoseID",
                                "TrialID",
                                "IndivID",
                                "USUBJID",
                                setdiff(c(names(PKPDparameters)),
                                        c("ID",
                                          "ScenID",
                                          "ExpID",
                                          "DoseID",
                                          "TrialID",
                                          "IndivID",
                                          "USUBJID"))))

  # Free some memory:
  rm(IndCovariate, ExpCovSampleExpanded)
  invisible(gc(verbose=FALSE))

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #---------------------------------------------------#
  # STEP 5: Create event Data ----
  #---------------------------------------------------#
  # Message:
  start_time <- Sys.time()
  cat("\tSTEP 5: Generate eventData table\n")

  # Generate Event Data:
  byCol <- intersect(names(dosing),names(PKPDparameters))
  byCol <- if (length(byCol)==0) NULL else byCol
  eventData <- merge(dosing,PKPDparameters, by=byCol, allow.cartesian=TRUE)

  # Set key:
  setkey(eventData, "ID", "ScenID", "ExpID", "DoseID", "TrialID", "IndivID", "USUBJID", "TIME")

  # Correct Dose with Fpediatric if necessary:
  if(!is.null(Fpediatric)){
    eventData[,AMT:=(AMT*Fpediatric)]
    eventData[     , Fpediatric:= NULL]
    PKPDparameters[, Fpediatric:= NULL]
  }

  # Remember dosing, experiment setting and ID of each Subject:
  dosID <- unique(eventData[, c("ID", "ScenID", "DoseID", "ExpID", "TrialID", "IndivID", "USUBJID")])

  # Check if absorption of 0-order:
  abs0inputs   <- NULL
  abs0Tk0param <- NULL
  if (!is.null(Tk0)){
    if (any(Tk0 %in% regressorName)){
      idx_Tk0      <- (Tk0 %in% regressorName)
      abs0inputs   <- as.numeric(gsub("INPUT","",names(Tk0)[idx_Tk0]))
      abs0Tk0param <- c(Tk0[idx_Tk0], use.names=FALSE)
    }
  }

  # Create input event table:
  eventTable <- eventData[,{
    # Dosing Information:
    L1 <- list(TIME = TIME,
               ADM  = ADM,
               AMT  = AMT)

    # Infusion Information
    L2 <- list()
    L2[["TINF"]] <- rep(1e-4,.N)
    if(!is.null(abs0inputs)){
      for(k in 1:length(abs0Tk0param)){
        Tk0_k   <- abs0Tk0param[k]
        input_k <- abs0inputs[k]
        L2[["TINF"]] <- ifelse(ADM==input_k,
                               get(Tk0_k),
                               L2[["TINF"]])
      }
    }

    # Regressors:
    L3 <- lapply(setNames(regressorName,regressorName),function(reg_k){
      ifelse(TIME==min(TIME) & ADM==min(ADM[TIME==min(TIME)]),get(reg_k), NA)
    })

    # Output:
    c(L1, L2, L3)

  },by=c("ScenID", "ExpID", "DoseID", "TrialID", "ID")]

  # Set keys:
  setkey(eventTable, "ID", "TIME", "ADM")

  # Split eventTable :
  eventTable_List <- split(eventTable, by= c("ScenID", "ExpID", "DoseID", "TrialID"), keep.by = TRUE)
  eventTable_List <- lapply(eventTable_List, function(eventTable_k){
    `class<-` (eventTable_k, c("IQReventTable", "data.frame"))
  })

  # save eventTable chunk for each trial into trial file.
  for(eventTable_k in eventTable_List) {
    saveToTrialFile(outputFolder, eventTable_k, listMemberName = "eventTable", setting$FLAGreuseStoredTrialSummaryFiles)
  }

  # free some memory, because eventTable_List is no more needed
  rm(eventTable_List)
  invisible(gc(verbose=FALSE))

  # Remove columns not acceptable for an eventTable
  eventTable[, ScenID:=NULL]
  eventTable[, ExpID:=NULL]
  eventTable[, DoseID:=NULL]
  eventTable[, TrialID:=NULL]
  eventTable        <- data.table::setDF(eventTable)
  class(eventTable) <- c("IQReventTable", "data.frame")

  # Measure Elapsed Time:
  end_time  <- Sys.time()
  Diff_time <- end_time - start_time
  cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")


  #---------------------------------------------------#
  # STEP 6: Simulations ----
  #---------------------------------------------------#

  if(FLAGrunSimulation){
    # Message:
    start_time <- Sys.time()
    cat("\tSTEP 6: Run Simulations\n")

    # Re-estimate Nparallel:
    maxCores     <- parallel::detectCores(logical = FALSE)
    NparallelNew <- min(Nparallel,maxCores)

    # Parallelisation:
    if (NparallelNew>1 && maxIDrun>=4) {
      # Flag:
      parallelFLAG <- TRUE

      # Start cluster:
      cluster_Sim <- parallel::makeCluster(NparallelNew)
      doParallel::registerDoParallel(cluster_Sim)

      # Try this line below to provide libPaths to all nodes:
      parallel::clusterCall(cl=cluster_Sim, ".libPaths", .libPaths())

      # Export needed objects from the global environment of the master to the global environment of the worker porcesses.
      if(!is.null(setting$.paropts$.export)) {
        parallel::clusterExport(cl = cluster_Sim, varlist = setting$.paropts$.export)
      }
      # Message:
      cat("\t\tParallelisation Activated:", NparallelNew, "cores are used in parallel.\n\n")

    } else {
      # Flag:
      parallelFLAG <- FALSE
      cat("\t\tNo Parallelisation.\n\n")
    }

    trialFilenames <- list.files(outputFolder, pattern = "trial.*.rds", full.names = TRUE)
    # Simulations:
    summary_simPKPD <- plyr::llply(
      trialFilenames,
      .parallel = parallelFLAG,
      .paropts = setting$.paropts,
      .fun = simulate_OneVirtualTrial)

    # Close clusters:
    if (parallelFLAG){
      parallel::stopCluster(cl = cluster_Sim)
    }
    invisible(gc(verbose=FALSE))

    # Collect summaries by trials into data.tables
    summaryPKPD.ByTrial    <- rbindlist(lapply(summary_simPKPD, function(sim_ScenIDTrialID) sim_ScenIDTrialID$summary_simPKPD$summaryPKPD.ByTrial))
    summaryClinEnd.ByIndiv <- rbindlist(lapply(summary_simPKPD, function(sim_ScenIDTrialID) sim_ScenIDTrialID$summary_simPKPD$summaryClinEnd.ByIndiv))
    summaryClinEnd.ByTrial <- rbindlist(lapply(summary_simPKPD, function(sim_ScenIDTrialID) sim_ScenIDTrialID$summary_simPKPD$summaryClinEnd.ByTrial))
    simPKPD                <- rbindlist(lapply(summary_simPKPD, function(sim_ScenIDTrialID) sim_ScenIDTrialID$simPKPD))

    # Measure Elapsed Time:
    end_time  <- Sys.time()
    Diff_time <- end_time - start_time
    cat("\t\tThis step took", ifelse(Diff_time<1, signif(Diff_time,2), round(Diff_time,1)), units(Diff_time),"\n")

  }else{
    cat("\tSTEP 6: Simulations not run\n")
    cat("\t\tAs requested\n")
    summaryPKPD.ByTrial    <- NULL
    summaryClinEnd.ByIndiv <- NULL
    summaryClinEnd.ByTrial <- NULL
    simPKPD                <- NULL
  }

  #---------------------------------------------------#
  # STEP 7: Output ----
  #---------------------------------------------------#

  # Message:
  cat("\tSTEP 7: Generate Output\n")
  cat("\t\tDone\n")
  cat("----------------------------------------------------\n\n")


  data.table::setDF(Covariates)
  data.table::setDF(PKPDparameters)
  # data.table::setDF(eventTable)  -- Already done
  data.table::setDF(eventData)

  # Create output:
  out <- list(summaryPKPD.ByTrial    = summaryPKPD.ByTrial,
              summaryClinEnd.ByIndiv = summaryClinEnd.ByIndiv,
              summaryClinEnd.ByTrial = summaryClinEnd.ByTrial,
              simPKPD                = simPKPD,
              Covariates             = Covariates,
              PKPDparameters         = PKPDparameters,
              eventTable             = eventTable,
              eventData              = eventData,
              regression             = regressorName,
              abs0inputs             = abs0inputs,
              abs0Tk0param           = abs0Tk0param,
              setting                = setting)

  # Output:
  return(out)
}

validateSetting_SimulateVirtualTrials <- function(setting = NULL) {

  if(is.character(setting)) {
    # Load setting from a file
    setting <- try(readRDS(setting), silent = TRUE)
    if(inherits(setting, "try-error")) {
      stop("validateSetting_SimulateVirtualTrials: error while loading setting file ", setting)
    }
  }
  if(!is.null(setting) && !is.list(setting)) {
    stop("validateSetting_SimulateVirtualTrials: setting should either be NULL or a list.")
  }
  if(is.null(setting)) {
    setting <- list()
  }

  validate <- function(setting, name, default, type = class(default), FLAGnullAllowed = is.null(default)) {
    if(name %in% names(setting)) {
      if( type == "logical" && !(is.logical(setting[[name]]) || (FLAGnullAllowed && is.null(setting[[name]]))) ) {
        stop(paste0("validateSetting_SimulateVirtualTrials: setting$", name, " should be logical", if(FLAGnullAllowed) " or NULL" else ""))
      }
      if(type == "numeric" && !(is.numeric(setting[[name]]) || (FLAGnullAllowed && is.null(setting[[name]]))) ) {
        stop(paste0("validateSetting_SimulateVirtualTrials: setting$", name, " should be numeric", if(FLAGnullAllowed) " or NULL" else ""))
      }
      if(type == "character" && !(is.character(setting[[name]]) || (FLAGnullAllowed && is.null(setting[[name]]))) ) {
        stop(paste0("validateSetting_SimulateVirtualTrials: setting$", name, " should be character", if(FLAGnullAllowed) " or NULL" else ""))
      }
      if(type == "list" && !(is.list(setting[[name]]) || (FLAGnullAllowed && is.null(setting[[name]]))) ) {
        stop(paste0("validateSetting_SimulateVirtualTrials: setting$", name, " should be list", if(FLAGnullAllowed) " or NULL" else ""))
      }
    } else {
      setting[[name]] <- default
    }
    setting
  }

  setting <- validate(setting, "simtime", seq(0,672,1))

  setting <- validate(setting, "FLAGsensitivity", FALSE)
  setting <- validate(setting, "FLAGuseSensEq", TRUE)
  setting <- validate(setting, "sensParams", NULL, "list")
  setting <- validate(setting, "FLAGoutputsOnly", FALSE)

  setting <- validate(setting, "FLAGreturnSimPKPD"               , FALSE)
  setting <- validate(setting, "FLAGsaveSimPKPDToFile"           , FALSE)
  setting <- validate(setting, "FLAGreuseStoredTrialSummaryFiles", FALSE)

  setting <- validate(setting, "fun_SummarizeTrial", "MMVmalaria::summarizeTrial_Basic")
  setting <- validate(setting, "args_SummarizeTrial", NULL, "list")

  setting <- validate(setting, "opt_eventTimes"           , NULL , "numeric")
  setting <- validate(setting, "opt_method_stiff"         , TRUE )
  setting <- validate(setting, "opt_abstol"               , 1e-09)
  setting <- validate(setting, "opt_reltol"               , 1e-06)
  setting <- validate(setting, "opt_minstep"              , 0    )
  setting <- validate(setting, "opt_maxstep"              , 0    )
  setting <- validate(setting, "opt_initstep"             , 0    )
  setting <- validate(setting, "opt_maxnumsteps"          , 1e+05)
  setting <- validate(setting, "opt_maxerrtestfails"      , 50   )
  setting <- validate(setting, "opt_maxorder_stiff"       , 5    )
  setting <- validate(setting, "opt_maxorder_nonstiff"    , 12   )
  setting <- validate(setting, "opt_maxconvfails"         , 10   )
  setting <- validate(setting, "opt_maxnonlineariter"     , 3    )
  setting <- validate(setting, "opt_usesymjac"            , TRUE )
  setting <- validate(setting, "opt_sens_simultaneous"    , FALSE)
  setting <- validate(setting, "opt_sens_errcon"          , FALSE)
  setting <- validate(setting, "opt_sens_maxnonlineariter", 3    )

  setting <- validate(setting, "verbose", FALSE)

  setting
}

#' Helper function: Load a named trial file-list for a given trial ScenID, ExpID, DoseID, TrialID
loadTrialFile <- function(trialFilename) {
  if(!file.exists(trialFilename)) {
    stop("loadTrialFile: trial file-name not found ", trialFilename)
  } else {
    readRDS(trialFilename)
  }
}

#' Helper function: Save a data.frame or a named list to a trial file
saveToTrialFile <- function(outputFolder, data, listMemberName, FLAGreuseStoredTrialSummaryFiles) {
  if(is.data.frame(data)) {
    ScenID  <- unique(data$ScenID)
    ExpID   <- unique(data$ExpID)
    DoseID  <- unique(data$DoseID)
    TrialID <- unique(data$TrialID)
  } else if(is.list(data)) {
    ScenID  <- data$ScenID
    ExpID   <- data$ExpID
    DoseID  <- data$DoseID
    TrialID <- data$TrialID
  } else {
    stop("saveToTrialFile: argument data should be either be a data.frame or a named list.")
  }

  stopifnot(
    "saveToTrialFile: failed to generate trialFilename because of non-unique ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)
  trialFilename <- file.path(outputFolder, paste0("trial_", ScenID, "_", ExpID, "_", DoseID, "_", TrialID, ".rds"))

  if(file.exists(trialFilename)) {
    lst <- as.list(readRDS(trialFilename))
  } else {
    lst <- list(
      ScenID  = ScenID,
      ExpID   = ExpID,
      DoseID  = DoseID,
      TrialID = TrialID)
  }

  if(is.data.frame(data)) {
    if(!is.character(listMemberName)) {
      stop("saveToTrialFile: if data is a data.frame listMemberName should be the name of this data.frame in a trial list .rds file.")
    }
    if(!is.null(lst[[listMemberName]])) {
      if(!isTRUE(all.equal(lst[[listMemberName]], data)) && FLAGreuseStoredTrialSummaryFiles) {
        stop("saveToTrialFile: a member of the same name is already stored, FLAGreuseStoredTrialSummaryFiles is TRUE, but stored data.frame is not all-equal to the passed data.frame; check: ", trialFilename, "$", listMemberName)
      } else {
        lst[[listMemberName]] <- data
        saveRDS(lst, file = trialFilename)
      }
    } else {
      lst[[listMemberName]] <- data
      saveRDS(lst, file = trialFilename)
    }
  } else {
    # data is a named list
    for(name in names(data)) {
      if(name %in% names(lst) && !all.equal(lst[[name]], data[[name]]) && FLAGreuseStoredTrialSummaryFiles) {
        stop("saveToTrialFile: a member of the same name is already stored, FLAGreuseStoredTrialSummaryFiles is TRUE, but stored data.frame is not all-equal to the passed data.frame; check: ", trialFilename, "$", name)
      } else {
        lst[[name]] <- data[[name]]
        saveRDS(lst, file = trialFilename)
      }
    }
  }
}

#' Simulate one virtual trial based on a input .rds file
#' @param trialFileName A character string denoting a file path to a trial .rds file. By convention, the general settings for the simulation
#' are loaded from a file setting.rds found in the parent directory of trialFileName.
#'
#' @return a named list
simulate_OneVirtualTrial <- function(trialFileName){

  # validate trial folder ----
  trialFolder <- dirname(trialFileName)
  if(!dir.exists(trialFolder)) {
    stop("simulate_OneVirtualTrial: Seemingly the directory of the specified trialFileName (", trialFilename, ") does not exist.")
  }

  # validate setting ----
  setting <- validateSetting_SimulateVirtualTrials(setting = file.path(trialFolder, "setting.rds"))

  # Get the function to summarize the simulation results ----
  if(!is.character(setting$fun_SummarizeTrial)) {
    stop("simulate_VirtualTrials: setting$fun_SummarizeTrial should be a character string denoting a function name.")
  }
  f_SummarizeTrial <- try(eval(parse(text = setting$fun_SummarizeTrial)), silent = TRUE)
  if(!is.function(f_SummarizeTrial)) {
    stop("simulate_VirtualTrials: ", setting$fun_SummarizeTrial, " should evaluate to a function, but evaluated to ", toString(f_SummarizeTrial))
  }


  # load trial file ----
  trial_k <- loadTrialFile(trialFileName)

  ScenID  <- trial_k$ScenID
  ExpID   <- trial_k$ExpID
  DoseID  <- trial_k$DoseID
  TrialID <- trial_k$TrialID

  stopifnot(
    "simulate_OneVirtualTrial: covariates data for each iteration should be for a single ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)


  if(setting$FLAGreuseStoredTrialSummaryFiles && !is.null(trial_k$summary_simPKPD)) {
    if(!is.null(trial_k$simPKPD) && setting$FLAGreturnSimPKPD) {
      return(list(ScenID          = ScenID,
                  ExpID           = ExpID,
                  DoseID          = DoseID,
                  TrialID         = TrialID,
                  summary_simPKPD = trial_k$summary_simPKPD,
                  simPKPD         = trial_k$simPKPD))
    } else {
      return(list(ScenID          = ScenID,
                  ExpID           = ExpID,
                  DoseID          = DoseID,
                  TrialID         = TrialID,
                  summary_simPKPD = trial_k$summary_simPKPD))
    }
  } else {
    # load model ----
    modelIQR <- try(IQRtools::IQRmodel(file.path(trialFolder, "model.txt")), silent = TRUE)
    if(inherits(modelIQR, "try-error")) {
      stop("simulate_OneVirtualTrial: error reloading IQRmodel object:", modelIQR)
    }

    eventTable_k <- data.table::setDT(trial_k$eventTable)
    covariates_k <- trial_k$covariates

    stopifnot(setequal(unique(eventTable_k$ID), unique(covariates_k$ID)))

    # Remove columns not acceptable for an eventTable
    eventTable_k[, ScenID:=NULL]
    eventTable_k[, ExpID:=NULL]
    eventTable_k[, DoseID:=NULL]
    eventTable_k[, TrialID:=NULL]
    eventTable_k        <- data.table::setDF(eventTable_k)
    class(eventTable_k) <- c("IQReventTable", "data.frame")

    # Run Simulations:
    simPKPD_k <- try(data.table::setDT(sim_IQRmodel(modelIQR,
                                                    simtime                   = setting$simtime,
                                                    eventTable                = eventTable_k,
                                                    FLAGsensitivity           = setting$FLAGsensitivity,
                                                    FLAGuseSensEq             = setting$FLAGuseSensEq,
                                                    sensParams                = setting$sensParams,
                                                    FLAGoutputsOnly           = setting$FLAGoutputsOnly,
                                                    opt_eventTimes            = setting$opt_eventTimes,
                                                    opt_method_stiff          = setting$opt_method_stiff,
                                                    opt_abstol                = setting$opt_abstol,
                                                    opt_reltol                = setting$opt_reltol,
                                                    opt_minstep               = setting$opt_minstep,
                                                    opt_maxstep               = setting$opt_maxstep,
                                                    opt_initstep              = setting$opt_initstep,
                                                    opt_maxnumsteps           = setting$opt_maxnumsteps,
                                                    opt_maxerrtestfails       = setting$opt_maxerrtestfails,
                                                    opt_maxorder_stiff        = setting$opt_maxorder_stiff,
                                                    opt_maxorder_nonstiff     = setting$opt_maxorder_nonstiff,
                                                    opt_maxconvfails          = setting$opt_maxconvfails,
                                                    opt_maxnonlineariter      = setting$opt_maxnonlineariter,
                                                    opt_usesymjac             = setting$opt_usesymjac,
                                                    opt_sens_simultaneous     = setting$opt_sens_simultaneous,
                                                    opt_sens_errcon           = setting$opt_sens_errcon,
                                                    opt_sens_maxnonlineariter = setting$opt_sens_maxnonlineariter,
                                                    verbose                   = setting$verbose)),
                     silent = TRUE)
    if(inherits(simPKPD_k, "try-error")) {
      stop("simulate_OneVirtualTrial: error in sim_IQRmodel:", simPKPD_k)
    }
    if(!("ID" %in% names(simPKPD_k))){
      simPKPD_k$ID <- covariates_k$ID[1]
    }

    # Add covariates and ID columns to simPKPD:
    # TODO to prevent the case of multiple rows due to time varying covariates
    simPKPD_k <- merge(simPKPD_k,
                       covariates_k, by = "ID")


    # summarize simulation result
    summary_simPKPD_k <- try(do.call(f_SummarizeTrial,
                                     c(list(simPKPD = simPKPD_k, trialFileName = trialFileName), setting$args_SummarizeTrial)),
                             silent = TRUE)

    if(inherits(summary_simPKPD_k, "try-error")) {
      stop("simulate_OneVirtualTrial: Summarizing simulation results for ScenID=", ScenID, " ExpID=", ExpID, ", DoseID=", DoseID, ", TrialID=", TrialID, " failed with error ", summary_simPKPD_k, ".")
    }

    saveToTrialFile(outputFolder = trialFolder,
                    data = if(setting$FLAGsaveSimPKPDToFile) {
                      list(ScenID          = ScenID,
                           ExpID           = ExpID,
                           DoseID          = DoseID,
                           TrialID         = TrialID,
                           summary_simPKPD = summary_simPKPD_k,
                           simPKPD         = simPKPD_k)
                    } else {
                      list(ScenID          = ScenID,
                           ExpID           = ExpID,
                           DoseID          = DoseID,
                           TrialID         = TrialID,
                           summary_simPKPD = summary_simPKPD_k)
                    },
                    FLAGreuseStoredTrialSummaryFiles = setting$FLAGreuseStoredTrialSummaryFiles)

    if(setting$FLAGreturnSimPKPD) {
      return(list(ScenID          = ScenID,
                  ExpID           = ExpID,
                  DoseID          = DoseID,
                  TrialID         = TrialID,
                  summary_simPKPD = summary_simPKPD_k,
                  simPKPD         = simPKPD_k))
    } else {
      return(list(ScenID          = ScenID,
                  ExpID           = ExpID,
                  DoseID          = DoseID,
                  TrialID         = TrialID,
                  summary_simPKPD = summary_simPKPD_k))
    }
  }
}
