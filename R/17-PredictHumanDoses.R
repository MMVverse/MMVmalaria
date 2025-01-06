#' Predict doses for a target time above MIC for a number of sampled virtual trial populations
#'
#' @param targetTimeAboveMIC a scalar specifying the target time for the compound concentration to stay above MIC.
#' @inheritParams predictDose_GenericTrial
#'
#' @return Data frame with coumns Variable, Metric, CI Low, CI Median, CI High, CI Level. The Variable column equals tMIC and
#' the Metric values correspond to the different values in the percentiles argument. The other columns correspond to the
#' predicted doses for each percentile in Metric.
#'
#' @author Venelin Mitov, IntiQuan
#' @export
#' @md
predictDose_TimeAboveMIC <- function(
  targetTimeAboveMIC,
  doseInterval,
  modelFile,

  projectPath,

  simtime,
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,
  args_IQRdosing = list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1, II = timeInterval),

  N_Trial,
  N_SubjperTrial,
  covariates     = NULL,
  DOSEcovariate  = NULL,
  Fpediatric     = NULL,
  regressorExtra = character(0),
  FLAG_SAMPLE    = 1,

  percentiles = c(5,50,95),
  CIlevel     = 90,

  addArgs_sim_IQRmodel = NULL,
  addArgs_uniroot = list(tol = 1),

  fun_EvaluateCriterion = "evaluateDoseCriterion_TimeAboveMIC",

  Nparallel = 1,
  seed = 1234567,
  FLAG_BROWSER1 = FALSE) {

  summaryPred <- predictDose_GenericTrial(fun_EvaluateCriterion = fun_EvaluateCriterion,
                                          targetCriterionValue  = targetTimeAboveMIC,
                                          doseInterval          = doseInterval,
                                          modelFile             = modelFile,

                                          projectPath           = projectPath,

                                          simtime               = simtime,
                                          Tk0                   = Tk0,
                                          nbrDoses              = nbrDoses,
                                          timeInterval          = timeInterval,
                                          args_IQRdosing        = args_IQRdosing,
                                          N_Trial               = N_Trial,
                                          N_SubjperTrial        = N_SubjperTrial,
                                          covariates            = covariates,
                                          DOSEcovariate         = DOSEcovariate,
                                          Fpediatric            = Fpediatric,
                                          regressorExtra        = regressorExtra,
                                          FLAG_SAMPLE           = FLAG_SAMPLE,

                                          percentiles           = percentiles,
                                          CIlevel               = CIlevel,

                                          addArgs_sim_IQRmodel  = addArgs_sim_IQRmodel,
                                          addArgs_uniroot       = addArgs_uniroot,

                                          Nparallel             = Nparallel,
                                          seed                  = seed,
                                          FLAG_BROWSER1         = FLAG_BROWSER1)
  summaryPred$Criteria <- paste0("tMIC=",targetTimeAboveMIC/24,"days")
  summaryPred
}


#' Predict doses for a target time above MPC90 for a number of sampled virtual trial populations
#'
#' @param targetTimeAboveMPC90 a scalar specifying the target time for the compound concentration to stay above MPC90.
#' @inheritParams predictDose_GenericTrial
#'
#' @return Data frame with coumns Variable, Metric, CI Low, CI Median, CI High, CI Level. The Variable column equals tMPC90 and
#' the Metric values correspond to the different values in the percentiles argument. The other columns correspond to the
#' predicted doses for each percentile in Metric.
#'
#' @author Venelin Mitov, IntiQuan
#' @export
#' @md
predictDose_TimeAboveMPC90 <- function(
  targetTimeAboveMPC90,
  doseInterval,
  modelFile,

  projectPath,

  simtime,
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,
  args_IQRdosing = list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1, II = timeInterval),

  N_Trial,
  N_SubjperTrial,
  covariates     = NULL,
  DOSEcovariate  = NULL,
  Fpediatric     = NULL,
  regressorExtra = character(0),
  FLAG_SAMPLE    = 1,

  percentiles = c(5,50,95),
  CIlevel     = 90,

  addArgs_sim_IQRmodel = NULL,
  addArgs_uniroot = list(tol = 1),

  fun_EvaluateCriterion = "evaluateDoseCriterion_TimeAboveMPC90",

  Nparallel = 1,
  seed = 1234567,
  FLAG_BROWSER1 = FALSE) {

  summaryPred <- predictDose_GenericTrial(fun_EvaluateCriterion = fun_EvaluateCriterion,
                                          targetCriterionValue  = targetTimeAboveMPC90,
                                          doseInterval          = doseInterval,
                                          modelFile             = modelFile,

                                          projectPath           = projectPath,

                                          simtime               = simtime,
                                          Tk0                   = Tk0,
                                          nbrDoses              = nbrDoses,
                                          timeInterval          = timeInterval,
                                          args_IQRdosing        = args_IQRdosing,
                                          N_Trial               = N_Trial,
                                          N_SubjperTrial        = N_SubjperTrial,
                                          covariates            = covariates,
                                          DOSEcovariate         = DOSEcovariate,
                                          Fpediatric            = Fpediatric,
                                          regressorExtra        = regressorExtra,
                                          FLAG_SAMPLE           = FLAG_SAMPLE,

                                          percentiles           = percentiles,
                                          CIlevel               = CIlevel,

                                          addArgs_sim_IQRmodel  = addArgs_sim_IQRmodel,
                                          addArgs_uniroot       = addArgs_uniroot,

                                          Nparallel             = Nparallel,
                                          seed                  = seed,
                                          FLAG_BROWSER1         = FLAG_BROWSER1)
  summaryPred$Criteria <- paste0("tMPC90=",targetTimeAboveMPC90/24,"days")
  summaryPred
}

#' Predict doses for a target PRRtot for a number of sampled virtual trial populations
#'
#' @param targetPRRtot a scalar specifying the target PRRtot.
#' @inheritParams predictDose_GenericTrial
#'
#' @return Data frame with coumns Variable, Metric, CI Low, CI Median, CI High, CI Level. The Variable column equals PRRtot and
#' the Metric values correspond to the different values in the percentiles argument. The other columns correspond to the
#' predicted doses for each percentile in Metric.
#'
#' @author Venelin Mitov, IntiQuan
#' @export
#' @md
predictDose_PRRtot <- function(
  targetPRRtot,
  doseInterval,
  modelFile,

  projectPath,

  simtime,
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,
  args_IQRdosing = list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1, II = timeInterval),

  N_Trial,
  N_SubjperTrial,
  covariates     = NULL,
  DOSEcovariate  = NULL,
  Fpediatric     = NULL,
  regressorExtra = character(0),
  FLAG_SAMPLE    = 1,

  percentiles = c(5,50,95),
  CIlevel     = 90,

  addArgs_sim_IQRmodel = NULL,
  addArgs_uniroot = list(tol = 1),

  fun_EvaluateCriterion = "evaluateDoseCriterion_PRRtot",

  Nparallel = 1,
  seed = 1234567,
  FLAG_BROWSER1 = FALSE) {

  summaryPred <- predictDose_GenericTrial(fun_EvaluateCriterion = fun_EvaluateCriterion,
                                          targetCriterionValue  = targetPRRtot,
                                          doseInterval          = doseInterval,
                                          modelFile             = modelFile,

                                          projectPath           = projectPath,

                                          simtime               = simtime,
                                          Tk0                   = Tk0,
                                          nbrDoses              = nbrDoses,
                                          timeInterval          = timeInterval,
                                          args_IQRdosing        = args_IQRdosing,
                                          N_Trial               = N_Trial,
                                          N_SubjperTrial        = N_SubjperTrial,
                                          covariates            = covariates,
                                          DOSEcovariate         = DOSEcovariate,
                                          Fpediatric            = Fpediatric,
                                          regressorExtra        = regressorExtra,
                                          FLAG_SAMPLE           = FLAG_SAMPLE,

                                          percentiles           = percentiles,
                                          CIlevel               = CIlevel,

                                          addArgs_sim_IQRmodel  = addArgs_sim_IQRmodel,
                                          addArgs_uniroot       = addArgs_uniroot,

                                          Nparallel             = Nparallel,
                                          seed                  = seed,
                                          FLAG_BROWSER1         = FLAG_BROWSER1)
  summaryPred$Criteria <- paste0("PRRtot=",targetPRRtot,"-log10")
  summaryPred
}


#' Predict doses reaching a minimum target criterion value for a number of sampled virtual trial populations
#'
#' @param fun_EvaluateCriterion a function object or a character string denoting the name of a function accessible
#' from the calling (parent) environment. This function should have the following synopsis:
#' \code{function(sim_results, parameters) { ... }}, where \code{sim_results} is a data.frame resulting from a
#' previous call to \code{\link{IQRtools::sim_IQRmodel}} and \code{parameters} is a data.frame with a single row
#' and numeric columns named as the model parameters; the function should return a numeric value. If
#' \code{fun_EvaluateCriterion} is a character string, the function is searched for via
#' \code{mget(fun_EvaluateCriterion, mode = "function", ifnotfound = list(NULL), inherits = TRUE)[[fun_EvaluateCriterion]]}.
#' An error will be thrown if no such function is found.
#'
#' @param targetCriterionValue a numeric value to compare the returned values from \code{fun_EvaluateCriterion} against.
#' @param doseInterval a numeric vector of length 2 denoting the dosing interval to search within.
#' @param modelFile       Path to a PKPD model to be used for the simulations or a parsed IQRmodel object.
#' @param projectPath     Path to a PK/PD `IQRprojectNLME` folder or a path to a GPF.xlsx file or a GPF object.
#'
#' @param simtime a sorted numeric vector denoting the simulation times (passed to \code{\link{IQRtools::sim_IQRmodel}}).
#' @param Tk0             A named character vector (default: NULL). Example : c(INPUT1 = "Tk0", INPUT2 = "Tk0_2"), meaning that
#' the IQRdosing object for each simulation would be adjusted to have TINF = set to the model parameter Tk0 for ADM = 1.
#' @param nbrDoses an integer (default 1), denoting the number of doses to simulate
#' @param timeInterval a numeric (default 24), denoting the time interval between doses (in the model time units).
#' @param args_IQRdosing a named list (default NULL) with arguments to be passed to \code{\link{IQRtools::IQRdosing}}. For flexibility, this
#' can be used to overwrite some of the otherwise set by default arguments to IQRdosing, such as \code{TIME}, \code{ADM}, \code{ADDL}, \code{II}
#' (advanced use only). Default: \code{list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1, II = timeInterval)}.
#'
#' @param N_Trial         Number of trials to be simulated.
#' @param N_SubjperTrial  Number of subjects per trial.
#' @param covariates      A dataframe with covariate values, and optionally regressors (see also arugment regressorExtra).
#' @param DOSEcovariate   Name of covariates describing dose. Example: \code{c("INPUT1"="DOSELEVEL")}. Default: NULL.
#' @param Fpediatric      Name of covariate describing pediatric dose. Example: \code{"Fpediatric"}. Default: NULL.
#' @param regressorExtra  A character vector indicating the names of regressors to take from the covariate data,
#' instead of sampling them as model parameters. Extra regressors not present in `projectPath` but present in `modelFile`.
#' @param FLAG_SAMPLE     An integer (by default 1) specifying how model parameters are sampled for each virtual trial population.
#' Possible values:
#'
#' |FLAG_SAMPLE | Sample from uncertainty distribution | Sample covariate data and apply formulas | Sample from IIV distribution  | Notes                                        |
#' |:-----------|:------------------------------------:|:----------------------------------------:|:-----------------------------:|:---------------------------------------------|
#' |0           |    No                                | Yes                                      |  Yes                          |Use point estimates as pop. parameter values  |
#' |1 (default) |   Yes                                | Yes                                      |  Yes                          |                                              |
#' |2           |   Yes                                | No                                       |   No                          |Forces N_SubjperTrial = 1                     |
#' |3           |    No                                | Yes                                      |   No                          |Use point estimates as pop. parameter values; |
#' |4           |   Yes                                | Yes                                      |   No                          |                                              |
#'
#'
#' @param addArgs_sim_IQRmodel a named list (default NULL) with additional arguments to be passed to \code{\link{IQRtools::sim_IQRmodel}}.
#' @param addArgs_uniroot a named list (default \code{list(tol = .1)}) with additional arguments to be passed to \code{\link{stats::uniroot}}.
#'
#' @param percentiles     Percentiles to estimate in each trial (Default: \code{c(5,50,95)}).
#' @param CIlevel         Confidence interval to estimate for each percentile (Default: 90),
#' based on the simulated trials
#'
#' @param seed An integer specifying a seed for the random generator (Default: 1234567).
#' @param Nparallel Number of cores to use for parallel simulations (Default: 1).
#' @param FLAG_BROWSER1 A logical (default: \code{FALSE}), indicating whether a browser should be launched for each individual simulation (debug purpose only).
#'
#' @return Data frame with columns Variable, Metric, CI Low, CI Median, CI High, CI Level. The Variable column equals e.g. tMIC and
#' the Metric values correspond to the different values in the percentiles argument. The other columns correspond to the
#' predicted doses for each percentile in Metric.
#'
#' @importFrom dplyr bind_rows select_if mutate_if
#' @importFrom foreach foreach %do% %dopar%
#'
#' @author Venelin Mitov, IntiQuan
#'
#' @export
#' @importFrom MMVbase summaryAcrossTrials
#' @md
predictDose_GenericTrial <- function(
  fun_EvaluateCriterion,
  targetCriterionValue,
  doseInterval,
  modelFile,

  projectPath,

  simtime,
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,
  args_IQRdosing = list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1, II = timeInterval),

  N_Trial,
  N_SubjperTrial,
  covariates     = NULL,
  DOSEcovariate  = NULL,
  Fpediatric     = NULL,
  regressorExtra = character(0),
  FLAG_SAMPLE    = 1,

  percentiles = c(5,50,95),
  CIlevel     = 90,

  addArgs_sim_IQRmodel = NULL,
  addArgs_uniroot = list(tol = 1),

  Nparallel = 1,
  seed = 1234567,
  FLAG_BROWSER1 = FALSE) {

  # . Create IQRmodel object ----
  if(is.character(modelFile)) {
    model <- IQRmodel(modelFile)
  }

  # . Check that all regressorExtra are parameters in model, present in covariates and non-na numeric values ----
  if( !(isTRUE(all(regressorExtra %in% names(model$parameters))) &&
        isTRUE(all(regressorExtra %in% names(covariates))) &&
        isTRUE(all(sapply(regressorExtra, function(rE) isTRUE(all(!is.na(as.numeric(covariates[[rE]]))))))) ) ) {
    stop("predictDose_GenericTrial: all regressorExtra should be parameters in the model, present as columns in covariates and non-NA values.")
  }

  # . Create GPF object ----
  if(is.character(projectPath)) {
    if(endsWith(projectPath, ".xlsx")) {
      # path to a GPF.xlsx file
      gpf <- try(IQRtools::GPF(filename = projectPath), silent = TRUE)
      if(inherits(gpf, "try-error")) {
        stop("predictDose_GenericTrial: error reading GPF file from projectPath ", projectPath)
      }
    } else {
      # projectPath is an IQRnlmerProject directory
      if(!dir.exists(projectPath)) {
        stop("predictDose_GenericTrial: Directory projectPath does not exist: ", projectPath)
      }
      gpf <- IQRtools::generate_GPFFromIQRnlmeProject(projectPath)
    }
  }

  # . Seed the random generator ----
  set.seed(seed)

  # . Sample parameters for the trials ----
  listSampledParamsTrials <- list()
  for(i in seq_len(N_Trial)) {
    if(!is.null(covariates) && ! "ID" %in% names(covariates)) {
      # The ID column will be stored as column ID0 in the member "sampledData" in the list
      # returned by sampleIndParamValues_MMVmalariaXLS
      covariates$ID <- seq_len(nrow(covariates))
    }
    listSampledParamsTrials[[i]] <- sampleIndParamValues_MMVmalariaXLS(spec = gpf,
                                                                       data = covariates,
                                                                       Nsamples = N_SubjperTrial,
                                                                       Npop = 1,
                                                                       FLAG_SAMPLE = FLAG_SAMPLE)
  }

  # . Predict dose for each subject in each trial ----
  `%op%` <- if(Nparallel > 1) `%dopar%` else `%do%`

  predictedDosesTrials <- foreach(i = seq_len(N_Trial), .combine = bind_rows) %op% {
    sampledParams <- listSampledParamsTrials[[i]]

    predictedDosesSubjs <- NULL
    for(j in seq_len(N_SubjperTrial)) {

      if(is.null(DOSEcovariate)) {
        indParams <- listSampledParamsTrials[[i]]$indParamValues[j, , drop = FALSE]
        vecIndParams <- unlist(indParams)
        # remove ID.POP and ID entries.
        parameterNames <- names(vecIndParams)
        parameterNames <- setdiff(parameterNames, c("ID.POP", "ID"))
        vecIndParams <- vecIndParams[parameterNames]

        # Regressor values stored as columns in covariates
        indRegressors <- listSampledParamsTrials[[i]]$sampledData[j, , drop = FALSE]
        # Extract only the columns from indRegressors that are convertible to numeric
        indRegressors <- mutate_if(indRegressors, is.character, as.numeric)
        indRegressors <- select_if(indRegressors, function(x) !is.na(x))
        vecIndRegressors <- unlist(indRegressors)

        # . . The regressors listed in regressorExtra overwrite sampled parameter values ----
        vecIndParams[intersect(regressorExtra, names(vecIndParams))] <- vecIndRegressors[intersect(regressorExtra, names(vecIndParams))]

        # . . Create a vector of the sampled model parameters and the regressors ----
        vecIndParamsRegs <- c(vecIndParams, vecIndRegressors[intersect(names(model$parameters), setdiff(regressorExtra, names(vecIndParams)))])

        predictedDose <- suppressWarnings(predictDose_Generic(fun_EvaluateCriterion = fun_EvaluateCriterion,
                                                              targetCriterionValue  = targetCriterionValue,
                                                              doseInterval          = doseInterval,
                                                              parameters            = vecIndParamsRegs,
                                                              model                 = model,
                                                              simtime               = simtime,
                                                              Tk0                   = Tk0,
                                                              nbrDoses              = nbrDoses,
                                                              timeInterval          = timeInterval,
                                                              args_IQRdosing        = args_IQRdosing,
                                                              addArgs_sim_IQRmodel  = addArgs_sim_IQRmodel,
                                                              addArgs_uniroot       = addArgs_uniroot,
                                                              FLAG_BROWSER          = FLAG_BROWSER1))
      } else {
        listParameters <- list(
          gpf = gpf,
          FLAG_SAMPLE = FLAG_SAMPLE,
          popParamValues = listSampledParamsTrials[[i]]$popParamValues,
          sampledData = listSampledParamsTrials[[i]]$sampledData[j, , drop = FALSE],
          randomEffects = listSampledParamsTrials[[i]]$randomEffects[j, , drop = FALSE],
          DOSEcovariate = DOSEcovariate,
          Fpediatric = Fpediatric,
          regressorExtra = regressorExtra)
        predictedDose <- suppressWarnings(predictDose_Generic(fun_EvaluateCriterion = fun_EvaluateCriterion,
                                                              targetCriterionValue  = targetCriterionValue,
                                                              doseInterval          = doseInterval,
                                                              parameters            = listParameters,
                                                              model                 = model,
                                                              simtime               = simtime,
                                                              Tk0                   = Tk0,
                                                              nbrDoses              = nbrDoses,
                                                              timeInterval          = timeInterval,
                                                              args_IQRdosing        = args_IQRdosing,
                                                              addArgs_sim_IQRmodel  = addArgs_sim_IQRmodel,
                                                              addArgs_uniroot       = addArgs_uniroot,
                                                              FLAG_BROWSER          = FLAG_BROWSER1))
      }
      predictedDosesSubjs <- bind_rows(predictedDosesSubjs, data.frame(IDTrial = i, IDSubj = j, Dose = predictedDose))
    }
    predictedDosesSubjs
  }

  # . Summarize predicted doses ----
  smmByTrial <- summaryByTrial(predictedDosesTrials, "Dose", percentiles = percentiles, usubjidCOL = "IDSubj", trialCOL = "IDTrial")
  smmAcrossTrials <- summaryAcrossTrials(smmByTrial, CIlevel = CIlevel)

  # . Output
  smmAcrossTrials
}

#' Find the minimal dose satisfying a criterion within a given dose interval for a single subject PKPD model simulation
#'
#' @param fun_EvaluateCriterion a function object or a character string denoting the name of a function accessible
#' from the calling (parent) environment. This function should have the following synopsis:
#' \code{function(sim_results, parameters) { ... }}, where \code{sim_results} is a data.frame resulting from a
#' previous call to \code{\link{IQRtools::sim_IQRmodel}} and \code{parameters} is a data.frame with a single row
#' and numeric columns named as the model parameters; the function should return a numeric value. If
#' \code{fun_EvaluateCriterion} is a character string, the function is searched for via
#' \code{mget(fun_EvaluateCriterion, mode = "function", ifnotfound = list(NULL), inherits = TRUE)[[fun_EvaluateCriterion]]}.
#' An error will be thrown if no such function is found.
#'
#' @param targetCriterionValue a numeric value to compare the returned values from \code{fun_EvaluateCriterion} against.
#' @param doseInterval a numeric vector of length 2 denoting the dosing interval to search within.
#' @param model an IQRmodel object denoting the structural model to simulate, i.e. an object returned by a previous call to \code{\link{IQRtools::IQRmodel}}.
#' @param parameters Either a named numeric vector denoting the model parameter values to simulate or, in the case of a dose covariate, a named list
#' with elements as follows:
#'   * popParamValues: a data.frame as the popParamValues member of the list returned by \code{\link{sampleIndParamValues_MMVmalariaXLS}} (only the 1st row will be used).
#'   * sampledData: a data.frame as the sampledData member of the list returned by \code{\link{sampleIndParamValues_MMVmalariaXLS}}.
#'     This data.frame should have an integer column 'ID'. Only the 1st row of this data.frame is going to be used.
#'   * randomEffects: a data.frame as the randomEffects member of the list returned by \code{\link{sampleIndParamValues_MMVmalariaXLS}} (only the 1st row will be used).
#'   * DOSEcovariate: a named character vector (see documentation of \code{\link{predictDose_GenericTrial}}).
#'   * Fpediatric: a character string (see documentation of \code{\link{predictDose_GenericTrial}}).
#'   * regressorExtra: a character vector (see documentation of \code{\link{predictDose_GenericTrial}}).
#'   * gpf: a GFP object (used to calculate typical individual values based on the covariates in sampledData and to calculate individual values based on
#'   randomEffects and transformation type).
#'   * Fpediatric: (optional) a character string denoting a column in sampledData to scale adult equivalent doses, before applying covariate formulae.
#'   * FLAG_SAMPLE: an integer (see documentation of \code{\link{predictDose_GenericTrial}} for possible values)
#'
#' @param simtime a sorted numeric vector denoting the simulation times (passed to \code{\link{IQRtools::sim_IQRmodel}}).
#' @param Tk0             A named character vector (default: NULL). Example : c(INPUT1 = "Tk0", INPUT2 = "Tk0_2"), meaning that
#' the IQRdosing object for each simulation would be adjusted to have TINF = set to the model parameter Tk0 for ADM = 1.
#' @param nbrDoses an integer (default 1), denoting the number of doses to simulate
#' @param timeInterval a numeric (default 24), denoting the time interval between doses (in the model time units).
#' @param args_IQRdosing a named list (default NULL) with arguments to be passed to \code{\link{IQRtools::IQRdosing}}. For flexibility, this
#' can be used to overwrite some of the otherwise set by default arguments to IQRdosing, such as \code{TIME}, \code{ADM}, \code{ADDL}, \code{II}
#' (advanced use only).
#' @param addArgs_sim_IQRmodel a named list (default NULL) with additional arguments to be passed to \code{\link{IQRtools::sim_IQRmodel}}.
#' @param addArgs_uniroot a named list (default \code{list(tol = .1)}) with additional arguments to be passed to \code{\link{stats::uniroot}}.
#' @param FLAG_BROWSER logical (default FALSE) indicating if a browser should be started before each simulation (debugging purpose only).
#' @return A (possibly infinite) numeric value. If the criterion value evaluated at \code{min(doseInteval) >= targetCriterionValue}, \code{min(doseInterval)} is returned;
#' if the criterion value evaluated at \code{max(doseInteval) < targetCriterionValue}, \code{Inf} is returned; otherwise, a value within the interval
#' \code{[min(doseInterval), max(doseInterval)]} is found using the function \code{\link{stats::uniroot}}. Note that the resulting value is undetermined if \code{criterion
#' function - targetCriterionValue} has more than 1 root within the interval \code{[min(doseInterval), max(doseInterval)]}.
#'
#' @details In the current implementation, the search for a dose is done using the \code{\link{stats::uniroot}} and in the log10-transformed dose interval.
#' This implies that the convergence tolerance gets bigger for higher doses.
#' @importFrom IQRtools sim_IQRmodel IQRdosing
#'
#' @author Venelin Mitov, IntiQuan
#'
#' @md
#' @export
predictDose_Generic <- function(fun_EvaluateCriterion,
                                targetCriterionValue,
                                doseInterval,
                                model,
                                parameters,
                                simtime,
                                Tk0 = NULL,
                                nbrDoses             = 1,
                                timeInterval         = 24,
                                args_IQRdosing       = list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1, II = timeInterval),
                                addArgs_sim_IQRmodel = NULL,
                                addArgs_uniroot      = list(tol = .1),
                                FLAG_BROWSER         = FALSE,
                                FLAG_RPROF           = FALSE,
                                FLAG_SAVECALLS_sim_IQRmodel = FALSE) {

  if(FLAG_RPROF) {
    utils::Rprof(filename = "Rprof_predictDose_Generic.out", append = TRUE, interval = 0.01, memory.profiling = TRUE, gc.profiling = TRUE, line.profiling = TRUE)
  }

  if(is.character(fun_EvaluateCriterion)) {
    f_EvaluateCriterion <- mget(fun_EvaluateCriterion, mode = "function", ifnotfound = list(NULL), inherits = TRUE)[[fun_EvaluateCriterion]]
    if(is.function(f_EvaluateCriterion)) {
      fun_EvaluateCriterion <- f_EvaluateCriterion
    } else {
      stop("predictDose_Generic: Cannot find a function ", fun_EvaluateCriterion, " in the parent frame")
    }
  } else if(!is.function(fun_EvaluateCriterion)) {
    stop("predictDose_Generic: The argument fun_EvaluateCriterion should be either a function or a character string denoting the name of a function accessible from the calling environment.")
  }
  if(! (is.list(args_IQRdosing) && isTRUE(all(c("TIME", "ADM", "AMT", "ADDL", "II") %in% names(args_IQRdosing)))) ) {
    stop("predictDose_Generic: The argument args_IQRdosing should be a named list with at least these members: TIME, ADM, AMT, ADDL, II. Please, check the documentation.")
  } else {
    force(args_IQRdosing)
  }


  # an environment to hold a memory for the previous dose simulation during uniroot call.
  envCacheSimulateCriterion                             <- new.env()
  # envCacheSimulateCriterion$model                       <- model
  # envCacheSimulateCriterion$parameters                  <- parameters
  # envCacheSimulateCriterion$Tk0                         <- Tk0
  # envCacheSimulateCriterion$args_IQRdosing              <- args_IQRdosing
  # envCacheSimulateCriterion$addArgs_sim_IQRmodel        <- addArgs_sim_IQRmodel
  # envCacheSimulateCriterion$FLAG_LOG10DOSE              <- TRUE
  # envCacheSimulateCriterion$FLAG_BROWSER                <- FLAG_BROWSER
  # envCacheSimulateCriterion$FLAG_SAVECALLS_sim_IQRmodel <- FLAG_SAVECALLS_sim_IQRmodel

  envCacheSimulateCriterion$prevDose    <- NA_real_
  envCacheSimulateCriterion$prevCritVal <- NA_real_

  # Call uniroot
  x <- try(do.call(uniroot,
                   c(list(simulateCriterion,
                          interval = log10(doseInterval),

                          envCache = envCacheSimulateCriterion,
                          model = model,
                          parameters = parameters,
                          simtime    = simtime,
                          Tk0 = Tk0,
                          args_IQRdosing = args_IQRdosing,
                          addArgs_sim_IQRmodel = addArgs_sim_IQRmodel,
                          fun_EvaluateCriterion = fun_EvaluateCriterion,
                          targetCriterionValue = targetCriterionValue,
                          FLAG_LOG10DOSE = TRUE,
                          FLAG_BROWSER = FLAG_BROWSER,
                          FLAG_SAVECALLS_sim_IQRmodel = FLAG_SAVECALLS_sim_IQRmodel),
                     addArgs_uniroot)), silent = TRUE)

  result <- if(inherits(x, "try-error")) {
    if(grepl('values at end points not of opposite sign', as.character(x))) {
      if(simulateCriterion(min(doseInterval),
                           envCache = envCacheSimulateCriterion,
                           model = model,
                           parameters = parameters,
                           simtime    = simtime,
                           Tk0 = Tk0,
                           args_IQRdosing = args_IQRdosing,
                           addArgs_sim_IQRmodel = addArgs_sim_IQRmodel,
                           fun_EvaluateCriterion = fun_EvaluateCriterion,
                           targetCriterionValue = targetCriterionValue,
                           FLAG_LOG10DOSE = FALSE,
                           FLAG_BROWSER = FLAG_BROWSER,
                           FLAG_SAVECALLS_sim_IQRmodel = FLAG_SAVECALLS_sim_IQRmodel) >= 0.0) {
        min(doseInterval)
      } else if(simulateCriterion(max(doseInterval),
                                  envCache = envCacheSimulateCriterion,
                                  model = model,
                                  parameters = parameters,
                                  simtime    = simtime,
                                  Tk0 = Tk0,
                                  args_IQRdosing = args_IQRdosing,
                                  addArgs_sim_IQRmodel = addArgs_sim_IQRmodel,
                                  fun_EvaluateCriterion = fun_EvaluateCriterion,
                                  targetCriterionValue = targetCriterionValue,
                                  FLAG_LOG10DOSE = FALSE,
                                  FLAG_BROWSER = FLAG_BROWSER,
                                  FLAG_SAVECALLS_sim_IQRmodel = FLAG_SAVECALLS_sim_IQRmodel) < 0.0) {
        Inf
      }
    } else {
      stop("predictDose_Generic: Unexpected error during call to uniroot: ", as.character(x))
    }
  } else {
    if(x$f.root >= 0) {
      10^x$root
    } else {
      envCacheSimulateCriterion$prevDose
    }
  }

  if(FLAG_RPROF) {
    # stop profiling
    utils::Rprof(NULL)
  }

  result
}

# Internal function to call from uniroot in predictDose_Generic
simulateCriterion <- function(dose, envCache, model,
                              parameters,
                              simtime,
                              Tk0,
                              args_IQRdosing,
                              addArgs_sim_IQRmodel,
                              fun_EvaluateCriterion,
                              targetCriterionValue,
                              FLAG_LOG10DOSE = TRUE,
                              FLAG_BROWSER = FALSE,
                              FLAG_SAVECALLS_sim_IQRmodel = FALSE) {

  if(FLAG_BROWSER == TRUE) {
    browser()
  }

  if(FLAG_LOG10DOSE) {
    dose <- 10^dose
  }

  # Set the dose to simulate in args_IQRdosing
  args_IQRdosing$AMT[is.na(args_IQRdosing$AMT)] <- dose

  # If parameters is a list, calculate individual simulation parameters
  if(is.list(parameters)) {
    # Check members of the parameters list
    if(!isTRUE(all(c("popParamValues", "gpf", "FLAG_SAMPLE") %in% names(parameters)))) {
      stop('predictDose_Generic: when the argument "parameters" is a list it should have at least the members "gpf", popParamValues, and "FLAG_SAMPLE" (check documentation).')
    }
    if(!is_GPF(parameters$gpf)) {
      stop('predictDose_Generic: when the argument parameters is a list parameters$gpf should be a GPF object.')
    }
    if(!is.data.frame(parameters$popParamValues)) {
      stop('predictDose_Generic: when the argument parameters is a list, parameters$popParamValeus should be a data.frame (with 1 row).')
    }
    if(!(is.numeric(parameters$FLAG_SAMPLE) && parameters$FLAG_SAMPLE %in% c(0, 1, 2, 3, 4))) {
      stop('predictDose_Generic: when the argument parameters is a list, parameters$FLAG_SAMPLE should be an integer within {0, 1, 2, 3, 4}.')
    }
    if(parameters$FLAG_SAMPLE %in% c(0, 1, 3, 4) && !is.data.frame(parameters$sampledData)) {
      stop('predictDose_Generic: when the argument parameters is a list, parameters$sampledData should be a data.frame (with 1 row) containing covariates.')
    }
    if(parameters$FLAG_SAMPLE %in% c(0, 1) && !is.data.frame(parameters$randomEffects)) {
      stop('predictDose_Generic: when the argument parameters is a list, parameters$randomEffects should be a data.frame (with 1 row).')
    }

    # We may change the covariate data if DOSEcovariate is to be applied, so we create a local scope copy that will be changed.
    sampledData <- parameters$sampledData
    # adjust DOSEcovariate according to args_IQRdosing$AMT
    if(!is.null(parameters$DOSEcovariate)) {
      if(!is.character(parameters$DOSEcovariate)) {
        stop('predictDose_Generic: when the argument parameters$DOSEcovariate is not null, it should be a named character vector.')
      }
      if(!isTRUE(all(names(parameters$DOSEcovariate) %in% paste0("INPUT", args_IQRdosing$ADM)))) {
        stop("predictDose_Generic: the argument parameters$DOSEcovariate should be either NULL or a named character vector with names starting with INPUT (check documentation).")
      }
      for(INPUTx in names(parameters$DOSEcovariate)) {
        x <- as.numeric(gsub("INPUT", "", INPUTx))
        sampledData[[parameters$DOSEcovariate[INPUTx]]] <- unique(args_IQRdosing$AMT[args_IQRdosing$ADM == x])
      }
    }

    # adjust DOSEcovariate according to Fpediatric
    if(!is.null(parameters$Fpediatric)) {
      if(!parameters$Fpediatric %in% names(parameters$sampledData)) {
        stop("predictDose_Generic: the argument parameters$Fpediatric should indicate a column in the covariate data (check documentation and columns in parameters$sampledData).")
      }
      if(!(is.character(parameters$DOSEcovariate) &&
           !is.null(names(parameters$DOSEcovariate)) &&
           isTRUE(all(names(parameters$DOSEcovariate) %in% paste0("INPUT", args_IQRdosing$ADM))))) {
        stop("predictDose_Generic: when parameters$Fpediatric is specified there should be at least one parameters$DOSEcovariate (check documentation).")
      }
      for(INPUTx in names(parameters$DOSEcovariate)) {
        x <- as.numeric(gsub("INPUT", "", INPUTx))
        sampledData[[parameters$DOSEcovariate[INPUTx]]] <- sampledData[[parameters$Fpediatric]] * sampledData[[parameters$DOSEcovariate[INPUTx]]]
      }
    }

    # Apply covariate formulae to population values to obtain typical individual values
    if(parameters$FLAG_SAMPLE %in% c(0, 1, 3, 4)) {
      typicalIndParamValues <- calcTypicalIndParamValues_MMVmalariaXLS(spec = parameters$gpf,
                                                                       referencePopParamValues = parameters$popParamValues[1, , drop = FALSE],
                                                                       data = sampledData[1, , drop = FALSE],
                                                                       doCartesian = FALSE)
    } else {
      typicalIndParamValues <- parameters$typicalIndParamValues
    }

    # remove ID.POP and ID entries.
    parameterNames <- names(typicalIndParamValues)
    parameterNames <- setdiff(parameterNames, c("ID.POP", "ID"))

    if(parameters$FLAG_SAMPLE %in% c(0, 1)) {
      indParams <- calcIndParamValues_MMVmalariaXLS(spec = parameters$gpf,
                                                    typicalIndParamValues = typicalIndParamValues[, parameterNames, drop = FALSE],
                                                    randomEffects = parameters$randomEffects[1, parameterNames, drop = FALSE])
    } else {
      indParams <- typicalIndParamValues
    }

    # combine individual parameter values and regressors into a vector
    vecIndParams <- unlist(indParams)
    vecIndParams <- vecIndParams[parameterNames]

    # Regressor values stored as columns in covariates
    indRegressors <- sampledData
    # Extract only the columns from indRegressors that are convertible to numeric
    indRegressors <- mutate_if(indRegressors, is.character, as.numeric)
    indRegressors <- select_if(indRegressors, function(x) !is.na(x))
    vecIndRegressors <- unlist(indRegressors)

    # The regressors listed in regressorExtra overwrite sampled parameter values
    vecIndParams[intersect(parameters$regressorExtra, names(vecIndParams))] <- vecIndRegressors[intersect(parameters$regressorExtra, names(vecIndParams))]

    # Create a vector of the sampled model parameters and the regressors
    vecIndParamsRegs <- c(vecIndParams, vecIndRegressors[intersect(names(model$parameters), setdiff(parameters$regressorExtra, names(vecIndParams)))])

    # Create local parameters vector that precedes the parameters list argument within this function scope
    parameters <- vecIndParamsRegs
  }

  # Set TINF according to Tk0
  if(!is.null(Tk0)) {
    if(!(is.character(Tk0) && !is.null(names(Tk0)) && isTRUE(all(names(Tk0) %in% paste0("INPUT", args_IQRdosing$ADM))))) {
      stop("predictDose_Generic: the argument Tk0 should be either NULL or a named character vector with names starting with INPUT (check documentation).")
    }
    if(!isTRUE(all(Tk0) %in% names(parameters))) {
      stop("predictDose_Generic: all values in Tk0 should be among the names(parameters).")
    }
    # Set argument TINF to the Tk0 parameter corresponding to ADM
    if(is.null(args_IQRdosing$TINF)) {
      if(!is.null(args_IQRdosing$RATE)) {
        # According to IQRdosing documentation RATE == 0 corresponds to BOLUS.
        args_IQRdosing$TINF <- ifelse(args_IQRdosing$RATE == 0, 0, dose / args_IQRdosing$RATE)
        # According to IQRdosing documentation specifying both RATE and TINF results in an error, so we set RATE to NULL
        args_IQRdosing$RATE <- NULL
      } else {
        # Default TINF value for BOLUS
        args_IQRdosing$TINF <- 0
      }
    }
    # Overwrite TINF values according to Tk0 parameters
    for(INPUTx in names(Tk0)) {
      x <- as.numeric(gsub("INPUT", "", INPUTx))
      args_IQRdosing$TINF[args_IQRdosing$ADM == x] <- parameters[Tk0[INPUTx]]
    }
  }

  # Create the dosing object to update the stored eventTable in envCache from it
  dosingTable  <- do.call(IQRdosing, args_IQRdosing)

  # Use a cached eventTable object to pass as eventTable argument to sim_IQRmodel.
  # This is to prevents expensive eventTable creation during each call to sim_IQRmodel.
  if(is.null(envCache$eventTable)) {
    envCache$eventTable <- IQRtools::IQReventTable(dosingTable)
  } else {
    envCache$eventTable$AMT <- dosingTable$AMT
    envCache$eventTable$TINF <- dosingTable$TINF
  }

  listArgs_sim_IQRmodel <- c(list(model = model, simtime = simtime, eventTable = envCache$eventTable, parameters = unlist(parameters)), addArgs_sim_IQRmodel)

  # Run the simulation
  sim_results <- do.call(sim_IQRmodel, listArgs_sim_IQRmodel)

  if(FLAG_SAVECALLS_sim_IQRmodel) {
    objCall <- list(listArgs = listArgs_sim_IQRmodel)
    objCall$listArgs$model <- NULL
    rdsFile <- "predictDose_Generic_calls_sim_IQRmodel.rds"
    if(file.exists(rdsFile)) {
      listCalls <- readRDS(rdsFile)
    } else {
      listCalls <- NULL
    }
    listCalls <- c(listCalls, list(objCall))
    saveRDS(listCalls, file = rdsFile)
  }

  # Evaluate whether the dose satisfies the criterion.
  # If the distance between this and the previously simulated dose is less than doseTol,
  # return 0.0. This will cause uniroot to stop the search.
  if(!is.null(sim_results)) {
    critVal <- fun_EvaluateCriterion(sim_results, parameters) - targetCriterionValue
    if(critVal >= 0.0) {
      envCache$prevDose    <- dose
      envCache$prevCritVal <- critVal
    }
  } else {
    warning("predictDose_Generic: failed simulation at dose ", dose, ". Returning NA for the criterion value.")
    critVal <- NA_real_
  }

  critVal
}


#' Calculate PRRtot for an individual parasitemia profile
#' @details This function can be passed as argument to \code{\link{predictDose_Generic}}.
#'
#' @param sim_results a data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}. It is assumed that
#' this data.frame has a column PL corresponding to log-transformed simulated parasitemia.
#' @param parameters currently not used.
#' @return a numeric value.
#'
#' @author Venelin Mitov, IntiQuan
#'
#' @export
evaluateDoseCriterion_PRRtot <- function(sim_results, parameters) {
  get_PRRtot(sim_results,paraCOL = "PL", Plog = TRUE)[["PRRtot"]]
}


#' Calculate the time above the minimum inhibitory concentration (MIC) for an individual parasitemia profile
#' @details This function can be passed as argument to \code{\link{predictDose_Generic}}.
#'
#' @param sim_results a data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}. It is assumed that
#' this data.frame has a column PL corresponding to log-transformed simulated parasitemia.
#' @param parameters a named numeric vector.
#' @return a numeric value.
#'
#' @seealso \code{\link{getKeysEMAX}}
#'
#' @author Venelin Mitov, IntiQuan
#'
#' @export
evaluateDoseCriterion_TimeAboveMIC <- function(sim_results, parameters) {
  dfParam <- as.data.frame(as.list(parameters))
  MIC     <- getKeysEMAX(x = dfParam)[["MIC"]]
  getTimeAboveMIC(sim_results, MIC = MIC)[["tMIC"]]
}


#' Calculate the time above the minimum parasiticidal concentration (MPC90) for an individual parasitemia profile
#' @details This function can be passed as argument to \code{\link{predictDose_Generic}}.
#'
#' @param sim_results a data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}. It is assumed that
#' this data.frame has a column PL corresponding to log-transformed simulated parasitemia.
#' @param parameters a named numeric vector.
#' @return a numeric value.
#'
#' @seealso \code{\link{getKeysEMAX}}
#'
#' @author Mohammed H. Cherkaoui, MMV
#'
#' @export
evaluateDoseCriterion_TimeAboveMPC90 <- function(sim_results, parameters) {
  dfParam <- as.data.frame(as.list(parameters))
  MPC90     <- getKeysEMAX(x = dfParam)[["MPC90"]]
  getTimeAboveMPC90(sim_results, MPC90 = MPC90)[["tMPC90"]]
}


