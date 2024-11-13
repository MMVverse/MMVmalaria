#' Summarize given outcome variable(s) for a single trial
#'
#' Summarize given outcome variable(s) for a single trial using a choice of metrics
#' (mean, percentile or geomean). If no metric is provided, percentile is used by default
#' for numeric variable, while mean is used for logical variable.
#'
#' @param dataInput data.table containing outcome variable(s) and selected other variables for each individual.
#' @param varCOL List containing metrics (percentiles, mean and/or geomean) to be calculated for each outcome variable - i.e. \code{list(Cmax=c("Percentile","Mean"),APR="Mean")}.
#' @param percentiles Percentiles of outcome variable(s) to be calculated if appropriate (Default: \code{c(5,50,90)}).
#' @param usubJIDCOL Name(s) of column(s) of `dataInput` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL Name(s) of column(s) of `dataInput` uniquely describing trial (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")}).
#'
#' @return data table containing values of outcome metrics for each trial
#'
#' @export
#' @author Sam Jones (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summaryByTrial <- function(dataInput,
                           varCOL,
                           percentiles = c(5,50,95),
                           usubjidCOL  = c("IndivID","USUBJID"),
                           trialCOL    = c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")){

  #---------------------------------------------------#
  # STEP 1: Run some checks and initialization ----
  #---------------------------------------------------#

  # Which metrics are accepted for calculation by this code:
  #   NOTE: Adding additional metrics will require additional coding
  accepted_metrics <- c("mean","geomean","percentiles")

  # Add/Adjust name to percentiles (Predictive percentiles):
  percentiles <- aux_addNamesToPercentiles(percentiles)

  # Generate colKey to use:
  colKey <- intersect(c(trialCOL,usubjidCOL),
                      names(dataInput))

  # Check if dataInput is a data.table:
  if(!data.table::is.data.table(dataInput)){
    data.table::setDT(dataInput, key = colKey)
  }

  # VarCOL checks:
  #   If varCOL is not provided
  if(is.null(varCOL)) {
    stop("varCOL is not provided - it must be provided in as a list with outcome measures as list names i.e: list(PRR=c('Percentiles','Mean')), please adjust")
  }
  #   If varCOL is provided but without correct structure, generate error
  if(is.null(names(varCOL))){
    varCOL <- sapply(varCOL, function(x){
      out <- "percentiles"
      names(out) <- x
      out
    },USE.NAMES = FALSE)
    # stop("varCOL must be provided as a list with outcome measures as list names i.e: list(PRR=c('Percentiles','Mean')), please adjust")
  }
  # If varCOL is provided, but not in the data:
  #   NOTE: This is calling names as varCOL is a list...
  idx_ToKeep <- (names(varCOL) %in% names(dataInput))
  if(!any(idx_ToKeep)){
    stop("None of the variables in 'varCOL' are present in dataInput: Please adjust.")
  }

  # Force varCOL to a list:
  varCOL <- as.list(varCOL[idx_ToKeep])

  # Create measure.vars:
  measure.vars <- names(varCOL)

  # Make sure that only mean is being used for logical variables:
  class.vars   <- sapply(measure.vars, function(x){class(dataInput[[x]])})
  logical.vars <- measure.vars[class.vars=="logical"]
  if(any(tolower(varCOL[logical.vars])!="mean")){
    varCOL[logical.vars] <- "mean"
    warning("Metrics other than arithmetic mean cannot be calculated for logical outcome variables: It was automatically adjusted.")
  }

  # Ensure that only mean, geomean and percentiles are chosen in varCOL
  metrics_diff <- lapply(varCOL, function(x){
    out <- setdiff(tolower(x), accepted_metrics)
    out
  })
  if(any(lapply(metrics_diff,length)>0)){
    stop("Metric provided in varCOL that cannot be calculated. Accepted metrics are 'mean','geomean','percentiles' ")
  }


  #---------------------------------------------------#
  # STEP 2: Identify what to apply for each variable in varCOL ----
  #---------------------------------------------------#
  idx_mean        <- sapply(varCOL, function(x){"mean"        %in% tolower(x)})
  idx_geomean     <- sapply(varCOL, function(x){"geomean"     %in% tolower(x)})
  idx_percentiles <- sapply(varCOL, function(x){"percentiles" %in% tolower(x)})


  #---------------------------------------------------#
  # STEP 3: Produce long format data tables to operate on ----
  #---------------------------------------------------#

  # Add variables and values to data-frame if tagged to calculate mean:
  if(any(idx_mean)){
    data.byIndiv.mean <- data.table::melt.data.table(dataInput,
                                                     id.vars       = colKey,
                                                     measure.vars  = measure.vars[idx_mean],
                                                     variable.name = "Variable",
                                                     value.name    = "Value")
  }else{
    data.byIndiv.mean <- NULL
  }
  # Add variables and values to data-frame if tagged to calculate geometric mean
  if(any(idx_geomean)){
    data.byIndiv.geomean <- data.table::melt.data.table(dataInput,
                                                        id.vars       = colKey,
                                                        measure.vars  = measure.vars[idx_geomean],
                                                        variable.name = "Variable",
                                                        value.name    = "Value")
  }else{
    data.byIndiv.geomean <- NULL
  }
  # Add variables and values to data-frame if tagged to calculate percentiles
  if(any(idx_percentiles)){
    data.byIndiv.percentiles <- data.table::melt.data.table(dataInput,
                                                            id.vars       = colKey,
                                                            measure.vars  = measure.vars[idx_percentiles],
                                                            variable.name = "Variable",
                                                            value.name    = "Value")
  }else{
    data.byIndiv.percentiles <- NULL
  }


  #---------------------------------------------------#
  # STEP 4: Operate on data tables to obtain outcome metrics described in varCOL ----
  #---------------------------------------------------#

  # Query: What's better here - re-using if(any(idx...)) or using if(length(data.by.Indiv... >0)) ?
  # Either way the intent is the same - process calculations for specific metrics if idx_ objects are TRUE for any of varCOL.

  # Mean:
  if(any(idx_mean)){
    colKey <- intersect(trialCOL,
                        names(data.byIndiv.mean))
    data.byTrial.mean <- data.byIndiv.mean[,{
      list(Metric = "Mean",
           Value  = mean(Value,na.rm=TRUE))
    },
    by = c(colKey,"Variable")]
  }else{
    data.byTrial.mean<-NULL
  }

  # Geometric Mean
  if(any(idx_geomean)){
    colKey <- intersect(trialCOL,
                        names(data.byIndiv.geomean))
    data.byTrial.geomean <- data.byIndiv.geomean[,{
      list(Metric = "Geomean",
           Value  = geomean(Value,na.rm=TRUE))
    },
    by = c(colKey,"Variable")]
  }else{
    data.byTrial.geomean<-NULL
  }

  # Percentiles
  if(any(idx_percentiles)){
    colKey <- intersect(trialCOL,
                        names(data.byIndiv.percentiles))
    data.byTrial.percentiles <- data.byIndiv.percentiles[,{
      list(Metric =  names(percentiles),
           Value  =  quantile(Value, probs = percentiles/100, na.rm = TRUE))
    },
    by = c(colKey,"Variable")]
  }else{
    data.byTrial.percentiles<-NULL
  }

  # Combine data tables of all metrics
  data.byTrial <- data.table::rbindlist(list(data.byTrial.mean,
                                             data.byTrial.geomean,
                                             data.byTrial.percentiles))

  # Output:
  data.byTrial
}


#' Summarize confidence intervals for provided metrics across trials.
#'
#' @param dataInput data.table containing summarized outcome variable(s) summarized by trials - produced by `summaryByTrial` or similar format.
#' @param CIlevel Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param metricCOL name of column of `dataInput` containing names of metrics on which CI is being estimated (Default: \code{"Metric"})
#' @param variableCOL name of column of `dataInput` containing names of variables (Default: \code{"Variable"})
#' @param valueCOL name of column of `dataInput` containing values of metrics for specific variables (Default: \code{"Value"})
#' @param scenCOL Name(s) of column(s) of `dataInput` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses","TIME")})
#'
#' @return data table containing confidence intervals for selected outcome variable(s) and summarized metric(s)
#'
#' @export
#' @author Sam Jones (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summaryAcrossTrials <- function(dataInput,
                                CIlevel     = 90,
                                metricCOL   = "Metric",
                                variableCOL = "Variable",
                                valueCOL    = "Value",
                                scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses","TIME")){

  #---------------------------------------------------#
  # STEP 1: Run some checks and initialization ----
  #---------------------------------------------------#

  # Define quantile of CI:
  CI.percentile <- aux_constructCIpercentiles(CIlevel)

  # Check if data is data.table:
  if(!data.table::is.data.table(dataInput)){
    colKey <- intersect(c(scenCOL,variableCOL,metricCOL,valueCOL),
                        names(dataInput))
    data.table::setDT(dataInput, key = colKey)
  }


  #---------------------------------------------------#
  # STEP 2: Calculate CI for metrics of variables ----
  #---------------------------------------------------#

  # Get CI:
  colKey <- intersect(c(scenCOL,variableCOL,metricCOL),
                      names(dataInput))
  data <- dataInput[,{
    Variable.Q        <- quantile(get(valueCOL), probs = CI.percentile/100, na.rm = TRUE)
    Variable.Q        <- as.list(Variable.Q)
    names(Variable.Q) <- names(CI.percentile)
    Variable.Q
  },
  by = colKey]

  # Add CI level:
  data$`CI Level` <- CIlevel

  # Output:
  data
}


#' Generic function to estimate PI and CI
#'
#' Generic function to estimate the confidence interval of a metric (mean, percentiles, and/or geoman)
#' of various given variables, when multiple trials have been simulated. It calls `summaryByTrial`
#' and `summmaryAcrossTrials` in a row. If no metric is provided, percentile is used by default
#' for numeric variable, while mean is used for logical variable.
#'
#' @param dataInput data.table containing outcome variable(s) and selected other variables for each individual.
#' @param varCOL List containing metrics (percentiles, mean and/or geomean) to be calculated for each outcome variable - i.e. \code{list(Cmax=c("Percentile","Mean"),APR="Mean")}.
#' @param percentiles Percentiles of outcome variable(s) to be calculated if appropriate (Default: \code{c(5,50,90)}).
#' @param CIlevel Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param usubJIDCOL Name(s) of column(s) of `dataInput` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL Name(s) of column(s) of `dataInput` uniquely describing trial (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")}).
#' @param scenCOL Name(s) of column(s) of `dataInput` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses","TIME")})
#'
#' @return
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summarize_PIandCIgeneric <- function(dataInput,
                                     varCOL,
                                     percentiles = c(5,50,90),
                                     CIlevel     = 90,
                                     usubjidCOL  = c("IndivID","USUBJID"),
                                     trialCOL    = c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME"),
                                     scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses","TIME")){

  # Summarize 'dataInput' by Trial:
  data.byTrial <- summaryByTrial(dataInput   = dataInput,
                                 percentiles = percentiles,
                                 varCOL      = varCOL,
                                 usubjidCOL  = usubjidCOL,
                                 trialCOL    = unique(c(scenCOL,trialCOL)))

  # Summarize 'data.byTrial' Across Trials:
  data <- summaryAcrossTrials(dataInput = data.byTrial,
                              CIlevel   = CIlevel,
                              scenCOL   = scenCOL)

  # Output:
  data
}


#' Summarize PK and PD profiles across trials
#'
#' Summarize PK and PD profiles entered in `outputCOL`, where the predictive interval
#' (defined in with the variable 'percentiles') is returned with the confidence
#' interval as given by `CIlevel`.
#'
#' @param simPKPD     Data table output of `simulate_virtualtrials()` (or similar format).
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param outputCOL   Vector containing names of desired columns for output (Default: \code{c("Cc", "PL")}).
#' @param percentiles Percentiles of outcome variable(s) to be calculated if appropriate (Default: \code{c(5,50,90)}).
#' @param CIlevel Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param usubjidCOL Name(s) of column(s) of `simPKPD` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL Name(s) of column(s) of `simPKPD` uniquely describing trial (Default: \code{c("TrialID")}).
#' @param scenCOL Name(s) of column(s) of `simPKPD` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses")})
#'
#'
#' @return
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summarize_PKPDprofilesFromSimulations <- function(simPKPD,
                                                  timeCOL     = "TIME",
                                                  outputCOL   = c("Cc", "PL"),
                                                  percentiles = c(5,50,95),
                                                  CIlevel     = 90,
                                                  usubjidCOL  = c("IndivID","USUBJID"),
                                                  trialCOL    = c("TrialID"),
                                                  scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses")){

  # Simply call summarize_PIandCIgeneric:
  summaryPKPD <- summarize_PIandCIgeneric(dataInput   = simPKPD,
                                          varCOL      = outputCOL,
                                          percentiles = percentiles,
                                          CIlevel     = CIlevel,
                                          usubjidCOL  = usubjidCOL,
                                          trialCOL    = unique(c(scenCOL,trialCOL,timeCOL)),
                                          scenCOL     = unique(c(scenCOL,timeCOL)))

  # Output:
  summaryPKPD
}


#' Summarize a selection of Key PK parameters from simulated data
#'
#' \code{"AUCinf"},\code{"Cmax"} and \code{"Tmax"} are estimated for each individual and summarized.
#' The predictive interval (defined in with the variable 'percentiles') is returned with the confidence
#' interval as given by `CIlevel`.
#'
#' @param simPKPD     Data-table output of `simulate_virtualtrials()` (or similar format).
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param concCOL     Name of column containing values for drug concentration (Default: \code{"Cc"}).
#' @param percentiles Percentiles of outcome variable(s) to be calculated if appropriate (Default: \code{c(5,50,90)}).
#' @param CIlevel     Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param usubjidCOL  Name(s) of column(s) of `simPKPD` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL    Name(s) of column(s) of `simPKPD` uniquely describing trial (Default: \code{c("TrialID")}).
#' @param scenCOL     Name(s) of column(s) of `simPKPD` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses")})
#'
#' @return            Data-table containing values of population PK parameters from simulated trials at specified percentiles
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summarize_KeyPKparametersFromSimulations <- function(simPKPD,
                                                     timeCOL     = "TIME",
                                                     concCOL     = "Cc",
                                                     percentiles = c(5,50,95),
                                                     CIlevel     = 90,
                                                     usubjidCOL  = c("IndivID","USUBJID"),
                                                     trialCOL    = c("TrialID"),
                                                     scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses")){


  #---------------------------------------------------#
  # STEP 1: Estimate the key PK parameters for each individuals----
  #---------------------------------------------------#

  # Check if simPKPD is data.table:
  if(!data.table::is.data.table(simPKPD)){
    colKey <- intersect(unique(c(scenCOL,trialCOL,usubjidCOL,timeCOL)),
                        names(simPKPD))
    data.table::setDT(simPKPD, key = colKey)
  }

  # Get Key PK parameters for each individual:
  colKey <- intersect(unique(c(scenCOL,trialCOL,usubjidCOL)),
                      names(simPKPD))
  summaryKeyPK.byIndiv <- simPKPD[,{
    list(AUCinf = estimate_AUC(Time    = get(timeCOL),
                               Conc    = get(concCOL),
                               AUCtype = "inf"),
         Cmax   = max(get(concCOL)),
         Tmax   = get(timeCOL)[which.max(get(concCOL))])
  },
  by = colKey]


  #---------------------------------------------------#
  # STEP 2: Summarize key PK parameters for each ScenID ----
  #---------------------------------------------------#

  # Simply call summarize_PIandCIgeneric:
  summaryKeyPK <- summarize_PIandCIgeneric(summaryKeyPK.byIndiv,
                                           varCOL = setdiff(names(summaryKeyPK.byIndiv),colKey),
                                           percentiles = percentiles,
                                           CIlevel     = CIlevel,
                                           usubjidCOL  = usubjidCOL,
                                           trialCOL    = trialCOL,
                                           scenCOL     = scenCOL)

  # Output:
  summaryKeyPK
}


#' Summarize a selection of Key PD parameters from simulated data
#'
#' \code{"MIC"},\code{"MPC90"}, \code{"PRR24"}, \code{"PRR48"}, \code{"PRR72"}, \code{"PRRtot"},
#' \code{"tMIC"} and \code{"tMPC90"} are estimated for each individual and summarized.
#' The predictive interval (defined in with the variable 'percentiles') is returned with the confidence
#' interval as given by `CIlevel`.
#'
#' @param simPKPD     Data-table output of `simulate_virtualtrials()` (or similar format).
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param concCOL     Name of column containing values for drug concentration (Default: \code{"Cc"}).
#' @param paraCOL     Name of column containing values for parasitaemia (Default:  \code{"PL"}).
#' @param Plog        Is paraistaemia provided in Log scale?  (Default: `TRUE`).
#' @param percentiles Percentiles of outcome variable(s) to be calculated if appropriate (Default: \code{c(5,50,90)}).
#' @param CIlevel     Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param usubjidCOL  Name(s) of column(s) of `simPKPD` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL    Name(s) of column(s) of `simPKPD` uniquely describing trial (Default: \code{c("TrialID")}).
#' @param scenCOL     Name(s) of column(s) of `simPKPD` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses")})
#'
#' @return            Data-table containing values of population PK parameters from simulated trials at specified percentiles
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summarize_KeyPDparametersFromSimulations <- function(simPKPD,
                                                     timeCOL     = "TIME",
                                                     concCOL     = "Cc",
                                                     paraCOL     = "PL",
                                                     Plog        = TRUE,
                                                     percentiles = c(5,50,95),
                                                     CIlevel     = 90,
                                                     usubjidCOL  = c("IndivID","USUBJID"),
                                                     trialCOL    = c("TrialID"),
                                                     scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses")){


  #---------------------------------------------------#
  # STEP 1: Estimate the key PD parameters for each individuals----
  #---------------------------------------------------#

  # Check if simPKPD is data.table:
  if(!data.table::is.data.table(simPKPD)){
    colKey <- intersect(unique(c(scenCOL,trialCOL,usubjidCOL,timeCOL)),
                        names(simPKPD))
    data.table::setDT(simPKPD, key = colKey)
  }

  # Get Key PD parameters for each individual:
  colKey <- intersect(unique(c(scenCOL,trialCOL,usubjidCOL)),
                      names(simPKPD))
  summaryKeyPD.byIndiv <- simPKPD[,{

    # Estimate the static key PD parameter using only the first row
    KeyPD <- getKeysEMAX(.SD[1,])

    # Output:
    list(MIC    = KeyPD[["MIC"]],
         MPC90  = KeyPD[["MPC90"]],
         PRR24  = KeyPD[["PRR24"]],
         PRR48  = KeyPD[["PRR48"]],
         PRR72  = KeyPD[["PRR72"]],
         PRRtot = get_PRRtot(.SD, paraCOL = paraCOL, Plog = Plog)[["PRRtot"]],
         tMIC   = getTimeAboveMIC(.SD,
                                  MIC      = KeyPD[["MIC"]],
                                  timeCOL  = timeCOL,
                                  # Venelin: the following code causes an error:
                                  # Error in `[.data.table`(simPKPD, , { :
                                  # Column 7 of j's result for the first group is NULL. We rely on the column types of the first result to decide the type expected for the remaining groups (and require consistency). NULL columns are acceptable for later groups (and those are replaced with NA of appropriate type and recycled) but not for the first. Please use a typed empty vector instead, such as integer() or numeric().
                                  concCOL  = concCOL)[["tMIC"]],
         tMPC90 = getTimeAboveMPC90(.SD,
                                    MPC90      = KeyPD[["MPC90"]],
                                  timeCOL  = timeCOL,
                                  # Error in `[.data.table`(simPKPD, , { :
                                  # Column 8 of j's result for the first group is NULL. We rely on the column types of the first result to decide the type expected for the remaining groups (and require consistency). NULL columns are acceptable for later groups (and those are replaced with NA of appropriate type and recycled) but not for the first. Please use a typed empty vector instead, such as integer() or numeric().
                                  concCOL  = concCOL)[["tMPC90"]])
  },
  by = colKey]


  #---------------------------------------------------#
  # STEP 2: Summarize key PD parameters for each ScenID ----
  #---------------------------------------------------#

  # Simply call summarize_PIandCIgeneric:
  summaryKeyPD <- summarize_PIandCIgeneric(summaryKeyPD.byIndiv,
                                           varCOL = setdiff(names(summaryKeyPD.byIndiv),colKey),
                                           percentiles = percentiles,
                                           CIlevel     = CIlevel,
                                           usubjidCOL  = usubjidCOL,
                                           trialCOL    = trialCOL,
                                           scenCOL     = scenCOL)

  # Output:
  summaryKeyPD
}


#' Summarizes the total parasite reduction ratio
#'
#' @param simPKPD     Data-table output of `simulate_virtualtrials()` (or similar format).
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param paraCOL     Name of column containing values for parasitaemia (Default:  \code{"PL"}).
#' @param Plog        Is paraistaemia provided in Log scale?  (Default: `TRUE`).
#' @param percentiles Percentiles of outcome variable(s) to be calculated if appropriate (Default: \code{c(5,50,90)}).
#' @param CIlevel     Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param usubjidCOL  Name(s) of column(s) of `simPKPD` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL    Name(s) of column(s) of `simPKPD` uniquely describing trial (Default: \code{c("TrialID")}).
#' @param scenCOL     Name(s) of column(s) of `simPKPD` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses")})
#'
#'
#' @return
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Summarize Functions
summarize_PRRtotFromSimulations <- function(simPKPD,
                                            timeCOL     = "TIME",
                                            paraCOL     = "PL",
                                            Plog        = TRUE,
                                            percentiles = c(5,50,95),
                                            CIlevel     = 90,
                                            usubjidCOL  = c("IndivID","USUBJID"),
                                            trialCOL    = c("TrialID"),
                                            scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses")){

  #---------------------------------------------------#
  # STEP 1: Estimate PRRtot for each individuals----
  #---------------------------------------------------#

  # Check if simPKPD is data.table:
  if(!data.table::is.data.table(simPKPD)){
    colKey <- intersect(unique(c(scenCOL,trialCOL,usubjidCOL,timeCOL)),
                        names(simPKPD))
    data.table::setDT(simPKPD, key = colKey)
  }

  # Get PRRtot for each individual:
  colKey <- intersect(unique(c(scenCOL,trialCOL,usubjidCOL)),
                      names(simPKPD))
  summaryPRRtot.byIndiv <- simPKPD[,{
    list(PRRtot = get_PRRtot(.SD, paraCOL = paraCOL, Plog = Plog)[["PRRtot"]])
  },
  by = colKey]


  #---------------------------------------------------#
  # STEP 2: Summarize PRRtot for each ScenID ----
  #---------------------------------------------------#

  # Simply call summarize_PIandCIgeneric:
  summaryPRRtot <- summarize_PIandCIgeneric(summaryPRRtot.byIndiv,
                                            varCOL = c("PRRtot"),
                                            percentiles = percentiles,
                                            CIlevel     = CIlevel,
                                            usubjidCOL  = usubjidCOL,
                                            trialCOL    = trialCOL,
                                            scenCOL     = scenCOL)

  # Output:
  summaryPRRtot
}

#' Summarize PKPD simulation results for a virtual trial of one experimental setting and one dose
#'
#' This function is called at runtime by \code{\link{simulate_OneVirtualTrial}}.
#' @md
#'
#' @param simPKPD A data.frame representing the merge (join) for data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}, and
#' individual and experimental covariates. This argument is set automatically by the calling function.
#' @param trialFileName A character string denoting the trial .rds file on the disk. This file can be loaded to get information about the
#' virtual subjects, such as individual and experimental covariates.
#' @param varPKPD a character vector indicating time varying PKPD variables of the simulated PKPD model. All of these
#' should be names of columns in \code{simPKPD}.
#' @param percentiles a numeric vector indicating percentiles (default: \code{c(5, 50, 95)}) for summarizing
#' individual PKPD variables and clinical end-points.
#' @details This function is designed to be passed as argument of \code{\link{simulate_VirtualTrials}}. Currently, the function
#' assumes that the \code{simPKPD} has a column PL denoting natural log parasitemia.
#'
#' @return a named list with the following members:
#'
#' * ScenID: an integer
#' * ExpID: an integer
#' * DoseID: an integer
#' * TrialID: an integer
#' * summaryPKPD.ByTrial: a data.frame
#' * summaryClinEnd.ByTrial: NULL
#' * summaryClinEnd.ByIndiv: NULL
#'
#' @export
summarizeTrial_Basic <- function(
    simPKPD,
    trialFileName,
    varPKPD             = c("Cc","PL"),
    percentiles         = c(5, 50, 95)) {

  ScenID  <- unique(simPKPD$ScenID)
  ExpID   <- unique(simPKPD$ExpID)
  DoseID  <- unique(simPKPD$DoseID)
  TrialID <- unique(simPKPD$TrialID)

  stopifnot(
    "summarizeTrial_Basic: function must be called on simPKPD for a single ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)
  stopifnot("summarizeTrial_Basic: All variables in varPKPD should be column names in simPKPD." = is.character(varPKPD) && isTRUE(all(varPKPD %in% names(simPKPD))))
  stopifnot('summarizeTrial_Basic: simPKPD should have at least the columns: ID, ScenID,ExpID,DoseID,TrialID,IndivID,USUBJID' =
              isTRUE(all(c("ID", "ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID") %in% names(simPKPD))))
  stopifnot("summarizeTrial_Basic: simPKPD must be a data.frame and must have a column ID" = is.data.frame(simPKPD) && "ID" %in% names(simPKPD))

  # Get PKPD summary:
  if(length(unique(simPKPD$ID))>1){
    summaryPKPD.ByTrial <- summaryByTrial(dataInput   = simPKPD,
                                          percentiles = percentiles,
                                          varCOL      = varPKPD,
                                          usubjidCOL  = c("ID","IndivID","USUBJID"),
                                          trialCOL    = unique(c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")))
  }else{
    summaryPKPD.ByTrial <- simPKPD
  }

  # Produce output
  result <- list(summaryPKPD.ByTrial    = summaryPKPD.ByTrial,
                 summaryClinEnd.ByIndiv = NULL,
                 summaryClinEnd.ByTrial = NULL)

  # Output:
  return(result)
}


#' Summarize PKPD simulation results for a virtual trial of one experimental setting and one dose
#'
#' This function is called at runtime by \code{\link{simulate_OneVirtualTrial}}.
#' @md
#'
#' @param simPKPD A data.frame representing the merge (join) for data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}, and
#' individual and experimental covariates. This argument is set automatically by the calling function.
#' @param trialFileName A character string denoting the trial .rds file on the disk. This file can be loaded to get information about the
#' virtual subjects, such as individual and experimental covariates.
#' @param varPKPD a character vector indicating time varying PKPD variables of the simulated PKPD model. All of these
#' should be names of columns in \code{simPKPD}.
#' @param percentiles a numeric vector indicating percentiles (default: \code{c(5, 50, 95)}) for summarizing
#' individual PKPD variables and clinical end-points.
#' @param samplePDtime a numeric vector (default \code{seq(0,28*24,24)}) of time-points at which PL is
#' needed for summarizing individual clinical endpoints. This should contain aprTime.
#' @param LLOQ.PD a numeric indicating log parasitemia LLOQ (default: \code{log(10000)}).
#' @param aprTime a numeric indicating time after beginning of treatment in days at which APR should be evaluated. Default: 28.
#' This argument is passed to function \code{\link{evaluate_APR}}.
#' @param extrapTol,FLAGextrapolateTime These arguments are passed to function \code{\link{evaluate_APR}}
#'
#' @details This function is designed to be passed as argument of \code{\link{simulate_VirtualTrials}}. Currently, the function
#' assumes that the \code{simPKPD} has a column PL denoting natural log parasitemia.
#'
#' @return a named list with the following members:
#'
#' * ScenID: an integer
#' * ExpID: an integer
#' * DoseID: an integer
#' * TrialID: an integer
#' * summaryPKPD.ByTrial: a data.frame
#' * summaryClinEnd.ByTrial: a data.frame of 1 row
#' * summaryClinEnd.ByIndiv: a data.frame of as many rows as there are virtual subjects per trial
#'
#' @export
summarizeTrial_APRLongitudinal <- function(
    simPKPD,
    trialFileName,
    varPKPD             = c("Ccx1","Ccx2","Effx1","Effx2","Kkillx1","Kkillx2","Kkill","PL"),
    percentiles         = c(5, 50, 95),
    samplePDtime        = seq(0,28*24,24),
    LLOQ.PD             = log(10000),
    aprTime             = 28,
    extrapTol           = 48,
    FLAGextrapolateTime = TRUE) {

  ScenID  <- unique(simPKPD$ScenID)
  ExpID   <- unique(simPKPD$ExpID)
  DoseID  <- unique(simPKPD$DoseID)
  TrialID <- unique(simPKPD$TrialID)

  stopifnot(
    "summarizeTrial_APRLongitudinal: function must be called on simPKPD for a single ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)
  stopifnot("summarizeTrial_APRLongitudinal: All variables in varPKPD should be column names in simPKPD." = is.character(varPKPD) && isTRUE(all(varPKPD %in% names(simPKPD))))
  stopifnot('summarizeTrial_APRLongitudinal: simPKPD should have at least the columns: ID, ScenID,ExpID,DoseID,TrialID,IndivID,USUBJID,WT0' =
              isTRUE(all(c("ID", "ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID","WT0") %in% names(simPKPD))))
  stopifnot("summarizeTrial_APRLongitudinal: simPKPD must be a data.frame and must have a column ID" = is.data.frame(simPKPD) && "ID" %in% names(simPKPD))
  stopifnot("summarizeTrial_APRLongitudinal: aprTime*24 should be in samplePDtime." = (aprTime*24) %in% samplePDtime)

  # Get PKPD summary:
  if(length(unique(simPKPD$ID))>1){
    summaryPKPD.ByTrial <- summaryByTrial(dataInput   = simPKPD,
                                          percentiles = percentiles,
                                          varCOL      = varPKPD,
                                          usubjidCOL  = c("ID","IndivID","USUBJID"),
                                          trialCOL    = unique(c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")))
  }else{
    summaryPKPD.ByTrial <- simPKPD
  }

  # Add curethreshold and PLcure:
  simPKPD$CureThreshold <- ifelse(simPKPD$WT0 <= 35,
                                  log(1/(simPKPD$WT0*80)),
                                  log(1/(simPKPD$WT0*70)))
  simPKPD <- simPKPD[,PLcure:={
    Tcure <- NA_real_
    if(any(PL <= CureThreshold )){
      Tcure <- min(TIME[PL < CureThreshold[1]])
    }
    PLcure <- PL
    PLcure[TIME>=Tcure] <- CureThreshold[1]
    PLcure
  }, by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Get individual clinical endpoint:
  summaryClinEnd.ByIndiv <- simPKPD[,{
    interpResult_i <- approx(x = .SD$TIME, y = .SD$PLcure, xout = samplePDtime, ties = list("ordered", mean), rule = 2)

    dat_i          <- data.frame(USUBJID          = USUBJID[1],
                                 TIME             = interpResult_i$x,
                                 PL               = interpResult_i$y,
                                 WT0              = WT0[1],
                                 CureThreshold    = CureThreshold[1],
                                 stringsAsFactors = FALSE)
    evaluate_APR(dat_i,
                 LLOQ                = LLOQ.PD,
                 timeCOL             = "TIME",
                 paraCOL             = "PL",
                 weightCOL           = "WT0",
                 curethreshCOL       = "CureThreshold",
                 aprTime             = aprTime,
                 extrapTol           = extrapTol,
                 FLAGextrapolateTime = FLAGextrapolateTime,
                 Plog                = TRUE,
                 FLAGverbose         = FALSE)

  },
  by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Summarize 'dataInput' by Trial:
  if(length(unique(summaryClinEnd.ByIndiv$ID))>1){
    summaryClinEnd.ByTrial <- summaryByTrial(dataInput   = summaryClinEnd.ByIndiv,
                                             percentiles = percentiles,
                                             varCOL      = setdiff(names(summaryClinEnd.ByIndiv),
                                                                   c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")),
                                             usubjidCOL  = c("ID","IndivID","USUBJID"),
                                             trialCOL    = c("ScenID", "ExpID", "DoseID", "TrialID"))
  }else{
    summaryClinEnd.ByTrial <- summaryClinEnd.ByIndiv
  }

  # Produce output
  result <- list(
    summaryPKPD.ByTrial    = summaryPKPD.ByTrial,
    summaryClinEnd.ByIndiv = summaryClinEnd.ByIndiv,
    summaryClinEnd.ByTrial = summaryClinEnd.ByTrial)

  # Output:
  return(result)
}



#' Summarize PKPD simulation results for a virtual trial of one experimental setting and one dose
#'
#' This function is called at runtime by \code{\link{simulate_OneVirtualTrial}}.
#' @md
#'
#' @param simPKPD A data.frame representing the merge (join) for data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}, and
#' individual and experimental covariates. This argument is set automatically by the calling function.
#' @param trialFileName A character string denoting the trial .rds file on the disk. This file can be loaded to get information about the
#' virtual subjects, such as individual and experimental covariates.
#' @param varPKPD a character vector indicating time varying PKPD variables of the simulated PKPD model. All of these
#' should be names of columns in \code{simPKPD}.
#' @param percentiles a numeric vector indicating percentiles (default: \code{c(5, 50, 95)}) for summarizing
#' individual PKPD variables and clinical end-points.
#' @param samplePDtime a numeric vector (default \code{seq(0,28*24,24)}) of time-points at which PL is
#' needed for summarizing individual clinical endpoints. This should contain aprTime.
#' @param LLOQ.PD a numeric indicating log parasitemia LLOQ (default: \code{log(10000)}).
#' @param aprTime a numeric indicating time after beginning of treatment in days at which APR should be evaluated. Default: 28.
#' This argument is passed to function \code{\link{evaluate_APR}}.
#' @param extrapTol,FLAGextrapolateTime These arguments are passed to function \code{\link{evaluate_APR}}
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param concCOL     Name of column containing values for drug concentration (Default: \code{"Cc"}).
#' @param CIlevel     Numeric containing the confidence interval in percent to be estimated (Default: 90)
#' @param usubjidCOL  Name(s) of column(s) of `simPKPD` uniquely describing individuals (Default: \code{c("IndivID","USUBJID")}).
#' @param trialCOL    Name(s) of column(s) of `simPKPD` uniquely describing trial (Default: \code{c("TrialID")}).
#' @param scenCOL     Name(s) of column(s) of `simPKPD` uniquely describing scenario (Default: \code{c("ScenID","ExpID","DoseID","Dose","nbrDoses")})
#'
#'
#' @details This function is designed to be passed as argument of \code{\link{simulate_VirtualTrials}}. Currently, the function
#' assumes that the \code{simPKPD} has a column PL denoting natural log parasitemia.
#' In addition to the primary parasitological clinical endpoints, this function also summarizes secondary endpoints for PK and PD.
#'
#' @return a named list with the following members:
#'
#' * ScenID: an integer
#' * ExpID: an integer
#' * DoseID: an integer
#' * TrialID: an integer
#' * summaryPKPD.ByTrial: a data.frame
#' * summaryClinEnd.ByTrial: a data.frame of 1 row
#' * summaryClinEnd.ByIndiv: a data.frame of as many rows as there are virtual subjects per trial
#' * summaryPK.ByTrial: a data.frame of 1 row
#' * summaryPK.ByIndiv: a data.frame of as many rows as there are virtual subjects per trial
#'
#' @export
summarizeTrial_APRLongitudinal_PK_secondary <- function(
    simPKPD,
    trialFileName,
    varPKPD             = c("Cc","Eff","Kkill","PL"),
    percentiles         = c(5, 50, 95),
    samplePDtime        = seq(0,28*24,24),
    LLOQ.PD             = log(10000),
    aprTime             = 28,
    extrapTol           = 48,
    FLAGextrapolateTime = TRUE,
    timeCOL     = "TIME",
    concCOL     = "Cc",
    CIlevel     = 90,
    usubjidCOL  = c("IndivID","USUBJID"),
    trialCOL    = c("TrialID"),
    scenCOL     = c("ScenID","ExpID","DoseID","Dose","nbrDoses")
) {

  ScenID  <- unique(simPKPD$ScenID)
  ExpID   <- unique(simPKPD$ExpID)
  DoseID  <- unique(simPKPD$DoseID)
  TrialID <- unique(simPKPD$TrialID)

  stopifnot(
    "summarizeTrial_APRLongitudinal_minPara: function must be called on simPKPD for a single ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)
  stopifnot("summarizeTrial_APRLongitudinal_minPara: All variables in varPKPD should be column names in simPKPD." = is.character(varPKPD) && isTRUE(all(varPKPD %in% names(simPKPD))))
  stopifnot('summarizeTrial_APRLongitudinal_minPara: simPKPD should have at least the columns: ID, ScenID,ExpID,DoseID,TrialID,IndivID,USUBJID,WT0' =
              isTRUE(all(c("ID", "ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID","WT0") %in% names(simPKPD))))
  stopifnot("summarizeTrial_APRLongitudinal_minPara: simPKPD must be a data.frame and must have a column ID" = is.data.frame(simPKPD) && "ID" %in% names(simPKPD))
  stopifnot("summarizeTrial_APRLongitudinal_minPara: aprTime*24 should be in samplePDtime." = (aprTime*24) %in% samplePDtime)

  # Get PKPD summary:
  if(length(unique(simPKPD$ID))>1){
    summaryPKPD.ByTrial <- summaryByTrial(dataInput   = simPKPD,
                                          percentiles = percentiles,
                                          varCOL      = varPKPD,
                                          usubjidCOL  = c("ID","IndivID","USUBJID"),
                                          trialCOL    = unique(c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")))
  }else{
    summaryPKPD.ByTrial <- simPKPD
  }

  # Add curethreshold and PLcure:
  simPKPD$CureThreshold <- ifelse(simPKPD$WT0 <= 35,
                                  log(1/(simPKPD$WT0*80)),
                                  log(1/(simPKPD$WT0*70)))
  simPKPD <- simPKPD[,PLcure:={
    Tcure <- NA_real_
    if(any(PL <= CureThreshold )){
      Tcure <- min(TIME[PL < CureThreshold[1]])
    }
    PLcure <- PL
    PLcure[TIME>=Tcure] <- CureThreshold[1]
    PLcure
  }, by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Get individual primary clinical endpoint:
  summaryClinEndPrimary.ByIndiv <- simPKPD[,{
    interpResult_i <- approx(x = .SD$TIME, y = .SD$PLcure, xout = samplePDtime, ties = list("ordered", mean), rule = 2)

    dat_i          <- data.frame(USUBJID          = USUBJID[1],
                                 TIME             = interpResult_i$x,
                                 PL               = interpResult_i$y,
                                 WT0              = WT0[1],
                                 CureThreshold    = CureThreshold[1],
                                 stringsAsFactors = FALSE)
    dAPR_i <- evaluate_APR(dat_i,
                           LLOQ                = LLOQ.PD,
                           timeCOL             = "TIME",
                           paraCOL             = "PL",
                           weightCOL           = "WT0",
                           curethreshCOL       = "CureThreshold",
                           aprTime             = aprTime,
                           extrapTol           = extrapTol,
                           FLAGextrapolateTime = FLAGextrapolateTime,
                           Plog                = TRUE,
                           FLAGverbose         = FALSE)
    dAPR_i
  },
  by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Summarize 'dataInput' by Trial:
  if(length(unique(summaryClinEndPrimary.ByIndiv$ID))>1){
    summaryClinEndPrimary.ByTrial <- summaryByTrial(dataInput   = summaryClinEndPrimary.ByIndiv,
                                                    percentiles = percentiles,
                                                    varCOL      = setdiff(names(summaryClinEndPrimary.ByIndiv),
                                                                          c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")),
                                                    usubjidCOL  = c("ID","IndivID","USUBJID"),
                                                    trialCOL    = c("ScenID", "ExpID", "DoseID", "TrialID"))
  }else{
    summaryClinEndPrimary.ByTrial <- summaryClinEndPrimary.ByIndiv
  }



  # Get individual secondary clinical endpoint:
  summaryClinEndSecondary.ByIndiv <- simPKPD[,{
    dat_i <- approx(x = .SD$TIME, y = .SD$PL, xout = unique(c(samplePDtime, 0, 48, 24*28)), ties = list("ordered", mean), rule = 2)

    #  min para
    minPara <- min(.SD$PL, na.rm = TRUE)

    # parasite clearnace half life
    # index of fist cross of PLbase -log(2)
    PLbase <- dat_i$y[dat_i$x == 0]
    PLhalfCross <- match(TRUE, .SD$PL < PLbase-log(2), nomatch = length(.SD$PL))
    approxIndex <- c(PLhalfCross-1,PLhalfCross)
    thalf <- approx(x = .SD$PL[approxIndex], y = .SD$TIME[approxIndex], xout = PLbase-log(2) , rule = 2)$y

    dScn_i <- data.frame(
      minPara = minPara,
      tminPara = .SD$TIME[which.min(.SD$PL)],
      # Parasite elimination half-life
      thalf = thalf,

      # PRRs
      PRR48 = (dat_i$y[dat_i$x == 0] - dat_i$y[dat_i$x == 48])/log(10),
      PRRtotal = (dat_i$y[dat_i$x == 0] - minPara)/log(10),
      PRRday28 = (dat_i$y[dat_i$x == 0] - dat_i$y[dat_i$x == 28*24])/log(10),
      PRRtminday28 = (minPara - dat_i$y[dat_i$x == 28*24])/log(10)
    )
    dScn_i
  },
  by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Summarize 'dataInput' by Trial:
  if(length(unique(summaryClinEndSecondary.ByIndiv$ID))>1){
    summaryClinEndSecondary.ByTrial <- summaryByTrial(dataInput   = summaryClinEndSecondary.ByIndiv,
                                                      percentiles = percentiles,
                                                      varCOL      = setdiff(names(summaryClinEndSecondary.ByIndiv),
                                                                            c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")),
                                                      usubjidCOL  = c("ID","IndivID","USUBJID"),
                                                      trialCOL    = c("ScenID", "ExpID", "DoseID", "TrialID"))
  }else{
    summaryClinEndSecondary.ByTrial <- summaryClinEndSecondary.ByIndiv
  }

  # Summarize PK for multiple PK loop over all compounds
  KeyPKlist <- lapply(seq_along(concCOL), function(k) {


    KeyPK <- simPKPD[,{
      thisconcCOL <- get(concCOL[[k]])
      # make all positive
      thisconcCOL[thisconcCOL <= 1e-09] <- 1e-09
      list(AUCinf = estimate_AUC(Time    = get(timeCOL),
                                 Conc    = thisconcCOL,
                                 AUCtype = "inf"),
           Cmax   = max(thisconcCOL),
           Tmax   = get(timeCOL)[which.max(thisconcCOL)])
    },
    by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

    # rename columns
    torename <- !(names(KeyPK) %in% c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID"))
    names(KeyPK)[torename] <- paste(names(KeyPK)[torename], concCOL[[k]], sep = "_")
    KeyPK
  })
  # merge into one data.frame
  mymerge <- function(x, y) merge(x, y, by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID"))
  summaryKeyPK.byIndiv <- Reduce(mymerge, KeyPKlist)

  # Summarize 'dataInput' by Trial:
  if(length(unique(summaryKeyPK.byIndiv$ID))>1){
    summaryKeyPK.ByTrial <- summaryByTrial(dataInput   = summaryKeyPK.byIndiv,
                                           percentiles = percentiles,
                                           varCOL      = setdiff(names(summaryKeyPK.byIndiv),
                                                                 c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")),
                                           usubjidCOL  = c("ID","IndivID","USUBJID"),
                                           trialCOL    = c("ScenID", "ExpID", "DoseID", "TrialID"))
  }else{
    summaryKeyPK.ByTrial <- summaryKeyPK.byIndiv
  }


  # Produce output
  result <- list(
    summaryPKPD.ByTrial    = summaryPKPD.ByTrial,
    summaryClinEnd.ByIndiv = Reduce(mymerge, list(summaryClinEndPrimary.ByIndiv, summaryClinEndSecondary.ByIndiv, summaryKeyPK.byIndiv)),
    summaryClinEnd.ByTrial = rbind(summaryClinEndPrimary.ByTrial, summaryClinEndSecondary.ByTrial ,summaryKeyPK.ByTrial)
  )

  # Output:
  return(result)
}

#' Summarize PKPD simulation results for a virtual trial of one experimental setting and one dose
#'
#' This function is called at runtime by \code{\link{simulate_OneVirtualTrial}}.
#' @md
#'
#' @param simPKPD A data.frame representing the merge (join) for data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}, and
#' individual and experimental covariates. This argument is set automatically by the calling function.
#' @param trialFileName A character string denoting the trial .rds file on the disk. This file can be loaded to get information about the
#' virtual subjects, such as individual and experimental covariates.
#' @param varPKPD a character vector indicating time varying PKPD variables of the simulated PKPD model. All of these
#' should be names of columns in \code{simPKPD}.
#' @param percentiles a numeric vector indicating percentiles (default: \code{c(5, 50, 95)}) for summarizing
#' individual PKPD variables and clinical end-points.
#' @param LLOQ.PD a numeric indicating log parasitemia LLOQ (default: \code{log(10000)}).
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param PbloodCOL Name of column containing values for blood stage parasitemia (Default: \code{"PBlood"}).
#' @param Plog Indicate if parasitemia is in the log or linear scale (Default: `TRUE` which means it is logged).
#' @param FLAGinterpolateTime A logical indicating if the PL measurements should be interpolated within the simulation period.
#' @param outputNames
#' @details This function is designed to be passed as argument of \code{\link{simulate_VirtualTrials}}.It performs
#' a summary of individual clinical endpoints using \code{\link{MMVmalaria:::evaluate_BreakthroughEvent}}. It does
#' not perform summaries by trial, or across trials.

#' @return a named list with the following members:
#' * summaryPKPD.ByTrial: a data.frame
#' * summaryClinEnd.ByIndiv: a data.frame of as many rows as there are virtual subjects per trial

#' @export
#' @author Catalina Barcelo (MMV), Sam Jones (MMV)
#' @family Summarize Functions
summarizeTrial_ChemoBreakthrough <- function(
    simPKPD,
    trialFileName,
    varPKPD             = c("Cc","PLiver","PBlood"),
    percentiles         = c(5, 50, 95),
    LLOQ.PD             = log(10000),
    timeCOL     = "TIME",
    PbloodCOL   = "PBlood",
    Plog                = TRUE,
    FLAGinterpolateTime = FALSE) {

  ScenID  <- unique(simPKPD$ScenID)
  ExpID   <- unique(simPKPD$ExpID)
  DoseID  <- unique(simPKPD$DoseID)
  TrialID <- unique(simPKPD$TrialID)


  stopifnot("summarizeTrial_ChemoBreakthrough: function must be called on simPKPD for a single ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)
  stopifnot("summarizeTrial_ChemoBreakthrough: All variables in varPKPD should be column names in simPKPD." = is.character(varPKPD) && isTRUE(all(varPKPD %in% names(simPKPD))))
  stopifnot('summarizeTrial_ChemoBreakthrough: simPKPD should have at least the columns: ID, ScenID,ExpID,DoseID,TrialID,IndivID,USUBJID' =
              isTRUE(all(c("ID", "ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID") %in% names(simPKPD))))
  stopifnot("summarizeTrial_ChemoBreakthrough: simPKPD must be a data.frame and must have a column ID" = is.data.frame(simPKPD) && "ID" %in% names(simPKPD))
  stopifnot("summarizeTrial_ChemoBreakthrough: PbloodCOL must be a valid column in simPKPD" = PbloodCOL %in% names(simPKPD))

  #   Get PKPD summary:
  if(length(unique(simPKPD$ID))>1){
    summaryPKPD.ByTrial <- summaryByTrial(dataInput   = simPKPD,
                                          percentiles = percentiles,
                                          varCOL      = varPKPD,
                                          # These are hard coded for the purpose of summaryPKPD.ByTrial
                                          usubjidCOL  = c("ID","IndivID","USUBJID"),
                                          trialCOL    = unique(c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")))
  }else{
    summaryPKPD.ByTrial <- simPKPD
  }

  #   Log-transform parasitemia
  #   Data.table call to PbloodCOL allows flexible naming of blood parasitemia compartment
  simPKPD <- simPKPD[, PL := {
    PL <- ifelse(.SD[[PbloodCOL]] < 0.0001, log(0.0001), log(.SD[[PbloodCOL]]))

    list(PL = PL)
  }, .SDcols = PbloodCOL]

  #   Evaluate events for each individual
  summaryClinEnd.ByIndiv <- simPKPD[, {
    events <- evaluate_BreakthroughEvent(dataSim = .SD,
                                         LLOQ    = LLOQ.PD,
                                         timeCOL = "timeCOL",
                                         paraCOL = "PL",
                                         Plog    = Plog,
                                         FLAGinterpolateTime = FLAGinterpolateTime)
    events
  },
  by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Produce output:
  result <- list(
    summaryPKPD.ByTrial    = summaryPKPD.ByTrial,
    summaryClinEnd.ByIndiv = summaryClinEnd.ByIndiv)

  # Output:
  return(result)
}

#' Summarize PKPD simulation results for a virtual trial of one experimental setting and one dose
#'
#' This function is called at runtime by \code{\link{simulate_OneVirtualTrial}}.
#' @md
#'
#' @param simPKPD A data.frame representing the merge (join) for data.frame returned by \code{\link{IQRtools::sim_IQRmodel}}, and
#' individual and experimental covariates. This argument is set automatically by the calling function.
#' @param trialFileName A character string denoting the trial .rds file on the disk. This file can be loaded to get information about the
#' virtual subjects, such as individual and experimental covariates.
#' @param varPKPD a character vector indicating time varying PKPD variables of the simulated PKPD model. All of these
#' should be names of columns in \code{simPKPD}.
#' @param percentiles a numeric vector indicating percentiles (default: \code{c(5, 50, 95)}) for summarizing
#' individual PKPD variables and clinical end-points.
#' @param LLOQ.PD a numeric indicating log parasitemia LLOQ (default: \code{log(10000)}).
#' @param timeCOL     Name of column containing values for simulation time (Default: \code{"TIME"}).
#' @param PbloodCOL Name of column containing values for blood stage parasitemia (Default: \code{"PBlood"}).
#' @param outputNames String of output names being passed in, should correspond to available output names in structural model 
#' @param Plog Indicate if parasitemia is in the log or linear scale (Default: `TRUE` which means it is logged).
#' @param FLAGinterpolateTime A logical indicating if the PL measurements should be interpolated within the simulation period.
#' @param IRdenominator Character string corresponding to the parameter to be used as the denominator of calculating
#' the incidence rate - corresponds to elements of survival objects generated by \code{\link{survival:::survfit}}.
#' (Default: \code{"n.start"})
#' @param TimePreDose Numeric, used when first dose is not given at the first time-step (i.e., to seed pre-infections). Allows
#' the summarize function to automatically disregard any events that occur before the first dose is given.

#' @details This function is designed to be passed as argument of \code{\link{simulate_VirtualTrials}}.It performs
#' a summary of individual clinical endpoints using \code{\link{MMVmalaria:::evaluate_BreakthroughEvent}}. Survival
#' objects are created for each simulated trial, and culmulative incidence calculated for each trial using the Kaplan-Meier
#' estimator from the \code{\link{survival:::survfit}}

#' @return a named list with the following members:
#' * summaryPKPD.ByTrial: a data.frame
#' * summaryClinEnd.ByIndiv: a data.frame of as many rows as there are virtual subjects per trial
#' * summaryClinEnd.Surv : a data.frame of as many rows as there are virtual subjects per trial showing survival time
#' * surv.ByTrial : survival object generated using survfit for each simulated trial
#' * summarySurv.ByTrial : a data frame of as many rows as there are simulated trials, containing culmulative incidence per trial
#' * removedIndv.Surv : a data frame containing USUBJID of individuals who were removed due to excluding pre-dose events

#' @export
#' @author Sam Jones (MMV)
#' @family Summarize Functions

summarizeTrial_ChemoSurvival <- function(
    simPKPD,
    trialFileName,
    varPKPD             = c("Cc","PLiver","PBlood"),
    percentiles         = c(5, 50, 95),
    LLOQ.PD             = log(10000),
    timeCOL     = "TIME",
    PbloodCOL   = "PBlood",
    outputNames = NULL,
    Plog                = TRUE,
    FLAGinterpolateTime = FALSE,
    IRdenominator       = "n.start",
    TimePreDose         = 0
) {

  ScenID  <- unique(simPKPD$ScenID)
  ExpID   <- unique(simPKPD$ExpID)
  DoseID  <- unique(simPKPD$DoseID)
  TrialID <- unique(simPKPD$TrialID)

  # Print warning for timePreDose
  # cat("TimePreDose is", TimePreDose, ". Survival time will be calculated by removing individuals with events occuring between timePreDose and time of first dose")
  
  # Change output names first:
  if(!is.null(outputNames)){
    data.table::setnames(simPKPD, old = names(outputNames), new = outputNames)
  }
  
  stopifnot(
    "summarizeTrial_ChemoSurvival: function must be called on simPKPD for a single ScenID, ExpID, DoseID, TrialID combination." =
      length(ScenID) == 1 && length(ExpID) == 1 && length(DoseID) == 1 && length(TrialID) == 1)
  stopifnot("summarizeTrial_ChemoSurvival: All variables in varPKPD should be column names in simPKPD." = is.character(varPKPD) && isTRUE(all(varPKPD %in% names(simPKPD))))
  stopifnot('summarizeTrial_ChemoSurvival: simPKPD should have at least the columns: ID, ScenID,ExpID,DoseID,TrialID,IndivID,USUBJID' =
              isTRUE(all(c("ID", "ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID") %in% names(simPKPD))))
  stopifnot("summarizeTrial_ChemoSurvival: simPKPD must be a data.frame and must have a column ID" = is.data.frame(simPKPD) && "ID" %in% names(simPKPD))
  stopifnot("summarizeTrial_ChemoSurvival: PbloodCOL must be a valid column in simPKPD" = PbloodCOL %in% names(simPKPD))

  # Get PKPD summary:
  if(length(unique(simPKPD$ID))>1){
    summaryPKPD.ByTrial <- summaryByTrial(dataInput   = simPKPD,
                                          percentiles = percentiles,
                                          varCOL      = varPKPD,
                                          usubjidCOL  = c("ID","IndivID","USUBJID"),
                                          trialCOL    = unique(c("ScenID","ExpID","DoseID","Dose","nbrDoses","TrialID","TIME")))
  }else{
    summaryPKPD.ByTrial <- simPKPD
  }

  #   Log-transform parasitemia
  #   Data.table call to PbloodCOL allows flexible naming of blood parasitemia compartment
  simPKPD <- simPKPD[, PL := {
    PL <- ifelse(.SD[[PbloodCOL]] < 0.0001, log(0.0001), log(.SD[[PbloodCOL]]))

    list(PL = PL)
  }, .SDcols = PbloodCOL]

  #   Evaluate events
  summaryClinEnd.ByIndiv <- simPKPD[, {
    events <- evaluate_BreakthroughEvent(dataSim = .SD,
                                         LLOQ    = LLOQ.PD,
                                         timeCOL = timeCOL,
                                         paraCOL = "PL",
                                         Plog    = Plog,
                                         FLAGinterpolateTime = FLAGinterpolateTime)
    events
  },
  by = c("ID","ScenID","ExpID","DoseID","TrialID","IndivID","USUBJID")]

  # Flag individuals to remove from analysis who suffer event in the time pre dose
  indvToRemove <- dplyr::filter(summaryClinEnd.ByIndiv, TIME<=TimePreDose)
  # Generate the ClindEnd.Surv and adjust TIME to reflect time observed without event after having removed the pre-dose time.
  summaryClinEnd.Surv <- dplyr::anti_join(summaryClinEnd.ByIndiv, indvToRemove)%>%
    mutate(TIME = TIME - TimePreDose)

  surv.ByTrial <- survival::survfit(survival::Surv(`TIME`, `VALUE`) ~ 1, data = summaryClinEnd.Surv)
  summarySurv.ByTrial <- summary(surv.ByTrial)$table
  # Add info
  summarySurv.ByTrial["ScenID"] <- ScenID
  summarySurv.ByTrial["ExpID"] <- ExpID
  summarySurv.ByTrial["DoseID"] <- DoseID
  summarySurv.ByTrial["TrialID"] <- TrialID
  # Calculate IR for trial, allowing denominator to vary in function argument if needed
  summarySurv.ByTrial["CulmulativeIncidence"] <- summarySurv.ByTrial["events"] / summarySurv.ByTrial[IRdenominator]

  # Produce output:
  result <- list(
    summaryPKPD.ByTrial    = summaryPKPD.ByTrial,
    summaryClinEnd.ByIndiv = summaryClinEnd.ByIndiv,
    summaryClinEnd.Surv = summaryClinEnd.Surv,
    surv.ByTrial = surv.ByTrial,
    summarySurv.ByTrial = summarySurv.ByTrial,
    removedIndv.Surv = indvToRemove)

  # Output:
  return(result)
}
