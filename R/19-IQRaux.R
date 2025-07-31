
# Workhorse for sim_IQRmodelConsistent
#
# This function is the original sim_IQRmodelConsistent function,
# which was hidden under another layer of abstraction to allow for parallelization.
# See the main function for documentation
doOneChunk_sim_IQRmodelConsistent <- function(simIQRmodelSpec) {
  
  # Prepare handling of special cases
  ## Special case 1: length(simtime) == 1
  FLAGTimeLength1 <- length(simIQRmodelSpec$simtime) == 1
  if (FLAGTimeLength1) {
    simIQRmodelSpec$simtime <- c(simIQRmodelSpec$simtime, simIQRmodelSpec$simtime + 1e-6)
  }
  
  # Simulate
  simres <- do.call(sim_IQRmodel, simIQRmodelSpec)
  if (is.null(simres)) {return(simres)} # Will be captured in sim_IQRmodelConsistent
  
  # Post-process
  ## Special case 1: length(simtime) == 1
  if (FLAGTimeLength1) {
    simres <- simres[simres$TIME == simIQRmodelSpec$simtime[1], , drop = FALSE]
  }
  ## Special case 2: No ID column
  if (! "ID" %in% names(simres)) {
    simres$ID <- unique(simIQRmodelSpec$eventTable$ID)
    # Move ID to front
    simres <- simres[,c("ID", setdiff(names(simres), c("ID")))]
  }
  ## Special case 3: No simNames attribute
  if (!"simNames" %in% names(attributes(simres))) {
    attr(simres, "simNames") <- setdiff(names(simres), c("ID", "TIME"))
  }
  
  simres
}



# Split a simIQRmodelSpec into simulation tasks with fewer IDs per eventTable
#
# @param simIQRmodelSpec simIQRmodelSpec
#
# @return simIQRmodelSpecList with smaller eventTables
# @export
#
# @examples
split_simIQRmodelSpecPerCore <- function(simIQRmodelSpec) {
  
  mc.cores <- attr(simIQRmodelSpec, "mc.cores")
  if (is.null(mc.cores)) {
    mc.cores <- 1
  }
  
  # Disregard the chunksize functionality for now - but in future it could be cool
  # to add the option of passing a post-processor function to quickly post-process
  # chunks - if the problem of memory overflow ever arises or if postprocessing of
  # the full simres gets too slow
  eTList <- split_eventTable(simIQRmodelSpec$eventTable, mc.cores = mc.cores, chunksize = 100000000)
  
  simIQRmodelSpecList <- lapply(eTList, function(x) {
    y <- simIQRmodelSpec
    y$eventTable <- x
    y
  })
  
  simIQRmodelSpecList
}

pred_IQRnlmeProject_Consistent <- function (project, model = NULL, data = NULL, simtime = NULL, 
                                            abs0inputs = NULL, abs0Tk0param = NULL, opt_abstol = 1e-06, 
                                            opt_reltol = 1e-06) 
{
  if (is.null(model)) {
    model <- getModel_IQRnlmeProject(project)
  }
  if (is.null(data)) {
    data <- getData_IQRnlmeProject(project)
  }
  indiv <- getIndivParameters_IQRnlmeProject(project)
  header <- getHeader_IQRnlmeProject(project)
  if (is.null(simtime)) {
    dataObs <- data[data$EVID == 0, ]
    simtime <- lapply(split(dataObs, dataObs$ID), function(d) unique(d$TIME))
  }
  if ("ABSORPTION0" %in% header$DOSINGTYPES) {
    if (is.null(abs0inputs) || is.null(abs0Tk0param)) 
      stopIQR("ABSORPTION0 present in model. Specify 'abs0inputs' and 'abs0Tk0param' arguments!")
  }
  dataparam <- dplyr::left_join(data, indiv, by = c("USUBJID", 
                                                    "ID"))
  if (header$REGRESSIONNAMES[1] == "") {
    regressors <- NULL
  }
  else {
    regressors <- header$REGRESSIONNAMES
  }
  regression <- c(regressors, setdiff(names(indiv), c("USUBJID", 
                                                      "ID")))
  eT <- IQReventTable(data = dataparam, regression = regression, 
                      abs0inputs = abs0inputs, abs0Tk0param = abs0Tk0param)
  simspec <- createSimIQRmodelSpec(model = model, eventTable = eT, simtime = simtime, 
                                   opt_abstol = opt_abstol, opt_reltol = opt_reltol)
  res <- sim_IQRmodelConsistent(simIQRmodelSpec = simspec)
  res
}

# -------------------------------------------------------------------------#
# Helper functions ----
# -------------------------------------------------------------------------#

# Run some consistency checks for a simIQRmodelSpec
#
# If min(simtime) differs from min(eventTable$TIME) for any ID, the simulation might do unintended things.
# Depending on the opt_timeHandling, emit a warning or an error so that the user can act accordingly.
#
# Please read the code for the implemented checks - I don't want to create too much redundancy and the code is straight forward.
#
# @param simIQRmodelSpec Result from [createSimIQRmodelSpec]
# @param opt_timeHandling See [createSimIQRmodelSpec]
#
# @return Nothing, or a warning or an error
# @author Daniel Lill (daniel.lill@intiquan.com)
# @md
checkSimSeqTimes <- function(simIQRmodelSpec,
                             opt_timeHandling = c("augmentEventTable", "warn", "stop")) {
  
  opt_timeHandling <- match.arg(opt_timeHandling)
  
  simMinTime <- sapply(simIQRmodelSpec$simtime, min)
  eTMinTimes <- vapply(split(simIQRmodelSpec$eventTable, simIQRmodelSpec$eventTable$ID),
                       function(x) min(x$TIME),
                       FUN.VALUE = 0)
  eTMinStartsEarlier <- eTMinTimes < simMinTime
  eTMinStartsLater   <- eTMinTimes > simMinTime
  
  if (opt_timeHandling %in% c("stop", "warn")) {
    error <- ""
    if (any(eTMinStartsEarlier)) {
      error <- c(error, "The following IDs have their first event before min(simtimes):",
                 paste0(names(eTMinStartsEarlier)[eTMinStartsEarlier], collapse = ", "))
    }
    if (any(eTMinStartsLater)) {
      error <- c(error, "The following IDs have their first event after min(simtimes):",
                 paste0(names(eTMinStartsLater)[eTMinStartsLater], collapse = ", "))
    }
    if (any(eTMinStartsEarlier) | any(eTMinStartsLater)) {
      error <- c("Time inconsistencies in the simulation specification:",
                 "-----------------------------------------------------",
                 error)
      if (opt_timeHandling == "stop") {
        stop(paste0(error, collapse = "\n"))
      } else if (opt_timeHandling == "warn") {
        error <- c(error, "The simulations might be wrong, because the initial conditions from the previous step are applied at min(eventTable$TIME)")
        warning(paste0(error, collapse = "\n"))
      }
    }
  }
  
}

# -------------------------------------------------------------------------#
# Note for Reviewer ----
# This function is from IQauxfn and can be considered validated.
# -------------------------------------------------------------------------#

#' Create an IQRdataGENEARL from a design specification
#'
#' Populate the data with [populateGDFWithSimValues()]
#'
#' @param design A list which is formatted like the `design`-argument of [IQRtools::IQRpopEDdb()]
#'
#' @return An [IQRtools::IQRdataGENERAL()]
#' @export
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @family IQRdataGENERAL
#' @family Simulation
#' @importFrom data.table data.table rbindlist
#'
#' @examples
#' design <- list(
#'   time = list(
#'     group1 = list(
#'       OUTPUT1 = c(1, 2, 5, 10),
#'       OUTPUT2 = c(1, 5, 20)
#'     ),
#'     group2 = list(
#'       OUTPUT1 = c(1, 2, 5, 10),
#'       OUTPUT2 = c(1, 5, 20)
#'     )
#'   ),
#'   dosing = list(
#'     group1 = IQRtools::IQRdosing(
#'       TIME = 0, ADM = 1, AMT = 40, TINF = 0
#'     ),
#'     group2 = IQRtools::IQRdosing(
#'       TIME = 0, ADM = 1, AMT = 40, TINF = 0
#'     )),
#'   groupsize = list(
#'     group1 = 12, group2 = 3
#'   )
#' )
#' create_GDF_fromDesignSpec(design)
create_GDF_fromDesignSpec <- function(design) {

  # Expand USUBJIDs
  USUBJIDs <- lapply(setNames(nm = names(design$groupsize)), function(nm) {
    data.table::data.table(USUBJID = paste0(nm, "_", seq_len(design$groupsize[[nm]])))
  })
  USUBJIDs <- data.table::rbindlist(USUBJIDs, idcol = "TRTNAME")

  # Observations with reasonable default values
  observations <- lapply(setNames(nm = names(design$time)), function(nmTRTNAME) {
    obsTimesList <- design$time[[nmTRTNAME]]
    x <- lapply(setNames(nm = names(obsTimesList)), function(nmOUTPUT) {
      data.table::data.table(
        TIME     = obsTimesList[[nmOUTPUT]],
        TIMEUNIT = "HOURS",
        NAME     = nmOUTPUT,
        YTYPE    = as.numeric(gsub("OUTPUT", "", nmOUTPUT)),
        DV       = 0,
        VALUE    = 0,
        UNIT     = "ug/mL",
        LLOQ = 0, CENS = 0, EVID = 0, MDV = 0,
        AMT = 0, ADM = NA, TINF = 0, II = 0, ADDL = 0, DURATION = 0, RATE = 0, ROUTE = NA_character_, DOSE = NA
      )
    })
    data.table::rbindlist(x)
  })
  observations <- data.table::rbindlist(observations, idcol = "TRTNAME")

  # Dosings with reasonable default values
  dosings <- lapply(setNames(nm = names(design$dosing)), function(nm) {
    d <- data.table::data.table(design$dosing[[nm]])
    d[,`:=`(
      TIMEUNIT = "HOURS",
      NAME     = paste0("INPUT", ADM),
      YTYPE = 0, DV = 0, VALUE = AMT, UNIT = "mg", ROUTE = "ORAL", DOSE = AMT,
      LLOQ = NA, CENS = 0, EVID = 1, MDV = 1,
      II = 0, ADDL = 0, DURATION = TINF, RATE = AMT/TINF)]
  })
  dosings <- data.table::rbindlist(dosings, idcol = "TRTNAME")

  # Combine everything and create last columns
  design0 <- data.table::rbindlist(list(dosings, observations), use.names = TRUE, fill = TRUE)
  design0 <- merge(design0,USUBJIDs, by = "TRTNAME", allow.cartesian = TRUE)
  design0 <- design0[order(USUBJID, TIME, YTYPE)]
  design0 <- data.frame(design0)
  IQRtools::IQRdataGENERAL(design0)
}

# -------------------------------------------------------------------------#
# Note for Reviewer ----
# This original function populateGDFWithSimValues is from IQauxfn and can be considered validated.
# This function was modified to use sim_IQRmodel_consistent from IQauxfn to allow for simulation at negative time, i.e. prior to first dose
# The modified function has not been reviewed
# -------------------------------------------------------------------------#

#' Simulate DV values and place them into a dataGENERAL
#'
#' The dataGEN most likely should come from [create_GDF_fromDesignSpec()].
#' You can add individual parameters to simulate based on individuals.
#' Simulated values will be placed into the respective VALUE and DV entries of the dataGEN.
#'
#' @param dataGEN An [IQRtools::IQRdataGENERAL()]
#' @param model An [IQRtools::IQRmodel()]
#' @param regression Character vector
#' @param abs0inputs Character vector of abs0inputs (see [IQRtools::IQReventTable()])
#' @param abs0Tk0param  Character vector of abs0 parameters (see [IQRtools::IQReventTable()])
#'
#' @return The `dataGEN` with modified DV and VALUE column
#' @export
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @family IQRdataGENERAL
#' @family Simulation
#' @importFrom data.table data.table melt
#'
#' @examples
#' library(ggplot2)
#' library(IQRtools)
#'
#' # Step 1: Get model
#'
#' model <-IQRtools::IQRmodel(file.path(get_MMVmalariaPath(subdir="inst"),"modelLibrary/PKmodels/model_2cpt_linsat_iv.txt"))
#' # Step 2: Create empty dataGEN
#' design <- list(
#'   time = list(
#'     group1 = list(
#'       OUTPUT1 = c(1, 2, 5, 10)
#'     ),
#'     group2 = list(
#'       OUTPUT1 = c(1, 2, 5, 10)
#'     )
#'   ),
#'   dosing = list(
#'     group1 = IQRtools::IQRdosing(
#'       TIME = 0, ADM = 1, AMT = 40, TINF = 0
#'     ),
#'     group2 = IQRtools::IQRdosing(
#'       TIME = 0, ADM = 1, AMT = 40, TINF = 0
#'     )),
#'   groupsize = list(
#'     group1 = 10, group2 = 10
#'   )
#' )
#' dataGEN <- create_GDF_fromDesignSpec(design)
#'
#' # Step 2.1: Add some individual parameters
#' gpf <- IQRtools::generate_GPFfromIQRmodel(model)
#' estimates <- gpf$estimates
#' estimates[VALUE > 0,`:=`(IIV = 0.1)] # Update IIV Value
#' gpf <- IQRtools::GPF(estimates = estimates)
#' indPars <- sample_GPF(gpf, Nsamples = length(unique(dataGEN$ID)), FLAGid = TRUE)
#' indPars <- indPars$indParamValues
#' dataGEN <- merge(dataGEN, indPars, by = "ID")
#' dataGEN <- IQRdataGENERAL(dataGEN)
#'
#' # Step 3: Populate dataGEN
#' dataGENsim <- populateGDFWithSimValues_consistent(dataGEN, model, regression = setdiff(names(indPars), "ID"), abs0inputs = NULL, abs0Tk0param = NULL)
#'
#' # Visualize
#' IQRtools::IQRggplot(dataGENsim, aes(TIME, DV, group = USUBJID)) +
#'   facet_wrap(~TRTNAME, scales = "free") +
#'   geom_line() +
#'   labs(title = "Simulated data (no residual error)")
populateGDFWithSimValues_consistent <- function(dataGEN, model, regression = NULL, abs0inputs = NULL, abs0Tk0param = NULL) {
  # For preserving attributes of dataGEN. Nothing should change except for the values in DV and VALUE,
  # so we can just take it as it is and replace the new attributes later.
  atr <- attributes(dataGEN)

  
  
  # Simulate with individual parameters and melt simres.
  eT <- IQRtools::IQReventTable(data = dataGEN,
                      regression = regression,
                      abs0inputs = abs0inputs,
                      abs0Tk0param = abs0Tk0param)
  
  # simres <- sim_IQRmodel(model, simtime = sort(unique(c(dataGEN$TIME))),
  #                        eventTable = eT, FLAGoutputsOnly = TRUE)
  # use sim_IQRmodelConsistent
  simIQRmodelSpec1 <- createSimIQRmodelSpec(model, simtime = sort(unique(c(dataGEN$TIME))),
                                            eventTable = eT, FLAGoutputsOnly = TRUE)
  simres <- IQRtools::sim_IQRmodel_strict(simIQRmodelSpec = simIQRmodelSpec1)
  
  simres <- data.table::data.table(simres)
  simres <- data.table::melt(simres, id.vars = intersect(names(simres), c("ID", "TIME")),
                             variable.name = "modelEntityId", variable.factor = FALSE, value.name = "simValue")
  simres[,`:=`(YTYPE = as.numeric(gsub("OUTPUT", "", modelEntityId)))]
  simres[,`:=`(modelEntityId = NULL)]

  # Populate dataGEN with simulated values: Create a copy and update
  dataGENNew <- data.table::data.table(dataGEN)
  dataGENNew <- simres[dataGENNew, on = .NATURAL] # merge but keep row order of dataGENNe
  dataGENNew[!is.na(simValue),`:=`(DV = simValue, VALUE = simValue)]
  dataGENNew[,`:=`(simValue = NULL)]
  dataGENNew <- data.frame(dataGENNew, stringsAsFactors = FALSE)

  # Preserve attributes
  dataGENNew <- IQRtools::IQRdataGENERAL(dataGENNew)
  attributes(dataGENNew) <- atr
  dataGENNew
}