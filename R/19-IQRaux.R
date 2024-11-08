# --- START IQauxfn function validation status ----
#
# Note to Reviewer: If all boxes in this section are ticked (marked with # [ ]),
#   the function can be considered validated for IQdesktop > 2022.11
#
# Function digest: 4fe7c0b0b59d292b46123d9ec8a6f9fe
#
# [x] The function is sufficiently documented according to the IQauxfn documentation standards
# [x] The function's examples cover the most frequent use cases
# [x] The function's examples run without error and are in sync with the expectations
# [x] The function returns informative error messages if user input is faulty.
#
# --- END IQauxfn function validation status ----


#' Verify that input arguments meet certain conditions
#'
#' @param x Any object
#' @param allowNull TRUE: x can be null. If TRUE and x is null, no additional tests are made.
#' @param expectedClass Missing or character denoting the class of x
#' @param expectedMode Missing or character denoting the mode of x
#' @param expectedLength Missing or integer denoting the length of x
#' @param expectedSign Missing or 1 or -1 denoting the sign of x
#' @param expectedNames Names which must at least be present in x
#' @param expectedTestFun A function for which the call `expectedTestFun(x)` returns TRUE (test is passed) or FALSE (test is failed)
#'
#' @return Called for side-effect. If all tests pass, nothing happens. If errors occur, they are collected in informative error messages.
#' @export
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @family Other functions
#'
#' @examples
#' \dontrun{
#' # Checking for integers is a bit tricky and should not necessarily be done with expectedClass
#'
#' # Doubles
#' verifyArg(1, expectedClass = "numeric")    # Gives no error
#' verifyArg(1, expectedClass = "double")     # Gives error
#' verifyArg(1, expectedMode = "numeric")     # Gives no error
#' verifyArg(1, expectedTestFun = is.numeric) # Gives no error, last resort
#'
#' # Integers
#' verifyArg(1L, expectedClass = "numeric") # Gives error
#' verifyArg(1L, expectedMode = "numeric")  # Gives no error
#' }
verifyArg <- function(x, allowNull = FALSE,
                      expectedClass,
                      expectedMode,
                      expectedLength,
                      expectedSign ,
                      expectedNames,
                      expectedTestFun) {
  subx <- substitute(x)
  charx <- deparse(subx)
  
  if (allowNull & is.null(x)) return()
  if (!allowNull & is.null(x)) stop(charx, " must not be NULL")
  
  errors <- c()
  if (!missing(expectedClass) && !expectedClass %in% class(x))
    errors <- c(errors,
                paste0("Class of ", charx, " -----------"),
                paste0("* Expected: '", expectedClass, "'"),
                paste0("* Supplied: '", paste0(class(x), collapse = ","), "'")
    )
  if (!missing(expectedMode) && !expectedMode %in% mode(x))
    errors <- c(errors,
                paste0("Mode of ", charx, " -----------"),
                paste0("* Expected: '", expectedMode, "'"),
                paste0("* Supplied: '", paste0(mode(x), collapse = ","), "'")
    )
  if (!missing(expectedLength) && length(x) != expectedLength)
    errors <- c(errors,
                paste0("Length of ", charx, " -----------"),
                paste0("* Expected: '", expectedLength, "'"),
                paste0("* Supplied: '", paste0(length(x), collapse = ","), "'")
    )
  if (!missing(expectedSign) && !expectedSign %in% c(1,-1))
    stop("expectedSign should be one of 1,-1. ",
         "Please update your call to verifyArg. ",
         "(This is an error in the call to verifyArg itself, not an argument check.)")
  if (!missing(expectedSign) && any(sign(x) != expectedSign))
    # [ ] Add support for sign(0)
    errors <- c(errors,
                paste0("Sign of ", charx, " -----------"),
                paste0("* Expected: '", expectedSign, "'"),
                paste0("* Elements with wrong sign: '", paste0(which(sign(x) != expectedSign), collapse = ","), "'")
    )
  if (!missing(expectedNames) && any(! expectedNames %in% names(x)))
    errors <- c(errors,
                paste0("Required names of ", charx, " -----------"),
                paste0("* Expected: '", paste0(expectedNames, collapse = ","), "'"),
                paste0("* Supplied:  '", paste0(intersect(names(x), expectedNames), collapse = ","), "'"),
                paste0("* Missing:  '", paste0(setdiff(expectedNames, names(x)), collapse = ","), "'"),
                paste0("* Additional:  '", paste0(setdiff(names(x), expectedNames), collapse = ","), "'")
    )
  if (!missing(expectedTestFun)) {
    subTest <- substitute(expectedTestFun)
    charTest <- as.character(subTest)
    test <- expectedTestFun(x)
    if (!test)
      errors <- c(errors,
                  paste0("The following test returned FALSE ----"),
                  paste0("  ", charTest, "(", charx, ")")
      )
  }
  
  if (length(errors))
    stop(paste0(errors, collapse = "\n"))
  
}


# --- START IQauxfn function validation status ----
#
# Note to Reviewer: If all boxes in this section are ticked (marked with # [ ]),
#   the function can be considered validated for IQdesktop > 2022.11
#
# Function digest: 53bde5531099e59603cf1639d70f2d80
#
# [x] The function is sufficiently documented according to the IQauxfn documentation standards
# [x] The function's examples cover the most frequent use cases
# [x] The function's examples run without error and are in sync with the expectations
# [x] The function returns informative error messages if user input is faulty.
#
# --- END IQauxfn function validation status ----

#' sim_IQRmodel with strict output format
#'
#' This function needs [createSimIQRmodelSpec()] and [getDefaultEventTable()] and [verifyArg()]
#'
#' [IQRtools::sim_IQRmodel()] output format depends on arguments of the sim_IQRmodel call.
#'
#' * TIME is very context dependent, eventfor explicit `simtime` values passed as argument.
#' * ID is not always present: If the eventTable contains one ID only, or sim_IQRmodel is called without eventTable, the output an ID column does not contain
#'
#' By contrast, **this** function has consistent output, suitable for programming.
#' * TIME corresponds exactly to the input given in `simtime`
#' * It always contains an ID column
#' * It always has a simNames attribute
#'
#'
#' @param simIQRmodelSpec Argument list for sim_IQRmodel, created by [createSimIQRmodelSpec()]
#'
#' @return `IQRsimres` with consistent format, see description
#' @export
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @family Simulation
#'
#' @examples
#' library(IQRtools)
#' library(ggplot2)
#' fileURL <- "http://192.168.1.12:7804/IQauxfn/inst/02-ExampleModels/ModelLibraryPK/model_2cpt_linsat_iv.txt"
#' fl <- file.path(tempdir(), basename(fileURL))
#' download.file(fileURL, fl)
#'
#' model <- IQRmodel(fl)
#' simIQRmodelSpec1 <- createSimIQRmodelSpec(model, simtime = c(0), parameters =  c(CL = 10), dosingTable = IQRdosing(0,1,1))
#' simres1 <- do.call(sim_IQRmodel, simIQRmodelSpec1)
#' simres1Consistent <- sim_IQRmodelConsistent(simIQRmodelSpec = simIQRmodelSpec1)
sim_IQRmodelConsistent <- function(simIQRmodelSpec) {
  verifyArg(simIQRmodelSpec, expectedClass = "simIQRmodelSpec")
  verifyArg(simIQRmodelSpec$eventTable, allowNull = FALSE, expectedClass = "IQReventTable")
  
  # Split into multi-core simulation
  simIQRmodelSpecList <- split_simIQRmodelSpecPerCore(simIQRmodelSpec)
  
  # Simulate parallelized
  mc.cores <- attr(simIQRmodelSpec, "mc.cores")
  if (is.null(mc.cores)) {
    mc.cores <- 1
  }
  simres <- parallel::mclapply(X = simIQRmodelSpecList,
                               FUN = function(x) {doOneChunk_sim_IQRmodelConsistent(x)},
                               mc.cores = mc.cores,
                               mc.silent = FALSE)
  simres <- do.call(rbind, simres) # Column order should be the same, and rbind conserves attributes
  
  # Handle failures
  failedIDs <- setdiff(simIQRmodelSpec$eventTable$ID, simres$ID)
  if (length(failedIDs)) warningIQR("The following IDs did not succeed - please do.call the original sim_IQRmodel function to get the error message if it was not printed yet.\n",
                                    paste0(failedIDs, collapse = ","))
  
  simres
}


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



# Split eventTable into smaller Tables
#
# If the simulation results get very big, e.g. because sim_IQRmodel is called
# with FLAGoutputsOnly = FALSE, it's faster to split the eventTable into pieces.
# This can also circumvent RAM limiations. Furthermore, it allows to parallelize the simulations.
#
# For few IDs, the function splits the eventTable in chunks of equal size
#
# @param eventTable Object of class IQReventTable
# @param mc.cores Number of cores for parallelization, determining the optimal split size for small eventTables
# @param chunksize How many IDs should be included in one slice? 1000 proved to be fast for the typical DDI models
#
# @importFrom data.table data.table copy ":="
# @importFrom ggplot2 cut_number
#
# @export
# @return list of smaller IQReventTable objects
# @author Daniel Lill (daniel.lill@intiquan.com)
# @md
# @family IQReventTable
#
# @examples
# eT <- data.frame(ID = 1:2000, TIME = 0, ADM = 1, AMT = 1, TINF = 1e-4)
# eT <- `class<-`(eT, c("IQReventTable", "data.frame"))
# split_eventTable(eT, ncores = 10)
# split_eventTable(eT, ncores = 10, chunksize = 50)
# split_eventTable(eT[eT$ID<=5,], ncores = 10)
split_eventTable <- function(eventTable, mc.cores = 1, chunksize = 1000) {
  
  # Helpers: See documentation in MMVIsoboles
  seqminmax <- function(from, to, by) {unique(sort(c(from, to, seq(from, to, by))))}
  cut_groupsize <- function(x, groupsize){cut(seq_along(x), seqminmax(0, length(x), groupsize))}
  
  eT <- data.table::copy(data.table::data.table(eventTable))
  eT_split <- eT[,list(ID = unique(ID))]
  if (length(eT_split$ID) >= mc.cores * chunksize){ # case 1: very many IDs, limit size of each chunk to chunksize
    eT_split[,`:=`(cuts = cut_groupsize(ID, chunksize))]
  } else if (length(eT_split$ID) > mc.cores) {      # case 2: distribute many IDs evenly across cores
    eT_split[,`:=`(cuts = ggplot2::cut_number(ID, n = mc.cores))]
  } else {                                        # case 3: one ID per core
    eT_split[,`:=`(cuts = factor(as.numeric(ID)))]
  }
  eT_split <- eT_split[eT, on = "ID"]
  eT_split <- as.data.frame(eT_split)
  eT <- split(eT,eT_split$cuts)
  eT <- lapply(eT, function(.x) `class<-`(.x, c("IQReventTable", "data.frame")))
  eT
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

# --- START IQauxfn function validation status ----
#
# Note to Reviewer: If all boxes in this section are ticked (marked with # [x]),
#   the function can be considered validated for IQdesktop > 2022.11
#
# Function digest: 2e2978e1792947d0462eeb234dc2916d
#
# [x] The function is sufficiently documented according to the IQauxfn documentation standards
# [x] The function's examples cover the most frequent use cases
# [x] The function's examples run without error and are in sync with the expectations
# [x] The function returns informative error messages if user input is faulty.
#
# --- END IQauxfn function validation status ----

#' Get a default IQReventTable
#'
#' @param model [IQRtools::IQRmodel()]
#' @param dosingTable NULL or [IQRtools::IQRdosing()]
#' @param parameters Named vector to override parameter values
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @family Simulation
#' @family eventTable
#'
#' @return [IQRtools::IQReventTable()]
#' @export
#'
#' @details
#' The two results are equivalent: \cr
#' ```
#' sim_IQRmodel(m)
#' sim_IQRmodel(m, eventTable = getDefaultEventTable(m))
#' ```
#'
#' Same with dosing: \cr
#' ```
#' sim_IQRmodel(m, dosingTable = IQRdosing(c(0,1), 1, 1))
#' sim_IQRmodel(m, eventTable = getDefaultEventTable(m, IQRdosing(c(0,1), 1, 1)))
#' ```
#'
#' @importFrom IQRtools IQRdosing create_IQReventTable
#'
#' @examples
#' library(IQRtools)
#' fileURL <- "http://192.168.1.12:7804/IQauxfn/inst/02-ExampleModels/ModelLibraryPK/model_2cpt_linsat_iv.txt"
#' fl <- file.path(tempdir(), basename(fileURL))
#' download.file(fileURL, fl)
#' m <- IQRmodel(fl)
#' d <- IQRdosing(c(0,1),1,1)
#' eT <- getDefaultEventTable(m, dosingTable = d)
#'
#' simres1 <- sim_IQRmodel(m, dosingTable = d)
#' simres2 <- sim_IQRmodel(m, eventTable = eT)
#'
#' p <- c(CL = 100)
#' eT <- getDefaultEventTable(m, parameters = p)
#'
#' simres1 <- sim_IQRmodel(m, parameters = p)
#' simres2 <- sim_IQRmodel(m, eventTable = eT)
getDefaultEventTable <- function(model, dosingTable = NULL, parameters = NULL) {
  # Collect all parameters in a data.frame called eventGrid ("eG")
  eG <- vapply(model$parameters, function(p) p$value, c("par" = 0.1))
  eG <- eG[grep("INPUT", names(eG), invert = TRUE)]
  eG <- as.data.frame(t(eG))
  eG <- cbind(ID = 1, eG)  # Dummy ID
  
  if (!is.null(parameters)) {
    parameters <- as.data.frame(as.list(parameters))
    eG <- cbind(eG[setdiff(names(eG), names(parameters))], parameters)
  }
  
  # Collect dosingTable ("dT"). If no dosingTable available use a default
  dT <- IQRtools::IQRdosing(TIME = 0, ADM = 1, AMT = 0)
  FLAGDosingProvided <- !is.null(dosingTable)
  if (FLAGDosingProvided)
    dT <- dosingTable
  
  
  # Create the eventTable
  eT <- IQRtools::create_IQReventTable(dT, eG)
  
  # Mini-Hack: - IQRtools::create_IQReventTable does not allow ADM = NA,
  # but for models without any INPUT, ADM needs to be NA for simulation
  if (!FLAGDosingProvided) eT$ADM <- NA
  
  eT
}

# --- START IQauxfn function validation status ----
#
# Note to Reviewer: If all boxes in this section are ticked (marked with # [ ]),
#   the function can be considered validated for IQdesktop > 2022.11
#
# Function digest: 4207eee540f594094d76e4629383677e
#
# [x] The function is sufficiently documented according to the IQauxfn documentation standards
# [x] The function's examples cover the most frequent use cases
# [x] The function's examples run without error and are in sync with the expectations
# [x] The function returns informative error messages if user input is faulty.
#
# --- END IQauxfn function validation status ----


#' Collect and standardize input arguments for sim_IQRmodel
#'
#' This function does four things
#'
#' * Collect input arguments differing from the default values for [IQRtools::sim_IQRmodel()] in a list
#' * Require that `simtime` is specified explicitly
#' * Replace `parameters` and `dosingTable` inputs by `eventTable`
#' * If min(simtime) < min(eventTable$TIME), add a dummy event at min(simtime), so no simtime is cut off ever.
#'
#'
#' @param model See [IQRtools::sim_IQRmodel()], also for other arguments
#' @param opt_timeHandling sanitize times in an eventTable (see description)
#' @param mc.cores For parallelization with parallel::mclapply
#'
#' @return List of its input arguments, with class simIQRmodelSpec
#' @export
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @family Simulation
#'
#' @examples
#' library(IQRtools)
#' library(ggplot2)
#' fileURL <- "http://192.168.1.12:7804/IQauxfn/inst/02-ExampleModels/ModelLibraryPK/model_2cpt_linsat_iv.txt"
#' fl <- file.path(tempdir(), basename(fileURL))
#' download.file(fileURL, fl)
#'
#' model <- IQRmodel(fl)
#' simIQRmodelSpec1 <- createSimIQRmodelSpec(model, simtime = c(0,1), parameters =  c(CL = 10), dosingTable = IQRdosing(0,1,1))
#' do.call(sim_IQRmodel, simIQRmodelSpec1)
createSimIQRmodelSpec <- function (model, simtime,
                                   IC = NULL, parameters = NULL,
                                   dosingTable = NULL, eventTable = NULL, FLAGsensitivity = FALSE,
                                   FLAGuseSensEq = TRUE, sensParams = NULL, FLAGoutputsOnly = FALSE,
                                   opt_eventTimes = NULL, opt_method_stiff = TRUE, opt_abstol = 1e-06,
                                   opt_reltol = 1e-06, opt_minstep = 0, opt_maxstep = 0, opt_initstep = 0,
                                   opt_maxnumsteps = 1e+05, opt_maxerrtestfails = 50, opt_maxorder_stiff = 5,
                                   opt_maxorder_nonstiff = 12, opt_maxconvfails = 10, opt_maxnonlineariter = 3,
                                   opt_usesymjac = TRUE, opt_sens_simultaneous = FALSE, opt_sens_errcon = FALSE,
                                   opt_sens_maxnonlineariter = 3, verbose = FALSE, FLAGprogressBar = FALSE,
                                   # simIQRmodelSpec specific
                                   opt_timeHandling = c("augmentEventTable", "warn", "stop"),
                                   mc.cores = 1
) {
  
  # Input verification
  verifyArg(simtime, allowNull = FALSE)
  opt_timeHandling <- match.arg(opt_timeHandling)
  
  # Collect all args in a list
  fx <- formals()
  l <- lapply(setNames(nm = names(fx)), function(x) eval(parse(text = x)))
  
  # Only return non-default values.
  # Could be implemented with formals(sim_IQRmodel), but let's keep it simple.
  # The values are not likely to change.
  defaultValues <- list(FLAGsensitivity = FALSE,
                        FLAGuseSensEq = TRUE, FLAGoutputsOnly = FALSE,
                        opt_method_stiff = TRUE, opt_abstol = 1e-06,
                        opt_reltol = 1e-06, opt_minstep = 0, opt_maxstep = 0, opt_initstep = 0,
                        opt_maxnumsteps = 1e+05, opt_maxerrtestfails = 50, opt_maxorder_stiff = 5,
                        opt_maxorder_nonstiff = 12, opt_maxconvfails = 10, opt_maxnonlineariter = 3,
                        opt_usesymjac = TRUE, opt_sens_simultaneous = FALSE, opt_sens_errcon = FALSE,
                        opt_sens_maxnonlineariter = 3, verbose = FALSE, FLAGprogressBar = FALSE)
  hasDefaultValue <- vapply(setNames(nm = names(defaultValues)), function(nm) {defaultValues[[nm]] == l[[nm]]}, FUN.VALUE = TRUE)
  hasDefaultValue <- names(hasDefaultValue)[hasDefaultValue]
  l <- l[setdiff(names(l), hasDefaultValue)]
  
  hasDefaultNull <- vapply(setNames(nm = c("sensParams", "opt_eventTimes")), function(nm) is.null(l[[nm]]), FUN.VALUE = TRUE)
  hasDefaultNull <- names(hasDefaultNull)[hasDefaultNull]
  l <- l[setdiff(names(l), hasDefaultNull)]
  
  
  # Move from parameters and dosingTable to eventTable
  if (is.null(l$eventTable)) {
    eventTable <- getDefaultEventTable(model = model, dosingTable = dosingTable, parameters = parameters)
    l$eventTable <- eventTable
    l$parameters <- NULL
    l$dosingTable <- NULL
  }
  
  # if simtime is vector with same simtime for all IDs then expand simtime to be a list with as many entries as IDs in the eventTable
  if (is.atomic(l$simtime)) {
    ids <- unique(l$eventTable$ID)
    l$simtime  <- replicate(length(ids), l$simtime, simplify = FALSE)
    names(l$simtime) <- ids
  }
  
  # Check consistency of times and remove opt_timeHandling from list for do.call(simIQRmodel, l)
  checkSimSeqTimes(simIQRmodelSpec = l, opt_timeHandling = opt_timeHandling)
  if (opt_timeHandling == "augmentEventTable") {
    
    l$eventTable <- addTminToEventTable(eT = l$eventTable, simtime = l$simtime)
  }
  l <- l[setdiff(names(l), "opt_timeHandling")]
  
  # Carry over mc.cores as attribute and remove mc.cores from list
  l <- l[setdiff(names(l), "mc.cores")]
  attr(l, "mc.cores") <- mc.cores
  
  class(l) <- c("simIQRmodelSpec", class(l))
  
  l
}

# -------------------------------------------------------------------------#
# Helper functions ----
# -------------------------------------------------------------------------#


# Add min(simtime) to eventTable
#
# sim_IQRmodel with eventTable discards simtimes < min(simtime).
# This is a reasonable choice for PK models with trivial initial steady state, but leads to
# problems in other scenarios, e.g., when you want to pre-equilibrate QSP models.
#
# Solution:
# For IDs with min(simtime) < min(eT$TIME), add a row with TIME = min(simtime),
# and carry backward the parameters
#
# @param eT [IQRtools::IQReventTable]
# @param simtime Vector of simulation times
#
# @return eT with min(simtimes added)
# @author Daniel Lill (daniel.lill@intiquan.com)
# @md
# @family IQReventTable
# @importFrom data.table data.table rbindlist copy
# @keywords Internal
#
# @examples
addTminToEventTable <- function(eT, simtime) {
  # Nothing works without data.table
  e0 <- data.table::data.table(eT)
  minSimTimes <- sapply(simtime, min)
  dt_minST <- data.table::data.table(ID = as.numeric(names(minSimTimes)), minSimTimes)
  
  # Work on first rows per subject if they are greater than min(simtime)
  e1 <- e0[,.SD[1], by = ID]
  
  # merge min simtimes
  e2 <- merge(e1, dt_minST)
  
  e2 <- e2[TIME > minSimTimes]
  
  # Check if anything needs to be done. If not, return original
  if (nrow(e2) == 0) {
    return(eT)
  }
  
  # Copy them, replace dosings by dummy dosing at min(simtime)
  eadd <- data.table::copy(e2)
  # eadd[,`:=`(AMT = 0)]
  eadd[,`:=`(TIME = minSimTimes)]
  # eadd[,`:=`(TINF = 1e-8)] # Should be small enough to avoid the error "Event table input argument seems to have same type of event happening at an identical time point in ID=1"
  eadd[,minSimTimes:=NULL]
  # Remove regressor values from eadd in the original eventtable -> they are set in eadd
  eSetNa <- data.table::copy(e2[,minSimTimes:=NULL])
  regnm <- setdiff(names(eSetNa), c("ID", "ID.POP", "TIME", "ADM", "AMT", "TINF", "RATE", "ADDL", "II"))
  eSetNa[,(regnm) := lapply(.SD, function(x) NA), .SDcols = regnm]
  e0[eSetNa, (regnm) := mget(sprintf("i.%s", regnm)), on = c("ID", "TIME")] # "Update-join"
  
  # Rowbind dummy doses and original eT
  e <- data.table::rbindlist(list(e0,eadd), use.names = TRUE)
  e <- e[order(ID, TIME)]
  
  # Turn into eventTable
  eT <- structure(data.frame(e), class = c("IQReventTable", "data.frame"))
  eT
}


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
#'     group1 = IQRdosing(
#'       TIME = 0, ADM = 1, AMT = 40, TINF = 0
#'     ),
#'     group2 = IQRdosing(
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
  IQRdataGENERAL(design0)
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
#' fileURL <- "http://192.168.1.12:7804/IQauxfn/inst/02-ExampleModels/ModelLibraryPK/model_2cpt_linsat_iv.txt"
#' fl <- file.path(tempdir(), basename(fileURL))
#' download.file(fileURL, fl)
#' model <- IQRmodel(fl)
#'
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
#'     group1 = IQRdosing(
#'       TIME = 0, ADM = 1, AMT = 40, TINF = 0
#'     ),
#'     group2 = IQRdosing(
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
#' gpf <- GPF(estimates = estimates)
#' indPars <- sample_GPF(gpf, Nsamples = length(unique(dataGEN$ID)), FLAGid = TRUE)
#' indPars <- indPars$indParamValues
#' dataGEN <- merge(dataGEN, indPars, by = "ID")
#' dataGEN <- IQRdataGENERAL(dataGEN)
#'
#' # Step 3: Populate dataGEN
#' dataGENsim <- populateGDFWithSimValues(dataGEN, model, regression = setdiff(names(indPars), "ID"), abs0inputs = NULL, abs0Tk0param = NULL)
#'
#' # Visualize
#' IQRggplot(dataGENsim, aes(TIME, DV, group = USUBJID)) +
#'   facet_wrap(~TRTNAME, scales = "free") +
#'   geom_line() +
#'   labs(title = "Simulated data (no residual error)")
populateGDFWithSimValues_consistent <- function(dataGEN, model, regression = NULL, abs0inputs = NULL, abs0Tk0param = NULL) {
  # For preserving attributes of dataGEN. Nothing should change except for the values in DV and VALUE,
  # so we can just take it as it is and replace the new attributes later.
  atr <- attributes(dataGEN)

  
  
  # Simulate with individual parameters and melt simres.
  eT <- IQReventTable(data = dataGEN,
                      regression = regression,
                      abs0inputs = abs0inputs,
                      abs0Tk0param = abs0Tk0param)
  
  # simres <- sim_IQRmodel(model, simtime = sort(unique(c(dataGEN$TIME))),
  #                        eventTable = eT, FLAGoutputsOnly = TRUE)
  # use sim_IQRmodelConsistent
  simIQRmodelSpec1 <- createSimIQRmodelSpec(model, simtime = sort(unique(c(dataGEN$TIME))),
                                            eventTable = eT, FLAGoutputsOnly = TRUE)
  simres <- sim_IQRmodelConsistent(simIQRmodelSpec = simIQRmodelSpec1)
  
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
  dataGENNew <- IQRdataGENERAL(dataGENNew)
  attributes(dataGENNew) <- atr
  dataGENNew
}