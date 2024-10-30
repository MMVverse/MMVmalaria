#' getPopIndPred
#'
#' @description
#' @param modelFolder
#' @param timeGrid Default: 'obsTimes'
#' @param abs0inputs Default: NULL
#' @param abs0Tk0param Default: NULL
#' @return
#' @author Anne Kümmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family Experimental
#' @importFrom dplyr full_join
getPopIndPred <- function(modelFolder, timeGrid = "obsTimes", abs0inputs = NULL, abs0Tk0param = NULL) {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check input arguments ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (!(timeGrid %in% c("fine","obsTimes"))) {
    stop("timeGrid input has to be 'fine' or 'obsTimes'.")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get project information and model results ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  modelRes <- getResults_IQRnlmeProjectMulti(as_IQRnlmeProjectMulti(modelFolder))[[1]]

  # Get the dataset for subject and compound information:
  data     <- IQRloadCSVdata(file.path(modelFolder,modelRes$projectHeader$DATA))
  subjInfo <- unique(data[,c("STUDY","USUBJID", "ID", "TRTNAME")])
  obsInfo  <- unique(data[data$EVID==0, c("YTYPE","NAME", "UNIT")])
  nsubj    <- dim(subjInfo)[1]

  # Get regression parameters if present
  regressionParameters <- modelRes$projectHeader$REGRESSIONNAMES

  # Get individual parameters
  indivParameters <- getIndivParameters_IQRnlmeProject(modelFolder)

  # Get information on covariates (needed to potential adjust population parameters)
  if (modelRes$projectHeader$COVARIATESUSED == "") {
    Covariates <- NULL
  } else {
    Covariates <- data[, c("ID", modelRes$projectHeader$COVARIATESUSED)]
  }

  # Get population parameters (that account for covariates)
  if (nsubj > 1) {
    popSample <- sample_IQRnlmeProject(modelFolder, Nsamples = nsubj, FLAG_SAMPLE = 3, covariates= Covariates)
  } else {
    popSample <- sample_IQRnlmeProject(modelFolder, Nsamples = 2, FLAG_SAMPLE = 3, covariates= Covariates)
  }
  popParameters <- cbind(subjInfo[,"ID", drop = FALSE], popSample$parameterValuesPopulation[1:nsubj,])

  # Names of model parameters
  modelParNames <- modelRes$rawParameterInfo$fixedEffects$names

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare info for simulations ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load model
  model <- IQRmodel(file.path(modelFolder, "model.txt"))

  # Define simulation times (either provided by user or )
  if (timeGrid == "fine") {
    simTimes <- lapply(unique(data$ID), function(ID) {
      timesID <- data$TIME[data$ID == ID]
      simtimeID <- seq(floor(min(timesID)), ceiling(max(timesID)))
      simtimeID
    })
  } else {
    simTimes <- lapply(unique(data$ID), function(ID) {
      timesID <- unique(data$TIME[data$ID == ID & data$EVID == 0])
      timesID
    })
  }
  names(simTimes) <- unique(data$ID)

  # Table with dosing and regression parameters
  eventData0 <- data[data$EVID == 1,c("ID", "TIME", "AMT", "ADM", "RATE", "II", "ADDL", "TINF", regressionParameters)]

  # Handle negative or observations before first dosing time
  FirstTime <- ddply(data, ~ID, function(x) c(FTIME = min(x$TIME)))
  # eventData0 <- merge(eventData0, FirstTime, by = "ID")
  eventData0 <- ddply(eventData0, ~ID, function(xx) {
    fTime <- FirstTime$FTIME[FirstTime$ID == xx$ID[1]]
    if (min(xx$TIME) > fTime) {
      xx <- xx[c(1,1:dim(xx)[1]),]
      xx$TIME[1] <- fTime
      xx$AMT[1] <- 0
      xx$ADDL[1] <- 0
      xx$II[1] <- 0
    }
    xx
  })
  # # Get information on individual time shift
  # ShiftTime <- plyr::ddply(eventData0, ~ID, function(xx) c(TSHIFT = xx$TIME[1] - xx$TIME0[1]))

  # Add population or individual model parameters
  eventDataPop <- merge(eventData0, popParameters)
  eventDataInd <- merge(eventData0, indivParameters)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Population predictions ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  eventTablePop <- IQReventTable(eventDataPop, regression = c(regressionParameters, modelParNames),
                                 abs0inputs = abs0inputs, abs0Tk0param = abs0Tk0param)
  simPop <- sim_IQRmodel(model, simtime = simTimes, eventTable = eventTablePop, FLAGoutputsOnly = TRUE)

  simPop <- gather(simPop, "YTYPE", "XPRED", grep("OUTPUT",names(simPop)))
  simPop$YTYPE <- as.numeric(gsub("OUTPUT","",simPop$YTYPE))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Individual predictions ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  eventTableInd <- IQReventTable(eventDataInd, regression = c(regressionParameters, modelParNames),
                                 abs0inputs = abs0inputs, abs0Tk0param = abs0Tk0param)
  simInd <- sim_IQRmodel(model, simtime = simTimes, eventTable = eventTableInd, FLAGoutputsOnly = TRUE)

  simInd <- gather(simInd, "YTYPE", "IPRED", grep("OUTPUT",names(simInd)))
  simInd$YTYPE <- as.numeric(gsub("OUTPUT","",simInd$YTYPE))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  simRes <- dplyr::full_join(simPop, simInd)

  return(simRes)

}

