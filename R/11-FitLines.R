#' define_tLLOQ
#' Defines the time when parasite clearance or growth crosses the LLOQ
#'
#' Determines the time when the LLOQ is reached based on linear segments
#' fitted to clearance phase (output from [fit_LinearSegments]) or re-growth
#' phase (output from [fit_LinearModel]). Needed to determine the time range within
#' which the MIC is reached assuming that it needs to in the time
#' when parasites are below LLOQ
#'
#' Input needs to contain Phase column indicating whether data is from
#' parasite clearance phase (Phase = "KILL") or from re-growth (Phase = "RECRUD").
#' Segment knots are defined by Xknot and Yknot columns
#'
#' @md
#'
#' @param xx data frame with parasitemia data
#' @param LLOQ LLOQ for parasite detection (log scale)
#'
#' @return input data frame amended by row indicating the time at which LLOQ is reached
#'         and Phase column indicating whether this is defining the lower or upper bound for the MIC.
#'         (MICmin for re-growth phase and MICmax for clearance phase)
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan)
define_tLLOQ <- function(xx, LLOQ) {


  Phase <- xx$Phase[1]
  if (!(Phase %in% c("KILL","RECRUD"))) stop("Phase in input needs to be 'KILL' or 'RECRUD'")

  if (Phase == "KILL") {
    xxx    <- xx[xx$Nknots[1]+c(-1,0),]
    PhaseXXX <- "MICmax"
  } else {
    xxx    <- xx[1:2,]
    PhaseXXX <- "MICmin"
  }

  a <- (xxx$Yknot[1]-xxx$Yknot[2])/(xxx$Xknot[1]-xxx$Xknot[2])
  b <- xxx$Yknot[1] - a*xxx$Xknot[1]
  tLLOQ <- (LLOQ - b)/a

  add <- data.frame(USUBJID=xx$USUBJID[1], Phase = PhaseXXX, Xknot = tLLOQ, Yknot = LLOQ)
  out <- rbind(xx[,c("USUBJID", "Phase", "Xknot", "Yknot")], add)
  out
}
#' define_tMIClinearSegments
#' Determine time of MIC based on fitted linear segments
#'
#' Determines intersection of last segment of clearance phase and line of re-growth
#' phase as the time of MIC.
#'
#' Input data frame contains the outputs of [fit_LinearModel]
#' and [fit_LinearSegments] and contain coordinates of knots of linear segments.
#' Phases are differentiated by 'Phase' column with KILL for clearance phase
#' and RECRUD for re-growth phase.
#'
#' @md
#'
#' @param FIT data frame defining linear segments fitted
#'         to clearance and regrowth data of individual
#'         (combined output of [fit_LinearModel] and [fit_LinearSegments]).
#'         Phase column needs to be defined (see Details)
#'
#' @return data frame with last knot of clearance phase,
#'         first knot of regrowth phase
#'         and interpolated time of MIC defined by Xknot, Yknot
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
define_tMIClinearSegments <- function(FIT)
{

  # ---
  # KILLing phase:
  # Take last segment and determine line to extrapolate
  FITkill <- subset(FIT, Phase == "KILL")
  FITkill <- FITkill[order(FITkill$Xknot),]
  FITkill <- FITkill[dim(FITkill)[1]+c(-1,0),]

  lineparKILL <- get_linepar(FITkill)

  # ---
  # RECRUDescence phase
  # Take first segment and determine line to extrapolate
  FITrecrud <- subset(FIT, Phase == "RECRUD")
  FITrecrud <- FITrecrud[order(FITrecrud$Xknot),]
  FITrecrud <- FITrecrud[c(1,2),]
  # browser()
  lineparRECRUD <- get_linepar(FITrecrud)

  # ---
  # Determine time at which MIC is reached
  tMIC  = (lineparRECRUD[["intercept"]]-lineparKILL[["intercept"]])/(lineparKILL[["slope"]]-lineparRECRUD[["slope"]])
  nadir = lineparRECRUD[["intercept"]] + lineparRECRUD[["slope"]] * tMIC
  EXTRAPOL <- rbind(
    FITkill[2,c("Phase","Xknot","Yknot")],
    data.frame(Phase = "MIC", Xknot=tMIC, Yknot=nadir),
    FITrecrud[1,c("Phase","Xknot","Yknot")]
  )

  return(EXTRAPOL = EXTRAPOL)

}
#' fit_LinearModel
#' Fits line to parasitemia data on log scale
#'
#' Fits a line using `lm` to parasitemia data (TIME and VALUE columns).
#'
#' This function is used in the model-free assessment of PD data by describing the
#' PD data by linear segments. The linear segments to describe the PD data
#' are characterized by knots which are the endpoint of the single segments.
#' This function is used to fit single line to growth data.
#' The output has the same format as fit_LinearSegments and gives the number
#' of knots, here 2, and the position of the knots, the sum of squares
#' for the fit and some status of the fit.
#' Knots are given by Xknot (time scale) and Yknot (logtransformed PD scale).
#' Uses only non-censored data.
#'
#'
#' @param x data frame with parasitemia records
#' @param FLAGlogtransY Flag whether to need to be log-transformed
#'                      (TRUE (default) in case data are provided on normal scale)
#'
#' @return data frame with locations of the knots (Xknot and Yknot), number of knots (Nknots),
#'         the sum of squares of the fit (SSresid) and fit information (fitinfo)
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
fit_LinearModel <- function(x, FLAGlogtransY = TRUE)
{
  # Get timepoints and values for uncensored data
  idxCens <-  x$CENS == 1
  x = x[!idxCens,]

  # Do fit only if there is more than 1 data point
  if (dim(x)[1] > 1) {
    if (FLAGlogtransY) x$VALUE = log(x$VALUE)

    # Make sure that ordered by TIME
    x <- x[order(x$TIME),]

    # Fit line
    linmod <- lm(VALUE~TIME, data=x)
    linpred <- predict(linmod, x[c(1,dim(x)[1]),])

    # Prepare output
    fit <- data.frame(
      Nknots      = 2,
      Xknot       = x[c(1,dim(x)[1]),"TIME"],
      Yknot       = linpred,
      SSresid     = sum(residuals(linmod)^2),
      fitinfo     = "Linear model fitted"
    )
  } else {
    fit <- data.frame(
      Nknots      = NA,
      Xknot       = NA,
      Yknot       = NA,
      SSresid     = NA,
      fitinfo     = "No line since only single data point."
    )
  }

  return(fit)
}


#' fit_LinearSegments
#' Fits line to parasitemia data on log scale
#'
#' Fits linear segments to parasitemia data (TIME and VALUE columns)
#' using `cobs` function from the cobs package.
#'
#' This function is used in the model-free assessment of PD data by describing the
#' PD data by linear segments. The linear segments to describe the PD data
#' are characterized by knots which are the endpoint of the single segments.
#' This function is used to data in the clearance phase.
#' The output has the same format as fit_LinearModel and gives the number
#' of knots, the position of the knots, the sum of squares
#' for the fit, and some status of the fit.
#' If there is only one observation, no fit is done. If there are two observations,
#' the line through these points is taken. More than two observations are fitted
#' using the cobs function (constrained linear splines). If cobs fit fails a line
#' is fitted trough first and last timepoint.
#' Knots are given by Xknot (time scale) and Yknot (logtransformed PD scale).
#' Uses only non-censored data.
#'
#' @param data data frame with parasitemia records
#' @param constraint character (string) specifying the kind of constraint
#'            (recommended: "none" or "decrease"). Decrease will only allow negatve slopes for segments.
#' @param Nsegments maximum number of knots, defaults to 3
#' @param FLAGlogtransY Flag whether to need to be log-transformed
#'                      (TRUE (default) in case data are provided on normal scale)
#'
#' @return data frame with locations of the knots (Xknot and Yknot), number of knots (Nknots),
#'         the sum of squares of the fit (SSresid) and fit information (fitinfo)
#' @export
#'
#' @seealso [cobs::cobs]
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
fit_LinearSegments <- function(data, constraint = "none", Nsegments = 3, FLAGlogtransY = TRUE)
{
  # Get timepoints and values for uncensored data
  idxCens <- data$CENS==1
  x <- data$TIME[!idxCens]
  if(FLAGlogtransY){
    y <- log(data$VALUE[!idxCens])
  }else{
    y <- data$VALUE[!idxCens]
  }

  # Get censoring information:
  #   NOTE: Take DV column as censoring value is set here
  censored <- cbind(-data$CENS[idxCens],data$NT[idxCens],data$DV[idxCens])

  if (length(x) == 1) {
    message("Only one datapoint.")
    fit <- data.frame(
      Nknots      = NA,
      Xknot       = NA,
      Yknot       = NA,
      SSresid     = NA,
      fitinfo     = "No line since only single data point",
      stringsAsFactors = FALSE)
  } else {

    if (length(x) == 2) {
      message("Only two datapoints. Extrapolate with line through these two.")
      fit <- data.frame(
        Nknots      = 2,
        Xknot       = x,
        Yknot       = y,
        SSresid     = 0,
        fitinfo     = "Line through two data points.",
        stringsAsFactors = FALSE)
    } else {
      message("Fitting constrained B-splines.")
      if (is.null(Nsegments)) {
        res <- suppressWarnings(try(cobs::cobs(x, y, constraint = constraint, pointwise = censored, degree = 1)))
      } else {
        res <- suppressWarnings(try(cobs::cobs(x, y, constraint = constraint, pointwise = censored, degree = 1, nknots = Nsegments+1)))
      }


      if (class(res) == "cobs"){
        message("... success!")
        fit <- data.frame(
          Nknots      = length(knots(res)),
          Xknot       = knots(res),
          Yknot       = res$coef,
          SSresid     = sum(res$resid^2),
          fitinfo     = "cobs fit",
          stringsAsFactors = FALSE)
      } else {
        message("Failed. Extrapolate with first and last.")

        fit <- data.frame(
          Nknots      = 2,
          Xknot       = x[1:length(x)],
          Yknot       = y[1:length(x)],
          SSresid     = 0,
          fitinfo     = "Line through first and last data point (cobs fit failed).",
          stringsAsFactors = FALSE)
      }

    }
  }
  return(fit)
}


#' get_CLmaxTlag
#' Determines maximum clearance and lag time of parasitemia
#'
#' Derives maximum rate of clearance and lag time before it is achieved
#' based on segments fitted to PD data in clearance phase.
#'
#' The input is a data frame with knots of the (at most 2) segments for an individual as
#' returned by fit_LinearSegments.
#' If only one segment present, lag time is zero and maximum clearance the slope of this segment.
#' If there are two segments, maximum clearance is the steeper slope of the two.
#' Tlag is the time of the middle knot if the first slope is smaller than 20% of the second and 0 otherwise.
#'
#' @md
#'
#' @param x data frame with one or two segments characterizing parasite clearance (output of [fit_LinearSegments])
#'
#' @return input data frame amended by CLmax and Tlag columns
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan)
get_CLmaxTlag <- function(x) {
  if(is.na(x$Nknots[1]) | is.na(x$Xknot[1]) | is.na(x$Yknot[1])){
    # No information
    x$CLmax <- NA
    x$Tlag  <- NA
  }else{
    if (x$Nknots[1]==2) {
      # Only single segment: take slope as max clearance and no Tlag
      x$CLmax <- -diff(x$Yknot)/diff(x$Xknot)
      x$Tlag  <- 0
    }else{
      # Two segments: Take larger slope as CLmax, and Tlag as middle timepoint
      # if first slope is smaller than 20% of second, otherwise Tlag = 0
      slopes  <- -diff(x$Yknot)/diff(x$Xknot)
      x$CLmax <- max(slopes)
      if(slopes[1] < 0.2*slopes[2]){
        x$Tlag <- x$Xknot[2]
      }else{
        x$Tlag <- 0
      }
    }
  }
  x
}
#' get_DataClearance
#' Gets PD subset from the parasite clearance phase
#'
#' Subsets parasitemia data to initial clearance phase
#' after dosing, i.e., until minimum in parasite count
#' or first occurence of LLOQ data point
#'
#' Expects that given data starts as drug administration
#' (no pre-dose growth data to be present).
#' Expects TIME, VALUE, and LLOQ columns.
#'
#' @param x data frame with parasitemia records for an individual
#'
#' @return subsetted input data frame
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan)
get_DataClearance <- function(x) {
  #
  # x is PD data from an indivdual
  x <- x[order(x$TIME),]
  x <- x[!is.na(x$VALUE),]
  if (dim(x)[1] == 0) return()
  idxLLOQ <- x$VALUE < x$LLOQ
  if (any(idxLLOQ)) {
    if (min(which(idxLLOQ)) == 1) {return()}
    out <- x[1:min(which(idxLLOQ)),]
  } else {
    idxMIN <- which.min(x$VALUE)
    if (idxMIN == 1) {return()}
    out <- x[1:idxMIN,]
  }
  out
}
#' get_DataGrowth
#' Gets PD subset from the parasite recrudescence phase
#'
#' Subsets parasitemia data to re-growth or recrudescence phase
#' , i.e., after minimum in parasite count
#' or last occurence of LLOQ data point
#'
#' Expects TIME, VALUE, and LLOQ columns.
#'
#' @param x data frame with parasitemia records for an individual
#'
#' @return subsetted input data frame
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan)
get_DataGrowth <- function(x) {
  #
  # x is PD data from an indivdual
  x <- x[order(x$TIME),]
  x <- x[!is.na(x$VALUE),]
  if (dim(x)[1] == 0) return()
  idxLLOQ <- x$VALUE < x$LLOQ
  if (any(idxLLOQ)) {
    if (max(which(idxLLOQ)) == dim(x)[1]) {return()}
    out <- x[max(which(idxLLOQ)):dim(x)[1],]
  } else {
    idxMIN <- which.min(x$VALUE)
    if (idxMIN == dim(x)[1]) {return()}
    out <- x[idxMIN:dim(x)[1],]
  }
  out
}
#' get_Growth
#' Determines growth rate from line through growth data
#'
#' Determines growth rate based on line fitted to PD data in growth phase.
#' As the line is on logscale the growth rate is the slope of the line.
#'
#' The input is a data frame with knots of the line for an individual as
#' returned by fit_LinearModel
#'
#' @md
#'
#' @param x data frame with two knots characterizing growth (output of [fit_LinearModel])
#'
#' @return input data frame amended by GR column with growth rate value
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan)
get_Growth <- function(x) {
  if (is.na(x$Nknots[1])) {
    # No information
    x$GR <- NA
  } else {
    x$GR <- diff(x$Yknot)/diff(x$Xknot)
  }
  x
}
#' Slope and intercept of line
#'
#' `get_linepar` Calculates the slope and intercept of a line defined by two points
#'
#' Given the coordinates of two points, the slope and the intercept
#' for the line through these points is determined.
#'
#' Each row of the data frame defines one point.
#'
#' @param twoknots data frame with Xknot and Yknot colums and two rows
#'
#' @return A name vector with intercept and slope
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
get_linepar <- function(twoknots)
{
  if (dim(twoknots)[1] != 2) stop("Number of rows of input is not 2.")

  # Difference:
  Dx <- diff(twoknots$Xknot)
  Dy <- diff(twoknots$Yknot)

  # Slope and Intercept:
  slope     <- Dy/Dx
  intercept <- twoknots$Yknot[1] - slope*twoknots$Xknot[1]

  # Output:
  return(c(intercept=intercept, slope=slope))
}
#' get_MICextrapolPK
#' Determine MIC based on linearily inter/extrapolated PK
#'
#' The MIC is determined as the concentration extrapolated
#' linearily based on the PK data at a given time of MIC
#' which was determined based on linear PD data interpolation.
#'
#' Linear segments are fitted to the individual PK data by [fit_LinSegments].
#' The MIC is taken as the value of the respective segment at the given time of MIC in 'tMICindiv'
#' If the time of MIC is later than the observed data, i.e., later than the time described by the segment,
#' the last segment is extrapolated.
#' Concentration data is expected to be given in ug/mL.
#' Note that only single dose data is handled appropriately.
#'
#' @md
#'
#' @param tMICindiv data frame with individual times of MIC (Xknot column)
#' @param dataPK data frame with PK data (expects columns USUBJID, TIME, VALUE, LLOQ, CENS, TIMEUNIT
#'              as provided by general dataset format in MMVmalaria or IQRtools)
#' @param plottitle character (string) to be printed as figure title
#'
#' @return list with data frame of individual MIC values and a plot (ggplot object) illustrating the MIC determination
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
get_MICextrapolPK <- function(tMICindiv,
                              dataPK,
                              plottitle = "MIC determined by linear extrapolation of parasitemia and PK data")
{

  # Get linear segments for extra/intrapolation
  LinSegPK <- plyr::ddply(dataPK, ~USUBJID+TRT+TRTNAME+STUDY, fit_LinearSegments, constraint = "none", Nsegments = NULL, FLAGlogtransY = TRUE)

  # Determine PK at tMIC
  LinSegPK <- merge(LinSegPK, tMICindiv[,c("USUBJID","Xknot")], by = c("USUBJID"), suffixes = c("","MIC"))
  MICextra <- plyr::ddply(LinSegPK, ~USUBJID+TRT+TRTNAME+STUDY, function(x) {
    x <- x[!is.na(x$Xknot),]
    if (dim(x)[1] > 1){
      x <- x[order(x$Xknot),]
      if (all(x$Xknot < x$XknotMIC)) { # do extrapolation
        support <- x[dim(x)[1]+c(-1,0),]
      } else { # interpolation between to observations
        support <- x[max(which(x$Xknot < x$XknotMIC))+c(0,1),]
      }
      YknotMIC <- support$Yknot[1] + (support$Yknot[2]-support$Yknot[1])/(support$Xknot[2]-support$Xknot[1]) * (support$XknotMIC[1]-support$Xknot[1])

      return(data.frame(
        XknotS   = support$Xknot[1],
        YknotS   = support$Yknot[1],
        XknotMIC = support$XknotMIC[1],
        YknotMIC = YknotMIC,
        timeMIC  = support$XknotMIC[1],
        MIC      = exp(YknotMIC)
      ))} else {
        NULL
      }
  })
  MICextra <- within(MICextra, {
    MICstr = paste0(sprintf("%.5f",MIC*1000), " ng/mL")
  } )

  # plot PK data and simulation
  xlabel <- paste0("Time (", dataPK$TIMEUNIT[1],")")
  ylabel <- paste0("Concentration (ng/mL)")

  gr <- IQRggplot(LinSegPK,aes(Xknot,exp(Yknot)*1000)) +
    # Extrapolation
    geom_line(size=0.8) +
    geom_segment(data = MICextra, aes(x=XknotS, y=exp(YknotS)*1000, xend=XknotMIC, yend=MIC*1000), linetype = "dotted", size = 0.6) +
    # Data
    geom_point(data=dataPK, aes(TIME,VALUE*1000, fill = as.character(CENS)), shape=21, size=3, color="white") +
    # Info on MIC
    geom_text(data=MICextra,
              aes(XknotMIC,MIC*1000, label=MICstr),
              hjust=1.1,vjust=-0.2, size = 3) +
    geom_hline(data=MICextra, aes(yintercept=MIC*1000), linetype = 2) +
    geom_vline(data=MICextra, aes(xintercept=XknotMIC), linetype = 2) +
    # Polising
    scale_fill_manual("BLQ value", values = c("navyblue","firebrick")) +
    scale_color_manual(values=c("olivedrab","firebrick","firebrick"),guide = FALSE) +
    scale_y_log10() +
    labs(x=xlabel, y=ylabel, title = plottitle) +
    facet_wrap(~USUBJID, scales="free") +
    theme(legend.position = "bottom")

  return(list(MICextra,gr))

}


#' get_MICsimPK
#' Determine MIC based on simulated PK
#'
#' The MIC is determined as the concentration simulated
#' at a given time of MIC which was determined based on
#' linear PD data interpolation.
#'
#' Individual PK timexourses are simulated by given model and individual PK parameters.
#' The MIC is taken as the simulated value given the time of MIC in 'tMICindiv'
#' Concentration data is expected to be given in ug/mL.
#'
#' @param dataDOS data frame with individual dosing information (defined by TIME, AMT and RATE)
#'                and individual PK parameters.
#'                Expects columns USUBJID, TRT, TRTNAME, and STUDY to identify and annotate individuals
#' @param model IQRmodel of the PK model to use for simulation
#' @param tMICindiv data frame with timepoints of MIC (Xknot) for individuals (USUBJID)
#' @param dataPK data frame with PK data (expects columns USUBJID, TIME, VALUE, LLOQ, CENS, TIMEUNIT
#'              as provided by general dataset format in MMVmalaria or IQRtools)
#' @param plottitle character (string) to be printed as figure title
#'
#' @return list with data frame of individual MIC values and a plot (ggplot object) illustrating the MIC determination
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)

get_MICsimPK <- function(dataDOS,
                         model,
                         tMICindiv,
                         dataPK,
                         plottitle = "MIC determined by linear extrapolation of parasitemia data and PK simulation")
{


  # Get MIC time to dosing data:
  dataDOS <- merge(dataDOS, tMICindiv[,c("USUBJID","Xknot")])

  # Get data time range and merge to dosing:
  timeRange <- ddply(dataPK, ~USUBJID, summarize, minTIME = min(TIME), maxTIME = max(TIME))
  dataDOS   <- merge(dataDOS, timeRange)


  # Simulate PK at MIC for each individual:
  simres <- ddply(dataDOS, ~USUBJID+TRT+TRTNAME+STUDY, simulate_indivPK, model=model)
  MICsim <- subset(simres, FLAGmic == 1, c("USUBJID","TRT","TRTNAME","STUDY","TIME","OUTPUT1"))
  names(MICsim) <- sub("OUTPUT1","MIC",names(MICsim))
  names(MICsim) <- sub("TIME","timeMIC",names(MICsim))
  MICsim <- within(MICsim, {
    MICstr = paste0(sprintf("%.5f",MIC*1000), " ng/mL")
  } )

  # plot PK data and simulation:
  dataPK$CENS <- (dataPK$VALUE < dataPK$LLOQ)
  xlabel <- paste0("Time (", dataPK$TIMEUNIT[1],")")
  ylabel <- paste0("Concentration (ng/mL)")

  gr <- IQRggplot(simres,aes(TIME,OUTPUT1*1000)) +
    # Simulation
    geom_line(size=0.8) +
    # Data
    geom_point(data = dataPK, aes(TIME, VALUE*1000, fill=CENS), shape=21, size=3, color="white") +
    # Info on MIC
    geom_text(data=MICsim,
              aes(timeMIC,MIC*1000, label=MICstr),
              hjust=1.1,vjust=-0.2, size=3) +
    geom_hline(data=MICsim, aes(yintercept=MIC*1000), linetype = 2) +
    geom_vline(data=MICsim, aes(xintercept=timeMIC), linetype = 2) +
    # Polising
    scale_fill_manual("BLQ value", values = c("navyblue","firebrick")) +
    scale_color_manual(values=c("olivedrab","firebrick","firebrick"),guide = FALSE) +
    scale_y_log10() +
    labs(x=xlabel, y=ylabel, title = plottitle) +
    facet_wrap(~USUBJID, scales="free") +
    theme(legend.position = "bottom")

  return(list(MICsim,gr))
}


#' plot_LinSegments
#' Plots linear segments fitted to clearance and recrudescence phase of PD data
#'
#' Parasitemia data is overlaid with linear segments for the clearance phase
#' and the (re-)growth phase separately.
#'
#' fitcl is the output of fit_LinearSegments and fitgr is the output of fit_LinearModel
#' Parasitemia data is assumed to be on linear scale and log-transformed before plotting.
#'
#' @param x   list with information for an individual with the following fields:
#'            data: PD data (linear scale)
#'            fitcl: segments (knots) fitted to clearance phase ().
#'            fitgr: segments (knots) fitted to growth phase.
#'
#' @return ggplot object of plot
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
plot_LinSegments <- function(x) {
  IQRggplot(x$data, aes(TIME, VALUE)) +
    geom_point(aes(shape = CENS==1)) +
    # add linear fits for growth
    geom_line(data=x$fitgr, aes(x=Xknot,y=exp(Yknot)), color = "olivedrab") +
    geom_text(data=x$fitgr, x=max(x$fitgr$Xknot),y=exp(max(x$fitgr$Yknot)), aes(label = paste0("Growth rate: ", round(GR,3), " 1/h")), color = "olivedrab", hjust=1, vjust=0.5) +
    # add linear segments for clearance
    geom_line(data=x$fitcl, aes(x=Xknot,y=exp(Yknot)), color = "firebrick") +
    geom_text(data=x$fitcl, x=min(x$fitcl$Xknot),y=exp(max(x$fitcl$Yknot)), aes(label = paste0("Max. clearance: ", round(CLmax,3), " 1/h\nLag time: ", round(Tlag,2), " h")), color = "firebrick", hjust=0, vjust=0.5) +

    scale_shape_manual("BLQ", values = c(19,4)) +
    scale_y_log10(breaks = 10^seq(-3,8), labels = 10^seq(-3,8)) +
    labs(x = "Time [hours]", y = "Parasitemia [%]", title = paste0(x$data$USUBJID[1], " in group ", x$data$TRTNAME[1]))
}
#' plot_tMICextrapolation
#' Visualize determination of time of MIC by linear segments
#'
#' Plots parasitemia data overlaid with linear segments fitted to
#' clearance and re-growth phase and the interpolation of the
#' time of MIC.
#'
#' Segments are expected to be defined on logtransformed parasitemia data.
#' Data may need to be transformed if given on linear scale by setting
#' FLAGlogtransData accordingly.
#' Phases are differentiated by 'Phase' column with KILL for clearance phase,
#' RECRUD for re-growth phase, and "MIC" for the interpolation between them.
#' Individuals are plotted to separate panels.
#'
#' @md
#'
#' @param EXTRA Data frame defining the  that define the time of MIC
#'         (output of [define_tMIClinearSegments]).
#'         Phase column need to be defined (see Details).
#' @param FIT Data frame defining linear segments fitted
#'         to clearance and re-growth data
#'         (combined output of [fit_LinearModel] and [fit_LinearSegments]).
#'         Phase column need to be defined (see Details).
#' @param data Data frame with parasitemia data (expects NT, DV, TIME and CENS columns).
#' @param xlabel Character (string) to be printed on x axis.
#' @param ylabel Character (string) to be printed on y axis.
#' @param title Character (string) to be printed as figure title.
#' @param FLAGlogtransData Flag whether data needs to be logtransformed before plotting.
#'
#' @return ggplot object
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @importFrom MMVbase get_ActivityPath
plot_tMICextrapolation <- function(EXTRA, FIT, data,
                                   xlabel = "x",
                                   ylabel = "y",
                                   title  = "Determination of time at MIC by linear extrapolation",
                                   FLAGlogtransData = TRUE,
                                   ActivityPath     = NULL)
{
  # get some dimensions
  dataxrange <- range(data$NT, na.rm = TRUE)
  datayrange <- range(data$DV, na.rm = TRUE)

  # log-transform data if required
  if (FLAGlogtransData) data$DV <- log(data$DV)


  gr <- IQRggplot(FIT, aes(Xknot,Yknot)) +
    # fitted lines and indicate knots
    geom_line(aes(group=Phase),size = 0.8) +
    geom_vline(aes(xintercept=Xknot), linetype = 3, size = 0.6) +
    # add extrapolated lines
    geom_line(data = EXTRA, linetype = 2, size = 0.8) +
    geom_vline(data = subset(EXTRA, Phase == "MIC"), aes(xintercept=Xknot), color = "olivedrab") +
    geom_text(data = subset(EXTRA, Phase == "MIC"),
              aes(label = paste0("time at MIC: ",sprintf("%.1f", Xknot))),
              color = "olivedrab", nudge_x = 0.1*diff(dataxrange), hjust = 0) +
    # add data
    geom_point(data=data, aes(TIME,DV, color = as.logical(CENS)), shape = 1, size = 3) +
    # info on fitting
    geom_text(
      data=subset(FIT, Phase == "KILL"),
      x=dataxrange[1],y=datayrange[2]+0.1*diff(datayrange),
      aes(label=paste0(' (',fitinfo,')')), hjust = 0, vjust=1
    ) +
    geom_text(
      data=subset(FIT, Phase == "RECRUD"),
      x=dataxrange[2],y=datayrange[2]+0.1*diff(datayrange),
      aes(label=paste0('(',fitinfo,')')), hjust = 1, vjust=1
    ) +
    scale_color_manual("BLQ value", values = c("navyblue","firebrick")) +
    facet_wrap(~USUBJID, scales = "free") +
    labs(x=xlabel,
         y=ylabel,
         title=title,
         caption = paste0("Activity: ", MMVbase::get_ActivityPath(ActivityPath))) +
    theme(legend.position = "bottom",
          plot.caption    = element_text(hjust=0))

  return(gr)
}
#' plot_tMICextrapolationRange
#' Visualize determination of possible time range for MIC
#'
#' Plots parasitemia data overlaid with linear segments fitted to
#' clearance and re-growth phase and the extrapolation to the times
#' at which the LLOQ is crosses while killing or re-growing of parasites
#'
#' Segments are expected to be defined on logtransformed parasitemia data.
#' Data may need to be transformed if given on linear scale by setting
#' FLAGlogtransData accordingly.
#' Phases are differentiated by 'Phase' column with KILL for clearance phase,
#' RECRUD for re-growth phase, and "MICmax" for the extrapolation of
#' the clearance phase to the LLOQ that is an upper limit for the MIC and
#' "MICmin" for the extrapolation of the re-growth phase to the LLOQ that
#' is a lower limit for the MIC.
#' Individuals are plotted to separate panels.
#'
#' @md
#'
#' @param EXTRA Data frame defining intersect that define the time of MIC
#'         (combined output of [define_tLLOQ] applied to clearance and regrowth data).
#'         Phase column need to be defined (see Details).
#' @param FIT Data frame defining linear segments fitted
#'         to clearance and re-growth data
#'         (combined output of [fit_LinearModel] and [fit_LinearSegments]).
#'         Phase column need to be defined (see Details).
#' @param data Data frame with parasitemia data (expects NT, DV, TIME and CENS columns).
#' @param xlabel Character (string) to be printed on x axis.
#' @param ylabel Character (string) to be printed on y axis.
#' @param title Character (string) to be printed as figure title.
#' @param FLAGlogtransData Flag whether data needs to be logtransformed before plotting.
#'
#' @return ggplot object
#' @export
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
plot_tMICextrapolationRange <- function(EXTRA, FIT, data,
                                        xlabel = "x",
                                        ylabel = "y",
                                        title  = "Determination of time at MIC min/max by linear extrapolation",
                                        FLAGlogtransData = TRUE,
                                        ActivityPath     = NULL)
{
  # Get some dimensions:
  dataxrange <- range(data$NT, na.rm = TRUE)
  datayrange <- range(data$DV, na.rm = TRUE)

  # log-transform data if required:
  if (FLAGlogtransData) data$DV <- log(data$DV)


  gr <- IQRggplot(FIT, aes(Xknot,Yknot)) +
    # fitted lines and indicate knots
    geom_line(aes(group=Phase),size = 0.8) +
    geom_vline(aes(xintercept=Xknot), linetype = 3, size = 0.6) +
    # add extrapolated lines
    geom_line(data = EXTRA, linetype = 2, size = 0.8) +
    geom_vline(data = subset(EXTRA, Phase %in% c("MICmin","MICmax")), aes(xintercept=Xknot), color = "olivedrab") +
    geom_text(data = subset(EXTRA, Phase == "MICmin"),
              aes(label = paste0("time at minimum MIC: ",sprintf("%.1f", Xknot))),
              color = "olivedrab", nudge_x = -0.01*diff(dataxrange),
              nudge_y = 1, hjust = 1, size = 3) +
    geom_text(data = subset(EXTRA, Phase == "MICmax"),
              aes(label = paste0("time at maximum MIC: ",sprintf("%.1f", Xknot))),
              color = "olivedrab", nudge_x = 0.01*diff(dataxrange),
              nudge_y = 2, hjust = 0, size = 3) +
    # add data
    geom_point(data=data, aes(TIME,DV, color = as.logical(CENS)), shape = 1, size = 3) +
    # # info on fitting
    # geom_text(
    #   data=subset(FIT, Phase == "KILL"),
    #   x=dataxrange[1],y=datayrange[2]+0.1*diff(datayrange),
    #   aes(label=paste0(' (',fitinfo,')')), hjust = 0, vjust=1
    # ) +
    # geom_text(
    #   data=subset(FIT, Phase == "RECRUD"),
    #   x=dataxrange[2],y=datayrange[2]+0.1*diff(datayrange),
    #   aes(label=paste0('(',fitinfo,')')), hjust = 1, vjust=1
    # ) +
  scale_color_manual("BLQ value", values = c("navyblue","firebrick")) +
    facet_wrap(~USUBJID, scales = "free") +
    labs(x=xlabel,
         y=ylabel,
         title=title,
         caption = paste0("Activity: ", MMVbase::get_ActivityPath(ActivityPath))) +
    theme(legend.position = "bottom",
          plot.caption    = element_text(hjust=0))

  return(gr)
}
#' simulate_indivPK
#' Simulation of PK for individual
#'
#' Concentrations are simulated with the given dosing, individual PK parameters
#' at 100 steps between minTIME and maxTIME and the time of MIC defined in the input data frame.
#' The simulation at the time of MIC is flagged in FLAGmic.
#'
#' @param x data frame with dosing records (defined by TIME, AMT and RATE), PK parameters and the time of MIC (Xknot) for one individual
#' @param model IQRmodel for the PK
#'
#' @return data frame with simulated PK timecourse (TIME and )
#'
#' @family Model-free PD assessment by linear segments
#' @author Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
simulate_indivPK <- function(x,
                             model = model)
{

  # Get model parameter names:
  parameterNames <- names(model$parameters)

  # Get regression parameters in dataset:
  regressionNames  <- intersect(names(x),parameterNames)
  regressionValues <- unlist(unique(x[,regressionNames]))
  names(regressionValues) <- regressionNames

  # Get dosing:
  dosing <- IQRdosing(
    TIME = x$TIME,
    ADM  = 1,
    AMT  = x$AMT,
    RATE = x$RATE
  )

  # Simulation time:
  tMIC         <- x$Xknot[1]
  simsteps     <- 100
  simtimerange <- c(x$minTIME[1],max(x$maxTIME[1],tMIC))
  simtimerange <- simtimerange * c(0.9,1.1)
  simtime      <- unique(sort(c(seq(simtimerange[1],simtimerange[2],length.out=simsteps),tMIC)))

  # Simulate:
  simres  <- sim_IQRmodel(
    model,
    simtime         = simtime,
    parameters      = regressionValues,
    dosingTable     = dosing,
    FLAGoutputsOnly = TRUE
  )

  # Indicate MIC:
  simres$FLAGmic <- (simres$TIME==tMIC)

  return(as.data.frame(simres))

}
