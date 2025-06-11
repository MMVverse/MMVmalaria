#' generate_InitialParametersIQRsys
#'
#' @description
#' @param projectPath Path of the project from which to retrieve modelSpec
#' @param IIVvalues0 Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family SYSFITmodeling
generate_InitialParametersIQRsys <- function(projectPath,
                                             IIVvalues0 = NULL) {

  # Get modelSpec of the model:
  estText   <- readRDS(file.path(projectPath,"project.est"))
  modelSpec <- estText$modelSpec

  # Get the estimated POP value:
  sysModel   <- load_IQRsysProject(projectPath, FLAGresultsOnly = TRUE)
  POPvalues0 <- getPars_IQRsysModel(sysModel)

  # Update POPvalues0 in modelspec:
  str_Paras  <- intersect(names(modelSpec$POPestimate), names(POPvalues0))
  modelSpec$POPvalues0[str_Paras] <- POPvalues0[str_Paras]

  # Manage IIVvalues0:
  if (!is.null(IIVvalues0)){
    if (is.numeric(IIVvalues0)){
      # Get parameters idx for which the IIV is gonna be estimated
      idx_Est <- which(as.numeric(modelSpec$IIVestimate)==1)

      # If IIVvalues has no names and length of one
      if (is.null(names(IIVvalues0)) && length(IIVvalues0)==1){
        modelSpec$IIVvalues0[idx_Est] <- IIVvalues0

      # If IIVvalues is a numeric vector
      }else{
        # Look for the parameters in POPvalues to change:
        idx_Name.POPvalues0 <- match(names(IIVvalues0), names(POPvalues0))

        # Look for the parameters in IIVvalues0 that can be changed and are estimated:
        idx_Name.IIVvalues0 <- match(idx_Est,idx_Name.POPvalues0)

        # Retain only the paramerter to be estimated:
        idx_Name.POPvalues0 <- idx_Name.POPvalues0[idx_Name.IIVvalues0]
        idx_Name.POPvalues0 <- idx_Name.POPvalues0[!is.na(idx_Name.POPvalues0)]
        idx_Name.IIVvalues0 <- idx_Name.IIVvalues0[!is.na(idx_Name.IIVvalues0)]

        # Update:
        modelSpec$IIVvalues0[idx_Name.POPvalues0] <- as.numeric(IIVvalues0[idx_Name.IIVvalues0])

      }

    # If the estimated values should be used as initial
    } else if(IIVvalues0=="Estimated"){
      # IIVvalues00          <- IQRnlmeResult$randomEffects$values
      # names(IIVvalues00)   <- names(POPvalues0)
      # modelSpec$IIVvalues0 <- IIVvalues00

      # TO BE IMPLEMENTED
    }
  }

  # Output:
  return(modelSpec)

}
#' generate_SYSalgorithmSettings
#'
#' @description
#' @param SIMOPT.method Default: 'lsodes'
#' @param SIMOPT.atol Default: 1e-06
#' @param SIMOPT.rtol Default: 1e-06
#' @param SIMOPT.hmin Default: 0
#' @param SIMOPT.hmax Default: NULL
#' @param SIMOPT.hini Default: 0
#' @param SIMOPT.maxsteps Default: 5000
#' @param SIMOPT.nauxtimes Default: 0
#' @param opt.nfits Default: 10
#' @param opt.sd Default: 1
#' @param opt.rinit Default: 1
#' @param opt.rmax Default: 10
#' @param opt.iterlim Default: 100
#' @param opt.prior_sigma Default: 10
#' @param opt.parlower Default: NULL
#' @param opt.parupper Default: NULL
#' @param algOpt.SEED Default: 123456
#' @param FLAGprofileLL Default: FALSE
#' @param FLAGkeepFits Default: FALSE
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family SYSFITmodeling
generate_SYSalgorithmSettings <- function(# Solver Settings:
                                          SIMOPT.method    = "lsodes",
                                          SIMOPT.atol	     = 1e-06,
                                          SIMOPT.rtol	     = 1e-06,
                                          SIMOPT.hmin	     = 0,
                                          SIMOPT.hmax	     = NULL,
                                          SIMOPT.hini	     = 0,
                                          SIMOPT.maxsteps	 = 5000,
                                          SIMOPT.nauxtimes = 0,

                                          # Algorithm Settings:
                                          opt.nfits        = 10,
                                          opt.sd           = 1,
                                          opt.rinit        = 1,
                                          opt.rmax         = 10,
                                          opt.iterlim      = 100,
                                          opt.prior_sigma  = 10,
                                          opt.parlower     = NULL,
                                          opt.parupper     = NULL,
                                          algOpt.SEED      = 123456,

                                          # Various Flags:
                                          FLAGprofileLL = FALSE,
                                          FLAGkeepFits  = FALSE) {

  # Create output:
  setting <- list()

  # Solver Settings:
  setting$SIMOPT.method    <- SIMOPT.method
  setting$SIMOPT.atol	     <- SIMOPT.atol
  setting$SIMOPT.rtol	     <- SIMOPT.rtol
  setting$SIMOPT.hmin	     <- SIMOPT.hmin
  setting$SIMOPT.hmax	     <- SIMOPT.hmax
  setting$SIMOPT.hini	     <- SIMOPT.hini
  setting$SIMOPT.maxsteps  <- SIMOPT.maxsteps
  setting$SIMOPT.nauxtimes <- SIMOPT.nauxtimes

  # Algorithm Settings:
  setting$opt.nfits       <- opt.nfits
  setting$opt.sd          <- opt.sd
  setting$opt.rinit       <- opt.rinit
  setting$opt.rmax        <- opt.rmax
  setting$opt.iterlim     <- opt.iterlim
  setting$opt.prior_sigma <- opt.prior_sigma
  setting$opt.parlower    <- opt.parlower
  setting$opt.parupper    <- opt.parupper
  setting$algOpt.SEED     <- algOpt.SEED

  # Various Flags:
  setting$FLAGprofileLL <- FLAGprofileLL
  setting$FLAGkeepFits  <- FLAGkeepFits

  # Return Dataset:
  return(setting)
}

#' table_BestFitEstimation
#'
#' @description
#' @param fit
#' @param filename Default: NULL
#' @param FLAGout Default: FALSE
#' @param title Default: NULL
#' @param footer Default: NULL
#' @return
#' @export
#' @importFrom MMVbase aux_formatErrorName
#' @author Aline Fuchs (MMV), Anne Kümmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family SYSFITmodeling
#' @importFrom plyr ldply
table_BestFitEstimation <- function(fit, filename = NULL, FLAGout = FALSE, title = NULL, footer = NULL) {
  # Input is a list of fitting objects produced by run_SysFitEstimation
  # For the IQRsysFit the field BestFit is used containing best fits the best fit based on OFV.
  # If interaction parameters are estimated (AlphaXX and Beta XX) these are converted to
  # interpretable values.

  # concatenate table with best models for all models:
  table <- plyr::ldply(fit, function(x) {

    FLAGci <- ("xopt_95CI_low" %in% names(x$BestFit))

    # Get estimates for hill, and estimated parameters and their estimates
    idxEST    <- x$BestFit$estObject$parameters$estimate==1
    idxHILL   <- tolower(names(x$BestFit$estObject$parameters$estimate))=="hill"
    estimates <- x$BestFit$xopt[idxHILL | idxEST]
    if (FLAGci) {
      CIlow  <- x$BestFit$xopt_95CI_low[idxHILL | idxEST]
      CIhigh <- x$BestFit$xopt_95CI_high[idxHILL | idxEST]
    } else {
      RSE <- x$BestFit$xopt_rse[idxHILL | idxEST]
    }
    OFV <- x$BestFit$fopt

    # # Translate interaction parameters
    # # - alpha: convert to maximum fold change of EC50 value
    # idxALPHA            <- grep("Alpha", names(estimates))
    # estimates[idxALPHA] <- exp(-estimates[idxALPHA])
    # if (FLAGci) {CIlow[idxALPHA] <- exp(-CIhigh[idxALPHA]); CIhigh[idxALPHA] <- exp(-CIlow[idxALPHA])} # note: swap of boundaries is correct
    # # - beta: convert to maximum percent change of Emax value
    # idxBETA            <- grep("Beta", names(estimates))
    # estimates[idxBETA] <- (exp(estimates[idxBETA]) - 1) *100
    # if (FLAGci) {CIlow[idxBETA] <- (exp(CIlow[idxBETA]) - 1) *100; CIhigh[idxBETA] <- (exp(CIhigh[idxBETA]) - 1) *100}


    # Translate interaction parameters
    # - alpha: convert to maximum fold change of EC50 value
    idxALPHA            <- grep("Alpha", names(estimates))
    estimates[idxALPHA] <- exp(-estimates[idxALPHA])
    if (FLAGci) {
      # Note: Boundaries might be in different order, due to +/- signs.
      #       It will be corrected later.
      CIlow[idxALPHA]  <- exp(-CIlow[idxALPHA])
      CIhigh[idxALPHA] <- exp(-CIhigh[idxALPHA])
    }
    # - beta: convert to maximum percent change of Emax value
    idxBETA            <- grep("Beta", names(estimates))
    estimates[idxBETA] <- (exp(estimates[idxBETA]) - 1) *100
    if (FLAGci) {
      CIlow[idxBETA]  <- (exp(CIlow[idxBETA]) - 1) *100
      CIhigh[idxBETA] <- (exp(CIhigh[idxBETA]) - 1) *100
    }

    # Adjust order of boundaries:
    if (FLAGci) {
      CIlow_temp  <- CIlow
      CIhigh_temp <- CIhigh
      CIlow       <- ifelse(CIlow_temp  < CIhigh_temp, CIlow_temp , CIhigh_temp)
      CIhigh      <- ifelse(CIhigh_temp > CIlow_temp , CIhigh_temp, CIlow_temp )
    }

    # Change name of Alpha and Beta:
    idxERROR                   <- grep("error", names(estimates))
    names(estimates)[idxALPHA] <- paste0(names(estimates)[idxALPHA]," (F.Ch.)")
    names(estimates)[idxBETA]  <- paste0(names(estimates)[idxBETA]," (% Ch.)")

    if (FLAGci) {
      y <- data.frame(
        row.names = c("OFV",names(estimates)),
        Estimate  = c(round(OFV), signif(estimates,3))
      )
      y$Estimate <- ifelse(c(TRUE, is.na(CIlow)), y$Estimate, paste0(y$Estimate, " [",sprintf("%.2g", c(NA,CIlow)),";",sprintf("%.2g", c(NA,CIhigh)),"]"))
    } else {
      y <- data.frame(
        row.names = c("OFV",names(estimates)),
        Estimate  = c(round(OFV), signif(estimates,3)),
        RSE       = c(NA,uncertain)
      )
      y <- within(y, {
        Estimate <- ifelse(is.na(RSE),Estimate,paste0(Estimate, " (",RSE, ")"))
        RSE      <- NULL
      })
    }

    y <- as.data.frame(t(y))

    # Nicer error parameter names
    names(y) <- MMVbase::aux_formatErrorName(names(y))

    y
  }, .id = "Model")

  # Put Error Column at the end:
  idxERROR <- grep("error", names(table))
  table    <- table[, c((1:ncol(table))[-idxERROR], idxERROR)]

  # Order by objective function value
  table$OFV <- as.numeric(as.character(table$OFV))
  table     <- table[order(table$OFV),]

  if (!is.null(filename))
    IQRoutputTable(table, filename = filename, xtitle = title, xfooter = footer, report = TRUE)

  if (FLAGout) table
}
