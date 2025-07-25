#' Estimate MIC, MPC90, PRR48 parameters
#'
#' @description Required: `Emax` `EC50`, `hill`. Emax, cleareance, delay, turnover and effect compartment models
#' @param modelFolder Relative or absolute path to the IQRnlmeProject PD model
#' @param dataPath Path to the dataset
#' @param outputFolder Folder where the plots of MIC, MPC90 and PRR48 will be saved
#' @param yaxis Default: `NULL`
#' @param yaxisunit Default: ''
#' @param GRmodelFolder Default: `NULL`
#' @param PKmodelFolder Default: '../04-Output/S131-01_PKmodel'
#' @param Apparent Default: `TRUE`
#' @param Tk0 Default: 'Tk0'
#' @param Indv Default: `FALSE`
#' @param ActivityPath Default: `NULL`
#' @param addArgs_sim_IQRmodel named list (default NULL) with additional arguments to be passed to sim_IQR_model. This
#' argument is relevant only if `Apparent == TRUE`.
#' @return Data frame or list of calculated parameters  
#' @export
#' @author Mohammed H. Cherkaoui (MMV), Sam Jones (MMV)
#' @family Key Parameters
assess_MICMPC90PRR48 <- function(modelFolder,
                                 dataPath,
                                 outputFolder,
                                 yaxis         = NULL,
                                 yaxisunit     = "",
                                 GRmodelFolder = NULL,
                                 PKmodelFolder = "../04-Output/S131-01_PKmodel",
                                 Apparent      = TRUE,
                                 Tk0           = "Tk0",
                                 Indv          = FALSE,
                                 ActivityPath  = NULL,
                                 addArgs_sim_IQRmodel = NULL) {

  #-------------------------------------#
  # Initial Settings ----
  #-------------------------------------#

  # Get individual and population parameters:
  IndParPD <- getIndivParameters_IQRnlmeProject(modelFolder)
  PopParPD <- getPopParameters_IQRnlmeProject(modelFolder)
  PopParPD <- as.data.frame(t(PopParPD))

  # Check which standard model:
  if ("CLPara" %in% names(PopParPD)) {FLAGclearance <- TRUE; cat("Assume Clearance model.\n")         } else FLAGclearance <- FALSE
  if ("ke" %in% names(PopParPD))     {FLAGeffectCpt <- TRUE; cat("Assume Effect Cpt model.\n")        } else FLAGeffectCpt <- FALSE
  if ("kin" %in% names(PopParPD))    {FLAGturnover  <- TRUE; cat("Assume Turnover model.\n")          } else FLAGturnover  <- FALSE

  # Check if Turnover Clearance model:
  if (FLAGturnover & FLAGclearance){
    FLAGturnclear <- TRUE
    FLAGturnover  <- FALSE
    FLAGclearance <- FALSE
    cat("Assume Turnover Clearance model.\n")
  } else {
    FLAGturnclear <- FALSE
  }

  if (all(!c(FLAGturnover,FLAGeffectCpt,FLAGclearance,FLAGturnclear)))
    cat("Assume Emax model.\n")

  if (sum(c(FLAGturnover,FLAGeffectCpt,FLAGclearance,FLAGturnclear)) > 1)
    stop("Standard model not detected unambigously.\n")


  # If growth model folder given get growth information from this model:
  if (!is.null(GRmodelFolder)) {
    IndParGR  <- getIndivParameters_IQRnlmeProject(GRmodelFolder)
    PopParGR  <- getPopParameters_IQRnlmeProject(GRmodelFolder)
    PopParGR  <- as.data.frame(t(PopParGR))
    IndParPD  <- merge(IndParPD, IndParGR[,c("ID","GR")])
    PopParPD  <- cbind(PopParPD, PopParGR)
  }

  # Get the dataset for subject information:
  data     <- IQRloadCSVdata(dataPath)
  dataAttr <- loadATRinfo_csvData(dataPath)
  subjInfo <- unique(data[data$TAD >= 0, c("STUDY","USUBJID", "ID", "TRTNAME","DOSE", dataAttr$regressorNames)])
  doseInfo <-unique(data[data$TYPENAME == "Dose", c("USUBJID", "ID", "TIME", "AMT", "ADM","TRTNAME")])
  IndPar   <- merge(IndParPD, subjInfo)


  #-------------------------------------#
  # Calculate MIC, MPC90 and PPR48 ----
  #-------------------------------------#

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Do everything in case of effect compartment and turnover model by simulations
  if ((FLAGeffectCpt|FLAGturnover|FLAGturnclear) & Apparent ) {

    # Define sim time:
    simtime <- seq(0,5000,0.5)

    # Load Model:
    modelResults <- getResults_IQRnlmeProjectMulti(as_IQRnlmeProjectMulti(modelFolder))[[1]]
    model        <- IQRmodel(file.path(modelResults$model, "model.txt"))

    # PK parameters:
    trtGroups <- unique(subjInfo[,c("TRTNAME","DOSE", "STUDY")])

    if (dim(trtGroups)[1] == 1) {
      PopParPK <- getPopParameters_IQRnlmeProject(PKmodelFolder, IndCovariates = trtGroups[c(1,1),])
      PopParPK <- PopParPK[1,]
    } else {
      PopParPK <- getPopParameters_IQRnlmeProject(PKmodelFolder, IndCovariates = trtGroups)
    }

    # Define ID:
    trtGroups$ID <- (1:dim(trtGroups)[1]) + 10000

    # Combined PK and PD parameters:
    PopPar <- within(cbind(trtGroups,PopParPK, PopParPD), {
      ID <- (1:dim(trtGroups)[1]) + 10000
      PLbase <- mean(IndPar$PLbase)
      PLerr  <- 0
      USUBJID <- "POP"
    })
    # PopPar <- PopPar[!duplicated(PopPar$DOSE), ]
    PopPar <- PopPar[!duplicated(PopPar$TRTNAME), ]

    # Event table for individual and population simulation:
    colNames <- intersect(names(PopPar), names(IndPar))
    eventTable <- rbind(PopPar[,colNames], IndPar[,colNames])
    eventTable <- dplyr::left_join(unique(doseInfo[c("TIME","AMT","ADM","TRTNAME")]), eventTable, by = "TRTNAME")
    eventTable <- eventTable[order(eventTable$ID, eventTable$TIME, eventTable$ADM, eventTable$AMT),]

    if (modelResults$projectHeader$DOSINGTYPES == "ABSORPTION0"){
      eventTable$TINF <- eventTable[,Tk0]
    }
    eventTable <- IQReventTable(eventTable, regression = unique(c(names(PopParPD), names(PopParPK), "PLbase")))

    # Run Simulation:
    simRes <- do.call(
      sim_IQRmodel,
      c(
        list(model = model,
             simtime = simtime,
             eventTable = eventTable),
        addArgs_sim_IQRmodel))

    # Adjust PL column for Turnover Clearance model:
    if (FLAGturnclear){
      simRes$PL <- simRes$PLliv
    }

    # Get the "Key" parameters for each individual:
    AppKeys <- plyr::ddply(simRes, ~ID, getApparentKeys)

    # Split Indiv and Pop values:
    IndPar  <- plyr::join(AppKeys[AppKeys$ID < 10000,], subjInfo[,c("STUDY","USUBJID", "ID", "TRTNAME", "DOSE")])
    PopPar  <- plyr::join(AppKeys[AppKeys$ID > 10000,], trtGroups)

    # Prepare plotting (MIC and MPC to be plotted in ng/mL, therefore multiply with 1000):
    IndParPlot <- reshape(
      within(IndPar[, c("STUDY","USUBJID", "ID", "TRTNAME", "DOSE", "MIC", "MPC90", "PRR48")], {MIC <- MIC*1000; MPC90 <- MPC90*1000}),
      direction = "long",
      varying   = c("MIC", "MPC90", "PRR48"),
      v.names   = "Value",
      times     = c("MIC App (ng/mL)", "MPC90 App (ng/mL)", "PRR48 App (log10)"),
      timevar   = "Parameter"
    )
    PopParPlot <- reshape(
      within(PopPar[, c("MIC", "MPC90", "PRR48", "DOSE")], {MIC <- MIC*1000; MPC90 <- MPC90*1000}),
      direction = "long",
      varying   = c("MIC", "MPC90", "PRR48"),
      v.names   = "Value",
      times     = c("MIC App (ng/mL)", "MPC90 App (ng/mL)", "PRR48 App (log10)"),
      timevar   = "Parameter"
    )

    if (!is.null(yaxis))
      warning("Y-axis of plots for turnover and clearance model always DOSE. 'yaxis' input argument ignored.")

    # Plot:
    gr <- IQRggplot(IndParPlot, aes(x=DOSE, y=Value)) +
      geom_point(position=position_jitter(width=0.2))  +
      geom_point(data = PopParPlot, color = "firebrick", shape = 3) +
      geom_text(data = PopParPlot, aes(label=round(Value,2)),
                color = "firebrick", vjust = 0, hjust=0) +
      facet_wrap(~Parameter, scales = "free_y") +
      scale_x_continuous("DOSE (mg)") +
      labs(x       = "",
           caption = paste0("Activity: ", MMVbase::get_ActivityPath(ActivityPath),
                            "\nModel: ", modelFolder)) +
      theme(strip.text   = element_text(size=12, face="bold"),
            plot.caption = element_text(hjust=0))

    # only return "PRR48","MPC90","MIC", but no other parameter values:
    #   POP values
    PopParRelevant <- data.frame(Model = aux_fileparts(modelFolder)$filename,
                                 stringsAsFactors = FALSE)
    PopParRelevant <- cbind(PopParRelevant, PopPar[c("DOSE","PRR48","MPC90","MIC")])
    #   INDIV values
    IndParRelevent <- data.frame(Model = aux_fileparts(modelFolder)$filename,
                                 stringsAsFactors = FALSE)
    IndParRelevent <- cbind(IndParRelevent, IndPar[, c("ID","USUBJID","DOSE","PRR48","MPC90","MIC")])


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Do everything in case of clearance and Emax model by formulas
  } else {

    if (FLAGclearance|FLAGturnclear) {
      IndPar <- getKeysEMAX(IndPar)
      PopPar <- getKeysEMAX(PopParPD)
    } else {
      IndPar <- getKeysEMAX(IndPar)
      PopPar <- getKeysEMAX(PopParPD)
    }
    # Prepare plotting (MIC and MPC to be plotted in ng/mL, therefore multiply with 1000)
    IndParPlot <- reshape(
      within(IndPar[, c("STUDY","USUBJID", "ID", "TRTNAME", "DOSE", "MIC", "MPC90", "PRR48")],{MIC <- MIC*1000; MPC90 <- MPC90*1000}),
      direction = "long",
      varying   = c("MIC", "MPC90", "PRR48"),
      v.names   = "Value",
      times     = c("MIC Static (ng/mL)", "MPC90 Static (ng/mL)", "PRR48 Static (log10)"),
      timevar   = "Parameter"
    )
    PopParPlot <- reshape(
      within(PopPar[, c("MIC", "MPC90", "PRR48")], {MIC <- MIC*1000; MPC90 <- MPC90*1000}),
      direction = "long",
      varying   = c("MIC", "MPC90", "PRR48"),
      v.names   = "Value",
      times     = c("MIC Static (ng/mL)", "MPC90 Static (ng/mL)", "PRR48 Static (log10)"),
      timevar   = "Parameter"
    )

    # Plot:
    if (is.null(yaxis)) {
      gr <- IQRggplot(IndParPlot, aes(x=0, y=Value)) +
        geom_point(position=position_jitter(width=0.2))  +
        geom_hline(data = PopParPlot, aes(yintercept=Value), color = "firebrick") +
        geom_text(data = PopParPlot, aes(x=-0.4,y=Value, label=paste0("Pop. estimate: ",round(Value,2))),
                  color = "firebrick", vjust = -0.2, hjust=0) +
        facet_wrap(~Parameter, scales = "free_y") +
        coord_cartesian(xlim = c(-0.5,0.5)) +
        scale_x_continuous(breaks = -1) +
        labs(x       = "",
             caption = paste0("Activity: ", MMVbase::get_ActivityPath(ActivityPath),
                              "\nModel: ", modelFolder)) +
        theme(strip.text   = element_text(size=12, face="bold"),
              plot.caption = element_text(hjust=0))

    } else {
      yaxisMin <- min(IndParPlot[[yaxis]])
      yaxisMax <- max(IndParPlot[[yaxis]])
      gr <- IQRggplot(IndParPlot, aes_string(x=yaxis, y="Value")) +
        geom_point(position=position_jitter(width=0.2))  +
        geom_hline(data = PopParPlot, aes(yintercept=Value), color = "firebrick") +
        geom_text(data = PopParPlot, aes(y=Value, label=paste0("Pop. estimate: ",round(Value,2))),
                  x=yaxisMin-0.4,
                  color = "firebrick", vjust = -0.2, hjust=0) +
        facet_wrap(~Parameter, scales = "free_y") +
        coord_cartesian(xlim = c(yaxisMin-1,yaxisMax+1)) +
        labs(x       = paste0(yaxis," (", yaxisunit,")"),
             y       = "",
             caption = paste0("Activity: ", MMVbase::get_ActivityPath(ActivityPath),
                              "\nModel: ", modelFolder)) +
        theme(strip.text   = element_text(size=12, face="bold"),
              plot.caption = element_text(hjust=0))
    }

    # only return "PRR48","MPC90","MIC", but no other parameter values:
    #   POP values
    PopParRelevant <- data.frame(Model = aux_fileparts(modelFolder)$filename,
                                 stringsAsFactors = FALSE)
    PopParRelevant <- cbind(PopParRelevant, PopPar[c("PRR48","MPC90","MIC")])
    #   INDIV values
    IndParRelevent <- data.frame(Model = aux_fileparts(modelFolder)$filename,
                                 stringsAsFactors = FALSE)
    IndParRelevent <- cbind(IndParRelevent, IndPar[, c("ID","USUBJID","PRR48","MPC90","MIC")])
  }



  #-------------------------------------#
  # Outputs ----
  #-------------------------------------#

  # Plot to file:
  if ((FLAGeffectCpt|FLAGturnover|FLAGturnclear) & Apparent ){
    IQRoutputPNG(gr, file.path(outputFolder, "MIC_MPC90_PRR48_Apparent"))
  } else{
    IQRoutputPNG(gr, file.path(outputFolder, "MIC_MPC90_PRR48_Static"))
  }

  # Check if individual parameters need to be returned:
  if (Indv){
    out <- list(PopParRelevant  = PopParRelevant,
                IndParRelevent  = IndParRelevent)
  } else{
    out <- PopParRelevant
  }

  # Attributes:
  ModelName <- "Emax"
  if (FLAGturnover) ModelName <- "Turnover"
  if (FLAGclearance) ModelName <- "Clearance"
  if (FLAGeffectCpt) ModelName <- "EffectCpt"
  if (FLAGturnclear) ModelName <- "TurnoverClearance"
  attributes(out)$ModelName <- ModelName

  # Output:
  return(out)
}

#' convert_SCIDconcentration
#'
#' @description convert_SCIDconcentration to selected matrix
#' @param ConcToConv Numeric vector containing concentration data to convert
#' @param HuRbp Numeric: value of human ratio blood:plasma
#' @param HuPPB Numeric: % value of human plasma protein binding as a proportion
#' @param MouseRbp Numeric: value of mice ratio blood:plasma
#' @param MousePPB Numeric: % value of mice plasma protein binding as a proportion
#' @param HuH Numeric: Value of % human hematocrit as a proportion. Default: 0.45 
#' @param MouseH Numeric: Value of % mouse hematocrit as a proportion. Default: 0.45
#' @param huRBCinSCID Numeric: value of % human red blood cells in SCID as proportion. Default: 0.6
#' @param scidH Numeric: Value of % SCID  hematocrit as a proportion. Default: 0.8
#' @param ConvertTo Character string defining matrix to convert to. One of 'Human Blood', 'Human Plasma', 'Mouse Blood', 'Mouse Plasma' or 'Free'. Default: Human Blood
#' @return Numeric vector containing converted concentration data 
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
convert_SCIDconcentration <- function(ConcToConv, HuRbp, HuPPB, MouseRbp, MousePPB,
                                      HuH = 0.45, MouseH = 0.45, huRBCinSCID = 0.6, scidH = 0.8,
                                      ConvertTo = c("Human Blood", "Human Plasma", "Mouse Blood", "Mouse Plasma", "Free")[1])
{
  if (!ConvertTo %in%  c("Human Blood", "Human Plasma", "Mouse Blood", "Mouse Plasma", "Free")) stop("Parameter 'ConvertTo' should be equal to 'Human Blood', 'Human Plasma', 'Mouse Blood', 'Mouse Plasma' or 'Free'.")

  # Estimate SCID Rbp:
  scidRbp <- formula_SCIDrbp(HuRbp       = HuRbp,
                             HuPPB       = HuPPB,
                             HuH         = HuH,
                             MouseRbp    = MouseRbp,
                             MousePPB    = MousePPB,
                             MouseH      = MouseH,
                             scidH       = scidH,
                             huRBCinSCID = huRBCinSCID)

  # Determine conversion target
  matrix_to <- dplyr::case_when(grepl("Blood", ConvertTo) ~ "blood",
                                grepl("Plasma", ConvertTo) ~ "plasma",
                                TRUE ~ "unbound")

  Rbp_to   = 1
  PPB_to   = 0
  if (grepl("Human", ConvertTo)) {
    Rbp_to   = HuRbp
    PPB_to   = HuPPB
  }
  if (grepl("Mouse", ConvertTo)) {
    Rbp_to   = MouseRbp
    PPB_to   = MousePPB
  }

  # Convert concentration
  ConvertedConc <- MMVbase::formula_convertConc(ConcToConv,
                                       from = "blood",
                                       to   = matrix_to,
                                       Rbp_from = scidRbp,
                                       PPB_from = MousePPB,
                                       Rbp_to   = Rbp_to,
                                       PPB_to   = PPB_to)


  # Output:
  return(ConvertedConc)
}

#' formula_SCIDrbp
#' @description Calculate ratio blood plasma in SCID
#'
#' @param HuRbp Human ratio blood plasma
#' @param HuPPB Human plasma protein binding
#' @param HuH Human hematocrit
#' @param MouseRbp Mouse ratio blood plasma
#' @param MousePPB Mouse plasma protein binding
#' @param MouseH Mouse hematocrit
#' @param scidH SCID hematocrit
#' @param huRBCinSCID Human red blood cells in SCID
#'
#' @return numeric, SCID ratio blood plasma
#' @author Anne Kuemmel (IntiQuan)
#' @family Key Parameters
#' @export
formula_SCIDrbp <- function(HuRbp,
                            HuPPB,
                            HuH = 0.45,
                            MouseRbp,
                            MousePPB,
                            MouseH = 0.45,
                            scidH  = 0.8,
                            huRBCinSCID = 0.6) {

  # Estimate the Kp,RBC:
  HuKpRBC    <- (HuRbp    - 1 + HuH   )/ HuH
  MouseKpRBC <- (MouseRbp - 1 + MouseH)/ MouseH

  # Take 0 if negative:
  HuKpRBC    <- max(HuKpRBC   ,0)
  MouseKpRBC <- max(MouseKpRBC,0)

  # Estimate Fup:
  HuFup    <- 1 - HuPPB
  MouseFup <- 1 - MousePPB

  # Estimate SCID Rbp:
  scidRbp <- 1 - scidH + MouseKpRBC*scidH*(1-huRBCinSCID)/(1+huRBCinSCID) +
    HuKpRBC*MouseFup/HuFup*scidH*(2*huRBCinSCID)/(1+huRBCinSCID)


  return(scidRbp)
}

#' estimate_AUC
#'
#' @description Estimate area under the curve (AUC) of provided time v concentration data
#' @param Time Numeric vector containing time data 
#' @param Conc Numeric vector containing concentration data
#' @param AUCtype Default: 'last'
#' @param intMethod Default: MMVbase::logLinTrapzMMV
#' @param LambdaZ Default: function(x, y) {
#'    estimate_LambdaZ(x, y)$LambdaZ
#'}
#' @return Estimated AUC for the provided `AUCtype`
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @importFrom MMVbase logLinTrapzMMV
#' @family Key Parameters
estimate_AUC <- function(Time,
                         Conc,
                         AUCtype   = "last",
                         intMethod = MMVbase::logLinTrapzMMV,
                         LambdaZ   = function(x,y){estimate_LambdaZ(x,y)$LambdaZ},
                         atol      = 1e-9){

  #--------------------------------------#
  # STEP 0: Checks ----
  #--------------------------------------#

  # Check if some concentration lower than -ATOL (ATOL should be defined within the solver tolerance)
  #   - Should not be the case
  #   => Generated Error
  if (any(Conc < (-atol))){
    stop("Negative concentration (i.e. 'Conc'<-", atol, "): Either the model used for simulation is not well coded, or the tolerance set on the solver are not matching the one is 'estimate_AUC'")
  }

  #--------------------------------------#
  # STEP 1: Index to keep for AUC ----
  #--------------------------------------#

  # Initialize idx_AUC:
  idx_AUC <- 1:length(Time)

  # Scenario I:
  #   - If AUCinf
  #   => Need to check to better estimate Lambda
  if (tolower(AUCtype)=="inf"){

    # Scenario I.1:
    #   - Some concentrations are between -ATOL and ATOL
    #   - Expected behavior if simulation time is long enough
    if(any(abs(Conc)<=atol)){

      # Number of concentrations above ATOL:
      n_Pos <- length(which(Conc>atol))

      # Scenario I.1.A:
      #   - At least two concentrations larger than ATOL
      if(n_Pos > 1){

        # Last position of concentration larger than ATOL:
        Last_Pos <- max(which(abs(Conc)>atol))

        # Correct idx_AUC:
        idx_AUC <- 1:Last_Pos


        # Scenario I.1.B:
        #   - All concentrations smaller than ATOL (in absolute value)
        #   - Could also only have one concentration larger than ATOL
        #   => Change AUCtype to last
      }else{
        AUCtype <- "last"
        warning("All concentrations are lower than 'atol'=", atol, "which would make the estimation of LambdaZ difficult: Therefore AUClast is returned")
      }


      # Scenario I.2:
      #   - All concentrations are larger than ATOL
      #   => Take all index
      #   => No need of "else if" or "else" statement
    }

    # Scenario II:
    #   - If AUClast
    #   => Take all index
    #   => No need of "else if" or "else" statement
  }

  # Subset of Time and Conc:
  Time <- Time[idx_AUC]
  Conc <- Conc[idx_AUC]


  #--------------------------------------#
  # STEP 2: Estimate AUC last ----
  #--------------------------------------#
  AUC <- intMethod(x = Time,y = Conc)


  #--------------------------------------#
  # STEP 3: Estimate AUCinf if needed ----
  #--------------------------------------#
  if (tolower(AUCtype)=="inf" & Conc[length(Conc)]!=0){
    # Estimate Lambda_Z:
    if (is.function(LambdaZ)){
      Lambda_Z <- tryCatch({
        out <- LambdaZ(Time,Conc)
        out
      }, error = function(e){
        out <- 1000
        out
      })
      # Lambda_Z <- LambdaZ(Time,Conc)
    }else if (is.numeric(LambdaZ)){
      Lambda_Z <- LambdaZ
    } else{
      stop("'LambdaZ' is neither a function or a scalar: Please update.")
    }

    # AUCinf:
    AUC <- AUC + Conc[length(Conc)]/Lambda_Z
  }

  # Return Results:
  return(AUC)
}


#' estimate_LambdaZ
#'
#' @description Estimates the slope (lambda-z) of provided time v concentration data
#' @param Time Numeric vector containing time data 
#' @param Conc Numeric vector containing concentration data 
#' @param R2.min Default: 0.95
#' @param atol Default: 1e-12
#' @return Estimated lambda-z
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
estimate_LambdaZ <- function(Time,
                             Conc,
                             R2.min = 0.95,
                             atol   = 1e-12){

  # Check that all Concentration values are positive:
  #   If not remove
  idx_Pos <- which(Conc>=0)
  if (length(idx_Pos)!=length(Conc)){
    Time <- Time[idx_Pos]
    Conc <- Conc[idx_Pos]
    warning("Some concentrations were negative and removed.")
  }

  # length of vectors:
  n <- length(Time)

  # Number of points to start the regression:
  k <- min(max(floor(n/10),4),n)

  # Loop to estimate the elimination rate:
  R2 <- 0
  LambdaZ <- 1e6
  while (k>2 & R2<R2.min){
    # Data to be used for the linear regression:
    dataLM <- data.frame(Time     = Time[(n-k):n],
                         Conc_log = log(Conc[(n-k):n]))

    # Linear Regression:
    fitLM <- lm(Conc_log ~ Time, data = dataLM)

    # Extract results:
    LambdaZ <- -fitLM$coefficients[[2]]
    if (sum(fitLM$residuals)<=atol){
      R2 <- 1
    }else{
      R2 <- summary(fitLM)$r.squared
    }

    # Reduce k:
    k <- k - 1
  }

  # Wargnings in the case tha tonly three points were used
  if (k==2 & R2<R2.min){
    warning(paste0("R-Squared is samller than ", R2.min, " and only the last three points were used for the regression."))
  }

  # Create Output:
  results_LambdaZ <- list(LambdaZ = LambdaZ,
                          R2      = R2)

  # Return Results:
  return(results_LambdaZ)
}


#' formulaPKpar_MacroToMicro
#'
#' @description Converts provided macro PK parameters to micro parameters
#' @param PKparameters_Macro Named vector containing macro PK parameters and their values
#' @param PKmacroName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    Alpha = "Alpha", Beta = "Beta", Gamma = "Gamma", A = "A",
#'    B = "B", C = "C")
#' @param PKmicroName Default: c(K = "K", Vc = "Vc", Kc1 = "Kc1", K1c = "K1c", Kc2 = "Kc2",
#'    K2c = "K2c")
#' @return Named vector containing micro PK parameters and their values
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_MacroToMicro <- function(PKparameters_Macro,
                                      PKmacroName = c("Fabs0" = "Fabs0",
                                                      "Tk0"   = "Tk0"  ,
                                                      "Fabs1" = "Fabs1",
                                                      "ka"    = "ka",
                                                      "Tlag1" = "Tlag1",
                                                      "Alpha" = "Alpha",
                                                      "Beta"  = "Beta",
                                                      "Gamma" = "Gamma",
                                                      "A"     = "A",
                                                      "B"     = "B",
                                                      "C"     = "C"),
                                      PKmicroName = c("K"     = "K",
                                                      "Vc"    = "Vc",
                                                      "Kc1"   = "Kc1",
                                                      "K1c"   = "K1c",
                                                      "Kc2"   = "Kc2",
                                                      "K2c"   = "K2c")) {

  # Adjust Name in PKparameters_Macro:
  #   NOTE: Just to simplify the code afterward
  idx_Ordered <- intersect(PKmacroName, names(PKparameters_Macro))
  PKparameters_Macro <- PKparameters_Macro[idx_Ordered]
  idx_1 <- which(names(PKparameters_Macro) %in% PKmacroName)
  idx_2 <- which(PKmacroName %in% names(PKparameters_Macro))
  names(PKparameters_Macro)[idx_1] <- names(PKmacroName)[idx_2]

  # Generate output:
  PKparameters_Micro <- c()

  # PK Parameters to keep:
  if ("Fabs0" %in% names(PKparameters_Macro)){
    PKparameters_Micro[["Fabs0"]] <- PKparameters_Macro[["Fabs0"]]
  }
  if ("Tk0" %in% names(PKparameters_Macro)){
    PKparameters_Micro[["Tk0"]]   <- PKparameters_Macro[["Tk0"]]
  }
  if ("Fabs1" %in% names(PKparameters_Macro)){
    PKparameters_Micro[["Fabs1"]] <- PKparameters_Macro[["Fabs1"]]
  }

  if ("ka" %in% names(PKparameters_Macro)){
    PKparameters_Micro[["ka"]]    <- PKparameters_Macro[["ka"]]
  }

  if ("Tlag1" %in% names(PKparameters_Macro)){
    PKparameters_Micro[["Tlag1"]] <- PKparameters_Macro[["Tlag1"]]
  }

  # Determine the number of compartments:
  n_CPT <- get_NumberOfPKcpt(PKparameters_Macro, PKparaType = "Macro")

  # Estimate the macroal parameters:
  if (n_CPT==1){
    PKparameters_Micro[["K"]]  <- PKparameters_Macro[["Alpha"]]
    PKparameters_Micro[["Vc"]] <- 1/PKparameters_Macro[["A"]]

  }else if(n_CPT==2){
    PKparameters_Micro[["Vc"]]  <- 1/(PKparameters_Macro[["A"]] + PKparameters_Macro[["B"]])
    PKparameters_Micro[["K1c"]] <- (PKparameters_Macro[["A"]]*PKparameters_Macro[["Beta"]] + PKparameters_Macro[["B"]]*PKparameters_Macro[["Alpha"]])*PKparameters_Micro[["Vc"]]
    PKparameters_Micro[["K"]]   <- PKparameters_Macro[["Alpha"]]*PKparameters_Macro[["Beta"]]/PKparameters_Micro[["K1c"]]
    PKparameters_Micro[["Kc1"]] <- PKparameters_Macro[["Alpha"]] + PKparameters_Macro[["Beta"]] - PKparameters_Micro[["K1c"]] - PKparameters_Micro[["K"]]

  }else if(n_CPT==3){
    PKparameters_Micro[["Vc"]] <- 1/(PKparameters_Macro[["A"]] + PKparameters_Macro[["B"]] + PKparameters_Macro[["C"]])

    # Estimate K1c & K2c:
    Sum   <- PKparameters_Macro[["Alpha"]] + PKparameters_Macro[["Beta"]] - (PKparameters_Macro[["A"]]*(PKparameters_Macro[["Alpha"]] - PKparameters_Macro[["Gamma"]]) + PKparameters_Macro[["B"]]*(PKparameters_Macro[["Beta"]] - PKparameters_Macro[["Gamma"]]))*PKparameters_Micro[["Vc"]]
    Prod  <- PKparameters_Macro[["Alpha"]]*Sum - PKparameters_Macro[["Alpha"]]^2 + PKparameters_Macro[["A"]]*(PKparameters_Macro[["Alpha"]] - PKparameters_Macro[["Beta"]])*(PKparameters_Macro[["Alpha"]] - PKparameters_Macro[["Gamma"]])*PKparameters_Micro[["Vc"]]
    delta <- Sum^2-4*Prod
    PKparameters_Micro[["K1c"]] <- (Sum + sqrt(delta))/2
    PKparameters_Micro[["K2c"]] <- (Sum - sqrt(delta))/2

    # Get Intermediate Constant:
    a2  <- PKparameters_Macro[["Alpha"]] + PKparameters_Macro[["Beta"]] + PKparameters_Macro[["Gamma"]]
    phi <- atan(sqrt(3)*(PKparameters_Macro[["Beta"]] - PKparameters_Macro[["Gamma"]])/(PKparameters_Macro[["Beta"]] + PKparameters_Macro[["Gamma"]] - 2*PKparameters_Macro[["Alpha"]]))
    r2  <- (PKparameters_Macro[["Beta"]] + PKparameters_Macro[["Gamma"]] - 2*PKparameters_Macro[["Alpha"]])/(3*cos(phi))
    r1  <- r2^3/8
    p   <- -(27*r1^2)^(1/3)
    a1  <- p + a2^2/3
    q   <- -cos(3*phi)*2*r1
    a0  <- q - 2*a2^3/27 + a1*a2/3

    # Get remaining Micro Constant:
    PKparameters_Micro[["K"]] <- a0/(PKparameters_Micro[["K1c"]]*PKparameters_Micro[["K2c"]])
    K1 <- a2 - PKparameters_Micro[["K"]]  - PKparameters_Micro[["K1c"]] - PKparameters_Micro[["K2c"]]
    K2 <- a1 - PKparameters_Micro[["K"]]*PKparameters_Micro[["K1c"]] - PKparameters_Micro[["K"]]*PKparameters_Micro[["K2c"]] - PKparameters_Micro[["K1c"]]*PKparameters_Micro[["K2c"]]
    PKparameters_Micro[["Kc1"]] <- (K2 - K1*PKparameters_Micro[["K1c"]])/(PKparameters_Micro[["K2c"]] - PKparameters_Micro[["K1c"]])
    PKparameters_Micro[["Kc2"]] <- K1 - PKparameters_Micro[["Kc1"]]
  }

  # For Larger Compartment Model:
  if(n_CPT>3){
    stop("'formulaPKpar_MacroToMicro' works only for 1, 2 and 3 compartments.")
  }

  # Adjust name in PKmicroName:
  #   NOTE: Add the missing name in PKmicroName from PKmacroName
  idx_Missing <- which(!(names(PKmacroName) %in% names(PKmicroName)))
  for (idx_k in idx_Missing){
    PKmicroName[[names(PKmacroName)[idx_k]]] <- as.character(PKmacroName[idx_k])
  }

  # Re-ordered PKmicroName & PKparameters_Micro:
  Ordered <- c("Fabs0", "Tk0", "Fabs1", "ka", "Tlag1", "K", "Vc", "Kc1", "K1c", "Kc2", "K2c")
  idx_Ordered_1 <- intersect(Ordered, names(PKmicroName))
  PKmicroName <- PKmicroName[idx_Ordered_1]
  idx_Ordered_2 <- intersect(Ordered, names(PKparameters_Micro))
  PKparameters_Micro <- PKparameters_Micro[idx_Ordered_2]

  # Adjust name of PKparameters_Micro:
  idx_1 <- which(names(PKparameters_Micro) %in% names(PKmicroName))
  idx_2 <- which(names(PKmicroName) %in% names(PKparameters_Micro))
  names(PKparameters_Micro)[idx_1] <- PKmicroName[idx_2]

  # Output:
  return(PKparameters_Micro)
}


#' formulaPKpar_MacroToPrimary
#'
#' @description Converts provided macro PK parameters to primary parameters
#' @param PKparameters_Macro Named vector containing macro PK parameters and their values
#' @param PKmacroName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    Alpha = "Alpha", Beta = "Beta", Gamma = "Gamma", A = "A",
#'    B = "B", C = "C")
#' @param PKprimaryName Default: c(CL = "CL", Vc = "Vc", Q1 = "Q1", Vp1 = "Vp1", Q2 = "Q2", Vp2 = "Vp2")
#' @return Named vector containing primary PK parameters and their values
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_MacroToPrimary <- function(PKparameters_Macro,
                                        PKmacroName = c("Fabs0" = "Fabs0",
                                                        "Tk0"   = "Tk0"  ,
                                                        "Fabs1" = "Fabs1",
                                                        "ka"    = "ka",
                                                        "Tlag1" = "Tlag1",
                                                        "Alpha" = "Alpha",
                                                        "Beta"  = "Beta",
                                                        "Gamma" = "Gamma",
                                                        "A"     = "A",
                                                        "B"     = "B",
                                                        "C"     = "C"),
                                        PKprimaryName = c("CL"    = "CL",
                                                          "Vc"    = "Vc",
                                                          "Q1"    = "Q1",
                                                          "Vp1"   = "Vp1",
                                                          "Q2"    = "Q2",
                                                          "Vp2"   = "Vp2")) {

  # Use formulaPKpar_MacroToMicro:
  PKparameters_Micro <- formulaPKpar_MacroToMicro(PKparameters_Macro,
                                                  PKmacroName = PKmacroName)

  # Transform the parameter to Primary:
  PKparameters_Primary <- formulaPKpar_MicroToPrimary(PKparameters_Micro,
                                                      PKprimaryName = PKprimaryName)

  # Output:
  return(PKparameters_Primary)
}


#' formulaPKpar_MicroToMacro
#'
#' @description Converts provided micro PK parameters to macro parameters
#' @param PKparameters_Micro Named vector containing micro PK parameters and their values
#' @param PKmicroName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    K = "K", Vc = "Vc", Kc1 = "Kc1", K1c = "K1c", Kc2 = "Kc2",
#'    K2c = "K2c")
#' @param PKmacroName Default: c(Alpha = "Alpha", Beta = "Beta", Gamma = "Gamma", A = "A", B = "B",
#'    C = "C")
#' @return Named vector containing macro PK parameters and their values
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_MicroToMacro <- function(PKparameters_Micro,
                                      PKmicroName = c("Fabs0" = "Fabs0",
                                                      "Tk0"   = "Tk0"  ,
                                                      "Fabs1" = "Fabs1",
                                                      "ka"    = "ka",
                                                      "Tlag1" = "Tlag1",
                                                      "K"     = "K",
                                                      "Vc"    = "Vc",
                                                      "Kc1"   = "Kc1",
                                                      "K1c"   = "K1c",
                                                      "Kc2"   = "Kc2",
                                                      "K2c"   = "K2c"),
                                      PKmacroName = c("Alpha" = "Alpha",
                                                      "Beta"  = "Beta",
                                                      "Gamma" = "Gamma",
                                                      "A"     = "A",
                                                      "B"     = "B",
                                                      "C"     = "C")) {

  # Adjust Name in PKparameters_Micro:
  #   NOTE: Just to simplify the code afterward
  idx_Ordered <- intersect(PKmicroName, names(PKparameters_Micro))
  PKparameters_Micro <- PKparameters_Micro[idx_Ordered]
  idx_1 <- which(names(PKparameters_Micro) %in% PKmicroName)
  idx_2 <- which(PKmicroName %in% names(PKparameters_Micro))
  names(PKparameters_Micro)[idx_1] <- names(PKmicroName)[idx_2]

  # Generate output:
  PKparameters_Macro <- c()

  # PK Parameters to keep:
  if ("Fabs0" %in% names(PKparameters_Micro)){
    PKparameters_Macro[["Fabs0"]] <- PKparameters_Micro[["Fabs0"]]
  }
  if ("Tk0" %in% names(PKparameters_Micro)){
    PKparameters_Macro[["Tk0"]]   <- PKparameters_Micro[["Tk0"]]
  }
  if ("Fabs1" %in% names(PKparameters_Micro)){
    PKparameters_Macro[["Fabs1"]] <- PKparameters_Micro[["Fabs1"]]
  }

  if ("ka" %in% names(PKparameters_Micro)){
    PKparameters_Macro[["ka"]]    <- PKparameters_Micro[["ka"]]
  }

  if ("Tlag1" %in% names(PKparameters_Micro)){
    PKparameters_Macro[["Tlag1"]] <- PKparameters_Micro[["Tlag1"]]
  }

  # Determine the number of compartments:
  n_CPT <- get_NumberOfPKcpt(PKparameters_Micro, PKparaType = "Micro")

  # Estimate the macro parameters:
  if (n_CPT==1){
    PKparameters_Macro[["Alpha"]] <- PKparameters_Micro[["K"]]
    PKparameters_Macro[["A"]]     <- 1/PKparameters_Micro[["Vc"]]

  }else if(n_CPT==2){
    PKparameters_Macro[["Beta"]]  <- 1/2*(PKparameters_Micro[["Kc1"]] + PKparameters_Micro[["K1c"]] + PKparameters_Micro[["K"]] - sqrt((PKparameters_Micro[["Kc1"]] + PKparameters_Micro[["K1c"]] + PKparameters_Micro[["K"]])^2 - 4*PKparameters_Micro[["K1c"]]*PKparameters_Micro[["K"]]))
    PKparameters_Macro[["Alpha"]] <- PKparameters_Micro[["K1c"]]*PKparameters_Micro[["K"]]/PKparameters_Macro[["Beta"]]
    PKparameters_Macro[["A"]]     <- 1/PKparameters_Micro[["Vc"]]*(PKparameters_Micro[["K1c"]]-PKparameters_Macro[["Alpha"]])/(PKparameters_Macro[["Beta"]] -PKparameters_Macro[["Alpha"]])
    PKparameters_Macro[["B"]]     <- 1/PKparameters_Micro[["Vc"]]*(PKparameters_Micro[["K1c"]]-PKparameters_Macro[["Beta"]] )/(PKparameters_Macro[["Alpha"]]-PKparameters_Macro[["Beta"]] )

  }else if(n_CPT==3){
    a0    <- PKparameters_Micro[["K"]]*PKparameters_Micro[["K1c"]]*PKparameters_Micro[["K2c"]]
    a1    <- PKparameters_Micro[["K"]]*PKparameters_Micro[["K2c"]] + PKparameters_Micro[["K1c"]]*PKparameters_Micro[["K2c"]] + PKparameters_Micro[["K1c"]]*PKparameters_Micro[["Kc2"]] + PKparameters_Micro[["K"]]*PKparameters_Micro[["K1c"]] + PKparameters_Micro[["K2c"]]*PKparameters_Micro[["Kc1"]]
    a2    <- PKparameters_Micro[["K"]] + PKparameters_Micro[["Kc1"]] + PKparameters_Micro[["Kc2"]] + PKparameters_Micro[["K1c"]] + PKparameters_Micro[["K2c"]]
    p     <- a1 - a2^2/3
    q     <- 2*a2^3/27 - a1*a2/3 + a0
    r1    <- sqrt(-p^3/27)
    r2    <- 2*r1^(1/3)
    phi   <- acos(-q/(2*r1))/3
    PKparameters_Macro[["Alpha"]] <- -(cos(phi       )*r2 - a2/3)
    PKparameters_Macro[["Beta"]]  <- -(cos(phi+2*pi/3)*r2 - a2/3)
    PKparameters_Macro[["Gamma"]] <- -(cos(phi+4*pi/3)*r2 - a2/3)
    PKparameters_Macro[["A"]]     <- 1/PKparameters_Micro[["Vc"]] * (PKparameters_Micro[["K1c"]]-PKparameters_Macro[["Alpha"]])/(PKparameters_Macro[["Beta"]] -PKparameters_Macro[["Alpha"]]) * (PKparameters_Micro[["K2c"]]-PKparameters_Macro[["Alpha"]])/(PKparameters_Macro[["Gamma"]]-PKparameters_Macro[["Alpha"]])
    PKparameters_Macro[["B"]]     <- 1/PKparameters_Micro[["Vc"]] * (PKparameters_Micro[["K1c"]]-PKparameters_Macro[["Beta"]] )/(PKparameters_Macro[["Alpha"]]-PKparameters_Macro[["Beta"]] ) * (PKparameters_Micro[["K2c"]]-PKparameters_Macro[["Beta"]] )/(PKparameters_Macro[["Gamma"]]-PKparameters_Macro[["Beta"]] )
    PKparameters_Macro[["C"]]     <- 1/PKparameters_Micro[["Vc"]] * (PKparameters_Micro[["K1c"]]-PKparameters_Macro[["Gamma"]])/(PKparameters_Macro[["Beta"]] -PKparameters_Macro[["Gamma"]]) * (PKparameters_Micro[["K2c"]]-PKparameters_Macro[["Gamma"]])/(PKparameters_Macro[["Alpha"]]-PKparameters_Macro[["Gamma"]])
  }

  # For Larger Compartment Model:
  if(n_CPT>3){
    stop("'formulaPKpar_MicroToMacro' works only for 1, 2 and 3 compartments.")
  }

  # Adjust name in PKmacroName:
  #   NOTE: Add the missing name in PKmacroName from PKmicroName
  idx_Missing <- which(!(names(PKmicroName) %in% names(PKmacroName)))
  for (idx_k in idx_Missing){
    PKmacroName[[names(PKmicroName)[idx_k]]] <- as.character(PKmicroName[idx_k])
  }

  # Re-ordered PKmacroName & PKparameters_Macro:
  Ordered <- c("Fabs0", "Tk0", "Fabs1", "ka", "Tlag1", "Alpha", "Beta", "Gamma", "A", "B", "C")
  idx_Ordered_1 <- intersect(Ordered, names(PKmacroName))
  PKmacroName <- PKmacroName[idx_Ordered_1]
  idx_Ordered_2 <- intersect(Ordered, names(PKparameters_Macro))
  PKparameters_Macro <- PKparameters_Macro[idx_Ordered_2]

  # Adjust name of PKparameters_Macro:
  idx_1 <- which(names(PKparameters_Macro) %in% names(PKmacroName))
  idx_2 <- which(names(PKmacroName) %in% names(PKparameters_Macro))
  names(PKparameters_Macro)[idx_1] <- PKmacroName[idx_2]

  # Output:
  return(PKparameters_Macro)
}


#' formulaPKpar_MicroToPrimary
#'
#' @description Converts provided micro PK parameters to primary parameters
#' @param PKparameters_Micro Named vector containing micro PK parameters and their values
#' @param PKmicroName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    K = "K", Vc = "Vc", Kc1 = "Kc1", K1c = "K1c", Kc2 = "Kc2",
#'    K2c = "K2c")
#' @param PKprimaryName Default: c(CL = "CL", Vc = "Vc", Q1 = "Q1", Vp1 = "Vp1", Q2 = "Q2", Vp2 = "Vp2")
#' @return Named vector containing primary PK parameters and their values
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_MicroToPrimary <- function(PKparameters_Micro,
                                        PKmicroName = c("Fabs0" = "Fabs0",
                                                        "Tk0"   = "Tk0"  ,
                                                        "Fabs1" = "Fabs1",
                                                        "ka"    = "ka",
                                                        "Tlag1" = "Tlag1",
                                                        "K"     = "K",
                                                        "Vc"    = "Vc",
                                                        "Kc1"   = "Kc1",
                                                        "K1c"   = "K1c",
                                                        "Kc2"   = "Kc2",
                                                        "K2c"   = "K2c"),
                                        PKprimaryName = c("CL"    = "CL",
                                                          "Vc"    = "Vc",
                                                          "Q1"    = "Q1",
                                                          "Vp1"   = "Vp1",
                                                          "Q2"    = "Q2",
                                                          "Vp2"   = "Vp2")) {

  # Adjust Name in PKparameters_Micro:
  #   NOTE: Just to simplify the code afterward
  idx_Ordered <- intersect(PKmicroName, names(PKparameters_Micro))
  PKparameters_Micro <- PKparameters_Micro[idx_Ordered]
  idx_1 <- which(names(PKparameters_Micro) %in% PKmicroName)
  idx_2 <- which(PKmicroName %in% names(PKparameters_Micro))
  names(PKparameters_Micro)[idx_1] <- names(PKmicroName)[idx_2]

  # Generate output:
  PKparameters_Primary <- c()

  # PK Parameters to keep:
  if ("Fabs0" %in% names(PKparameters_Micro)){
    PKparameters_Primary[["Fabs0"]] <- PKparameters_Micro[["Fabs0"]]
  }
  if ("Tk0" %in% names(PKparameters_Micro)){
    PKparameters_Primary[["Tk0"]]   <- PKparameters_Micro[["Tk0"]]
  }
  if ("Fabs1" %in% names(PKparameters_Micro)){
    PKparameters_Primary[["Fabs1"]] <- PKparameters_Micro[["Fabs1"]]
  }

  if ("ka" %in% names(PKparameters_Micro)){
    PKparameters_Primary[["ka"]]    <- PKparameters_Micro[["ka"]]
  }

  if ("Tlag1" %in% names(PKparameters_Micro)){
    PKparameters_Primary[["Tlag1"]] <- PKparameters_Micro[["Tlag1"]]
  }

  # Determine the number of compartments:
  n_CPT <- get_NumberOfPKcpt(PKparameters_Micro, PKparaType = "Micro")

  # Estimate the Primary parameters:
  if (n_CPT>=1){
    PKparameters_Primary[["CL"]] <- PKparameters_Micro[["K"]]*PKparameters_Micro[["Vc"]]
    PKparameters_Primary[["Vc"]] <- PKparameters_Micro[["Vc"]]
  }
  if(n_CPT>=2){
    PKparameters_Primary[["Q1"]]  <- PKparameters_Micro[["Vc"]]*PKparameters_Micro[["Kc1"]]
    PKparameters_Primary[["Vp1"]] <- PKparameters_Micro[["Vc"]]*PKparameters_Micro[["Kc1"]]/PKparameters_Micro[["K1c"]]
  }
  if(n_CPT>=3){
    PKparameters_Primary[["Q2"]]  <- PKparameters_Micro[["Vc"]]*PKparameters_Micro[["Kc2"]]
    PKparameters_Primary[["Vp2"]] <- PKparameters_Micro[["Vc"]]*PKparameters_Micro[["Kc2"]]/PKparameters_Micro[["K2c"]]
  }

  # For Larger Compartment Model:
  if(n_CPT>3){
    stop("'formulaPKpar_MicroToPrimary' works only for 1, 2 and 3 compartments.")
  }

  # Adjust name in PKprimaryName:
  #   NOTE: Add the missing name in PKprimaryName from PKmicroName
  idx_Missing <- which(!(names(PKmicroName) %in% names(PKprimaryName)))
  for (idx_k in idx_Missing){
    PKprimaryName[[names(PKmicroName)[idx_k]]] <- as.character(PKmicroName[idx_k])
  }

  # Re-ordered PKprimaryName & PKparameters_Primary:
  Ordered <- c("Fabs0", "Tk0", "Fabs1", "ka", "Tlag1", "CL", "Vc", "Q1", "Vp1", "Q2", "Vp2")
  idx_Ordered_1 <- intersect(Ordered, names(PKprimaryName))
  PKprimaryName <- PKprimaryName[idx_Ordered_1]
  idx_Ordered_2 <- intersect(Ordered, names(PKparameters_Primary))
  PKparameters_Primary <- PKparameters_Primary[idx_Ordered_2]

  # Adjust name of PKparameters_Primary:
  idx_1 <- which(names(PKparameters_Primary) %in% names(PKprimaryName))
  idx_2 <- which(names(PKprimaryName) %in% names(PKparameters_Primary))
  names(PKparameters_Primary)[idx_1] <- PKprimaryName[idx_2]

  # Output:
  return(PKparameters_Primary)
}


#' formulaPKpar_PrimaryToMacro
#'
#' @description Converts provided primary PK parameters to macro parameters
#' @param PKparameters_Primary Named vector containing primary PK parameters and their values
#' @param PKprimaryName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    CL = "CL", Vc = "Vc", Q1 = "Q1", Vp1 = "Vp1", Q2 = "Q2",
#'    Vp2 = "Vp2")
#' @param PKmacroName Default: c(Alpha = "Alpha", Beta = "Beta", Gamma = "Gamma", A = "A", B = "B",
#'    C = "C")
#' @return Named vector containing macro PK parameters and their values
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_PrimaryToMacro <- function(PKparameters_Primary,
                                        PKprimaryName = c("Fabs0" = "Fabs0",
                                                          "Tk0"   = "Tk0"  ,
                                                          "Fabs1" = "Fabs1",
                                                          "ka"    = "ka",
                                                          "Tlag1" = "Tlag1",
                                                          "CL"    = "CL",
                                                          "Vc"    = "Vc",
                                                          "Q1"    = "Q1",
                                                          "Vp1"   = "Vp1",
                                                          "Q2"    = "Q2",
                                                          "Vp2"   = "Vp2"),
                                        PKmacroName = c("Alpha" = "Alpha",
                                                        "Beta"  = "Beta",
                                                        "Gamma" = "Gamma",
                                                        "A"     = "A",
                                                        "B"     = "B",
                                                        "C"     = "C")) {

  # Transform the parameter to Micro:
  PKparameters_Micro <- formulaPKpar_PrimaryToMicro(PKparameters_Primary,
                                                    PKprimaryName = PKprimaryName)

  # Use formulaPKpar_MicroToMacro:
  PKparameters_Macro <- formulaPKpar_MicroToMacro(PKparameters_Micro,
                                                  PKmacroName = PKmacroName)

  # Output:
  return(PKparameters_Macro)
}


#' formulaPKpar_PrimaryToMicro
#'
#' @description Converts provided primary PK parameters to micro parameters
#' @param PKparameters_Primary Named vector containing primary PK parameters and their values
#' @param PKprimaryName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    CL = "CL", Vc = "Vc", Q1 = "Q1", Vp1 = "Vp1", Q2 = "Q2",
#'    Vp2 = "Vp2")
#' @param PKmicroName Default: c(K = "K", Vc = "Vc", Kc1 = "Kc1", K1c = "K1c", Kc2 = "Kc2",
#'    K2c = "K2c")
#' @return Named vector containing micro PK parameters and their values
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_PrimaryToMicro <- function(PKparameters_Primary,
                                        PKprimaryName = c("Fabs0" = "Fabs0",
                                                          "Tk0"   = "Tk0"  ,
                                                          "Fabs1" = "Fabs1",
                                                          "ka"    = "ka",
                                                          "Tlag1" = "Tlag1",
                                                          "CL"    = "CL",
                                                          "Vc"    = "Vc",
                                                          "Q1"    = "Q1",
                                                          "Vp1"   = "Vp1",
                                                          "Q2"    = "Q2",
                                                          "Vp2"   = "Vp2"),
                                        PKmicroName = c("K"     = "K",
                                                        "Vc"    = "Vc",
                                                        "Kc1"   = "Kc1",
                                                        "K1c"   = "K1c",
                                                        "Kc2"   = "Kc2",
                                                        "K2c"   = "K2c")) {

  # Adjust Name in PKparameters_Primary:
  #   NOTE: Just to simplify the code afterward
  idx_Ordered <- intersect(PKprimaryName, names(PKparameters_Primary))
  PKparameters_Primary <- PKparameters_Primary[idx_Ordered]
  idx_1 <- which(names(PKparameters_Primary) %in% PKprimaryName)
  idx_2 <- which(PKprimaryName %in% names(PKparameters_Primary))
  names(PKparameters_Primary)[idx_1] <- names(PKprimaryName)[idx_2]

  # Generate output:
  PKparameters_Micro <- c()

  # PK Parameters to keep:
  if ("Fabs0" %in% names(PKparameters_Primary)){
    PKparameters_Micro[["Fabs0"]] <- PKparameters_Primary[["Fabs0"]]
  }
  if ("Tk0" %in% names(PKparameters_Primary)){
    PKparameters_Micro[["Tk0"]]   <- PKparameters_Primary[["Tk0"]]
  }
  if ("Fabs1" %in% names(PKparameters_Primary)){
    PKparameters_Micro[["Fabs1"]] <- PKparameters_Primary[["Fabs1"]]
  }

  if ("ka" %in% names(PKparameters_Primary)){
    PKparameters_Micro[["ka"]]    <- PKparameters_Primary[["ka"]]
  }

  if ("Tlag1" %in% names(PKparameters_Primary)){
    PKparameters_Micro[["Tlag1"]] <- PKparameters_Primary[["Tlag1"]]
  }

  # Determine the number of compartments:
  n_CPT <- get_NumberOfPKcpt(PKparameters_Primary, PKparaType = "Primary")

  # Estimate the Primary parameters:
  if (n_CPT>=1){
    PKparameters_Micro[["K"]]  <- PKparameters_Primary[["CL"]]/PKparameters_Primary[["Vc"]]
    PKparameters_Micro[["Vc"]] <- PKparameters_Primary[["Vc"]]
  }
  if(n_CPT>=2){
    PKparameters_Micro[["Kc1"]] <- PKparameters_Primary[["Q1"]]/PKparameters_Primary[["Vc"]]
    PKparameters_Micro[["K1c"]] <- PKparameters_Primary[["Q1"]]/PKparameters_Primary[["Vp1"]]
  }
  if(n_CPT>=3){
    PKparameters_Micro[["Kc2"]] <- PKparameters_Primary[["Q2"]]/PKparameters_Primary[["Vc"]]
    PKparameters_Micro[["K2c"]] <- PKparameters_Primary[["Q2"]]/PKparameters_Primary[["Vp2"]]
  }

  # For Larger Compartment Model:
  if(n_CPT>3){
    stop("'formulaPKpar_PrimaryToMicro' works only for 1, 2 and 3 compartments.")
  }

  # Adjust name in PKmicroName:
  #   NOTE: Add the missing name in PKmicroName from PKprimaryName
  idx_Missing <- which(!(names(PKprimaryName) %in% names(PKmicroName)))
  for (idx_k in idx_Missing){
    PKmicroName[[names(PKprimaryName)[idx_k]]] <- as.character(PKprimaryName[idx_k])
  }

  # Re-ordered PKmicroName & PKparameters_Micro:
  Ordered <- c("Fabs0", "Tk0", "Fabs1", "ka", "Tlag1", "K", "Vc", "Kc1", "K1c", "Kc2", "K2c")
  idx_Ordered_1 <- intersect(Ordered, names(PKmicroName))
  PKmicroName   <- PKmicroName[idx_Ordered_1]
  idx_Ordered_2 <- intersect(Ordered, names(PKparameters_Micro))
  PKparameters_Micro <- PKparameters_Micro[idx_Ordered_2]

  # Adjust name of PKparameters_Micro:
  idx_1 <- which(names(PKparameters_Micro) %in% names(PKmicroName))
  idx_2 <- which(names(PKmicroName) %in% names(PKparameters_Micro))
  names(PKparameters_Micro)[idx_1] <- PKmicroName[idx_2]

  # Output:
  return(PKparameters_Micro)
}


#' formulaPKpar_PrimaryToSecondary
#'
#' @description Converts provided primary PK parameters to secondary parameters
#' @param PKparameters_Primary Named vector containing primary PK parameters and their values
#' @param PKPrimaryName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    CL = "CL", Vc = "Vc", Q1 = "Q1", Vp1 = "Vp1", Q2 = "Q2",
#'    Vp2 = "Vp2")
#' @param PKsecondaryName Default: c(Vss = "Vss", Thalf0 = "Thalf0", Thalf1 = "Thalf1", Thalf2 = "Thalf2",
#'    Thalf3 = "Thalf3", MRT = "MRT")
#' @return Named vector containing secondary PK parameters and their values
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
formulaPKpar_PrimaryToSecondary <- function(PKparameters_Primary,
                                            PKPrimaryName = c("Fabs0" = "Fabs0",
                                                              "Tk0"   = "Tk0"  ,
                                                              "Fabs1" = "Fabs1",
                                                              "ka"    = "ka",
                                                              "Tlag1" = "Tlag1",
                                                              "CL"    = "CL",
                                                              "Vc"    = "Vc",
                                                              "Q1"    = "Q1",
                                                              "Vp1"   = "Vp1",
                                                              "Q2"    = "Q2",
                                                              "Vp2"   = "Vp2"),
                                            PKsecondaryName = c("Vss"    = "Vss",
                                                                "Thalf0" = "Thalf0",
                                                                "Thalf1" = "Thalf1",
                                                                "Thalf2" = "Thalf2",
                                                                "Thalf3" = "Thalf3",
                                                                "MRT"    = "MRT")) {

  # Adjust Name in PKparameters_Primary:
  #   NOTE: Just to simplify the code afterward
  idx_Ordered <- intersect(PKPrimaryName, names(PKparameters_Primary))
  PKparameters_Primary <- PKparameters_Primary[idx_Ordered]
  idx_1 <- which(names(PKparameters_Primary) %in% PKPrimaryName)
  idx_2 <- which(PKPrimaryName %in% names(PKparameters_Primary))
  names(PKparameters_Primary)[idx_1] <- names(PKPrimaryName)[idx_2]

  # Get Macro Parameters:
  PKparameters_Macro <- formulaPKpar_PrimaryToMacro(PKparameters_Primary)

  # Determine the number of compartments:
  n_CPT <- get_NumberOfPKcpt(PKparameters_Primary, PKparaType = "Primary")

  # Generate output:
  PKparameters_Secondary <- c()

  # Estimate Vss & Thalf1/2/3:
  if (n_CPT>=1){
    PKparameters_Secondary[["Vss"]]    <- PKparameters_Primary[["Vc"]]
    PKparameters_Secondary[["Thalf1"]] <- log(2)/PKparameters_Macro[["Alpha"]]

  }
  if(n_CPT>=2){
    PKparameters_Secondary[["Vss"]] <- PKparameters_Secondary[["Vss"]] + PKparameters_Primary[["Vp1"]]
    PKparameters_Secondary[["Thalf2"]] <- log(2)/PKparameters_Macro[["Beta"]]

  }
  if(n_CPT>=3){
    PKparameters_Secondary[["Vss"]] <- PKparameters_Secondary[["Vss"]] + PKparameters_Primary[["Vp2"]]
    PKparameters_Secondary[["Thalf3"]] <- log(2)/PKparameters_Macro[["Gamma"]]

  }

  # For Larger Compartment Model:
  if(n_CPT>3){
    stop("'formulaPKpar_PrimaryToSecondary' works only for 1, 2 and 3 compartments.")
  }

  # Estimate MRT:
  PKparameters_Secondary[["MRT"]] <- PKparameters_Secondary[["Vss"]]/PKparameters_Primary[["CL"]]

  # Estimate Apparent Half-Life:
  PKparameters_Secondary[["Thalf0"]] <- PKparameters_Secondary[["MRT"]]*log(2)

  # Adjust name in PKsecondaryName:
  #   NOTE: Add the missing name in PKsecondaryName from PKPrimaryName
  idx_Missing <- which(!(names(PKPrimaryName) %in% names(PKsecondaryName)))
  for (idx_k in idx_Missing){
    PKsecondaryName[[names(PKPrimaryName)[idx_k]]] <- as.character(PKPrimaryName[idx_k])
  }

  # Re-ordered PKsecondaryName & PKparameters_Secondary:
  Ordered <- c("Vss", "Thalf0", "Thalf1", "Thalf2", "Thalf3", "MRT")
  idx_Ordered_1 <- intersect(Ordered, names(PKsecondaryName))
  PKsecondaryName   <- PKsecondaryName[idx_Ordered_1]
  idx_Ordered_2 <- intersect(Ordered, names(PKparameters_Secondary))
  PKparameters_Secondary <- PKparameters_Secondary[idx_Ordered_2]

  # Adjust name of PKparameters_Secondary:
  idx_1 <- which(names(PKparameters_Secondary) %in% names(PKsecondaryName))
  idx_2 <- which(names(PKsecondaryName) %in% names(PKparameters_Secondary))
  names(PKparameters_Secondary)[idx_1] <- PKsecondaryName[idx_2]

  # Output:
  return(PKparameters_Secondary)
}
#' get_LagTime
#'
#' @description Estimates lag-time of drug effect from simulated time v parasitaemia data 
#' @param dataSim Simulation data generated by `sim_IQRmodel()` or equivalent
#' @param EMAX Numeric, value of EMAX. 
#' @param GR Numeric, value of growth rate. 
#' @param effCOL Character string denoting name of eff column in `dataSim`. Default: 'Eff'
#' @param Parasitemia Character string denoting name of parasitaemia column in `dataSim`.Default: "PL"
#' @param Plog Logical, is parasitaemia in `dataSim` in log-scale? Default: `TRUE`
#' @param lagTimeTYPE Character string: Which threshold to assess lag time against? Must be One of `Baseline`, `MIC`, `MPC90` or `1Log`. Default: "MIC" 
#' @return Estimate of lag time. 
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
get_LagTime <- function(dataSim,
                        GR,
                        EMAX,
                        effCOL      = "Eff",
                        Parasitemia = "PL",
                        Plog      = TRUE,
                        lagTimeTYPE = "MIC"
) {

  # Get GR from dataSim:
  if (!is.null(GR)) {
    dataSim$GR <- GR
  }else{
    if (is.null(dataSim$GR)) {
      stop("No GR value given.")
    }
  }

  # Get EMAX from dataSim:
  if (!is.null(EMAX)) {
    dataSim$EMAX <- EMAX
  }else{
    if (is.null(dataSim$EMAX)) {
      stop("No EMAX value given.")
    }
  }

  out <- data.frame(ID = unique(dataSim$ID), LagTime = NA)
  for (k in seq_along(out$ID)) {
    # Get subset of dataSim:
    dataSim_k <- subset(dataSim, ID == out$ID[k])

    # Check for times that are above the threshold:
    if (lagTimeTYPE=="MIC"){
      threshold <- dataSim_k$GR[1]/dataSim_k$EMAX[1]
      colToCheck <- effCOL

    }else if(lagTimeTYPE=="MPC90"){
      threshold <- 0.9
      colToCheck <- effCOL

    }else if (lagTimeTYPE=="1Log"){
      if (Plog){
        dataSim_k[,Parasitemia] <- exp(dataSim_k[,Parasitemia])
      }
      threshold              <- - dataSim_k[1,Parasitemia]/10
      dataSim_k[,Parasitemia] <- -dataSim_k[,Parasitemia]
      colToCheck             <- Parasitemia

    }else if (lagTimeTYPE=="Baseline"){
      idx_BelowBaseline <- which(dataSim_k[[Parasitemia]][dataSim_k$TIME<24*4]<dataSim_k[[Parasitemia]][1])
      lagtime_k <- dataSim_k$TIME[min(idx_BelowBaseline)]

    }else{
      stop("lagTimeTYPE' is not valid. It should be 'MIC', 'MPC90' or '1Log'.")
    }

    # Get time at which it is above the threshold:
    if (lagTimeTYPE!="Baseline"){
      idx <- which(dataSim_k[,colToCheck]>threshold)
      if (length(idx)==0) {
        warning(lagTimeTYPE, " not reached in simulations, lag time not defined.")
        out$LagTime[k] <- NA

      } else {
        idx_min        <- min(idx)
        lagtime_k      <- approx(dataSim_k[,colToCheck][idx_min+c(-1,0)], dataSim_k$TIME[idx_min+c(-1,0)], threshold)$y
        out$LagTime[k] <- lagtime_k
      }
    }else{
      out$LagTime[k] <- lagtime_k
    }

  }

  # Output:
  return(out)
}


#' get_NumberOfPKcpt
#'
#' @description From provided PK parameters, determine the number of compartments of a model. 
#' @param PKparameters Named vector of PK parameters and their values 
#' @param PKparaType Default: 'Primary'
#' @param PKprimaryName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    CL = "CL", Vc = "Vc", Q1 = "Q1", Vp1 = "Vp1", Q2 = "Q2",
#'    Vp2 = "Vp2")
#' @param PKmacroName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    Alpha = "Alpha", Beta = "Beta", Gamma = "Gamma", A = "A",
#'    B = "B", C = "C")
#' @param PKmicroName Default: c(Fabs0 = "Fabs0", Tk0 = "Tk0", Fabs1 = "Fabs1", ka = "ka", Tlag1 = "Tlag1",
#'    K = "K", Vc = "Vc", Kc1 = "Kc1", K1c = "K1c", Kc2 = "Kc2",
#'    K2c = "K2c")
#' @return Number of compartments of the model 
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
get_NumberOfPKcpt <- function(PKparameters,
                              PKparaType  = "Primary",     # Micro, Macro or Primary
                              PKprimaryName = c("Fabs0" = "Fabs0",
                                                "Tk0"   = "Tk0"  ,
                                                "Fabs1" = "Fabs1",
                                                "ka"    = "ka",
                                                "Tlag1" = "Tlag1",
                                                "CL"    = "CL",
                                                "Vc"    = "Vc",
                                                "Q1"    = "Q1",
                                                "Vp1"   = "Vp1",
                                                "Q2"    = "Q2",
                                                "Vp2"   = "Vp2"),
                              PKmacroName = c("Fabs0" = "Fabs0",
                                              "Tk0"   = "Tk0"  ,
                                              "Fabs1" = "Fabs1",
                                              "ka"    = "ka",
                                              "Tlag1" = "Tlag1",
                                              "Alpha" = "Alpha",
                                              "Beta"  = "Beta",
                                              "Gamma" = "Gamma",
                                              "A"     = "A",
                                              "B"     = "B",
                                              "C"     = "C"),
                              PKmicroName = c("Fabs0" = "Fabs0",
                                              "Tk0"   = "Tk0"  ,
                                              "Fabs1" = "Fabs1",
                                              "ka"    = "ka",
                                              "Tlag1" = "Tlag1",
                                              "K"     = "K",
                                              "Vc"    = "Vc",
                                              "Kc1"   = "Kc1",
                                              "K1c"   = "K1c",
                                              "Kc2"   = "Kc2",
                                              "K2c"   = "K2c")) {

  # Initialize the number of CPT:
  n_CPT <- 0

  if (toupper(PKparaType)=="PRIMARY"){

    # Adjust Name in PKparameters with PKprimaryName:
    #   NOTE: Just to simplify the code afterward
    idx_1 <- which(names(PKparameters) %in% PKprimaryName)
    idx_2 <- which(PKprimaryName %in% names(PKparameters))
    names(PKparameters)[idx_1] <- names(PKprimaryName)[idx_2]

    #   - One Compartment at least?
    if (("CL" %in% names(PKparameters)) && ("Vc" %in% names(PKparameters))){
      n_CPT <- 1
    }else{
      stop("'CL' and/or 'Vc' is missing in 'PKparameters_Primary': Please adjust 'PKparameters_Primary'.")
    }
    #   - Two Compartment at least?
    is_CPT2 <- FALSE
    if (("Q1" %in% names(PKparameters)) && ("Vp1" %in% names(PKparameters))){
      n_CPT   <- 2
      is_CPT2 <- TRUE
    }else if (xor(("Q1" %in% names(PKparameters)), ("Vp1" %in% names(PKparameters)))){
      stop("'Q1' or 'Vp1' is missing in 'PKparameters_Primary': Please adjust 'PKparameters_Primary'.")
    }
    #   - Three Compartment at least?
    is_CPT3 <- FALSE
    if (is_CPT2 && ("Q2" %in% names(PKparameters)) && ("Vp2" %in% names(PKparameters))){
      n_CPT   <- 3
      is_CPT3 <- TRUE
    }else if (xor(("Q2" %in% names(PKparameters)), ("Vp2" %in% names(PKparameters)))){
      stop("'Q2' or 'Vp2' is missing in 'PKparameters_Primary': Please adjust 'PKparameters_Primary'.")
    }else if(!is_CPT2 && ("Q2" %in% names(PKparameters)) && ("Vp2" %in% names(PKparameters))){
      stop("It seems that 'Q2' and 'Vp2' are present, but that 'Q1' and 'Vp1' are absent: Please add them to 'PKparameters_Primary'.")
    }


    # Check if some of the parameters are equal to 0:
    if (n_CPT==3){
      if (PKparameters[["Q2"]]==0){
        n_CPT <- 2
        warning("The parameters 'Q2' is equal to 0, therefore, it is assumed to be a 2-cpt model.")
      }
    }
    if (n_CPT==2){
      if (PKparameters[["Q1"]]==0){
        n_CPT <- 1
        warning("The parameters 'Q1' is equal to 0, therefore, it is assumed to be a 1-cpt model.")
      }
    }


  }else if (toupper(PKparaType)=="MACRO"){
    # Adjust Name in PKparameters with PKmacroName:
    #   NOTE: Just to simplify the code afterward
    idx_1 <- which(names(PKparameters) %in% PKmacroName)
    idx_2 <- which(PKmacroName %in% names(PKparameters))
    names(PKparameters)[idx_1] <- names(PKmacroName)[idx_2]

    #   - One Compartment at least?
    if (("Alpha" %in% names(PKparameters)) && ("A" %in% names(PKparameters))){
      n_CPT <- 1
    }else{
      stop("'Alpha' and/or 'A' is missing in 'PKparameters_Macro': Please adjust 'PKparameters_Macro'.")
    }
    #   - Two Compartment at least?
    is_CPT2 <- FALSE
    if (("Beta" %in% names(PKparameters)) && ("B" %in% names(PKparameters))){
      n_CPT   <- 2
      is_CPT2 <- TRUE
    }else if (xor(("Beta" %in% names(PKparameters)), ("B" %in% names(PKparameters)))){
      stop("'Beta' or 'B' is missing in 'PKparameters_Macro': Please adjust 'PKparameters_Macro'.")
    }
    #   - Three Compartment at least?
    is_CPT3 <- FALSE
    if (is_CPT2 && ("Gamma" %in% names(PKparameters)) && ("C" %in% names(PKparameters))){
      n_CPT   <- 3
      is_CPT3 <- TRUE
    }else if (xor(("Gamma" %in% names(PKparameters)), ("C" %in% names(PKparameters)))){
      stop("'Gamma' or 'C' is missing in 'PKparameters_Macro': Please adjust 'PKparameters_Macro'.")
    }else if(!is_CPT2 && ("Gamma" %in% names(PKparameters)) && ("C" %in% names(PKparameters))){
      stop("It seems that 'Gamma' and 'C' are present, but that 'Beta' and 'B' are absent: Please add them to 'PKparameters_Macro'.")
    }


  }else if (toupper(PKparaType)=="MICRO"){
    # Adjust Name in PKparameters with PKmicroName:
    #   NOTE: Just to simplify the code afterward
    idx_1 <- which(names(PKparameters) %in% PKmicroName)
    idx_2 <- which(PKmicroName %in% names(PKparameters))
    names(PKparameters)[idx_1] <- names(PKmicroName)[idx_2]

    #   - One Compartment at least?
    if (("K" %in% names(PKparameters)) && ("Vc" %in% names(PKparameters))){
      n_CPT <- 1
    }else{
      stop("'K' and/or 'Vc' is missing in 'PKparameters_Micro': Please adjust 'PKparameters_Micro'.")
    }
    #   - Two Compartment at least?
    is_CPT2 <- FALSE
    if (("Kc1" %in% names(PKparameters)) && ("K1c" %in% names(PKparameters))){
      n_CPT   <- 2
      is_CPT2 <- TRUE
    }else if (xor(("Kc1" %in% names(PKparameters)), ("K1c" %in% names(PKparameters)))){
      stop("'Kc1' or 'K1c' is missing in 'PKparameters_Micro': Please adjust 'PKparameters_Micro'.")
    }
    #   - Three Compartment at least?
    is_CPT3 <- FALSE
    if (is_CPT2 && ("Kc2" %in% names(PKparameters)) && ("K2c" %in% names(PKparameters))){
      n_CPT   <- 3
      is_CPT3 <- TRUE
    }else if (xor(("Kc2" %in% names(PKparameters)), ("K2c" %in% names(PKparameters)))){
      stop("'Kc2' or 'K2c' is missing in 'PKparameters_Micro': Please adjust 'PKparameters_Micro'.")
    }else if(!is_CPT2 && ("Kc2" %in% names(PKparameters)) && ("K2c" %in% names(PKparameters))){
      stop("It seems that 'Kc2' and 'K2c' are present, but that 'Kc1' and 'K1c' are absent: Please add them to 'PKparameters_Micro'.")
    }

    # Check if some of the parameters are equal to 0:
    if (n_CPT==3){
      if (PKparameters[["Kc2"]]==0){
        n_CPT <- 2
        warning("The parameters 'Kc2' is equal to 0, therefore, it is assumed to be a 2-cpt model.")
      }
    }
    if (n_CPT==2){
      if (PKparameters[["Kc1"]]==0){
        n_CPT <- 1
        warning("The parameters 'Kc1' is equal to 0, therefore, it is assumed to be a 1-cpt model.")
      }
    }


  }else{
    stop("'PKparaType' should be equal to 'Macro', 'Micro' or 'Primary'.")
  }

  # Output:
  return(n_CPT)
}


#' get_PRRtot
#'
#' @description Estimate total parasite reduction ratio from provided simulation data
#' @param dataSim Simulation data generated by `sim_IQRmodel()` or equivalent, containing at least a column named TIME and a column defined by `paraCOL`
#' @param paraCOL Character string denoting name of parasitaemia column in `dataSim`.Default: "PL"
#' @param Plog Logical, is parasitaemia in `dataSim` in log-scale? Default: `TRUE`
#' @param MethodMin Default: 'CubicSpline'
#' @return Estimated log10 total PRR 
#' @export
#' @importFrom MMVbase find_MinMMV
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
get_PRRtot <- function(dataSim,
                       paraCOL    = "PL",
                       Plog       = TRUE,
                       MethodMin  = "CubicSpline") {

  # Get Initial Parasitemia:
  Pstart <- dataSim[[paraCOL]][1]

  # Get Minimum Parasitemia:
  Pmin <- MMVbase::find_MinMMV(X  = dataSim$TIME,
                      Y  = dataSim[[paraCOL]],
                      dx = 0.1,
                      Method = MethodMin)$Ymin

  # Estimate PRRtotal:
  if (Plog) {
    PRRtot <- (Pstart-Pmin)/log(10)
  } else {
    PRRtot <- log10(Pstart/Pmin)
  }

  # Output:
  data.frame(PRRtot=PRRtot) # in log10 difference
}


#' getApparentKeys
#'
#' @description Estimate apparent MIC, MPC90, PRR24, PRR48 and PPR72 from provided simulation data
#' @param dataSim Simulation data generated by `sim_IQRmodel()` or equivalent
#' @param concCOL Character string denoting name of concentration column in `dataSim. `Default: 'Cc'
#' @param effCOL Character string denoting name of eff column in `dataSim`Default: 'Eff'
#' @param paraCOL Character string denoting name of parasitaemia column in `dataSim` Default: "PL" 
#' @param convConc Numeric: scalar for concentration data. Default: 1
#' @param Plog Logical, is parasitaemia in `dataSim` in log-scale? Default: `TRUE`
#' @return Apparent MIC, MPC90, PRR24, PRR48 and PPR72 
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
getApparentKeys <- function(dataSim,
                            concCOL  = "Cc",
                            effCOL   = "Eff",
                            paraCOL  = "PL",
                            convConc = 1,
                            Plog     = TRUE){

  # Apparent MIC: concentration when dP/dt=0, i.e., when P minimal
  idxMIC <- which.min(dataSim[[paraCOL]])
  if (idxMIC==length(dataSim[[paraCOL]])) {
    warning("Simulation time not long enough to determine MIC.")
    MIC   <- NA
  } else {
    if (idxMIC==1) {
      warning("First point minimum PL, skip MIC.")
      MIC <- NA
    } else {

      # FIRST ---
      # Get Time and Minimum Parasitemia happens:
      dMin <- MMVbase::find_MinMMV(X  = dataSim$TIME,
                          Y  = dataSim[[paraCOL]],
                          dx = 0.1,
                          Method = "CubicSpline")
      Tmin  <- dMin$Xmin
      PLmin <- dMin$Ymin

      # SECOND ---
      # Estimate the concentration at Tmin:
      if (Tmin<=dataSim$TIME[idxMIC]){
        slope <- (log(dataSim[[concCOL]][idxMIC-1]) - log(dataSim[[concCOL]][idxMIC]))/(dataSim$TIME[idxMIC-1] - dataSim$TIME[idxMIC])
      } else{
        slope <- (log(dataSim[[concCOL]][idxMIC+1]) - log(dataSim[[concCOL]][idxMIC]))/(dataSim$TIME[idxMIC+1] - dataSim$TIME[idxMIC])
      }
      logMIC <- slope*(Tmin - dataSim$TIME[idxMIC]) + log(dataSim[[concCOL]][idxMIC])
      MIC    <- exp(logMIC) * convConc
      tMIC   <- Tmin
    }
  }

  # Apparent MPC: concentration when Eff = 90% during washout
  if (all(dataSim[[effCOL]] < 0.9)) {
    MPC90 <- NA
  } else {
    idxAbove <- dataSim[[effCOL]] > 0.9
    idxLast  <- max(which(idxAbove))
    if (idxLast == length(dataSim[[effCOL]])) {
      warning("Simulation time not long enough to determine MPC90.")
      MPC90 <- NA
    } else {
      tMPC90 <- approx(x=dataSim[[effCOL]][c(idxLast,idxLast+1)], y=dataSim$TIME[c(idxLast,idxLast+1)], xout = 0.9)$y
      MPC90  <- approx(x=dataSim$TIME[c(idxLast,idxLast+1)]     , y=dataSim[[concCOL]][c(idxLast,idxLast+1)], xout = tMPC90)$y
      MPC90  <- MPC90*convConc
    }
  }

  # Apparent or rather effective PRR48:
  idx24 <- dataSim$TIME == 24
  idx48 <- dataSim$TIME == 48
  idx72 <- dataSim$TIME == 72
  idx0  <- dataSim$TIME ==  0
  if (Plog) {
    PRR24 <- (dataSim[[paraCOL]][idx0] - dataSim[[paraCOL]][idx24])/log(10)
    PRR48 <- (dataSim[[paraCOL]][idx0] - dataSim[[paraCOL]][idx48])/log(10)
    PRR72 <- (dataSim[[paraCOL]][idx0] - dataSim[[paraCOL]][idx72])/log(10)

  } else {
    PRR24 <- log(dataSim[[paraCOL]][idx0] / dataSim[[paraCOL]][idx24])/log(10)
    PRR48 <- log(dataSim[[paraCOL]][idx0] / dataSim[[paraCOL]][idx48])/log(10)
    PRR72 <- log(dataSim[[paraCOL]][idx0] / dataSim[[paraCOL]][idx72])/log(10)

  }

  # Output:
  out <- data.frame(MIC = MIC, MPC90 = MPC90, PRR24 = PRR24, PRR48 = PRR48, PRR72 = PRR72)
  return(out)
}

#' getCure
#'
#' @description Evaluate simulation data to see if cure has occured. 
#' @param dataSim Simulation data generated by `sim_IQRmodel()` or equivalent
#' @param paraCOL Character string denoting name of parasitaemia column in `dataSim` Default: "PL" 
#' @param LLOQ Numeric: LLOQ of parasitaemia measurements. Default: 10
#' @param CureThreshold Numeric: Cure threshold Default: 1/5000
#' @param Day Numeric: Which day to determine if cure has occured by. Default: 28
#' @return Vector of cure events or time at which <LLOQ occurred or maximum time, if cure does not occur. 
#' @export
#' @author Aline Fuchs (MMV)
#' @family Key Parameters
getCure <- function(dataSim,
                    paraCOL = "Parasitemia",
                    LLOQ = 10,
                    CureThreshold = 1/5000,
                    Day = 28) {
  # expects single time course of parasitemia contained in dataSim

  if (any(dataSim[[paraCOL]]<CureThreshold)) {
    Cure = TRUE
  } else {
    idxT <- dataSim$TIME == 24*Day
    Cure <- dataSim[[paraCOL]][idxT] < LLOQ
  }

  out <- data.frame(Cure=Cure)
  names(out) <- paste0("CureDay", Day)
  out

}
#' getKeysEMAX
#'
#' @description Estimate MIC, MPC90, PRR24, PRR48 and PRR72 from EMAX model 
#' @param x data.frame containing the columns `GR`, `EMAX`. `EC50` and `hill`
#' @param convConc Default: 1
#' @return Static MIC, MPC90, PRR24, PRR48 and PPR72 
#' @export
#' @author Aline Fuchs (MMV)
#' @importFrom MMVbase formula_EMAXmodelParsToMIC formula_EC50toMPC90 formula_EMAXtoPRR
#' 
#' @family Key Parameters
getKeysEMAX <- function(x, convConc = 1) {
  within(x, {
    MIC   <- MMVbase::formula_EMAXmodelParsToMIC(GR, EMAX, EC50, hill) * convConc
    MPC90 <- MMVbase::formula_EC50toMPC90(EC50, hill) * convConc
    PRR24 <- MMVbase::formula_EMAXtoPRR(GR, EMAX, timePRR = 24)
    PRR48 <- MMVbase::formula_EMAXtoPRR(GR, EMAX, timePRR = 48)
    PRR72 <- MMVbase::formula_EMAXtoPRR(GR, EMAX, timePRR = 72)
  })
}

#' getTimeAboveMIC
#'
#' @description Calculates the total time during which the concentration is above the MIC 
#' @param dataSim A data frame or data.table containing the simulation data.
#' @param MIC Numeric: Value of MIC. 
#' @param timeCOL A string specifying the name of the time column in `dataSim`. Default is `"TIME"`.
#' @param concCOL A string specifying the name of the concentration column in `dataSim`. Default is `"Cc"`.
#' @param convConc Numeric: scalar for concentration data. Default: 1
#' @return data.frame containing Time above MIC
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @importFrom MMVbase rectintMMV
#' @family Key Parameters
getTimeAboveMIC <- function(dataSim,
                            MIC, 
                            timeCOL  = "TIME",
                            concCOL  = "Cc",
                            convConc = 1) {

  # Transform Units if required ( eg unit, huRBP, DBScorrection)
  dataSim[[concCOL]] <- dataSim[[concCOL]] * convConc

  # Check whether MIC value given as argument or within dataset:
  if (is.null(MIC)) {
    if (is.null(dataSim$MIC)) {
      stop("No MIC value given.")
    } else {
      MIC <- unique(dataSim$MIC)
    }
  }

  # Check whether unique value given:
  if (length(MIC) != 1){
    stop("MIC not unique.")
  }

  # If MIC is Inf tMIC should be 0:
  if (MIC==Inf){
    out <- data.frame(tMIC = 0,
                      stringsAsFactors = FALSE)

  }else{
    # Check if Concentration is above MIC or not:
    dtMIC <- ifelse(dataSim[[concCOL]]>MIC, 1, 0)

    # Calculate Time Above MIC:
    tMIC <- MMVbase::rectintMMV(dataSim[[timeCOL]],dtMIC)

    #------------------------------------------------#
    # When dtMIC goes from 0 to 1 (or 1 to 0), a piece
    # of time between t[k-1] and t[k] is missing (Or added as left integration).
    # If the time step is very small, this
    # numerical estimation is not an issue. But if
    # the time scale is large, this could change
    # dramatically the value of tMIC. THEREFORE,
    # A correction is added assuming linearity when
    # going up, and exponential decrease when going
    # down between concentration points.
    #------------------------------------------------#

    # Estimate when Concentration goes above and below MIC:
    #   - BtoA stands for Below to Above
    #   - AtoB stands for Above to Below
    ddtMIC  <- c(0,dtMIC[2:(length(dtMIC))]-dtMIC[1:(length(dtMIC)-1)])
    idxBtoA <- which(ddtMIC==1)
    idxAtoB <- which(ddtMIC==-1)

    # Correction when it goes above:
    dt_BtoA <- (dataSim[[timeCOL]][idxBtoA] - dataSim[[timeCOL]][idxBtoA-1]) * (dataSim[[concCOL]][idxBtoA] - MIC) / (dataSim[[concCOL]][idxBtoA]-dataSim[[concCOL]][idxBtoA-1])
    tMIC    <- tMIC + sum(dt_BtoA)

    # Correction when it goes below:
    dt_AtoB <- (dataSim[[timeCOL]][idxAtoB] - dataSim[[timeCOL]][idxAtoB-1]) * log(MIC/dataSim[[concCOL]][idxAtoB]) / log(dataSim[[concCOL]][idxAtoB-1]/dataSim[[concCOL]][idxAtoB])
    tMIC    <- tMIC - sum(dt_AtoB)

    # Prepare output:
    out <- data.frame(tMIC = tMIC,
                      stringsAsFactors = FALSE)
  }

  # Output:
  out
}
#' getTimeAboveMICsim
#'
#' @description Calculates the total time during which the concentration is above the MIC from simulation data
#' @param dataSim A data frame or data.table containing the simulation data.
#' @param timeCOL A string specifying the name of the time column in `dataSim`. Default is `"TIME"`.
#' @param effCOL A string specifying the name of the eff column in `dataSim`. Default is `"Eff"`.
#' @param GR Numeric, value of growth rate. 
#' @param EMAX Numeric, value of EMAX. 
#' @param returnLag Logical: Return Lag time? Default: `FALSE`
#' @return data.frame containing time above MIC, and lag time if flagged
#' @export
#' @importFrom MMVbase trapzMMV
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
getTimeAboveMICsim <- function(dataSim,
                               timeCOL = "TIME",
                               effCOL  = "Eff",
                               GR      = NULL,
                               EMAX    = NULL,
                               returnLag = FALSE) {

  # Get GR from dataset:
  if (is.null(GR)) {
    if (is.null(dataSim$GR)) {
      stop("No GR value given.")
    } else {
      GR <- unique(dataSim$GR)
    }
  }

  # Get EMAX from dataset:
  if (is.null(EMAX)) {
    if (is.null(dataSim$EMAX)) {
      stop("No EMAX value given.")
    } else {
      EMAX <- unique(dataSim$EMAX)
    }
  }

  # Estimate Effect at which there is no growth or decrease in parasitemia:
  Eff.MIC <- GR/EMAX

  # Check if Concentration is above Eff.MIC or not:
  dtMIC <- ifelse(dataSim[[effCOL]]>Eff.MIC, 1, 0)

  # Calculate Time Above MIC:
  tMIC <- MMVbase::trapzMMV(dataSim[[timeCOL]],dtMIC)

  #------------------------------------------------#
  # When dtMIC goes from 0 to 1 (or 1 to 0), the
  # integration will add 0.5*(t[k+1]-t[k]), which
  # means that the changes happen in the middle
  # of the time interval. Which is not necessarly
  # the case. If the time step is very small, this
  # numerical estimation is not an issue. But if
  # the time scale is large, this could change
  # dramatically the value of tMIC. THEREFORE,
  # A correction is added.
  #------------------------------------------------#

  # Estimate when Effect goes above and below MIC:
  #   - BtoA stands for Below to Above
  #   - AtoB stands for Above to Below
  ddtMIC  <- c(0,dtMIC[2:(length(dtMIC))]-dtMIC[1:(length(dtMIC)-1)])
  idxBtoA <- which(ddtMIC==1)
  idxAtoB <- which(ddtMIC==-1)

  # Correction when it goes above:
  for (k in idxBtoA){
    dt_BtoA  <- dataSim[[timeCOL]][k] - approx(dataSim[[effCOL]][(k-1):k], dataSim[[timeCOL]][(k-1):k], Eff.MIC)$y
    dt_trapz <- 0.5*(dataSim[[timeCOL]][k]-dataSim[[timeCOL]][k-1])
    tMIC     <- tMIC + dt_BtoA - dt_trapz
  }

  # Correction when it goes below:
  for (k in idxAtoB){
    dt_AtoB  <- approx(dataSim[[effCOL]][(k-1):k], dataSim[[timeCOL]][(k-1):k], Eff.MIC)$y - dataSim[[timeCOL]][k-1]
    dt_trapz <- 0.5*(dataSim[[timeCOL]][k]-dataSim[[timeCOL]][k-1])
    tMIC     <- tMIC + dt_AtoB - dt_trapz
  }

  # Prepare Output:
  resultMIC <- data.frame(tMIC=tMIC)

  # Add Lag:
  if (returnLag){
    indMIC <- which((dataSim[[effCOL]]*EMAX-GR)>0)
    if (length(indMIC)==0){
      tLAG <- 0
      tMIClast <- 0
    }else{
      tLAG <- dataSim$TIME[indMIC[1]]
      tMIClast <- dataSim$TIME[indMIC[length(indMIC)]]
    }
    resultMIC$tLAG     <- tLAG
    resultMIC$tMIClast <- tMIClast
  }

  # Output:
  return(resultMIC)
}

#' getTimeAboveMPC90
#'
#' @description Calculates the total time during which the concentration is above the MPC90
#' @param dataSim A data frame or data.table containing the simulation data.
#' @param MPC90 Numeric: Value of MPC90
#' @param timeCOL A string specifying the name of the time column in `dataSim`. Default is `"TIME"`.
#' @param concCOL A string specifying the name of the concentration column in `dataSim`. Default is `"Cc"`.
#' @param convConc Numeric: scalar for concentration data. Default: 1
#' @return data.frame containing time above MPC90
#' @export
#' @importFrom MMVbase rectintMMV
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
getTimeAboveMPC90 <- function(dataSim,
                              MPC90,
                              timeCOL  = "TIME",
                              concCOL  = "Cc",
                              convConc = 1) {

  # Transform Units if required (e.g. unit, huRBP, DBScorrection)
  dataSim[[concCOL]] <- dataSim[[concCOL]] * convConc

  # Check whether MPC90 value given as argument or within dataset:
  if (is.null(MPC90)) {
    if (is.null(dataSim$MPC90)) {
      stop("No MPC90 value given.")
    } else {
      MPC90 <- unique(dataSim$MPC90)
    }
  }

  # Check whether unique value given:
  if (length(MPC90) != 1){
    stop("MPC90 not unique.")
  }

  # If MPC90 is Inf tMPC90 should be 0:
  if (MPC90==Inf){
    out <- data.frame(tMPC90 = 0,
                      stringsAsFactors = FALSE)

  }else{
    # Check if Concentration is above MPC90 or not:
    dtMPC90 <- ifelse(dataSim[[concCOL]]>MPC90, 1, 0)

    # Calculate Time Above MPC90:
    tMPC90 <- MMVbase::rectintMMV(dataSim[[timeCOL]],dtMPC90)

    #------------------------------------------------#
    # When dtMPC goes from 0 to 1 (or 1 to 0), a piece
    # of time between t[k-1] and t[k] is missing (Or added as left integration).
    # If the time step is very small, this
    # numerical estimation is not an issue. But if
    # the time scale is large, this could change
    # dramatically the value of tMIC. THEREFORE,
    # A correction is added assuming linearity when
    # going up, and exponential decrease when going
    # down between concentration points.
    #------------------------------------------------#

    # Estimate when Concentration goes above and below MPC90:
    #   - BtoA stands for Below to Above
    #   - AtoB stands for Above to Below
    ddtMPC90 <- c(0,dtMPC90[2:(length(dtMPC90))]-dtMPC90[1:(length(dtMPC90)-1)])
    idxBtoA  <- which(ddtMPC90==1)
    idxAtoB  <- which(ddtMPC90==-1)

    # Correction when it goes above:
    dt_BtoA <- (dataSim[[timeCOL]][idxBtoA] - dataSim[[timeCOL]][idxBtoA-1]) * (dataSim[[concCOL]][idxBtoA] - MPC90) / (dataSim[[concCOL]][idxBtoA]-dataSim[[concCOL]][idxBtoA-1])
    tMPC90  <- tMPC90 + sum(dt_BtoA)

    # Correction when it goes below:
    dt_AtoB <- (dataSim[[timeCOL]][idxAtoB] - dataSim[[timeCOL]][idxAtoB-1]) * log(MPC90/dataSim[[concCOL]][idxAtoB]) / log(dataSim[[concCOL]][idxAtoB-1]/dataSim[[concCOL]][idxAtoB])
    tMPC90  <- tMPC90 - sum(dt_AtoB)

    # Prepare output:
    out <- data.frame(tMPC90=tMPC90)
  }

  # Output:
  out
}


#' getTimeAboveMPC90sim
#'
#' @description Calculates the total time during which the concentration is above the MPC90 from simulation data
#' @param dataSim A data frame or data.table containing the simulation data.
#' @param timeCOL A string specifying the name of the time column in `dataSim`. Default is `"TIME"`.
#' @param effCOL A string specifying the name of the eff column in `dataSim`. Default is `"Eff"`.
#' @param returnLag Logical: Return Lag time? Default: `FALSE`
#' @return data.frame containing time above MIC, and lag time if flagged
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
getTimeAboveMPC90sim <- function(dataSim,
                                 timeCOL   = "TIME",
                                 effCOL    = "Eff",
                                 returnLag = FALSE) {

  # Check if Effect is above 0.9 or not:
  dtMPC90 <- ifelse(dataSim[[effCOL]]>0.9, 1, 0)

  # Calculate Time Above MPC90:
  tMPC90 <- MMVbase::trapzMMV(dataSim[[timeCOL]],dtMPC90)

  #------------------------------------------------#
  # When dtMPC90 goes from 0 to 1 (or 1 to 0), the
  # integration will add 0.5*(t[k+1]-t[k]), which
  # means that the changes happen in the middle
  # of the time interval. Which is not necessarly
  # the case. If the time step is very small, this
  # numerical estimation is not an issue. But if
  # the time scale is large, this could change
  # dramatically the value of tMPC90. THEREFORE,
  # A correction is added.
  #------------------------------------------------#

  # Estimate when Effect goes above and below 0.9:
  #   - BtoA stands for Below to Above
  #   - AtoB stands for Above to Below
  ddtMPC90 <- c(0,dtMPC90[2:(length(dtMPC90))]-dtMPC90[1:(length(dtMPC90)-1)])
  idxBtoA  <- which(ddtMPC90==1)
  idxAtoB  <- which(ddtMPC90==-1)

  # Correction when it goes above:
  for (k in idxBtoA){
    dt_BtoA  <- dataSim[[timeCOL]][k] - approx(dataSim[[effCOL]][(k-1):k], dataSim[[timeCOL]][(k-1):k], 0.9)$y
    dt_trapz <- 0.5*(dataSim[[timeCOL]][k]-dataSim[[timeCOL]][k-1])
    tMPC90   <- tMPC90 + dt_BtoA - dt_trapz
  }

  # Correction when it goes below:
  for (k in idxAtoB){
    dt_AtoB  <- approx(dataSim[[effCOL]][(k-1):k], dataSim[[timeCOL]][(k-1):k], 0.9)$y - dataSim[[timeCOL]][k-1]
    dt_trapz <- 0.5*(dataSim[[timeCOL]][k]-dataSim[[timeCOL]][k-1])
    tMPC90   <- tMPC90 + dt_AtoB - dt_trapz
  }

  # Prepare Output:
  resultMPC90 <- data.frame(tMPC90=tMPC90)

  # Add Lag:
  if (returnLag){
    indMPC <- which(dataSim[[effCOL]]>0.9)
    if (length(indMPC)==0){
      tLAG     <- 0
      tMPClast <- 0
    }else{
      tLAG     <- dataSim[[timeCOL]][indMPC[1]]
      tMPClast <- dataSim[[timeCOL]][indMPC[length(indMPC)]]
    }
    resultMPC90$tLAG     <- tLAG
    resultMPC90$tMPClast <- tMPClast
  }

  # Output:
  return(resultMPC90)
}


#' @title Get Time When Kill Rate is Above Growth Rate (tKRGR)
#'
#' @description
#' Calculates the total time during which the kill rate exceeds the growth rate in a simulation dataset.
#'
#' @param dataSim A data frame or data.table containing the simulation data.
#' @param timeCOL A string specifying the name of the time column in `dataSim`. Default is `"TIME"`.
#' @param killCOL A string specifying the name of the kill rate column in `dataSim`. Default is `"KillBlood"`.
#' @param GR A numeric value representing the growth rate threshold. Default is 0.07. Expects a single value, but if growth rate
#' is a column in `dataSim`, arguments such as `dataSim$GR[1]` or `unique(dataSim$GR)` are suitable.
#'
#' @details
#' The function computes the total time (`tKRGR`) during which the kill rate exceeds the growth rate in the provided simulation data.
#'
#' It handles cases where the kill rate crosses the growth rate threshold between time points by applying corrections using linear interpolation. This ensures a more accurate estimation of `tKRGR`, especially when time steps are large.
#'
#' @return
#' A data frame with one row and one column:
#' \describe{
#'   \item{`tKRGR`}{Numeric value representing the total time during which the kill rate exceeds the growth rate.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example simulation data
#' dataSim <- data.frame(
#'   TIME = seq(0, 24, by = 1),
#'   KillBlood = runif(25, min = 0, max = 1),
#'   GR = 0.5
#' )
#'
#' # Calculate tKRGR
#' resultKRGR <- getTimeKRaboveGR(
#'   dataSim = dataSim,
#'   timeCOL = "TIME",
#'   killCOL = "KillBlood",
#'   GR = 0.5
#' )
#'
#' print(resultKRGR)
#' }
#'
#' @export
#' @author Sam Jones (MMV)
#' @family Key Parameters
getTimeKRaboveGR<- function(dataSim,
                            timeCOL = "TIME",
                            killCOL  = "KillBlood",
                            GR      = 0.07) {

  # Get GR from dataset:
  if (is.null(GR)) {
    if (is.null(dataSim$GR)) {
      stop("No GR value given.")
    } else {
      GR <- unique(dataSim$GR)
    }
  }

  # Check if Concentration is above Eff.MIC or not:
  dtKRGR <- ifelse(dataSim[[killCOL]]>GR, 1, 0)
  # Calculate Time Above MIC:
  tKRGR <- MMVbase::trapzMMV(dataSim[[timeCOL]],dtKRGR)

  #------------------------------------------------#
  # When dtKRGR goes from 0 to 1 (or 1 to 0), the
  # integration will add 0.5*(t[k+1]-t[k]), which
  # means that the changes happen in the middle
  # of the time interval. Which is not necessarly
  # the case. If the time step is very small, this
  # numerical estimation is not an issue. But if
  # the time scale is large, this could change
  # dramatically the value of tKRGR THEREFORE,
  # A correction is added.
  #------------------------------------------------#

  # Estimate when Effect goes above and below dtKRGR:
  #   - BtoA stands for Below to Above
  #   - AtoB stands for Above to Below
  ddtKRGR  <- c(0,dtKRGR[2:(length(dtKRGR))]-dtKRGR[1:(length(dtKRGR)-1)])
  idxBtoA <- which(ddtKRGR==1)
  idxAtoB <- which(ddtKRGR==-1)

  # Correction when it goes above:
  for (k in idxBtoA){
    dt_BtoA  <- dataSim[[timeCOL]][k] - approx(dataSim[[killCOL]][(k-1):k], dataSim[[timeCOL]][(k-1):k], GR)$y
    dt_trapz <- 0.5*(dataSim[[timeCOL]][k]-dataSim[[timeCOL]][k-1])
    tKRGR     <- tKRGR + dt_BtoA - dt_trapz
  }

  # Correction when it goes below:
  for (k in idxAtoB){
    dt_AtoB  <- approx(dataSim[[killCOL]][(k-1):k], dataSim[[timeCOL]][(k-1):k], GR)$y - dataSim[[timeCOL]][k-1]
    dt_trapz <- 0.5*(dataSim[[timeCOL]][k]-dataSim[[timeCOL]][k-1])
    tKRGR     <- tKRGR + dt_AtoB - dt_trapz
  }

  # Prepare Output:
  resultKRGR <- data.frame(tKRGR=tKRGR)

  # Output:
  return(resultKRGR)
}

#' getTimeRecrudescence
#'
#' @description Estimate time of recrudescence 
#' @param x data.frame containing `"TIME"` and `valCol` where `valCol` is a measurement of parasitaemia
#' @param GR Numeric, value of growth rate. 
#' @param LLOQ Numeric: LLOQ of parasitaemia measurements. 
#' @param valCol A string specifying the name of the column containing parasitaemia measurements. Default is `"TIME"`. Default: 'DV'
#' @param log Logical, is parasitaemia in log-scale? Default: `TRUE`
#' @param Time2RecDef String cenoting methodology with which to determine recrudescence: One of `MMV` or `TAD` MMV = once min Para or LOQ is achieved. TAD = return to Para BL value: Default: 'MMV'
#' @return data.frame containing time of recrudesence 
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
getTimeRecrudescence <- function(x,
                                 GR,
                                 LLOQ,
                                 valCol      = "DV",
                                 log         = TRUE,
                                 Time2RecDef = "MMV")  # MMV = once min Para or LOQ is achieved # TAD = return to Para BL value
{
  # Expects time course of individual (x) with log-transformed parasite levels

  if (log) {
    cat("Input assumed to be log-transformed parasite levels.\n")
  } else {
    cat("Input assumed to be parasite levels on linear scale.\n -> LLOQ and parasite levels will be log-transformed.\n")
    LLOQ <- log(LLOQ)
    x[[valCol]] <- log(x[[valCol]])
  }

  # Remove NA values
  x <- x[!is.na(x[[valCol]]),]

  # If MMV definition for Time 2 Rec:
  if (Time2RecDef=="MMV"){
    # Determine time 2 recrudescence and whether that is only minimum estimate
    if (dim(x)[1]>0) { # Do only if sufficient data points available
      idxBLOQ <- x[[valCol]]<=LLOQ + abs(LLOQ)*1e-6
      if (sum(idxBLOQ)==0) {
        # if never below limit of quantification, take time of parasitemia nadir
        if (which.min(x[[valCol]])==length(x$TIME)){
          timeR   <- NA
          flagMin <- TRUE
        }else{
          timeR   <- x$TIME[which.min(x[[valCol]])]
          flagMin <- FALSE
        }
      } else {
        # if below limit of quantification, extrapolate time of recrudescence or maximum time if no recrudescence observed
        idxTR <- max(which(idxBLOQ))
        if (idxTR==length(x$TIME)){
          # A distinction is made between experiment with no recrudescence and with:
          #   A cut off of 10 days is used for the disctinction and experiments
          #   with no recrudescence are ignored
          if (max(x$TIME)<24*10){
            timeR   <- NA
            flagMin <- TRUE
          }else{
            timeR   <- max(x$TIME)
            flagMin <- TRUE
          }
        } else {
          timeR   <- x$TIME[idxTR+1] - (x[[valCol]][idxTR+1] - LLOQ)/GR
          flagMin <- FALSE
        }
      }
    } else { # Set to NA if not sufficient data points
      timeR   <- NA
      flagMin <- NA
    }

    # If TAD definition for Time 2 Rec:
  }else if (Time2RecDef=="TAD"){
    # Determine time 2 recrudescence and whether that is only minimum estimate
    if (dim(x)[1]>0) { # Do only if sufficient data points available
      # Data below LOQ:
      idxBLOQ <- x[[valCol]]<=LLOQ + abs(LLOQ)*1e-6

      # Data below Initial Baseline:
      idxBPL0 <- max(which(x[[valCol]]<=x[[valCol]][1]))

      # If no points under LLOQ
      if (sum(idxBLOQ)==0){
        # if never below limit of quantification, take time of parasitemia nadir
        if (which.min(x[[valCol]])==length(x$TIME)){
          timeR   <- NA
          flagMin <- TRUE
        }else if (which.min(x[[valCol]])==1){
          timeR   <- 0
          flagMin <- FALSE
        }else{
          if (idxBPL0==length(x$TIME)){
            nLast   <- length(x$TIME)
            timeR   <- x$TIME[nLast] + (x[[valCol]][1] - x[[valCol]][nLast])/GR
            flagMin <- FALSE
          }else{
            GR_Inter <- (x[[valCol]][idxBPL0+1] - x[[valCol]][idxBPL0])/(x$TIME[idxBPL0+1] - x$TIME[idxBPL0])
            timeR    <- x$TIME[idxBPL0] + (x[[valCol]][1] - x[[valCol]][idxBPL0])/GR_Inter
            flagMin  <- FALSE
          }
        }

      }else{
        # if the last point is above the LLOQ
        idxTR <- max(which(idxBLOQ))
        if (idxTR==length(x$TIME)){
          # A distinction is made between experiment with no recrudescence and with:
          #   A cut off of 10 days is used for the disctinction and experiments
          #   with no recrudescence are ignored
          if (max(x$TIME)<24*10){
            timeR   <- NA
            flagMin <- TRUE
          }else{
            timeR   <- max(x$TIME)
            flagMin <- TRUE
          }

        }else{
          if (idxBPL0==length(x$TIME)){
            nLast   <- length(x$TIME)
            timeR   <- x$TIME[nLast] + (x[[valCol]][1] - x[[valCol]][nLast])/GR
            flagMin <- FALSE
          }else{
            GR_Inter <- (x[[valCol]][idxBPL0+1] - x[[valCol]][idxBPL0])/(x$TIME[idxBPL0+1] - x$TIME[idxBPL0])
            timeR    <- x$TIME[idxBPL0] + (x[[valCol]][1] - x[[valCol]][idxBPL0])/GR_Inter
            flagMin  <- FALSE
          }
        }
      }
    }

  }else{
    stop("'Time2RecDef' should be equal to 'MMV' or 'TAD'.")
  }

  # Output
  out <- data.frame(Time2Recrud = unname(unlist(timeR)), flagMin=flagMin)
  out
}


#' plot_keyPDparametersMMV
#'
#' @description Plot key PD parameters from an IQRnlmeProjectMulti
#' @param object String denoting path to an IQRnlmeProjectMulti object
#' @param outputFolder String denoting folder in which to write plots 
#' @param ActivityPath Default: `NULL`
#' 
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
#' @importFrom plyr ldply
plot_keyPDparametersMMV <- function(object,
                                    outputFolder,
                                    ActivityPath = NULL) {

  # check inputs:
  if (!is_IQRnlmeProjectMulti(object))
    stop("Need to provide IQRnlmeProjectMulti object")
  if (!is.character(outputFolder) | length(outputFolder) != 1)
    stop("Need to provide single string as output folder")

  # Parse IQRnlmeProject Results for projects that actually have been run ... sort by "ORDER":
  RESULTS <- getResults_IQRnlmeProjectMulti(object)

  # General information:
  folders <- sapply(RESULTS,function(x) x$model)

  # Model information:
  foldersSplit1 <- aux_fileparts(folders)
  foldersSplit2 <- aux_fileparts(foldersSplit1$pathname)
  foldersSplit3 <- aux_fileparts(foldersSplit2$pathname)
  modelInfo <- data.frame(
    Folder= folders,
    Model = foldersSplit3$filename,
    Hill  = as.numeric(gsub("^[[:alpha:]_]*","",foldersSplit2$filename)),
    Run   = gsub("^[[:alpha:]_]*","",foldersSplit1$filename)
  )

  # Get model metrics:
  OBJ <- round(sapply(RESULTS,function(x) x$OBJ))
  AIC <- round(sapply(RESULTS,function(x) x$AIC))
  BIC <- round(sapply(RESULTS,function(x) x$BIC))

  # Get model parameters:
  ParsSelect <- c("GR","EMAX","EC50","CLPara", "kin", "ke", "IC50")
  estimates <- plyr::ldply(RESULTS, function(x) {
    idx <- x$parameternames %in% ParsSelect
    pars <- x$parametervalues[idx]
    se   <- x$stderrors[idx]
    data.frame(Folder = x$model, Parameter= c(names(pars),"BIC"), Value = c(pars,x$BIC), SE = c(se,NA))
  })

  # Merge model information to estimates:
  estimates <- merge(estimates, modelInfo)

  dodge <- position_dodge(width=0.5)

  # Plot:
  gr <- MMVbase::MMVggplot(estimates, aes(Hill, Value, color=Run, shape = !is.na(SE)), ActivityPath = ActivityPath) +
    geom_point(position=dodge) +
    geom_linerange(aes(ymin=Value-SE,ymax=Value+SE), position=dodge) +
    facet_wrap(~Parameter, scales="free_y") +
    scale_x_continuous(breaks = unique(estimates$Hill)) +
    # scale_color_manual(values=MMVbase::MMVcolors) +
    scale_shape_manual("SE Available", breaks=c(TRUE,FALSE), values=c(17,16))

  # Save Plot:
  IQRoutputPNG(gr, filename = file.path(outputFolder, "ModelParameterComparison"))
}


#' Evaluate APR
#'
#' @description Evaluate APR at Day `aprTime` given an individual time-course of parasitemia measurements,
#' baseline parasitemia, LLOQ parasitemia and weight at TIME 0. Two APRs are return, an optimistic one (CSoff - clinical symptoms off),
#' which assumes that no patient has fever, and a pessimistic (CSon - clinical symptoms on) one, which assumes that all patients have fever.
#' Note that times up to 7 days are hard-coded in and the function will crash (currently with no error message)
#' if shorter trials are evaluated
#'
#' @param dataSim data.table produced by `simulate_virtualtrials()`, specifically `$simPKPD`. `WT0` and `curethreshCOL` should also be included in this data.table if desired
#' @param LLOQ LLOQ parasitemia in the same scale as in `paraCOL` (i.e. linear or ln - numeric).
#' @param weightCOL optional column name of `dataSim` containing body weight (Default: `WT0`).
#' @param curethreshCOL optional column name of `dataSim` containing cure threshold (Default: "Cure_Thresh_col")
#' @param aprTime day that denotes final day of follow-up - metrics will be calculated based on this day
#' @param defaultWT default weight (kg) to be assumed if neither cure threshold or patient weight are provided (Default: 55)
#' @param extrapTol Numeric indicating the maximum allowed tolerance in hours for the interpolation/extrapolation in case
#' of last measurement not exactly at `aprTime`*24hr (Default: 48hr).
#' @param Plog Indicate if parasitemia is in the log or linear scale (Default: `TRUE` which means it is logged).
#' @param FLAGtruncateTime A logical indicating if the PL measurements should be truncated at `aprTime`*24hr,
#' for subjects with measurements later than `aprTime`*24hr.
#' If `TRUE`, all measurements after `aprTime`*24hr+`extrapTol` are removed (Default: `TRUE`).
#' @param FLAGextrapolateTime A logical indicating if the PL measurements should be extrapolated up to `aprTime`*24hr,
#' for subjects which measurements end before `aprTime`*24hr. If `TRUE` the last measured value is repeated up to
#' `aprTime`*24hr (Default: `FALSE`).
#' @param FLAGverbose Logical (Default: `FALSE`).
#'
#' @return A data.frame of a single row and columns USUBJID and binary (0 or 1) columns as follows:
#' \itemize{
#' \item USUBJID
#' \item TrialID
#' \item DOSE
#' \item PCloq: Time to first drop below LLOQ
#' \item Tcure: Time at which cure threshold is reached
#' \item P2BL : Presence of parasitaemia higher than baseline at 48 hours
#' \item P3BL : Presence of parasitaemia higher than 25% of baseline at 72 hours
#' \item P3 :  Presence of parasitaemia higher than LLOQ at 72 hours
#' \item ETF_CSon: Parasitaemia on 48h higher than baseline or parasitaemia > LLOQ on day 3 or parasitaemia > 25% baseline on day 3
#' \item ETF_CSoff : Parasitaemia on 48h higher than baseline or parasitaemia >25% baseline on 48h
#' \item simIPFabs : Absolute presence of parasites >LLOQ between day 4 and 7, irrespective of any other criteria being met
#' \item simLPFabs : Absolute presence of parasites >LLOQ between day 7 and aprTime, irrespective of any other criteria being met
#' \item LCFabs_CSon : Absolute presence of parasites >LLOQ between day 4 and aprTime, irrespective of any other criteria being met
#' \item LCF_CSon : Presence of parasitaemia between day 4 and aprTime, with respect to ETF_CSon
#' \item LPF_CSon : Presence of parasitaemia between day 7 and aprTime, with respect to ETF_CSon and LCF_CSon
#' \item LPF_CSoff : Presence of parasitaemia between day 7 and aprTime, with respect to ETF_CSoff
#' \item APR_CSon : Absence of parasitaemia on aprTime, if the patient did not meet criteria of ETF_CSon, LCF_CSon or LPF_CSon
#' \item APR_CSoff : Absence of parasitaemia on aprTime, if the patient did not meet criteria of ETF_CSoff or LPF_CSoff
#' \item APRinf_CSon : Attainment of the cure threshold within aprTime, in the absence of ETF_CSon
#' \item APRinf_CSoff : Attainment of the cure threshold within aprTime, in the absence of ETF_CSoff
#' }
#'
#' @export
#' @author Sam Jones (MMV), Venelin Mitov (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
evaluate_APR <- function(dataSim,
                         LLOQ,
                         timeCOL       = "TIME",
                         paraCOL       = "PL",
                         weightCOL     = "WT0",
                         curethreshCOL = "CureThreshold",
                         aprTime       = 28,
                         defaultWT     = 55,
                         extrapTol     = 48,
                         Plog          = TRUE,
                         FLAGtruncateTime    = TRUE,
                         FLAGextrapolateTime = FALSE,
                         FLAGverbose         = FALSE){

  #--------------------------------------#
  # STEP 1: Checks ----
  #--------------------------------------#

  # Make sure we have the right column:
  if(!(timeCOL %in% names(dataSim)) || !(paraCOL %in% names(dataSim))){
    stop("'", timeCOL, "', '", paraCOL, " 'are not columns available in 'dataSim': Please Adjust.")
  }

  # Check that WT0 exists:
  if(!is.null(weightCOL) && !(weightCOL %in% names(dataSim))){
    warning("'", weightCOL, "' is either not a column available in 'datasim' or NULL: Column 'WT0' was added with a default value of ", defaultWT)
    weightCOL   <- "WT0"
    dataSim$WT0 <- defaultWT
  }

  # Check that curethreshCOL exists:
  #   If cure threshold is note provided in dataSim,
  #   use calculations as below:.
  #   Cure Threshold: defined as 1 parasite in the body
  #   Cure Threshold (log) = log(1/( WT0* Vblood_mLperKg))
  #   With Vblood_mLperKg = 80 up to 35 kg and then 70 for BW > 35 kg
  #   if WT0 is also not provided, default WT will be used.
  if(!is.null(curethreshCOL) && !(curethreshCOL %in% names(dataSim))){
    warning("'", curethreshCOL, "' is either not a column available in 'datasim' or NULL: Columns 'CureThreshold' was added with default formula (see....).")
    curethreshCOL <- "CureThreshold"
    dataSim$CureThreshold <- ifelse(dataSim$WT0 <= 35,
                                    log(1/(dataSim$WT0*80)),
                                    log(1/(dataSim$WT0*70)))
  }

  # Convert to data.table and change column names:
  data.table::setDT(dataSim)
  data.table::setnames(dataSim,
                       c(timeCOL, paraCOL, weightCOL, curethreshCOL),
                       c("TIME" , "PL"   , "WT0"    , "CureThreshold"),
                       skip_absent = TRUE)
  dataSim <- dataSim[,list(TIME, PL, WT0, CureThreshold)]

  # Make sure you only have one body weight:
  if(length(unique(dataSim$WT0))>1){
    stop("'evaluate_APR' - More than one '", weightCOL, "' detected for subject ", dataSim$ID[1], ": Please Adjust.")
  }

  #   Log Parasitemia if necessary:
  #   NOTE: It is assumed that same transformation is necessary for LLOQ
  if(!Plog){
    dataSim$PL <- log(dataSim$PL)
    LLOQ       <- log(LLOQ)
  }

  # Make sure we have TIME ordered:
  if(is.unsorted(dataSim$TIME)) {
    dataSim <- dataSim[order(dataSim$TIME),]
  }

  # Check there are positive time:
  if (nrow(dataSim)==0 || !any(dataSim$TIME>0)) {
    if("ID" %in% names(dataSim)){
      msg <- paste0("'evaluate_APR' - ", dataSim$ID[1], " has no measurements after start of treatment: Returning NULL.")
    }else{
      msg <- "'evaluate_APR' has no measurements after start of treatment: Returning NULL."
    }
    warning(msg)
    return(NULL)
  }


  #--------------------------------------#
  # STEP 2: Define PLbase, Curethreshold and Cleaning ----
  #--------------------------------------#

  # Make sure that PLbase is set:
  T0 <- 0
  idx_TIMEleq0 <- which(dataSim$TIME <= 0)
  if(length(idx_TIMEleq0) > 0) {
    # Averaging done if several measurements at the latest time before 0:
    T0     <- max(dataSim$TIME[idx_TIMEleq0])
    PLbase <- mean(dataSim$PL[dataSim$TIME == T0], na.rm = TRUE)
  } else {
    warning("PLbase not specified and no PL measurements before or at TIME 0. Setting PLbase to first PL-value.")
    # Averaging done if several measurements at the earliest time.
    PLbase <- mean(dataSim$PL[dataSim$TIME == min(dataSim$TIME)], na.rm = TRUE)
  }

  # Remove pre-treatment data
  dataSim  <- dataSim[dataSim$TIME>=T0,]
  dataSim$TIME[dataSim$TIME == T0] <- 0
  LastTime <- max(dataSim$TIME)
  PLlast   <- dataSim$PL[which.max(dataSim$TIME)]

  # Curethreshold should be unique:
  CureThreshold <- unique(dataSim$CureThreshold)
  if(length(CureThreshold)!=1){
    stop(length(CureThreshold), " curethresholds are defined for the enetered individual: Only one value per individual is allowed")
  }

  # Set '<LLOQ'-infinite parasitemia to LLOQ - 0.0001
  dataSim$PL[dataSim$PL < LLOQ & is.infinite(dataSim$PL)] <- LLOQ - 0.0001

  # Take average at ties in TIME, i.e. time-points for which more than one measurement has been done.
  dataSim[, list(PL = mean(PL, na.rm = TRUE)), keyby = TIME]


  # Remove any datapoint after aprTime*24+extrapTol:
  if(FLAGtruncateTime){
    dataSim <- dataSim[dataSim$TIME<=(aprTime*24+extrapTol),]
  }


  #--------------------------------------#
  # STEP 3: Interpolation (+ Extrapolate) parasitemia at point of interest ----
  #--------------------------------------#

  if (length(dataSim$TIME)>=2) {
    # If any of the time points 48, 72, 96 and 'acprTime'*24 are missing, perform linear
    # interpolation at these time-points.
    interpTimes <- setdiff(c(48, 72, 96, aprTime*24), dataSim$TIME)
    if (!FLAGextrapolateTime) {
      # Do not do extrapolation of missing points after max(TIME) (Accounting for tolerence after 168hr).
      interpTimes <- interpTimes[interpTimes <= ifelse(max(dataSim$TIME)<=168, max(dataSim$TIME), max(dataSim$TIME) + extrapTol)]
    }
    if(length(interpTimes) > 0L) {
      interpResult <- approx(
        x = dataSim$TIME, y = dataSim$PL, xout = interpTimes, ties = list("ordered", mean), rule = 2)

      # Add interpolation result available data:
      dataSim <- data.table::rbindlist(list(dataSim,
                                            data.frame(TIME    = interpResult$x,
                                                       PL      = interpResult$y,
                                                       WT0     = dataSim$WT0[1],
                                                       CureThreshold = dataSim$CureThreshold[1],
                                                       stringsAsFactors = FALSE)))

      # Re-order needed:
      dataSim <- dataSim[order(dataSim$TIME),]
    }
  }else{
    warning("evaluate_APR - The dataset has less than two measurement points:",
            "Impossible to interpolate.")
  }

  # Remove data after 'acprTime' days:
  #		NOTE: LastTime might be modified for subjects for which the initial
  #			  LastTime > 'acprTime' days, it will return LastTime='acprTime'days ('acprTime'*24hours)
  dataSim  <- dataSim[dataSim$TIME <= aprTime*24,]
  LastTime <- max(dataSim$TIME)
  PLlast   <- dataSim$PL[which.max(dataSim$TIME)]


  #--------------------------------------#
  # STEP 4: Calculate criteria and adjust PL ----
  #--------------------------------------#

  # Establish whether LLOQ is reached, and time if so:
  #   NOTE: PCloq is time at which LLOQ is achieved

  if(all(dataSim$PL > LLOQ)){
    PCloq <- NA
  }else{
    PCloq <- min(dataSim$TIME[dataSim$PL <= LLOQ])
  }

  # Establish what time parasites reach cure threshold | Tcure <- NA if they do not.
  if(any(dataSim$PL <= CureThreshold )){
    Tcure <- min(dataSim$TIME[dataSim$PL < CureThreshold])
  }else{
    Tcure <- NA
  }

  # Duplicate PL data to modify based on CureThreshold.
  dataSim$PLcure <- as.double(dataSim$PL)

  # Modify PLcure based on TCure.
  if(!is.na(Tcure)){
    dataSim$PLcure[dataSim$TIME>=Tcure] <- CureThreshold
  }

  #+++++++#
  # 4-a: Initialize all conditions
  #+++++++#
  # Initialize all conditions as FALSE:
  # These could all be the same line, but this way helps clarity when writing this function.
  P2BL <- P3BL <- P3 <- FALSE
  ETF_CSon <- ETF_CSoff <- FALSE
  simIPFabs <- simLPFabs <- LCFabs_CSon <- FALSE
  LCF_CSon <- LPF_CSon <- LPF_CSoff <- FALSE
  APR_CSon <- APR_CSoff <- FALSE
  APRinf_CSon <- APRinf_CSoff <- FALSE

  # ETF, IPF, LPF, APR calculations can then occur using dataSim$PLcure

  #+++++++#
  # 4-b: Classifications based on parasitaemia within first 3 days
  #+++++++#

  # P2BL : Is parasitaemia at 48h higher than at baseline?
  if(dataSim$PLcure[dataSim$TIME==48] > dataSim$PL[1]){
    P2BL <- TRUE
  }

  # P3BL : Is parasitaemia at 72h higher than 25% baseline?
  # Note : PL (and PLcure) are log-transformed.
  # .: log(PLbase) + log(0.25) is = PL*0.25
  if(dataSim$PLcure[dataSim$TIME==72] >= dataSim$PL[1]+log(0.25)){
    P3BL <- TRUE
  }

  # P3 : Is parasitaemia at 72 higher than LLOQ?
  if(dataSim$PLcure[dataSim$TIME==72] >= LLOQ){
    P3 <- TRUE
  }

  #+++++++#
  # 4-c: ETF classifications
  #+++++++#

  # ETF_CSon : Less permissive criteria - TRUE if any of P2BL, P3BL or P3 are met.
  # Note: In practice, we would assume that P3BL must always be TRUE if P3 is TRUE?
  if(P2BL == TRUE | P3BL == TRUE | P3 == TRUE){
    ETF_CSon <- TRUE
  }

  # ETF_CSoff : More permissive criteria - TRUE if any of P2BL, P3BL are met, but not P3.
  if(P2BL == TRUE | P3BL == TRUE){
    ETF_CSoff <- TRUE
  }

  #+++++++#
  # 4-d: "Absolute classifications" - based on simulated parasitaemia (PLcure) without respect to any other classification being declared
  #+++++++#

  # Intermediate parasitological failure
  # IPFabs : Presence of parasites day 4-7, irrespective of other criteria
  if(any(dataSim$PLcure[dataSim$TIME >=96 & dataSim$TIME < 168] >= LLOQ)){
    simIPFabs <- TRUE
  }else if(max(dataSim$TIME) < 168){
    simIPFabs <- NA
  }

  # Late parasitological failure
  # LPFabs : Presence of parasites day 7-aprTime*24, irrespective of other criteria
  if(any(dataSim$PLcure[dataSim$TIME >=168 & dataSim$TIME <= aprTime*24] >= LLOQ)){
    simLPFabs <- TRUE
  }else if(max(dataSim$TIME) < aprTime*24){
    simLPFabs <- NA
  }

  # Late clinical failure CSon scenario LCFabs_CSon:
  #   If simIPF or simLPFabs == TRUE
  #   i.e. presence of parasites between day 4 and final day.
  if(!is.na(simIPFabs) && simIPFabs == TRUE){
    LCFabs_CSon <- TRUE
  }else if(!is.na(simLPFabs) && simLPFabs == TRUE){
    LCFabs_CSon <- TRUE
  }else if(is.na(simIPFabs) || is.na(simLPFabs)){
    LCFabs_CSon <- NA
  }

  #+++++++#
  # 4-e: Late failure (clinical or parasitological)
  #+++++++#

  # LCF_CSon : Late clinical failure in CSon scenario - presence of parasites between day 4 and final day, with respect to ETF.
  if(ETF_CSon == FALSE && !is.na(LCFabs_CSon) && LCFabs_CSon == TRUE){
    LCF_CSon <- TRUE
  }else if(is.na(LCFabs_CSon)){
    LCF_CSon <- NA
  }

  # Late parasitological failure with fever LPF_Fever:
  #   Presence of parasites between day 7 and final day, with respect to ETF_CSon and LCF.
  #   The difference between LCF and LPF is somewhat confusing at first glance; a failure day 4-7
  #   is LCF, not LPF and precludes an LPF.
  if(ETF_CSon == FALSE && !is.na(LCF_CSon) && LCF_CSon == FALSE && !is.na(simLPFabs) && simLPFabs == TRUE){
    LPF_CSon <- TRUE
  }else if(ETF_CSon == FALSE && !is.na(LCF_CSon) && LCF_CSon == FALSE && is.na(simLPFabs)){
    LPF_CSon <- NA
  }else if(ETF_CSon == FALSE && is.na(LCF_CSon)){
    LPF_CSon <- NA
  }

  # Late parasitological failure in CSoff scenario LPF_CSoff:
  #   Presence of parasites between day 7 and final day, with respect to ETF_noFever
  #   Also without respect to LCF_CSoff (unlike LPF_CSon!)
  if(ETF_CSoff == FALSE && !is.na(simLPFabs) && simLPFabs == TRUE){
    LPF_CSoff <- TRUE
  }else if(ETF_CSoff == FALSE && is.na(simLPFabs)){
    LPF_CSoff <- NA
  }

  #+++++++#
  # 4-f: APR
  #+++++++#

  # Adequate parasitological response:
  # APR_CSon - Absence of parasitaemia on final day, with respect to ETF_CSon, LCF_CSon or LPF_CSon
  if(ETF_CSon == FALSE && !is.na(LCF_CSon) && LCF_CSon == FALSE && !is.na(LPF_CSon) && LPF_CSon == FALSE){
    APR_CSon <- TRUE
  }else if(ETF_CSon == FALSE && !is.na(LPF_CSon) && LCF_CSon == FALSE && is.na(LPF_CSon)){
    APR_CSon <- NA
  }else if(ETF_CSon == FALSE && is.na(LCF_CSon)){
    APR_CSon <- NA
  }

  # APR_CSoff - Absence of parasitaemia on final day, with respect to ETF_CSoff and LPF_CSoff
  if(ETF_CSoff == FALSE && !is.na(LPF_CSoff) && LPF_CSoff == FALSE){
    APR_CSoff <- TRUE
  }else if(ETF_CSoff == FALSE && is.na(LPF_CSoff)){
    APR_CSoff <- NA
  }

  #+++++++#
  # 4-g: APRinf
  #+++++++#

  # Use Tcure to establish if parasites reach the threshold.
  # Tcure is NA if they don't reach the threshold within aprTime....
  # So in principle, this is simple. But we may want to later change the code to call specific days?
  # We also may later want to see if Parasites at aprTime are lower than aprTime-1 (i.e, decreasing)

  #APRinf_CSon
  if(!is.na(Tcure) && ETF_CSon == FALSE){
    APRinf_CSon <- TRUE
  }

  #APRinf_CSoff
  if(!is.na(Tcure) && ETF_CSoff == FALSE){
    APRinf_CSoff <- TRUE
  }

  #--------------------------------------#
  # STEP 5: Prepare output ----
  #--------------------------------------#

  # Update names after creating data frame to paste in aprTime*24 for all LPF, LCF, APR criteria.
  # Coercing numerical to as.double( ) appears to be needed to stop data.table operations (i.e., [ , evaluate_apr(), ]) from crashing in the workflow?
  # See this stack overflow:
  # https://stackoverflow.com/questions/12125364/why-does-median-trip-up-data-table-integer-versus-double
  # Although I am not sure which particular part of the output is causing issues.

  res <- data.frame(PLbase    = as.double(PLbase),
                    PCloq     = as.double(PCloq),
                    Tcure     = as.double(Tcure),
                    P2BL      = P2BL,
                    P3BL      = P3BL,
                    P3        = P3,
                    ETF_CSon  = ETF_CSon,
                    ETF_CSoff = ETF_CSoff,
                    simIPFabs = simIPFabs,
                    simLPFabs = simLPFabs,
                    LCFabs_CSon  = LCFabs_CSon,
                    LCF_CSon     = LCF_CSon,
                    LPF_CSon     = LPF_CSon,
                    LPF_CSoff    = LPF_CSoff,
                    APR_CSon     = APR_CSon,
                    APR_CSoff    = APR_CSoff,
                    APRinf_CSon  = APRinf_CSon,
                    APRinf_CSoff = APRinf_CSoff,
                    stringsAsFactors   = FALSE)

  # Update names to include, where appropriate, aprTime
  names(res) <- c("PLbase","PCloq","Tcure",
                  "P2BL","P3BL","P3","ETF_CSon","ETF_CSoff","simIPFabs",
                  paste0("simLPF",aprTime,"abs"),
                  paste0("LCF",aprTime,"abs_CSon"),
                  paste0("LCF",aprTime,"_CSon") ,
                  paste0("LPF",aprTime,"_CSon") ,
                  paste0("LPF",aprTime,"_CSoff"),
                  paste0("APR",aprTime,"_CSon"),
                  paste0("APR",aprTime,"_CSoff"),
                  "APRinf_CSon", "APRinf_CSoff")

  # Add attributes:
  attr(res,'aprTime')             <- aprTime
  attr(res,'extrapTol')           <- extrapTol
  attr(res,'FLAGtruncateTime')    <- FLAGtruncateTime
  attr(res,'FLAGextrapolateTime') <- FLAGextrapolateTime

  # Output:
  res
}


#' Evaluate APR at Day 28
#'
#' Evaluate APR at Day 28 given an individual time-course of parasitemia measurements,
#' baseline parasitemia, LLOQ parasitemia and weight at TIME 0. Two APRs are return, an Optimistic (CSoff - clinical symptoms off) one,
#' which assumes that no patient has fever, and a pessmistic (CSon - clinical symptons on) one, which assumes that all patients have fever.
#' Note that times up to 7 days are hard-coded in and the function will crash (currently with no error message)
#' if shorter trials are evaluated
#'
#' @param dataSim data.table produced by `simulate_virtualtrials()`, specifically `$simPKPD`. `WT0` and `curethreshCOL` should also be included in this data.table if desired
#' @param LLOQ LLOQ parasitemia in the same scale as in `paraCOL` (i.e. linear or ln - numeric).
#' @param weightCOL optional column name of `dataSim` containing body weight (Default: `WT0`).
#' @param curethreshCOL optional column name of `dataSim` containing cure threshold (Default: "Cure_Thresh_col")
#' @param aprTime day that denotes final day of follow-up - metrics will be calculated based on this day
#' @param defaultWT default weight (kg) to be assumed if neither cure threshold or patient weight are provided (Default: 55)
#' @param extrapTol Numeric indicating the maximum allowed tolerance in hours for the interpolation/extrapolation in case
#' of last measurement not exactly at 28*24hr (Default: 48hr).
#' @param Plog Indicate if parasitemia is in the log or linear scale (Default: `TRUE` which means it is logged).
#' @param FLAGtruncateTime A logical indicating if the PL measurements should be truncated at 28*24hr,
#' for subjects with measurements later than 28*24hr.
#' If `TRUE`, all measurements after 28*24hr+`extrapTol` are removed (Default: `TRUE`).
#' @param FLAGextrapolateTime A logical indicating if the PL measurements should be extrapolated up to 28*24hr,
#' for subjects which measurements end before 28*24hr. If `TRUE` the last measured value is repeated up to
#' 28*24hr (Default: `FALSE`).
#' @param FLAGverbose Logical (Default: `FALSE`).
#'
#' @return A data.frame of a single row and columns USUBJID and binary (0 or 1) columns as follows:
#' \itemize{
#' \item USUBJID
#' \item TrialID
#' \item DOSE
#' \item PCloq: Time to first drop below LLOQ
#' \item Tcure: Time at which cure threshold is reached
#' \item P2BL : Presence of parasitaemia higher than baseline at 48 hours
#' \item P3BL : Presence of parasitaemia higher than 25% of baseline at 72 hours
#' \item P3 :  Presence of parasitaemia higher than LLOQ at 72 hours
#' \item ETF_CSon : Parasitaemia on 48h higher than baseline or parasitaemia > LLOQ on day 3 or parasitaemia > 25% baseline on day 3
#' \item ETF_CSoff : Parasitaemia on 48h higher than baseline or parasitaemia >25% baseline on 48h
#' \item simIPFabs : Absolute presence of parasites >LLOQ between day 4 and 7, irrespective of any other criteria being met
#' \item simLPF28abs : Absolute presence of parasites >LLOQ between day 7 and 28, irrespective of any other criteria being met
#' \item LCF28abs_CSon : Absolute presence of parasites >LLOQ between day 4 and 28, irrespective of any other criteria being met
#' \item LCF28_CSon : Presence of parasitaemia between day 4 and 28, with respect to ETF_CSon
#' \item LPF28_CSon : Presence of parasitaemia between day 7 and 28, with respect to ETF_CSon and LCF28_CSon
#' \item LPF28_CSoff : Presence of parasitaemia between day 7 and 28, with respect to ETF_CSoff
#' \item APR28_CSon : Absence of parasitaemia on day 28, if the patient did not meet criteria of ETF_CSon, LCF28_CSon or LPF28_CSon
#' \item APR28_CSoff : Absence of parasitaemia on day 28, if the patient did not meet criteria of ETF_CSoff or LPF28_CSoff
#' }
#'
#' @export
#' @author Sam Jones (MMV), Venelin Mitov (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family Key Parameters
evaluate_APR28 <- function(dataSim,
                           LLOQ,
                           timeCOL       = "TIME",
                           paraCOL       = "PL",
                           weightCOL     = "WT0",
                           curethreshCOL = "CureThreshold",
                           defaultWT     = 55,
                           extrapTol     = 48,
                           Plog          = TRUE,
                           FLAGtruncateTime    = TRUE,
                           FLAGextrapolateTime = FALSE,
                           FLAGverbose         = FALSE){

  # use 'evaluate_APR' with 'aprTime'=28:
  res <- evaluate_APR(dataSim       = dataSim,
                      LLOQ          = LLOQ,
                      timeCOL       = timeCOL,
                      paraCOL       = paraCOL,
                      weightCOL     = weightCOL,
                      curethreshCOL = curethreshCOL,
                      aprTime       = 28,
                      defaultWT     = defaultWT,
                      extrapTol     = extrapTol,
                      Plog          = Plog,
                      FLAGtruncateTime    = FLAGtruncateTime,
                      FLAGextrapolateTime = FLAGextrapolateTime,
                      FLAGverbose         = FLAGverbose)

  # Output:
  res
}
#' ConfusionMatrix
#' @description Calculate confusion matrix and derived statistics from a predicted and reference binary data
#' @param data a data.frame with at least two columns of logical values named as the arguments colPred and colRef.
#' @param colPred,colRef character strings denoting names of logical (TRUE/FALSE) columns in data.
#' @param positiveClass character string denoting the human readable meaning of TRUE values in data[[colPred]] and data[[colRef]].
#' @param negativeClass character string (default paste0("No ", positiveClass)), denoting the human readable meaning of FALSE
#' values in data[[colPred]] and data[[colRef]].
#' @return a named list as described in documentation of `[caret::confusionMatrix()]`.
#' @note This function is using `[caret::confusionMatrix()]`
#' @export
ConfusionMatrix <- function(data, colPred, colRef, positiveClass, negativeClass = paste0('No ', positiveClass)) {

  if(!is.logical(data[[colPred]]) || !is.logical(data[[colRef]])) {
    stop("colPred and colRef should be logical columns in data.")
  }

  # map 'TRUE' or 'FALSE' values of colPred and colRef to negative and positive class respectively.
  logicalToClass = c('TRUE' = positiveClass, 'FALSE' = negativeClass)

  # predicted classes
  predicted <- logicalToClass[as.character(as.logical(data[[colPred]]))]
  # observed or other type of reference
  reference <- logicalToClass[as.character(as.logical(data[[colRef]]))]

  res <- caret::confusionMatrix(data      = factor(predicted, levels = unname(logicalToClass), labels = unname(logicalToClass)),
                                reference = factor(reference, levels = unname(logicalToClass), labels = unname(logicalToClass)),
                                positive  = positiveClass,
                                mode = "everything")
  res$negative <- negativeClass

  res
}

#' SummarizeConfusionMatrices
#' @description Summarize a list of confusion matrices into a format suitable for the IQRoutputTable function.
#' If the list contains one confusion matrix, the matrix TP, FN, FP, TN values are returned in the form of a data.frame.
#' If the list contains more than one confusion matrix a summary matrix is returned containing the median, Q05 and Q95 quantiles of
#' TP, FN, FP, TN.
#' @param cmList a list of one or more confusionMatrix objects created by calling the function ConfusionMatrix.
#' @return a data.frame of 3 columns named according to colNames.
#' @export
SummarizeConfusionMatrices <- function(cmList) {
  if(is.null(cmList) || length(cmList) == 0) {
    NULL
  } else if(!is.list(cmList)) {
    stop("SummarizeConfusionMatrices: argument cmList should be a list of confusionMatrix objects created by calling ConfusionMatrix().")
  } else {

    # . Confusion matrix summary ----
    dt <- data.table::rbindlist(lapply(cmList, function(x) {
      if(!inherits(x, "confusionMatrix")) {
        stop("SummarizeConfusionMatrices: Argument cmList is a list but not all of its elements are confusionMatrix objects.")
      }
      as.data.frame(x$table)
    }))
    positive <- cmList[[1]]$positive
    negative <- cmList[[1]]$negative

    # Add the commonly used names for the statistics in the confusion matrix
    dtNames <- data.table(
      Prediction = c(positive, negative, positive, negative),
      Reference  = c(positive, positive, negative, negative),
      Name       = c("TP",     "FN",     "FP",     "TN"    ))
    dt <- dplyr::left_join(dtNames, dt, by = c("Prediction", "Reference"))

    if(length(cmList) > 1) {
      # more than one confusionMatrix object
      matrixSummary <- dt[, list(Median = median(Freq), Q05 = quantile(Freq, probs = 0.05), Q95 = quantile(Freq, probs = 0.95)), by = list(Prediction, Reference, Name)]
      matrixSummary[, SummaryType := " Median [Q05,Q95]"]

      matrixSummary[, Freq := sapply(.I, function(i) sprintf("%.5g [%.5g, %.5g]", Median[i], Q05[i], Q95[i]))]
    } else {
      matrixSummary <- dt[, list(Prediction, Reference, Name, Freq)]
      matrixSummary[, SummaryType := ""]
    }

    matrixSummary[, Text:=paste0(Name, SummaryType, ": ", Freq)]
    matrixSummary <- matrixSummary %>% tidyr::pivot_wider(id_cols = c(Prediction), names_from = c("Reference"), values_from = "Text")
    names(matrixSummary) <- c("Predicted", paste0("Observed ", positive), paste0("Observed ", negative))

    # . Confusion matrix statistics summary ----

    dtStats <- data.table::rbindlist(lapply(cmList, function(cm) {
      dt <- data.table(
        Statistic = c("Prevalence", "Sensitivity (Recall)", "Specificity", "Precision", "Accuracy", "Balanced Accuracy", "F1 Score"),
        Value = c(cm$byClass['Prevalence'], cm$byClass['Sensitivity'], cm$byClass['Specificity'], cm$byClass['Precision'],
                  cm$overall['Accuracy'], cm$byClass['Balanced Accuracy'], cm$byClass['F1']),
        Interpretation = c(
          "Fraction observed POSITIVECLASS of all the subjects",
          "Fraction correctly predicted observed POSITIVECLASS",
          "Fraction correctly predicted observed NEGATIVECLASS",
          "Fraction observed POSITIVECLASS of the subjects predicted as POSITIVECLASS",
          "Fraction true predictions",
          "Mean of sensitivity and specificity",
          "Harmonic mean between precision and recall"),
        Formula = c(
          "(TP+FN)/(TP+FN+FP+TN)",
          "TP/(TP+FN)",
          "TN/(FP+TN)",
          "TP/TP+FP)",
          "(TP+TN)/(TP+FN+FP+TN)",
          "0.5*(Sensitivity+Specificity)",
          "2*Precision*Recall/(Precision+Recall)"))
      dt$Interpretation <- gsub("NEGATIVECLASS", cm$negative, dt$Interpretation, fixed = TRUE)
      dt$Interpretation <- gsub("POSITIVECLASS", cm$positive, dt$Interpretation, fixed = TRUE)
      dt
    }))

    if(length(cmList) > 1) {
      # more than one confusionMatrix object
      statsSummary <- dtStats[, list(Median = median(Value, na.rm = TRUE),
                                     Q05 = quantile(Value, probs = 0.05, na.rm = TRUE),
                                     Q95 = quantile(Value, probs = 0.95, na.rm = TRUE),
                                     NnonNA = sum(!is.na(Value))),
                              by = list(Statistic, Interpretation, Formula)]

      statsSummary[, `Median [Q05,Q95] (# not NA)` := sapply(.I, function(i) sprintf("%.2f [%.2f, %.2f] (%d)", Median[i], Q05[i], Q95[i], NnonNA[i]))]
      statsSummary <- statsSummary[, list(Statistic, `Median [Q05,Q95] (# not NA)`, Interpretation, Formula)]
    } else {
      statsSummary <- dtStats[, list(Statistic, Value = sprintf("%.2f", Value), Interpretation, Formula)]
    }

    list(matrixSummary = as.data.table(matrixSummary), statsSummary = statsSummary)
  }
}


