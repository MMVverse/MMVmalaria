#' display_ComboTreatments
#'
#' @description
#' @param data
#' @param stratify Default: 'STUDY'
#' @param fileplot Default: NULL
#' @param filetable Default: NULL
#' @param CompoundList Default: NULL
#' @param FLAGreturnObject Default: FALSE
#' @return
#' @export
#' @author Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom MMVbase MMVcolors
display_ComboTreatments <- function(data,
                                    stratify     = "STUDY",
                                    fileplot     = NULL,
                                    filetable    = NULL,
                                    CompoundList = NULL,
                                    FLAGreturnObject = FALSE

) {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Input checks ----

  # Check whether output is defined:
  if (is.null(fileplot) & is.null(filetable) & !FLAGreturnObject) stop("No output (either fileplot, filetable, or return graph object) is defined.\nPlease do so to get me doing something.")

  # Check whether data contained combination study data:
  if (!all(c("DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2") %in% names(data))) stop("Data do not seem to be from combination study.\n(DOSELEVEL1, DOSELEVEL2, DOSEMULT1, and/or DOSEMULT2 missing)")


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Gather information ----

  # Update Compound Name if wanted:
  if (!is.null(CompoundList)){
    data <- swapName_MMVnameToName(data         = data,
                                   CompoundList = CompoundList)
  }

  # Get compound names:
  tmp <- strsplit(grep("+" ,data$COMPOUND, value = TRUE, fixed = TRUE)[1], split = "+", fixed = TRUE)[[1]]
  Compound1 <- tmp[1]
  Compound2 <- tmp[2]

  # Get dosing units:
  DOSEUNIT1 <- data$UNIT[data$NAME == paste0(Compound1," Dose")][1]
  DOSEUNIT2 <- data$UNIT[data$NAME == paste0(Compound2," Dose")][1]

  # Handle stratification:
  if (is.null(stratify)) {
    data$STRAT <- "dummy"
  } else {
    data$STRAT <- data[[stratify]]
  }
  stratLevels <- sort(unique(data$STRAT))
  nStratLevels <- length(stratLevels)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract treatment information and number of individuals ----

  # Get subject-treatment overview:
  subjInfo <- unique(as.data.frame(data)[,c("USUBJID","TRTNAME", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2", "STRAT")])

  # Set all number of doses for vehicles to 0:
  idxVehicle <- (subjInfo$DOSELEVEL1+subjInfo$DOSELEVEL2) <= 2.0e-12
  subjInfo <- within(subjInfo, {
    TRTNAME[idxVehicle] <- "Vehicle"
    DOSEMULT1[idxVehicle] <- 0
    DOSEMULT2[idxVehicle] <- 0
  })

  # Get count of individuals per unique treatment:
  trtInfo <- ddply(subjInfo, ~TRTNAME+DOSELEVEL1+DOSEMULT1+DOSELEVEL2+DOSEMULT2+STRAT, function(xxx) {
    out <- data.frame(
      N = length(xxx$TRTNAME)
    )
    out
  })

  # Combine DOSELEVEL and DOSEMULT for classifying treatment of one drug:
  trtInfo <- within(trtInfo, {
    TRTNAME1 <- paste0(ifelse(DOSEMULT1==1, "", paste0(DOSEMULT1, "x")), DOSELEVEL1, " ", DOSEUNIT1)
    TRTNAME2 <- paste0(ifelse(DOSEMULT2==1, "", paste0(DOSEMULT2, "x")), DOSELEVEL2, " ", DOSEUNIT2)
    TRTNAME1[DOSELEVEL1<=1e-12] <- paste0("0 ", DOSEUNIT1)
    TRTNAME2[DOSELEVEL2<=1e-12] <- paste0("0 ", DOSEUNIT2)
  })

  # Order treatments per drug
  trtInfo1 <- unique(trtInfo[order(trtInfo$DOSEMULT1,trtInfo$DOSELEVEL1),c("TRTNAME1", "DOSELEVEL1", "DOSEMULT1")])
  trtLevels1 <- trtInfo1$TRTNAME1
  nLevels1 <- length(trtLevels1)
  trtInfo$TRTNAME1 <- factor(trtInfo$TRTNAME1, levels = trtLevels1)
  trtInfo2 <- unique(trtInfo[order(trtInfo$DOSEMULT2,trtInfo$DOSELEVEL2),c("TRTNAME2", "DOSELEVEL2", "DOSEMULT2")])
  trtLevels2 <- trtInfo2$TRTNAME2
  nLevels2 <- length(trtLevels2)
  trtInfo$TRTNAME2 <- factor(trtInfo$TRTNAME2, levels = trtLevels2)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot ----

  # Spread numbers from statification:
  dodgeD <- 0.3
  if (nStratLevels == 1) stratD <- 0 else stratD <- seq(-dodgeD, dodgeD, length.out = nStratLevels)
  trtInfoPlot <- within(trtInfo, {
    trt1pos <- as.numeric(TRTNAME1)
    trt2pos <- as.numeric(TRTNAME2)
    stratNo <- as.numeric(factor(STRAT, levels = stratLevels))
    strat1D <- stratD[stratNo]
  })

  # Do plot:
  gr <- IQRggplot(trtInfoPlot, aes(trt1pos+strat1D, TRTNAME2)) +
    geom_hline(yintercept=seq(1.5,nLevels2), color = "grey", size = 0.5, linetype = 2)+
    geom_vline(xintercept=seq(1.5,nLevels1), color = "grey", size = 0.5, linetype = 2)+
    geom_label(aes(label=N, color = STRAT), label.padding = unit(0.1, "lines"), label.r = unit(0.05, "lines")) +
    scale_y_discrete(limits = rev(trtLevels2)) +
    scale_x_continuous(breaks = 1:nLevels1, labels = trtLevels1, position = "top") +
    labs(x = Compound1,
         y = Compound2) +
    theme(panel.grid      = element_blank(),
          legend.position = "bottom",
          axis.text.x     = element_text(angle=0, hjust=0.5),
          axis.ticks      = element_blank()) +
    coord_cartesian(xlim=c(0.5,nLevels1+0.5), ylim=c(0.5,nLevels2+0.5), expand = FALSE)
  if (is.null(stratify)) {
    gr <- gr + scale_color_manual(values=MMVcolors, guide=FALSE)
  } else {
    gr <- gr + scale_color_manual(stratify, values=MMVcolors) +
      guides(color=guide_legend(override.aes = list(size = 5, label = "\u25FB")))
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Table ----

  # Function to concatenate number over stratification variable:
  concatenateNs <- function(x, stratLevels) {
    idxCol <- which(names(x) %in% stratLevels)
    nstr <- x[,idxCol[1]]
    for (k in idxCol[-1]) nstr <- paste(nstr, x[, k], sep="/")
    out <- x[,-idxCol]
    out$N <- nstr
    out
  }

  # Create contingency table:
  trtTable <- select(trtInfo,TRTNAME1, TRTNAME2, STRAT, N) # column selection
  trtTable <- spread(trtTable, STRAT,N, fill = "-") # spread by stratifier
  trtTable <- select(trtTable, TRTNAME1, TRTNAME2, stratLevels) # make sure that order is correct
  trtTable <- concatenateNs(trtTable, stratLevels)
  trtTable <- spread(trtTable, TRTNAME2, N)
  trtTable[is.na(trtTable)] <- ""
  names(trtTable) <- sub("TRTNAME1", paste0(Compound1, " (row)/", Compound2, " (column)"), names(trtTable))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output ----

  if (!is.null(fileplot)) {
    # Print to PNG:
    IQRoutputPNG(gr, fileplot, width = 3+0.5*nLevels1+0.5*nStratLevels, height = 2+0.3*nLevels2)
  }

  if (!is.null(filetable)) {
    filetable <- ifelse(grepl(".txt$",filetable, fixed=TRUE),filetable, paste0(filetable, ".txt"))
    # Print to txt file:
    if (nStratLevels == 1) stratInfo <- NULL else stratInfo <- paste0(stratLevels, collapse = "/")
    IQRoutputTable(trtTable, filename = filetable, xfooter = paste0("Order of stratified numbers: ", stratInfo), report = TRUE)
  }

  if (FLAGreturnObject) {
    return(list(Table = trtTable, Plot = gr, StratifyOrder = stratLevels))
  }
}
#' plot_PKdataMMV
#'
#' @description
#' @param dataGen
#' @param filePath
#' @param PKname Default: NULL
#' @param scale_y Default: 'log'
#' @param facet_scales Default: 'fixed'
#' @return
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join
plot_PKdataMMV <- function(dataGen,
                           filePath,
                           PKname       = NULL,
                           scale_y      = "log",
                           facet_scales = "fixed") {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle input arguments ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (length(unique(dataGen$NAME[dataGen$TYPENAME == "Dose"])) > 1)
    stop("Function only handles monotherapy data. Two types of doses found.")

  if (is.null(PKname))
    PKname <- unique(dataGen$NAME[dataGen$TYPENAME == "PK"])

  if (length(PKname) == 0)
    stop("No PK data found.")


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load dataset ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dataPlot       <- as.data.frame(dataGen)
  dataPlot$Label <- with(dataPlot, paste0(NAME, " (", UNIT,")"))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get nicer treatment name ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TRTann <- within(unique(dataPlot[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME
    Treatment <- ifelse(DOSELEVEL<=1e-12, "Vehicle", Treatment)
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]
  if ("Vehicle" %in% TRTann$Treatment){
    TRTann$Treatment <- factor(TRTann$Treatment, levels = c("Vehicle",setdiff(TRTann$Treatment, "Vehicle")))
  } else{
    TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment)
  }

  dataPlot <- plyr::join(dataPlot, TRTann[,c("TRTNAME","Treatment")])

  doseunit <- dataPlot[dataPlot$NAME %in% attr(dataGen, "doseNAMES"),"UNIT"][1]


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PK ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  k = 1
  for (PKnamek in PKname){
    # Reduce dataset to observation of interest:
    dataPlotPK <- dataPlot[dataPlot$NAME %in% c(PKnamek) & dataPlot$MDV==0 & dataPlot$Treatment!="Vehicle", ]
    UnitPK     <- dataPlotPK$UNIT[1]

    # Get LLOQ information:
    lloq <- unique(dataPlotPK[,c("Label","LLOQ","STUDY")])

    # Plot with PK and PD data:
    gr <- IQRggplot(dataPlotPK, aes(TIME, VALUE)) +
      geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
      geom_line(aes(group = USUBJID), color = "grey") +
      scale_color_IQRtools() +
      facet_wrap(~STUDY, scales = facet_scales) +
      labs(x = "Time [h]",
           y = paste0(PKnamek," [",UnitPK,"]"))

    # To avoid to have a legend if only one treatment:
    if (length(unique(dataPlotPK$Treatment))>1){
      gr <- gr + geom_point(aes(color = Treatment))
    } else{
      gr <- gr + geom_point(color = "dodgerblue3")
    }

    # log10 scale?
    if (scale_y=="log"){
      gr <- gr + scale_y_log10()
    }

    # Save Plot:
    if (length(PKname)==1){
      IQRoutputPNG(gr, filename = file.path(filePath, "01_PKlinePlot.png"))
    } else{
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("01_PKlinePlot_DRUG_",k,".png")))
      k = k+1
    }
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PK per treatment ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  k = 1
  for (PKnamek in PKname){
    # Reduce dataset to observation of interest:
    dataPlotPK <- dataPlot[dataPlot$NAME %in% PKnamek & dataPlot$Treatment != "Vehicle", ]
    UnitPK     <- dataPlotPK$UNIT[1]

    # Get LLOQ information:
    lloq <- unique(dataPlotPK[,c("Label","LLOQ","STUDY")])

    # Plot with PK and PD data:
    gr <- IQRggplot(dataPlotPK, aes(TIME, VALUE)) +
      geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
      geom_line(aes(group = USUBJID), color = "grey") +
      scale_color_IQRtools() +
      facet_wrap(~Treatment, scales = facet_scales) +
      labs(x = "Time [h]",
           y = paste0(PKnamek," [",UnitPK,"]"))

    # To avoid to have a legend if only one study:
    if (length(unique(dataPlotPK$STUDY))>1){
      gr <- gr + geom_point(aes(color = STUDY))
    } else{
      gr <- gr + geom_point(color = "dodgerblue3")
    }

    # log10 scale?
    if (scale_y=="log"){
      gr <- gr + scale_y_log10()
    }

    # Save Plot:
    if (length(PKname)==1){
      IQRoutputPNG(gr, filename = file.path(filePath, "02_PKlinePlot_ByTRT.png"))
    } else{
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("02_PKlinePlot_ByTRT_DRUG_",k,".png")))
      k = k+1
    }
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot dose normalized PK ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Unit Dose:
  dataDose = dataGen[dataGen$TYPENAME == "Dose",]
  UnitDose = dataDose$UNIT[1]

  k = 1
  for (PKnamek in PKname){
    # Reduce dataset to observation of interest:
    dataPlotPK           <- dataPlot[dataPlot$NAME %in% PKnamek & dataPlot$Treatment != "Vehicle", ]
    dataPlotPK$VALUEnorm <- dataPlotPK$VALUE / dataPlotPK$DOSELEVEL
    UnitPK               <- dataPlotPK$UNIT[1]

    if (dim(dataPlotPK)[1] > 0) {
      gr <- IQRggplot(dataPlotPK, aes(TIME, VALUEnorm)) +
        geom_line(aes(group = USUBJID), color = "grey") +
        scale_color_IQRtools() +
        facet_wrap(~STUDY, scales = facet_scales) +
        labs(x = "Time [h]",
             y = paste0("Dose-Normalized Concentrations (",dataPlotPK$UNIT[1],"/",doseunit,")"))

      # To avoid to have a legend if only one treatment:
      if (length(unique(dataPlotPK$Treatment))>1){
        gr <- gr + geom_point(aes(color = Treatment))
      } else{
        gr <- gr + geom_point(color = "dodgerblue3")
      }

      # log10 scale?
      if (scale_y=="log"){
        gr <- gr + scale_y_log10()
      }

      # Save Plot:
      if (length(PKname)==1){
        IQRoutputPNG(gr, filename = file.path(filePath, "03_PKnormLinePlot.png"))
      } else{
        IQRoutputPNG(gr, filename = file.path(filePath, paste0("03_PKnormLinePlot_DRUG_",k,".png")))
        k = k+1
      }
    }
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot dose normalized PK  per treatment ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Unit Dose:
  dataDose = dataGen[dataGen$TYPENAME == "Dose",]
  UnitDose = dataDose$UNIT[1]

  k = 1
  for (PKnamek in PKname){
    # Reduce dataset to observation of interest:
    dataPlotPK           <- dataPlot[dataPlot$NAME %in% PKnamek & dataPlot$Treatment != "Vehicle", ]
    dataPlotPK$VALUEnorm <- dataPlotPK$VALUE / dataPlotPK$DOSELEVEL
    dataPlotPK$LLOQnorm <- dataPlotPK$LLOQ / dataPlotPK$DOSELEVEL
    UnitPK               <- dataPlotPK$UNIT[1]

    # Get LLOQ information:
    lloq <- unique(dataPlotPK[,c("Label","LLOQnorm","STUDY","Treatment")])

    if (dim(dataPlotPK)[1] > 0) {
      gr <- IQRggplot(dataPlotPK, aes(TIME, VALUEnorm)) +
        geom_hline(data=lloq, aes(yintercept = LLOQnorm), linetype = 2, color = "black") +
        geom_line(aes(group = USUBJID), color = "grey") +
        scale_color_IQRtools() +
        facet_wrap(~Treatment, scales = facet_scales) +
        labs(x = "Time [h]",
             y = paste0("Dose-Normalized Concentrations (",dataPlotPK$UNIT[1],"/",doseunit,")"))

      # To avoid to have a legend if only one study:
      if (length(unique(dataPlotPK$STUDY))>1){
        gr <- gr + geom_point(aes(color = STUDY))
      } else{
        gr <- gr + geom_point(color = "dodgerblue3")
      }

      # log10 scale?
      if (scale_y=="log"){
        gr <- gr + scale_y_log10()
      }

      # Save Plot:
      if (length(PKname)==1){
        IQRoutputPNG(gr, filename = file.path(filePath, "04_PKnormLinePlot_ByTRT.png"))
      } else{
        IQRoutputPNG(gr, filename = file.path(filePath, paste0("04_PKnormLinePlot_ByTRT_DRUG_",k,".png")))
        k = k+1
      }
    }
  }
}
#' plot_PKPDdataMMV
#'
#' @description
#' @param dataGen
#' @param filePath
#' @param PDname Default: 'Parasitemia'
#' @param PKname Default: NULL
#' @param HuErys Default: NULL
#' @param doseUnit Default: NULL
#' @param ActivityPath Default: NULL
#' @return
#' @export
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join
plot_PKPDdataMMV <- function(dataGen,
                             filePath,
                             PDname   = "Parasitemia",
                             PKname   = NULL,
                             HuErys   = NULL,
                             doseUnit = NULL,
                             ActivityPath = NULL) {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle input arguments ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (length(unique(dataGen$NAME[dataGen$TYPENAME=="Dose"])) > 1){
    stop("Function only handles monotherapy data. Two types of doses found.")
  }
  if (is.null(PKname)){
    PKname <- unique(dataGen$NAME[dataGen$TYPENAME=="PK"])
  }

  if (length(PKname)>1){
    stop("Two types of PK data found.")
  }

  if (length(PKname)==0){
    stop("No PK data found.")
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load dataset ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dataPlot       <- as.data.frame(dataGen)
  dataPlot$Label <- with(dataPlot, paste0(NAME, " [", UNIT,"]"))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get nicer treatment name ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TRTann <- within(unique(dataPlot[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME
    Treatment <- ifelse(DOSELEVEL <= 1e-12, "Vehicle", Treatment)
    Treatment <- ifelse(is.na(DOSELEVEL), "Unknown", Treatment)
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]
  # Get nice order for treatment
  # Note that if multiple 0 dose treatments exist, there might be multiple groups with "Vehicle",
  # Therefore set level order carefully regarding the vehicle treatment
  TRTann$Treatment <- factor(TRTann$Treatment, levels = c("Vehicle",setdiff(TRTann$Treatment, "Vehicle")))

  dataPlot <- plyr::join(dataPlot, TRTann[,c("TRTNAME","Treatment")])

  if (is.null(doseUnit)) {
    doseunit <- dataPlot[dataPlot$NAME %in% attr(dataGen, "doseNAMES"),"UNIT"][1]
  } else {
    doseunit <- doseUnit
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check that STUDY is not numeric ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (!class(dataPlot$STUDY) %in% c("factor", "character")) {
    dataPlot$STUDY <- as.character(dataPlot$STUDY)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PK and PD observations ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Reduce dataset to observation of interest:
  dataPlotPKPD <- dataPlot[dataPlot$NAME %in% c(PDname, PKname) & dataPlot$MDV == 0, ]

  # Get LLOQ information:
  lloq <- unique(dataPlotPKPD[,c("Label","LLOQ","STUDY")])

  # Plot with PK and PD data:
  gr <- MMVggplot(dataPlotPKPD, aes(TIME, VALUE), ActivityPath = ActivityPath) +
    geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
    geom_line(aes(group = USUBJID), color = "grey") +
    geom_point(aes(color = Treatment)) +
    scale_color_IQRtools() +
    scale_y_log10() +
    facet_grid(Label~STUDY, scales = "free") +
    labs(x = "Time [hr]",
         y = "")

  # Save Plot:
  IQRoutputPNG(gr, filename = file.path(filePath, "01-PKPDlinePlot.png"))

  # Plot with PK and PD data:
  gr <- MMVggplot(dataPlotPKPD, aes(TIME, VALUE), ActivityPath = ActivityPath) +
    geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
    geom_line(aes(group = USUBJID), color = "grey") +
    geom_point(aes(color = Treatment, shape = STUDY)) +
    scale_color_IQRtools() +
    scale_y_log10() +
    facet_grid(Label~., scales = "free") +
    labs(x = "Time [hr]",
         y = "")

  # Save Plot:
  IQRoutputPNG(gr, filename = file.path(filePath, "01-PKPDlinePlot_Unstratified.png"))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PK ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (k in 1:length(PKname)){
    # PK name:
    PKname_k <- PKname[k]

    # Reduce dataset to observation of interest:
    dataPlotPK      <- dataPlot[dataPlot$NAME %in% PKname_k & dataPlot$Treatment != "Vehicle", ]
    dataPlotPK.LLOQ <- dataPlot[dataPlot$NAME %in% PKname_k & dataPlot$Treatment != "Vehicle" & dataPlot$CENS == 1, ]

    if (dim(dataPlotPK)[1] > 0) {
      # Get LLOQ information:
      lloq <- unique(dataPlotPK[,c("Label","LLOQ","STUDY")])


      #~~~~~~~
      # Plot PK Data by STUDY:
      gr <- MMVggplot(dataPlotPK, aes(TIME, VALUE), ActivityPath = ActivityPath) +
        geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = Treatment)) +
        geom_point(data = dataPlotPK.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        facet_wrap(~STUDY, scales = "free") +
        labs(x = "Time [hr]",
             y = dataPlotPK$Label[1])
      # Save Plot:
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("02-", letters[k], "-PKlinePlot_BySTUDY.png")))


      #~~~~~~~
      # Plot PK Data by TRT:
      gr <- MMVggplot(dataPlotPK, aes(TIME, VALUE), ActivityPath = ActivityPath) +
        geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = STUDY)) +
        geom_point(data = dataPlotPK.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        facet_wrap(~Treatment, scales = "free") +
        labs(x = "Time [hr]",
             y = dataPlotPK$Label[1])
      # Save Plot:
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("02-", letters[k], "-PKlinePlot_ByTRT.png")))
    }
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot dose normalized PK ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (k in 1:length(PKname)){
    # PK name:
    PKname_k <- PKname[k]

    # Reduce dataset to observation of interest:
    dataPlotPK <- dataPlot[dataPlot$NAME %in% PKname_k & dataPlot$Treatment != "Vehicle", ]
    dataPlotPK$VALUEnorm <- dataPlotPK$VALUE / dataPlotPK$DOSELEVEL
    dataPlotPK.LLOQ <- dataPlot[dataPlot$NAME %in% PKname_k & dataPlot$Treatment != "Vehicle" & dataPlot$CENS == 1, ]
    dataPlotPK.LLOQ$VALUEnorm <- dataPlotPK.LLOQ$VALUE / dataPlotPK.LLOQ$DOSELEVEL

    if (dim(dataPlotPK)[1] > 0) {

      #~~~~~~~
      # Plot PK Data Unstratified:
      gr <- MMVggplot(dataPlotPK, aes(TIME, VALUEnorm), ActivityPath = ActivityPath) +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = Treatment, shape = STUDY)) +
        geom_point(data = dataPlotPK.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        labs(x = "Time [hr]",
             y = paste0("Dose-Normalized Concentrations [",dataPlotPK$UNIT[1],"/",doseunit,"]"))
      # Save Plot:
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("03-", letters[k], "-PKnormlizedPlot_Unstratified.png")))

      #~~~~~~~
      # Plot PK Data by STUDY:
      gr <- MMVggplot(dataPlotPK, aes(TIME, VALUEnorm), ActivityPath = ActivityPath) +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = Treatment)) +
        geom_point(data = dataPlotPK.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        facet_wrap(~STUDY, scales = "free") +
        labs(x = "Time [hr]",
             y = paste0("Dose-Normalized Concentrations [",dataPlotPK$UNIT[1],"/",doseunit,"]"))
      # Save Plot:
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("03-", letters[k], "-PKnormlizedPlot_BySTUDY.png")))


      #~~~~~~~
      # Plot PK Data by TRT:
      gr <- MMVggplot(dataPlotPK, aes(TIME, VALUEnorm), ActivityPath = ActivityPath) +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = STUDY)) +
        geom_point(data = dataPlotPK.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        facet_wrap(~Treatment, scales = "free") +
        labs(x = "Time [hr]",
             y = paste0("Dose-Normalized Concentrations [",dataPlotPK$UNIT[1],"/",doseunit,"]"))
      # Save Plot:
      IQRoutputPNG(gr, filename = file.path(filePath, paste0("03-", letters[k], "-PKnormlizedPlot_ByTRT.png")))
    }
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PD ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  replaceSpecChars <- function(name) {
    name <- gsub("+", "_plus_", name, fixed = TRUE)
    name <- gsub(" ", "_", name, fixed = TRUE)
    name
  }

  for(PDname.i in PDname) {
    # Reduce dataset to observation of interest:
    dataPlotPD <- dataPlot[dataPlot$NAME == PDname.i, ]
    dataPlotPD.LLOQ <- dataPlot[dataPlot$NAME %in% PDname.i & dataPlot$Treatment != "Vehicle" & dataPlot$CENS == 1, ]

    if (dim(dataPlotPD)[1] > 0) {
      # Get LLOQ information:
      lloq <- unique(dataPlotPD[,c("Label","LLOQ","STUDY")])
      # Get PD label
      PDlabel <- dataPlotPD$Label[1]

      #~~~~~~~
      # Plot PD Data by STUDY:
      gr <- MMVggplot(dataPlotPD, aes(TIME, VALUE), ActivityPath = ActivityPath) +
        geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = Treatment)) +
        geom_point(data = dataPlotPD.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        facet_wrap(~STUDY, scales = "free") +
        labs(x = "Time [hr]",
             y = dataPlotPD$Label[1])
      # Save Plot:
      IQRoutputPNG(
        gr,
        filename = file.path(filePath, "04-PDlinePlots", replaceSpecChars(paste0(PDname.i, "_BySTUDY.png")))
        )

      #~~~~~~~
      # Plot PD Data by TRT:
      gr <- MMVggplot(dataPlotPD, aes(TIME, VALUE), ActivityPath = ActivityPath) +
        geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = STUDY)) +
        geom_point(data = dataPlotPD.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        facet_wrap(~Treatment, scales = "free") +
        labs(x = "Time [hr]",
             y = dataPlotPD$Label[1])
      # Save Plot:
      IQRoutputPNG(
        gr,
        filename = file.path(filePath, "04-PDlinePlots", replaceSpecChars(paste0(PDname.i, "_ByTRT.png")))
      )
    }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PD of untreated if exist ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Reduce dataset to observation of interest:
    dataPlotVehicle      <- dataPlot[dataPlot$NAME %in% PDname.i & dataPlot$Treatment == "Vehicle", ]
    dataPlotVehicle.LLOQ <- dataPlot[dataPlot$NAME %in% PDname.i & dataPlot$Treatment == "Vehicle" & dataPlot$CENS == 1, ]

    #if ("Vehicle" %in% dataPlot$Treatment) {
    #  dataPlotVehicle <- dataPlot[dataPlot$Treatment == "Vehicle" & dataPlot$NAME == PDname.i, ]
    if (dim(dataPlotVehicle)[1] > 0) {

      # Get LLOQ information
      lloq <- unique(dataPlotVehicle[,c("Label","LLOQ","STUDY")])

      # Plot with PK and PD data
      gr <- MMVggplot(dataPlotVehicle, aes(TIME, VALUE), ActivityPath = ActivityPath) +
        geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
        geom_line(aes(group = USUBJID), color = "grey") +
        geom_point(aes(color = STUDY)) +
        geom_point(data = dataPlotVehicle.LLOQ, shape = 4, color = "red") +
        scale_color_IQRtools() +
        scale_y_log10() +
        labs(x = "Time [hr]",
             y = dataPlotVehicle$Label[1])

      # Save Plot:
      IQRoutputPNG(
        gr,
        filename = file.path(filePath, "05-UntreatedPDlinePlots", replaceSpecChars(paste0(PDname.i, ".png")))
      )
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot with erythrocytes if given ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (!is.null(HuErys)) {
    dataPlotHuEry <- dataPlot[dataPlot$NAME %in% HuErys, ]

    gr <- MMVggplot(dataPlotHuEry, aes(TIME, VALUE), ActivityPath = ActivityPath) +
      geom_line(aes(group = USUBJID), color = "grey") +
      geom_point(aes(color = Treatment, shape = STUDY)) +
      scale_color_IQRtools() +
      # scale_y_log10() +
      # coord_cartesian(ylim = c(1,140)) +
      # facet_grid(.~STUDY, scales = "free") +
      labs(x = "Time [hr]",
           y = dataPlotHuEry$Label[1])

    # Adjust shape id to many studies:
    if (length(unique(dataPlotHuEry$STUDY)) > 6){
      gr <- gr + scale_shape_manual(values = LETTERS[1:26])
    }

    # Save Plot:
    IQRoutputPNG(gr, filename = file.path(filePath, "06-HuErythrocyteslinePlot.png"))
  }

}
#' plot_PKPDhuChCombo
#'
#' @description
#' @param dataGen
#' @param filePath
#' @param PDname Default: 'Parasitemia'
#' @param Gametos Default: 'Parasitemia Gametocytes'
#' @param PKname Default: NULL
#' @param TRTNAME Column name of treatment group to use
#' @param ByUSUBJID If the plot should be split by subject (Default: TRUE)
#' @param ActivityPath Path of the current activity (Default: NULL)
#' @return
#' @export
#' @author Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join dlply
plot_PKPDhuChCombo <- function(dataGen,
                               filePath,
                               PDname  = "Parasitemia",
                               Gametos = "Parasitemia Gametocytes",
                               PKname  = NULL,
                               TRTNAME      = "TRTNAME",
                               ByUSUBJID    = TRUE,
                               ActivityPath = NULL) {


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # STEP 1: Load dataset ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Transform to data.frame:
  dataPlot       <- as.data.frame(dataGen)

  # Adjust TRTNAME column:
  if(TRTNAME!="TRTNAME"){
    dataPlot$TRTNAME <- NULL
    names(dataPlot)[names(dataPlot)==TRTNAME] <- "TRTNAME"
  }

  # Add Unit Label:
  dataPlot$Label <- with(dataPlot, paste0(NAME, " [", UNIT,"]"))

  # Get Compound Name:
  Compounds <- with(dataGen, COMPOUND[DOSELEVEL1 > 0.001 & DOSELEVEL2 > 0.001][1])
  Compound1 <- gsub("[+][[:alnum:]]*$","",Compounds)
  Compound2 <- gsub("^[[:alnum:]]*[+]","",Compounds)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # STEP 2: Handle input arguments ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Adjust PK name:
  if (is.null(PKname)){
    PKname <- paste0(c(Compound1, Compound2), " Plasma Concentration")
  }

  # Make sure that at least one PK name is defined:
  if (length(PKname) == 0){
    stop("No PK data found.")
  }

  # Get Dose and Resuce Name:
  drugDose   <- unique(dataGen$NAME[dataGen$TYPENAME=="Dose"])
  rescueDose <- grep("Riamet",doseNAMES(dataGen), value = TRUE)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # STEP 3: Get nicer treatment name ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Change Names to look better:
  TRTann <- within(unique(dataPlot[,c("TRTNAME", "STUDY", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")]), {
    Treatment <- TRTNAME
    Treatment <- gsub("[+]", " + ", Treatment)
    Treatment <- gsub("mg", "mg " , Treatment)
	Treatment <- gsub("mg/kg", "mg/kg " , Treatment, fixed = TRUE)
	Treatment <- gsub("mg /kg", "mg/kg ", Treatment, fixed = TRUE)
    idx <- 0
    idx[TRTNAME == "Vehicle"] <- 1
    idx[DOSEMULT2 == 0 & DOSEMULT1 > 0] <- 2
    idx[DOSEMULT2 > 0 & DOSEMULT1 == 0] <- 3
    idx[DOSEMULT2 > 0 & DOSEMULT1 > 0] <- 4
  })

  # Order nicely:
  TRTann           <- TRTann[with(TRTann, order(idx, DOSEMULT1, DOSEMULT2, DOSELEVEL1, DOSELEVEL2)),]
  TRTann$idx       <- NULL
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])
  dataPlot         <- plyr::join(dataPlot, TRTann[,c("TRTNAME","Treatment")])
  doseunit         <- dataPlot[dataPlot$NAME %in% attr(dataGen, "doseNAMES"),"UNIT"][1]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # STEP 4: Plot PK ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Patients with PD data of interest
  PatientsOfInterestPK <- unique(dataPlot[dataPlot$NAME %in% PKname,"USUBJID"])

  # Reduce dataset to observation of interest:
  dataPlotPK <- dataPlot[dataPlot$NAME %in% c(PKname, drugDose, rescueDose), ]
  dataPlotPK <- dataPlotPK[dataPlotPK$USUBJID %in% PatientsOfInterestPK, ]

  # Get LLOQ information:
  lloq <- unique(subset(dataPlotPK, TYPENAME=="PK", c("Label","LLOQ","STUDY")))

  # Get Breaks:
  dx      <- 24*7   # Correspond to a week
  x_break <- seq(floor(min(dataPlotPK$TIME/dx))*dx, ceiling(max(dataPlotPK$TIME/dx))*dx, dx)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4.a. All PK ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  # Plots: dlply instead of a for loop
  grList <- plyr::dlply(dataPlotPK, ~Treatment, function(x) {

    # Datasubet for plot:
    xPK  <- x[x$TYPENAME=="PK",]
    xPKc <- x[x$TYPENAME=="PK" & x$MDV == 1,]
    xD   <- x[x$TYPENAME=="Dose",]
    xR   <- x[x$TYPENAME=="Rescue",]

    # Caption:
    Caption <- "Experimental drug and recue medication doses indicated by black and purple lines respectively."

    # Plot:
    gr_k <- MMVggplot(xPK, aes(TIME, VALUE), ActivityPath = ActivityPath, Caption = Caption) +
      geom_vline(data=xD, aes(xintercept = TIME), linetype = 2) +
      geom_vline(data=xR, aes(xintercept = TIME), linetype = 4, color = "purple") +
      geom_line(aes(group = interaction(USUBJID,NAME)), color = "grey") +
      geom_point(aes(color = NAME)) +
      geom_point(data=xPKc, shape = 4) +
      scale_y_log10() +
      scale_x_continuous(breaks = x_break) +
      labs(x = "Time [hr]",
           y = "Concentration [ug/mL]") +
      theme(legend.position = "bottom")

    # Adjust Colors:
    base::suppressMessages(gr_k <- gr_k +
                             scale_color_manual("Compound", values = MMVcolors[2:length(MMVcolors)]))

    # By Subject:
    if(ByUSUBJID){
      gr_k <- gr_k +
        facet_wrap(~SUBJECT) +
        labs(title = paste0("PK Plot for ", x$Treatment[1], " by Subject"))
    }else{
      gr_k <- gr_k +
        labs(title = paste0("PK Spaghetti Plot for ", x$Treatment[1]))
    }

    # Output:
    return(gr_k)
  })

  # Save Plots:
  IQRoutputPDF(grList, filename = file.path(filePath, "07-PKplot_All.pdf"))


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4.b. PK by PKname ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (k in seq_along(PKname)){
    # Get PK name:
    PKname_k <- PKname[[k]]

    # Get Compound Name:
    CPD_k <-  gsub("Concentration", "", PKname_k)
    # CPD_k <-  gsub("Blood"        , "", CPD_k)
    # CPD_k <-  gsub("Plasma"       , "", CPD_k)
    # CPD_k <-  gsub("DBS"          , "", CPD_k)
    # CPD_k <-  gsub("DPS"          , "", CPD_k)
    CPD_k <-  gsub("Blood"        , "blood" , CPD_k)
    CPD_k <-  gsub("Plasma"       , "plasma", CPD_k)
    CPD_k <-  gsub("DBS"          , "dbs"   , CPD_k)
    CPD_k <-  gsub("DPS"          , "dps"   , CPD_k)
    CPD_k <-  gsub(" "            , ""      , CPD_k)

    # Patients with PK data of interest
    PatientsOfInterestPK_k <- unique(dataPlotPK[dataPlotPK$NAME %in% PKname_k,"USUBJID"])

    # Reduce dataset to observation of interest:
    dataPlotPK_k <- dataPlotPK[dataPlotPK$NAME %in% c(PKname_k, drugDose, rescueDose), ]
    dataPlotPK_k <- dataPlotPK_k[dataPlotPK_k$USUBJID %in% PatientsOfInterestPK_k, ]

    # Get LLOQ information:
    lloq <- unique(subset(dataPlotPK_k, TYPENAME=="PK", c("Label","LLOQ","STUDY")))

    # Plots: dlply instead of a for loop
    grList_k <- plyr::dlply(dataPlotPK_k, ~Treatment, function(x) {

      # Datasubet for plot:
      xPK  <- x[x$TYPENAME=="PK",]
      xPKc <- x[x$TYPENAME=="PK" & x$MDV==1,]
      xD   <- x[x$TYPENAME=="Dose",]
      xR   <- x[x$TYPENAME=="Rescue",]

      # Caption:
      Caption <- "Experimental drug and recue medication doses indicated by black and purple lines respectively."

      # Plot:
      gr_i <- MMVggplot(xPK, aes(TIME, VALUE), ActivityPath = ActivityPath, Caption = Caption) +
        geom_vline(data=xD, aes(xintercept = TIME), linetype = 2) +
        geom_vline(data=xR, aes(xintercept = TIME), linetype = 4, color = "purple") +
        geom_line(aes(group = interaction(USUBJID,NAME)), color = "grey") +
        #geom_point(aes(color = NAME)) +
        geom_point(color = MMVcolors[2]) +
        geom_point(data=xPKc, shape = 4) +
        scale_y_log10() +
        scale_x_continuous(breaks = x_break) +
        labs(x = "Time [hr]",
             y = xPK$Label[1]) +
        theme(legend.position = "bottom")

      # Adjust Colors:
      base::suppressMessages(gr_i <- gr_i +
                               scale_color_manual("Compound", values = MMVcolors[2:length(MMVcolors)]))

      # By Subject:
      if(ByUSUBJID){
        gr_i <- gr_i +
          facet_wrap(~SUBJECT) +
          labs(title = paste0("PK Plot for ", x$Treatment[1], " by Subject"))
      }else{
        gr_i <- gr_i +
          labs(title = paste0("PK Spaghetti Plot for ", x$Treatment[1]))
      }

      # Output:
      return(gr_i)
    })

    # Save Plots:
    IQRoutputPDF(grList_k, filename = file.path(filePath, paste0("07-PKplot_", CPD_k, ".pdf")))
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # STEP 5: Plot parasites total and gametocytes if given ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Patients with PD data of interest
  PatientsOfInterestPD <- unique(dataPlot[dataPlot$NAME %in% PDname,"USUBJID"])

  # Reduce dataset to observation of interest:
  dataPlotPD <- dataPlot[dataPlot$NAME %in% c(PDname, Gametos, drugDose, rescueDose), ]
  dataPlotPD <- dataPlotPD[dataPlotPD$USUBJID %in% PatientsOfInterestPD, ]

  # Get LLOQ information:
  lloq <- unique(subset(dataPlotPD, TYPENAME=="Efficacy", c("Label","LLOQ","STUDY")))

  # Get Breaks:
  dx      <- 24*7   # Correspond to a week
  x_break <- seq(floor(min(dataPlotPD$TIME/dx))*dx, ceiling(max(dataPlotPD$TIME/dx))*dx, dx)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 5.a. All D ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Plots: dlply instead of a for loop
  grList <- plyr::dlply(dataPlotPD, ~Treatment, function(x) {

    # Datasubet for plot:
    xP  <- x[x$TYPENAME=="Efficacy",]
    xPc <- x[x$TYPENAME=="Efficacy" & x$MDV == 1,]
    xD  <- x[x$TYPENAME=="Dose",]
    xR  <- x[x$TYPENAME=="Rescue",]

    # Caption:
    Caption <- "Experimental drug and recue medication doses indicated by black and purple lines respectively."

    # Plot:
    gr_k <- MMVggplot(xP, aes(TIME, VALUE), ActivityPath = ActivityPath, Caption = Caption) +
      geom_vline(data=xD, aes(xintercept = TIME), linetype = 2) +
      geom_vline(data=xR, aes(xintercept = TIME), linetype = 4, color = "purple") +
      geom_line(aes(group = interaction(USUBJID,NAME)), color = "grey") +
      geom_point(aes(color = NAME)) +
      geom_point(data=xPc, shape = 4) +
      scale_y_log10() +
      scale_x_continuous(breaks = x_break) +
      labs(x = "Time [hr]",
           y = "Parasite Count [p/mL]") +
      theme(legend.position = "bottom")

    # Adjust Colors:
    base::suppressMessages(gr_k <- gr_k +
                             scale_color_manual("Parasite Type", values = MMVcolors[2:length(MMVcolors)]))

    # By Subject:
    if(ByUSUBJID){
      gr_k <- gr_k +
        facet_wrap(~SUBJECT) +
        labs(title = paste0("Parasitemia Plot for ", x$Treatment[1], " by Subject"))
    }else{
      gr_k <- gr_k +
        labs(title = paste0("Parasitemia Spaghetti Plot for ", x$Treatment[1]))
    }

    # Output:
    return(gr_k)
  })

  # Save Plots:
  IQRoutputPDF(grList, filename = file.path(filePath, "08-PDplot.pdf"))


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 5.b. PD by PDname ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (k in seq_along(PDname)){
    # Get PK name:
    PDname_k <- PDname[[k]]

    # Get Parasitemia Name:
    Para_k <-  gsub(" ", "", PDname_k)

    # Patients with PD data of interest
    PatientsOfInterestPD_k <- unique(dataPlotPD[dataPlotPD$NAME %in% PDname_k,"USUBJID"])

    # Reduce dataset to observation of interest:
    dataPlotPD_k <- dataPlotPD[dataPlotPD$NAME %in% c(PDname_k, drugDose, rescueDose), ]
    dataPlotPD_k <- dataPlotPD_k[dataPlotPD_k$USUBJID %in% PatientsOfInterestPD_k, ]

    # Get LLOQ information:
    lloq <- unique(subset(dataPlotPD_k, TYPENAME=="Efficacy", c("Label","LLOQ","STUDY")))

    # Plots: dlply instead of a for loop
    grList_k <- plyr::dlply(dataPlotPD_k, ~Treatment, function(x) {

      # Datasubet for plot:
      xPD  <- x[x$TYPENAME=="Efficacy",]
      xPDc <- x[x$TYPENAME=="Efficacy" & x$MDV==1,]
      xD   <- x[x$TYPENAME=="Dose",]
      xR   <- x[x$TYPENAME=="Rescue",]

      # Caption:
      Caption <- "Experimental drug and recue medication doses indicated by black and purple lines respectively."

      # Plot:
      gr_i <- MMVggplot(xPD, aes(TIME, VALUE), ActivityPath = ActivityPath, Caption = Caption) +
        geom_vline(data=xD, aes(xintercept = TIME), linetype = 2) +
        geom_vline(data=xR, aes(xintercept = TIME), linetype = 4, color = "purple") +
        geom_line(aes(group = interaction(USUBJID,NAME)), color = "grey") +
        #geom_point(aes(color = NAME)) +
        geom_point(color = MMVcolors[2]) +
        geom_point(data=xPDc, shape = 4) +
        scale_y_log10() +
        scale_x_continuous(breaks = x_break) +
        labs(x = "Time [hr]",
             y = xPD$Label[1]) +
        theme(legend.position = "bottom")

      # Adjust Colors:
      base::suppressMessages(gr_i <- gr_i +
                               scale_color_manual("Parasite Type", values = MMVcolors[2:length(MMVcolors)]))

      # By Subject:
      if(ByUSUBJID){
        gr_i <- gr_i +
          facet_wrap(~SUBJECT) +
          labs(title = paste0("Parasitemia Plot for ", x$Treatment[1], " by Subject"))
      }else{
        gr_i <- gr_i +
          labs(title = paste0("Parasitemia Spaghetti Plot for ", x$Treatment[1]))
      }

      # Output:
      return(gr_i)
    })

    # Save Plots:
    IQRoutputPDF(grList_k, filename = file.path(filePath, paste0("08-PDplot_", Para_k, ".pdf")))
  }
}
#' plot_PKPDscidCombo
#'
#' @description
#' @param dataGen
#' @param filePath
#' @param PDname Default: 'Parasitemia'
#' @param PKname Default: NULL
#' @param HuErys Default: NULL
#' @param CompoundList Default: NULL
#' @param ActivityPath Path of the current activity (Default: NULL)
#' @return
#' @export
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join dlply
plot_PKPDscidCombo <- function(dataGen,
                               filePath,
                               PDname = "Parasitemia",
                               PKname = NULL,
                               HuErys = NULL,
                               CompoundList = NULL,
                               ActivityPath = NULL) {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load dataset ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Update Compound Name if wanted:
  if (!is.null(CompoundList)){
    dataGen <- swapName_MMVnameToName(data         = dataGen,
                                      CompoundList = CompoundList)
  }

  # Generate Label:
  dataPlot       <- as.data.frame(dataGen)
  dataPlot$Label <- with(dataPlot, paste0(NAME, " [", UNIT,"]"))

  # Get Compound Names:
  Compounds <- with(dataGen, COMPOUND[DOSELEVEL1>1e-12 & DOSELEVEL2>1e-12][1])
  Compound1 <- gsub("[+][[:alnum:]]*$","",Compounds)
  Compound2 <- gsub("^[[:alnum:]]*[+]","",Compounds)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle input arguments ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Adjust PK name:
  if (is.null(PKname)){
    PKname <- paste0(c(Compound1, Compound2), " Blood Concentration")
  }

  # Make sure that at least one PK name is defined:
  if (length(PKname) == 0){
    stop("No PK data found.")
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get nicer treatment name ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Change Names to look better:
  TRTann <- within(unique(dataPlot[,c("TRTNAME", "STUDY", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")]), {
    Treatment <- TRTNAME
    Treatment <- gsub("[+]", " + ", Treatment)
    Treatment <- gsub("mg", "mg ", Treatment)
    Treatment <- gsub("mg /kg", "mg/kg ", Treatment, fixed = TRUE)
    idx <- 0
    idx[TRTNAME == "Vehicle"] <- 1
    idx[DOSEMULT2 == 0 & DOSEMULT1 > 0] <- 2
    idx[DOSEMULT2 > 0 & DOSEMULT1 == 0] <- 3
    idx[DOSEMULT2 > 0 & DOSEMULT1 > 0] <- 4
  })

  # Order nicely:
  TRTann           <- TRTann[with(TRTann, order(idx, DOSEMULT1, DOSEMULT2, DOSELEVEL1, DOSELEVEL2)),]
  TRTann$idx       <- NULL
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])
  dataPlot         <- plyr::join(dataPlot, TRTann[,c("TRTNAME","Treatment")])
  doseunit         <- dataPlot[dataPlot$NAME %in% attr(dataGen, "doseNAMES"),"UNIT"][1]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PK and PD observations ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Reduce dataset to observation of interest:
  dataPlotPKPD <- dataPlot[dataPlot$NAME %in% c(PDname, PKname) & dataPlot$MDV == 0, ]

  # Get LLOQ information:
  lloq <- unique(dataPlotPKPD[,c("Label","LLOQ","STUDY", "TYPENAME")])

  # Plot with PK and PD data:
  gr <- MMVggplot(dataPlotPKPD, aes(TIME, VALUE), ActivityPath = ActivityPath) +
    geom_hline(data=lloq, aes(yintercept = LLOQ, color = Label), linetype = 2) +
    geom_line(aes(group = interaction(USUBJID, Label)), color = "grey") +
    geom_point(aes(color = Label, shape = STUDY)) +
    #scale_color_IQRtools() +
    scale_y_log10() +
    scale_shape_manual(values = LETTERS[1:26]) +
    facet_grid(TYPENAME~Treatment,
               scales = "free") +
    labs(x = "Time [hr]",
         y = "") +
    theme(legend.position  = "bottom",
          legend.direction = "vertical")

  # Save Plot:
  IQRoutputPNG(gr, filename = file.path(filePath, "01-PKPDindivPlot.png"), width = 2+dim(TRTann)[1]*1)


  # Instead of having then in one PNG file, have them in a PDF:
  #IQRoutputPDFstart(file=file.path(filePath,"01-PKPDindivPlot.pdf"))
  grList <- list()
  for (TRT_k in unique(dataPlotPKPD$Treatment)){
    # Data subset:
    dataPlotPKPDk <- dataPlotPKPD[dataPlotPKPD$Treatment==TRT_k,]

    # Plot:
    gr <- MMVggplot(dataPlotPKPDk, aes(TIME, VALUE), ActivityPath = ActivityPath) +
      geom_hline(data=lloq, aes(yintercept = LLOQ, color = Label), linetype = 2) +
      geom_line(aes(group = interaction(USUBJID, Label)), color = "grey") +
      geom_point(aes(color = Label, shape = STUDY)) +
      #scale_color_IQRtools() +
      scale_y_log10() +
      scale_shape_manual(values = LETTERS[1:26]) +
      facet_grid(TYPENAME~Treatment,scales = "free") +
      labs(x = "Time [hr]",
           y = "") +
      theme(legend.position  = "bottom",
            legend.direction = "vertical")

    # Add plot to gList:
    grList[[TRT_k]] <- gr
  }

  # Export PDF:
  IQRoutputPDF(grList,
               filename=file.path(filePath,"01-PKPDindivPlot.pdf"))


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PK per compound per combi treatment with corresponding mono data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (k in 1:2) {
    if (k == 1) kk <- 2 else kk <- 1

    # Reduce dataset to observation of interest:
    dataPlotPK <- dataPlot[dataPlot$NAME %in% PKname[k] & dataPlot$Treatment!="Vehicle", ]

    # Treatment labels separate for drug 1 and drug 2
    dataPlotPK$TreatmentA <- with(dataPlotPK, gsub(paste0("[[:alnum:][:punct:]]*[/]?[[:alnum:]]*[ ]", get(paste0("Compound",kk))),"",Treatment))
    dataPlotPK$TreatmentA <- with(dataPlotPK, gsub("[ ][+][ ]","",TreatmentA))

    dataPlotPK$TreatmentB <- with(dataPlotPK, gsub(paste0("[[:alnum:][:punct:]]*[/]?[[:alnum:]]*[ ]", get(paste0("Compound",k)) ),"",Treatment))
    dataPlotPK$TreatmentB <- with(dataPlotPK, gsub("[ ][+][ ]","",TreatmentB))
    dataPlotPK$TreatmentB <- with(dataPlotPK, ifelse(TreatmentB %in% c("", " (bid)"),"-",TreatmentB))


    # Get LLOQ information:
    lloq <- unique(dataPlotPK[,c("Label","LLOQ","STUDY")])

    # Get better order of treatment groups
    trtInfo <- unique(dataPlotPK[, c(paste0("DOSELEVEL", k), paste0("DOSEMULT", k), "TreatmentA")])
    trtInfo <- trtInfo[order(trtInfo[,2], trtInfo[,1]),]
    dataPlotPK$TreatmentA <- factor(dataPlotPK$TreatmentA, levels = trtInfo$TreatmentA)

    # Plot with PK and PD data:
    gr <- MMVggplot(dataPlotPK, aes(TIME, VALUE), ActivityPath = ActivityPath) +
      geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
      geom_line(aes(group = interaction(USUBJID, TreatmentB), color = TreatmentB)) +
      geom_point(aes(color = TreatmentB, shape = STUDY)) +
      scale_y_log10(breaks = 10^seq(-3,4)) +
      facet_wrap(~TreatmentA) +
      labs(x = "Time [hr]",
           y = dataPlotPK$Label[1]) +
      guides(color = guide_legend(ncol=3)) +
      theme(legend.position  = "bottom",
            legend.direction = "vertical")

    # Adjust Colors:
    base::suppressMessages(gr <- gr +
                             scale_color_manual("", values = MMVcolors[2:length(MMVcolors)]))

    # Save Plot:
    IQRoutputPNG(gr,
                 filename = file.path(filePath, paste0(sprintf("%02.0f",1+k),"-PKlinePlot_",get(paste0("Compound",k)),".png")),
                 width = 9.5, height = 7)
    # Only monotherapy data:
    if (sum(dataPlotPK$TreatmentB == "-")>0) {
    grM <- gr %+% subset(dataPlotPK, TreatmentB == "-")
    IQRoutputPNG(grM,
                 filename = file.path(filePath, paste0(sprintf("%02.0f",1+k),"-a-PKlinePlotMono_",get(paste0("Compound",k)),".png")),
                 width = 9.5, height = 7)
    }
    # Only combotherapy data:
    grC <- gr %+% subset(dataPlotPK, TreatmentA %in% unique(TreatmentA[TreatmentB != "-"]))
    IQRoutputPNG(grC,
                 filename = file.path(filePath, paste0(sprintf("%02.0f",1+k),"-b-PKlinePlotCombo_",get(paste0("Compound",k)),".png")),
                 width = 9.5, height = 7)

    # Do it dose normalized
    dataPlotPK$VALUEnorm <- dataPlotPK$VALUE / dataPlotPK[[paste0("DOSELEVEL",k)]]

    # Plot:
    gr <- MMVggplot(dataPlotPK, aes(TIME, VALUEnorm), ActivityPath = ActivityPath) +
      geom_line(aes(group = interaction(USUBJID, TreatmentB), color = TreatmentB)) +
      geom_point(aes(color = TreatmentB, shape = STUDY)) +
      scale_y_log10(breaks = 10^seq(-3,4)) +
      facet_wrap(~TreatmentA) +
      scale_shape_manual(values = LETTERS[1:26]) +
      guides(color = guide_legend(ncol=3)) +
      labs(x = "Time [hr]",
           y = paste0("Dose-normalized\n",gsub("]$",paste0("/(",doseunit,")]"),dataPlotPK$Label[1]))) +
      theme(legend.position  = "bottom",
            legend.direction = "vertical")

    # Adjust Colors:
    base::suppressMessages(gr <- gr +
                             scale_color_manual("", values = MMVcolors[2:length(MMVcolors)]))

    # Save Plot:
    IQRoutputPNG(gr,
                 filename = file.path(filePath, paste0(sprintf("%02.0f",3+k),"-PKlinePlotDoseNorm_",get(paste0("Compound",k)),".png")),
                 width = 9.5, height = 7)
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PD per combo treatment ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Reduce dataset to observation of interest:
  dataPlotPD <- dataPlot[with(dataPlot, NAME == PDname & (DOSELEVEL1>1e-12 | DOSELEVEL2>1e-12))  , ]

  # Panels by combo treatment groups and add sorresponding mono therapy
  panels <- unique(dataPlotPD[with(dataPlotPD, DOSELEVEL1>1e-12 & DOSELEVEL2>1e-12) ,c("Treatment", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")])
  panels$Treatment <- droplevels(panels$Treatment)
  panels$Panel <- as.numeric(panels$Treatment)

  # Create data set by concatenating data per panel,
  # note that mono therapies may appear in different panels, such that they are concatenated multiple times
  dataPlotPD2 <- data.frame()
  for (k in seq_along(panels$Panel)) {
    dataPlotPD2 <- rbind(
      dataPlotPD2,
      within(
        rbind(
          subset(dataPlotPD, DOSELEVEL1==panels$DOSELEVEL1[k] & DOSEMULT1 == panels$DOSEMULT1[k] & DOSELEVEL2<=1e-12                 & DOSEMULT2== 0),
          subset(dataPlotPD, DOSELEVEL1<=1e-12                & DOSEMULT1 == 0                   & DOSELEVEL2==panels$DOSELEVEL2[k] & DOSEMULT2== panels$DOSEMULT2[k]),
          subset(dataPlotPD, DOSELEVEL1==panels$DOSELEVEL1[k] & DOSEMULT1 == panels$DOSEMULT1[k] & DOSELEVEL2==panels$DOSELEVEL2[k] & DOSEMULT2== panels$DOSEMULT2[k])
        )
        ,TreatmentCombo <- panels$Treatment[k])
    )
  }

  # Flag combo therapy & Compound used for TRT:
  dataPlotPD2$isCOMBO <- ifelse(with(dataPlotPD2, DOSELEVEL1>1e-12 & DOSELEVEL2>1e-12), TRUE, FALSE)
  dataPlotPD2$TRTCPD  <- ifelse(with(dataPlotPD2, DOSELEVEL1>1e-12 & DOSELEVEL2>1e-12), paste0(Compound1,"+",Compound2),
                                ifelse(with(dataPlotPD2, DOSELEVEL1>1e-12 & DOSELEVEL2<=1e-12), Compound1,
                                       ifelse(with(dataPlotPD2, DOSELEVEL1<=1e-12 & DOSELEVEL2>1e-12), Compound2,
                                              "Vehicle")))

  # Adapt Factor:
  dataPlotPD2$TRTCPD  <- factor(dataPlotPD2$TRTCPD, levels = c(Compound1,
                                                               Compound2,
                                                               paste0(Compound1,"+",Compound2)))

  # Get LLOQ information:
  lloq <- unique(dataPlotPD[,c("Label","LLOQ","STUDY")])

  # Plot PD data: With Treatment as color
  gr <- MMVggplot(dataPlotPD2, aes(TIME, VALUE), ActivityPath = ActivityPath) +
    geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
    geom_line(aes(color = Treatment, group = interaction(USUBJID, Treatment, isCOMBO), linetype = isCOMBO)) +
    geom_point(aes(color = Treatment, shape = STUDY)) +
    #scale_color_IQRtools() +
    scale_y_log10() +
    facet_wrap(~TreatmentCombo, scales = "free") +
    labs(x = "Time [hr]",
         y = dataPlotPD$Label[1]) +
    scale_shape_manual(values = LETTERS[1:26]) +
    guides(color = guide_legend(ncol=3)) +
    theme(legend.position  = "bottom",
          legend.direction = "vertical",
          legend.text      = element_text(size=7),
          strip.text       = element_text(size = 8))
  #   Save Plot:
  IQRoutputPNG(gr, filename = file.path(filePath, "06-a-PDlinePlotCombo_TRT.png"), width = 12, height = 7)

  # Plot PD data: With Compound used as color:
  gr <- MMVggplot(dataPlotPD2, aes(TIME, VALUE), ActivityPath = ActivityPath) +
    geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
    geom_line(aes(color = TRTCPD, group = interaction(USUBJID, Treatment)), linetype = 1) +
    geom_point(aes(color = TRTCPD, shape = STUDY)) +
    scale_y_log10() +
    facet_wrap(~TreatmentCombo, scales = "free") +
    labs(x = "Time [hr]",
         y = dataPlotPD$Label[1]) +
    scale_shape_manual(values = LETTERS[1:26]) +
    guides(color = guide_legend(ncol=1),
           shape = guide_legend(ncol=2)) +
    theme(legend.position  = "bottom",
          legend.direction = "vertical",
          legend.text      = element_text(size=7),
          strip.text       = element_text(size = 8))

  # Adjust Colors:
  base::suppressMessages(gr <- gr +
                           scale_color_manual(values = c("dodgerblue3", "green4", "darkorange3")))

  #   Save Plot:
  IQRoutputPNG(gr, filename = file.path(filePath, "06-b-PDlinePlotCombo_CPD.png"), width = 12, height = 7)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PD for mono treatment ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Get monotherapies
  dataPlotPD3 <- subset(dataPlotPD, (DOSELEVEL1>1e-12 & DOSELEVEL2<=1e-12) | (DOSELEVEL1<=1e-12 & DOSELEVEL2>1e-12))
  dataPlotPD3$CompoundLabel <- with(dataPlotPD3, ifelse(DOSELEVEL1>1e-12, Compound1, Compound2))

  # Do plot
  grList <- plyr::dlply(dataPlotPD3, ~CompoundLabel, function(x) {
    gr <- MMVggplot(x, aes(TIME, VALUE), ActivityPath = ActivityPath) +
      geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
      geom_line(aes(color = Treatment, group = interaction(USUBJID, Treatment))) +
      geom_point(aes(color = Treatment, shape = STUDY)) +
      #scale_color_IQRtools() +
      scale_y_log10() +
      labs(x = "Time [hr]",
           y = x$Label[1],
           title = paste0(x$CompoundLabel[1], " monotherapy")) +
      scale_shape_manual(values = LETTERS[1:26]) +
      guides(color = guide_legend(ncol=3)) +
      theme(legend.position  = "bottom",
            legend.direction = "vertical",
            legend.text      = element_text(size=7),
            strip.text       = element_text(size = 8))

    return(gr)
  })
  for (k in 1:length(grList))
    IQRoutputPNG(grList[[k]], filename = file.path(filePath, paste0("07-PDlinePlotMono_",get(paste0("Compound",k)),".png")))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot PD of untreated if exist ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Reduce dataset to observation of interest:
  dataPlotVehicle <- dataPlot[dataPlot$Treatment=="Vehicle" & dataPlot$NAME == PDname, ]

  #if ("Vehicle" %in% dataPlot$Treatment) {
  #  dataPlotVehicle <- dataPlot[dataPlot$Treatment == "Vehicle" & dataPlot$NAME == PDname, ]
  if (nrow(dataPlotVehicle)!=0) {


    # Get LLOQ information
    lloq <- unique(dataPlotVehicle[,c("Label","LLOQ","STUDY")])

    # Plot with PK and PD data
    gr <- MMVggplot(dataPlotVehicle, aes(TIME, VALUE), ActivityPath = ActivityPath) +
      geom_hline(data=lloq, aes(yintercept = LLOQ), linetype = 2, color = "black") +
      geom_line(aes(group = USUBJID), color = "grey") +
      geom_point(aes(color = Treatment, shape = STUDY)) +
      #scale_color_IQRtools() +
      scale_shape_manual(values = LETTERS[1:26]) +
      scale_y_log10() +
      labs(x = "Time [hr]",
           y = dataPlotVehicle$Label[1])

    # Save Plot:
    IQRoutputPNG(gr, filename = file.path(filePath, "08-UntreatedPDlinePlot.png"))
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot with erythrocytes if given ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (!is.null(HuErys)) {
    dataPlotHuEry <- dataPlot[dataPlot$NAME %in% HuErys, ]

    # Plot:
    gr <- MMVggplot(dataPlotHuEry, aes(TIME, VALUE), ActivityPath = ActivityPath) +
      geom_line(aes(group = USUBJID), color = "grey") +
      geom_point(aes(color = Treatment, shape = STUDY)) +
      # scale_y_log10() +
      # coord_cartesian(ylim = c(1,140)) +
      # facet_grid(.~STUDY, scales = "free") +
      labs(x = "Time [hr]",
           y = dataPlotHuEry$Label[1])

    # Adjust Colors:
    base::suppressMessages(gr <- gr +
                             scale_color_manual(limits = as.character(TRTann$Treatment),
                                                values = MMVcolors[2:length(MMVcolors)]))

    # Adjust Shapes:
    if (length(unique(dataPlotHuEry$Treatment)) > 6)
      gr <- gr + scale_shape_manual(values = LETTERS[1:26])

    # Save Plot:
    IQRoutputPNG(gr, filename = file.path(filePath, "09-HuErythrocyteslinePlot.png"))
  }
}

#' summary_PDdata_byStudyTRT
#'
#' @description
#' @param dataGen
#' @param filename Default: NULL
#' @param EfficacyNAME Default: 'Parasitemia Asexual + Gametocyte Mean'
#' @param NewTRTNAME Default: NULL
#' @param CovToTRTNAME Default: NULL
#' @param CovNewTRTNAME Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join
#' @importFrom tidyr spread
summary_PDdata_byStudyTRT <- function(dataGen,
                                      filename     = NULL,
                                      EfficacyNAME = "Parasitemia Asexual + Gametocyte Mean",
                                      NewTRTNAME   = NULL,
                                      CovToTRTNAME  = NULL,
                                      CovNewTRTNAME = NULL) {

  # Change TRTNAME by adding Covariable(s):
  if (!is.null(CovToTRTNAME)){
    for (k in 1:length(CovToTRTNAME)){
      dataGen$TRTNAME <- paste0(dataGen$TRTNAME, "_", dataGen[,CovToTRTNAME[k]])
    }
  }

  # Transform dataset to dataFrame:
  dataStat <- as.data.frame(dataGen)

  # Get nicer treatment name ----
  TRTann <- within(unique(dataStat[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME

    # Change name for DOSE=0: Vehicle was chosen as it was first developped for SCID
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)

    # Change the name of other Treatment if desired:
    for (TRT_k in Treatment){
      idx_k = which(Treatment==TRT_k)
      if (!is.null(NewTRTNAME[[TRT_k]])){
        Treatment[idx_k] = NewTRTNAME[[TRT_k]]
      }
    }
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])

  # Add Treatment Name:
  dataStat <- plyr::join(dataStat, TRTann)

  # PD data:
  dataStatPD <- dataStat[dataStat$TYPENAME=="Efficacy" & dataStat$NAME==EfficacyNAME,]

  # Table on summarizing PD Data ----
  #   Nbr of Subj. per trial and Dose Group:
  tmp             <- as.data.frame(with(unique(dataStatPD[,c("USUBJID","STUDY","Treatment")]), table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableGR         <- tidyr::spread(tmp, Study, Freq)

  #   Nbr of PD Ob. per trial and Dose Group:
  tmp             <- as.data.frame(with(dataStatPD[,c("STUDY","Treatment")], table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tablePD         <- tidyr::spread(tmp, Study, Freq)

  #   Nbr of PD Ob. below LLOQ per trial and Dose Group:
  tmp             <- as.data.frame(with(dataStatPD[dataStatPD$VALUE<dataStatPD$LLOQ,c("STUDY","Treatment")], table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableLLOQ         <- tidyr::spread(tmp, Study, Freq)

  #   Look for the LLOQ of eqch study, qnd convert to ng/mL is possible:
  LLOQ <- unique(dataStatPD[dataStatPD$TYPENAME=="Efficacy" & dataStatPD$NAME==EfficacyNAME,c("STUDY","LLOQ","UNIT")])

  #   Create dataframe to be saved:
  n_TRT   = length(unique(dataStatPD$Treatment))
  n_STUDY = length(unique(dataStatPD$STUDY))
  n_ROW   = n_STUDY*n_TRT+2
  PDsummary <- data.frame(STUDY     = character(n_ROW),
                          TRTNAME   = character(n_ROW),
                          Nsubj     = character(n_ROW),
                          PDobs     = character(n_ROW),
                          PDobsLLOQ = character(n_ROW),
                          LLOQ      = character(n_ROW),
                          stringsAsFactors = FALSE
  )

  #   Define Study and Treatment for looping:
  STUDY     <- unique(dataStatPD[order(dataStatPD$STUDY),]$STUDY)
  Treatment <- unique(dataStatPD[order(dataStatPD$DOSELEVEL),]$Treatment)

  # Loop over Study:
  for (i_STUDY in 1:length(STUDY)){
    # Get Study Name:
    STUDY_i <- STUDY[i_STUDY]

    # Loop over Treatment
    HeadStudy = TRUE
    for (j_TRT in 1:length(Treatment)){
      #Get Treatment Name:
      Treatment_j <- Treatment[j_TRT]

      # Save information into PDsummary:
      PDsummary$TRTNAME[n_TRT*(i_STUDY-1)+j_TRT] <- as.character(Treatment_j)
      PDsummary$Nsubj[n_TRT*(i_STUDY-1)+j_TRT] <- tableGR[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      PDsummary$PDobs[n_TRT*(i_STUDY-1)+j_TRT] <- tablePD[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      PDsummary$PDobsLLOQ[n_TRT*(i_STUDY-1)+j_TRT] <- tableLLOQ[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      if (PDsummary$Nsubj[n_TRT*(i_STUDY-1)+j_TRT]!=0 & HeadStudy){
        PDsummary$STUDY[n_TRT*(i_STUDY-1)+j_TRT] <- STUDY_i
        PDsummary$LLOQ[n_TRT*(i_STUDY-1)+j_TRT] <- paste0(LLOQ$LLOQ[LLOQ$STUDY==STUDY_i]," ",LLOQ$UNIT[LLOQ$STUDY==STUDY_i])
        HeadStudy = FALSE
      }


    }
  }

  #   Adjust last row:
  PDsummary$STUDY[n_ROW] = "Total"
  PDsummary$Nsubj[n_ROW] = sum(as.numeric(PDsummary$Nsubj[1:(n_ROW-2)]))
  PDsummary$PDobs[n_ROW] = sum(as.numeric(PDsummary$PDobs[1:(n_ROW-2)]))
  PDsummary$PDobsLLOQ[n_ROW] = sum(as.numeric(PDsummary$PDobsLLOQ[1:(n_ROW-2)]))

  #   Adjust PDobs with LLOQ information:
  PDsummary$PDobs = paste0(PDsummary$PDobs, " (",PDsummary$PDobsLLOQ,")")
  PDsummary$PDobs[n_ROW-1] = ""

  #   Keep only columns and lines of interest
  PDsummary <- PDsummary[,c("STUDY","TRTNAME","Nsubj","PDobs","LLOQ")]
  PDsummary <- PDsummary[PDsummary$Nsubj!="0",]

  # Change Column Name:
  names(PDsummary) <- c("Study", "Treatment Groups", "N subjects", "Total PD Observation (<LLOQ)","LLOQ")

  # Generate Extra columns to split the treatment group, if needed:
  if (!is.null(CovToTRTNAME)){
    # Dose Group Column:
    PDsummary$DOSE <- ""

    # Extra Column:
    for (k in 1:length(CovToTRTNAME)){
      PDsummary[,CovToTRTNAME[k]] <- ""
    }

    # Fill out columns
    for (TRT_k in PDsummary[,'Treatment Groups']){
      idx <- which((names(CovNewTRTNAME)==TRT_k))
      if (length(idx)!=0){
        PDsummary$DOSE[PDsummary[,'Treatment Groups']==TRT_k] <- CovNewTRTNAME[[idx]][1]
        for (i in 1:length(CovToTRTNAME)){
          PDsummary[PDsummary[,'Treatment Groups']==TRT_k,CovToTRTNAME[i]] <- CovNewTRTNAME[[idx]][i+1]
        }
      }
    }

    #   Re-organize table
    PDsummary        <- PDsummary[,c("Study", "DOSE", CovToTRTNAME, "N subjects", "Total PD Observation (<LLOQ)", "LLOQ")]
    if (!is.null(names(CovToTRTNAME))){
      names(PDsummary) <- c("Study", "Dose Groups", names(CovToTRTNAME), "N subjects", "Total PD Observation (<LLOQ)", "LLOQ")
    }else{
      names(PDsummary) <- c("Study", "Dose Groups", CovToTRTNAME, "N subjects", "Total PD Observation (<LLOQ)", "LLOQ")
    }
  }

  # Save Table:
  if (!is.null(filename)){
    IQRoutputTable(PDsummary,
                   xtitle   = "PD data summarized by study and treatment group included in the PD analysis",
                   xfooter  = paste0("N: Number of<br>PD Data: ",EfficacyNAME),
                   filename = filename,
                   report   = TRUE)
  }

  # Output:
  return(PDsummary)
}

#' summary_PKdata_byStudyTRT
#'
#' @description
#' @param dataGen
#' @param filename Default: NULL
#' @param NewTRTNAME Default: NULL
#' @param CovToTRTNAME Default: NULL
#' @param CovNewTRTNAME Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join
#' @importFrom tidyr spread
summary_PKdata_byStudyTRT <- function(dataGen,
                                      filename      = NULL,
                                      NewTRTNAME    = NULL,
                                      CovToTRTNAME  = NULL,
                                      CovNewTRTNAME = NULL) {

  # Change TRTNAME by adding Covariable(s):
  if (!is.null(CovToTRTNAME)){
    for (k in 1:length(CovToTRTNAME)){
      dataGen$TRTNAME <- paste0(dataGen$TRTNAME, "_", dataGen[,CovToTRTNAME[k]])
    }
  }

  # Transform dataset to dataFrame:
  dataStat <- as.data.frame(dataGen)

  # Get nicer treatment name ----
  TRTann <- within(unique(dataStat[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME

    # Change name for DOSE=0: Vehicle was chosen as it was first developped for SCID
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)

    # Change the name of other Treatment if desired:
    for (TRT_k in Treatment){
      idx_k = which(Treatment==TRT_k)
      if (!is.null(NewTRTNAME[[TRT_k]])){
        Treatment[idx_k] = NewTRTNAME[[TRT_k]]
      }
    }
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])

  # Add Treatment Name:
  dataStat <- plyr::join(dataStat, TRTann)

  # Check that STUDY is not numeric ----
  if (!class(dataStat$STUDY) %in% c("factor", "character")) {
    dataStat$STUDY <- as.character(dataStat$STUDY)
  }

  # Table on summarizing PK Data ----
  #   Nbr of Subj. per trial and Dose Group:
  tmp             <- as.data.frame(with(unique(dataStat[,c("USUBJID","STUDY","Treatment")]), table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableGR         <- tidyr::spread(tmp, Study, Freq)

  #   Nbr of PK Ob. per trial and Dose Group:
  tmp             <- as.data.frame(with(dataStat[dataStat$TYPENAME=="PK",c("STUDY","Treatment")], table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tablePK         <- tidyr::spread(tmp, Study, Freq)

  #   Nbr of PK Ob. below LLOQ per trial and Dose Group:
  tmp             <- as.data.frame(with(dataStat[dataStat$TYPENAME=="PK" & dataStat$VALUE<dataStat$LLOQ,c("STUDY","Treatment")], table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableLLOQ         <- tidyr::spread(tmp, Study, Freq)

  #   Look for the LLOQ of eqch study, qnd convert to ng/mL is possible:
  LLOQ <- unique(dataStat[dataStat$TYPENAME=="PK",c("STUDY","LLOQ","UNIT")])
  if (LLOQ$UNIT[1]=="ug/mL"){
    LLOQ$LLOQ = LLOQ$LLOQ*1000
    LLOQ$UNIT = "ng/mL"
  }

  #   Create dataframe to be saved:
  n_TRT   = length(unique(dataStat$Treatment))
  n_STUDY = length(unique(dataStat$STUDY))
  n_ROW   = n_STUDY*n_TRT+2
  PKsummary <- data.frame(STUDY     = character(n_ROW),
                          TRTNAME   = character(n_ROW),
                          Nsubj     = character(n_ROW),
                          PKobs     = character(n_ROW),
                          PKobsLLOQ = character(n_ROW),
                          LLOQ      = character(n_ROW),
                          stringsAsFactors = FALSE
  )

  # Define Study and Treatment for looping:
  STUDY     <- unique(dataStat[order(dataStat$STUDY),]$STUDY)
  Treatment <- unique(dataStat[order(dataStat$DOSELEVEL),]$Treatment)

  # Loop over Study:
  for (i_STUDY in 1:length(STUDY)){
    # Get Study Name:
    STUDY_i <- STUDY[i_STUDY]

    # Loop over Treatment
    HeadStudy = TRUE
    for (j_TRT in 1:length(Treatment)){
      #Get Treatment Name:
      Treatment_j <- Treatment[j_TRT]

      # Save information into PKsummary:
      PKsummary$TRTNAME[n_TRT*(i_STUDY-1)+j_TRT] <- as.character(Treatment_j)
      PKsummary$Nsubj[n_TRT*(i_STUDY-1)+j_TRT] <- tableGR[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      PKsummary$PKobs[n_TRT*(i_STUDY-1)+j_TRT] <- tablePK[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      PKsummary$PKobsLLOQ[n_TRT*(i_STUDY-1)+j_TRT] <- tableLLOQ[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      if (PKsummary$Nsubj[n_TRT*(i_STUDY-1)+j_TRT]!=0 & HeadStudy){
        PKsummary$STUDY[n_TRT*(i_STUDY-1)+j_TRT] <- STUDY_i
        PKsummary$LLOQ[n_TRT*(i_STUDY-1)+j_TRT] <- paste0(LLOQ$LLOQ[LLOQ$STUDY==STUDY_i]," ",LLOQ$UNIT[LLOQ$STUDY==STUDY_i])
        HeadStudy = FALSE
      }


    }
  }

  #   Adjust last row:
  PKsummary$STUDY[n_ROW] = "Total"
  PKsummary$Nsubj[n_ROW] = sum(as.numeric(PKsummary$Nsubj[1:(n_ROW-2)]))
  PKsummary$PKobs[n_ROW] = sum(as.numeric(PKsummary$PKobs[1:(n_ROW-2)]))
  PKsummary$PKobsLLOQ[n_ROW] = sum(as.numeric(PKsummary$PKobsLLOQ[1:(n_ROW-2)]))

  #   Adjust PKobs with LLOQ information:
  PKsummary$PKobs = paste0(PKsummary$PKobs, " (",PKsummary$PKobsLLOQ,")")
  PKsummary$PKobs[n_ROW-1] = ""

  #   Keep only columns and lines of interest
  PKsummary <- PKsummary[,c("STUDY","TRTNAME","Nsubj","PKobs","LLOQ")]
  PKsummary <- PKsummary[PKsummary$Nsubj!="0",]

  # Change Column Name:
  names(PKsummary) <- c("Study", "Treatment Groups", "N subjects", "Total PK Observation (<LLOQ)","LLOQ")


  # Generate Extra columns to split the treatment group, if needed:
  if (!is.null(CovToTRTNAME)){
    # Dose Group Column:
    PKsummary$DOSE <- ""

    # Extra Column:
    for (k in 1:length(CovToTRTNAME)){
      PKsummary[,CovToTRTNAME[k]] <- ""
    }

    # Fill out columns
    for (TRT_k in PKsummary[,'Treatment Groups']){
      idx <- which((names(CovNewTRTNAME)==TRT_k))
      if (length(idx)!=0){
        PKsummary$DOSE[PKsummary[,'Treatment Groups']==TRT_k] <- CovNewTRTNAME[[idx]][1]
        for (i in 1:length(CovToTRTNAME)){
          PKsummary[PKsummary[,'Treatment Groups']==TRT_k,CovToTRTNAME[i]] <- CovNewTRTNAME[[idx]][i+1]
        }
      }
    }

    #   Re-organize table
    PKsummary        <- PKsummary[,c("Study", "DOSE", CovToTRTNAME, "N subjects", "Total PK Observation (<LLOQ)", "LLOQ")]
    if (!is.null(names(CovToTRTNAME))){
      names(PKsummary) <- c("Study", "Dose Groups", names(CovToTRTNAME), "N subjects", "Total PK Observation (<LLOQ)", "LLOQ")
    }else{
      names(PKsummary) <- c("Study", "Dose Groups", CovToTRTNAME, "N subjects", "Total PK Observation (<LLOQ)", "LLOQ")
    }
  }

  # Save Table:
  if (!is.null(filename)){
    IQRoutputTable(PKsummary,
                   xtitle   = "PK data summarized by study and treatment group included in the PK analysis",
                   xfooter  = "N: Number of",
                   filename = filename,
                   report   = TRUE)
  }

  # Output:
  return(PKsummary)
}

#' summary_PKdataMMV
#'
#' @description
#' @param dataGen
#' @param filePath Default: NULL
#' @param NewTRTNAME Default: NULL
#' @return
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join ddply
#' @importFrom tidyr spread
summary_PKdataMMV <- function(dataGen,
                              filePath   = NULL,
                              NewTRTNAME = NULL) {

  # Transform dataset to dataFrame:
  dataStat <- as.data.frame(dataGen)

  # Get nicer treatment name ----
  TRTann <- within(unique(dataStat[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME

    # Change name for DOSE=0: Vehicle was chosen as it was first developped for SCID
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)

    # Change the name of other Treatment if desired:
    for (TRT_k in Treatment){
      idx_k = which(Treatment==TRT_k)
      if (!is.null(NewTRTNAME[[TRT_k]])){
        Treatment[idx_k] = NewTRTNAME[[TRT_k]]
      }
    }
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])

  # Add Treatment Name:
  dataStat <- plyr::join(dataStat, TRTann)


  # Count # of subjects per treatment group and study ----
  tmp             <- as.data.frame(with(unique(dataStat[,c("USUBJID","STUDY","Treatment")]), table(Treatment,STUDY)), stringsAsFactors = FALSE)
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableGR         <- tidyr::spread(tmp, Study, Freq)
  tableGR$Total   <- rowSums(as.data.frame(tableGR[,as.character(unique(tmp$Study))]))
  TotRow          <- c("Total",as.numeric(colSums(as.data.frame(tableGR[,as.character(unique(tmp$Study))]))),sum(tableGR$Total))
  names(TotRow)   <- names(tableGR)
  TotRow          <- as.data.frame(t(TotRow))
  EmptyRow        <- character(length=length(TotRow))
  names(EmptyRow) <- names(tableGR)
  EmptyRow        <- as.data.frame(t(EmptyRow))
  tableGR         <- rbind(tableGR,EmptyRow,TotRow)
  tableGR[tableGR == 0] <- "."
  IQRoutputTable(tableGR,
                 xtitle   = "Number of subjects per study and treatment group",
                 xfooter  = "",
                 filename = file.path(filePath, "01_NumSUBJ_byTRT.txt"),
                 report   = TRUE)


  # Table on sampling/observations ----
  tmp <- as.data.frame(with(subset(dataStat, NAME %in% attr(dataGen, "obsNAMES")), table(NAME,NT)), stringsAsFactors = FALSE)
  names(tmp) <- c("Observation", "NT", "Freq")
  tmp <- tmp[tmp$Freq != 0,]
  tableOBS <- plyr::ddply(tmp, ~Observation, function(x) {
    x$xx <- paste0(x$NT, "h (", x$Freq, ")")
    # x$xx <- paste0(x$NT, "h")
    out <- data.frame('Sampling Times' =  paste0(x$xx, collapse = ", "))
    out
  })

  names(tableOBS)[2] <- "Sampling Times"
  IQRoutputTable(tableOBS,
                 xtitle = "Nominal sampling times per observation type",
                 xfooter = "Number of samples given in parenthesis",
                 filename = file.path(filePath, "02_NumObs_NT.txt"),
                 report = TRUE)


  # Table on dosing ----
  tableDOS = as.data.frame(levels(dataStat$Treatment))
  names(tableDOS) = c("Treatment.Group")
  nameTRT  = tableDOS$Treatment.Group
  for (listNAME in attr(dataGen, "doseNAMES")){
    tmp <- as.data.frame(with(subset(dataStat,NAME %in% listNAME), table(Treatment,NT)), stringsAsFactors = FALSE)
    names(tmp) <- c("Treatment.Group", "NT", "Freq")
    tmp <- tmp[tmp$Freq != 0,]
    tmp2 <- plyr::ddply(tmp, ~Treatment.Group, function(x) {
      x$xx <- paste0(x$NT, "h (", x$Freq, ")")
      # x$xx <- paste0(x$NT, "h")
      out <- data.frame('Dosing times' =  paste0(x$xx, collapse = ", "))
      out
    })
    names(tmp2) <- c("Treatment.Group", paste0(listNAME, " - Time"))

    tableDOS <- merge(tableDOS,tmp2,by="Treatment.Group",all.x=TRUE)
  }
  tableDOS <- tableDOS[nameTRT,]
  for (col in 1:ncol(tableDOS)){
    tableDOS[,col]=as.character(tableDOS[,col])
  }
  tableDOS[is.na(tableDOS)] <- "."

  names(tableDOS)[1] <- "Treatment Group"

  # Save Table:
  if (!is.null(filePath)){
    IQRoutputTable(tableDOS,
                   xtitle = "Dosing time in NT per treatment group and dosing/inoculation type",
                   xfooter = "Number of subjects in parenthesis",
                   filename = file.path(filePath, "03_DosingTime_NT.txt"),
                   report = TRUE)
  }

  # Output:
  return(tableDOS)

}
#' summary_PKPDdataCombo
#'
#' @description
#' @param dataGen
#' @param filePath Default: NULL
#' @param CompoundList Default: NULL
#' @param positiveQC Default: NULL
#' @return
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join ddply
#' @importFrom tidyr spread
summary_PKPDdataCombo <- function(dataGen,
                                  filePath     = NULL,
                                  CompoundList = NULL,
                                  positiveQC   = NULL,
                                  TRTNAME      = "TRTNAME") {

  # Update Compound Name if wanted:
  if (!is.null(CompoundList)){
    dataGen        <- swapName_MMVnameToName(dataGen,
                                             CompoundList)
  }

  # As dataFrame:
  dataStat <- as.data.frame(dataGen)

  # Adjust TRTNAME column:
  if(TRTNAME!="TRTNAME"){
    dataStat$TRTNAME <- NULL
    names(dataStat)[names(dataStat)==TRTNAME] <- "TRTNAME"
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get nicer treatment names and order ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Change Names to look better:
  TRTann <- within(unique(dataStat[,c("TRTNAME", "STUDY", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")]), {
    Treatment <- TRTNAME
    Treatment <- gsub("[+]", " + ", Treatment)
    Treatment <- gsub("mg/kg", "mg/kg ", Treatment)
    idx <- 0
    idx[TRTNAME=="Vehicle"] <- 1
    # idx[DOSELEVEL1<=1e-12 & DOSELEVEL2<=1e-12 ] <- 1
    idx[DOSEMULT2==0 & DOSEMULT1>0 ] <- 3
    idx[DOSEMULT2>0  & DOSEMULT1==0] <- 4
    idx[DOSEMULT2>0  & DOSEMULT1>0 ] <- 5
    idx <- ifelse(grepl(paste0(positiveQC,collapse = " | "), TRTNAME), 2, idx)
  })

  # Order nicely:
  TRTann           <- TRTann[with(TRTann, order(idx, DOSEMULT1, DOSEMULT2, DOSELEVEL1, DOSELEVEL2)),]
  TRTann$idx       <- NULL
  if(is.factor(TRTann$TRTNAME)){
    levels_TRT <- unique(TRTann$Treatment)[match(levels(TRTann$TRTNAME), unique(TRTann$TRTNAME))]
    TRTann$Treatment <- factor(TRTann$Treatment, levels = levels_TRT)
  }else{
    TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])
  }
  dataStat         <- plyr::join(dataStat, TRTann)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # TABLE 1: Count # of subjects per treatment group and study ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp             <- as.data.frame(with(unique(dataStat[,c("USUBJID","STUDY","Treatment")]), table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableGR         <- tidyr::spread(tmp, Study, Freq)
  tableGR$Total   <- rowSums(as.data.frame(tableGR[,as.character(unique(tmp$Study))]))
  TotRow          <- c("Total",as.numeric(colSums(as.data.frame(tableGR[,as.character(unique(tmp$Study))]))),sum(tableGR$Total))
  names(TotRow)   <- names(tableGR)
  TotRow          <- as.data.frame(t(TotRow))
  EmptyRow        <- character(length=length(TotRow))
  names(EmptyRow) <- names(tableGR)
  EmptyRow        <- as.data.frame(t(EmptyRow))
  tableGR         <- rbind(tableGR,EmptyRow,TotRow)
  tableGR[tableGR == 0] <- "."
  IQRoutputTable(tableGR,
                 xtitle   = "Number of subjects per study and treatment group",
                 xfooter  = "",
                 filename = file.path(filePath, "01-NumSUBJ_byTRT.txt"),
                 report   = TRUE)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # TABLE 2: Dosing Table ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tableDOS        <- as.data.frame(levels(dataStat$Treatment))
  names(tableDOS) <- c("Treatment.Group")
  nameTRT         <- tableDOS$Treatment.Group
  for (listNAME in attr(dataGen, "doseNAMES")){
    tmp        <- as.data.frame(with(subset(dataStat,NAME %in% listNAME), table(Treatment,NT)))
    names(tmp) <- c("Treatment.Group", "NT", "Freq")
    tmp        <- tmp[tmp$Freq != 0,]
    tmp2       <- plyr::ddply(tmp, ~Treatment.Group, function(x) {
      x$xx <- paste0(x$NT, "h (", x$Freq, ")")
      # x$xx <- paste0(x$NT, "h")
      out  <- data.frame("Dosing Times"=paste0(x$xx, collapse = ", "))
      return(out)
    }
    )
    names(tmp2) <- c("Treatment.Group", paste0(listNAME," Time"))
    tableDOS    <- merge(tableDOS,tmp2,by="Treatment.Group",all.x=TRUE)
  }

  tableDOS <- tableDOS[nameTRT,]

  for (col in 1:ncol(tableDOS)){
    tableDOS[,col] <- as.character(tableDOS[,col])
  }

  # Replace NA cell with ".":
  tableDOS[is.na(tableDOS)] <- "."

  # Update Treatment Column Name:
  names(tableDOS)[1] <- "Treatment Group"

  # Replace Administration Time of Vehicle woth ".":
  # If in the data set there is a dose event with a dose of 0 it would appear
  idx_Vehicle <- (tableDOS[,"Treatment Group"]=="Vehicle")
  tableDOS[idx_Vehicle,2:length(tableDOS)] <- "."

  # Save Table:
  if (!is.null(filePath)){
    IQRoutputTable(tableDOS,
                   xtitle = "Dosing/Inoculation time in NT per treatment group and dosing/inoculation type",
                   xfooter = "Number of subjects in parenthesis",
                   filename = file.path(filePath, "02-DosingInoculationTime_NT.txt"),
                   report = TRUE)
  }

  # Output:
  return(tableDOS)
}
#' summary_PKPDdataMMV
#' Create summary tables for pharmacokinetic data
#'
#' Summarizes pharmacokinetic data in 3 different tables: (1) number of subjects per study and treatment group, (2) nominal sampling times per observation type and (3) dosing/inoculation nominal time per treatment group and dosing/inoculation type.
#'
#' @param dataGen
#' @param filePath
#' @param NewTRTNAME
#' @param positiveQC
#'
#' @return
#' @export
#'
#' @examples
#' @author Aline Fuchs (MMV), catalinbarcelo (MMV), Mohammed H. Cherkaoui (MMV)
summary_PKPDdataMMV <- function(dataGen,
                                filePath   = NULL,
                                NewTRTNAME = NULL,
                                positiveQC = NULL) {

  # Transform dataset to dataFrame:
  dataStat <- as.data.frame(dataGen)

  # Get nicer treatment name ----
  TRTann <- within(unique(dataStat[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME

    # Change name for DOSE=0: Vehicle was chosen as it was first developped for SCID
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)

    # Change the name of other Treatment if desired:
    for (TRT_k in Treatment){
      idx_k <- which(Treatment==TRT_k)
      if (!is.null(NewTRTNAME[[TRT_k]])){
        Treatment[idx_k] <- NewTRTNAME[[TRT_k]]
      }
    }
  })

  # Order
  TRTann$idx <- 0
  TRTann$idx[TRTann$DOSEMULT == 0 ] <- 3
  TRTann$idx[TRTann$DOSEMULT > 0] <- 4
  TRTann$idx[TRTann$TRTNAME == "Vehicle"] <- 1
  TRTann$idx <- ifelse(grepl(paste0(positiveQC,collapse = " | "), TRTann$TRTNAME), 2, TRTann$idx)

  TRTann <- TRTann[with(TRTann, order(idx,DOSEMULT, DOSELEVEL)),]
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])

  # Add Treatment Name:
  dataStat <- plyr::join(dataStat, TRTann)

  # Count # of subjects per treatment group and study ----
  tmp             <- as.data.frame(with(unique(dataStat[,c("USUBJID","STUDY","Treatment")]), table(Treatment,STUDY)), stringsAsFactors = FALSE)
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableTRT        <- tidyr::spread(tmp, Study, Freq)
  tableTRT$Total  <- rowSums(tableTRT[,as.character(unique(tmp$Study)), drop = FALSE], na.rm = TRUE)
  TotRow          <- c("Total",as.numeric(colSums(tableTRT[,as.character(unique(tmp$Study)), drop = FALSE])),sum(tableTRT$Total))
  names(TotRow)   <- names(tableTRT)
  TotRow          <- as.data.frame(t(TotRow))
  EmptyRow        <- character(length=length(TotRow))
  names(EmptyRow) <- names(tableTRT)
  EmptyRow        <- as.data.frame(t(EmptyRow))
  tableTRT[tableTRT == 0] <- "."
  tableTRT        <- tableTRT[with(tableTRT, order(levels(TRTann$Treatment))),]
  tableTRT        <- rbind(tableTRT,EmptyRow,TotRow)
  rownames(tableTRT) <- NULL

  IQRoutputTable(tableTRT,
                 xtitle   = "Number of subjects per study and treatment group",
                 xfooter  = "",
                 filename = file.path(filePath, "01-NumSUBJ_byTRT.txt"),
                 report   = TRUE)


  # Table on sampling/observations ----
  tmp <- as.data.frame(with(subset(dataStat, NAME %in% attr(dataGen, "obsNAMES")), table(NAME,NT)), stringsAsFactors = FALSE)
  names(tmp) <- c("Observation", "NT", "Freq")
  tmp <- tmp[tmp$Freq != 0,]
  tableOBS <- plyr::ddply(tmp, ~Observation, function(x) {
    x$xx <- paste0(x$NT, "h (", x$Freq, ")")
    # x$xx <- paste0(x$NT, "h")
    out <- data.frame('Sampling Times' =  paste0(x$xx, collapse = ", "))
    out
  })

  names(tableOBS)[2] <- "Sampling Times"
  IQRoutputTable(tableOBS,
                 xtitle = "Nominal sampling times per observation type",
                 xfooter = "Number of samples given in parenthesis",
                 filename = file.path(filePath, "02-NumObs_NT.txt"),
                 report = TRUE)


  # Table on dosing ----
  tableDOSE = as.data.frame(levels(dataStat$Treatment))
  names(tableDOSE) = c("Treatment.Group")
  nameTRT  = tableDOSE$Treatment.Group
  for (listNAME in attr(dataGen, "doseNAMES")){
    tmp <- as.data.frame(with(subset(dataStat,NAME %in% listNAME), table(Treatment,NT)), stringsAsFactors = FALSE)
    names(tmp) <- c("Treatment.Group", "NT", "Freq")
    tmp <- tmp[tmp$Freq != 0,]
    tmp2 <- plyr::ddply(tmp, ~Treatment.Group, function(x) {
      x$xx <- paste0(x$NT, "h (", x$Freq, ")")
      # x$xx <- paste0(x$NT, "h")
      out <- data.frame('Dosing times' =  paste0(x$xx, collapse = ", "))
      out
    })
    names(tmp2) <- c("Treatment.Group", paste0(listNAME, " - Time"))

    tableDOSE <- merge(tableDOSE,tmp2,by="Treatment.Group",all.x=TRUE)
  }
  tableDOSE <- tableDOSE[nameTRT,]
  for (col in 1:ncol(tableDOSE)){
    tableDOSE[,col]=as.character(tableDOSE[,col])
  }
  tableDOSE[is.na(tableDOSE)] <- "."

  names(tableDOSE)[1] <- "Treatment Group"

  # Save Table:
  if (!is.null(filePath)){
    IQRoutputTable(tableDOSE,
                   xtitle = "Dosing/Inoculation time in NT per treatment group and dosing/inoculation type",
                   xfooter = "Number of subjects in parenthesis",
                   filename = file.path(filePath, "03-DosingInoculationTime_NT.txt"),
                   report = TRUE)
  }

  # Output:
  return(list(tableTRT  = tableTRT,
              tableOBS  = tableOBS,
              tableDOSE = tableDOSE))
}

#' summary_PKPDdataSimplified
#'
#' @description
#' @param dataGen
#' @param filename Default: NULL
#' @param PKname Default: NULL
#' @param PDname Default: NULL
#' @param StudyTYPE
#' @param Center
#' @param NewTRTNAME Default: NULL
#' @param CovToTRTNAME Default: NULL
#' @param CovNewTRTNAME Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
#' @importFrom plyr join
#' @importFrom tidyr spread
summary_PKPDdataSimplified <- function(dataGen,
                                       filename      = NULL,
                                       PKname        = NULL,
                                       PDname        = NULL,
                                       StudyTYPE,
                                       Center,
                                       NewTRTNAME    = NULL,
                                       CovToTRTNAME  = NULL,
                                       CovNewTRTNAME = NULL) {

  # Change TRTNAME by adding Covariable(s):
  if (!is.null(CovToTRTNAME)){
    for (k in 1:length(CovToTRTNAME)){
      dataGen$TRTNAME <- paste0(dataGen$TRTNAME, "_", dataGen[,CovToTRTNAME[k]])
    }
  }

  # Transform dataset to dataFrame:
  dataStat <- as.data.frame(dataGen)

  # Get nicer treatment name ----
  TRTann <- within(unique(dataStat[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME

    # Change name for DOSE=0: Vehicle was chosen as it was first developped for SCID
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)

    # Change the name of other Treatment if desired:
    for (TRT_k in Treatment){
      idx_k = which(Treatment==TRT_k)
      if (!is.null(NewTRTNAME[[TRT_k]])){
        Treatment[idx_k] = NewTRTNAME[[TRT_k]]
      }
    }
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]
  TRTann$Treatment <- factor(TRTann$Treatment, levels = TRTann$Treatment[!duplicated(TRTann$Treatment)])

  # Add Treatment Name:
  dataStat <- plyr::join(dataStat, TRTann)


  # Table on summarizing PK Data ----
  #   Nbr of Subj. per trial and Dose Group:
  tmp             <- as.data.frame(with(unique(dataStat[,c("USUBJID","STUDY","Treatment")]), table(Treatment,STUDY)))
  names(tmp)      <- c("Treatment Group", "Study", "Freq")
  tableGR         <- tidyr::spread(tmp, Study, Freq)

  #   Create dataframe to be saved:
  n_TRT   = length(unique(dataStat$Treatment))
  n_STUDY = length(unique(dataStat$STUDY))
  n_ROW   = n_STUDY*n_TRT+2
  PKPDsummary <- data.frame(STUDY     = character(n_ROW),
                            StudyTYPE = character(n_ROW),
                            TRTNAME   = character(n_ROW),
                            Nsubj     = character(n_ROW),
                            Analysis  = character(n_ROW),
                            stringsAsFactors = FALSE
  )

  #   Define Study and Treatment for looping:
  STUDY     <- unique(dataStat[order(dataStat$STUDY),]$STUDY)
  Treatment <- unique(dataStat[order(dataStat$DOSELEVEL),]$Treatment)


  # Loop over Study:
  for (i_STUDY in 1:length(STUDY)){
    # Get Study Name:
    STUDY_i <- STUDY[i_STUDY]

    #PK or/and PD study?
    dataStat_STUDY_i <- dataStat[dataStat$STUDY==STUDY_i,]
    #   PK study?
    PKstudy = FALSE
    if (!is.null(PKname)){
      if (PKname %in% unique(dataStat_STUDY_i$NAME)){
        PKstudy = TRUE
      }
    } else{
      if ("PK" %in% unique(dataStat_STUDY_i$TYPENAME)){
        PKstudy = TRUE
      }
    }
    # PD study?
    PDstudy = FALSE
    if (!is.null(PDname)){
      if (PDname %in% unique(dataStat_STUDY_i$NAME)){
        PDstudy = TRUE
      }
    } else{
      if ("Efficacy" %in% unique(dataStat_STUDY_i$TYPENAME)){
        PDstudy = TRUE
      }
    }


    # Loop over Treatment
    HeadStudy = TRUE
    for (j_TRT in 1:length(Treatment)){
      #Get Treatment Name:
      Treatment_j <- Treatment[j_TRT]

      # Save information into PKPDsummary:
      PKPDsummary$TRTNAME[n_TRT*(i_STUDY-1)+j_TRT] <- as.character(Treatment_j)
      PKPDsummary$Nsubj[n_TRT*(i_STUDY-1)+j_TRT] <- tableGR[STUDY_i][[1]][which(tableGR[[1]]==Treatment_j)]
      if (PKPDsummary$Nsubj[n_TRT*(i_STUDY-1)+j_TRT]!=0 & HeadStudy){
        PKPDsummary$STUDY[n_TRT*(i_STUDY-1)+j_TRT]     <- paste0(STUDY_i, " (", Center[[which(names(Center)==STUDY_i)]], ")")
        PKPDsummary$StudyTYPE[n_TRT*(i_STUDY-1)+j_TRT] <- StudyTYPE[[which(names(StudyTYPE)==STUDY_i)]]
        PKPDsummary$Analysis[n_TRT*(i_STUDY-1)+j_TRT]  <- ifelse((PKstudy & PDstudy),
                                                               "PK & PD",
                                                               ifelse(PKstudy,
                                                                      "PK",
                                                                      ifelse(PDstudy,
                                                                             "PD",
                                                                             "-")))
        HeadStudy = FALSE
      }
    }
  }

  #   Adjust last row:
  PKPDsummary$STUDY[n_ROW] = "Total"
  PKPDsummary$Nsubj[n_ROW] = sum(as.numeric(PKPDsummary$Nsubj[1:(n_ROW-2)]))

  #   Keep only columns and lines of interest
  PKPDsummary <- PKPDsummary[,c("STUDY", "StudyTYPE", "TRTNAME", "Nsubj", "Analysis")]
  PKPDsummary <- PKPDsummary[PKPDsummary$Nsubj!="0",]

  # Change Column Name:
  names(PKPDsummary) <- c("Study (Center)", "Type", "Treatment Groups", "N subjects", "Analysis")

  # Generate Extra columns to split the treatment group, if needed:
  if (!is.null(CovToTRTNAME)){
    # Dose Group Column:
    PKPDsummary$DOSE <- ""

    # Extra Column:
    for (k in 1:length(CovToTRTNAME)){
      PKPDsummary[,CovToTRTNAME[k]] <- ""
    }

    # Fill out columns
    for (TRT_k in PKPDsummary[,'Treatment Groups']){
      idx <- which((names(CovNewTRTNAME)==TRT_k))
      if (length(idx)!=0){
        PKPDsummary$DOSE[PKPDsummary[,'Treatment Groups']==TRT_k] <- CovNewTRTNAME[[idx]][1]
        for (i in 1:length(CovToTRTNAME)){
          PKPDsummary[PKPDsummary[,'Treatment Groups']==TRT_k,CovToTRTNAME[i]] <- CovNewTRTNAME[[idx]][i+1]
        }
      }
    }

    #   Re-organize table
    PKPDsummary        <- PKPDsummary[,c("Study (Center)", "Type", "DOSE", CovToTRTNAME, "N subjects", "Analysis")]
    if (!is.null(names(CovToTRTNAME))){
      names(PKPDsummary) <- c("Study (Center)", "Type", "Dose Groups", names(CovToTRTNAME), "N subjects", "Analysis")
    }else{
      names(PKPDsummary) <- c("Study (Center)", "Type", "Dose Groups", CovToTRTNAME, "N subjects", "Analysis")
    }
  }

  # Save Table:
  if (!is.null(filename)){
    IQRoutputTable(PKPDsummary,
                   xtitle   = "Summary information of the studies included in the dataset.",
                   xfooter  = "N: Number of",
                   filename = filename,
                   report   = TRUE)
  }

  # Output:
  return(PKPDsummary)
}

#' symmarybyGroup_MMVnca
#'
#' @description
#' @param resultNCA
#' @param stratify Default: NULL
#' @param stats Default: 'geomean'
#' @param colSummary Default: NULL
#' @param SIGNIF Default: 4
#' @param report Default: TRUE
#' 
#' @importFrom MMVbase GeometricMean
#' 
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
summarybyGroup_MMVnca <- function(resultNCA,
                                  stratify   = NULL,
                                  stats      = "geomean",
                                  colSummary = NULL,
                                  SIGNIF     = 4,
                                  report     = TRUE)
{

  # Check that resultNCA is a IQRnca object:
  if (!("IQRnca" %in% class(resultNCA))){
    stop("'resultNCA' should be an 'IQRnca' object.")
  }

  # Get attirbutes:
  attNCA <- attributes(resultNCA)

  # Get default colSummary if not specified:
  if (is.null(colSummary)){
    if ("CLO" %in% attNCA$names){
      colSummary <- c("AUCIFOD", "CLO" , "CMAXD", "LAMZHL", "TMAX", "VZO")

    }else if ("CLFO" %in% attNCA$names){
      colSummary <- c("AUCIFOD", "CLFO" , "CMAXD", "LAMZHL", "TMAX", "VZFO")

    }
  }

  # Check if VSS is desired and absent from resultNCA:
  #   Observed VSS - IV
  if (("VSSO" %in% colSummary) & !("VSSO" %in% attNCA$names) & ("CLO" %in% attNCA$names)){
    resultNCA$VSSO <- resultNCA$MRTIVIFO * resultNCA$CLO
    attNCA$units[["VSSO"]]         <- attNCA$units[["VZO"]]
    attNCA$IQRWNLmapping[["VSSO"]] <- "Vss_obs"
    attNCA$outputInfo[["VSSO"]]    <- "Observed Steady State Volume of Distribution"
  }
  #   Observed VSS - ORAL
  if (("VSSFO" %in% colSummary) & !("VSSFO" %in% attNCA$names) & ("CLFO" %in% attNCA$names)){
    resultNCA$VSSFO <- resultNCA$MRTEVIFO * resultNCA$CLFO
    attNCA$units[["VSSFO"]]         <- attNCA$units[["VZFO"]]
    attNCA$IQRWNLmapping[["VSSFO"]] <- "Vss_F_obs"
    attNCA$outputInfo[["VSSFO"]]    <- "Observed Steady State Volume of Distribution of Absorbed Fraction"
  }
  #   Predicted VSS - IV
  if (("VSSP" %in% colSummary) & !("VSSP" %in% attNCA$names) & ("CLP" %in% attNCA$names)){
    resultNCA$VSSP <- resultNCA$MRTIVIFP * resultNCA$CLP
    attNCA$units[["VSSP"]]         <- attNCA$units[["VZP"]]
    attNCA$IQRWNLmapping[["VSSP"]] <- "Vss_pred"
    attNCA$outputInfo[["VSSP"]]    <- "Predicted Steady State Volume of Distribution"
  }
  #   Predicted VSS - ORAL
  if (("VSSFP" %in% colSummary) & !("VSSFP" %in% attNCA$names) & ("CLFP" %in% attNCA$names)){
    resultNCA$VSSFP <- resultNCA$MRTEVIFP * resultNCA$CLFP
    attNCA$units[["VSSFP"]]         <- attNCA$units[["VZFP"]]
    attNCA$IQRWNLmapping[["VSSFP"]] <- "Vss_F_pred"
    attNCA$outputInfo[["VSSFP"]]    <- "Predicted Steady State Volume of Distribution of Absorbed Fraction"
  }

  # Choose the appropriate function to do the mean and sd/CV
  if(tolower(stats)=="geomean"){
    # Mean:
    FuncMean <- function(x){
      out <- MMVbase::GeometricMean(x, na.rm = TRUE)
    }
    # sd/CV:
    FuncSD <- function (x){
      out <- sqrt(exp(sd(log(x[!is.na(x)]))^2)-1)*100

    }
    # Text:
    FuncText <- function (){
      return("Values represent Geometric Mean (CV%)")
    }

  }else if(tolower(stats)=="mean"){
    # Mean:
    FuncMean <- function(x){
      out <- mean(x, na.rm = TRUE) # MMVbase::GeometricMean(x)
    }
    # sd/CV:
    FuncSD <- function (x){
      out <- sd(x, na.rm = TRUE)
    }
    # Text:
    FuncText <- function (){
      return("Values represent Mean (SD)")
    }

  }else{
    stop("'stats' should be 'mean' or 'geomean'.")
  }

  # Generate output:
  tabSumNCA <- list(xtable  = NULL,
                    xtitle  = "Comparison of NCA results",
                    xfooter = FuncText(),
                    report  = report)

  # Loop around colSummary:
  for (col_k in colSummary){

    # Check if it needs to be stratified:
    if (is.null(stratify)){
      # Summarize:
      tabSumNCA_k <- data.frame(N     = nrow(resultNCA),
                                Value = paste0(signif(FuncMean(resultNCA[,col_k]), SIGNIF), " (",
                                               signif(FuncSD(resultNCA[,col_k])  , SIGNIF), ")"),
                                stringsAsFactors = FALSE)

      # Adjust column name:
      names(tabSumNCA_k) <- c("N", paste0(attNCA$IQRWNLmapping[[col_k]], " [", attNCA$units[[col_k]], "]"))

    }else{
      # Summarize:
      tabSumNCA_k <- ddply(resultNCA, stratify, function(x){
        # Summarize:
        out <- data.frame(N     = nrow(x),
                          Value = paste0(signif(FuncMean(x[,col_k]), SIGNIF), " (",
                                         signif(FuncSD(x[,col_k]), SIGNIF), ")"),
                          stringsAsFactors = FALSE)

        # Adjust column name:
        names(out) <- c("N", paste0(attNCA$IQRWNLmapping[[col_k]], " [", attNCA$units[[col_k]], "]"))

        # Output:
        return(out)
      })
    }

    # Concatenate tabSumNCA_k to tabSumNCA$xtable:
    if (is.null(tabSumNCA$xtable)){
      tabSumNCA$xtable <- tabSumNCA_k
    }else{
      tabSumNCA$xtable <- left_join(tabSumNCA$xtable, tabSumNCA_k)
    }

    # Adjust footer:
    tabSumNCA$xfooter <- paste0(tabSumNCA$xfooter,
                                "\n", attNCA$IQRWNLmapping[[col_k]], ": ", attNCA$outputInfo[[col_k]])
  }

  # Adjust class of tabSumNCA:
  class(tabSumNCA) <- c("IQRoutputTable", class(tabSumNCA))

  # Ouput:
  return(tabSumNCA)
}

#' Create stratified barplot of ACPR
#'
#' @param data general dataset
#' @param type 'Crude' or 'PCR-Adjusted'
#' @param acprDay ACPR assessment day, can only be 28, 42, or 63
#' @param stratifyBy Column to stratify by
#'
#' @return
#' @export
#' @author Anne Kuemmel (IntiQuan)
#' @family DataExploration
barplotACPR <- function(data,
                        type       = "Crude",
                        acprDay    = 28,
                        stratifyBy = "BWB") {

  # Input check:
  if (!(type    %in% c("Crude", "PCR-Adjusted"))) stop("Wrong type input: Only 'Crude' or 'PCR-Adjusted' allowed.")
  if (!(acprDay %in% c(14, 21, 28, 42, 63)))      stop("Wrong day input: Only 14, 21, 28, 42 and 63 allowed.")

  # Get Categorical information:
  catInfo <- catInfo(data)

  # Stratification:
  if (!(stratifyBy %in% catInfo$COLNAME)) {
    warning("Stratification is not a categorical covariate column.")
    levels <- sort(unique(data[[stratifyBy]]))
    labels <- levels
    xtitle <- paste0("Column ", stratifyBy)
  } else {
    levels <- as.numeric(IQRtools::aux_explode(catInfo$VALUES[catInfo$COLNAME==stratifyBy]))
    labels <- IQRtools::aux_explode(catInfo$VALUETXT[catInfo$COLNAME==stratifyBy])
    xtitle <- catInfo$NAME[catInfo$COLNAME==stratifyBy]
  }

  # Determine which ACPR:
  varname <- paste0(case_when(type == "Crude"       ~ "c",
                              type == "PCR-Adjusted" ~ "pcr"),
                    "ACPR", acprDay)
  varsym  <- sym(varname)
  ytitle  <- paste0(type, " ACPR", acprDay)

  # Get Numeric and TXT value of Cure:
  VALUETXT <- strsplit(catInfo[catInfo$COLNAME==varsym,"VALUETXT"],split=",", fixed=TRUE)[[1]]
  VALUES   <- as.numeric(strsplit(catInfo[catInfo$COLNAME==varsym,"VALUES"],split=",", fixed=TRUE)[[1]])

  # Generate the data to analyze/plot:
  dataNumbers <- as.data.frame(data) %>%
    # Select one row per subject
    select(USUBJID, ID, one_of(c(stratifyBy, varname))) %>%
    unique() %>%
    # Calculate numbers subject in ACPR categories
    group_by_at(stratifyBy) %>%
    summarise(Ncure    = sum(!!varsym == VALUES[which(VALUETXT=="Cure")   ], na.rm = TRUE),
              Nfailure = sum(!!varsym == VALUES[which(VALUETXT=="Failure")], na.rm = TRUE),
              Nmissing = sum(!!varsym == VALUES[which(VALUETXT=="Missing")], na.rm = TRUE)) %>%
    ungroup()

  # Summary of ACPR by stratifier:
  dataSummary <- dataNumbers %>%
    group_by_at(stratifyBy) %>%
    tidyr::nest() %>%
    mutate(ACPRstats = purrr::map(data, ~switch(type,
                                                "Crude"        = Hmisc::binconf(.$Ncure, .$Ncure+.$Nfailure),
                                                "PCR-Adjusted" = Hmisc::binconf(.$Ncure, .$Ncure+.$Nfailure)))) %>%
    mutate(ACPR   = purrr::map_dbl(ACPRstats, ~.[1]),
           ACPRlb = purrr::map_dbl(ACPRstats, ~.[2]),
           ACPRub = purrr::map_dbl(ACPRstats, ~.[3])) %>%
    select(-ACPRstats) %>%
    tidyr::unnest() %>%
    mutate(N = switch(type,
                      "Crude"        = Ncure+Nfailure,
                      "PCR-Adjusted" = Ncure+Nfailure))

  # Nicer x mark text:
  stratsym    <- sym(stratifyBy)
  dataSummary <- mutate(dataSummary, STRAT = factor(!!stratsym, levels = levels, labels = labels))

  # Bar plot
  myplot <- dataSummary %>%
    MMVggplot(aes(STRAT, y = ACPR, ymin = ACPRlb, ymax = ACPRub)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_errorbar(width = 0.3) +
    geom_text(aes(y=0,label = paste0("N=",N)), vjust=1.2, size = 4) +
    labs(x=xtitle, y=ytitle, caption = "Errorbar indicates 95% confidence interval.") +
    scale_y_continuous(labels = scales::percent)

  # Create table
  mytable <- dataSummary
  names(mytable) <- gsub("STRAT", xtitle,  names(mytable))
  mytable <- select(mytable, -!!stratsym) %>%
    select(one_of(xtitle), N, Ncure, Nfailure, Nmissing, ACPR, ACPRlb, ACPRub) %>%
    mutate('ACPR (%)' = signif(100*ACPR,3),
           'ACPR 95%-CI' = paste0(signif(ACPRlb*100, 3), " - ", signif(ACPRub*100,3))) %>%
    select(-ACPR,-ACPRlb, -ACPRub)

  mytable <- IQRoutputTable(mytable,
                            xtitle = paste0(type, "ACPR stratified by ", xtitle),
                            xfooter = "Values rounded to 3 significant digits.")

  # Output:
  return(list(plot=myplot, table=mytable))
}

#' Scatterplot concentrationand ACPR
#'
#' @param data general data set
#' @param conc NAME of the concentration to plot
#' @param concDay VISNAME of the timepoint to take the concentrations from
#' @param type 'Crude' or 'PCR-Adjusted'
#' @param acprDay ACPR assessment day, can only be 28, 42, or 63
#'
#' @return
#' @export
#' @author Anne Kuemmel (IntiQuan)
#' @family DataExploration
scatterConcACPR <- function(data,
                            conc,
                            concDay = "Day 7",
                            type    = "Crude",
                            acprDay = 28) {

  # Input check:
  if (!(type     %in% c("Crude", "PCR-Adjusted"))) stop("Wrong type input: only 'Crude' or 'PCR-Adjusted' allowed.")
  if (!(acprDay %in% c(14, 21, 28, 42, 63)))       stop("Wrong day input: Only 14, 21, 28, 42 and 63 allowed.")
  if (!(conc     %in% data$NAME))                  stop("Given conc does not exist as NAME in data set.")

  # Day Name:
  if (is.numeric(concDay)){
    concDayName <- paste0("Day ", concDay)
  }else if(is.character(concDay)){
    concDayName <- concDay
  }else{
    stop("'concDay' should be numeric or a character corresponding to a 'VISNAME'")
  }

  # Get Categorical information:
  catInfo <- catInfo(data)

  # Determine which ACPR:
  varname <- paste0(case_when(type == "Crude" ~ "c",
                              type == "PCR-Adjusted" ~ "pcr"),
                    "ACPR", acprDay)
  varsym  <- sym(varname)
  ytitle  <- paste0(type, " ACPR", acprDay)

  # Get Numeric and TXT value of Cure:
  VALUETXT <- strsplit(catInfo[catInfo$COLNAME==varsym,"VALUETXT"],split=",", fixed=TRUE)[[1]]
  VALUES   <- as.numeric(strsplit(catInfo[catInfo$COLNAME==varsym,"VALUES"],split=",", fixed=TRUE)[[1]])

  # Select ACPR categories depending on ACPR type:
  if (type == "Crude")        selectedCat <- VALUES[which(VALUETXT %in% c("Cure", "Failure"))]
  if (type == "PCR-Adjusted") selectedCat <- VALUES[which(VALUETXT %in% c("Cure", "Failure"))]

  # Generate the data to analyze/plot:
  dataPlot <- as.data.frame(data) %>%
    # Select only Day7 concentrations
    filter(NAME==conc & (if (is.numeric(concDay)) NT==concDay*24 else VISNAME==concDay)) %>%
    # Filter ACPR categories
    filter(!!varsym %in% selectedCat) %>%
    mutate(ACPR = as.numeric(!!varsym==1))

  # x axis label:
  xtitle <- paste0(conc, " ", concDayName, " (",data$UNIT[1],")")

  # Scatter plot:
  myplot <- dataPlot %>%
    MMVggplot(aes(x = DV, y = ACPR)) +
    geom_point(shape = 1) +
    geom_smooth(method = "loess") +
    labs(x=xtitle, y=ytitle) +
    scale_x_log10() +
    scale_y_continuous(breaks = c(0,1)) +
    coord_cartesian(ylim = c(0,1))

  # Output:
  return(myplot)
}


#' plot_2DconcVsACPR
#'
#' @param data general data set
#' @param conc1 NAME of the concentration of Drug 1
#' @param conc2 NAME of the concentration of Drug 2
#' @param concDay VISNAME of the timepoint to take the concentrations from
#' @param type 'Crude' or 'PCR-Adjusted'
#' @param acprDay ACPR assessment day, can only be 28, 42, or 63
#'
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family DataExploration
plot_2DconcVsACPR <- function(data,
                              conc1,
                              conc2,
                              concDay = "Day 7",
                              type    = "Crude",
                              acprDay = 28) {

  # Input check:
  if (!(type    %in% c("Crude", "PCR-Adjusted"))){stop("Wrong type input: only 'Crude' or 'PCR-Adjusted' allowed.")}
  if (!(acprDay %in% c(14, 21, 28, 42, 63)))     {stop("Wrong day input: Only 14, 21, 28, 42 and 63 allowed.")}
  if (!(conc1   %in% data$NAME))                 {stop("Given conc does not exist as NAME in data set.")}
  if (!(conc2   %in% data$NAME))                 {stop("Given conc does not exist as NAME in data set.")}

  # Day Name:
  if (is.numeric(concDay)){
    concDayName <- paste0("Day ", concDay)
  }else if(is.character(concDay)){
    concDayName <- concDay
  }else{
    stop("'concDay' should be numeric or a character corresponding to a 'VISNAME'")
  }

  # Get Categorical information:
  catInfo <- catInfo(data)

  # Determine which ACPR:
  varname <- paste0(case_when(type == "Crude" ~ "c",
                              type == "PCR-Adjusted" ~ "pcr"),
                    "ACPR", acprDay)
  varsym  <- sym(varname)
  title  <- paste0(type, " ACPR", acprDay)

  # Get Numeric and TXT value of Cure:
  VALUETXT <- strsplit(catInfo[catInfo$COLNAME==varsym,"VALUETXT"],split=",", fixed=TRUE)[[1]]
  VALUES   <- as.numeric(strsplit(catInfo[catInfo$COLNAME==varsym,"VALUES"],split=",", fixed=TRUE)[[1]])

  # Select ACPR categories depending on ACPR type:
  if (type == "Crude")        {selectedCat <- VALUES[which(VALUETXT %in% c("Cure", "Failure"))]}
  if (type == "PCR-Adjusted") {selectedCat <- VALUES[which(VALUETXT %in% c("Cure", "Failure"))]}

  # Generate the data to analyze/plot:
  dataPlot <- as.data.frame(data) %>%
    select(USUBJID, VISNAME, NT, TIME, TYPENAME, NAME, VALUE, VALUETXT, LLOQ, !!varsym) %>%
    # Select only Day7 concentrations
    filter((NAME %in% c(conc1,conc2)) & (if (is.numeric(concDay)) NT==concDay*24 else VISNAME==concDay)) %>%
    # Filter ACPR categories
    filter(!!varsym %in% selectedCat)

  # Long to Wide:
  dataPlot <- transform_dataFrame_LongToWide(dataPlot,
                                             key   = c("NAME"),
                                             value = c("VALUE", "LLOQ"),
                                             colMaster = c("USUBJID", "VISNAME", "NT", "TIME", "TYPENAME", as.character(varsym)))

  # Rename column:
  names(dataPlot) <- c("USUBJID", "VISNAME", "NT", "TIME", "TYPENAME", "ACPR", "VALUE1", "LLOQ1", "VALUE2", "LLOQ2")

  # Get LLOQ:
  LLOQ1 <- unique(dataPlot$LLOQ1[!is.na(dataPlot$LLOQ1)])[1]
  LLOQ2 <- unique(dataPlot$LLOQ2[!is.na(dataPlot$LLOQ2)])[1]

  # Adjust na value or below LLOQ:
  dataPlot$VALUE1 <- ifelse(is.na(dataPlot$VALUE1) | dataPlot$VALUE1<=LLOQ1,
                            LLOQ1,
                            dataPlot$VALUE1)
  dataPlot$VALUE2 <- ifelse(is.na(dataPlot$VALUE2) | dataPlot$VALUE2<=LLOQ2,
                            LLOQ2,
                            dataPlot$VALUE2)

  # Rename ACPR:
  dataPlot$ACPR <- VALUETXT[match(dataPlot$ACPR,VALUES)]

  # Get Median Concentration:
  C1median <- quantile(dataPlot$VALUE1, 0.50)
  C2median <- quantile(dataPlot$VALUE2, 0.50)

  # Some stats:
  dataStat <- data.frame(Corner     = c("BottomLeft", "TopLeft", "BottomRight", "TopRight"),
                         NbrSucces  = c(length(dataPlot$ACPR[dataPlot$VALUE1< C1median & dataPlot$VALUE2< C2median & dataPlot$ACPR=="Cure"]),
                                        length(dataPlot$ACPR[dataPlot$VALUE1< C1median & dataPlot$VALUE2>=C2median & dataPlot$ACPR=="Cure"]),
                                        length(dataPlot$ACPR[dataPlot$VALUE1>=C1median & dataPlot$VALUE2< C2median & dataPlot$ACPR=="Cure"]),
                                        length(dataPlot$ACPR[dataPlot$VALUE1>=C1median & dataPlot$VALUE2>=C2median & dataPlot$ACPR=="Cure"])),
                         NbrSubject = c(length(dataPlot$ACPR[dataPlot$VALUE1< C1median & dataPlot$VALUE2< C2median]),
                                        length(dataPlot$ACPR[dataPlot$VALUE1< C1median & dataPlot$VALUE2>=C2median]),
                                        length(dataPlot$ACPR[dataPlot$VALUE1>=C1median & dataPlot$VALUE2< C2median]),
                                        length(dataPlot$ACPR[dataPlot$VALUE1>=C1median & dataPlot$VALUE2>=C2median])),
                         stringsAsFactors = FALSE)
  dataStat <- within(dataStat,{
    Sucess <- paste0("ACPR=", NbrSucces, "/", NbrSubject, "=", signif(NbrSucces/NbrSubject*100,2), "%")
  })

  # Scatter plot:
  myplot <- dataPlot %>%
    MMVggplot(aes(x = VALUE1*1000, y = VALUE2*1000, color = ACPR)) +
    #   LLOQ lines
    geom_hline(yintercept=LLOQ2*1000, linetype=4, size=1) +
    geom_text(x = log10(max(dataPlot$VALUE1)*1000), y = log10(LLOQ2*1000), label=paste0("LLOQ=",LLOQ2*1000,"ng/mL"), color="black", hjust=1, vjust=1.2, size=3) +
    geom_vline(xintercept=LLOQ1*1000, linetype=4, size=1) +
    geom_text(x = log10(LLOQ1*1000), y = log10(max(dataPlot$VALUE2)*1000), label=paste0("LLOQ=",LLOQ1*1000,"ng/mL"), color="black", hjust=0, vjust=1.2, size=3, angle = -90) +
    #   Data
    geom_point(shape = 16, size=2.5) +
    scale_color_manual("ACPR", values = MMVcolors[3:2]) +
    labs(title = title,
         x=paste0(conc1, " ", concDayName, " [ng/mL]"),
         y=paste0(conc2, " ", concDayName, " [ng/mL]")) +
    scale_x_log10() +
    scale_y_log10() +
    #   Median Concentration
    geom_hline(yintercept=C2median*1000, linetype=1, size=1.25) +
    geom_text(x = log10(LLOQ1*1000)   , y = log10(C2median*1000)            , label="Median Concentration", color="black", hjust=-0.05, vjust=-0.4, size=3) +
    geom_vline(xintercept=C1median*1000, linetype=1, size=1.25) +
    geom_text(x = log10(C1median*1000), y = log10(max(dataPlot$VALUE2)*1000), label="Median Concentration", color="black", hjust= 0 , vjust=-0.5, size=3, angle = -90) +
    #   Success by Corner
    geom_text(label=dataStat$Sucess[dataStat$Corner=="BottomLeft"] , x = log10(LLOQ1*1000)               , y = log10(LLOQ2*1000)               , color="navyblue", hjust=-0.05, vjust=-0.4, size=4) +
    geom_text(label=dataStat$Sucess[dataStat$Corner=="TopLeft"]    , x = log10(LLOQ1*1000)               , y = log10(max(dataPlot$VALUE2)*1000), color="navyblue", hjust=-0.05, vjust=-0.4, size=4) +
    geom_text(label=dataStat$Sucess[dataStat$Corner=="BottomRight"], x = log10(max(dataPlot$VALUE1)*1000), y = log10(LLOQ2*1000)               , color="navyblue", hjust=1    , vjust=-0.4, size=4) +
    geom_text(label=dataStat$Sucess[dataStat$Corner=="TopRight"]   , x = log10(max(dataPlot$VALUE1)*1000), y = log10(max(dataPlot$VALUE2)*1000), color="navyblue", hjust=1    , vjust=-0.4, size=4)

  # Output:
  return(myplot)
}


#' Plot PKPD data available from a Human Challenge study in "standard" MMV format, suitable for exploration and reporting
#'
#' @param dataGen : Object containing data to plot in `IQRDataGeneral` format
#' @param filePath : Filepath to which plots will be saved
#' @param PDname : Name of PD observations in the NAME column of IQRDataGeneral. Default: 'Parasitemia'
#' @param Gametos : Name of Gametocyte observations in the NAME column of IQRDataGeneral. Default: 'Parasitemia Gametocytes'
#' @param rescuename : Name of observations in the NAME column that correspond to Rescue medication doses.
#' @param PercentPD : Proportion of total parasitaemia at which to adjust naming of datapoints i.e. for gametocytes: Default: 0.1
#'
#' @return Figures as .png to destination specified in `filePath`
#'
#' @export
#' @md
#'
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV), Karsten Kuritz (IntiQuan)
#' @family DataExploration
#' @importFrom plyr join dlply
plot_PKPDdataMMVhuCh <- function(dataGen,
                                     filePath,
                                     PDname  = "Parasitemia",
                                     Gametos = "Parasitemia Gametocytes",
                                     rescuename = "Riamet",
                                     PercentPD = 0.1) {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle input arguments ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (length(unique(dataGen$NAME[dataGen$TYPENAME == "Dose"])) > 1)
    stop("Function only handles monotherapy data. Two types of doses found.")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load dataset ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dataPlot       <- as.data.frame(dataGen)
  dataPlot$Label <- with(dataPlot, paste0(NAME, " (", UNIT,")"))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get nicer treatment name ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TRTann <- within(unique(dataPlot[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)
    DOSEMULT  <- ifelse(DOSELEVEL == 0, 0, DOSEMULT)
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]

  dataPlot <- plyr::join(dataPlot, TRTann[,c("TRTNAME","Treatment")])

  doseunit <- dataPlot[dataPlot$NAME %in% attr(dataGen, "doseNAMES"),"UNIT"][1]


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot parasites total and gametocytes if given ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Patients with PD data of interest
  PatientsOfInterest <- unique(dataPlot[dataPlot$NAME %in% PDname,"USUBJID"])

  # Reduce dataset to observation of interest:
  drugDose     <- unique(dataGen$NAME[dataGen$TYPENAME == "Dose"])
  rescueDose   <- grep(rescuename,doseNAMES(dataGen), value = TRUE)
  dataPlotPDPD <- dataPlot[dataPlot$NAME %in% c(PDname, Gametos, drugDose, rescueDose), ]
  dataPlotPDPD <- dataPlotPDPD[dataPlotPDPD$USUBJID %in% PatientsOfInterest, ]


  # Get LLOQ information
  lloq <- unique(subset(dataPlotPDPD, TYPENAME == "Efficacy", c("Label","LLOQ","STUDY")))

  # Plots: dlply instead of a for loop
  gr <- plyr::dlply(dataPlotPDPD, ~Treatment, function(x) {
    xP <- x[x$TYPENAME == "Efficacy",]
    xPD <- x[x$NAME == PDname,]
    xPD$VALUE <- xPD$VALUE*PercentPD
    xPD$NAME <- paste(round(100*PercentPD), "% of", PDname)
    xPx <- rbind(xP, xPD)
    xPc<- x[x$TYPENAME == "Efficacy" & x$MDV == 1,]
    xD <- x[x$TYPENAME == "Dose",]
    xR <- x[x$TYPENAME == "Rescue",]
    # Keep only the first time for each:
    xR_temp = NULL
    for (SUBJECT_k in unique(xR$SUBJECT)){
      # Sub-dataset:
      xR_k <- xR[xR$SUBJECT==SUBJECT_k,]

      # Index of the minimum:
      k = which.min(xR_k$TIME)

      # Add to dataFrame:
      xR_temp <- rbind(xR_temp, xR_k[k,])
    }
    xR <- xR_temp

    # Create Plot:
    MMVggplot(xPx, aes(TIME, VALUE)) +
      geom_vline(data=xD, aes(xintercept = TIME), linetype = 2) +
      geom_vline(data=xR, aes(xintercept = TIME), linetype = 4, color = "purple") +
      geom_line(aes(group = interaction(USUBJID,NAME)), color = "grey") +
      geom_point(aes(color = NAME)) +
      geom_point(data=xPc, shape = 4) +
      scale_color_manual("Parasite type",values=MMVcolors[2:10]) +
      scale_y_log10() +
      scale_x_continuous(breaks = seq(-120,1000,48)) +
      facet_wrap(~SUBJECT) +
      labs(
        x = "Time (hours)",
        y = "Parasite Counts (p/mL)",
        title = paste0("Individual parasite data for ", x$Treatment[1]),
        caption = "Experimental drug and recue medication doses indicated by black and purple lines, respectively."
      ) +
      theme(legend.position = "bottom", legend.direction = "vertical")
  })

  # Save Plots:
  for (k in 1:length(gr))
    IQRoutputPNG(gr[[k]], filename = file.path(filePath, sprintf("08-%02.0f_GameTotPDlinePlot.png", k)))
}

#' Plot PKPD data available from a Human Challenge study in "standard" MMV format, suitable for exploration and reporting.
#' plots ratio of Individual parasite data over gametocyte counts, stratified by treatment.
#'
#' @param dataGen : Object containing data to plot in `IQRDataGeneral` format
#' @param filePath : Filepath to which plots will be saved
#' @param PDname : Name of PD observations in the NAME column of IQRDataGeneral. Default: 'Parasitemia'
#' @param Gametos : Name of Gametocyte observations in the NAME column of IQRDataGeneral. Default: 'Parasitemia Gametocytes'
#' @param rescuename : Name of observations in the NAME column that correspond to Rescue medication doses.
#'
#' @return
#' @export
#' @md
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV), Karsten Kuritz #' (IntiQuan)
#' @family DataExploration
#' @importFrom plyr join dlply
plot_PKPDdataMMVhuCh_Ratio_Gam <- function(dataGen,
                                           filePath,
                                           PDname  = "Parasitemia",
                                           Gametos = "Parasitemia Gametocytes",
                                           rescuename = "Riamet") {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle input arguments ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (length(unique(dataGen$NAME[dataGen$TYPENAME == "Dose"])) > 1)
    stop("Function only handles monotherapy data. Two types of doses found.")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load dataset ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dataPlot       <- as.data.frame(dataGen)
  dataPlot$Label <- with(dataPlot, paste0(NAME, " (", UNIT,")"))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get nicer treatment name ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TRTann <- within(unique(dataPlot[,c("TRTNAME", "DOSELEVEL", "DOSEMULT")]), {
    Treatment <- TRTNAME
    Treatment <- ifelse(DOSELEVEL == 0, "Vehicle", Treatment)
    DOSEMULT  <- ifelse(DOSELEVEL == 0, 0, DOSEMULT)
  })
  TRTann <- TRTann[with(TRTann, order(DOSEMULT, DOSELEVEL)),]

  dataPlot <- plyr::join(dataPlot, TRTann[,c("TRTNAME","Treatment")])

  doseunit <- dataPlot[dataPlot$NAME %in% attr(dataGen, "doseNAMES"),"UNIT"][1]


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot parasites total and gametocytes if given ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Patients with PD data of interest
  PatientsOfInterest <- unique(dataPlot[dataPlot$NAME %in% PDname,"USUBJID"])

  # Reduce dataset to observation of interest:
  drugDose     <- unique(dataGen$NAME[dataGen$TYPENAME == "Dose"])
  rescueDose   <- grep(rescuename,doseNAMES(dataGen), value = TRUE)
  dataPlotPDPD <- dataPlot[dataPlot$NAME %in% c(PDname, Gametos, drugDose, rescueDose), ]
  dataPlotPDPD <- dataPlotPDPD[dataPlotPDPD$USUBJID %in% PatientsOfInterest, ]


  # Get LLOQ information
  lloq <- unique(subset(dataPlotPDPD, TYPENAME == "Efficacy", c("Label","LLOQ","STUDY")))

  # Plots: dlply instead of a for loop
  gr <- plyr::dlply(dataPlotPDPD, ~Treatment, function(x) {
    xP <- x[x$TYPENAME == "Efficacy",]
    xRatio <- xP %>% tidyr::pivot_wider(values_from = "VALUE", names_from = "NAME")
    xRatio$VALUE <- xRatio[[Gametos]]/xRatio[[PDname]]
    xRatio$NAME <- paste(Gametos, "over", PDname)
    xPc<- x[x$TYPENAME == "Efficacy" & x$MDV == 1,]
    xD <- x[x$TYPENAME == "Dose",]
    xR <- x[x$TYPENAME == "Rescue",]

    # Keep only the first time for each:
    xR_temp = NULL
    for (SUBJECT_k in unique(xR$SUBJECT)){
      # Sub-dataset:
      xR_k <- xR[xR$SUBJECT==SUBJECT_k,]

      # Index of the minimum:
      k = which.min(xR_k$TIME)

      # Add to dataFrame:
      xR_temp <- rbind(xR_temp, xR_k[k,])
    }
    xR <- xR_temp

    # Create Plot:
    MMVggplot(xRatio, aes(TIME, VALUE)) +
      geom_vline(data=xD, aes(xintercept = TIME), linetype = 2) +
      geom_vline(data=xR, aes(xintercept = TIME), linetype = 4, color = "purple") +
      geom_line(aes(group = interaction(USUBJID,NAME)), color = "grey") +
      geom_point(aes(color = NAME)) +
      geom_point(data=xPc, shape = 4) +
      scale_color_manual("Parasite type",values=MMVcolors[2:10]) +
      scale_y_log10() +
      scale_x_continuous(breaks = seq(-120,1000,48)) +
      facet_wrap(~SUBJECT) +
      labs(
        x = "Time (hours)",
        y = "Ratio Parasite Counts (p/mL) / Gametocytes",
        title = paste0("Ratio of Individual parasite data over gametocyte counts for ", x$Treatment[1]),
        caption = "Experimental drug and recue medication doses indicated by black and purple lines, respectively."
      ) +
      theme(legend.position = "bottom")
  })

  # Save Plots:
  for (k in 1:length(gr))
    IQRoutputPNG(gr[[k]], filename = file.path(filePath, sprintf("08-%02.0f_RatioGameTotPDlinePlot.png", k)))
}
