
#' censorGametocyte
#'
#' @description
#' @param data
#' @param threshold
#' @param nameAll Default: 'Parasitemia asexual + gametocyte'
#' @param nameGam Default: 'Parasitemia female gametocyte'
#' @param FLAGfollowing Default: TRUE
#' @return
#' @export
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
#' @family Data Preparation
censorGametocyte <- function(
  data,
  threshold,                                    # Percentage of gametocytes of total. If this is exceeded total will be censored.
  nameAll = "Parasitemia asexual + gametocyte",
  nameGam = "Parasitemia female gametocyte",
  FLAGfollowing = TRUE
) {

  cat(paste0(data$USUBJID[1], ": "))
  # Check whether both readouts are in dataset, if not do nothing.
  if ( !any(data$NAME == nameAll) | !any(data$NAME == nameGam) ) {cat("Nothing to do.\n"); return()}

  # Set BLOQ values to 0
  dataPara <- as.data.frame(data)[data$NAME %in% c(nameAll, nameGam),c("USUBJID", "TIME", "VALUE", "NAME", "LLOQ")]
  dataPara$VALUE <- with(dataPara, ifelse(VALUE<LLOQ, 0, VALUE))
  dataPara$LLOQ  <- NULL

  # Reshape to align total and gametocyte counts
  dataPara <- tidyr::spread(
    dataPara,
    key = "NAME",
    value = "VALUE"
  )

  # Calculate percent gametocytes
  dataPara$GametosPerc <- dataPara[[nameGam]] / dataPara[[nameAll]] * 100

  # Handle BLQ values which were set to 0:
  # - if both 0, NA is produced and record discarded as not useful for censoring
  # - if gametocytes below LLOQ, but total counts not, percent is zero and therefore not censoring
  # - if gametocytes are above LLOQ, but total counts below, NA would be produced but it is assumed that gametocyte levels are high and percent set to 100
  dataPara$GametosPerc <- ifelse(dataPara[[nameGam]] > 0 & dataPara[[nameAll]] == 0, 100, dataPara$GametosPerc)

  # Discard NAs as they cannot be used for censoring
  dataPara <- dataPara[!is.na(dataPara$GametosPerc),]

  # Make sure it is ordered by time
  dataPara <- dataPara[order(dataPara$TIME),]

  # Get timepoints for which gametocytes above threshold
  idxCensT <- dataPara$GametosPerc >= threshold
  if (all(!idxCensT)) {
    # case there is no gametocyte level above threshold
    cat(paste0("No gametocyte levels above threshold.\n"))
    return(NULL)
  } else {
    if (FLAGfollowing) {
      TimeCensT <- dataPara$TIME[min(which(idxCensT))]
      IXGDFcens <- data$IXGDF[data$TIME >= TimeCensT & data$NAME == nameAll]
      cat(paste0("Censoring from ",TimeCensT, " hours on.\n"))
    } else {
      IXGDFcens <- data$IXGDF[data$TIME %in% dataPara$TIME[idxCensT] & data$NAME == nameAll]
      cat(paste0("Censoring of ",sum(idxCensT), " timepoints.\n"))
    }
    return((IXGDFcens))
  }
}

#' Convert Gametocytes
#'
#' Converts gametocyte/parasitemia variable units from laboratory (generally copies/mL) to IQRdataset format units (p/mL).
#' Examples of gametocyte/parasitemia variables that might need unit conversion are Parasitemia female gametocytes, Parasitemia male gametocytes and Parasitemia Trophozoite.
#' Conversion factor (GamFactor) from copies/mL to p/mL will depend on the gametocyte/parasitemia variable and the laboratory.
#' oldStdCurve and oldNewStdFactor might be needed with old data from QIMR laboratory.
#'
#' @param data Data frame or IQRdataGENERAL object with gametocyte data.
#' @param gamName Variable name to convert (Default: `"Parasitemia female gametocytes"`).
#' @param convUnit Conversion factor from copies/ul to copies/mL (Default: 1).
#' @param GamFactor Conversion factor from copies/mL to p/mL (p being e.g. female gametocytes).
#' @param oldStdCurve TRUE or FALSE indicating whether old standard curve was used (Default: `"FALSE"`).
#' @param oldNewStdFactor Conversion factor from old to new standard curve (Default: 62).
#'
#' @return
#'
#' @examples
#' dat <- data.frame(NAME  = c( rep("Parasitemia female gametocytes", 4)),
#'                   TIME  = c( 100,  120,  200, 240),
#'                   VALUE = c(9400, 1500, 2100, 3700),
#'                   UNIT  = c(rep("copies/mL", 4)),
#'                   LLOQ  = c(rep(1089, 4)),
#'                   stringsAsFactors = FALSE)
#'
#' dat <- convert_Gametocytes(data            = dat,
#'                            gamName         = "Parasitemia female gametocytes",
#'                            GamFactor       = 279.3,
#'                            convUnit        = 1,
#'                            oldStdCurve     = FALSE
#' )
#'
#' @export
#' @author Aline Fuchs (MMV)
#' @family Data Preparation
convert_Gametocytes <- function(
  data,
  gamName      = "Parasitemia female gametocytes",
  GamFactor,
  convUnit        = 1,
  oldStdCurve     = FALSE,
  oldNewStdFactor = 62
) {

  # Initialize factor by which to convert current gametocyte data to p/mL in whole blood:
  Fact <- 1

  # ~~~~~~~~~~~~~~~~~~~~
  # Factors for converting data to copies/mL in whole blood:

  # Convert from old to new standard curve if old was used:
  if (oldStdCurve) {
    Fact <- Fact / oldNewStdFactor
  }

  # ~~~~~~~~~~~~~~~~~~~~
  # Conversion from copies to parasites (meaning e.g. female gametocytes):
  Fact <- Fact * convUnit/ GamFactor

  # ~~~~~~~~~~~~~~~~~~~~
  # Do conversion:
  idxGam <- data$NAME == gamName
  data$VALUE[idxGam] <- data$VALUE[idxGam] * Fact
  data$LLOQ[idxGam]  <- data$LLOQ[idxGam] * Fact
  data$UNIT[idxGam]  <- "p/mL"

  # Output:
  data
}



#' import_SCIDpkpdData
#'
#' @description
#' @param dataFile
#' @param compound
#' @param PKsheet Default: NULL
#' @param PKrange Default: NULL
#' @param PKmapping
#' @param PKlloq Default: 1
#' @param PKlloqIdentifier Default: c("BLQ", "<LLOQ", "< LLOQ")
#' @param PKfactor Default: 0.001
#' @param PDsheet Default: NULL
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01", "< 0,01", "< 0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param study Default: NULL
#' @param DayOfFirstDrugAdmin Default: NULL
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @param nrdoses Default: NULL
#' @param intervaldose Default: NULL
#' @param route Default: NULL
#' @param DoseTimeSheet Default: 'DosingTimeTable'
#' @param DoseTimeRange Default: NULL
#' @param DoseTimeMapping
#' @param DoseSheet Default: 'IndividualMetadataDrug1'
#' @param DoseRange Default: NULL
#' @param DoseMapping Default: NULL
#' @param positiveQC
#' @param importVersion Default: centerName
#' @return
#' @export
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Catalina Barcelo (MMV), Mohammed H. Cherkaoui (MMV), \email{cherkaouim@@mmv.org}
#' @family Data Preparation
import_SCIDpkpdData <- function(dataFile,
                                compound,
                                PKsheet             = NULL,
                                PKrange             = NULL,
                                PKmapping,
                                PKlloq              = 1,
                                PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ"),
                                PKfactor            = 1e-3,  # From ng/mL to ug/mL
                                PDsheet             = NULL,
                                PDrange             = NULL,
                                PDmapping,
                                PDlloq              = 0.01,
                                PDlloqIdentifier    = c("BLQ", "<0.01", "< 0,01", "< 0.01"),
                                PDname              = "Parasitemia Total",
                                study               = NULL,
                                DayOfFirstDrugAdmin = NULL,
                                centerNumber        = -1,
                                centerName          = "",
                                visitNumber         = -1,
                                nrdoses             = NULL,
                                intervaldose        = NULL,
                                route               = NULL,
                                DoseTimeSheet       = "DosingTimeTable",
                                DoseTimeRange       = NULL,
                                DoseTimeMapping,
                                DoseSheet           = "IndividualMetadataDrug1",
                                DoseRange           = NULL,
                                DoseMapping         = NULL,
                                positiveQC,
                                importVersion       = centerName
) {

  # Use the appropriate import function according to the center name:
  if (importVersion=="TAD"){
    # PK sheet:
    if (is.null(PKsheet)){
      PKsheet <- "DrugConcentrationData"
    }

    # PD sheet:
    if (is.null(PDsheet)){
      PDsheet <- "ParasitemiaData"
    }

    # DayOfFirstDrugAdmin
    if (is.null(DayOfFirstDrugAdmin)){
      DayOfFirstDrugAdmin <- 1
    }

    # Import Data:
    data <- import_SCIDpkpdData_TAD(dataFile            = dataFile,
                                    compound            = compound,
                                    PKsheet             = PKsheet,
                                    PKrange             = PKrange,
                                    PKmapping           = PKmapping,
                                    PKlloq              = PKlloq,
                                    PKlloqIdentifier    = PKlloqIdentifier,
                                    PKfactor            = PKfactor,
                                    PDsheet             = PDsheet,
                                    PDrange             = PDrange,
                                    PDmapping           = PDmapping,
                                    PDlloq              = PDlloq,
                                    PDlloqIdentifier    = PDlloqIdentifier,
                                    PDname              = PDname,
                                    DayOfFirstDrugAdmin = DayOfFirstDrugAdmin,
                                    centerNumber        = centerNumber,
                                    centerName          = centerName,
                                    visitNumber         = visitNumber,
                                    nrdoses             = nrdoses,
                                    intervaldose        = intervaldose,
                                    route               = route,
                                    DoseTimeSheet       = DoseTimeSheet,
                                    DoseTimeRange       = DoseTimeRange,
                                    DoseTimeMapping     = DoseTimeMapping,
                                    DoseSheet           = DoseSheet,
                                    DoseRange           = DoseRange,
                                    DoseMapping         = DoseMapping,
                                    positiveQC          = positiveQC)

  }else if(importVersion=="oldTAD"){
    # PK sheet:
    if (is.null(PKsheet)){
      PKsheet <- "BloodLevels"
    }

    # PD sheet:
    if (is.null(PDsheet)){
      PDsheet <- "Parasitemia"
    }

    # DayOfFirstDrugAdmin
    if (is.null(DayOfFirstDrugAdmin)){
      DayOfFirstDrugAdmin <- 1
    }

    # Import Data:
    data <- import_SCIDpkpdData_oldTAD(dataFile            = dataFile,
                                       compound            = compound,
                                       PKsheet             = PKsheet,
                                       PKrange             = PKrange,
                                       PKmapping           = PKmapping,
                                       PKlloq              = PKlloq,
                                       PKlloqIdentifier    = PKlloqIdentifier,
                                       PKfactor            = PKfactor,
                                       PDsheet             = PDsheet,
                                       PDrange             = PDrange,
                                       PDmapping           = PDmapping,
                                       PDlloq              = PDlloq,
                                       PDlloqIdentifier    = PDlloqIdentifier,
                                       PDname              = PDname,
                                       study               = study,
                                       DayOfFirstDrugAdmin = DayOfFirstDrugAdmin,
                                       centerNumber        = centerNumber,
                                       centerName          = centerName,
                                       visitNumber         = visitNumber,
                                       nrdoses             = nrdoses,
                                       intervaldose        = intervaldose,
                                       positiveQC          = positiveQC)

  }else if(importVersion=="GSK-TresCantos"){
    # PK sheet:
    if (is.null(PKsheet)){
      PKsheet <- "BloodLevels"
    }

    # PD sheet:
    if (is.null(PDsheet)){
      PDsheet <- "Parasitemia"
    }

    # Dose Interval:
    if (is.null(intervaldose)){
      intervaldose <- "qd"
    }

    # DayOfFirstDrugAdmin
    if (is.null(DayOfFirstDrugAdmin)){
      DayOfFirstDrugAdmin <- 3
    }

    # Import Data:
    data <- import_SCIDpkpdData_GSK(dataFile            = dataFile,
                                    compound            = compound,
                                    PKsheet             = PKsheet,
                                    PKrange             = PKrange,
                                    PKmapping           = PKmapping,
                                    PKlloq              = PKlloq,
                                    PKlloqIdentifier    = PKlloqIdentifier,
                                    PKfactor            = PKfactor,
                                    PDsheet             = PDsheet,
                                    PDrange             = PDrange,
                                    PDmapping           = PDmapping,
                                    PDlloq              = PDlloq,
                                    PDlloqIdentifier    = PDlloqIdentifier,
                                    PDname              = PDname,
                                    study               = study,
                                    DayOfFirstDrugAdmin = DayOfFirstDrugAdmin,
                                    centerNumber        = centerNumber,
                                    centerName          = centerName,
                                    visitNumber         = visitNumber,
                                    nrdoses             = nrdoses,
                                    intervaldose        = intervaldose,
                                    positiveQC          = positiveQC)

  }else{
    stop("'", importVersion, "' is not a valid 'importVersion'; Please choose between 'TAD', 'oldTAD' and 'GSK-TresCantos'")
  }

  # Output:
  return(data)
}

#' import_SCIDpkpdData_GSK
#'
#' @description
#' @param dataFile
#' @param compound Default: NULL
#' @param PKsheet Default: 'BloodLevels'
#' @param PKrange Default: NULL
#' @param PKmapping
#' @param PKlloq Default: 1
#' @param PKlloqIdentifier Default: c("BLQ", "<LLOQ", "< LLOQ")
#' @param PKfactor Default: 0.001
#' @param PDsheet Default: 'Parasitemia'
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01", "< 0,01", "< 0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param study Default: NULL
#' @param DayOfFirstDrugAdmin Default: 3
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @param nrdoses Default: NULL
#' @param intervaldose Default: 'qd'
#' @param positiveQC Default: NULL
#' @return
#' @export
#' @author catalinbarcelo (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
#' @importFrom plyr rename join ddply
#' @importFrom MMVbase aux_createUSUBJID aux_removeEscapeChar
import_SCIDpkpdData_GSK <- function(dataFile,
                                    compound            = NULL,
                                    PKsheet             = "BloodLevels",
                                    PKrange             = NULL,
                                    PKmapping,
                                    PKlloq              = 1,
                                    PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ"),
                                    PKfactor            = 1e-3,  # From ng/mL to ug/mL
                                    PDsheet             = "Parasitemia",
                                    PDrange             = NULL,
                                    PDmapping,
                                    PDlloq              = 0.01,
                                    PDlloqIdentifier    = c("BLQ", "<0.01", "< 0,01", "< 0.01"),
                                    PDname              = "Parasitemia Total",
                                    study               = NULL,
                                    DayOfFirstDrugAdmin = 3,
                                    centerNumber        = -1,
                                    centerName          = "",
                                    visitNumber         = -1,
                                    nrdoses             = NULL,
                                    intervaldose        = "qd",   # only qd and bid allowed
                                    positiveQC          = NULL
) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Final Column Names ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  colNames <- c("IGNORE"   , "STUDY"     , "USUBJID" , "GROUP"   , "SUBJECT",
                "CENTER"   , "CENTERNAME", "COMPOUND", "TRTNAME" , "VISIT"  ,
                "TIME"     , "NT"        , "TIMEUNIT", "TYPENAME", "NAME"   ,
                "VALUE"    , "VALUETXT"  , "UNIT"    , "LLOQ"    , "ROUTE"  ,
                "DOSELEVEL", "DOSEMULT")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Input Variables ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check is the PK column have the appropriate information:
  reqColsPK <- c("SUBJECT", "DOSELEVEL", "NT", "VALUE")
  if (!all(reqColsPK %in% PKmapping))
    stop("Not all required PK columns supplied in mapping: see import_SCIDpkpdData.R")

  # Check is the number of dose is defined:
  if (is.null(nrdoses) & !("DOSEMULT" %in% PKmapping))
    stop("If number of doses not supplied in PK data set, needs to be defined by input argument.")

  # Check is the PD column have the appropriate information:
  reqColsPD <- c("SUBJECT", "DOSELEVEL", "NT", "VALUE", "huErythro")
  if (!all(reqColsPD %in% PDmapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdData.R")

  # Check is the number of dose is defined:
  if (is.null(nrdoses) & !("DOSEMULT" %in% PDmapping))
    stop("If no of doses not supplied in PD data set, needs to be defined by input argument.")

  # Check whether dosing interval either bid or qd
  if (!intervaldose %in% c("bid", "qd")) stop("Only once or twice daily dosing implemented.")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PK data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load PK dataset:
  dataPKraw <- readxl::read_excel(dataFile, PKsheet, range = PKrange)

  # Remove escape characters from column names
  names(dataPKraw) <- aux_removeEscapeChar(names(dataPKraw))

  # Apply mapping:
  dataPK <- plyr::rename(dataPKraw, replace = PKmapping)
  dataPK <- as.data.frame(dataPK)

  # Add number of doses:
  if (!is.null(nrdoses)){
    dataPK$DOSEMULT <- nrdoses
  }

  # Make sure DOSELEVEL and DOSEMULT are numeric:
  dataPK$DOSELEVEL <- as.numeric(dataPK$DOSELEVEL)
  dataPK$DOSEMULT  <- as.numeric(dataPK$DOSEMULT)

  # Adjust DOSELMULT:
  dataPK$DOSEMULT <- ifelse(dataPK$DOSELEVEL==0, 0, dataPK$DOSEMULT)



  # Remove NA rows:
  dataPK <- dataPK[!is.na(dataPK$VALUE),]

  # use STUDY column if exist, otherwise "study":
  if ("STUDY" %in% names(dataPK)) {
    if (is.character(study)) cat("\nStudy taken from STUDY column, given study entry not used.\n")
  } else {
    if (is.character(study)) dataPK$STUDY <- study
    else stop("STUDY needs to be defined either by column or by study argument as character.")
  }

  # use COMPOUND column if exist, otherwise "compound":
  if ("COMPOUND" %in% names(dataPK)) {
    if (is.character(compound)) cat("\nCompound taken from COMPOUND column, given compound entry not used.\n")
  } else {
    if (is.character(compound)) dataPK$COMPOUND <- compound
    else stop("COMPOUND needs to be defined either by column or by compound argument as character.")
  }

  #Check if only one compound was measured for PK:
  FLAGcpdGiven <- FALSE
  if (!is.character(compound)) {compound <- unique(dataPK$COMPOUND); FLAGcpdGiven <- TRUE}
  if (!length(compound) == 1) stop("More than one compounds measured for PK? This is not supported.")

  # Define extra columns for PK:
  dataPK <- within(dataPK, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    TIME      <- NT
    TYPENAME  <- "PK"
    NAME      <- paste0(compound, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq/2/PKfactor, VALUE)) * PKfactor
    VALUETXT  <- NA
    UNIT      <- "ug/mL"
    ROUTE     <- NA
    LLOQ      <- PKlloq
    TRTNAME   <- paste0(ifelse(DOSEMULT == 1, paste0(DOSELEVEL, "mg/kg"), paste0(DOSEMULT, "x", DOSELEVEL, "mg/kg")), compound)
    IGNORE    <- NA
  })

  # Adjust treatment name if dosing twice per day instead of once a day
  if (intervaldose == "bid")
    dataPK$TRTNAME <- paste0(dataPK$TRTNAME, " (bid)")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PD data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load PD dataset:
  dataPDraw <- readxl::read_excel(dataFile, PDsheet, range = PDrange)

  # Remove escape characters from column names
  names(dataPDraw) <- aux_removeEscapeChar(names(dataPDraw))

  # Apply mapping:
  dataPD <- plyr::rename(dataPDraw, replace = PDmapping)
  dataPD <- as.data.frame(dataPD)

  # Add number of doses:
  if (!is.null(nrdoses)){
    dataPD$DOSEMULT <- nrdoses
  }

  # Make sure DOSELEVEL and DOSEMULT are numeric:
  dataPD$DOSELEVEL <- as.numeric(dataPD$DOSELEVEL)
  dataPD$DOSEMULT  <- as.numeric(dataPD$DOSEMULT)

  # Adjust DOSELMULT:
  dataPD$DOSEMULT <- ifelse(dataPD$DOSELEVEL==0, 0, dataPD$DOSEMULT)

  # Remove NA rows:
  dataPD <- dataPD[!is.na(dataPD$VALUE),]

  # use STUDY column if exist, otherwise "study":
  if ("STUDY" %in% names(dataPD)) {
    if (is.character(study)) cat("\nStudy taken from STUDY column, given study entry not used.\n")
  } else {
    if (is.character(study)) dataPD$STUDY <- study
    else stop("STUDY needs to be defined either by column or by study argument as character.")
  }

  # use COMPOUND column if exist, otherwise "compound":
  if ("COMPOUND" %in% names(dataPD)) {
    if (is.character(compound))
      if (FLAGcpdGiven) {
        cat("\nCompound taken from COMPOUND column, given compound entry not used.\n")
      }
  } else {
    if (is.character(compound)) {
      dataPD$COMPOUND <- compound
      if (!FLAGcpdGiven) cat("\nCompound for PD data taken from PK data.\n")
    } else stop("COMPOUND needs to be defined either by column or by compound argument as character.")
  }

  # Define Group TRT:
  if (!("GROUP" %in% names(dataPD))) {
    groupInfo <- unique(dataPD[, c("DOSELEVEL", "DOSEMULT")])
    groupInfo <- rbind(
      groupInfo[groupInfo$DOSELEVEL == 0,][with(groupInfo[groupInfo$DOSELEVEL == 0,], order(DOSEMULT)),],
      groupInfo[groupInfo$DOSELEVEL >  0,][with(groupInfo[groupInfo$DOSELEVEL >  0,], order(DOSEMULT, DOSELEVEL)),]
    )
    groupInfo$GROUP <- 1:dim(groupInfo)[1]
    dataPD <- plyr::join(dataPD, groupInfo)
  } else {
    groupInfo <- unique(dataPD[, c("DOSELEVEL", "DOSEMULT", "GROUP")])
  }

  # Remove all controls except vehicle
  # idx_keep <- (grepl(compound, dataPD$COMPOUND) |dataPD$DOSELEVEL==0)
  # dataPD   <- dataPD[idx_keep,]
  idx_keep <- (grepl(paste0(c(compound, positiveQC), collapse = "|"), dataPD$COMPOUND) |dataPD$DOSELEVEL==0)
  dataPD   <- dataPD[idx_keep,]

  # Define extra columns for PD:
  dataPD <- within(dataPD, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    NT        <- (NT-DayOfFirstDrugAdmin)*24
    TIME      <- NT
    TYPENAME  <- "Efficacy"
    NAME      <- PDname
    IGNORE    <- ifelse(VALUE %in% c("Dead"), "Dead", NA)
    VALUE     <- ifelse(VALUE %in% c("Dead"), NA, VALUE)
    VALUE     <- ifelse(VALUE %in% PDlloqIdentifier, 0, VALUE)
    VALUE     <- as.numeric(VALUE)
    VALUE     <- ifelse(VALUE <= PDlloq, PDlloq/2, VALUE)
    VALUETXT  <- NA
    UNIT      <- "percent"
    ROUTE     <- NA
    LLOQ      <- PDlloq
    TRTNAME   <- ifelse(DOSELEVEL == 0, "Vehicle", paste0(ifelse(DOSEMULT == 1, paste0(DOSELEVEL, "mg/kg"), paste0(DOSEMULT, "x", DOSELEVEL, "mg/kg")), COMPOUND))
    huErythro <- ifelse(huErythro %in% c("Dead"), NA, huErythro)
    huErythro <- as.numeric(huErythro)
  })

  # Adjust treatment name if dosing twice per day instead of once a day
  if (intervaldose == "bid")
    dataPD$TRTNAME <- ifelse(dataPD$TRTNAME=="Vehicle",dataPD$TRTNAME,paste0(dataPD$TRTNAME, " (bid)"))

  # Normalize parasitemia to human erythrocyte percent:
  dataPD$VALUE <- ifelse(dataPD$VALUE == PDlloq/2,PDlloq/2, dataPD$VALUE /dataPD$huErythro * 100)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Erythrocyte Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dataHE   <- within(dataPD, {VALUE    <- huErythro
  TYPENAME <- "Vital Signs"
  NAME     <- "Human Erythrocytes"})


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create USUBJID ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Add group information to PK data
  if (!("GROUP" %in% names(dataPK)))
    dataPK <- plyr::join(dataPK, groupInfo)
  if (!("GROUP" %in% names(dataPD)))
    dataPD <- plyr::join(dataPD, groupInfo)
  if (!("GROUP" %in% names(dataHE)))
    dataHE <- plyr::join(dataHE, groupInfo)

  # Add UBSUBJID
  dataPK$USUBJID <- aux_createUSUBJID(dataPK)
  dataPD$USUBJID <- aux_createUSUBJID(dataPD)
  dataHE$USUBJID <- aux_createUSUBJID(dataHE)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Dosing Records ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Define columns to retain from PK and PD dataset: Only keep one line per subject:
  colDos  <- c("USUBJID", "STUDY", "GROUP", "TRTNAME", "SUBJECT", "COMPOUND", "DOSELEVEL", "DOSEMULT", "TIMEUNIT")
  dataDos <- unique(rbind(dataPK[, colDos],dataPD[, colDos]))

  # Check that there is only one line per subject:
  if (any(duplicated(dataDos$USUBJID)))
    stop("Dosing information not unique for individuals: Check dataset")

  if (intervaldose == "bid") II <- 12 else II <- 24
  # Create Dose events for each subject:
  dataDos <- plyr::ddply(dataDos, ~USUBJID, function(x) {
    if (x$DOSEMULT > 0){
      out <- x[rep(1,x$DOSEMULT),]
      out <- within(out,{TIME      <- seq(0, by = II, length.out = x$DOSEMULT)
      NT        <- TIME
      TYPENAME  <- "Dose"
      NAME      <- ifelse(COMPOUND==compound,
                          paste0(compound, " Dose"),
                          paste0(positiveQC, " Dose"))
      VALUE     <- x$DOSELEVEL
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- "ORAL"
      IGNORE    <- NA
      }
      )
    }else{
      out <- x
      out <- within(out,{TIME <- 0
      NT        <- TIME
      TYPENAME  <- "Dose"
      NAME      <- paste0(compound, " Dose")
      VALUE     <- x$DOSELEVEL
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- "ORAL"
      IGNORE    <- NA
      }
      )
    }
  })

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Remove positive controls ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # # Identify mice with measured PK and with doses of 0:
  # IDwPK <- unique(dataPK$USUBJID)
  # IDw0  <- unique(dataDos$USUBJID[dataDos$DOSELEVEL == 0])
  #
  # # Identify positive controls as the ones having no PK measurement, but a dose > 0:
  # IDposCtrl <- setdiff(dataPD$USUBJID, c(IDwPK, IDw0))
  # dataPD    <- dataPD[!(dataPD$USUBJID %in% IDposCtrl), ]
  # dataHE    <- dataHE[!(dataHE$USUBJID %in% IDposCtrl), ]
  # dataDos   <- dataDos[!(dataDos$USUBJID %in% IDposCtrl), ]
  # if (length(IDposCtrl) > 0)
  #   cat(paste0("\nRemove ", length(IDposCtrl), " positive control mice from data (no PK and dose > 0): ", paste0(IDposCtrl, collapse = ", "), "\n"))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Concatenate All Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Add CENTER & VISIT column:
  dataPK$CENTER      <- centerNumber
  dataPK$CENTERNAME  <- centerName
  dataPK$VISIT       <- visitNumber
  dataPD$CENTER      <- centerNumber
  dataPD$CENTERNAME  <- centerName
  dataPD$VISIT       <- visitNumber
  dataHE$CENTER      <- centerNumber
  dataHE$CENTERNAME  <- centerName
  dataHE$VISIT       <- visitNumber
  dataDos$CENTER     <- centerNumber
  dataDos$CENTERNAME <- centerName
  dataDos$VISIT      <- visitNumber

  # Bind all data:
  data <- rbind(dataPK[,colNames],
                dataPD[,colNames],
                dataHE[,colNames],
                dataDos[,colNames])

  # Order:
  data <- data[order(data$USUBJID, data$TIME),]

  # Check whether mice with dose > 0 have only one compound assigned and assign this also to untreated mice:
  # cpdCheck <- unique(data$COMPOUND[data$DOSELEVEL > 0])
  # if (length(cpdCheck) > 1) stop("Compound information for subjects with dose > 0 not unique.")
  # data$COMPOUND <- cpdCheck
  cpdCheck <- unique(data$COMPOUND[data$DOSELEVEL > 0])
  if (length(setdiff(cpdCheck,positiveQC)) > 1) stop("Compound information for subjects with dose > 0 not unique.")
  if (length(cpdCheck)>1) data$COMPOUND <- setdiff(cpdCheck,positiveQC)

  # Output:
  return(data)
}
#' import_SCIDpkpdData_oldTAD
#'
#' @description
#' @param dataFile
#' @param compound Default: NULL
#' @param PKsheet Default: 'BloodLevels'
#' @param PKrange Default: NULL
#' @param PKmapping
#' @param PKlloq Default: 1
#' @param PKlloqIdentifier Default: c("BLQ", "<LLOQ", "< LLOQ")
#' @param PKfactor Default: 0.001
#' @param PDsheet Default: 'Parasitemia'
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01", "< 0,01", "< 0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param study Default: NULL
#' @param DayOfFirstDrugAdmin Default: 3
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @param nrdoses Default: NULL
#' @param intervaldose Default: 'qd'
#' @param positiveQC Default: NULL
#' @return
#' @export
#' @author catalinbarcelo (MMV)
#' @family Data Preparation
#' @importFrom plyr rename join ddply
import_SCIDpkpdData_oldTAD <- function(dataFile,
                                       compound            = NULL,
                                       PKsheet             = "BloodLevels",
                                       PKrange             = NULL,
                                       PKmapping,
                                       PKlloq              = 1,
                                       PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ"),
                                       PKfactor            = 1e-3,  # From ng/mL to ug/mL
                                       PDsheet             = "Parasitemia",
                                       PDrange             = NULL,
                                       PDmapping,
                                       PDlloq              = 0.01,
                                       PDlloqIdentifier    = c("BLQ", "<0.01", "< 0,01", "< 0.01"),
                                       PDname              = "Parasitemia Total",
                                       study               = NULL,
                                       DayOfFirstDrugAdmin = 3,
                                       centerNumber        = -1,
                                       centerName          = "",
                                       visitNumber         = -1,
                                       nrdoses             = NULL,
                                       intervaldose        = "qd",   # only qd and bid allowed
                                       positiveQC          = NULL
) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Final Column Names ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  colNames <- c("IGNORE"   , "STUDY"     , "USUBJID" , "GROUP"   , "SUBJECT",
                "CENTER"   , "CENTERNAME", "COMPOUND", "TRTNAME" , "VISIT"  ,
                "TIME"     , "NT"        , "TIMEUNIT", "TYPENAME", "NAME"   ,
                "VALUE"    , "VALUETXT"  , "UNIT"    , "LLOQ"    , "ROUTE"  ,
                "DOSELEVEL", "DOSEMULT")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Input Variables ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check is the PK column have the appropriate information:
  reqColsPK <- c("SUBJECT", "DOSELEVEL", "NT", "VALUE")
  if (!all(reqColsPK %in% PKmapping))
    stop("Not all required PK columns supplied in mapping: see import_SCIDpkpdData.R")

  # Check is the number of dose is defined:
  if (is.null(nrdoses) & !("DOSEMULT" %in% PKmapping))
    stop("If number of doses not supplied in PK data set, needs to be defined by input argument.")

  # Check is the PD column have the appropriate information:
  reqColsPD <- c("SUBJECT", "DOSELEVEL", "NT", "VALUE", "huErythro")
  if (!all(reqColsPD %in% PDmapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdData.R")

  # Check is the number of dose is defined:
  if (is.null(nrdoses) & !("DOSEMULT" %in% PDmapping))
    stop("If no of doses not supplied in PD data set, needs to be defined by input argument.")

  # Check whether dosing interval either bid or qd
  if (!intervaldose %in% c("bid", "qd")) stop("Only once or twice daily dosing implemented.")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PK data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load PK dataset:
  dataPKraw <- readxl::read_excel(dataFile, PKsheet, range = PKrange)

  # Remove escape characters from column names
  names(dataPKraw) <- aux_removeEscapeChar(names(dataPKraw))

  # Apply mapping:
  dataPK <- plyr::rename(dataPKraw, replace = PKmapping)
  dataPK <- as.data.frame(dataPK)

  # Add number of doses:
  if (!is.null(nrdoses)){
    dataPK$DOSEMULT <- nrdoses
  }
  dataPK$DOSEMULT <- ifelse(dataPK$DOSELEVEL==0, 0, dataPK$DOSEMULT)

  # Remove NA rows:
  dataPK <- dataPK[!is.na(dataPK$VALUE),]

  # use STUDY column if exist, otherwise "study":
  if ("STUDY" %in% names(dataPK)) {
    if (is.character(study)) cat("\nStudy taken from STUDY column, given study entry not used.\n")
  } else {
    if (is.character(study)) dataPK$STUDY <- study
    else stop("STUDY needs to be defined either by column or by study argument as character.")
  }

  # use COMPOUND column if exist, otherwise "compound":
  if ("COMPOUND" %in% names(dataPK)) {
    if (is.character(compound)) cat("\nCompound taken from COMPOUND column, given compound entry not used.\n")
  } else {
    if (is.character(compound)) dataPK$COMPOUND <- compound
    else stop("COMPOUND needs to be defined either by column or by compound argument as character.")
  }

  #Check if only one compound was measured for PK:
  FLAGcpdGiven <- FALSE
  if (!is.character(compound)) {compound <- unique(dataPK$COMPOUND); FLAGcpdGiven <- TRUE}
  if (!length(compound) == 1) stop("More than one compounds measured for PK? This is not supported.")

  # Define extra columns for PK:
  dataPK <- within(dataPK, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    TIME      <- NT
    TYPENAME  <- "PK"
    NAME      <- paste0(COMPOUND, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq/2/PKfactor, VALUE)) * PKfactor
    VALUETXT  <- NA
    UNIT      <- "ug/mL"
    ROUTE     <- NA
    LLOQ      <- PKlloq
    TRTNAME   <- paste0(ifelse(DOSEMULT == 1, paste0(DOSELEVEL, "mg/kg"), paste0(DOSEMULT, "x", DOSELEVEL, "mg/kg")), COMPOUND)
    IGNORE    <- NA
  })

  # Adjust treatment name if dosing twice per day instead of once a day
  if (intervaldose == "bid")
    dataPK$TRTNAME <- paste0(dataPK$TRTNAME, " (bid)")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PD data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load PD dataset:
  dataPDraw <- readxl::read_excel(dataFile, PDsheet, range = PDrange)

  # Remove escape characters from column names
  names(dataPDraw) <- aux_removeEscapeChar(names(dataPDraw))

  # Apply mapping:
  dataPD <- plyr::rename(dataPDraw, replace = PDmapping)
  dataPD <- as.data.frame(dataPD)

  # Add number of doses:
  if (!is.null(nrdoses)){
    dataPD$DOSEMULT <- nrdoses
  }
  dataPD$DOSEMULT <- ifelse(dataPD$DOSELEVEL==0, 0, dataPD$DOSEMULT)

  # Remove NA rows:
  dataPD <- dataPD[!is.na(dataPD$VALUE),]

  # use STUDY column if exist, otherwise "study":
  if ("STUDY" %in% names(dataPD)) {
    if (is.character(study)) cat("\nStudy taken from STUDY column, given study entry not used.\n")
  } else {
    if (is.character(study)) dataPD$STUDY <- study
    else stop("STUDY needs to be defined either by column or by study argument as character.")
  }

  # use COMPOUND column if exist, otherwise "compound":
  if ("COMPOUND" %in% names(dataPD)) {
    if (is.character(compound))
      if (FLAGcpdGiven) {
        cat("\nCompound taken from COMPOUND column, given compound entry not used.\n")
      }
  } else {
    if (is.character(compound)) {
      dataPD$COMPOUND <- compound
      if (!FLAGcpdGiven) cat("\nCompound for PD data taken from PK data.\n")
    } else stop("COMPOUND needs to be defined either by column or by compound argument as character.")
  }

  # Define Group TRT:
  if (!("GROUP" %in% names(dataPD))) {
    groupInfo <- unique(dataPD[, c("DOSELEVEL", "DOSEMULT")])
    groupInfo <- rbind(
      groupInfo[groupInfo$DOSELEVEL == 0,][with(groupInfo[groupInfo$DOSELEVEL == 0,], order(DOSEMULT)),],
      groupInfo[groupInfo$DOSELEVEL >  0,][with(groupInfo[groupInfo$DOSELEVEL >  0,], order(DOSEMULT, DOSELEVEL)),]
    )
    groupInfo$GROUP <- 1:dim(groupInfo)[1]
    dataPD <- plyr::join(dataPD, groupInfo)
  } else {
    groupInfo <- unique(dataPD[, c("DOSELEVEL", "DOSEMULT", "GROUP")])
  }
  # Remove all controls except vehicle
  idx_keep <- (grepl(paste0(c(compound, positiveQC), collapse = "|"), dataPD$COMPOUND) |dataPD$DOSELEVEL==0)
  dataPD   <- dataPD[idx_keep,]

  # Define extra columns for PD:
  dataPD <- within(dataPD, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    NT        <- (NT-DayOfFirstDrugAdmin)*24
    TIME      <- NT
    TYPENAME  <- "Efficacy"
    NAME      <- PDname
    IGNORE    <- ifelse(VALUE %in% c("Dead"), "Dead", NA)
    VALUE     <- ifelse(VALUE %in% c("Dead"), NA, VALUE)
    VALUE     <- ifelse(VALUE %in% PDlloqIdentifier, 0, VALUE)
    VALUE     <- as.numeric(VALUE)
    VALUE     <- ifelse(VALUE <= PDlloq, PDlloq/2, VALUE)
    VALUETXT  <- NA
    UNIT      <- "percent"
    ROUTE     <- NA
    LLOQ      <- PDlloq
    TRTNAME   <- ifelse(DOSELEVEL == 0, "Vehicle", paste0(ifelse(DOSEMULT == 1, paste0(DOSELEVEL, "mg/kg"), paste0(DOSEMULT, "x", DOSELEVEL, "mg/kg")), COMPOUND))
    huErythro <- ifelse(huErythro %in% c("Dead"), NA, huErythro)
    huErythro <- as.numeric(huErythro)
  })

  # Adjust treatment name if dosing twice per day instead of once a day
  if (intervaldose == "bid")
    dataPD$TRTNAME <- ifelse(dataPD$TRTNAME=="Vehicle",dataPD$TRTNAME,paste0(dataPD$TRTNAME, " (bid)"))

  # Normalize parasitemia to human erythrocyte percent:
  dataPD$VALUE <- ifelse(dataPD$VALUE == PDlloq/2,PDlloq/2, dataPD$VALUE /dataPD$huErythro * 100)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Erythrocyte Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dataHE   <- within(dataPD, {VALUE    <- huErythro
  TYPENAME <- "Vital Signs"
  NAME     <- "Human Erythrocytes"})


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create USUBJID ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Add group information to PK data
  if (!("GROUP" %in% names(dataPK)))
    dataPK <- plyr::join(dataPK, groupInfo)
  if (!("GROUP" %in% names(dataPD)))
    dataPD <- plyr::join(dataPD, groupInfo)
  if (!("GROUP" %in% names(dataHE)))
    dataHE <- plyr::join(dataHE, groupInfo)

  # Add UBSUBJID
  dataPK$USUBJID <- aux_createUSUBJID(dataPK)
  dataPD$USUBJID <- aux_createUSUBJID(dataPD)
  dataHE$USUBJID <- aux_createUSUBJID(dataHE)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Dosing Records ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Define columns to retain from PK and PD dataset: Only keep one line per subject:
  colDos  <- c("USUBJID", "STUDY", "GROUP", "TRTNAME", "SUBJECT", "COMPOUND", "DOSELEVEL", "DOSEMULT", "TIMEUNIT")
  dataDos <- unique(rbind(dataPK[, colDos],dataPD[, colDos]))

  # Check that there is only one line per subject:
  if (any(duplicated(dataDos$USUBJID)))
    stop("Dosing information not unique for individuals: Check dataset")

  if (intervaldose == "bid") II <- 12 else II <- 24
  # Create Dose events for each subject:
  dataDos <- plyr::ddply(dataDos, ~USUBJID, function(x) {
    if (x$DOSEMULT > 0){
      out <- x[rep(1,x$DOSEMULT),]
      out <- within(out,{TIME      <- seq(0, by = II, length.out = x$DOSEMULT)
      NT        <- TIME
      TYPENAME  <- "Dose"
      NAME      <- ifelse(COMPOUND==compound,
                          paste0(compound, " Dose"),
                          paste0(positiveQC, " Dose"))
      VALUE     <- x$DOSELEVEL
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- "ORAL"
      IGNORE    <- NA
      }
      )
    }else{
      out <- x
      out <- within(out,{TIME <- 0
      NT        <- TIME
      TYPENAME  <- "Dose"
      NAME      <- paste0(compound, " Dose")
      VALUE     <- x$DOSELEVEL
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- "ORAL"
      IGNORE    <- NA
      }
      )
    }
  })

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Remove positive controls ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # # Identify mice with measured PK and with doses of 0:
  # IDwPK <- unique(dataPK$USUBJID)
  # IDw0  <- unique(dataDos$USUBJID[dataDos$DOSELEVEL == 0])
  #
  # # Identify positive controls as the ones having no PK measurement, but a dose > 0:
  # IDposCtrl <- setdiff(dataPD$USUBJID, c(IDwPK, IDw0))
  # dataPD    <- dataPD[!(dataPD$USUBJID %in% IDposCtrl), ]
  # dataHE    <- dataHE[!(dataHE$USUBJID %in% IDposCtrl), ]
  # dataDos   <- dataDos[!(dataDos$USUBJID %in% IDposCtrl), ]
  # if (length(IDposCtrl) > 0)
  #   cat(paste0("\nRemove ", length(IDposCtrl), " positive control mice from data (no PK and dose > 0): ", paste0(IDposCtrl, collapse = ", "), "\n"))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Concatenate All Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Add CENTER & VISIT column:
  dataPK$CENTER      <- centerNumber
  dataPK$CENTERNAME  <- centerName
  dataPK$VISIT       <- visitNumber
  dataPD$CENTER      <- centerNumber
  dataPD$CENTERNAME  <- centerName
  dataPD$VISIT       <- visitNumber
  dataHE$CENTER      <- centerNumber
  dataHE$CENTERNAME  <- centerName
  dataHE$VISIT       <- visitNumber
  dataDos$CENTER     <- centerNumber
  dataDos$CENTERNAME <- centerName
  dataDos$VISIT      <- visitNumber

  # Bind all data:
  data <- rbind(dataPK[,colNames],
                dataPD[,colNames],
                dataHE[,colNames],
                dataDos[,colNames])

  # Order:
  data <- data[order(data$USUBJID, data$TIME),]

  # Check whether mice with dose > 0 have only one compound assigned and assign this also to untreated mice:
  cpdCheck <- unique(data$COMPOUND[data$DOSELEVEL > 0])
  if (length(setdiff(cpdCheck,positiveQC)) > 1) stop("Compound information for subjects with dose > 0 not unique.")
  if (length(cpdCheck)>1) data$COMPOUND <- setdiff(cpdCheck,positiveQC)

  # Output:
  return(data)
}
#' Import TAD Data to MMV format
#'
#' @description
#' @param dataFile
#' @param compound
#' @param PKsheet Default: 'DrugConcentrationData'
#' @param PKrange Default: NULL
#' @param PKmapping
#' @param PKlloq Default: 1
#' @param PKlloqIdentifier Default: c("BLQ", "<LLOQ", "< LLOQ")
#' @param PKfactor Default: 0.001
#' @param PDsheet Default: 'ParasitemiaData'
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01", "< 0,01", "< 0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param DayOfFirstDrugAdmin Default: 1
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @param nrdoses Default: NULL
#' @param intervaldose Default: NULL
#' @param route Default: NULL
#' @param DoseTimeSheet Default: 'DosingTimeTable'
#' @param DoseTimeRange Default: NULL
#' @param DoseTimeMapping
#' @param DoseSheet Default: 'DrugTreatmentTable'
#' @param DoseRange Default: NULL
#' @param DoseMapping Default: NULL
#' @param positiveQC
#' @return
#' @export
#' @author Aline Fuchs (MMV), Catalin Barcelo (MMV), Mohammed H. Cherkaoui (MMV), Nathalie Gobeau (MMV)
#' @family Data Preparation
#' @importFrom plyr rename join ddply
import_SCIDpkpdData_TAD <- function(dataFile,
                                    compound,
                                    PKsheet             = "DrugConcentrationData",
                                    PKrange             = NULL,
                                    PKmapping,
                                    PKlloq              = 1,
                                    PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ"),
                                    PKfactor            = 1e-3,  # From ng/mL to ug/mL
                                    PDsheet             = "ParasitemiaData",
                                    PDrange             = NULL,
                                    PDmapping,
                                    PDlloq              = 0.01,
                                    PDlloqIdentifier    = c("BLQ", "<0.01", "< 0,01", "< 0.01"),
                                    PDname              = "Parasitemia Total",
                                    DayOfFirstDrugAdmin = 1,
                                    centerNumber        = -1,
                                    centerName          = "",
                                    visitNumber         = -1,
                                    nrdoses             = NULL,
                                    intervaldose        = NULL,   # only id different from once a day; for name purpose
                                    route               = NULL,
                                    DoseTimeSheet       = "DosingTimeTable",
                                    DoseTimeRange       = NULL,
                                    DoseTimeMapping,
                                    DoseSheet           = "DrugTreatmentTable",
                                    DoseRange           = NULL,
                                    DoseMapping         = NULL,
                                    positiveQC
) {



  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Input Variables ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check if the PK column have the appropriate information:
  reqColsPK <- c("SUBJECT", "STUDY", "COMPOUND", "NT", "VALUE")
  if (!all(reqColsPK %in% PKmapping))
    stop("Not all required PK columns supplied in mapping: see import_SCIDpkpdData_TAD.R")

  # Check if the PD column have the appropriate information:
  reqColsPD <- c("SUBJECT", "STUDY", "NT", "VALUE", "huErythro")
  if (!all(reqColsPD %in% PDmapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdData_TAD.R")

  # Check if the Dose Time column have the appropriate information:
  reqColsDoseTime <- c("SUBJECT", "STUDY", "COMPOUND", "NT")
  if (!all(reqColsDoseTime %in% DoseTimeMapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdData_TAD.R")

  # Check if the Dose column have the appropriate information:
  reqColsDose <- c("SUBJECT", "STUDY", "COMPOUND", "DOSELEVEL")
  if (!all(reqColsDose %in% DoseMapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdData_TAD.R")

  # Check is the number of dose is defined:
  if (is.null(nrdoses) & !("DOSEMULT" %in% DoseMapping))
    stop("If number of doses not supplied in Dose data set, needs to be defined by input argument.")

  # Check is route of administration is defined:
  if (is.null(route) & !("ROUTE" %in% DoseMapping))
    stop("If route of administration is not supplied in Dose data set, needs to be defined by input argument.")

  # Check whether dosing interval either bid or qd
  if (is.null(DoseTimeSheet))
    warning("Dosing time sheet not provided, once daily dosing is assumed.")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Dosing data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # -- Load DoseTime dataset:
  dataDoseTimeraw <- readxl::read_excel(dataFile, DoseTimeSheet, range = DoseTimeRange)

  # Remove escape characters from column names
  names(dataDoseTimeraw) <- aux_removeEscapeChar(names(dataDoseTimeraw))

  # Apply mapping:
  dataDoseTime <- plyr::rename(dataDoseTimeraw, replace = DoseTimeMapping)
  dataDoseTime <- as.data.frame(dataDoseTime)
  
  # Remove dosing data for compounds other than QC, NONE and compound
  dataDoseTime <- dataDoseTime[dataDoseTime$COMPOUND %in% c(positiveQC, compound, "NONE"), ]

  # Load Dose dataset:
  dataDoseraw <- readxl::read_excel(dataFile, DoseSheet, range = DoseRange)

  # Remove escape characters from column names
  names(dataDoseraw) <- aux_removeEscapeChar(names(dataDoseraw))

  # Apply mapping:
  dataDose <- plyr::rename(dataDoseraw, replace = DoseMapping)
  dataDose <- as.data.frame(dataDose)

  # Remove data for compounds other than QC, NONE and compound
  dataDose <- dataDose[dataDose$COMPOUND %in% c(positiveQC, compound, "NONE"), ]

  # -- Join dosing info:
  dataDos <- plyr::join(dataDoseTime, dataDose)

  # handle multiple dose
  if (!is.null(nrdoses)){
    dataDos$DOSEMULT <- nrdoses
  }
  dataDos$DOSEMULT <- ifelse(dataDos$DOSELEVEL==0, 0, dataDos$DOSEMULT)

  # Remove potential NA entries in doselevel and number of doses: Set to 0
  dataDos$DOSELEVEL <- as.numeric(dataDos$DOSELEVEL)
  dataDos$DOSEMULT <- as.numeric(dataDos$DOSEMULT)
  dataDos$NT <- as.numeric(dataDos$NT)
  dataDos$ROUTE <- as.character(dataDos$ROUTE)

  dataDos <- within(dataDos, {
    DOSELEVEL <- ifelse(is.na(DOSELEVEL), 0, DOSELEVEL)
    DOSEMULT  <- ifelse(DOSELEVEL==0    , 0, DOSEMULT)
    NT <- ifelse(is.na(NT), 0, NT)
    ROUTE <- ifelse(is.na(ROUTE) | ROUTE == "NA", unique(ROUTE[!is.na(ROUTE) & ROUTE != "NA"])[1], ROUTE)
})

  # Define Group TRT: (include positive control; in the new format group should never be a column in any sheet actually)
  if (!("GROUP" %in% names(dataDos))) {
    groupInfo <- unique(dataDos[, c("COMPOUND","DOSELEVEL", "DOSEMULT")])
    groupInfo <- rbind(
      groupInfo[groupInfo$DOSELEVEL == 0,][with(groupInfo[groupInfo$DOSELEVEL == 0,], order(DOSEMULT)),],
      groupInfo[groupInfo$DOSELEVEL >  0,][with(groupInfo[groupInfo$DOSELEVEL >  0,], order(DOSEMULT, DOSELEVEL)),]
    )
    groupInfo$GROUP <- 1:dim(groupInfo)[1]
    dataDos <- plyr::join(dataDos, groupInfo)
  } else {
    groupInfo <- unique(dataDos[, c("COMPOUND","DOSELEVEL", "DOSEMULT", "GROUP")])
  }

  # Add UBSUBJID
  dataDos$USUBJID <- aux_createUSUBJID(dataDos)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Dosing Records ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check that there is only one line per subject:
  if (any(duplicated(dataDos[c("USUBJID","NT")])))
    stop("Dosing information not unique for individuals: Check dataset")

  # Create Dose events for each subject:
  dataDos <- plyr::ddply(dataDos, ~USUBJID, function(x) {
    out <- within(x,{
      TIMEUNIT  <- rep("hours", length(STUDY))
      TIME      <- NT
      TYPENAME  <- "Dose"
      NAME      <- ifelse(COMPOUND!="NONE",paste0(COMPOUND, " Dose"),paste0(compound, " Dose"))
      VALUE     <- x$DOSELEVEL
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- ifelse("ROUTE" %in% names(dataDos),  toupper(ROUTE), toupper(route))
      IGNORE    <- NA
    }
    )
  })

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PK data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load PK dataset:
  dataPKraw <- readxl::read_excel(dataFile, PKsheet, range = PKrange)

  # Remove escape characters from column names
  names(dataPKraw) <- aux_removeEscapeChar(names(dataPKraw))

  # Apply mapping:
  dataPK <- plyr::rename(dataPKraw, replace = PKmapping)
  dataPK <- as.data.frame(dataPK)

  # Remove data for compounds other than QC, NONE and compound
  dataPK <- dataPK[dataPK$COMPOUND %in% c(positiveQC, compound, "NONE"), ]

  # Remove NA rows:
  dataPK <- dataPK[!is.na(dataPK$VALUE),]

  # Get rid of concentration at TIME = 0; no measurement
  dataPK <- dataPK[dataPK$NT!=0,]

  # Define extra columns for PK:
  dataPK <- within(dataPK, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    TIME      <- NT
    TYPENAME  <- "PK"
    NAME      <- paste0(COMPOUND, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq/2/PKfactor, VALUE)) * PKfactor
    VALUETXT  <- NA
    UNIT      <- "ug/mL"
    ROUTE     <- NA
    LLOQ      <- PKlloq
    IGNORE    <- NA
  })


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PD data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load PD dataset:
  dataPDraw <- readxl::read_excel(dataFile, PDsheet, range = PDrange)

  # Remove escape characters from column names
  names(dataPDraw) <- aux_removeEscapeChar(names(dataPDraw))

  # Apply mapping:
  dataPD <- plyr::rename(dataPDraw, replace = PDmapping)
  dataPD <- as.data.frame(dataPD)

  # Remove NA rows:
  dataPD <- dataPD[!is.na(dataPD$VALUE),]

  # Define extra columns for PD:
  dataPD <- within(dataPD, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    NT        <- (NT-DayOfFirstDrugAdmin)*24
    TIME      <- NT
    TYPENAME  <- "Efficacy"
    NAME      <- PDname
    IGNORE    <- ifelse(VALUE %in% c("Dead"), "Dead", NA)
    VALUE     <- ifelse(VALUE %in% c("Dead"), NA, VALUE)
    VALUE     <- ifelse(VALUE %in% PDlloqIdentifier, 0, VALUE)
    VALUE     <- as.numeric(VALUE)
    VALUE     <- ifelse(VALUE <= PDlloq, PDlloq/2, VALUE)
    VALUETXT  <- NA
    UNIT      <- "percent"
    ROUTE     <- NA
    LLOQ      <- PDlloq
    huErythro <- ifelse(huErythro %in% c("Dead"), NA, huErythro)
    huErythro <- as.numeric(huErythro)
  })


  # Normalize parasitemia to human erythrocyte percent:
  dataPD$VALUE <- ifelse(dataPD$VALUE == PDlloq/2,PDlloq/2, dataPD$VALUE /dataPD$huErythro * 100)

  # Add compound to PD
  dataPD <- merge(dataPD,unique(dataDos[,c("SUBJECT","STUDY","COMPOUND")]))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Erythrocyte Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dataHE <- within(dataPD, {VALUE    <- huErythro
  TYPENAME <- "Vital Signs"
  NAME     <- "Human Erythrocytes"})


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Concatenate All Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  colKeep <- intersect(names(dataPK),
                       intersect(names(dataPD),
                                 intersect(names(dataHE),names(dataDos))))
  dataPK <- dataPK[,colKeep]
  dataPD <- dataPD[,colKeep]
  dataHE <- dataHE[,colKeep]
  dataDO <- dataDos[,colKeep]

  data <- rbind(dataPK,
                dataPD,
                dataHE,
                dataDO)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # finalization ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Add GROUP to be further used for USUBJID
  data <- plyr::join(data, unique(dataDos[,c("SUBJECT", "STUDY","COMPOUND", "GROUP", "USUBJID", "DOSELEVEL",  "DOSEMULT" )]))

  # Create TRTNAME
  data$TRTNAME <- with(data, ifelse(DOSELEVEL == 0, "Vehicle", paste0(ifelse(DOSEMULT==1, paste0(DOSELEVEL, "mg/kg"), paste0(DOSEMULT, "x", DOSELEVEL, "mg/kg")), COMPOUND)))

  # Adjust treatment name if dosing different of once a day
  if (!is.null(intervaldose) && !is.na(intervaldose)){
    data$TRTNAME <- ifelse(data$TRTNAME=="Vehicle",data$TRTNAME,paste0(data$TRTNAME, " ",intervaldose))
  }

  # Add CENTER & VISIT column:
  data$CENTER      <- centerNumber
  data$CENTERNAME  <- centerName
  data$VISIT       <- visitNumber

  # Create UBSUBJID
  data$USUBJID <- aux_createUSUBJID(data)

  # Order:
  data <- data[order(data$USUBJID, data$TIME),]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Final Column Names ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  colNames <- c("IGNORE"   , "STUDY"     , "USUBJID" , "GROUP"   , "SUBJECT",
                "CENTER"   , "CENTERNAME", "COMPOUND", "TRTNAME" , "VISIT"  ,
                "TIME"     , "NT"        , "TIMEUNIT", "TYPENAME", "NAME"   ,
                "VALUE"    , "VALUETXT"  , "UNIT"    , "LLOQ"    , "ROUTE"  ,
                "DOSELEVEL", "DOSEMULT")

  data <- data[,colNames]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get data only for compound of interest ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Get STUDY of the Compound of interest ad subset on appropriate study

  StudyID <- unique(data$STUDY[data$COMPOUND %in% compound])
  data    <- data[data$STUDY %in% StudyID,]

  # Select only for compound of interest + vehicle + positive control
  data <- rbind(data[data$COMPOUND==compound,],
                data[data$COMPOUND==positiveQC,],
                data[data$TRTNAME =="Vehicle",])

  # Check whether mice with dose > 0 have only one compound assigned and assign this also to untreated mice:
  cpdCheck <- unique(data$COMPOUND[data$DOSELEVEL > 0 & data$COMPOUND != positiveQC])
  if (length(cpdCheck)>1) stop("Compound information for subjects with dose > 0 & not positive control is not unique.")

  data$COMPOUND <- cpdCheck

  # Output:
  return(data)
}

#' import_SCIDpkpdDataCombo
#'
#' @description
#' @param dataFile
#' @param Compound1 Default: NULL
#' @param Compound2 Default: NULL
#' @param PKsheets Default: NULL
#' @param PKranges Default: NULL
#' @param PKmapping
#' @param PKlloq Default: c(0.001, 0.001)
#' @param PKlloqIdentifier Default: c("BLQ", "<LLOQ", "< LLOQ", "BLOQ")
#' @param PKfactor Default: c(0.001, 0.001)
#' @param PDsheet Default: NULL
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param study Default: NULL
#' @param studyID Default: NULL
#' @param DayOfFirstDrugAdmin Default: NULL
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @param intervaldose Default: NULL
#' @param DoseTimeSheet Default: 'DosingTimeTable'
#' @param DoseTimeRange Default: NULL
#' @param DoseTimeMapping
#' @param DoseSheet Default: 'DrugTreatmentTable'
#' @param DoseRange Default: NULL
#' @param DoseMapping Default: NULL
#' @param positiveQC Default: NULL
#' @param route Default: NULL
#' @param importVersion Default: 'oldTAD'
#' @return
#' @export
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
import_SCIDpkpdDataCombo <- function(dataFile,
                                     Compound1           = NULL,
                                     Compound2           = NULL,
                                     PKsheets            = NULL,
                                     PKranges            = NULL,
                                     PKmapping,
                                     PKlloq              = c(0.001,0.001),   # Needs to be in ug/mL.
                                     PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ", "BLOQ"),
                                     PKfactor            = c(1e-3,1e-3),
                                     PDsheet             = NULL,
                                     PDrange             = NULL,
                                     PDmapping,
                                     PDlloq              = 0.01,
                                     PDlloqIdentifier    = c("BLQ", "<0.01"),
                                     PDname              = "Parasitemia Total",
                                     study               = NULL,
                                     studyID             = NULL,
                                     DayOfFirstDrugAdmin = NULL,
                                     centerNumber        = -1,
                                     centerName          = "",
                                     visitNumber         = -1,
                                     intervaldose        = NULL,
                                     DoseTimeSheet       = "DosingTimeTable",
                                     DoseTimeRange       = NULL,
                                     DoseTimeMapping,
                                     DoseSheet           = "DrugTreatmentTable",
                                     DoseRange           = NULL,
                                     DoseMapping         = NULL,
                                     positiveQC          = NULL,
                                     route               = NULL,
                                     importVersion       = "oldTAD"
) {

  # Use the appropriate import function according to the center name:
  if (importVersion=="TAD"){

    #stop("Function 'import_SCIDpkpdDataCombo_TAD' not ready yet.")

    # PK sheet:
    if (is.null(PKsheets)){
      PKsheets = "DrugConcentrationData"
    }

    # PD sheet:
    if (is.null(PDsheet)){
      PDsheet = "ParasitemiaData"
    }

    # DayOfFirstDrugAdmin
    if (is.null(DayOfFirstDrugAdmin)){
      DayOfFirstDrugAdmin = 1
    }

    # Import Data:
    data <- import_SCIDpkpdDataCombo_TAD(dataFile            = dataFile,
                                         Compound1           = Compound1,
                                         Compound2           = Compound2,
                                         PKsheets             = PKsheets,
                                         PKranges             = PKranges,
                                         PKmapping           = PKmapping,
                                         PKlloq              = PKlloq,
                                         PKlloqIdentifier    = PKlloqIdentifier,
                                         PKfactor            = PKfactor,
                                         PDsheet             = PDsheet,
                                         PDrange             = PDrange,
                                         PDmapping           = PDmapping,
                                         PDlloq              = PDlloq,
                                         PDlloqIdentifier    = PDlloqIdentifier,
                                         PDname              = PDname,
                                         DayOfFirstDrugAdmin = DayOfFirstDrugAdmin,
                                         centerNumber        = centerNumber,
                                         centerName          = centerName,
                                         visitNumber         = visitNumber,
                                         intervaldose        = intervaldose,
                                         route               = route,
                                         DoseTimeSheet       = DoseTimeSheet,
                                         DoseTimeRange       = DoseTimeRange,
                                         DoseTimeMapping     = DoseTimeMapping,
                                         DoseSheet           = DoseSheet,
                                         DoseRange           = DoseRange,
                                         DoseMapping         = DoseMapping,
                                         positiveQC          = positiveQC)


  }else if(importVersion=="GSK-TresCantos"){
    # PK sheet:
    if (is.null(PKsheets)){
      PKsheets = c("Blood levels 1", "Blood levels 2")
    }

    # PD sheet:
    if (is.null(PDsheet)){
      PDsheet = "Parasitemia"
    }

    # DayOfFirstDrugAdmin
    if (is.null(DayOfFirstDrugAdmin)){
      DayOfFirstDrugAdmin = 3
    }

    # Import Data:
    data <- import_SCIDpkpdDataCombo_GSK(dataFile            = dataFile,
                                         Compound1           = Compound1,
                                         Compound2           = Compound2,
                                         PKsheets            = PKsheets,
                                         PKranges            = PKranges,
                                         PKmapping           = PKmapping,
                                         PKlloq              = PKlloq,
                                         PKlloqIdentifier    = PKlloqIdentifier,
                                         PKfactor            = PKfactor,
                                         PDsheet             = PDsheet,
                                         PDrange             = PDrange,
                                         PDmapping           = PDmapping,
                                         PDlloq              = PDlloq,
                                         PDlloqIdentifier    = PDlloqIdentifier,
                                         PDname              = PDname,
                                         study               = study,
                                         studyID             = studyID,
                                         DayOfFirstDrugAdmin = DayOfFirstDrugAdmin,
                                         centerNumber        = centerNumber,
                                         centerName          = centerName,
                                         visitNumber         = visitNumber)


  }else if(importVersion=="oldTAD"){
    # PK sheet:
    if (is.null(PKsheets)){
      PKsheets <- c("Blood levels 1", "Blood levels 2")
    }

    # PD sheet:
    if (is.null(PDsheet)){
      PDsheet <- "Parasitemia"
    }

    # DayOfFirstDrugAdmin
    if (is.null(DayOfFirstDrugAdmin)){
      DayOfFirstDrugAdmin <- 1
    }

    # Import Data:
    data <- import_SCIDpkpdDataCombo_oldTAD(dataFile            = dataFile,
                                            PKsheets            = PKsheets,
                                            PKranges            = PKranges,
                                            PKmapping           = PKmapping,
                                            PKlloq              = PKlloq,
                                            PKlloqIdentifier    = PKlloqIdentifier,
                                            PKfactor            = PKfactor,
                                            PDsheet             = PDsheet,
                                            PDrange             = PDrange,
                                            PDmapping           = PDmapping,
                                            PDlloq              = PDlloq,
                                            PDlloqIdentifier    = PDlloqIdentifier,
                                            PDname              = PDname,
                                            study               = study,
                                            studyID             = studyID,
                                            centerNumber        = centerNumber,
                                            centerName          = centerName,
                                            visitNumber         = visitNumber,
                                            DayOfFirstDrugAdmin = DayOfFirstDrugAdmin)

  }else{
    stop("'", importVersion, "' is not a valid 'importVersion'; Please choose between 'TAD' and 'GSK-TresCantos'")
  }

  # Output:
  return(data)
}


#' import_SCIDpkpdDataCombo_GSK
#'
#' @description
#' @param dataFile
#' @param Compound1 Default: NULL
#' @param Compound2 Default: NULL
#' @param PKsheets Default: c("Blood levels 1", "Blood levels 2")
#' @param PKranges Default: NULL
#' @param PKmapping
#' @param PKlloq Default: c(1, 1)
#' @param PKlloqIdentifier Default: '< LLOQ'
#' @param PKfactor Default: c(0.001, 0.001)
#' @param PKfact Default: PKfactor
#' @param PDsheet Default: 'Parasitemia'
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param study Default: NULL
#' @param studyID Default: NULL
#' @param DayOfFirstDrugAdmin Default: 3
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @return
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
#' @importFrom plyr rename ddply
#' @importFrom dplyr left_join
import_SCIDpkpdDataCombo_GSK <- function(dataFile,
                                         Compound1 = NULL,
                                         Compound2 = NULL,
                                         PKsheets = c("Blood levels 1", "Blood levels 2"),
                                         PKranges = NULL,
                                         PKmapping,
                                         PKlloq  = c(1,1),
                                         PKlloqIdentifier  = "< LLOQ",
                                         PKfactor = c(1e-3,1e-3),
                                         PKfact   = PKfactor,
                                         PDsheet  = "Parasitemia",
                                         PDrange  = NULL,
                                         PDmapping,
                                         PDlloq  = 0.01,
                                         PDlloqIdentifier = c("BLQ", "<0.01"),
                                         PDname           = "Parasitemia Total",
                                         study            = NULL,
                                         studyID          = NULL,
                                         DayOfFirstDrugAdmin = 3,
                                         centerNumber     = -1,
                                         centerName       = "",
                                         visitNumber      = -1
) {

  colNames <- c("IGNORE"    , "STUDY"    , "USUBJID"   , "GROUP"    , "SUBJECT"   ,
                "COMPOUND"  , "TRTNAME"  , "TIME"      , "NT"       , "TIMEUNIT"  ,
                "TYPENAME"  , "NAME"     , "VALUE"     , "VALUETXT", "UNIT"      ,
                "LLOQ"      , "ROUTE"    , "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")


  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  #        PK data

  reqColsPK <- c("SUBJECT", "GROUP", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2", "NT", "VALUE", "COMPOUND")
  if (!all(reqColsPK %in% PKmapping))
    stop("Not all required PK columns supplied.")

  # Load both PK data sets
  dataPKraw1 <- readxl::read_excel(dataFile, PKsheets[1], range = PKranges[1])
  dataPKraw2 <- readxl::read_excel(dataFile, PKsheets[2], range = PKranges[2])

  # Remove escape characters from column names
  names(dataPKraw1) <- aux_removeEscapeChar(names(dataPKraw1))
  names(dataPKraw2) <- aux_removeEscapeChar(names(dataPKraw2))

  # Apply mapping
  dataPK1 <- plyr::rename(dataPKraw1, replace = PKmapping, warn_missing = FALSE)
  dataPK2 <- plyr::rename(dataPKraw2, replace = PKmapping, warn_missing = FALSE)

  # Remove not matched columns
  dataPK1 <- dataPK1[, intersect(names(dataPK1),PKmapping)]
  dataPK2 <- dataPK2[, intersect(names(dataPK2),PKmapping)]

  # Indicate compound number  for DOSELEVEL and DOSEMULT
  #dataPK1 <- plyr::rename(dataPK1, replace = c("DOSELEVEL" = "DOSELEVEL1", "DOSEMULT" = "DOSEMULT1"))
  #dataPK2 <- plyr::rename(dataPK2, replace = c("DOSELEVEL" = "DOSELEVEL2", "DOSEMULT" = "DOSEMULT2"))
  dataPK1 <- dataPK1[,setdiff(names(dataPK1), c("DOSELEVEL2", "DOSEMULT2"))]
  dataPK2 <- dataPK2[,setdiff(names(dataPK2), c("DOSELEVEL1", "DOSEMULT1"))]

  # Remove NA rows
  dataPK1 <- dataPK1[!is.na(dataPK1$VALUE), ]
  dataPK2 <- dataPK2[!is.na(dataPK2$VALUE), ]

  # Get compound names
  compound1 <- grep("+",unique(dataPK1$COMPOUND), invert = TRUE, value = TRUE, fixed = TRUE)
  compound2 <- grep("+",unique(dataPK2$COMPOUND), invert = TRUE, value = TRUE, fixed = TRUE)

  # Define compound specific columns
  dataPK1 <- within(dataPK1, {
    NAME      <- paste0(compound1, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq[1]/2/PKfact[1], VALUE)) * PKfact[1]
    LLOQ      <- PKlloq[1]
  })
  dataPK2 <- within(dataPK2, {
    NAME      <- paste0(compound2, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq[2]/2/PKfact[2], VALUE)) * PKfact[2]
    LLOQ      <- PKlloq[2]
  })

  # Merge treatment information
  dataPK1 <- dplyr::left_join(dataPK1, dataPK2[,c("SUBJECT", "GROUP", "NT", "COMPOUND", "DOSELEVEL2", "DOSEMULT2")])
  dataPK1$DOSELEVEL2 <- with(dataPK1, ifelse(is.na(DOSELEVEL2), 0, DOSELEVEL2))
  dataPK1$DOSEMULT2 <- with(dataPK1, ifelse(is.na(DOSEMULT2), 0, DOSEMULT2))
  dataPK2 <- dplyr::left_join(dataPK2, dataPK1[,c("SUBJECT", "GROUP", "NT", "COMPOUND", "DOSELEVEL1", "DOSEMULT1")])
  dataPK2$DOSELEVEL1 <- with(dataPK2, ifelse(is.na(DOSELEVEL1), 0, DOSELEVEL1))
  dataPK2$DOSEMULT1 <- with(dataPK2, ifelse(is.na(DOSEMULT1), 0, DOSEMULT1))

  # Concatenate PK datasets
  dataPK <- rbind(dataPK1, dataPK2)

  # (Re-)define compound (order as defined by 1 and 2 in this script) and treatment
  dataPK$COMPOUND <- paste0(compound1,"+",compound2)
  dataPK <- within(dataPK, {
    str1 <- ifelse(DOSEMULT1 == 1, paste0(DOSELEVEL1, "mg/kg",compound1), paste0(DOSEMULT1,"x",DOSELEVEL1, "mg/kg",compound1))
    str2 <- ifelse(DOSEMULT2 == 1, paste0(DOSELEVEL2, "mg/kg",compound2), paste0(DOSEMULT2,"x",DOSELEVEL2, "mg/kg",compound2))
    TRTNAME <- ifelse(DOSELEVEL1 > 0 & DOSELEVEL2 == 0, str1, NA)
    TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 > 0, str2, TRTNAME)
    TRTNAME <- ifelse(DOSELEVEL1 > 0 & DOSELEVEL2 > 0, paste0(str1, "+", str2), TRTNAME)
    str1 <- str2 <- NULL
  })

  # use "studyID" if possible, otherwise "study"
  dataPK$STUDYID <- NaN
  if (!is.null(studyID)){
    dataPK$STUDYID <- dataPK[[studyID]]
  } else{
    dataPK$STUDYID <- study
  }

  # Define columns common for both PK
  dataPK <- within(dataPK, {
    STUDY     <- STUDYID
    TIMEUNIT  <- rep("hours", length(STUDY))
    TIME      <- NT
    TYPENAME  <- "PK"
    VALUETXT  <- NA
    UNIT      <- "ug/mL"
    ROUTE     <- NA
    IGNORE    <- NA
  })
  dataPK$USUBJID   <- aux_createUSUBJID(dataPK)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  #        PD data
  reqColsPD <- c("SUBJECT", "GROUP", "NT", "VALUE", "huErythro")
  if (!all(reqColsPD %in% PDmapping))
    stop("Not all required PD columns supplied in mapping.")

  # Load data
  dataPDraw <- readxl::read_excel(dataFile, PDsheet, range = PDrange)

  # Remove escape characters from column names
  names(dataPDraw) <- aux_removeEscapeChar(names(dataPDraw))

  # Apply matich
  dataPD <- plyr::rename(dataPDraw, replace = PDmapping)

  # Remove not matched columns
  dataPD <- dataPD[, intersect(names(dataPD),PDmapping)]

  # use "studyID" if possible, otherwise "study"
  dataPD$STUDYID <- NaN
  if (!is.null(studyID)){
    dataPD$STUDYID <- dataPD[[studyID]]
  } else{
    dataPD$STUDYID <- study
  }

  # Remove records with no entry for value (should remove also the gap between experiment parts)
  dataPD <- dataPD[!is.na(dataPD$VALUE), ]

  dataPD <- within(dataPD, {
    STUDY     <- STUDYID
    TIMEUNIT  <- rep("hours", length(STUDY))
    NT        <- (NT-DayOfFirstDrugAdmin)*24
    TIME      <- NT
    TYPENAME  <- "Efficacy"
    NAME      <- PDname
    VALUE     <- ifelse(VALUE %in% c("Dead"), NA, VALUE)
    IGNORE    <- ifelse(VALUE %in% c("Dead"), "Dead", NA)
    VALUE     <- ifelse(VALUE %in% PDlloqIdentifier, 0, VALUE)
    VALUE     <- as.numeric(VALUE)
    VALUE     <- ifelse(VALUE <= PDlloq, PDlloq/2, VALUE)
    VALUETXT  <- NA
    UNIT      <- "percent"
    ROUTE     <- NA
    LLOQ      <- PDlloq
    huErythro <- ifelse(huErythro %in% c("Dead"), NA, huErythro)
    huErythro <- as.numeric(huErythro)
    DOSELEVEL1 <- as.numeric(DOSELEVEL1)
    DOSEMULT1  <- as.numeric(DOSEMULT1)
    DOSELEVEL2 <- as.numeric(DOSELEVEL2)
    DOSEMULT2  <- as.numeric(DOSEMULT2)
    str1 <- ifelse(DOSEMULT1 == 1, paste0(DOSELEVEL1, "mg/kg",compound1), paste0(DOSEMULT1,"x",DOSELEVEL1, "mg/kg",compound1))
    str2 <- ifelse(DOSEMULT2 == 1, paste0(DOSELEVEL2, "mg/kg",compound2), paste0(DOSEMULT2,"x",DOSELEVEL2, "mg/kg",compound2))
    TRTNAME <- ifelse(DOSELEVEL1 > 0 & DOSELEVEL2 == 0, str1, NA)
    TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 > 0, str2, TRTNAME)
    TRTNAME <- ifelse(DOSELEVEL1 > 0 & DOSELEVEL2 > 0, paste0(str1, "+", str2), TRTNAME)
    str1 <- str2 <- NULL
  })
  dataPD$USUBJID   <- aux_createUSUBJID(dataPD)

  # Normalize parasitemia to human erythrocyte percent
  dataPD$VALUE <- ifelse(dataPD$VALUE == PDlloq/2,PDlloq/2, dataPD$VALUE /dataPD$huErythro * 100)

  # Get treatment information from PK data
  dataPD <- dplyr::left_join(dataPD, unique(dataPK[,c("USUBJID", "TRTNAME", "COMPOUND", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")]))
  dataPD <- within(dataPD, {
    COMPOUND = paste0(compound1, "+", compound2)
    TRTNAME  = ifelse(is.na(TRTNAME), "Vehicle", TRTNAME)
    DOSELEVEL1 <- ifelse(TRTNAME == "Vehicle", 0, DOSELEVEL1)
    DOSELEVEL2 <- ifelse(TRTNAME == "Vehicle", 0, DOSELEVEL2)
    DOSEMULT1 <- ifelse(TRTNAME == "Vehicle", 0, DOSEMULT1)
    DOSEMULT2 <- ifelse(TRTNAME == "Vehicle", 0, DOSEMULT2)
  })

  # Human erythrocyte data
  dataHE <- within(dataPD, {
    VALUE <- huErythro
    TYPENAME <- "Vital signs"
    NAME  <- "Human erythrocytes"
  })


  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  #     Dosing records
  colDos  <- c("USUBJID", "STUDY", "GROUP", "TRTNAME", "SUBJECT", "COMPOUND", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2", "TIMEUNIT")
  dataDos <- unique(rbind(dataPK[, colDos],dataPD[, colDos]))

  if (any(duplicated(dataDos$USUBJID)))
    stop("Dosing information not unique for individuals.")

  dataDos <- plyr::ddply(dataDos, ~USUBJID, function(x) {
    if (x$DOSEMULT1 > 0) {
      out1 <- x[rep(1,x$DOSEMULT1),]
      out1 <- within(out1, {
        TIME     <- seq(0, by = 24, length.out = x$DOSEMULT1)
        NT       <- TIME
        TYPENAME <- "Dose"
        NAME     <- paste0(compound1, " Dose")
        VALUE    <- x$DOSELEVEL1
        VALUETXT <- NA
        UNIT     <- "mg/kg"
        LLOQ     <- NA
        ROUTE    <- "ORAL"
        IGNORE   <- NA
      })
    } else out1 <- NULL

    if (x$DOSEMULT2 > 0) {
      out2 <- x[rep(1,x$DOSEMULT2),]

      out2 <- within(out2, {
        TIME <- seq(0, by = 24, length.out = x$DOSEMULT2)
        NT   <- TIME
        TYPENAME <- "Dose"
        NAME  <- paste0(compound2, " Dose")
        VALUE <- x$DOSELEVEL2
        VALUETXT  <- NA
        UNIT <- "mg/kg"
        LLOQ <- NA
        ROUTE <- "ORAL"
        IGNORE <- NA
      })
    } else out2 <- NULL
    out <- rbind(out1,out2)
  })

  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  #     Concatenate

  data <- rbind(
    dataPK[,colNames],
    dataPD[,colNames],
    dataHE[,colNames],
    dataDos[,colNames]
  )
  # Add CENTER column
  data$CENTER      <- centerNumber
  data$CENTERNAME  <- centerName
  data$VISIT       <- visitNumber

  # order
  data <- data[order(data$USUBJID, data$TIME),]

  # Swap Names:
  if(!is.null(Compound1)){
    data <- within(data,
                   {COMPOUND <- gsub(Compound1$Name,
                                     Compound1$MMVname,
                                     COMPOUND)
                   TRTNAME  <- gsub(Compound1$Name,
                                    Compound1$MMVname,
                                    TRTNAME)
                   NAME     <- gsub(Compound1$Name,
                                    Compound1$MMVname,
                                    NAME)
                   })
  }
  if(!is.null(Compound2)){
    data <- within(data,
                   {COMPOUND <- gsub(Compound2$Name,
                                     Compound2$MMVname,
                                     COMPOUND)
                   TRTNAME  <- gsub(Compound2$Name,
                                    Compound2$MMVname,
                                    TRTNAME)
                   NAME     <- gsub(Compound2$Name,
                                    Compound2$MMVname,
                                    NAME)
                   })
  }

  # output
  data

}

#' import_SCIDpkpdDataCombo_oldTAD
#'
#' @description
#' @param dataFile
#' @param PKsheets Default: c("Blood levels 1", "Blood levels 2")
#' @param PKranges Default: NULL
#' @param PKmapping
#' @param PKlloq Default: c(0.001, 0.001)
#' @param PKlloqIdentifier Default: c("BLQ", "<LLOQ", "< LLOQ", "BLOQ")
#' @param PKfactor Default: c(0.001, 0.001)
#' @param PDsheet Default: 'Parasitemia'
#' @param PDrange Default: NULL
#' @param PDmapping
#' @param PDlloq Default: 0.01
#' @param PDlloqIdentifier Default: c("BLQ", "<0.01")
#' @param PDname Default: 'Parasitemia Total'
#' @param study Default: NULL
#' @param studyID Default: NULL
#' @param DayOfFirstDrugAdmin Default: 1
#' @param centerNumber Default: -1
#' @param centerName Default: ''
#' @param visitNumber Default: -1
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
#' @importFrom plyr rename ddply
import_SCIDpkpdDataCombo_oldTAD <- function(dataFile,
                                            PKsheets            = c("Blood levels 1", "Blood levels 2"),
                                            PKranges            = NULL,
                                            PKmapping,
                                            PKlloq              = c(0.001,0.001),   # Needs to be in ug/mL.
                                            PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ", "BLOQ"),
                                            PKfactor            = c(1e-3,1e-3),
                                            PDsheet             = "Parasitemia",
                                            PDrange             = NULL,
                                            PDmapping,
                                            PDlloq              = 0.01,
                                            PDlloqIdentifier    = c("BLQ", "<0.01"),
                                            PDname              = "Parasitemia Total",
                                            study               = NULL,
                                            studyID             = NULL,
                                            DayOfFirstDrugAdmin = 1,
                                            centerNumber        = -1,
                                            centerName          = "",
                                            visitNumber         = -1
) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Final Column Names ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  colNames <- c("IGNORE"    , "STUDY"     , "USUBJID"   , "GROUP"    , "SUBJECT",
                "CENTER"    , "CENTERNAME", "COMPOUND"  , "TRTNAME"  , "VISIT"  ,
                "TIME"      , "NT"        , "TIMEUNIT"  , "TYPENAME" , "NAME"   ,
                "VALUE"     , "VALUETXT"  , "UNIT"      , "LLOQ"     , "ROUTE"  ,
                "DOSELEVEL1", "DOSEMULT1" , "DOSELEVEL2", "DOSEMULT2")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Input Variables ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check is the PK column have the appropriate information:
  reqColsPK <- c("SUBJECT", "GROUP", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2", "NT", "VALUE", "COMPOUND")
  if (!all(reqColsPK %in% PKmapping))
    stop("Not all required PK columns supplied in mapping: see import_SCIDpkpdDataCombo.R")

  # Check is the PD column have the appropriate information:
  reqColsPD <- c("SUBJECT", "GROUP", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2", "NT", "VALUE", "huErythro", "COMPOUND")
  if (!all(reqColsPD %in% PDmapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdDataCombo.R")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PK data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load both PK datasets:
  dataPKraw1 <- readxl::read_excel(dataFile, PKsheets[1], range = PKranges[1], na = c("NA", "", "."))
  dataPKraw2 <- readxl::read_excel(dataFile, PKsheets[2], range = PKranges[2], na = c("NA", "", "."))

  # Remove escape characters from column names:
  names(dataPKraw1) <- aux_removeEscapeChar(names(dataPKraw1))
  names(dataPKraw2) <- aux_removeEscapeChar(names(dataPKraw2))

  # Apply mapping:
  dataPK1 <- plyr::rename(dataPKraw1, replace = PKmapping, warn_missing = FALSE)
  dataPK2 <- plyr::rename(dataPKraw2, replace = PKmapping, warn_missing = FALSE)

  # Remove not matched columns:
  dataPK1 <- dataPK1[, intersect(names(dataPK1),PKmapping)]
  dataPK2 <- dataPK2[, intersect(names(dataPK2),PKmapping)]

  # Remove NA rows: NOT SURE IT IS A GOOD IDEA????? as sometimes they forget to add "<LLOQ"
  dataPK1 <- dataPK1[!is.na(dataPK1$VALUE),]
  dataPK2 <- dataPK2[!is.na(dataPK2$VALUE),]

  # Get compound names:
  compound1 <- grep("+",unique(dataPK1$COMPOUND), invert = TRUE, value = TRUE, fixed = TRUE)
  compound2 <- grep("+",unique(dataPK2$COMPOUND), invert = TRUE, value = TRUE, fixed = TRUE)

  # Define compound specific columns:
  dataPK1 <- within(dataPK1, {
    NAME      <- paste0(compound1, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq[1]/2/PKfactor[1], VALUE)) * PKfactor[1]
    LLOQ      <- PKlloq[1]
  })
  dataPK2 <- within(dataPK2, {
    NAME      <- paste0(compound2, " Blood Concentration")
    VALUE     <- as.numeric(ifelse(VALUE %in% PKlloqIdentifier, PKlloq[2]/2/PKfactor[2], VALUE)) * PKfactor[2]
    LLOQ      <- PKlloq[2]
  })

  # Check that numeric columns are numeric:
  numcols <- c("DOSELEVEL1", "DOSELEVEL2", "DOSEMULT1", "DOSEMULT2", "TIME", "NT", "LLOQ")
  dataPK1 <- cbind(dataPK1[, setdiff(names(dataPK1), numcols)], sapply(dataPK1[, intersect(names(dataPK1), numcols)], as.numeric))
  dataPK2 <- cbind(dataPK2[, setdiff(names(dataPK2), numcols)], sapply(dataPK2[, intersect(names(dataPK2), numcols)], as.numeric))

  # Remove potential NA entries in doselevel and number of doses: Set to 0
  dataPK1 <- within(dataPK1, {
    DOSELEVEL2 <- ifelse(is.na(DOSELEVEL2), 0, DOSELEVEL2)
    DOSEMULT2  <- ifelse(is.na(DOSEMULT2), 0, DOSEMULT2)
  })
  dataPK2 <- within(dataPK2, {
    DOSELEVEL1 <- ifelse(is.na(DOSELEVEL1), 0, DOSELEVEL1)
    DOSEMULT1  <- ifelse(is.na(DOSEMULT1), 0, DOSEMULT1)
  })

  # Concatenate PK datasets:
  dataPK <- rbind(dataPK1, dataPK2)

  # (Re-)define compound (order as defined by 1 and 2 in this script) and treatment:
  dataPK$COMPOUND <- paste0(compound1,"+",compound2)
  dataPK <- within(dataPK,{str1    <- ifelse(DOSEMULT1 == 1, paste0(DOSELEVEL1, "mg/kg",compound1), paste0(DOSEMULT1,"x",DOSELEVEL1, "mg/kg",compound1))
  str2    <- ifelse(DOSEMULT2 == 1, paste0(DOSELEVEL2, "mg/kg",compound2), paste0(DOSEMULT2,"x",DOSELEVEL2, "mg/kg",compound2))
  TRTNAME <- ifelse(DOSELEVEL1 >  0 & DOSELEVEL2 == 0, str1, NA)
  TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 > 0, str2, TRTNAME)
  TRTNAME <- ifelse(DOSELEVEL1 >  0 & DOSELEVEL2 > 0, paste0(str1, "+", str2), TRTNAME)
  str1    <- str2 <- NULL
  }
  )

  # use "studyID" if possible, otherwise "study":
  dataPK$STUDYID <- NaN
  if (!is.null(studyID)){
    dataPK$STUDYID <- dataPK[[studyID]]
  } else{
    dataPK$STUDYID <- study
  }

  # Define common columns for both PK:
  dataPK <- within(dataPK, {
    STUDY     <- STUDYID
    TIMEUNIT  <- rep("hours", length(STUDY))
    TIME      <- NT
    TYPENAME  <- "PK"
    VALUETXT  <- NA
    UNIT      <- "ug/mL"
    ROUTE     <- NA
    IGNORE    <- NA
  })
  dataPK$USUBJID <- aux_createUSUBJID(dataPK)

  # Adjust DOSEMULT1 & DOSEMULT2 if DOSELEVEL=0:
  dataPK$DOSEMULT1 <- ifelse(dataPK$DOSELEVEL1==0, 0, dataPK$DOSEMULT1)
  dataPK$DOSEMULT2 <- ifelse(dataPK$DOSELEVEL2==0, 0, dataPK$DOSEMULT2)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PD data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load data:
  dataPDraw <- readxl::read_excel(dataFile, PDsheet, range = PDrange, na = c("NA", "", "."))

  # Remove escape characters from column names:
  names(dataPDraw) <- aux_removeEscapeChar(names(dataPDraw))

  # Apply matching:
  dataPD <- plyr::rename(dataPDraw, replace = PDmapping)

  # Remove not matched columns:
  dataPD <- dataPD[, intersect(names(dataPD),PDmapping)]

  # Remove all controls except vehicle:
  idx_keep = (grepl(compound1, dataPD$COMPOUND) |
                grepl(compound2, dataPD$COMPOUND) |
                (dataPD$DOSELEVEL1==0 & dataPD$DOSELEVEL2==0))
  dataPD <- dataPD[idx_keep,]

  # use "studyID" if possible, otherwise "study":
  dataPD$STUDYID <- NaN
  if (!is.null(studyID)){
    dataPD$STUDYID <- dataPD[[studyID]]
  } else{
    dataPD$STUDYID <- study
  }

  # Check that numeric columns are numeric:
  numcols <- c("DOSELEVEL1", "DOSELEVEL2", "DOSEMULT1", "DOSEMULT2", "TIME", "NT", "LLOQ")
  dataPD <- cbind(dataPD[, setdiff(names(dataPD), numcols)], sapply(dataPD[, intersect(names(dataPD), numcols)], as.numeric))

  # Remove potential NA entries in doselevel and number of doses: Set to 0
  dataPD <- within(dataPD, {
    DOSELEVEL2 <- ifelse(is.na(DOSELEVEL2), 0, DOSELEVEL2)
    DOSEMULT2  <- ifelse(is.na(DOSEMULT2), 0, DOSEMULT2)
    DOSELEVEL1 <- ifelse(is.na(DOSELEVEL1), 0, DOSELEVEL1)
    DOSEMULT1  <- ifelse(is.na(DOSEMULT1), 0, DOSEMULT1)
  })

  # Remove records with no entry for value:
  # Should remove also the gap between experiment parts
  # NOT SURE IT IS A GOOD IDEA TO REMOVE THE NA?????
  dataPD <- dataPD[!is.na(dataPD$VALUE), ]

  # Add Common Clumn and treat data<LLOQ and dead mices:
  dataPD <- within(dataPD, {
    STUDY     <- STUDYID
    TIMEUNIT  <- rep("hours", length(STUDY))
    NT        <- (NT-DayOfFirstDrugAdmin)*24
    TIME      <- NT
    TYPENAME  <- "Efficacy"
    NAME      <- PDname
    IGNORE    <- ifelse(VALUE %in% c("Dead"), "Dead", NA)
    VALUE     <- ifelse(VALUE %in% c("Dead"), NA, VALUE)
    VALUE     <- ifelse(VALUE %in% PDlloqIdentifier, 0, VALUE)
    VALUE     <- as.numeric(VALUE)
    VALUE     <- ifelse(VALUE <= PDlloq, PDlloq/2, VALUE)
    VALUETXT  <- NA
    UNIT      <- "percent"
    ROUTE     <- NA
    LLOQ      <- PDlloq
    huErythro <- ifelse(huErythro %in% c("Dead"), NA, huErythro)
    huErythro <- as.numeric(huErythro)
  })
  dataPD$USUBJID   <- aux_createUSUBJID(dataPD)

  # Normalize parasitemia to human erythrocyte percent:
  dataPD$VALUE <- ifelse(dataPD$VALUE == PDlloq/2,PDlloq/2, dataPD$VALUE /dataPD$huErythro * 100)


  # (Re-)define compound and treatment:
  dataPD$COMPOUND <- paste0(compound1,"+",compound2)
  dataPD <- within(dataPD,{str1 <- ifelse(DOSEMULT1 == 1, paste0(DOSELEVEL1, "mg/kg",compound1), paste0(DOSEMULT1,"x",DOSELEVEL1, "mg/kg",compound1))
  str2 <- ifelse(DOSEMULT2 == 1, paste0(DOSELEVEL2, "mg/kg",compound2), paste0(DOSEMULT2,"x",DOSELEVEL2, "mg/kg",compound2))
  TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 == 0, "Vehicle", NA)
  TRTNAME <- ifelse(DOSELEVEL1 >  0 & DOSELEVEL2 == 0, str1, TRTNAME)
  TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 >  0, str2, TRTNAME)
  TRTNAME <- ifelse(DOSELEVEL1 >  0 & DOSELEVEL2 >  0, paste0(str1, "+", str2), TRTNAME)
  str1    <- str2 <- NULL
  }
  )

  # Adjust DOSEMULT1 & DOSEMULT2 if DOSELEVEL=0:
  dataPD$DOSEMULT1 <- ifelse(dataPD$DOSELEVEL1==0, 0, dataPD$DOSEMULT1)
  dataPD$DOSEMULT2 <- ifelse(dataPD$DOSELEVEL2==0, 0, dataPD$DOSEMULT2)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Erythrocyte Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dataHE <- within(dataPD,{VALUE    <- huErythro
  TYPENAME <- "Vital Signs"
  NAME     <- "Human Erythrocytes"
  }
  )


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Dosing Records ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Define columns to retain from PK and PD dataset: Only keep one line per subject:
  colDos  <- c("USUBJID", "STUDY", "GROUP", "TRTNAME", "SUBJECT", "COMPOUND", "DOSELEVEL1", "DOSEMULT1", "DOSELEVEL2", "DOSEMULT2", "TIMEUNIT")
  dataDos <- unique(rbind(dataPK[, colDos],dataPD[, colDos]))

  # Check that there is only one line per subject:
  if (any(duplicated(dataDos$USUBJID)))
    stop("Dosing information not unique for individuals: Check dataset")

  # Create Dose events for each subject:
  dataDos <- plyr::ddply(dataDos, ~USUBJID, function(x){
    if (x$DOSEMULT1 > 0){
      out1 <- x[rep(1,x$DOSEMULT1),]
      out1 <- within(out1,{TIME      <- seq(0, by = 24, length.out = as.numeric(as.character(x$DOSEMULT1)))
      NT        <- TIME
      TYPENAME  <- "Dose"
      NAME      <- paste0(compound1, " Dose")
      VALUE     <- as.numeric(as.character(x$DOSELEVEL1))
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- "ORAL"
      IGNORE    <- NA
      }
      )
    }else{
      out1 <- NULL
    }

    if (x$DOSEMULT2 > 0){
      out2 <- x[rep(1,x$DOSEMULT2),]
      out2 <- within(out2, {TIME      <- seq(0, by = 24, length.out = as.numeric(as.character(x$DOSEMULT2)))
      NT        <- TIME
      TYPENAME  <- "Dose"
      NAME      <- paste0(compound2, " Dose")
      VALUE     <- as.numeric(as.character(x$DOSELEVEL2))
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- "ORAL"
      IGNORE    <- NA
      }
      )
    }else{
      out2 <- NULL
    }

    out <- rbind(out1,out2)
  })


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Concatenate All Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Add CENTER & VISIT column:
  dataPK$CENTER      <- centerNumber
  dataPK$CENTERNAME  <- centerName
  dataPK$VISIT       <- visitNumber
  dataPD$CENTER      <- centerNumber
  dataPD$CENTERNAME  <- centerName
  dataPD$VISIT       <- visitNumber
  dataHE$CENTER      <- centerNumber
  dataHE$CENTERNAME  <- centerName
  dataHE$VISIT       <- visitNumber
  dataDos$CENTER     <- centerNumber
  dataDos$CENTERNAME <- centerName
  dataDos$VISIT      <- visitNumber

  # Bind all data:
  data <- rbind(dataPK[,colNames],
                dataPD[,colNames],
                dataHE[,colNames],
                dataDos[,colNames]
  )

  # Order:
  data <- data[order(data$USUBJID, data$TIME),]

  # Output:
  data
}

#' import_SCIDpkpdDataCombo_TAD
#' Import SCID Combo data from TAD
#'
#' Import original SCID combination data from TAD and convert it to a single structured data format
#'
#' Import original SCID combination data from TAD and convert it to a single structured data format
#' ready to be imported by IQRTool (IQRdataGENERAL) or pre-ready used for MONOLIX or NONMEM
#' It includes positive controls data
#' Uses function [readxl::read_excel]
#'
#' @param dataFile Original TAD excel datafile
#' @param Compound1 List containining user prefered compound name and MMV compound name. One of the two must be the one used in the originla data provided by TAD (cpd1  and cpd2 by alphabetic order)
#' @param Compound2 List containining user prefered compound name and MMV compound name. One of the two must be the one used in the originla data provided by TAD (cpd1  and cpd2 by alphabetic order)
#' @param PKsheets Name of the single sheet recording PK data
#' @param PKranges Cells range of the PKsheet to read (eg A1:C20)
#' @param PKmapping Character vector mapping PKsheet headers (left) to desired name in converted data format (right) (eg c("StudyID"= "STUDY", "Concentration" = "VALUE"))
#' @param PKlloq Numeric vector of PK LLOQ for compound 1 and compound 2 respectively
#' @param PKlloqIdentifier Character vector for PKlloq identifier (eg c("BLQ", "<0.01"))
#' @param PKfactor Numeric vector for compound 1 and compound 2 respectively, converting factor for PK data to get desired unit (eg c(1e-3,1e-3) from ng/mL to ug/mL for both compounds) )
#' @param PDsheet Name of the sheet recording PD data
#' @param PDrangeDoseTimeRange
#' @param PDmapping Character vector mapping PDsheet headers (left) to desired name in converted final dataset (right) (eg c("StudyID"= "STUDY", "Parasitemia" = "VALUE"))
#' @param PDlloq Numeric value of PD LLOQ
#' @param PDlloqIdentifier Numeric value of PD LLOQ
#' @param PDname Character, define the name of the PD variable to be recorded in the converted final dataset
#' @param DayOfFirstDrugAdmin Day of first drug administration as recorded in the original dataset
#' @param centerNumber Number corresponding to center name, see 'list_Center()'
#' @param centerName Center name, see 'list_Center()'
#' @param visitNumber Numeric Value, visit Number. Default = -1 when not relevant
#' @param intervaldose Character providing dosing interval. only to provide if dosing different of once a day.Used to define TRTNAME (eg "bid")
#' @param route Character, route of administration. Only useful if not provided in the original dataset
#' @param DoseTimeSheet Name of the sheet recording dosing time of drug administration for both compounds
#' @param DoseTimeRange Cells range of the DoseTimeSheet to read (eg A1:C20)
#' @param DoseTimeMapping Character vector mapping DoseTimesheet headers (left) to desired name in converted final dataset (right) (eg c("StudyID"="STUDY", "DosingTime"="NT"))
#' @param DoseSheet Name of the sheet recording drugs dosing for both compounds
#' @param DoseRange Cells range of the DoseSheet to read (eg A1:C20)
#' @param DoseMapping Character vector mapping DoseSheet headers (left) to desired name in converted final dataset (right) (eg c("StudyID"="STUDY", "Route" = "ROUTE"))
#' @param positiveQC Character vector identifying the name of positive control (must matchupper or lower case as in original dataset)
#'
#' @importFram reshape2 dcast melt
#' @export
#' @seealso [readxl], [plyr], [dplyr], [reshape2]
#' @family Data Preparation
#' @author Aline Fuchs MMV
import_SCIDpkpdDataCombo_TAD <- function(dataFile,
                                         Compound1           = NULL,
                                         Compound2           = NULL,
                                         PKsheets             = "DrugConcentrationData",
                                         PKranges             = NULL,
                                         PKmapping,
                                         PKlloq              = c(1,1),
                                         PKlloqIdentifier    = c("BLQ","<LLOQ","< LLOQ", "BLOQ"),
                                         PKfactor            = c(1e-3,1e-3),
                                         PDsheet             = "ParasitemiaData",
                                         PDrange             = NULL,
                                         PDmapping,
                                         PDlloq              = 0.01,
                                         PDlloqIdentifier    = c("BLQ", "<0.01", "< 0,01", "< 0.01"),
                                         PDname              = "Parasitemia Total",
                                         DayOfFirstDrugAdmin = 1,
                                         centerNumber        = -1,
                                         centerName          = "",
                                         visitNumber         = -1,
                                         intervaldose        = NULL, # only to spell TRTNAME; only to provide if dosing different of once a day
                                         route               = NULL,
                                         DoseTimeSheet       = "DosingTimeTable",
                                         DoseTimeRange       = NULL,
                                         DoseTimeMapping,
                                         DoseSheet           = "DrugTreatmentTable",
                                         DoseRange           = NULL,
                                         DoseMapping,
                                         positiveQC

) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Final Column Names ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  colNames <- c("IGNORE"    , "STUDY"     , "USUBJID"   , "GROUP"    , "SUBJECT",
                "CENTER"    , "CENTERNAME", "COMPOUND"  , "TRTNAME"  , "VISIT"  ,
                "TIME"      , "NT"        , "TIMEUNIT"  , "TYPENAME" , "NAME"   ,
                "VALUE"     , "VALUETXT"  , "UNIT"      , "LLOQ"     , "ROUTE"  ,
                "DOSELEVEL1", "DOSEMULT1" , "DOSELEVEL2", "DOSEMULT2")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Input Variables ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check if the PK column have the appropriate information:
  reqColsPK <- c("SUBJECT", "STUDY", "COMPOUND", "NT", "VALUE")
  if (!all(reqColsPK %in% PKmapping))
    stop("Not all required PK columns supplied in mapping: see import_SCIDpkpdDataCombo.R")

  # Check if the PD column have the appropriate information:
  reqColsPD <- c("SUBJECT", "STUDY", "NT", "VALUE", "huErythro")
  if (!all(reqColsPD %in% PDmapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdDataCombo.R")

  # Check if the Dose Time column have the appropriate information:
  reqColsDoseTime <- c("SUBJECT", "STUDY", "COMPOUND", "NT")
  if (!all(reqColsDoseTime %in% DoseTimeMapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdDataCombo.R")

  # Check if the Dose column have the appropriate information:
  reqColsDose <- c("SUBJECT", "STUDY", "COMPOUND", "DOSELEVEL")
  if (!all(reqColsDose %in% DoseMapping))
    stop("Not all required PD columns supplied in mapping: see import_SCIDpkpdDataCombo.R")

  # Check is the number of dose is defined:
  if (!("DOSEMULT" %in% DoseMapping))
    stop("Number of doses not to be supplied in Dose data set.")

  # Check is route of administration is defined:
  if (is.null(route) & !("ROUTE" %in% DoseMapping))
    stop("If route of administration is not supplied in Dose data set, needs to be defined by input argument.")


  # Check whether dosing interval is provided
  if (is.null(DoseTimeSheet))
    warning("Dosing time sheet not provided, once daily dosing is assumed.")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle Dosing data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # -- Load Dose dataset:
  dataDoseraw <- readxl::read_excel(dataFile, DoseSheet, range = DoseRange)

  # Remove escape characters from column names
  names(dataDoseraw) <- aux_removeEscapeChar(names(dataDoseraw))

  # Apply mapping:
  dataDose <- plyr::rename(dataDoseraw, replace = DoseMapping)
  dataDose <- as.data.frame(dataDose)

  # Remove not matched columns:
  dataDose <- dataDose[, intersect(names(dataDose),DoseMapping)]

  # Identify Treatment received by SUBJECT
  cpdAll <- unique(dataDose$COMPOUND)

  # Get Compounds as names in the original data
  cpd1 <- cpdAll[cpdAll %in% c(Compound1$Name,Compound1$MMVname)]
  cpd2 <- cpdAll[cpdAll %in% c(Compound2$Name,Compound2$MMVname)]

  # Identify Mice that receive compounds which we are  NOT interested
  idx_rm <-  dataDose$SUBJECT[!dataDose$COMPOUND %in% c("NONE",positiveQC,cpd1,cpd2)]
  # Identify Subjects which received only compound of interest and select them only
  idx_keep <-  setdiff(dataDose$SUBJECT,idx_rm)
  dataDose <- dataDose[dataDose$SUBJECT %in% idx_keep,]

  # Get Appropriate DoseLevel
  DL <- reshape2::dcast(dataDose,SUBJECT+STUDY+ROUTE~COMPOUND,value.var = "DOSELEVEL")
  # Doselevel 1
  DL1 <- reshape2::melt(DL,
                        id.vars=c(names(DL)[!names(DL) %in% cpdAll]),
                        measure.vars = c("NONE",positiveQC,cpd1))
  names(DL1) <- c(names(DL1)[!names(DL1) %in% c("variable","value")],"COMPOUND","DOSELEVEL1")
  DL1<-DL1[complete.cases(DL1), ]
  # Doselevel 2
  DL2 <- reshape2::melt(DL,
                        id.vars=c(names(DL)[!names(DL) %in% cpdAll]),
                        measure.vars = c(cpd2))
  names(DL2) <- c(names(DL2)[!names(DL2) %in% c("variable","value")],"COMPOUND","DOSELEVEL2")
  DL2 <- DL2[complete.cases(DL2), ]

  # Get Appropriate DoseMulti:
  DM <- reshape2::dcast(dataDose,SUBJECT+STUDY+ROUTE~COMPOUND,value.var = "DOSEMULT")
  # DoseMulti 1
  DM1 <- reshape2::melt(DM,
                        id.vars=c(names(DM)[!names(DM) %in% cpdAll]),
                        measure.vars = c("NONE",positiveQC,cpd1))
  names(DM1) <- c(names(DM1)[!names(DM1) %in% c("variable","value")],"COMPOUND","DOSEMULT1")
  DM1<-DM1[complete.cases(DM1), ]
  # DoseMulti 2
  DM2 <- reshape2::melt(DM,
                        id.vars=c(names(DM)[!names(DM) %in% cpdAll]),
                        measure.vars = c(cpd2))
  names(DM2) <- c(names(DM2)[!names(DM2) %in% c("variable","value")],"COMPOUND","DOSEMULT2")
  DM2<-DM2[complete.cases(DM2), ]

  # Join Doselevel1 and DoseMult1
  cpd1Dose <- full_join(DL1,DM1,by=c("SUBJECT","STUDY","ROUTE","COMPOUND"))
  cpd2Dose <- full_join(DL2,DM2,by=c("SUBJECT","STUDY","ROUTE","COMPOUND"))

  # -- Join dosing info:
  dataDos <- dplyr::full_join(cpd1Dose, cpd2Dose,by=c("SUBJECT","STUDY","ROUTE"),suffix=c("1","2"))
  dataDos$COMPOUND1 <- as.character(dataDos$COMPOUND1)
  dataDos$COMPOUND2 <- as.character(dataDos$COMPOUND2)


  # Check that numeric columns are numeric:
  numcols <- c("DOSELEVEL1", "DOSEMULT1","DOSELEVEL2", "DOSEMULT2","TIME", "NT")
  dataDos <- cbind(dataDos[, setdiff(names(dataDos), numcols)], sapply(dataDos[, intersect(names(dataDos), numcols)], as.numeric))

  # Remove potential NA entries in doselevel and number of doses: Set to 0
  dataDos <- within(dataDos, {
    DOSELEVEL1 <- ifelse(is.na(DOSELEVEL1), 0, DOSELEVEL1)
    DOSELEVEL2 <- ifelse(is.na(DOSELEVEL2), 0, DOSELEVEL2)
  })

  # handle multiple dose
  dataDos$DOSEMULT1 <- ifelse(dataDos$DOSELEVEL1==0, 0, dataDos$DOSEMULT1)
  dataDos$DOSEMULT2 <- ifelse(dataDos$DOSELEVEL2==0, 0, dataDos$DOSEMULT2)

  # (Re-)define compound (order as defined by 1 and 2 in this script) and treatment:
  dataDos$COMPOUND <- paste0(cpd1,"+",cpd2)
  dataDos$COMPOUND <- ifelse(dataDos$COMPOUND1 == positiveQC & !is.na(dataDos$COMPOUND1), positiveQC,dataDos$COMPOUND)


  dataDos <- within(dataDos,{
    str1    <- ifelse(DOSEMULT1 == 1 & COMPOUND!=positiveQC , paste0(DOSELEVEL1, "mg/kg",cpd1), paste0(DOSEMULT1,"x",DOSELEVEL1, "mg/kg",cpd1))
    str2    <- ifelse(DOSEMULT2 == 1 & COMPOUND!=positiveQC , paste0(DOSELEVEL2, "mg/kg",cpd2), paste0(DOSEMULT2,"x",DOSELEVEL2, "mg/kg",cpd2))
    str3    <- ifelse(DOSEMULT1 == 1 & COMPOUND==positiveQC , paste0(DOSELEVEL1, "mg/kg",positiveQC), paste0(DOSEMULT1,"x",DOSELEVEL1, "mg/kg",positiveQC))
    TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 == 0 & COMPOUND!=positiveQC, "Vehicle", NA)
    TRTNAME <- ifelse(DOSELEVEL1 >  0 & DOSELEVEL2 == 0 & COMPOUND!=positiveQC, str1, TRTNAME)
    TRTNAME <- ifelse(DOSELEVEL1 == 0 & DOSELEVEL2 > 0 & COMPOUND!=positiveQC, str2, TRTNAME)
    TRTNAME <- ifelse(DOSELEVEL1 >  0 & DOSELEVEL2 > 0 & COMPOUND!=positiveQC, paste0(str1, "+", str2), TRTNAME)
    TRTNAME <- ifelse(COMPOUND==positiveQC ,  str3, TRTNAME)
    str1    <- str2 <- str3 <- NULL
  }
  )

  # Define Group TRT: (include positive control; in the new format group should never be a column in any sheet actually)
  if (!("GROUP" %in% names(dataDos))) {
    groupInfo <- unique(dataDos[, c("COMPOUND","DOSELEVEL1", "DOSEMULT1","DOSELEVEL2", "DOSEMULT2")])
    groupInfo <- rbind(
      groupInfo[groupInfo$DOSELEVEL1 == 0 & groupInfo$DOSELEVEL2 == 0,][with(groupInfo[groupInfo$DOSELEVEL1 == 0 & groupInfo$DOSELEVEL2 == 0,], order(DOSEMULT1,DOSEMULT2)),],
      groupInfo[groupInfo$DOSELEVEL1 == 0 & groupInfo$DOSELEVEL2 > 0,][with(groupInfo[groupInfo$DOSELEVEL1 == 0 & groupInfo$DOSELEVEL2 > 0,]  , order(DOSEMULT1,DOSELEVEL1,DOSEMULT2,DOSELEVEL2)),],
      groupInfo[groupInfo$DOSELEVEL1 > 0 & groupInfo$DOSELEVEL2  == 0,][with(groupInfo[groupInfo$DOSELEVEL1 > 0 & groupInfo$DOSELEVEL2 == 0,] , order(DOSEMULT1,DOSELEVEL1,DOSEMULT2,DOSELEVEL2)),],
      groupInfo[groupInfo$DOSELEVEL1 > 0 & groupInfo$DOSELEVEL2  > 0,][with(groupInfo[groupInfo$DOSELEVEL1 > 0 & groupInfo$DOSELEVEL2 > 0,] , order(DOSEMULT1,DOSELEVEL1,DOSEMULT2,DOSELEVEL2)),]
    )
    groupInfo$GROUP <- 1:dim(groupInfo)[1]
    dataDos <- plyr::join(dataDos, groupInfo)
  } else {
    groupInfo <- unique(dataDos[, c("COMPOUND","DOSELEVEL1", "DOSEMULT1","DOSELEVEL2", "DOSEMULT2", "GROUP")])
  }

  # Add UBSUBJID
  dataDos$USUBJID <- aux_createUSUBJID(dataDos)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Time Dosing Records ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # -- Load DoseTime dataset:
  dataDoseTimeraw <- readxl::read_excel(dataFile, DoseTimeSheet, range = DoseTimeRange)

  # Remove escape characters from column names
  names(dataDoseTimeraw) <- aux_removeEscapeChar(names(dataDoseTimeraw))

  # Apply mapping:
  dataDoseTime <- plyr::rename(dataDoseTimeraw, replace = DoseTimeMapping)
  dataDoseTime = as.data.frame(dataDoseTime)

  # Remove not matched columns:
  dataDoseTime <- dataDoseTime[, intersect(names(dataDoseTime),DoseTimeMapping)]

  # Keep only IDs of interest:
  dataDoseTime <- dataDoseTime[dataDoseTime$SUBJECT %in% idx_keep,]

  # -- Join dosing info:
  dataDosTime1 <- dplyr::left_join(dataDoseTime[dataDoseTime$COMPOUND %in% c("NONE",positiveQC,cpd1),], dataDos,
                                   by = c("SUBJECT" = "SUBJECT", "STUDY" = "STUDY", "COMPOUND" = "COMPOUND1"))
  dataDosTime1$NAME      <- with(dataDosTime1,ifelse(COMPOUND!="NONE"&COMPOUND!=positiveQC,paste0(cpd1, " Dose"),paste0(COMPOUND, " Dose")))
  dataDosTime1$VALUE     <- as.numeric(as.character(dataDosTime1$DOSELEVEL1))

  dataDosTime2 <- dplyr::left_join(dataDoseTime[dataDoseTime$COMPOUND %in% c(cpd2),], dataDos,
                                   by = c("SUBJECT" = "SUBJECT", "STUDY" = "STUDY", "COMPOUND" = "COMPOUND2"))
  dataDosTime2$NAME      <- with(dataDosTime2,ifelse(COMPOUND!="NONE"&COMPOUND!=positiveQC,paste0(cpd2, " Dose"),paste0(COMPOUND, " Dose")))
  dataDosTime2$VALUE     <- as.numeric(as.character(dataDosTime2$DOSELEVEL2))

  # Bind
  dataDosTime <- rbind.fill(dataDosTime1,dataDosTime2)

  # Reduce and rename
  colToKeep             <- c("SUBJECT","STUDY","USUBJID","GROUP","NT","ROUTE","COMPOUND.y","TRTNAME","VALUE" ,"NAME","DOSELEVEL1","DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")
  dataDosTime           <- dataDosTime[,colToKeep]
  colnames(dataDosTime) <- c("SUBJECT","STUDY","USUBJID","GROUP","NT","ROUTE","COMPOUND" ,"TRTNAME","VALUE" ,"NAME","DOSELEVEL1","DOSEMULT1", "DOSELEVEL2", "DOSEMULT2")

  # Check that there is only one line per subject:
  if (any(duplicated(dataDosTime[c("USUBJID","NT","TRTNAME","NAME")])))
    stop("Dosing information not unique for individuals: Check dataset")


  # Add Dose columns:
  dataDosTime <- plyr::ddply(dataDosTime, ~USUBJID, function(x) {
    out <- within(x,{
      TIMEUNIT  <- rep("hours", length(STUDY))
      TIME      <- NT
      TYPENAME  <- "Dose"
      VALUETXT  <- NA
      UNIT      <- "mg/kg"
      LLOQ      <- NA
      ROUTE     <- ifelse("ROUTE" %in% names(dataDos),  toupper(ROUTE), toupper(route))
      IGNORE    <- NA
    }
    )
  })


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PK data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load both PK datasets:
  dataPKraw <- readxl::read_excel(dataFile, PKsheets, range = PKranges, na = c("NA", "", "."))

  # Remove escape characters from column names:
  names(dataPKraw) <- aux_removeEscapeChar(names(dataPKraw))

  # Apply mapping:
  dataPK <- plyr::rename(dataPKraw, replace = PKmapping, warn_missing = FALSE)

  # Get only compound of interest
  # dataPK <- dataPK[dataPK$COMPOUND %in% c(cpd1,cpd2),]
  # Get only subject of interest
  dataPK <- dataPK[dataPK$SUBJECT %in% idx_keep,]

  # Remove not matched columns:
  dataPK <- dataPK[, intersect(names(dataPK),PKmapping)]

  # Remove NA rows: NOT SURE IT IS A GOOD IDEA????? as sometimes they forget to add "<LLOQ"
  dataPK <- dataPK[!is.na(dataPK$VALUE),]

  # Get rid of concentration at TIME = 0; no measurement
  dataPK <- dataPK[dataPK$NT!=0,]

  # Define extra columns for PK:
  dataPK <- mutate(dataPK,
                   LLOQ      = ifelse(COMPOUND==cpd1, PKlloq[1], PKlloq[2]),
                   PKFACTOR  = ifelse(COMPOUND==cpd1,PKfactor[1], PKfactor[2]),
                   TIMEUNIT  = rep("hours", length(STUDY)),
                   TIME      = NT,
                   TYPENAME  = "PK",
                   NAME      = paste0(COMPOUND, " Blood Concentration"),
                   VALUE     = as.numeric(as.character(ifelse(VALUE %in% PKlloqIdentifier, LLOQ/2/PKFACTOR, VALUE))) * PKFACTOR,
                   VALUETXT  = NA,
                   UNIT      = "ug/mL",
                   ROUTE     = NA,
                   IGNORE    = NA
  )

  # Check that numeric columns are numeric:
  numcols <- c("TIME", "NT", "LLOQ","VALUE")
  dataPK <- cbind(dataPK[, setdiff(names(dataPK), numcols)], sapply(dataPK[, intersect(names(dataPK), numcols)], as.numeric))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle PD data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Load data:
  dataPDraw <- readxl::read_excel(dataFile, PDsheet, range = PDrange, na = c("NA", "", "."))

  # Remove escape characters from column names:
  names(dataPDraw) <- aux_removeEscapeChar(names(dataPDraw))

  # Apply matching:
  dataPD <- plyr::rename(dataPDraw, replace = PDmapping)

  # Remove not matched columns:
  dataPD <- dataPD[, intersect(names(dataPD),PDmapping)]

  # Get only subject of interest
  dataPD <- dataPD[dataPD$SUBJECT %in% idx_keep,]

  # Remove records with no entry for value:
  # Should remove also the gap between experiment parts
  # NOT SURE IT IS A GOOD IDEA TO REMOVE THE NA?????
  dataPD <- dataPD[!is.na(dataPD$VALUE), ]

  # Add Common column and treat data<LLOQ and dead mices:
  dataPD <- within(dataPD, {
    TIMEUNIT  <- rep("hours", length(STUDY))
    NT        <- (NT-DayOfFirstDrugAdmin)*24
    TIME      <- NT
    TYPENAME  <- "Efficacy"
    NAME      <- PDname
    IGNORE    <- ifelse(VALUE %in% c("Dead"), "Dead", NA)
    VALUE     <- ifelse(VALUE %in% c("Dead"), NA, VALUE)
    VALUE     <- ifelse(VALUE %in% PDlloqIdentifier, 0, VALUE)
    VALUE     <- as.numeric(as.character(VALUE))
    VALUE     <- ifelse(VALUE <= PDlloq, PDlloq/2, VALUE)
    VALUETXT  <- NA
    UNIT      <- "percent"
    ROUTE     <- NA
    LLOQ      <- PDlloq
    huErythro <- ifelse(huErythro %in% c("Dead"), NA, huErythro)
    huErythro <- as.numeric(huErythro)
  })

  # Check that numeric columns are numeric:
  numcols <- c("TIME", "NT", "LLOQ", "VALUE")
  dataPD <- cbind(dataPD[, setdiff(names(dataPD), numcols)], sapply(dataPD[, intersect(names(dataPD), numcols)], as.numeric))

  # Normalize parasitemia to human erythrocyte percent:
  dataPD$VALUE <- ifelse(dataPD$VALUE == PDlloq/2,PDlloq/2, dataPD$VALUE /dataPD$huErythro * 100)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define Erythrocyte Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dataHE <- within(dataPD,{VALUE    <- huErythro
  TYPENAME <- "Vital Signs"
  NAME     <- "Human Erythrocytes"
  }
  )


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Concatenate All Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dataPK <- dataPK[,-which(names(dataPK)          %in% setdiff(names(dataPK),names(dataPD)))]
  dataPD <- dataPD[,-which(names(dataPD)          %in% setdiff(names(dataPD),names(dataPK)))]
  dataHE <- dataHE[,-which(names(dataHE)          %in% setdiff(names(dataHE),names(dataPK)))]
  dataDO <- dataDosTime[,-which(names(dataDosTime)%in% setdiff(names(dataDosTime),names(dataPK)))]

  data <- rbind(dataPK,
                dataPD,
                dataHE,
                dataDO)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Concatenate All Data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  data <- dplyr::left_join(data, unique(dataDosTime[,c("SUBJECT", "STUDY","USUBJID","GROUP","COMPOUND",
                                                       "DOSELEVEL1",  "DOSEMULT1", "DOSELEVEL2",  "DOSEMULT2","TRTNAME")]))

  # Adjust treatment name if dosing different of once a day
  if (!is.null(intervaldose))
    data$TRTNAME <- ifelse(data$TRTNAME=="Vehicle",data$TRTNAME,paste0(data$TRTNAME, " ",intervaldose))

  # Add CENTER & VISIT column:
  data$CENTER      <- centerNumber
  data$CENTERNAME  <- centerName
  data$VISIT       <- visitNumber

  # Order:
  data <- data[order(data$USUBJID, data$TIME),]

  # Select columns
  data <- data[,colNames]

  # Remove Dose for vehicle
  data<-data[!(data$TYPENAME=="Dose" & data$TRTNAME=="Vehicle"),]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get data only for compoun of interest ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Select only for compound of interest + vehicle + positive control
  data <- rbind(data[data$COMPOUND==paste0(cpd1,"+",cpd2),],
                data[data$COMPOUND==positiveQC,])
  data<- data[complete.cases(data[ , "COMPOUND"]),]

  # Replace positiveQC by appropriate compound NAme (positiveQC info is kept in TRTNAME)
  data$COMPOUND <- paste0(cpd1,"+",cpd2)

  # Rename to MMVName
  data$TRTNAME <- gsub(Compound1$Name, Compound1$MMVname, data$TRTNAME)
  data$NAME    <- gsub(Compound1$Name, Compound1$MMVname, data$NAME)
  data$COMPOUND<- gsub(Compound1$Name, Compound1$MMVname, data$COMPOUND)
  data$TRTNAME <- gsub(Compound2$Name, Compound2$MMVname, data$TRTNAME)
  data$NAME    <- gsub(Compound2$Name, Compound2$MMVname, data$NAME)
  data$COMPOUND<- gsub(Compound2$Name, Compound2$MMVname, data$COMPOUND)

  # Get only Study on which the loop is performed
  data <- data[data$STUDY==StudyName_k,]
  # data <- data.frame(lapply(data, function(x) {
  #                     gsub(Compound1$Name, Compound1$MMVname, x)
  #                 }))
  # data <- data.frame(lapply(data, function(x) {
  #   gsub(Compound2$Name, Compound2$MMVname, x)
  # }))


  # Output:
  return(data)

}

#' load_AfricanPediatricMalariaPopulation2to5
#'
#' @description load Malaria population demographic data only for African children 2 to 5 yrs old
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
load_AfricanPediatricMalariaPopulation2to5 <- function(){

  # If not installed as a package:
  if(!("MMVmalaria" %in% .packages())){
    load(file.path(get_MMVmalariaPath(), "data/MalariaPopulation.RData"))
    data <- MalariaPopulation
  }else{
    data <- MMVmalaria::MalariaPopulation
  }

  # Subset of itnerest:
  data <- subset(data, REGION=="Africa" & AGE_years<=5 & AGE_years>2)

  # Output:
  return(data)
}


#' load_MalariaPopulation
#'
#' @description load Malaria population demographic data
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
load_MalariaPopulation <- function(){

  # If not installed as a package:
  if(!("MMVmalaria" %in% .packages())){
    load(file.path(get_MMVmalariaPath(), "data/MalariaPopulation.RData"))
    data <- MalariaPopulation
  }else{
    data <- MMVmalaria::MalariaPopulation
  }

  # Output:
  return(data)
}


#' Load Pediatric Parameters
#'
#' Load default pediatric parameters such as dose adjustement (`Fpediatric`) and
#' volume of blood per kg, per body weight band.
#'
#' @return
#'
#' @export
#'
#' @author Mohammed H. Cherkaoui (MMV)
#' @family Data Preparation
load_PediatricParameters <- function(){

  # Load dataset:
  data <- IQRloadCSVdata(file.path(get_MMVmalariaPath(subdir="inst"),"dataLibrary/pediatricPopulation/Population.csv"))

  # Output:
  return(data)
}

