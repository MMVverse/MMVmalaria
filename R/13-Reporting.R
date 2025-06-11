#' parameter_strToLatex
#'
#' @description

#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Reporting
parameter_strToLaTeX <- function() {

  strToLaTeX <- list(#PK Parameters:
    Fabs0 = "F_{abs,0}",
    Fabs1 = "F_{abs,1}",
    ka    = "k_{a}",
    CL    = "CL",
    CLpm1 = "CL_{p,m1}",
    CLpm2 = "CL_{p,m2}",
    CLpx  = "CL_{p,x}",
    CLm1x = "CL_{m1,x}",
    CLm2x = "CL_{m2,x}",
    CLm3x = "CL_{m3,x}",
    VMAX  = "V_{max}",
    KM    = "K_{m}",
    Vcm1  = "V_{c,m1}",
    Vcm2  = "V_{c,m2}",
    Vcm3  = "V_{c,m3}",
    Vc    = "V_{c}",
    Q1m1  = "Q_{1,m1}",
    Q1m2  = "Q_{1,m2}",
    Q1    = "Q_{1}",
    Vp1m1 = "V_{p1,m1}",
    Vp1m2 = "V_{p1,m2}",
    Vp1   = "V_{p,1}",
    Q2    = "Q_{2}",
    Vp2   = "V_{p,2}",
    Tk0   = "T_{k,0}",
    Tlag1 = "T_{lag,1}",

    # PD Parameters:
    GR     = "GR",
    kGR    = "k_{GR}",
    PLbase = "PL_{base}",
    PLerr  = "PL_{err}",
    EMAX   = "E_{max}",
    EC50   = "EC_{50}",
    hill   = "Hill",
    Hill   = "Hill",
    CLPara = "CL_{Para}",
    kin    = "k_{in}",
    ke     = "k_{e}",

    # Error Parameters:
    error_ADD1  = "error_{ADD1}",
    error_PROP1 = "error_{PROP1}",
    error_ADD2  = "error_{ADD2}",
    error_PROP2 = "error_{PROP2}",
    error_ADD3  = "error_{ADD3}",
    error_PROP3 = "error_{PROP3}",
    error_ADD4  = "error_{ADD4}",
    error_PROP4 = "error_{PROP4}"
  )

  return(strToLaTeX)
}


#' Convert equations/parameters to LaTeX mode
#'
#' Parameters and/or equations within IQRtable can be automatically converted to a LaTeX equation
#' for nicer reporting. `colToChange` specify to which column the changes should be made, while
#' `ListStrToLatex` is a list providing the characters to detect and their respective LaTeX notation.
#'
#' @description
#' @param tableInputPath Path of the IQR table to adjust
#' @param tableOutputPath Output path to save the updated table (Default: `NULL` which will save it to tableInputPath"_LaTeX").
#' @param colToChange Column in which to apply the changes (Default: `PARAMETER`)
#' @param ListStrToLatex list the character to change and their LaTeX expression; e.g. \code{list(Fabs0 = "F_{abs,0}", ka = "k_{a}", CL = "CL")} (Default: `MMVmalaria::parameter_strToLaTeX()`)
#' @param TitleRow To change adjust the title of each column (Default: `FALSE`).
#' @param keepCol List of column of th etable to keep (Default: `NULL` which will keep all columns).
#' @param report Logical to choose report mode or not (Default: `TRUE`).
#' @param FLAGcompliance Logical to choose compliance mode or not (Default: `TRUE`).
#'
#' @return An IQR table with the specified changes
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Reporting
table_strToLaTeX <- function(tableInputPath,
                             tableOutputPath = NULL,
                             colToChange     = "PARAMETER",
                             ListStrToLatex  = parameter_strToLaTeX(),
                             TitleRow        = FALSE,
                             keepCol         = NULL,
                             report          = TRUE,
                             FLAGcompliance  = TRUE) {

  # Define tableOutputPath if NULL:
  if (is.null(tableOutputPath)){tableOutputPath = paste0(gsub(".txt", "", tableInputPath, fixed=TRUE),"_LaTeX.txt")}

  # Import table:
  TableImport <- MMVbase::IQRtableToDataFrame(tableInputPath)
  dataFrame   <- TableImport$dataFrame
  xTitle      <- TableImport$xTitle
  xFooter     <- TableImport$xFooter

  # Check that 'colToChange' is a column of dataFrame:
  if (!(colToChange %in% names(dataFrame)) & !TitleRow)
  {
    stop("The value of the variable 'colToChange needs to be a column of the table to change.")
  }

  # Change along a column:
  if (colToChange %in% names(dataFrame)){
    # As Character:
    dataFrame[,colToChange] <- as.character(dataFrame[,colToChange])

    for (k in 1:nrow(dataFrame)){

      # Check if we need to add $: By Default FALSE
      AddDollars = FALSE

      # Check if there is omega, and change to LaTeX if present:
      if (grepl("omega", dataFrame[k,colToChange])){
        dataFrame[k,colToChange] <- gsub("omega(", "\\omega_{", dataFrame[k,colToChange], fixed=TRUE)
        dataFrame[k,colToChange] <- gsub(")", "}", dataFrame[k,colToChange]             , fixed=TRUE)

        # $s will need to be added:
        AddDollars = TRUE
      }

      # Check if there is beta, and change to LaTeX if present:
      if (grepl("beta", dataFrame[k,colToChange])){
        dataFrame[k,colToChange] <- gsub("beta_", "\\beta_{", dataFrame[k,colToChange], fixed=TRUE)
        dataFrame[k,colToChange] <- gsub("(", ",", dataFrame[k,colToChange]           , fixed=TRUE)
        dataFrame[k,colToChange] <- gsub(")", "}", dataFrame[k,colToChange]           , fixed=TRUE)

        # $s will need to be added:
        AddDollars = TRUE
      }

      # Check if there is rho, and change to LaTeX if present:
      if (grepl("corr_(", dataFrame[k,colToChange], fixed=TRUE)){
        dataFrame[k,colToChange] <- gsub("corr_(", "\\rho_{", dataFrame[k,colToChange], fixed=TRUE)
        dataFrame[k,colToChange] <- gsub(")", "}", dataFrame[k,colToChange]          , fixed=TRUE)

        # $s will need to be added:
        AddDollars = TRUE
      }

      # Change parameter to LaTeX:
      for (parLaTeX in names(ListStrToLatex)){
        if(grepl(parLaTeX, dataFrame[k,colToChange])){
          # Change:
          idx                      <- which(names(ListStrToLatex)==parLaTeX)
          dataFrame[k,colToChange] <- gsub(names(ListStrToLatex)[idx], ListStrToLatex[[parLaTeX]], dataFrame[k,colToChange], fixed=TRUE)

          # $s will need to be added:
          AddDollars = TRUE
        }
      }

      # Add $s:
      if (AddDollars){
        dataFrame[k,colToChange] <- paste0("$", dataFrame[k,colToChange], "$")
      }
    }
  }


  # Change Name of DataFrame:
  if (TitleRow){
    for (k in 1:length(names(dataFrame))){
      # Check if we need to add $: By Default FALSE
      AddDollars = FALSE

      # Check if there is omega, and change to LaTeX if present:
      if (grepl("omega", names(dataFrame)[k])){
        names(dataFrame)[k] <- gsub("omega(", "\\omega_{", names(dataFrame)[k], fixed=TRUE)
        names(dataFrame)[k] <- gsub(")", "}", names(dataFrame)[k]             , fixed=TRUE)

        # $s will need to be added:
        AddDollars = TRUE
      }

      # Check if there is beta, and change to LaTeX if present:
      if (grepl("beta", names(dataFrame)[k])){
        names(dataFrame)[k] <- gsub("beta_", "\\beta_{", names(dataFrame)[k], fixed=TRUE)
        names(dataFrame)[k] <- gsub("(", ",", names(dataFrame)[k]           , fixed=TRUE)
        names(dataFrame)[k] <- gsub(")", "}", names(dataFrame)[k]           , fixed=TRUE)

        # $s will need to be added:
        AddDollars <- TRUE
      }

      # Check if there is rho, and change to LaTeX if present:
      if (grepl("corr_(", dataFrame[k,colToChange], fixed=TRUE)){
        dataFrame[k,colToChange] <- gsub("corr_(", "\\rho_{", dataFrame[k,colToChange], fixed=TRUE)
        dataFrame[k,colToChange] <- gsub(")", "}", dataFrame[k,colToChange]          , fixed=TRUE)

        # $s will need to be added:
        AddDollars = TRUE
      }

      # Change parameter to LaTeX:
      for (parLaTeX in names(ListStrToLatex)){
        if(grepl(parLaTeX, names(dataFrame)[k])){
          # Change:
          idx                 <- which(names(ListStrToLatex)==parLaTeX)
          names(dataFrame)[k] <- gsub(names(ListStrToLatex)[idx], ListStrToLatex[[parLaTeX]], names(dataFrame)[k], fixed=TRUE)

          # $s will need to be added:
          AddDollars = TRUE
        }
      }

      # Add $s:
      if (AddDollars){
        names(dataFrame)[k] <- paste0("$", names(dataFrame)[k], "$")
      }
    }
  }

  # Keep columns of interest:
  if (!is.null(keepCol) & !is.na(keepCol)){
    dataFrame <- dataFrame[,keepCol]
  }

  # Generat output table:
  IQRoutputTable(dataFrame, xfooter=xFooter, xtitle=xTitle,
                 filename=tableOutputPath, report=report)

  # Remove log file if Compliance not needed:
  #   e.g. for publication
  if (!FLAGcompliance && file.exists(paste0(tableOutputPath, ".log"))){
    file.remove(paste0(tableOutputPath, ".log"))
  }
}


#' Generate Estimate Table for Report
#'
#' Generate a table of estimate from IQR or GPF output in a summarized format
#'
#' @param projectPath IQR pr GPF project path
#' @param tableOutputPath
#' @param FLAGLaTeX Flag to convert parameters name into LaTeX format (Default: `TRUE`)
#' @param setting Optional settings for `table_strToLaTeX`; `ListStrToLatex`, `TitleRow`, `report` and `FLAGcompliance`
#'
#' @return List containing the table as data.frame, title and footer
#'
#' @examples
#'
#' @export
#' @importFrom MMVbase get_fileNameExtension
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Reporting
generate_EstimateTableForReport <- function(projectPath,
                                            tableOutputPath = NULL,
                                            FLAGLaTeX       = TRUE,
                                            setting         = NULL)
{
  #---------------------------------------------------#
  # STEP 0: Define default setting ----
  #---------------------------------------------------#

  # Define tableOutputPath:
  if(is.null(tableOutputPath)){
    if(is_IQRnlmeProject(projectPath = projectPath) || is_IQRsysProject(projectPath = projectPath)){
      tableOutputPath <- file.path(projectPath, "project_parameters_table_Report.txt")
    }else if(MMVbase::get_fileNameExtension(projectPath) %in% c("xlsx","xls","csv")){
      tableOutputPath <- paste0(gsub(paste0(".",MMVbase::get_fileNameExtension(projectPath)),"",projectPath),"_Report.txt")
    }else{
      stop("Please define 'tableOutputPath'")
    }
  }

  # Define variable for LaTeX output:
  #if ("tableOutputPath" %in% names(setting)){tableOutputPath = setting$tableOutputPath} else{tableOutputPath = NULL                  }
  if ("ListStrToLatex" %in% names(setting)){ListStrToLatex = setting$ListStrToLatex} else{ListStrToLatex = parameter_strToLaTeX()}
  if ("TitleRow"       %in% names(setting)){TitleRow       = setting$TitleRow      } else{TitleRow       = FALSE                 }
  if ("report"         %in% names(setting)){report         = setting$report        } else{report         = TRUE                  }
  if ("FLAGcompliance" %in% names(setting)){FLAGcompliance = setting$FLAGcompliance} else{FLAGcompliance = TRUE                  }


  #---------------------------------------------------#
  # STEP 1: Prepare Table ----
  #---------------------------------------------------#

  # Load into parameters into GPF (ex-XLS format):
  if(is_IQRnlmeProject(projectPath = projectPath) || is_IQRsysProject(projectPath = projectPath) || MMVbase::get_fileNameExtension(projectPath) %in% c("xlsx","xls","csv")){
    estimates <- load_GPF(projectPath)$estimates

  }else if(is_GPF(projectPath)){
    estimates <- projectPath$estimates

  }else{
    stop("'projectPath' is not a valid argument. It should be a IQRnlmeProject, IQRsysProject or a GPF file: Please adjust the argument.")
  }

  # Split 'estimates' by category:
  #   Pop Parameters
  estimates.Pop <- subset(estimates,TYPE=="MODEL PARAMETER")
  #   Continuous Covariates
  estimates.ConCov <- subset(estimates,TYPE=="CONTINUOUS COVARIATE")
  #   Categorical Covariates
  estimates.CatCov <- subset(estimates,TYPE=="CATEGORICAL COVARIATE")
  #   Correlation
  estimates.Cor <- subset(estimates,TYPE=="IIV CORRELATION")
  #   Error
  estimates.Error <- subset(estimates,grepl("ERROR",TYPE))


  # Detect PARA and COV:
  #   Continuous Covariates
  if(nrow(estimates.ConCov)>0){
    estimates.ConCov$PARA <- sapply(estimates.ConCov$PARAMETER, function(x){
      regmatches(x, regexec("beta_(.*?)\\(",x))[[1]][2]
    },
    USE.NAMES = FALSE)
    estimates.ConCov$COV <- sapply(estimates.ConCov$PARAMETER, function(x){
      regmatches(x, regexec("\\((.*?)\\)",x))[[1]][2]
    },
    USE.NAMES = FALSE)
  }
  #   Categorical Covariates
  if(nrow(estimates.CatCov)>0){
    estimates.CatCov$PARA <- sapply(estimates.CatCov$PARAMETER, function(x){
      regmatches(x, regexec("beta_(.*?)\\(",x))[[1]][2]
    },
    USE.NAMES = FALSE)
    estimates.CatCov$COV <- sapply(estimates.CatCov$PARAMETER, function(x){
      regmatches(x, regexec("\\((.*?)\\)",x))[[1]][2]
    },
    USE.NAMES = FALSE)
    estimates.CatCov$COV.VALUE <- as.numeric(strsplit(estimates.CatCov$COV,"_")[[1]][2])
    estimates.CatCov$COV       <- as.character(strsplit(estimates.CatCov$COV,"_")[[1]][1])
  }


  # Construct Output for Population Parameters:
  estimates.Pop$VALUETXT <- paste0(signif(estimates.Pop$VALUE,3),
                                   ifelse(estimates.Pop$VALUE.RSE.PERCENT==0,
                                          "",
                                          paste0("(",round(estimates.Pop$VALUE.RSE.PERCENT,1),")")))
  estimates.Pop$IIVTXT <- paste0(round(estimates.Pop$IIV,2), "(", round(estimates.Pop$IIV.RSE.PERCENT,1), ")")
  #   Add Continuous Covariate
  if(nrow(estimates.ConCov)>0){
    estimates.ConCov$VALUETXT <- gsub("X=",
                                      "",
                                      estimates.ConCov$COV.FORMULA)
    for(k in 1:nrow(estimates.ConCov)){
      idx_k <- which(estimates.Pop$PARAMETER==estimates.ConCov$PARA[k])
      estimates.ConCov$VALUETXT[k] <- gsub("Beta",
                                           paste0("(",
                                                  estimates.ConCov$VALUE[k],
                                                  ifelse(estimates.ConCov$VALUE.RSE.PERCENT[k]==0,
                                                         "",
                                                         paste0("(",round(estimates.ConCov$VALUE.RSE.PERCENT[k],1),")")),
                                                  ")"),
                                           estimates.ConCov$VALUETXT[k])
      estimates.ConCov$VALUETXT[k] <- gsub("REF",
                                           estimates.ConCov$COV.REFERENCE[k],
                                           estimates.ConCov$VALUETXT[k])
      estimates.ConCov$VALUETXT[k] <- gsub("X_ref",
                                           estimates.Pop$VALUETXT[idx_k],
                                           estimates.ConCov$VALUETXT[k])
      estimates.Pop$VALUETXT[idx_k] <- estimates.ConCov$VALUETXT[k]
    }
  }
  #   Add Categorical Covariate
  if(nrow(estimates.CatCov)>0){
    estimates.CatCov$VALUETXT <- gsub("X=",
                                      "",
                                      estimates.CatCov$COV.FORMULA)
    for(k in 1:nrow(estimates.CatCov)){
      idx_k <- which(estimates.Pop$PARAMETER==estimates.CatCov$PARA[k])
      estimates.CatCov$VALUETXT[k] <- gsub("Beta",
                                           paste0("(",
                                                  signif(estimates.CatCov$VALUE[k],2),
                                                  ifelse(estimates.CatCov$VALUE.RSE.PERCENT[k]==0,
                                                         "",
                                                         paste0("(",round(estimates.CatCov$VALUE.RSE.PERCENT[k],1),")")),
                                                  "[If ",
                                                  estimates.CatCov$COV[k],
                                                  "=",
                                                  estimates.CatCov$COV.VALUE[k],
                                                  "])"),
                                           estimates.CatCov$VALUETXT[k])
      estimates.CatCov$VALUETXT[k] <- gsub("REF",
                                           estimates.CatCov$COV.REFERENCE[k],
                                           estimates.CatCov$VALUETXT[k])
      estimates.CatCov$VALUETXT[k] <- gsub("X_ref",
                                           estimates.Pop$VALUETXT[idx_k],
                                           estimates.CatCov$VALUETXT[k])
      estimates.Pop$VALUETXT[idx_k] <- estimates.CatCov$VALUETXT[k]
    }
  }


  # Add Unit to NAME:
  estimates.Pop$NAME <- ifelse(estimates.Pop$UNIT==estimates.Pop$NAME,
                               paste0(estimates.Pop$NAME,
                                      " [-]"),
                               paste0(estimates.Pop$NAME,
                                      " [",
                                      estimates.Pop$UNIT,
                                      "]"))


  # Prepare Output:
  tableOut <- estimates.Pop[,c("PARAMETER","VALUETXT","IIVTXT","NAME")]
  names(tableOut) <- c("Parameter", "Estimate", "IIV", "Description")


  # Convert into LaTeX format:
  if (FLAGLaTeX){
    # Save as text:
    IQRoutputTable(tableOut,
                   filename = "temp.txt")

    # Convert to LaTeX:
    table_strToLaTeX(tableInputPath  = "temp.txt",
                     tableOutputPath = NULL,
                     colToChange     = "Parameter",
                     ListStrToLatex  = ListStrToLatex,
                     TitleRow        = TitleRow,
                     keepCol         = names(tableOut),
                     report          = report,
                     FLAGcompliance  = FLAGcompliance)

    # Load LaTeX table:
    tableOut <- as.data.frame(as.matrix(MMVbase::IQRtableToDataFrame("temp_LaTeX.txt")$dataFrame),stringsAsFactors = FALSE)

    # Remove temporary tables:
    if(file.exists("temp.txt")){
      file.remove("temp.txt")
    }
    if(file.exists("temp.txt.log")){
      file.remove("temp.txt.log")
    }
    if(file.exists("temp_LaTeX.txt")){
      file.remove("temp_LaTeX.txt")
    }
    if(file.exists("temp_LaTeX.txt.log")){
      file.remove("temp_LaTeX.txt.log")
    }
  }

  # Save as text:
  xtitle  <- "Parameter estimate table."
  xfooter <- "Estimate (Relative Standard Error %); IIV = Inter-Individual Variability"
  IQRoutputTable(tableOut,
                 xtitle  = xtitle,
                 xfooter = xfooter,
                 filename = tableOutputPath)

  # Output:
  out <- list(dataFrame = tableOut,
              xtitle    = xtitle,
              xfooter   = xfooter)
  out
}


#' table_covariateSummaryResults
#'
#' @description
#' @param tableInputPath
#' @param tableOutputPath (Default: `NULL`).
#' @param report Default: TRUE
#' @param FUNnameModel Default: function(x) {
#'    gsub("MODEL_", "", x)
#'}
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Reporting
table_covariateSummaryResults <- function(tableInputPath,
                                          tableOutputPath = NULL,
                                          report          = TRUE,
                                          FUNnameModel = function(x){gsub("MODEL_","",x)}) {

  # Define tableOutputPath if NULL:
  if (is.null(tableOutputPath)){tableOutputPath = paste0(gsub(".txt", "", tableInputPath),"_Summary.txt")}

  # Import table:
  TableImport <- MMVbase::IQRtableToDataFrame(tableInputPath)
  xTitle      <- TableImport$xTitle
  xFooter     <- TableImport$xFooter
  TableImport <- TableImport$dataFrame

  # Generate dataframe to be exported:
  TableExport <- data.frame(MODEL       = character(0),
                            Metric      = character(0),
                            Covariate   = character(0),
                            Value       = character(0),
                            CI95        = character(0),
                            STAT.SIGNIF = character(0),
                            stringsAsFactors = FALSE)
  NewModel = TRUE
  for (k in 1:nrow(TableImport)){

    if (NewModel){
      # Get Metric:
      Metric.Type  <- as.character(strsplit(TableImport$MODEL[k+1],": ")[[1]][1])
      Metric.Value <- as.character(strsplit(TableImport$MODEL[k+1],": ")[[1]][2])

      # New Row:
      NewRow <- data.frame(MODEL       = TableImport$MODEL[k],
                           Metric      = Metric.Value,
                           Covariate   = TableImport$NAME[k],
                           Value       = as.character(TableImport$VALUE[k]),
                           CI95        = paste0("[",as.character(TableImport$LOW.95.CI[k]),";",as.character(TableImport$HIGH.95.CI[k]),"]"),
                           STAT.SIGNIF = TableImport$STAT.SIGNIF.05[k],
                           stringsAsFactors = FALSE)

      # Add New Row
      TableExport <- rbind(TableExport,NewRow)

      # Keep Reading the result of this model:
      NewModel = FALSE

    } else if(grepl("^\\s*$", TableImport$NAME[k])){
      # Do Nothing

    } else if(grepl("NAME", TableImport$NAME[k])){

      # Next Row should be a new model:
      NewModel = TRUE

    }else{

      # New Row:
      NewRow <- data.frame(MODEL       = "",
                           Metric      = "",
                           Covariate   = TableImport$NAME[k],
                           Value       = as.character(TableImport$VALUE[k]),
                           CI95        = paste0("[",as.character(TableImport$LOW.95.CI[k]),";",as.character(TableImport$HIGH.95.CI[k]),"]"),
                           STAT.SIGNIF = TableImport$STAT.SIGNIF.05[k],
                           stringsAsFactors = FALSE)
      TableExport <- rbind(TableExport,NewRow)
    }
  }

  # Rename MODELs:
  if (!is.null(FUNnameModel)){
    TableExport$MODEL <- FUNnameModel(TableExport$MODEL)
  }

  # Rename columns:
  names(TableExport) <- c("MODEL", Metric.Type, "Covariate", "Value", "95% C.I.", "Signif.")

  # Generat output table:
  IQRoutputTable(TableExport, xfooter=xFooter, xtitle=xTitle,
                 filename=tableOutputPath, report=report)

}
#' table_PKsummaryCompare
#'
#' @description
#' @param tableInputPath1
#' @param tableInputPath2
#' @param tableOutputPath (Default: `NULL`).
#' @param renameMODEL (Default: `NULL`).
#' @param addDescription (Default: `NULL`).
#' @param Metric Default: 'BIC'
#' @param xTitle (Default: `NULL`).
#' @param xFooter (Default: `NULL`).
#' @param report Default: TRUE
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Reporting
table_PKsummaryCompare <- function(tableInputPath1,
                                   tableInputPath2,
                                   tableOutputPath = NULL,
                                   renameMODEL     = NULL,
                                   addDescription  = NULL,
                                   Metric          = "BIC",
                                   xTitle          = NULL,
                                   xFooter         = NULL,
                                   report          = TRUE) {

  # Define tableOutputPath if NULL:
  if(is.null(tableOutputPath)){
    tableOutputPath1 <- paste0(gsub("_metrics", "", gsub("txt", "", gsub(".txt", "", tableInputPath1),fixed = TRUE)),"_PKsummaryCompare.txt")
    tableOutputPath2 <- paste0(gsub("_metrics", "", gsub("txt", "", gsub(".txt", "", tableInputPath2),fixed = TRUE)),"_PKsummaryCompare.txt")
  } else{
    tableOutputPath1 <- tableOutputPath
    tableOutputPath2 <- tableOutputPath
  }

  # Import tables:
  #   Table 1
  TableImport1 <- MMVbase::IQRtableToDataFrame(tableInputPath1)
  TableImport1 <- TableImport1$dataFrame
  #   Table 2
  TableImport2 <- MMVbase::IQRtableToDataFrame(tableInputPath2)
  TableImport2 <- TableImport2$dataFrame
  #   Merge Table:
  TableImport <- rbind(TableImport1,TableImport2)

  # Get ID Code of the PK model:
  if (!is.null(renameMODEL)){
    TableImport$MODEL <- renameMODEL(TableImport)
  }

  # Table to Export:
  TableExport <- data.frame(MODEL  = TableImport$MODEL,
                            Metric = as.numeric(as.character(TableImport[,Metric])),
                            stringsAsFactors = FALSE)

  # Add Description:
  if (!is.null(addDescription)){
    TableExport$Description <- addDescription(TableExport)
  }

  # Rename columns:
  if (!is.null(addDescription)){
    names(TableExport) <- c("MODEL", Metric, "Description")
  }else{
    names(TableExport) <- c("MODEL", Metric)
  }

  # Order by Metric:
  TableExport <- TableExport[order(TableExport[,Metric]),]


  # Define Title:
  if (is.null(xTitle))
    xTitle = paste0("Comparison of PK models ordered by ", Metric)

  # Define Footer:
  if (is.null(xFooter))
    xFooter = paste0("Models ordered by ", Metric, "<br>", Metric, " values are rounded", "<br> Models in ", tableInputPath1," and in ", tableInputPath2," are compared")

  # Generat output table:
  #   It puts it in two locations if tableOutput was not specified
  IQRoutputTable(TableExport, xfooter=xFooter, xtitle=xTitle,
                 filename=tableOutputPath1, report=report)
  IQRoutputTable(TableExport, xfooter=xFooter, xtitle=xTitle,
                 filename=tableOutputPath2, report=report)
}
#' table_PKsummaryResults
#'
#' @description
#' @param tableInputPath
#' @param tableOutputPath (Default: `NULL`).
#' @param colToAdd Default: c("Nr. of Compartments", "Clearance", "Absorption", "Time Lag",
#'    "Error Model")
#' @param Metric Default: 'BIC'
#' @param xTitle (Default: `NULL`).
#' @param xFooter (Default: `NULL`).
#' @param report Default: TRUE
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Reporting
table_PKsummaryResults <- function(tableInputPath,
                                   tableOutputPath = NULL,
                                   colToAdd        = c("Nr. of Compartments", "Clearance", "Absorption", "Time Lag", "Error Model"),
                                   Metric          = "BIC",
                                   xTitle          = NULL,
                                   xFooter         = NULL,
                                   report          = TRUE) {

  # Define tableOutputPath if NULL:
  if (is.null(tableOutputPath)){tableOutputPath = paste0(gsub("_metrics", "", gsub("txt", "", gsub(".txt", "", tableInputPath),fixed = TRUE)),"_PKsummaryResults.txt")}

  # Import table:
  TableImport <- MMVbase::IQRtableToDataFrame(tableInputPath)
  TableImport <- TableImport$dataFrame

  # Get ID Code of the PK model:
  TableImport$IDcode <- as.character(t(as.data.frame(strsplit(TableImport$MODEL,"_")))[,2])

  # Table to Export:
  TableExport <- data.frame(MODEL  = TableImport$MODEL,
                            Metric = as.numeric(as.character(TableImport[,Metric])),
                            Ncpt   = substr(TableImport$IDcode,1,1),
                            CL     = ifelse(substr(TableImport$IDcode,2,2)=="1","Linear",ifelse(substr(TableImport$IDcode,2,2)=="2","Saturable","Linear+Saturable")),
                            Abs    = ifelse(substr(TableImport$IDcode,3,3)=="0","0th order","1st order"),
                            Lag    = ifelse(substr(TableImport$IDcode,4,4)=="0","No Lag","With Lag"),
                            Err    = ifelse(substr(TableImport$IDcode,5,5)=="0","Additive",ifelse(substr(TableImport$IDcode,5,5)=="1","Proportional","Combined")),
                            stringsAsFactors = FALSE)

  # Rename columns:
  names(TableExport) <- c("MODEL", Metric, "Nr. of Compartments", "Clearance", "Absorption", "Time Lag", "Error Model")

  # Order by Metric:
  TableExport <- TableExport[order(TableExport[,Metric]),]

  # Export only the column of interest:
  TableExport <- TableExport[,c("MODEL", Metric, colToAdd)]

  # Define Title:
  if (is.null(xTitle))
    xTitle = paste0("Results of the tested structural PK models ordered by ", Metric)

  # Define Footer:
  if (is.null(xFooter))
    xFooter = paste0("Models ordered by ", Metric,"<br>", Metric, " values are rounded")

  # Generat output table:
  IQRoutputTable(TableExport, xfooter=xFooter, xtitle=xTitle,
                 filename=tableOutputPath, report=report)

}

#' @description Convert list input to text
#' @param list list to convert
#' @param sepWord Character string to separate elements of list

#' @return
#' @export
#' @author Unknown
#' @family Reporting
reportText_convertListToText <- function(list, sepWord){
  if (length(list) == 0){
    warning("Empty list cannot be converted to text")
  } else if (length(list) == 1) {
    return(list)
  } else if (length(list) == 2) {
    text <- paste(list[1], sepWord, list[2])
    return(text)
  } else {
    text <- paste(sepWord, list[length(list)])
    text2 <- paste(list[1:(length(list)-1)], collapse = ", ")
    return(paste(text2, text))
  }
}

#' @description Generate report text from list of PK models to test
#' @param PKmodelTest list containing description of PK models to test

#' @return
#' @export
#' @author Unknown
#' @family Reporting
reportText_PKmodelHypos <- function(PKmodelTest){
  # General text
  text <- "PK models with"

  # Compartments
  text <- paste(text,
                paste(reportText_convertListToText(PKmodelTest$Compartments, "or"),
                      "compartments"))

  # Absorption
  text <- append(text,
                 paste(reportText_convertListToText(PKmodelTest$Absorption, "or"),
                       "absorption"))

  # Elimination
  reptext <- lapply(X = PKmodelTest$Elimination,
                    FUN = function(t) gsub(pattern = "+", replacement = " and ", x = t, fixed = TRUE))
  text <- append(text,
                 paste(reportText_convertListToText(reptext, "or"),
                       "elimination"))

  # LagTime
  reptext <- ifelse(PKmodelTest$LagTime == FALSE, "no", "a")

  text <- append(text,
                 paste(reptext,
                       "lag time"))

  # ErrorModels
  reptext <- lapply(X = PKmodelTest$ErrorModels,
                    FUN = function(t) gsub(pattern = "rel", replacement = "relative", x = t, fixed = TRUE))
  reptext <- lapply(X = reptext,
                    FUN = function(t) gsub(pattern = "abs", replacement = "absolute", x = t, fixed = TRUE))
  reptext <- lapply(X = reptext,
                    FUN = function(t) gsub(pattern = "absoluterelative", replacement = "absolute and relative", x = t, fixed = TRUE))

  text <- append(text,
                 paste(reportText_convertListToText(reptext, "or"),
                       "error models"))

  # General text
  text[length(text)] <- paste(text[length(text)], "were tested")

  # Combine to single string
  return(reportText_convertListToText(text, "and"))
}
