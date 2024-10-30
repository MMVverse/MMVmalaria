#' summary_ParametersDDI
#' Summarize selected parameters over multiple SysFit projects
#'
#' Reads the parameters from SysFit projects and creates an
#' overview table with the different selected parameters (typically
#' the DDI parameters) and fit metrics (e.g. BIC). One column
#' per parameter.
#'
#' @param projfolder character vector with the project folders
#' @param pars Parameter names or metric names (OBJ, AIC, BIC) to be
#' exatracted from the parameter table.
#'
#' @return IQRtable object
#'
#' @seealso \code{\link{summary_ParametersMulti}}
#'
#' @export
#' @author Daniel Kaschek (IntiQuan), Mohammed H. Cherkaoui (MMV)
summary_ParametersDDI <- function(projfolders,
                                  pars = c("BIC", "Alpha", "Beta", "Gamma")) {

  # Transform to list:
  multi <- as.list(seq_along(projfolders))
  names(multi) <- projfolders

  # Get model names from folder names:
  modelnames__ <- sapply(strsplit(names(multi), "/", fixed = TRUE), function(x__) rev(x__)[1])

  # Produce list of tables (models times purpose):
  tables__ <- lapply(seq_along(multi), function(i__) {

    # Read single table from project:
    res__ <- load_IQRsysProject(names(multi)[i__], FLAGresultsOnly = TRUE)
    M__ <- IQRtools::tablePars_IQRsysModel(res__)

    # Merge value and RSE:
    M__[,3] <- ifelse(M__[, 3] == "", "", paste0("(", M__[,3], ")"))
    M__[,2] <- paste(M__[,2], M__[, 3])

    # Reduce to selected parameters and transpose the table:
    #M__ <- t(M__[c(1, unlist(lapply(pars, function(x__) grep(x__, M__[,1])))),1:2])
    M__ <- t(M__[(M__[,1] %in% pars),1:2])

    # Reorder column:
    M__ <- M__[,match(pars, M__[1,],nomatch = FALSE)]

    # Fill first column with model name information:
    #M__[,1] <- c("Model", modelnames__[i__])
    M__ <- cbind(c("Model", modelnames__[i__]),M__)


    # Output:
    return(M__)
  })

  # Reduce table:
  outtable__ <- Reduce(rjoin_IQRtable, tables__)

  # Output:
  return(outtable__)
}
#' summary_ParametersMulti
#' Summarize parameters over multiple SysFit projects
#'
#' Reads the parameter results from different project folders and
#' creates an output table with separate columns for each model/folder.
#'
#' @param projfolder Character vector with the project folders
#' @param sortBy Either NULL (no sorting), character denoting a model metric (BIC, OBJ, AIC)
#'               or a parameter by which to sort the output columns.
#'
#' @return IQRtable object
#'
#' @seealso \code{\link{summary_ParametersDDI}}
#'
#' @export
#' @author Daniel Kaschek (IntiQuan), Mohammed H. Cherkaoui (MMV)
summary_ParametersMulti <- function(projfolders,
                                    sortBy = NULL) {

  # Get all models in 'projfolders':
  multi        <- as.list(seq_along(projfolders))
  names(multi) <- projfolders

  # Get model names from folder names:
  modelnames__ <- sapply(strsplit(names(multi), "/", fixed = TRUE), function(x__) rev(x__)[1])

  # Function to be applied to a list of lists which swaps inner and outer index:
  swaplist <- function(mylist) {
    nouter__ <- length(mylist)
    ninner__ <- max(sapply(mylist, length))

    # Get names from first row with all non-NA entries:
    names_inner__ <- do.call(rbind, lapply(mylist, function(x__) names(x__)[1:ninner__]))
    nonNA__       <- apply(names_inner__, 1, function(line__) !any(is.na(line__)))
    names_inner__ <- names_inner__[which(nonNA__)[1],]

    # Construct output list:
    out_outer__ <- lapply(names_inner__, function(i__) {
      out_inner__ <- lapply(mylist, function(l__) l__[[i__]])
      names(out_inner__) <- names(mylist)

      return(out_inner__)
    })

    names(out_outer__) <- names_inner__

    return(out_outer__)
  }

  # Produce list of list of tables (models times purpose):
  tables__ <- lapply(seq_along(multi), function(i__) {

    # Read single table from project:
    res__   <- load_IQRsysProject(names(multi)[i__], FLAGresultsOnly = TRUE)
    M__     <- IQRtools::tablePars_IQRsysModel(res__)
    M__[,3] <- ifelse(M__[, 3]=="", "", paste0("(", M__[,3], ")"))
    M__[,2] <- paste(M__[,2], M__[, 3])

    splitidx__ <- grep("^\\*\\*", M__[,1])

    # Split the table into the parameter purposes (fixed, random, metric, etc.):
    out__ <- lapply(seq_along(splitidx__), function(j__) {
      start__ <- splitidx__[j__] + 1

      if (j__==length(splitidx__))
        end__ <- length(M__[,1])

      else
        end__ <- splitidx__[j__ + 1] - 1

      M__[start__:end__, 1:2]
    })

    names(out__) <- M__[splitidx__, 1]

    return(out__)
  })

  tables__ <- swaplist(tables__)

  # Join tables first for each purpose then compose tables of different purpose:
  tables__ <- lapply(seq_along(tables__), function(i__) {
    t__ <- Reduce(cjoin_IQRtable, tables__[[i__]])
    t__ <- t__[apply(t__, 1, function(line__) any(gsub(" ", "", line__) != "")),]
    compose_IQRtable(H = names(tables__)[i__], M = t__, pattern = "H\n M")
  })
  tables__ <- composeByRule_IQRtable(tables__, rule = "x \n \n y")

  # Provide header line to output table
  header__ <- matrix(c("Parameter", modelnames__), nrow = 1)
  outtable__ <- compose_IQRtable(H = header__, M = tables__, pattern = "H\nM")

  # Sort by:
  if (!is.null(sortBy)) {
    values__          <- outtable__[match(sortBy, trimws(outtable__[, 1])), -1]
    values__          <- utils::type.convert(values__)
    outtable__        <- outtable__[, c(1, (2:ncol(outtable__))[order(values__)])]
    class(outtable__) <- union("IQRtable", class(outtable__))
  }

  # Output:
  return(outtable__)
}

#' table_ModelEstimates
#'
#' @description
#' @param modelFolderList
#' @param filename Default: NULL
#' @param xtitle Default: NULL
#' @param xfooter Default: NULL
#' @param FLAGreturnOBJECT Default: FALSE
#' @param FLAGreport Default: TRUE
#' @return
#' @export
#' @author Anne Kümmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family ModelComparison
#' @importFrom plyr llply join_all
table_ModelEstimates <- function(modelFolderList,
                                 filename = NULL,
                                 xtitle   = NULL,
                                 xfooter  = NULL,
                                 FLAGreturnOBJECT = FALSE,
                                 FLAGreport = TRUE) {

  # Get estimates for all models in the list:
  estimateTables <- plyr::llply(modelFolderList, getEstimateTable)

  # Rename estimate column to differentiate models:
  if (is.null(names(estimateTables))) {
    warning(paste0("No model names given. Just generic names Model1, Model2, ..."))
    names(estimateTables) <- paste0("Model", 1:length(estimateTables))
  }
  for (namek in names(estimateTables)) {
    names(estimateTables[[namek]]) <- gsub("Estimate", namek, names(estimateTables[[namek]]) )
  }

  # Merge to one table:
  estimateTableAll <- plyr::join_all(estimateTables, type = "full")

  if (!is.null(filename))
    dummy <- IQRoutputTable(estimateTableAll, filename = filename, report = FLAGreport, xtitle = xtitle, xfooter = xfooter)

  if (FLAGreturnOBJECT) return(estimateTableAll)

}

#' getEstimateTable
#'
#' @description
#' @param modelFolder
#' @return
#' @export
#' @author Anne Kümmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family ModelComparison
getEstimateTable <- function(modelFolder) {
  # Detect IQR version:
  IQRversion <- sessionInfo()$otherPkgs$IQRtools$Version

  # Check whether NLME or SysFit results:
  if (exists("is_IQRsysProject", mode="function") && is_IQRsysProject(modelFolder)){
    modelType <- "sysFIT"
  }else if(is_IQRnlmeProject(modelFolder)){
    modelType <- "NLME"
  }else if(is_IQRnlmeProject(modelFolder)){
    folderContent <- list.files(modelFolder)
    if (any(grepl(".RData$",folderContent)) & all(!grepl("^RESULTS$",folderContent, ignore.case = FALSE))){
      modelType <- "sysFITold"
    }
  }else{
    stop(paste0("Model type (naive pool or NLME) not detected. Please check input model folder:\n", modelFolder,"\n"))
  }
  # modelFolder
  #
  #
  # folderContent <- list.files(modelFolder)
  # modelType     <- NA
  # if (any(grepl(".RData$",folderContent)) & all(!grepl("^RESULTS$",folderContent, ignore.case = FALSE))) modelType <- "sysfit"
  # if (all(!grepl(".RData$",folderContent)) & any(grepl("^RESULTS$",folderContent, ignore.case = FALSE))) modelType <- "NLME"
  # if (is.na(modelType)) stop(paste0("Model type (naive pool or NLME) not detected. Please check input model folder:\n", modelFolder,"\n"))

  # Read results
  if (modelType=="sysFITold") {
    # Load results:
    resultsFile <- grep(".RData$",folderContent, value = TRUE)
    load(file.path(modelFolder, resultsFile))
    resultsSysfit <- out$BestFit

    # Get estimates and their confidence interval:
    #   First, check whether confidence intervals for parameter estimates are available in model output.
    #   If not, an old version of IQRtools was used and confidence intervals will not be displayed.
    flagCI <- "xopt_95CI_high" %in% names(resultsSysfit)
    if (!flagCI) resultsSysfit$xopt_95CI_low <- resultsSysfit$xopt_95CI_high <- NA
    EstimatesTable <- data.frame(
      Parameter = names(resultsSysfit$xopt),
      Estimate  = sprintf("%.3f",resultsSysfit$xopt),
      stringsAsFactors = FALSE
    )

    # Append estimates with confidence intervals if estimated:
    estPars <- names(resultsSysfit$estObject$parameters$estimate)[resultsSysfit$estObject$parameters$estimate == 1]
    for (estk in estPars) {
      idxk <- EstimatesTable$Parameter == estk
      EstimatesTable$Estimate[idxk] <- paste0(EstimatesTable$Estimate[idxk], " (",sprintf("%.3f - %.3f", resultsSysfit$xopt_95CI_low[idxk], resultsSysfit$xopt_95CI_high[idxk]),")")
    }

    # Nicer error parameter names:
    EstimatesTable$Parameter <- aux_formatErrorName(EstimatesTable$Parameter)

  }else if (modelType=="NLME"){
    # Load results:
    resultsNLME <- getResults_IQRnlmeProjectMulti(as_IQRnlmeProjectMulti(modelFolder))[[1]]

    EstimatesTable <- data.frame(
      Parameter = resultsNLME$parameternames,
      Estimate  = ifelse(resultsNLME$stderrors==0 & !is.na(resultsNLME$stderrors),
                         sprintf("%.3f",resultsNLME$parametervalues),
                         ifelse(is.na(resultsNLME$stderrors),
                                paste0(sprintf("%.3f",resultsNLME$parametervalues), " (NaN%)"),
                                paste0(sprintf("%.3f",resultsNLME$parametervalues), " (",sprintf("%.1f",resultsNLME$stderrors/resultsNLME$parametervalues*100),"%)"))
      ),
      stringsAsFactors = FALSE
    )

  }else if (modelType=="sysFIT"){
    # Load results:
    resultsSysFIT <- getResults_IQRnlmeProject(as_IQRnlmeProject(modelFolder))
    #   Fixed Effect
    resFixedEffects <- resultsSysFIT[["fixedEffects"]]
    resFixedEffects$stderr[is.na(resFixedEffects$stderr) & resFixedEffects$estimated==0] <- 0
    EstimatesTable.fixedEffects <- data.frame(
      Parameter = resFixedEffects$names,
      Estimate  = ifelse(resFixedEffects$stderr==0 & !is.na(resFixedEffects$stderr),
                         sprintf("%.3f",resFixedEffects$values),
                         ifelse(is.na(resFixedEffects$stderr),
                                paste0(sprintf("%.3f",resFixedEffects$values), " (NaN%)"),
                                paste0(sprintf("%.3f",resFixedEffects$values), " (",sprintf("%.1f",resFixedEffects$stderr/resFixedEffects$values*100),"%)"))
      ),
      stringsAsFactors = FALSE
    )
    #   Random effect:
    resRandomEffects <- resultsSysFIT[["randomEffects"]]
    resRandomEffects$stderr[is.na(resRandomEffects$stderr) & resRandomEffects$estimated==0] <- 0
    EstimatesTable.randomEffects <- data.frame(
      Parameter = resRandomEffects$names,
      Estimate  = ifelse(resRandomEffects$stderr==0 & !is.na(resRandomEffects$stderr),
                         sprintf("%.3f",resRandomEffects$values),
                         ifelse(is.na(resRandomEffects$stderr),
                                paste0(sprintf("%.3f",resRandomEffects$values), " (NaN%)"),
                                paste0(sprintf("%.3f",resRandomEffects$values), " (",sprintf("%.1f",resRandomEffects$stderr/resRandomEffects$values*100),"%)"))
      ),
      stringsAsFactors = FALSE
    )
    #   Error Parameters:
    resErrorPar <- resultsSysFIT[["errorParameter"]]
    resErrorPar$stderr[is.na(resErrorPar$stderr) & resErrorPar$estimated==0] <- 0
    EstimatesTable.errorPar <- data.frame(
      Parameter = resErrorPar$names,
      Estimate  = ifelse(resErrorPar$stderr==0 & !is.na(resErrorPar$stderr),
                         sprintf("%.3f",resErrorPar$values),
                         ifelse(is.na(resErrorPar$stderr),
                                paste0(sprintf("%.3f",resErrorPar$values), " (NaN%)"),
                                paste0(sprintf("%.3f",resErrorPar$values), " (",sprintf("%.1f",resErrorPar$stderr/resErrorPar$values*100),"%)"))
      ),
      stringsAsFactors = FALSE
    )

    # Bind EstimatesTableS:
    EstimatesTable <- rbind(EstimatesTable.fixedEffects,
                            EstimatesTable.randomEffects,
                            EstimatesTable.errorPar)
  }

  # Output:
  return(EstimatesTable)
}
