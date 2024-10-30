#' Check Activity Informatiom
#'
#' Automatic check that the information provided by the modeler are valid.
#'
#' @param ActivityInfoFile Path to the activity information file (Default: \code{'../01-Data/S001-ActivityInfo'}).
#' @param tagFile          Path to the tag file  saved from PiNK (Default: \code{'../tags.RData'}).
#' @param FLAGMMVazure     To indicate if MMVazure, therefore PiNK, is being used or not (Default: \code{"MMVazure" %in% .packages(all.available = TRUE)})
#'
#' @return Message indicating if valid or not.
#'
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
check_ActivityInfo <- function(ActivityInfoFile = "../01-Data/S001-ActivityInfo",
                               tagFile          = "../tags.RData",
                               FLAGMMVazure     = ("MMVazure" %in% .packages(all.available = TRUE))) {

  # The purpose of check_ActivityInfo.R function is to gurantee that the ActivityInfo was:
  #   - Generated
  #   - Contain the appropriate variables for its activity type

  if(FLAGMMVazure){
    check_ActivityInfo_PiNK(ActivityInfoFile = ActivityInfoFile,
                            tagFile          = tagFile)
  }else{
    check_ActivityInfo_Sdrive(ActivityInfoFile = ActivityInfoFile)
  }
}


#' Check Activity Information on PiNK
#'
#' Automatic check that the information provided by the modeler are valid for the activities on PiNK.
#'
#' @param ActivityInfoFile Path to the activity information file (Default: \code{'../01-Data/S001-ActivityInfo'}).
#' @param tagFile          Path to the tag file  saved from PiNK (Default: \code{'../tags.RData'}).

#'
#' @return Message indicating if valid or not.
#'
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
check_ActivityInfo_PiNK <- function(ActivityInfoFile = "../01-Data/S001-ActivityInfo",
                                    tagFile          = "../tags.RData") {

  # warning("It seems that 'MMVazure', therefore PiNK, are being used. The function 'check_ActivityInfo' has not been updated for PiNK yet.")
  # return("Activity info not check as PiNK is being used")

  #------------------------------------------------------------#
  # STEP 1: Check ActivityInfo & tags exists ----
  #------------------------------------------------------------#
  if (!file.exists(ActivityInfoFile) || !file.exists(tagFile)){
    stop("Please start initializing the activity with the script S001_Initialization.R.")
  }


  #------------------------------------------------------------#
  # STEP 2: Load Activity Info & tags ----
  #------------------------------------------------------------#
  load(ActivityInfoFile)  # This one is still needed for the compounds
  load(tagFile)


  #------------------------------------------------------------#
  # STEP 3: Get info on the variable to check ----
  #------------------------------------------------------------#
  # NOTE: Done by PiNK now


  #------------------------------------------------------------#
  # STEP 4: Check if the existence of variables ----
  #------------------------------------------------------------#
  # NOTE: Done by PiNK now


  #------------------------------------------------------------#
  # STEP 5: Check Compound ----
  #------------------------------------------------------------#
  # Check for Compound/Compound1/Compound2/etc...
  # As explained above, the number of compound can't be know in advance, therefore,
  # the process of checking the number of CPD is different and will depend on
  # what is defined in ActivityType in ActivityProperties.xlsx

  # If 0, it means the user is free to use as many Compound name as he/she wants:
  ActivityTypeList <- list_ActivityType()
  ActivityType     <- tags[["Activity Type"]]
  if (0 %in% ActivityTypeList[[ActivityType]]$nbrCPD){
    if (length(ls(pattern="Compound"))==0){
      stop("It seems that you haven't defined at least one compound. Please do, e.g. Compound <- list(Name='MMVxxx'")
    }

    # ...Otherwise the number of compound has to respect the activity type properties
  } else{
    # Count the number of compound defined by the user:
    n_CPD <- length(ls(pattern = "Compound"))

    # Check that the number of compound is valid:
    if (n_CPD %in% ActivityTypeList[[ActivityType]]$nbrCPD){

      # if 1:
      if (n_CPD==1){
        if (!(exists("Compound"))){
          stop("Please define 'Compound'.")
        }

        # Check if it is well defined:
        #   Is List
        if (!is.list(Compound)){
          stop("'Compound' should be a list that contains at least the variable 'MMVname'.")
        }
        #   Does it contain 'MMVname'
        if (!("MMVname" %in% names(Compound))){
          stop("'MMVname' is not present in the list 'Compound': please add 'MMVname'.")
        }

        # if greater than 1:
      } else{
        for (k in 1:n_CPD){
          obj_k <- paste0("Compound", k)
          if (!(exists(obj_k))){
            stop("'", obj_k, "' is not defined.")
          }

          # Check if it is well defined:
          #   Is List
          if (!is.list(get(obj_k))){
            stop("'", obj_k, "' should be a list that contains at least the variable 'MMVname'.")
          }
          #   Does it contain 'MMVname'
          if (!("MMVname" %in% names(get(obj_k)))){
            stop("'MMVname' is not present in the list '", obj_k, "': please add 'MMVname'.")
          }
        }
      }

    }else{
      stop("The number of compound for activity type '",  ActivityType, "' should be ", collapseMMV(x = ActivityTypeList[[ActivityType]]$nbrCPD, andSymbole = " or "), ".")
    }
  }

  return("All activity info are valid.")
}


#' Check Activity Information on S Drive
#'
#' Automatic check that the information provided by the modeler are valid for the activities on the S Drive.
#'
#' @param ActivityInfoFile Path to the activity information file (Default: \code{'../01-Data/S001-ActivityInfo'}).
#'
#' @return Message indicating if valid or not.
#'
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
check_ActivityInfo_Sdrive <- function(ActivityInfoFile = "../01-Data/S001-ActivityInfo") {

  #------------------------------------------------------------#
  # STEP 1: Check ActivityInfo exists ----
  #------------------------------------------------------------#
  if (!file.exists(ActivityInfoFile)){
    stop("Please start initializing the activity with the script S001_Initialization.R.")
  }


  #------------------------------------------------------------#
  # STEP 2: Load Activity Info ----
  #------------------------------------------------------------#
  load(ActivityInfoFile)


  #------------------------------------------------------------#
  # STEP 3: Get info on the variable to check ----
  #------------------------------------------------------------#
  # NOTE:
  # 'Compound' will be check but treated differently as the number
  # of compounds can be different for each project

  # Path of the ActivityProperties.xlsx file:
  ActivityProperties.Path <- file.path(get_MMVmalariaPath(subdir = "inst"),"activityProperties/ActivityProperties.xlsx")

  # Get the names of all variables to be defined in the initialization:
  #   NOTE: Compound is not included as it will be an exception
  variableToCheck <- read_xlsx(ActivityProperties.Path, sheet = "Variables")$Initialization

  # Get the variable names with restricted values:
  variableToCheck.strictValue <- setdiff(excel_sheets(path = ActivityProperties.Path),
                                         c("Variables","ResultType"))

  # Check that all variables in 'variableToCheck.strictValue' are
  #   defined in 'variableToCheck':
  idx <- which(!(variableToCheck.strictValue %in% variableToCheck))
  if (length(idx)>0){
    stop("Add the variable(s) '", collapseMMV(x = variableToCheck.strictValue[idx]),"' to the sheet 'Variables$Initialization' in the file ActivityProperties.xlsx.")
  }


  #------------------------------------------------------------#
  # STEP 4: Check if the existence of variables ----
  #------------------------------------------------------------#
  # Variables defined in variableToCheck
  for (variable_k in variableToCheck){
    if (!(exists(variable_k))){
      stop("Please define '", variable_k, "' in the initialization script.")


      # Further check if variable_k is in 'variableToCheck.strictValue'
    }else if (variable_k %in% variableToCheck.strictValue){
      # Generate list of acceptable value for varaible_k:
      GenericTypeList <- list_GenericType(GenericType = variable_k)

      # Split Entery:
      GenericType.tmp <- strsplit(get(variable_k), "+", fixed = TRUE)[[1]]
      # Check if valid:
      if (!all(GenericType.tmp %in% names(GenericTypeList))){
        stop("'", variable_k, "' is not well defined. It should be equal to one of the following element: \n\t\t", collapseMMV(x = names(GenericTypeList), collapseSymbole="\n\t\t", andSymbole = "\n\t\t"))
      }
    }
  }


  #------------------------------------------------------------#
  # STEP 5: Check Compound ----
  #------------------------------------------------------------#
  # Check for Compound/Compound1/Compound2/etc...
  # As explained above, the number of compound can't be know in advance, therefore,
  # the process of checking the number of CPD is different and will depend on
  # what is defined in ActivityType in ActivityProperties.xlsx

  # If 0, it means the user is free to use as many Compound name as he/she wants:
  ActivityTypeList <- list_ActivityType()
  if (0 %in% ActivityTypeList[[ActivityType]]$nbrCPD){
    if (length(ls(pattern="Compound"))==0){
      stop("It seems that you haven't defined at least one compound. Please do, e.g. Compound<-list(Name='MMVxxx'")
    }

    # ...Otherwise the number of compound has to respect the activity type properties
  } else{
    # Count the number of compound defined by the user:
    n_CPD <- length(ls(pattern = "Compound"))

    # Check that the number of compound is valid:
    if (n_CPD %in% ActivityTypeList[[ActivityType]]$nbrCPD){

      # if 1:
      if (n_CPD==1){
        if (!(exists("Compound"))){
          stop("Please define 'Compound'.")
        }

        # Check if it is well defined:
        #   Is List
        if (!is.list(Compound)){
          stop("'Compound' should be a list that contains at least the variable 'MMVname'.")
        }
        #   Does it contain 'MMVname'
        if (!("MMVname" %in% names(Compound))){
          stop("'MMVname' is not present in the list 'Compound': please add 'MMVname'.")
        }

        # if greater than 1:
      } else{
        for (k in 1:n_CPD){
          obj_k <- paste0("Compound", k)
          if (!(exists(obj_k))){
            stop("'", obj_k, "' is not defined.")
          }

          # Check if it is well defined:
          #   Is List
          if (!is.list(get(obj_k))){
            stop("'", obj_k, "' should be a list that contains at least the variable 'MMVname'.")
          }
          #   Does it contain 'MMVname'
          if (!("MMVname" %in% names(get(obj_k)))){
            stop("'MMVname' is not present in the list '", obj_k, "': please add 'MMVname'.")
          }
        }
      }

    }else{
      stop("The number of compound for activity type '",  ActivityType, "' should be ", collapseMMV(x = ActivityTypeList[[ActivityType]]$nbrCPD, andSymbole = " or "), ".")
    }
  }

  return("All activity info are valid.")
}


#' Check Activity Results
#'
#' The function check if the results for the activity of interest are valid,
#' i.e. that the expected results for a given Activity Type are furnished.
#'
#' @param ActivityInfoFile    Path to an R object containing all the information of the activity of interest.
#' @param ActivityResultsFile Path to an R object containing all the results of the activity of interest.
#' @param tagFile             Path to the tag file  saved from PiNK (Default: \code{'../tags.RData'}).
#' @param FLAGMMVazure        To indicate if MMVazure, therefore PiNK, is being used or not (Default: \code{"MMVazure" %in% .packages(all.available = TRUE)})
#'
#' @return Message indicating if valid or not.
#'
#' @export
#' @family Activity Monitoring
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
check_ActivityResults <- function(ActivityInfoFile = "../01-Data/S001-ActivityInfo",
                                  ActivityResultsFile = "../01-Data/S999-ActivityResults",
                                  tagFile             = "../tags.RData",
                                  FLAGMMVazure        = ("MMVazure" %in% .packages(all.available = TRUE))) {

  #------------------------------------------------------------#
  # STEP 1: Check if ActivityInfo exists ----
  #------------------------------------------------------------#
  if (!file.exists(ActivityInfoFile)){
    stop("Please start initializing the activity with the script S001_Initialization.R.")

  }else{
    # Check Activity Info:
    #   Just in case that something changed since the beginning of the activity
    check_ActivityInfo(ActivityInfoFile=ActivityInfoFile)

    # Load Activity Info:
    load(ActivityInfoFile)
    if(FLAGMMVazure){
      load(tagFile)
      ActivityType <- tags[["Activity Type"]]
    }
  }


  #------------------------------------------------------------#
  # STEP 2: Check if ActivityResultsFile exists ----
  #------------------------------------------------------------#
  if (!file.exists(ActivityResultsFile)){
    stop("The file '", ActivityResultsFile, "' does not seem to exists: Please save ActivityResults in S999_ActivityFinalization.R prior using 'check_ActivityResults'.")

  }else{
    # Load Activity Results:
    load(ActivityResultsFile)
  }


  #------------------------------------------------------------#
  # STEP 3: Check if the activity was aborted ----
  #------------------------------------------------------------#
  if (any(names(ActivityResults)=="Aborted")){
    # Keep only the 'Aborted' result:
    ActivityResults <- ActivityResults["Aborted"]
    save(list = c("ActivityResults", "FinalDate"),
         file = ActivityResultsFile)
    return("Activity Aborted")
  }


  #------------------------------------------------------------#
  # STEP 4: Check if results are needed and missing ----
  #------------------------------------------------------------#

  # Get list_ActivityType:
  ActivityTypeList <- list_ActivityType()
  # Check what type of results is expected:
  if (all(ActivityTypeList[[ActivityType]]$results=="free")){
    # Results=free => any type of results is allowed and no check needed
    message(paste0("The activity type '", ActivityType,"' allows to define any type of results."))


    # If activity type results were not defined yet:
  } else if(all(ActivityTypeList[[ActivityType]]$results=="TBD")){
    warning(paste0("The results for the activity type '", ActivityType,"' were not defined yet. Please update 'ActivityProperties.xlsx'."))


    # If Results!=free & !=TBD, therefore, specific results are expected:
  } else{

    # Check that all the minimum results corresponding to the activity type are well defined:
    if (!(all(ActivityTypeList[[ActivityType]]$results %in% names(ActivityResults)))){
      idx_FALSE <- which(!(ActivityTypeList[[ActivityType]]$results %in% names(ActivityResults)))
      if (length(idx_FALSE)>0){
        stop(cat(paste0("Some results related to the activity type: '", ActivityType,"' are missing.\nPlease add the following variables:"), ActivityTypeList[[ActivityType]]$results[idx_FALSE],sep="\n\t"))
      }
    }
  }


  #------------------------------------------------------------#
  # STEP 5: Check if the results are properly defined ----
  #------------------------------------------------------------#

  # Path of the ActivityProperties.xlsx file:
  ActivityProperties.Path <- file.path(get_MMVmalariaPath(subdir = "inst"),"activityProperties/ActivityProperties.xlsx")

  # Get the function to test the result validity:
  resultToCheck <- read_xlsx(ActivityProperties.Path, sheet = "ResultType")

  # Check if the path of each result exists and if a valid result:
  for (result_k in names(ActivityResults)){

    # Check if file(s) exist(s):
    #   Test each element of the list
    TestExistence <- sapply(ActivityResults[[result_k]], function(x){
      if(!(dir.exists(x) | file.exists(x))){
        FALSE
      }else{
        TRUE
      }
    })
    #   Read outcome
    if(any(!TestExistence)){
      stop("The path '",
           collapseMMV(x = ActivityResults[[result_k]][!TestExistence], collapseSymbole = "', '", andSymbole = "' & '"),
           "' for the '", result_k, "'  does not exist: Please define the appropriate path.")
    }


    # Check if valid result for reserved ResultType:
    if(result_k %in% resultToCheck$ResultType){
      # Get index of result_k in resultToCheck$ResultType data frame:
      idx_result <- which(resultToCheck$ResultType==result_k)

      # Test each element of the list:
      TestResultFormat <- sapply(ActivityResults[[result_k]], function(x){
        # Generate the text to test:
        TestText <- gsub("ResultX",
                         x,
                         resultToCheck$FunctionTest[[idx_result]])

        # Test if valid result format:
        if(!eval(parse(text=TestText))){
          FALSE
        }else{
          TRUE
        }
      })

      # Read outcome:
      if(any(!TestResultFormat)){
        stop("The path '",
             collapseMMV(x = ActivityResults[[result_k]][!TestResultFormat], collapseSymbole = "', '", andSymbole = "' & '"),
             "' for the '", result_k, "' was expected to pass the test: '", resultToCheck$FunctionTest[[idx_result]], "'")
      }
    }
  }

  return("All activity results are valid")
}


#' create_RprofileMMV
#'
#' @description
#' @param IQRversion Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
create_RprofileMMV <- function(IQRversion = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Initialize ----
  #-----------------------------------------------------------------------------#

  # Check if IQRversion was set up:
  if (is.null(IQRversion)){
    stop("Please define 'IQRversion'.")
  }

  # If Unix server:
  if (.Platform$OS.type=="unix" & IQRversion<"0.9.2"){
    stop("Only IQRtools above version '0.9.2' should be used with the unix server to avoid messing up with 'RProfile'.")
  } else if (.Platform$OS.type=="unix" & IQRversion>="0.9.2"){
    # Stop function:
    return(cat("\n------------\n   The IQRversion is >='0.9.2', therefore, there is no need in modifying 'RProfile'.\nFurthremore, no control file was added as it is running on the unix server.\n------------\n"))
  }

  # Get current library Path:
  currentLibPaths <- .libPaths()

  # Reset .libPaths():
  invisible(.libPaths(c(unlist(strsplit(Sys.getenv("R_LIBS"), ";")),
                        unlist(strsplit(Sys.getenv("R_LIBS_USER"), ";"))
  )))

  # Session Information:
  UserName     <- Sys.info()[["user"]]
  thisFileName <- thisFile()
  thisFileName <- gsub(".R", "", thisFileName, fixed=TRUE)

  # Locate the file RProfile:
  baseR    <- file.path(.libPaths(),"base/R")
  RProfile <- file.path(baseR,"Rprofile")
  RProfile_Original <- paste0(RProfile,"_Original")


  #-----------------------------------------------------------------------------#
  # STEP 1: IQRversion>="0.9.2" ----
  #-----------------------------------------------------------------------------#
  if (IQRversion>="0.9.2"){
    if (file.exists(RProfile_Original)){
      stop("'RProfile' was modified by a previous version of IQRtools or R is running with a different 'IQRtools' version. Please reset 'RProfile' to its original form and delete 'RProfile_Original'.")
    } else{

      # Add Controle to avoid having other IQRtools running at the same time and modifying RProfile:
      RunningID   <- paste0("Running_IQRtools_v", IQRversion, "_by_", UserName, "_from_", thisFileName)
      TextRunning <- paste0(UserName, " is using ", IQRversion, ". Please do not delete RProfile while ", UserName, " is still using R.")
      writeLines(TextRunning, file.path(baseR, RunningID))

      # Reset library Path:
      .libPaths(currentLibPaths)

      # Stop function:
      return(cat("\n------------\n   The IQRversion is >='0.9.2', therefore, there is no need in modifying 'RProfile'.\n------------\n"))
    }
  }


  #-----------------------------------------------------------------------------#
  # STEP 3: Check if there is another user with a different IQRtools version ----
  #-----------------------------------------------------------------------------#

  # Running file:
  idx_runningFile <- which(grepl("Running_IQRtools_v",list.files(baseR)))
  runningFiles    <- list.files(baseR)[idx_runningFile]

  # various IQRtools version:
  allIQRtoolsVersion <- unique(substr(runningFiles,19,23))
  if (length(allIQRtoolsVersion)>1){
    stop("There are more than one IQRtools version defined. Please delete the 'Running files manually or talk to Mohammed")
  }

  # if there is a different version, stop:
  if (any(allIQRtoolsVersion!=IQRversion)){
    stop("Another IQRtools version is running. Please wait before running your script and/or talk to your colleagues.")
  }


  #-----------------------------------------------------------------------------#
  # STEP 4: Create new Rprofile ----
  #-----------------------------------------------------------------------------#

  # Locate the file RProfile_Original:
  #   Here its assumes that if it does not exist
  #   RProfile will be the original file
  RProfile_Original <- paste0(RProfile,"_Original")
  if (!file.exists(RProfile_Original)){
    file.copy(RProfile, RProfile_Original)
  }

  # Open RProfile_Original:
  Text <- readLines(RProfile_Original)

  # Update Text in RProfile_Original:
  Text <- c(Text,
            "",
            "#~~~~~~~~~",
            paste0("# Lines Added by ", UserName, ":"),
            "#   Reset .libPath()",
            "invisible(.libPaths(c(unlist(strsplit(Sys.getenv(\"R_LIBS\"), \";\")),",
            "                      unlist(strsplit(Sys.getenv(\"R_LIBS_USER\"), \";\"))",
            ")))",
            "#   Update .libPath() for all R instances",
            paste0(".libPaths(\"", libraryPath, "\")"))

  # Save Updates:
  writeLines(Text, RProfile)


  #-----------------------------------------------------------------------------#
  # STEP 5: Create new Rprofile ----
  #-----------------------------------------------------------------------------#

  # Add Control if multiple user:
  RunningID   <- paste0("Running_IQRtools_v", IQRversion, "_by_", UserName, "_from_", thisFileName)
  TextRunning <- paste0(UserName, " is using ", IQRversion, ". Please do not delete RProfile while ", UserName, " is still using R.")
  writeLines(TextRunning, file.path(baseR, RunningID))

  # Reset library Path:
  .libPaths(currentLibPaths)
}
#' generate_ActivityContent
#'
#' @description
#' @param ActivityList
#' @param colADD Default: c("ProjectName", "CompoundSummary", "ActivityName", "ActivityType",
#'    "ActivityPath", "ModellerName", "Status", "StartDate", "LastRunDay",
#'    "FinalDate", "StudyType", "Species", "SubjectType", "ParasiteType",
#'    "InfectionType", "IndicationType", "ActivityResultsSummary",
#'    "reviewDate", "revierName", "workflowName", "workflowVersion",
#'    "Rversion", "IQRversion", "tool", "toolVersion", "ActivityDescription")
#' @param filename Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
generate_ActivityContent <- function(ActivityList,  # Always use get_ActivityInfo to have it in the right format
                                     colADD = c("ProjectName", "CompoundSummary", "ActivityName", "ActivityType",
                                                "ActivityPath", "ModellerName", "Status", "StartDate",
                                                "LastRunDay", "FinalDate", "StudyType", "Species",
                                                "SubjectType", "ParasiteType", "InfectionType", "IndicationType",
                                                "ActivityResultsSummary", "reviewDate", "revierName", "workflowName",
                                                "workflowVersion", "Rversion", "IQRversion", "tool",
                                                "toolVersion", "ActivityDescription"),
                                     filename = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Variable to use for the data Frame ----
  #-----------------------------------------------------------------------------#

  # Variable that should always be present
  colADD_Default <-  c("ProjectName", "CompoundSummary", "ActivityName", "ActivityPath")

  # Add colADD_Default to colADD
  colADD <- unique(c(colADD_Default, colADD))

  # Define List Name:
  colADD_List <- paste0(colADD, "_List")


  #-----------------------------------------------------------------------------#
  # STEP 2: Generate Activity Summary as a data frame ----
  #-----------------------------------------------------------------------------#

  # Loop over colADD to add information
  ActivitiesContent <- NULL
  for (k in seq_along(colADD)){
    # Variable name to add:
    colADD_k      <- colADD[k]
    colADD_List_k <- colADD_List[k]

    # Temporary Data Frame:
    data_temp           <- data.frame(t(as.data.frame(ActivityInfo[[colADD_List_k]])),
                                      stringsAsFactors = FALSE)
    names(data_temp)    <- colADD_k
    rownames(data_temp) <- NULL

    # Concatenate:
    if (is.null(ActivitiesContent)){
      ActivitiesContent <- data_temp
    }else{
      ActivitiesContent <- cbind(ActivitiesContent, data_temp)
    }
  }

  # Adjust Column Name:
  names(ActivitiesContent) <- gsub("Summary", "", names(ActivitiesContent))


  #-----------------------------------------------------------------------------#
  # STEP 3: Output ----
  #-----------------------------------------------------------------------------#

  # Generate CSV file
  if (!is.null(filename)){
    IQRoutputCSV(ActivitiesContent,
                 filename = filename,
                 replaceComma = ";")
  }

  # Return output:
  return(ActivitiesContent)
}
#' generate_ActivityResults
#'
#' @description
#' @param ActivityInfo
#' @param colADD Default: c("ProjectName", "CompoundSummary", "ActivityName", "ActivityType",
#'    "ActivityPath", "FinalDate", "Status", "StudyType", "Species",
#'    "SubjectType", "ParasiteType", "InfectionType", "IndicationType",
#'    "workflowName", "workflowVersion", "Rversion", "IQRversion",
#'    "tool", "toolVersion")
#' @param filename Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
generate_ActivityResults <- function(ActivityInfo,  # Always use get_ActivityInfo to have it in the right format
                                     colADD = c("ProjectName", "CompoundSummary", "ActivityName", "ActivityType",
                                                "ActivityPath", "FinalDate", "Status", "StudyType", "Species",
                                                "SubjectType", "ParasiteType", "InfectionType", "IndicationType",
                                                "workflowName", "workflowVersion", "Rversion", "IQRversion", "tool",
                                                "toolVersion"),
                                     filename = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Variable to use for the data Frame ----
  #-----------------------------------------------------------------------------#

  # Variable that should always be present
  colADD_Default <-  c("ProjectName", "CompoundSummary", "ActivityName", "ActivityType", "ActivityPath", "FinalDate", "Status")

  # Add colADD_Default to colADD
  colADD <- unique(c(colADD_Default, colADD))

  # Define List Name:
  colADD_List <- paste0(colADD, "_List")


  #-----------------------------------------------------------------------------#
  # STEP 2: Generate Activity Summary as a data frame ----
  #-----------------------------------------------------------------------------#

  # Loop over colADD to add information
  ActivitiesContent <- data.frame(Name = names(ActivityInfo[["ProjectName_List"]]),
                                  stringsAsFactors = FALSE)
  for (k in seq_along(colADD)){
    # Variable name to add:
    colADD_k      <- colADD[k]
    colADD_List_k <- colADD_List[k]

    # Temporary Data Frame:
    data_temp          <- data.frame(t(as.data.frame(ActivityInfo[[colADD_List_k]])),
                                     stringsAsFactors = FALSE)
    names(data_temp)    <- colADD_k
    rownames(data_temp) <- NULL

    # Concatenate:
    ActivitiesContent <- cbind(ActivitiesContent, data_temp)
  }

  # Adjust Column Name:
  names(ActivitiesContent) <- gsub("Summary", "", names(ActivitiesContent))

  # Keep only finalized activities:
  ActivitiesContent <- ActivitiesContent[ActivitiesContent$Status=="Finalized",]


  #-----------------------------------------------------------------------------#
  # STEP 3: Generate Activity Results as a data frame ----
  #-----------------------------------------------------------------------------#

  ActivitiesResult <- NULL
  for (Name_k in ActivitiesContent$Name){
    # Get Results:
    Result_k <- ActivityInfo$ActivityResults_List[[Name_k]]

    # Temporary Data Frame:
    data_temp <- data.frame(Name       = Name_k,
                            ResultType = names(Result_k),
                            Result     = unlist(Result_k),
                            stringsAsFactors = FALSE)
    row.names(data_temp) <- NULL

    # Concatenate:
    if (is.null(ActivitiesResult)){
      ActivitiesResult <- data_temp
    }else{
      ActivitiesResult <- rbind(ActivitiesResult, data_temp)
    }
  }

  # Concatenate ActivitiesContent & ActivitiesResult:
  ActivitiesResult <- merge(ActivitiesContent, ActivitiesResult)

  # Remove Name Column:
  ActivitiesResult$Name <- NULL


  #-----------------------------------------------------------------------------#
  # STEP 4: Output ----
  #-----------------------------------------------------------------------------#

  # Generate CSV file:
  if (!is.null(filename)){
    IQRoutputCSV(ActivitiesResult,
                 filename = filename)
  }

  # Return output:
  return(ActivitiesResult)
}

#' generate_DuplicatedActivityReport
#'
#' @description
#' @param Activity_List
#' @param filename Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
generate_DuplicatedActivityReport <- function(Activity_List,
                                              filename = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Generate Report ----
  #-----------------------------------------------------------------------------#

  if (length(Activity_List$DuplicatedActivity)>0) {
    # Generate message to save
    Msg <- c()
    for (Activity_k in Activity_List$DuplicatedActivity){
      Msg <- c(Msg, sprintf("The activity overview file %s seems to be a duplicate.", Activity_k))
      Msg <- c(Msg, "")
    }

  } else{
    Msg <- "No activity overview file is duplicated."
  }


  #-----------------------------------------------------------------------------#
  # STEP 2: Generate Output ----
  #-----------------------------------------------------------------------------#

  # Generate Report:
  if (!is.null(filename)){
    IQRoutputFile(Msg,
                  filename = filename)
  }

  # Return output:
  return(Msg)
}
#' generate_MissingActivityReport
#'
#' @description
#' @param Activity_List
#' @param filename Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
generate_MissingActivityReport <- function(Activity_List,
                                           filename = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Generate Report ----
  #-----------------------------------------------------------------------------#

  if (length(Activity_List$MissingActivity)>0) {
    # Generate message to save
    Msg <- c()
    for (Activity_k in names(Activity_List$MissingActivity)){
      Msg <- c(Msg, sprintf("Activity %s is missing.", Activity_k))
      Msg <- c(Msg, sprintf("     It is supposed to be located in %s.", Activity_List$MissingActivity[[Activity_k]]))
      Msg <- c(Msg, "")
    }

  } else{
    Msg <- "No activity is missing."
  }


  #-----------------------------------------------------------------------------#
  # STEP 2: Generate Output ----
  #-----------------------------------------------------------------------------#

  # Generate Report:
  if (!is.null(filename)){
    IQRoutputFile(Msg,
                  filename = filename)
  }

  # Return output:
  return(Msg)
}
#' generate_OverviewFile
#'
#' @description
#' @param outputOverview Default: '../../../../../Content/00-ActivityList'
#' @param ActivityName Default: NULL
#' @param ProjectName Default: NULL
#' @param ActivityPath Default: NULL
#' @param NewOverviewFolder Default: FALSE
#' @param FLAGforceOverview Default: FALSE
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
generate_OverviewFile <- function(outputOverview = "../../../../../Content/00-ActivityList",
                                  ActivityName   = NULL,
                                  ProjectName    = NULL,
                                  ActivityPath   = NULL,
                                  NewOverviewFolder = FALSE,
                                  FLAGforceOverview = FALSE) {



  # Generate Activity Information:
  if (!dir.exists("../00-Toolbox/MMVmalaria") | FLAGforceOverview) {
    #----------------------------------------------#
    # STEP 0: Check Variables: ----
    #----------------------------------------------#
    if (is.null(outputOverview)){
      stop("Please define 'outputOverview'.")
    }

    #----------------------------------------------#
    # STEP 1: Get Current Working Directory ----
    #----------------------------------------------#
    WorkDir       <- getwd()
    WorkDir_Split <- strsplit(WorkDir, "/", fixed = TRUE)[[1]]
    n_str         <- length(WorkDir_Split)


    #----------------------------------------------#
    # STEP 2: Define ActivityName if NULL ----
    #----------------------------------------------#
    if (is.null(ActivityName) | WorkDir_Split[n_str-4]=="Projects"){
      # Activity Name:
      ActivityName  <- WorkDir_Split[n_str-1]
    }


    #---------------------------------------------------------------#
    # STEP 3: Over-write ProjectName for activities in Projects----
    #---------------------------------------------------------------#
    # If it is an activity defined in the project folder it is straight forward
    if (WorkDir_Split[n_str-4]=="Projects" && WorkDir_Split[n_str-5]!="InVitroCombo"){
      ProjectName   <- WorkDir_Split[n_str-3]

    } else if(WorkDir_Split[n_str-4]=="Projects" && WorkDir_Split[n_str-5]=="InVitroCombo"){
      ProjectName   <- paste0("InVitroCombo_", WorkDir_Split[n_str-3])

    } else if(is.null(ProjectName)){
      stop("Please define 'ProjectName'.")
    }


    #----------------------------------------------#
    # STEP 4: Define ActivityPath if NULL ----
    #----------------------------------------------#
    if (is.null(ActivityPath)){
      # Count the number of '..' in outputOverview:
      n_2dots     <- (nchar(outputOverview) - nchar(gsub("..","",outputOverview,fixed=TRUE)))/2

      # From the current working directory, only keep the last n_2dots folder but last:
      WorkDir_Split_Short <- WorkDir_Split[(n_str-(n_2dots-1)):(n_str-1)]

      # Count the number of sub folder in outputOverview: neglecting the numbe rof '..'
      n_temp   <- (nchar(outputOverview) - nchar(gsub("/","",outputOverview,fixed=TRUE)))
      n_folder <- n_temp+1-n_2dots

      # Generate the relative path of the activity relative the location of outputOverview
      ActivityPath <- ""
      for (k in 1:(n_folder)){
        ActivityPath <- paste0(ActivityPath, "../")
      }
      ActivityPath <- paste0(ActivityPath, paste0(WorkDir_Split_Short, collapse = '/'))
    }


    #----------------------------------------------#
    # STEP 5: Save Activity's Info ----
    #----------------------------------------------#

    # Check if outputOverview exist:
    if (NewOverviewFolder){
      if (!dir.exists(outputOverview)) dir.create(outputOverview,recursive=TRUE)
    } else{
      if (!dir.exists(outputOverview)) stop("The path of 'outputOverview' does not exist. Create the folder manually or check that the path 'outputOverview' is correct.")
    }

    # Get Start Date:
    ActivityInfoFile <- "../01-Data/S001-ActivityInfo"
    if (file.exists(ActivityInfoFile)){
      StartDate <- paste0(substr(file.info(ActivityInfoFile)$ctime, 1, 4), "/", substr(file.info(ActivityInfoFile)$ctime, 1, 7))
    } else{
      StartDate <- paste0(substr(Sys.Date(), 1, 4), "/", substr(Sys.Date(), 1, 7))
    }

    # Create Sub-Folder:
    if (!dir.exists(file.path(outputOverview, StartDate))) dir.create(file.path(outputOverview, StartDate),recursive=TRUE)


    # Save information of interest for Overview:
    save(list = c("ActivityName","ProjectName","ActivityPath"),
         file = file.path(outputOverview, StartDate, paste0(ProjectName,"_",ActivityName)))
  }
}
#' get_ActivityInfo
#'
#' @description
#' @param Activity_List
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
get_ActivityInfo <- function(Activity_List) {

  #-----------------------------------------------------------------------------#
  # STEP 1-a: Variables defined in ActivityProperties.xlsx ----
  #-----------------------------------------------------------------------------#

  # Path of the ActivityProperties.xlsx file:
  ActivityProperties.Path <- file.path(get_MMVmalariaPath(subdir = "inst"),"activityProperties/ActivityProperties.xlsx")

  # Get the names of all variables defined in the initialization:
  #   NOTE: Compound is not included as it will be an exception
  variableToLoad <- read_xlsx(ActivityProperties.Path, sheet = "Variables")$Initialization
  variableToLoad_List <- paste0(variableToLoad, "_List")
  for (variableToLoad_List_k in variableToLoad_List){
    assign(variableToLoad_List_k, list())
  }

  # TO ADD: Finalization & Review variables
  #.........


  #-----------------------------------------------------------------------------#
  # STEP 1-b: Variables NOT defined in ActivityProperties.xlsx YET----
  #-----------------------------------------------------------------------------#

  # Defined variable to be loaded not defined in ActivityProperties.xlsx:
  manualVariableToLoad <- c("Status", "LastRunDay", "FinalDate", "reviewDate", "revierName")
  manualVariableToLoad_List <- paste0(manualVariableToLoad, "_List")
  for (manualVariableToLoad_List_k in manualVariableToLoad_List){
    assign(manualVariableToLoad_List_k, list())
  }

  # Compound List as it is treated differently:
  COMPOUND_List        <- list()
  COMPOUNDSummary_List <- list()

  # Activity Result List as it is treated differently:
  ActivityResults_List     <- list()
  ActivityResultsSummary_List <- list()


  #-----------------------------------------------------------------------------#
  # STEP 2: Load All Variables ----
  #-----------------------------------------------------------------------------#

  for (Activity_k in Activity_List$ValidActivity){
    # Activity Path:
    ActivityPath_List_k <- Activity_List$ActivityPath_List[[Activity_k]]

    # Load Activity Info:
    load(file.path(ActivityPath_List_k,"01-Data/S000-LatestUpdate"))
    load(file.path(ActivityPath_List_k,"01-Data/S001-ActivityInfo"))

    # Load Activity Results if exist
    if (file.exists(file.path(ActivityPath_List_k,"01-Data/S999-ActivityResults"))){
      load(file.path(ActivityPath_List_k,"01-Data/S999-ActivityResults"))
    }


    # Load variables defined in 'variableToLoad':
    for (i in seq_along(variableToLoad)){
      # Variable & List i:
      variableToLoad_i      <- variableToLoad[i]
      variableToLoad_List_i <- variableToLoad_List[i]

      # Update List:
      ls_Temp <- get(variableToLoad_List_i, inherits = FALSE)
      if (exists(variableToLoad_i, inherits = FALSE)){
        ls_Temp[[Activity_k]] <- get(variableToLoad_i, inherits = FALSE)
      }else{
        ls_Temp[[Activity_k]] <- "NA"
      }
      assign(variableToLoad_List_i, ls_Temp)
    }

    # Adjust Activity Description:
    ActivityDescription_List[[Activity_k]] <- gsub("\n", "", ActivityDescription_List[[Activity_k]])


    # Load variables defined in 'manualVariableToLoad':
    for (i in seq_along(manualVariableToLoad)){
      # Variable & List i:
      manualVariableToLoad_i      <- manualVariableToLoad[i]
      manualVariableToLoad_List_i <- manualVariableToLoad_List[i]

      # Update List:
      ls_Temp <- get(manualVariableToLoad_List_i, inherits = FALSE)
      if (exists(manualVariableToLoad_i, inherits = FALSE)){
        ls_Temp[[Activity_k]] <- get(manualVariableToLoad_i, inherits = FALSE)
      }else{
        ls_Temp[[Activity_k]] <- ""
      }
      assign(manualVariableToLoad_List_i, ls_Temp)
    }


    # Compound:
    #   Get list of the compounds:
    COMPOUND_List[[Activity_k]] <- list()
    COMPOUNDSummary_tmp <- c()
    for (i in 1:length(ls(pattern = "Compound"))){
      if ("Name" %in% names(get(ls(pattern = "Compound")[i], inherits = FALSE))){
        COMPOUND_List[[Activity_k]][[i]] <- get(ls(pattern = "Compound")[i], inherits = FALSE)
        COMPOUNDSummary_tmp              <- c(COMPOUNDSummary_tmp,
                                              get(ls(pattern = "Compound")[i])$Name)
      }else{
        COMPOUND_List[[Activity_k]][[i]] <- get(ls(pattern = "Compound")[i], inherits = FALSE)
        COMPOUNDSummary_tmp              <- c(COMPOUNDSummary_tmp,
                                              get(ls(pattern = "Compound")[i], inherits = FALSE))
      }
    }
    #   Adjust list name:
    if(length(ls(pattern = "Compound"))==1){
      names(COMPOUND_List[[Activity_k]]) <- "Compound"
    }else{
      names(COMPOUND_List[[Activity_k]]) <- paste0("Compound", seq_along(ls(pattern = "Compound")))
    }
    #   Get a summary of the name
    COMPOUNDSummary_List[[Activity_k]] <- collapseMMV(x               = COMPOUNDSummary_tmp,
                                                      collapseSymbole = "+",
                                                      andSymbole      = "+")


    # Activity Results:
    if (exists("ActivityResults", inherits = FALSE)){
      # Get List of Results:
      ActivityResults_List[[Activity_k]] <- ActivityResults

      # Create Result Summary:
      ActivityResults_tmp <- ActivityResults
      for (k in 1:length(ActivityResults)){
        ActivityResults_tmp[[k]] <- paste0(names(ActivityResults)[k], ": ", ActivityResults[[k]])
      }

      # Assign Summary Result:
      ActivityResultsSummary_List[[Activity_k]] <- collapseMMV(x               = ActivityResults_tmp,
                                                               collapseSymbole = "+",
                                                               andSymbole      = "+")
    }else{
      ActivityResults_List[[Activity_k]]        <- NA
      ActivityResultsSummary_List[[Activity_k]] <- ""
    }


    # Remove some variables:
    obj <- c(variableToLoad, manualVariableToLoad, "ls_Temp",
             "COMPOUNDSummary_tmp", "ActivityResults", "ActivityResults_tmp")
    obj <- c(obj, ls(pattern = "Compound"))
    for (obj_i in obj){
      if(exists(obj_i, inherits = FALSE)) {
        rm(list = obj_i, inherits = FALSE)
      }
    }
  }


  #-----------------------------------------------------------------------------#
  # STEP 3: Output ----
  #-----------------------------------------------------------------------------#

  # Create Output:
  out <- list()

  # Add Activity List Info:
  out[["ActivityPath_List"]] <- Activity_List$ActivityPath_List
  out[["ProjectName_List"]]  <- Activity_List$ProjectName_List
  out[["ActivityName_List"]] <- Activity_List$ActivityName_List
  out[["ValidActivity"]]     <- Activity_List$ValidActivity

  # Add variables define in 'variableToLoad_List':
  for (variableToLoad_List_k in variableToLoad_List){
    out[[variableToLoad_List_k]] <- get(variableToLoad_List_k, inherits = FALSE)
  }

  # Add variables define in 'manualVariableToLoad_List':
  for (manualVariableToLoad_List_k in manualVariableToLoad_List){
    out[[manualVariableToLoad_List_k]] <- get(manualVariableToLoad_List_k, inherits = FALSE)
  }

  # Add Compound Lists:
  out[["Compound_List"]]        <- COMPOUND_List
  out[["CompoundSummary_List"]] <- COMPOUNDSummary_List

  # Add Results Lists:
  out[["ActivityResults_List"]]        <- ActivityResults_List
  out[["ActivityResultsSummary_List"]] <- ActivityResultsSummary_List

  # Return output:
  return(out)
}
#' get_ActivityList
#'
#' @description
#' @param folder_List Default: '../00-ActivityList'
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
get_ActivityList <- function(folder_List = "../00-ActivityList") {

  # List of all activities
  Activity_List <- list.files(folder_List, recursive = TRUE)

  # Load the path/names of all activities:
  ActivityPath_List <- list()
  ProjectName_List  <- list()
  ActivityName_List <- list()
  ValidActivity     <- c()
  MissingActivity   <- list()
  for (Activity_k in Activity_List){

    # Load Activity Path:
    load(file.path(folder_List,Activity_k))

    # Check that the activity folder exist:
    if (file.exists(ActivityPath)){
      # Add info to lists:
      ActivityPath_List[[Activity_k]] <- ActivityPath
      ProjectName_List[[Activity_k]]  <- ProjectName
      ActivityName_List[[Activity_k]] <- ActivityName

      # Add to valid list of activities:
      ValidActivity <- c(ValidActivity, Activity_k)

    }else{
      # Add to missing list of activities:
      MissingActivity[[Activity_k]] <- ActivityPath
    }

    # Remove some variables:
    rm(list=c("ActivityPath", "ProjectName", "ActivityName"))
  }

  # Order Valid Activity with Project Name, then Activity Name:
  idx_Order <- order(unlist(ProjectName_List), unlist(ActivityName_List))
  ActivityPath_List <- ActivityPath_List[idx_Order]
  ProjectName_List  <- ProjectName_List[idx_Order]
  ActivityName_List <- ActivityName_List[idx_Order]
  ValidActivity     <- ValidActivity[idx_Order]

  # Check for duplicates:
  idx_Duplicate <- duplicated(ActivityPath_List)

  # Get Duplicated Activities:
  DuplicatedActivity <- ValidActivity[idx_Duplicate]

  # Update Valid Activities:
  ActivityPath_List <- ActivityPath_List[!idx_Duplicate]
  ProjectName_List  <- ProjectName_List[!idx_Duplicate]
  ActivityName_List <- ActivityName_List[!idx_Duplicate]
  ValidActivity     <- ValidActivity[!idx_Duplicate]

  # Create Output:
  out <- list(ActivityPath_List  = ActivityPath_List,
              ProjectName_List   = ProjectName_List,
              ActivityName_List  = ActivityName_List,
              ValidActivity      = ValidActivity,
              DuplicatedActivity = DuplicatedActivity,
              MissingActivity    = MissingActivity)

  # Return output:
  return(out)
}
#' list_ActivityType
#'
#' @description

#' @return
#' @export
#' @author Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
list_ActivityType <- function() {

  # Locate the file 'ActivityProperties.xlsx within MMVmalaria:
  ActivityProperties.Path <- file.path(get_MMVmalariaPath(subdir = "inst"),"activityProperties/ActivityProperties.xlsx")

  # Obtain the ActivityType sheet:
  ActivityType.Excel <- read_xlsx(ActivityProperties.Path, sheet = "ActivityType")

  # Generate ActivityType list:
  ActivityType <- list()
  for (ActivityType_k in ActivityType.Excel$ActivityType){
    # Get Index Number:
    k <- which(ActivityType.Excel$ActivityType==ActivityType_k)

    # Get nbrCPD:
    nbrCPD_k <- as.numeric(strsplit(ActivityType.Excel$nbrCPD[k], "+", fixed = TRUE)[[1]])

    # Get result:
    results_k <- strsplit(ActivityType.Excel$results[k], "+", fixed = TRUE)[[1]]

    # Add nbrCPD_k and result_K to Activity Type:
    ActivityType[[ActivityType_k]] <- list(nbrCPD  = nbrCPD_k,
                                           results = results_k)
  }

  # Return list:
  return(ActivityType)
}
#' list_Center
#' Loads center list
#'
#' Loads up-to-date center list to perform center name and number verification during data preparation, see [check_dataGenetalMMV_Center] and [check_dataGeneralMMV].
#' This list should be updated for each new center, adding the center name and number.
#'
#' @md
#'
#' @return valid list of center names and numbers
#' @export
#' @seealso [check_dataGenetalMMV_Center], [check_dataGeneralMMV]
#' @family Data Preparation
#' @author To be defined

list_Center <- function() {

  CenterList <- list("TBD"                   = -1,
                     "Publication"           = 0,
                     "UCT"                   = 101001,
                     "Q2Sol-Centurion"       = 101002,
                     "CREC"                  = 102001,
                     "CERPAGE"               = 102002,
                     "IRSS-Nanoro"           = 103001,
                     "CNRFP-Banfora"         = 103002,
                     "CNRFP-Niangoloko"      = 103003,
                     "MURAZ"                 = 103004,
                     "CRSN-Nouna"            = 103005,
                     "CNRFP-Ouagadougou"     = 103006,
                     "USS-Libreville"        = 104001,
                     "CERMEL"                = 104002,
                     "KEMRI-Kisumu"          = 105001,
                     "KEMRI-Kondele"         = 105002,
                     "KEMRI-Kilifi"          = 105003,
                     "KEMRI-WRP"             = 105004,
                     "KEMRI-Nairobi"         = 105005,
                     "KEMRI-Siaya"           = 105006,
                     "CITSC"                 = 106001,
                     "CISM"                  = 106002,
                     "IDRC-Tororo"           = 107001,
                     "XXX-Yaounde"           = 108001,
                     "NIMR-Korogwe"          = 109001,
                     "CRCHMA"                = 110001,
                     "Faculte de Medecine Kinshasa" = 110002,
                     "Clinical Pharmacology and Pharmacovigilance Unit Kinshasa" = 110003,
                     "Institut M�dicale �vang�lique de Kimpes�" = 110004,
                     "Centre de sant� FCRM - Massisia" = 110005,
                     "MRCU-Gambia"           = 111001,
                     "MRC Laboratory Farafeni" = 111002,
                     "Albert Schweitzer Hospital Libreville" = 111003,
                     "Nchelenge-Zambia"      = 112001,
                     "Service de Parasitologie Dakar" = 113001,
                     "Institut Pasteur ABIDJAN" = 114001,
                     "IPR/INSP, Bouak�"      = 114002,
                     "Komfo Anokye Teaching Hospital Kumasi, Ghana" = 115001,
                     "Malaria Research and Training Center Bamako" = 116001,
                     "CNFSR Maferenya - Guinee" = 117001,
                     "UW-Seattle"            = 201001,
                     "Abbvie"                = 201002,
                     "Celgene"               = 201003,
                     "EISAI-Andover"         = 201004,
                     "SRI-MenloPark"         = 201005,
                     "Abbott-AbbottPark"     = 201006,
                     "St-Jude"               = 201007,
                     "UoK"                   = 201008,
                     "OHSU"                  = 201009,
                     "Covance-SaltLakeCity"  = 201010,
                     "ACSA-Iquitos"          = 301001,
                     "TCGLS"                 = 501001,
                     "Zydus"                 = 501002,
                     "Aurigene"              = 501003,
                     "Wenlock District Hospital Mangalore" = 501004,
                     "EISAI-Tsukuba"         = 502001,
                     "Wuxi-City"             = 503001,
                     "DH-KhanhVinh"          = 504001,
                     "DH-HuongHoa"           = 504002,
                     "DH-BuDang"             = 504003,
                     "DH-PhuThien"           = 504004,
                     "Choray Hospital Ho Chi Minh City" = 504005,
                     "NIMPE Hanoi"           = 504006,
                     "Q2Sol-Singapore"       = 505001,
                     "Mahidol-Bangok"        = 506001,
                     "Mahidol-Tak"           = 506002,
                     "Clinic-MawkerThai"     = 506003,
                     "Hospital-MaeRamat"     = 506004,
                     "Hospital-Mae Sot"      = 506005,
                     "National Centre for Parasitology Phnom Penh" = 507001,
                     "Bethesda Hospital Tomohon" = 508001,
                     "RSUD Indonesia"        = 508002,
                     "Jayapura General Hospital Indonesia" = 508003,
                     "Ospital ng Palawan"    = 509001,
                     "Eulji General Hospital 280-1 Seoul" = 510001,
                     "Ilsan Paik Hospital"   = 510002,
                     "QIMR"                  = 601001,
                     "CDCO"                  = 601002,
                     "Nucleus"               = 601003,
                     "USC Moreton Bay"       = 601004,
                     "USC South Bank"        = 601005,
                     "SydPath Clinical Trials" = 601006,
                     "Sullivan Nicolaides Pathology" = 601007,
                     "StGeorge"              = 701001,
                     "Sequani"               = 701002,
                     "Richmond Pharmacology Ltd" = 701003,
                     "The Doctors Laboratory" = 701004,
                     "Quotient Sciences"     = 701005,
                     "GSK-TresCantos"        = 702001,
                     "TAD"                   = 702002,
                     "STPH"                  = 703001,
                     "SwissBQ"               = 703002,
                     "ITM-Tubingen"          = 704001,
                     "Sanofi-Aventis-Chilly" = 705001,
                     "SGS-Antwerpen"         = 706001,
                     "ITM-Antwerpen"         = 706002,
                     "ZNA-Antwerpen"         = 706003)
}


#' list_GenericType
#'
#' @description
#' @param GenericType
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
list_GenericType <- function(GenericType) {

  # Locate the file 'ActivityProperties.xlsx within MMVmalaria:
  ActivityProperties.Path <- file.path(get_MMVmalariaPath(subdir = "inst"),"activityProperties/ActivityProperties.xlsx")

  # Obtain the GenericType sheet:
  GenericType.Excel <- read_xlsx(ActivityProperties.Path, sheet = GenericType)

  # Generate GenericType list:
  GenericTypeList <- list()
  for (GenericType_k in GenericType.Excel[[GenericType]]){
    GenericTypeList[[GenericType_k]] <- GenericType_k
  }

  # Return list:
  return(GenericTypeList)
}

#' saveActivityInfo
#' Save Activity Information
#'
#' [saveActivityInfo] save a 'list' activity information into an R object
#' to the specified file 'ActivityInfoFile'.
#'
#' @param list A character vector containing the names of objects to be saved.
#' @param ActivityInfoFile The name of the file where the activity information will be saved.
#'
#' @md
#'
#' @export
#' @seealso [base::save], [saveMMV]
#' @family General Functions
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
saveActivityInfo <- function(list,
                             ActivityInfoFile = "../01-Data/S001-ActivityInfo"){

  # Start Date:
  if (file.exists(ActivityInfoFile)){
    StartDate <- substr(file.info(ActivityInfoFile)$ctime, 1, 10)
  } else{
    StartDate <- Sys.Date()
  }

  # Save Activity Info:
  saveMMV(list = list,
          file = ActivityInfoFile)
}
#' set_MMVsettings
#'
#' @description
#' @param IQRversion Default: NULL
#' @param ModifyRprofile Default: TRUE
#' @param UnloadMORpackage Default: TRUE
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Activity Monitoring
set_MMVsettings <- function(IQRversion       = NULL,
                            ModifyRprofile   = TRUE,
                            UnloadMORpackage = TRUE) {

  # Modify R Profile:
  if (ModifyRprofile){
    create_RprofileMMV(IQRversion)
  }

  # Unload MOR package RevoUtils & RevoUtilsMath:
  UnloadedPackage <- list("RevoUtils"     = FALSE,
                          "RevoUtilsMath" = FALSE)
  if (UnloadMORpackage){
    sessInfo <- sessionInfo()
    if ("RevoUtils" %in% c(sessInfo$basePkgs, names(sessInfo$otherPkgs))){
      detach("package:RevoUtils", unload=TRUE)
      UnloadedPackage[["RevoUtils"]] <- TRUE
    }
    if ("RevoUtilsMath" %in% c(sessInfo$basePkgs, names(sessInfo$otherPkgs))){
      detach("package:RevoUtilsMath", unload=TRUE)
      UnloadedPackage[["RevoUtilsMath"]] <- TRUE
    }
  }

  # Genereate output:
  settingsInfo <- list(UnloadedPackage = UnloadedPackage,
                       ModifyRprofile = ModifyRprofile)

  # Output:
  return(settingsInfo)
}

