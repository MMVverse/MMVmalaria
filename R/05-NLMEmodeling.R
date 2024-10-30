#' generate_InitialParametersIQRnlme
#'
#' @description
#' @param projectPath
#' @param IIVvalues0 Default: \code{NULL}
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family NLME Modeling
generate_InitialParametersIQRnlme <- function(projectPath,
                                              IIVvalues0 = NULL) {

  # Get modelSpec of the model:
  estText   <- readRDS(file.path(projectPath,"project.est"))
  modelSpec <- estText$modelSpec

  # Get the estimated POP value:
  IQRnlmeResult <- getResults_IQRnlmeProject(projectPath)
  POPvalues0    <- IQRnlmeResult$fixedEffects$values

  # Update POPvalues0 in modelspec:
  modelSpec$POPvalues0 <- POPvalues0

  # Manage IIVvalues0:
  if (!is.null(IIVvalues0)){
    if (is.numeric(IIVvalues0)){
      # Get parameters idx for which the IIV is gonna be estimated
      idx_Est <- which(as.numeric(IQRnlmeResult$randomEffects$estimated)==1)

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
      IIVvalues00          <- IQRnlmeResult$randomEffects$values
      names(IIVvalues00)   <- names(POPvalues0)
      modelSpec$IIVvalues0 <- IIVvalues00
    }
  }

  # Output:
  return(modelSpec)

}
#' generate_NLMEalgorithmSettings
#'
#' @description
#' @param multiTestN Default: 1
#' @param multiTestSD Default: 0.5
#' @param FLAGanalytic Default: \code{TRUE}
#' @param keepProjectFolder Default: \code{FALSE}
#' @param algOpt.SEED Default: 123456
#' @param algOpt.K1 Default: 500
#' @param algOpt.K2 Default: 200
#' @param algOpt.NRCHAINS Default: 5
#' @param algOpt.NONMEM.METHOD Default: \code{"SAEM"}
#' @param algOpt.NONMEM.MAXEVAL Default: 9999
#' @param algOpt.NONMEM.SIGDIGITS Default: 3
#' @param algOpt.NONMEM.PRINT Default: 1
#' @param algOpt.NONMEM.COVSTEP_MATRIX Default: \code{"S"}
#' @param algOpt.NONMEM.ADVAN7 Default: \code{TRUE}
#' @param algOpt.NONMEM.N1 Default: 1000
#' @param algOpt.NONMEM.TOL Default: 6
#' @param algOpt.NONMEM.SIGL Default: \code{NULL}
#' @param algOpt.NONMEM.M4 Default: \code{FALSE}
#' @param algOpt.NONMEM.FOCEIOFV Default: \code{FALSE}
#' @param algOpt.NONMEM.IMPORTANCESAMPLING Default: \code{TRUE}
#' @param algOpt.NONMEM.IMP_ITERATIONS Default: 10
#' @param algOpt.NONMEM.ITS Default: \code{TRUE}
#' @param algOpt.NONMEM.ITS_ITERATIONS Default: 10
#' @param algOpt.MONOLIX.individualParameters Default: \code{"conditionalMode"}
#' @param algOpt.MONOLIX.logLikelihood Default: \code{"Linearization"}
#' @param algOpt.MONOLIX.fim Default: \code{"Linearization"}
#' @param algOpt.MONOLIX.variability Default: \code{"Decreasing"}
#' @param algOpt.MONOLIX.startTime Default: \code{NULL}
#' @param algOpt.MONOLIX.STIFF Default: \code{TRUE}
#' @param algOpt.NLMIXR.method Default: \code{"SAEM"}
#' @param algOpt.NLMIXR.control Default: \code{NULL}
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family NLME Modeling
generate_NLMEalgorithmSettings <- function(# General setting:
                                           multiTestN        = 1    ,
                                           multiTestSD       = 0.5  ,
                                           FLAGanalytic      = TRUE ,
                                           keepProjectFolder = FALSE,

                                           # General algo. options:
                                           algOpt.SEED     = 123456,
                                           algOpt.K1       = 500   ,
                                           algOpt.K2       = 200   ,
                                           algOpt.NRCHAINS = 5     ,

                                           # NONMEM algo. options:
                                           algOpt.NONMEM.METHOD             = "SAEM",
                                           algOpt.NONMEM.MAXEVAL            = 9999  ,
                                           algOpt.NONMEM.SIGDIGITS          = 3     ,
                                           algOpt.NONMEM.PRINT              = 1     ,
                                           algOpt.NONMEM.COVSTEP_MATRIX     = "S"   ,
                                           algOpt.NONMEM.ADVAN7             = TRUE  ,
                                           algOpt.NONMEM.N1                 = 1000  ,
                                           algOpt.NONMEM.TOL                = 6     ,
                                           algOpt.NONMEM.SIGL               = NULL  ,
                                           algOpt.NONMEM.M4                 = FALSE ,
                                           algOpt.NONMEM.FOCEIOFV           = FALSE ,
                                           algOpt.NONMEM.IMPORTANCESAMPLING = TRUE  ,
                                           algOpt.NONMEM.IMP_ITERATIONS     = 10    ,
                                           algOpt.NONMEM.ITS                = TRUE  ,
                                           algOpt.NONMEM.ITS_ITERATIONS     = 10    ,
                                           algOpt.NONMEM.WRES               = NULL  ,
                                           algOpt.NONMEM.PRED               = NULL  ,
                                           algOpt.NONMEM.RES                = NULL  ,

                                           # Monolix algo. options:
                                           algOpt.MONOLIX.individualParameters = "conditionalMode",
                                           algOpt.MONOLIX.indivMCMClength      = 50               ,
                                           algOpt.MONOLIX.indivNsim            = 10               ,
                                           algOpt.MONOLIX.indivRatio           = 0.05             ,
                                           algOpt.MONOLIX.logLikelihood        = "Linearization"  ,
                                           algOpt.MONOLIX.fim                  = "Linearization"  ,
                                           algOpt.MONOLIX.variability          = "Decreasing"     ,
                                           algOpt.MONOLIX.startTime            = NULL             ,
                                           algOpt.MONOLIX.STIFF                = TRUE             ,

                                           # NLMIXR algo. options:
                                           algOpt.NLMIXR.method                      = "SAEM",
                                           algOpt.NLMIXR.control                     = NULL
                                           ) {

  # # Get IQRtools version
  # IQRversion     <- sessionInfo()$otherPkgs$IQRtools$Version

  # Create output:
  setting <- list()

  # Set "setting":
  setting$multiTestN        = multiTestN
  setting$multiTestSD       = multiTestSD
  setting$FLAGanalytic      = FLAGanalytic
  setting$keepProjectFolder = keepProjectFolder

  # General algo. options:
  setting$algOpt.SEED     = algOpt.SEED
  setting$algOpt.K1       = algOpt.K1
  setting$algOpt.K2       = algOpt.K2
  setting$algOpt.NRCHAINS = algOpt.NRCHAINS

  # NONMEM algo. options:
  setting$algOpt.NONMEM.METHOD             = algOpt.NONMEM.METHOD
  setting$algOpt.NONMEM.MAXEVAL            = algOpt.NONMEM.MAXEVAL
  setting$algOpt.NONMEM.SIGDIGITS          = algOpt.NONMEM.SIGDIGITS
  setting$algOpt.NONMEM.PRINT              = algOpt.NONMEM.PRINT
  setting$algOpt.NONMEM.COVSTEP_MATRIX     = algOpt.NONMEM.COVSTEP_MATRIX
  setting$algOpt.NONMEM.ADVAN7             = algOpt.NONMEM.ADVAN7
  setting$algOpt.NONMEM.N1                 = algOpt.NONMEM.N1
  setting$algOpt.NONMEM.TOL                = algOpt.NONMEM.TOL
  setting$algOpt.NONMEM.SIGL               = algOpt.NONMEM.SIGL
  setting$algOpt.NONMEM.M4                 = algOpt.NONMEM.M4
  setting$algOpt.NONMEM.FOCEIOFV           = algOpt.NONMEM.FOCEIOFV
  setting$algOpt.NONMEM.IMPORTANCESAMPLING = algOpt.NONMEM.IMPORTANCESAMPLING
  setting$algOpt.NONMEM.IMP_ITERATIONS     = algOpt.NONMEM.IMP_ITERATIONS
  setting$algOpt.NONMEM.ITS                = algOpt.NONMEM.ITS
  setting$algOpt.NONMEM.ITS_ITERATIONS     = algOpt.NONMEM.ITS_ITERATIONS
  setting$algOpt.NONMEM.WRES               = algOpt.NONMEM.WRES
  setting$algOpt.NONMEM.PRED               = algOpt.NONMEM.PRED
  setting$algOpt.NONMEM.RES                = algOpt.NONMEM.RES

  # Monolix algo. options:
  setting$algOpt.MONOLIX.individualParameters = algOpt.MONOLIX.individualParameters
  setting$algOpt.MONOLIX.indivMCMClength      = algOpt.MONOLIX.indivMCMClength
  setting$algOpt.MONOLIX.indivNsim            = algOpt.MONOLIX.indivNsim
  setting$algOpt.MONOLIX.indivRatio           = algOpt.MONOLIX.indivRatio
  setting$algOpt.MONOLIX.logLikelihood        = algOpt.MONOLIX.logLikelihood
  setting$algOpt.MONOLIX.fim                  = algOpt.MONOLIX.fim
  setting$algOpt.MONOLIX.variability          = algOpt.MONOLIX.variability
  setting$algOpt.MONOLIX.startTime            = algOpt.MONOLIX.startTime
  setting$algOpt.MONOLIX.STIFF                = algOpt.MONOLIX.STIFF

  # NLMIXR algo. options:
  setting$algOpt.NLMIXR.method  = algOpt.NLMIXR.method
  setting$algOpt.NLMIXR.control = algOpt.NLMIXR.control

  # Return Dataset:
  return(setting)
}
#' testPKmodelsIQR
#'
#' @description
#' @param testSetup
#' @param data
#' @param parOptions Default: list()
#' @param projectPath Default: 'PKmodels'
#' @param FLAGrun Default: \code{TRUE}
#' @param tool Default: 'MONOLIX'
#' @param toolVersion Default: \code{NULL}
#' @param ncores Default: 1
#' @param Nparallel Default: 1
#' @param NPROCESSORS Default: 1
#' @param NPROCESSORSpar Default: 1
#' @param setting Default: \code{NULL}
#' @return
#' @export
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Daniel Kaschek (IntiQuan), Mohammed H. Cherkaoui (MMV)
#' @family NLME Modeling
testPKmodelsIQR <- function(testSetup, data, parOptions = list(),
                            projectPath = "PKmodels", FLAGrun = TRUE,
                            tool = "MONOLIX", toolVersion = NULL,
                            ncores = 1, Nparallel = 1,
                            NPROCESSORS = 1, NPROCESSORSpar = 1, setting = NULL){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function Description ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Model ID:
  # The model ID is automatically generated by accounting for the number of compartment,
  # the elimination type, the absorption type, lag/no lag time, the error model, the
  # covariance model and the covariate model. The ID is generated acccording to:
  #    - Number of Compartment: 1, 2 or 3
  #    - Elimination type: 1->Linear, 2->Saturation or 3->Linear+Saturation
  #    - Absorption: 0->Zero Order, 1->First Order or 2->I.V.
  #    - Lag Time: 0->No Lag Time accounted for or 1->With Lag Time
  #    - Error Model: 0->Additive, 1->Proportional, 2-> Combined 1
  #    - Covariance Model: ID of the model (Between 1 and number of Covariance Model)
  #    - Covariate Model: ID of the model (Between 1 and number of Covariate Model)
  # e.g.:
  #  For example, Model_ID=2110111 means
  #  2-> 2 compartments
  #  1-> Linear Elimination
  #  1-> First order absorption
  #  0-> No lag Time
  #  1-> Prop error model
  #  1-> First covariance model
  #  1-> First covariate model


  # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # # Detect IQR version ----
  # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # IQRversion <- sessionInfo()$otherPkgs$IQRtools$Version


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Handle input arguments ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check testSetup: Model Assumptions to test
  if (!is.list(testSetup)) {
    stop("Input testSetup needs to be a list with fields 'Compartments', 'Absorption', 'Elimination', 'LagTime', and 'ErrorModels'")
  } else {
    if (!all(c('Compartments', 'Absorption', 'Elimination', 'LagTime', 'ErrorModels') %in% names(testSetup) ))
      stop("Input testSetup needs to be a list with fields 'Compartments', 'Absorption', 'Elimination', 'LagTime', and 'ErrorModels'")
  }

  # Clear project directory:
  if (dir.exists(projectPath)) {
    unlink(projectPath)
  }

  # Check subject number in dataset
  Nsubj <- length(unique(data$USUBJID))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define default setting (e.g algorithm's option) ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # General setting:
  if ("multiTestN"        %in% names(setting)){multiTestN        = setting$multiTestN       } else{multiTestN        = 1    }
  if ("multiTestSD"       %in% names(setting)){multiTestSD       = setting$multiTestSD      } else{multiTestSD       = 0.5  }
  if ("FLAGanalytic"      %in% names(setting)){FLAGanalytic      = setting$FLAGanalytic     } else{FLAGanalytic      = TRUE }
  if ("keepProjectFolder" %in% names(setting)){keepProjectFolder = setting$keepProjectFolder} else{keepProjectFolder = FALSE}

  # General algo. options:
  if ("algOpt.SEED"     %in% names(setting)){algOpt.SEED     = setting$algOpt.SEED    } else{algOpt.SEED     = 123456                    }
  if ("algOpt.K1"       %in% names(setting)){algOpt.K1       = setting$algOpt.K1      } else{algOpt.K1       = 500                       }
  if ("algOpt.K2"       %in% names(setting)){algOpt.K2       = setting$algOpt.K2      } else{algOpt.K2       = 200                       }
  if ("algOpt.NRCHAINS" %in% names(setting)){algOpt.NRCHAINS = setting$algOpt.NRCHAINS} else{algOpt.NRCHAINS = min(ceiling(50/Nsubj),10) }

  # NONMEM algo. options:
  if ("algOpt.NONMEM.METHOD"             %in% names(setting)){algOpt.NONMEM.METHOD             = setting$algOpt.NONMEM.METHOD            } else{algOpt.NONMEM.METHOD             = "SAEM"}
  if ("algOpt.NONMEM.MAXEVAL"            %in% names(setting)){algOpt.NONMEM.MAXEVAL            = setting$algOpt.NONMEM.MAXEVAL           } else{algOpt.NONMEM.MAXEVAL            = 9999  }
  if ("algOpt.NONMEM.SIGDIGITS"          %in% names(setting)){algOpt.NONMEM.SIGDIGITS          = setting$algOpt.NONMEM.SIGDIGITS         } else{algOpt.NONMEM.SIGDIGITS          = 3     }
  if ("algOpt.NONMEM.PRINT"              %in% names(setting)){algOpt.NONMEM.PRINT              = setting$algOpt.NONMEM.PRINT             } else{algOpt.NONMEM.PRINT              = 1     }
  if ("algOpt.NONMEM.COVSTEP_MATRIX"     %in% names(setting)){algOpt.NONMEM.COVSTEP_MATRIX     = setting$algOpt.NONMEM.COVSTEP_MATRIX    } else{algOpt.NONMEM.COVSTEP_MATRIX     = "S"   }
  if ("algOpt.NONMEM.ADVAN7"             %in% names(setting)){algOpt.NONMEM.ADVAN7             = setting$algOpt.NONMEM.ADVAN7            } else{algOpt.NONMEM.ADVAN7             = TRUE  }
  if ("algOpt.NONMEM.N1"                 %in% names(setting)){algOpt.NONMEM.N1                 = setting$algOpt.NONMEM.N1                } else{algOpt.NONMEM.N1                 = 1000  }
  if ("algOpt.NONMEM.TOL"                %in% names(setting)){algOpt.NONMEM.TOL                = setting$algOpt.NONMEM.TOL               } else{algOpt.NONMEM.TOL                = 6     }
  if ("algOpt.NONMEM.SIGL"               %in% names(setting)){algOpt.NONMEM.SIGL               = setting$algOpt.NONMEM.SIGL              } else{algOpt.NONMEM.SIGL               = NULL  }
  if ("algOpt.NONMEM.M4"                 %in% names(setting)){algOpt.NONMEM.M4                 = setting$algOpt.NONMEM.M4                } else{algOpt.NONMEM.M4                 = FALSE }
  if ("algOpt.NONMEM.FOCEIOFV"           %in% names(setting)){algOpt.NONMEM.FOCEIOFV           = setting$algOpt.NONMEM.FOCEIOFV          } else{algOpt.NONMEM.FOCEIOFV           = FALSE }
  if ("algOpt.NONMEM.IMPORTANCESAMPLING" %in% names(setting)){algOpt.NONMEM.IMPORTANCESAMPLING = setting$algOpt.NONMEM.IMPORTANCESAMPLING} else{algOpt.NONMEM.IMPORTANCESAMPLING = TRUE  }
  if ("algOpt.NONMEM.IMP_ITERATIONS"     %in% names(setting)){algOpt.NONMEM.IMP_ITERATIONS     = setting$algOpt.NONMEM.IMP_ITERATIONS    } else{algOpt.NONMEM.IMP_ITERATIONS     = 10    }
  if ("algOpt.NONMEM.ITS"                %in% names(setting)){algOpt.NONMEM.ITS                = setting$algOpt.NONMEM.ITS               } else{algOpt.NONMEM.ITS                = TRUE  }
  if ("algOpt.NONMEM.ITS_ITERATIONS"     %in% names(setting)){algOpt.NONMEM.ITS_ITERATIONS     = setting$algOpt.NONMEM.ITS_ITERATIONS    } else{algOpt.NONMEM.ITS_ITERATIONS     = 10    }
  if ("algOpt.NONMEM.WRES"               %in% names(setting)){algOpt.NONMEM.WRES               = setting$algOpt.NONMEM.WRES              } else{algOpt.NONMEM.WRES               = NULL  }
  if ("algOpt.NONMEM.PRED"               %in% names(setting)){algOpt.NONMEM.PRED               = setting$algOpt.NONMEM.PRED              } else{algOpt.NONMEM.PRED               = NULL  }
  if ("algOpt.NONMEM.RES"                %in% names(setting)){algOpt.NONMEM.RES                = setting$algOpt.NONMEM.RES               } else{algOpt.NONMEM.RES                = NULL  }

  # Monolix algo. options:
  if ("algOpt.MONOLIX.individualParameters" %in% names(setting)){algOpt.MONOLIX.individualParameters = setting$algOpt.MONOLIX.individualParameters} else{algOpt.MONOLIX.individualParameters = "conditionalMode"}
  if ("algOpt.MONOLIX.indivMCMClength"      %in% names(setting)){algOpt.MONOLIX.indivMCMClength      = setting$algOpt.MONOLIX.indivMCMClength     } else{algOpt.MONOLIX.indivMCMClength      = 50               }
  if ("algOpt.MONOLIX.indivNsim"            %in% names(setting)){algOpt.MONOLIX.indivNsim            = setting$algOpt.MONOLIX.indivNsim           } else{algOpt.MONOLIX.indivNsim            = 10               }
  if ("algOpt.MONOLIX.indivRatio"           %in% names(setting)){algOpt.MONOLIX.indivRatio           = setting$algOpt.MONOLIX.indivRatio          } else{algOpt.MONOLIX.indivRatio           = 0.05             }
  if ("algOpt.MONOLIX.logLikelihood"        %in% names(setting)){algOpt.MONOLIX.logLikelihood        = setting$algOpt.MONOLIX.logLikelihood       } else{algOpt.MONOLIX.logLikelihood        = "Linearization"  }
  if ("algOpt.MONOLIX.fim"                  %in% names(setting)){algOpt.MONOLIX.fim                  = setting$algOpt.MONOLIX.fim                 } else{algOpt.MONOLIX.fim                  = "Linearization"  }
  if ("algOpt.MONOLIX.variability"          %in% names(setting)){algOpt.MONOLIX.variability          = setting$algOpt.MONOLIX.variability         } else{algOpt.MONOLIX.variability          = "FirstStage"     }
  if ("algOpt.MONOLIX.startTime"            %in% names(setting)){algOpt.MONOLIX.startTime            = setting$algOpt.MONOLIX.startTime           } else{algOpt.MONOLIX.startTime            = NULL             }
  if ("algOpt.MONOLIX.STIFF"                %in% names(setting)){algOpt.MONOLIX.STIFF                = setting$algOpt.MONOLIX.STIFF               } else{algOpt.MONOLIX.STIFF                = TRUE             }

  # NLMIXR algo. options:
  if ("algOpt.NLMIXR.method"  %in% names(setting)){algOpt.NLMIXR.method  = setting$algOpt.NLMIXR.method } else{algOpt.NLMIXR.method  = "SAEM"}
  if ("algOpt.NLMIXR.control" %in% names(setting)){algOpt.NLMIXR.control = setting$algOpt.NLMIXR.control} else{algOpt.NLMIXR.control = NULL  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Establish models to test ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Get model files in library:
  modelLibraryPath <- file.path(get_MMVmalariaPath(subdir = "inst"), "modelLibrary/PKmodels")
  modelLibrary     <- IQRloadCSVdata(file.path(modelLibraryPath, "ModelLibrary.csv"))

  # Model selection based on setup options:
    # Filter compartments
  idxCOMP <- modelLibrary$Compartments %in% testSetup$Compartments
    # Filter absorption
  idxABS  <- modelLibrary$Absorption %in% testSetup$Absorption
    # Filter elimination
  idxELIM <- modelLibrary$Elimination %in% testSetup$Elimination
    # Intersection of all filters
  idxKEEP <- idxCOMP & idxABS & idxELIM

  # Check if the models are in the library:
  if (sum(idxKEEP) == 0){
    stop("No model in library found with required options.")
  }

  # List of Models Files to use:
  modelsTest <- paste0(modelLibrary$ModelFile[idxKEEP], ".txt")

  # Define List for Folder Name:
  Compartments <- modelLibrary$Compartments[idxKEEP]
  Elimination  <- ifelse(modelLibrary$Elimination[idxKEEP]=="linear"    , 1, ifelse(modelLibrary$Elimination[idxKEEP]=="saturable"  ,2, 3))
  Absorption   <- ifelse(modelLibrary$Absorption[idxKEEP] =="zero order", 0, ifelse(modelLibrary$Absorption[idxKEEP] =="first order",1, 2))
  Lag          <- ifelse(testSetup$LagTim                               , 1, 0)
  Error        <- ifelse(testSetup$ErrorModels=="abs"                   , 0, ifelse(testSetup$ErrorModels=="rel"  ,1, 2))


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define dosing ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dosingBOLUS <- list(  # to be used for IV and first order absorption
    INPUT1 = c(type="BOLUS")
  )
  dosingABS0  <- list(  # to be used for zero order absorption
    INPUT1 = c(type="ABSORPTION0", Tk0 = "Tk0")
  )


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Default settings for parameter specification ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  POPvalues0      <- c(Fabs0=1  , Fabs1= 1 , ka=10 , CL= 1 , Vc=10 , VMAX=10 , KM= 1 , Q1= 1 , Vp1=10 , Q2= 1  , Vp2=10 , Tlag1=0.5, Tk0=0.2 )
  POPestimate     <- c(Fabs0=0  , Fabs1= 0 , ka= 1 , CL= 1 , Vc= 1 , VMAX= 1 , KM= 1 , Q1= 1 , Vp1= 1 , Q2= 1  , Vp2= 1 , Tlag1= 1 , Tk0= 1  )
  IIVdistribution <- c(Fabs0='N', Fabs1='N', ka='L', CL='L', Vc='L', VMAX='L', KM='L', Q1='L', Vp1='L', Q2='L' , Vp2='L', Tlag1='L', Tk0='L' )
  IIVvalues0      <- c(Fabs0=0.5, Fabs1=0.5, ka=0.5, CL=0.5, Vc=0.5, VMAX=0.5, KM=0.5, Q1=0.5, Vp1=0.5, Q2= 0.5, Vp2=0.5, Tlag1=0.5, Tk0=0.5 )
  IIVestimate     <- c(Fabs0=0  , Fabs1= 0 , ka= 1 , CL= 1 , Vc= 1 , VMAX= 1 , KM= 1 , Q1= 1 , Vp1= 1 , Q2= 1  , Vp2= 1 , Tlag1= 1 , Tk0= 1  )


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Default settings for covariance model: none ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  covarianceModel <- list("Diagonal" = "diagonal")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Default settings for covariance model: none ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  covariateModel       <- list("NoCovariates" = NULL)
  covariateModelValues <- list("NoCovariates" = NULL)
  COVestimate          <- list("NoCovariates" = NULL)
  COVcentering         <- list("NoCovariates" = NULL)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Initialize error models ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  errorModels = list(
    abs    = list(
      OUTPUT1 = c("abs", c(abs0=1))
      ),
    rel    = list(
      OUTPUT1 = c("rel", c(rel0=0.3))
      ),
    absrel = list(
      OUTPUT1 = c("absrel", c(abs0=1,rel0=0.3))
      )
  )


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Adjust parameter specification based on parOptions ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (k in seq_along(parOptions$POPvalues0)){
    POPvalues0[[names(parOptions$POPvalues0[k])]] <- parOptions$POPvalues0[k]
  }
  for (k in seq_along(parOptions$POPestimate)){
    POPestimate[[names(parOptions$POPestimate[k])]] <- parOptions$POPestimate[k]
  }
  for (k in seq_along(parOptions$IIVdistribution)){
    IIVdistribution[[names(parOptions$IIVdistribution[k])]] <- parOptions$IIVdistribution[k]
  }
  for (k in seq_along(parOptions$IIVvalues0)){
    IIVvalues0[[names(parOptions$IIVvalues0[k])]] <- parOptions$IIVvalues0[k]
  }
  for (k in seq_along(parOptions$IIVestimate)){
    IIVestimate[[names(parOptions$IIVestimate[k])]] <- parOptions$IIVestimate[k]
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Adjust covariance specification based on parOptions ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ("covarianceModel" %in% names(parOptions)) {
    covarianceModel <- parOptions$covarianceModel
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Adjust covariate specification based on parOptions ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ("covariateModel" %in% names(parOptions)) {
    # Update covariateModel:
    covariateModel <- parOptions$covariateModel

    # Check covariateModelValues:
    if ("covariateModelValues" %in% names(parOptions)) {
      if (length(covariateModel) == length(parOptions$covariateModelValues)) {
        covariateModelValues <- parOptions$covariateModelValues
      } else {
        stop("Number of covariateModelValues sets needs to be same as covariateModel list.")
      }
    }

    # Check COVestimate:
    if ("COVestimate" %in% names(parOptions)) {
      if (length(covariateModel) == length(parOptions$COVestimate)) {
        COVestimate <- parOptions$COVestimate
      } else {
        stop("Number of COVestimate sets needs to be same as covariateModel list.")
      }
    }

    # Check COVcentering:
    if ("COVcentering" %in% names(parOptions)) {
      COVcentering <- parOptions$COVcentering
    }

  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Adjust error model specification based on parOptions ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ("abs0" %in% names(parOptions$errorValues0)) {
    errorModels$abs$OUTPUT1[["abs0"]] <- parOptions$errorValues0[["abs0"]]
    errorModels$absrel$OUTPUT1[["abs0"]] <- parOptions$errorValues0[["abs0"]]
  }
  if ("rel0" %in% names(parOptions$errorValues0)) {
    errorModels$rel$OUTPUT1[["rel0"]] <- parOptions$errorValues0[["rel0"]]
    errorModels$absrel$OUTPUT1[["rel0"]] <- parOptions$errorValues0[["rel0"]]
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Loop over models to test  and create projects ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ModelInfo <- data.frame(
    # Creation of model IDs. Mind the order of settings:
    # structure (Cpt, Eli, Abs), lagtime, error, covariance, covariate
    ModelID = suppressWarnings(
      levels(
        interaction(
          Compartments,
          Elimination,
          Absorption,
          Lag,
          Error,
          seq_along(covarianceModel),
          seq_along(covariateModel), sep = ""))
    ),
    Model = NA,
    LagTime = NA,
    ErrorModel = NA,
    CovarianceModel = NA,
    CovariateModel = NA
  )

  # Loop over structural models:
  for (kmod in 1:length(modelsTest)) {

    # Load the IQRmodel:
    model          <- IQRmodel(file.path(modelLibraryPath, modelsTest[kmod]))
    modelParameter <- grep("INPUT1",names(model$parameters), invert = TRUE, value = TRUE)

    # Select appropriate dosing:
    dosing <- NULL
    if (grepl("_iv|_abs1",modelsTest[kmod])){
      dosing <- dosingBOLUS
    }
    if (grepl("_abs0",modelsTest[kmod])){
      dosing <- dosingABS0
    }
    if (is.null(dosing)){
      stop("No dosing defined for ", modelsTest[kmod],".")
    }

    # Select Parameters:
    POPvalues0k      <- POPvalues0[modelParameter]
    POPestimatek     <- POPestimate[modelParameter]
    IIVdistributionk <- IIVdistribution[modelParameter]
    IIVvalues0k      <- IIVvalues0[modelParameter]
    IIVestimatek     <- IIVestimate[modelParameter]

    # Loop over LagTime to test:
    for (klt in 1:length(testSetup$LagTime)) {

      # Switch off estmation of lag time in case ...
      if (!testSetup$LagTime[klt]) {
        POPvalues0k[["Tlag1"]]      <-  0
        POPestimatek[["Tlag1"]]     <-  0
        IIVvalues0k[["Tlag1"]]      <-  0
        IIVestimatek[["Tlag1"]]     <-  0
        IIVdistributionk[["Tlag1"]] <- 'N'
      } else {
        POPvalues0k[["Tlag1"]]      <- POPvalues0[["Tlag1"]]
        POPestimatek[["Tlag1"]]     <- 1
        IIVvalues0k[["Tlag1"]]      <- IIVvalues0[["Tlag1"]]
        IIVestimatek[["Tlag1"]]     <- IIVestimate[["Tlag1"]]
        IIVdistributionk[["Tlag1"]] <- IIVdistribution[["Tlag1"]]
      }

      # Loop over error models to test:
      for (kerr in 1:length(testSetup$ErrorModels)) {

        # Loop over covariance models to test
        for (kcovn in 1:length(covarianceModel)) {

          # Loop over covariate models to test
          for (kcovtn in 1:length(covariateModel)) {
            # Define path for project:
            kCpt  <- Compartments[kmod]
            kEli  <- Elimination[kmod]
            kAbs  <- Absorption[kmod]
            kLag  <- Lag[klt]
            kErr  <- Error[kerr]
            IDk   <- paste0(kCpt, kEli, kAbs, kLag, kErr, kcovn, kcovtn)
            projk <- paste0("MODEL_", IDk)

            # Model identification to table
            modelIdx  <- (ModelInfo$ModelID==IDk)
            ModelInfo <- within(ModelInfo, {
              Model[modelIdx]           <- gsub(".txt","",modelsTest[kmod], fixed = TRUE)
              LagTime[modelIdx]         <- ifelse(POPestimatek[["Tlag1"]] == 0, "no lag time", "with lag time")
              ErrorModel[modelIdx]      <- testSetup$ErrorModels[[kerr]]
              CovarianceModel[modelIdx] <- names(covarianceModel)[kcovn]
              CovariateModel[modelIdx]  <- names(covariateModel)[kcovtn]
            })

            # Display message:
            cat("\n******\n   Creating project", projk, "\n")

            # Model specification
            modelSpec     <- list(
              # Parameters:
              POPvalues0      = POPvalues0k,
              POPestimate     = POPestimatek,
              IIVdistribution = IIVdistributionk,
              IIVvalues0      = IIVvalues0k,
              IIVestimate     = IIVestimatek,

              # Define error model and the initial guesses:
              errorModel      = errorModels[[testSetup$ErrorModels[[kerr]]]],

              # Define covariance model:
              covarianceModel = covarianceModel[[kcovn]],

              # Define covariate model:
              covariateModel       = covariateModel[[kcovtn]][intersect(modelParameter,names(covariateModel[[kcovtn]]))],
              covariateModelValues = covariateModelValues[[kcovtn]][intersect(modelParameter,names(covariateModelValues[[kcovtn]]))],
              COVestimate          = COVestimate[[kcovtn]][intersect(modelParameter,names(COVestimate[[kcovtn]]))],
              COVcentering         = COVcentering[[kcovtn]]
            )

            # Create IQRnlmeEst object:
            est <- IQRnlmeEst(model         = model,
                              dosing        = dosing,
                              data          = data,
                              modelSpec     = modelSpec)

            # Create NMLE project (monolix):
            proj <- IQRnlmeProject(est,projectPath=file.path(projectPath, projk),tool=tool,toolVersion = toolVersion,
                                   # General setting :
                                   multiTestN        = multiTestN,
                                   multiTestSD       = multiTestSD,
                                   FLAGanalytic      = FLAGanalytic,
                                   keepProjectFolder = keepProjectFolder,

                                   # General algo. options:
                                   algOpt.SEED     = algOpt.SEED,
                                   algOpt.K1       = algOpt.K1,
                                   algOpt.K2       = algOpt.K2,
                                   algOpt.NRCHAINS = algOpt.NRCHAINS,

                                   # NONMEM algo. options:
                                   algOpt.NONMEM.METHOD             = algOpt.NONMEM.METHOD,
                                   algOpt.NONMEM.MAXEVAL            = algOpt.NONMEM.MAXEVAL,
                                   algOpt.NONMEM.SIGDIGITS          = algOpt.NONMEM.SIGDIGITS,
                                   algOpt.NONMEM.PRINT              = algOpt.NONMEM.PRINT,
                                   algOpt.NONMEM.COVSTEP_MATRIX     = algOpt.NONMEM.COVSTEP_MATRIX,
                                   algOpt.NONMEM.ADVAN7             = algOpt.NONMEM.ADVAN7,
                                   algOpt.NONMEM.N1                 = algOpt.NONMEM.N1,
                                   algOpt.NONMEM.TOL                = algOpt.NONMEM.TOL,
                                   algOpt.NONMEM.SIGL               = algOpt.NONMEM.SIGL,
                                   algOpt.NONMEM.M4                 = algOpt.NONMEM.M4,
                                   algOpt.NONMEM.FOCEIOFV           = algOpt.NONMEM.FOCEIOFV,
                                   algOpt.NONMEM.IMPORTANCESAMPLING = algOpt.NONMEM.IMPORTANCESAMPLING,
                                   algOpt.NONMEM.IMP_ITERATIONS     = algOpt.NONMEM.IMP_ITERATIONS,
                                   algOpt.NONMEM.ITS                = algOpt.NONMEM.ITS,
                                   algOpt.NONMEM.ITS_ITERATIONS     = algOpt.NONMEM.ITS_ITERATIONS,
                                   algOpt.NONMEM.WRES               = algOpt.NONMEM.WRES,
                                   algOpt.NONMEM.PRED               = algOpt.NONMEM.PRED,
                                   algOpt.NONMEM.RES                = algOpt.NONMEM.RES,

                                   # Monolix algo. options:
                                   algOpt.MONOLIX.individualParameters = algOpt.MONOLIX.individualParameters,
                                   algOpt.MONOLIX.indivMCMClength      = algOpt.MONOLIX.indivMCMClength,
                                   algOpt.MONOLIX.indivNsim            = algOpt.MONOLIX.indivNsim,
                                   algOpt.MONOLIX.indivRatio           = algOpt.MONOLIX.indivRatio,
                                   algOpt.MONOLIX.logLikelihood        = algOpt.MONOLIX.logLikelihood,
                                   algOpt.MONOLIX.fim                  = algOpt.MONOLIX.fim,
                                   algOpt.MONOLIX.variability          = algOpt.MONOLIX.variability,
                                   algOpt.MONOLIX.startTime            = algOpt.MONOLIX.startTime,
                                   algOpt.MONOLIX.STIFF                = algOpt.MONOLIX.STIFF,

                                   # NLMIXR algo. options:
                                   algOpt.NLMIXR.method                = algOpt.NLMIXR.method,
                                   algOpt.NLMIXR.control               = algOpt.NLMIXR.control)
          } # End of loop covariate models
        } # End of loop covariance models
      } # End loop error models
    } # End loop lag times
  } # End loop structural models

  # Remove models from information that were not run
  ModelInfo <- subset(ModelInfo, !is.na(Model))
  # Re-order DataFrame:
  ModelInfo <- ModelInfo[order(ModelInfo$ModelID),]
  # Print Model information table
  IQRoutputTable(ModelInfo, xtitle = "Settings for tested models", filename = file.path(projectPath, "ModelInfo.txt"), report = TRUE)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Run projects----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (FLAGrun) {
    # Run estimation:
    cat("\n******\n   Running models in", projectPath, "\n")
    allproj <- as_IQRnlmeProjectMulti(projectPath, FLAGrecursive = TRUE)
    run_IQRnlmeProjectMulti(allproj, ncores=ncores, Nparallel=Nparallel)
    cat(" .. Finished!\n")

    # Summary tables:
    if (length(list.dirs(projectPath,recursive = FALSE))>1){
      summary(allproj, pathname=projectPath, FLAGreport=TRUE, FLAGremovePath=TRUE, order="BIC")
    }else{
      projectPath_1 <- list.dirs(projectPath,recursive = FALSE)[[1]]
      allproj_1    <- as_IQRnlmeProjectMulti(projectPath_1)
      summaryComments_IQRnlmeProjectMulti(allproj_1,
                                          order = "BIC",
                                          FLAGreport = TRUE,
                                          FLAGremovePath = TRUE,
                                          filename = file.path(projectPath,"01_model_overview"))
      summaryParameters_IQRnlmeProjectMulti(allproj_1,
                                          order = "BIC",
                                          FLAGreport = TRUE,
                                          FLAGremovePath = TRUE,
                                          filename = file.path(projectPath,"02_model_parameters"))
      summaryCorrelations_IQRnlmeProjectMulti(allproj_1,
                                          order = "BIC",
                                          FLAGreport = TRUE,
                                          FLAGremovePath = TRUE,
                                          filename = file.path(projectPath,"03_model_correlations"))
      summaryCovariates_IQRnlmeProjectMulti(allproj_1,
                                          order = "BIC",
                                          FLAGreport = TRUE,
                                          FLAGremovePath = TRUE,
                                          filename = file.path(projectPath,"04_model_covariates"))
    }

    # Add summary table in sub folder:
    if (multiTestN>1){
	  # Overall with Path:
	  summary(allproj, pathname=projectPath, FLAGreport=TRUE, FLAGremovePath=FALSE, order="BIC")

	  # Sub-Folders:
      projectPath_List <- list.dirs(path = projectPath, full.names = TRUE, recursive = FALSE)
      for (projectPath_k in projectPath_List){
        allproj_k <- as_IQRnlmeProjectMulti(projectPath_k, FLAGrecursive = TRUE)
        summary(allproj_k, pathname=projectPath_k, FLAGreport=TRUE, FLAGremovePath=TRUE, order="BIC")
      }
    }
  }
}


#' Get ModelSpec from a MMVmalaria Project
#'
#' Determines model specification based on parameter estimation results of a MMVmalaria
#' project stored in \code{projectPath}. Estimated parameter are defined as initial guess.
#'
#' @param projectPath A character string with the path to an \code{IQRsysProject} or \code{IQRnlmeProject} folder, or a GPF file.
#' @param fixParameters Flag whether to fix the estimated parameters
#'
#' @return Model specification list
#'
#' @export
#'
#' @author Mohammed H. Cherkaoui (MMV)
#' @family NLME Modeling
modelSpec_MMVmalariaProject <- function(projectPath, fixParameters = FALSE) {

  x <- getModelParameters_MMVmalariaProject(projectPath)

  # Determine estimation setting
  x$Parameters <- dplyr::mutate(x$Parameters,
                                POPestimate = ifelse(VALUE.RSE > 0, 1, 0),
                                IIVestimate = dplyr::case_when(
                                  IIV > 0 & IIV.RSE  > 0 ~ 1,
                                  IIV > 0 & IIV.RSE == 0 ~ 2,
                                  IIV == 0               ~ 0)
  )
  x$Beta <- dplyr::mutate(
    dplyr::rowwise(x$Beta),
    estimate  = ifelse(VALUE.RSE > 0, 1, 0),
    reference = ifelse(cov0, COV.REFERENCE, catN)
  )

  # Change estimation settings if parameters should be fixed
  if (fixParameters){
    x$Parameters$POPestimate <- 0
    x$Parameters$IIVestimate[x$Parameters$IIVestimate == 1] <- 2
    x$Beta$estimate <- 0
  }

  # Population values
  POPvalues0  <- structure(x$Parameters$VALUE,       names = x$Parameters$NAME)
  POPestimate <- structure(x$Parameters$POPestimate, names = x$Parameters$NAME)

  # Extract IIV values
  IIVvalues0      <- structure(x$Parameters$IIV,         names = x$Parameters$NAME)
  IIVdistribution <- structure(x$Parameters$DIST,        names = x$Parameters$NAME)
  IIVestimate     <- structure(x$Parameters$IIVestimate, names = x$Parameters$NAME)

  # Extract covariate model
  splitPar <- split(x$Beta, x$Beta$PARAMETER)
  covSpec <- lapply(splitPar, function(.p) {
    covariateModel       <- structure(list(.p$COVARIATE), names = unique(.p$PARAMETER))
    covariateModelValues <- structure(list(structure(.p$VALUE, names = .p$COVARIATE)),
                                      names = unique(.p$PARAMETER))
    COVestimate <- structure(list(structure(.p$estimate, names = .p$COVARIATE)),
                             names = unique(.p$PARAMETER))
    COVcentering <- structure(list(structure(.p$reference, names = .p$COVARIATE)),
                              names = unique(.p$PARAMETER))
    list(covariateModel=covariateModel,
         covariateModelValues = covariateModelValues,
         COVestimate = COVestimate,
         COVcentering = COVcentering)
  })
  covariateModel       <- purrr::flatten(purrr::map(covSpec, purrr::pluck("covariateModel")))
  covariateModelValues <- purrr::flatten(purrr::map(covSpec, purrr::pluck("covariateModelValues")))
  COVestimate          <- purrr::flatten(purrr::map(covSpec, purrr::pluck("COVestimate")))
  COVcentering         <- purrr::flatten(purrr::map(covSpec, purrr::pluck("COVcentering")))
  COVcentering <- purrr::flatten(COVcentering)
  COVcentering <- COVcentering[!duplicated(names(COVcentering))]

  # Covariance not yet implemented

  # Error model
  x$ResidualError <- dplyr::mutate(x$ResidualError,
                                   output = paste0("OUTPUT",gsub("[[:alpha:]]|[_]", "", x$ResidualError$NAME)),
                                   type   = dplyr::case_when(grepl("PROP", NAME) ~ "rel",
                                                             grepl("ADD", NAME) ~ "abs",
                                                             TRUE ~ NA_character_))
  splitErr <- split(x$ResidualError, x$ResidualError$output)
  errorModel <- lapply(splitErr, function(.e) {
    .e <- dplyr::arrange(.e, type)
    c(type = paste0(.e$type, collapse = ""), .e$VALUE)
  })

  out <- modelSpec_IQRest(
    POPvalues0 = POPvalues0,
    POPestimate = POPestimate,
    IIVdistribution = IIVdistribution,
    IIVvalues0 = IIVvalues0,
    IIVestimate = IIVestimate,
    covariateModel = covariateModel,
    covariateModelValues = covariateModelValues,
    COVestimate = COVestimate,
    COVcentering = COVcentering,
    errorModel = errorModel
  )
  out
}



#' Helper function to extract continuous covariate reference value
#' from FORMULA entry of MMVmalariaCSV model estimate table
#'
#' Looks basically for numerical value in string. If none or more than one
#' is found NA is return and a warung issued.
#'
#' @param formula Character string of the covariate equation used
#'
#' @return Reference covariate value (as character)
#' @export
#'
#' @examples
#' getCovRefValue("X=X_Pop*(WT0/65)^Beta")
getCovRefValue <- function(formula) {

  splitchar <- unlist(strsplit(formula, split = "[[:punct:]]"))
  value <- as.numeric(splitchar[grep("^[0-9 ]+$", splitchar)])

  if (length(value) != 1) {warning("Could not guess reference value."); value = NA}

  value
}
