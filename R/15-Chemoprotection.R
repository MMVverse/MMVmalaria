#' Generate Chemoprevention Model (Upwind First Order scheme)
#'
#' Generates a chemoprevention model script using a Upwind First Order scheme to characterize the liver stage
#'
#' @param nCpt Integer representing the number of liver compartments to be generated in the model
#' @param Finv Fraction of parasites invading the liver per infectious bite (Default: 1)
#' @param VBlood Blood volume (mL) (Default: 5000)
#' @param LiverDuration Duration of liver stage (days) (Default: 6)
#' @param LiverGrowth List of GR (referring to the PD model, e.g. exponential growth), parameters and extraStates if applicable
#' @param LiverActivity List of Kill (referring to the PD model, e.g. Emax equation), parameters and extraStates if applicable
#' @param BloodGrowth List of GR (referring to the PD model, e.g. exponential growth), parameters and extraStates if applicable
#' @param BloodActivity List of Kill (referring to the PD model, e.g. Emax equation), parameters and extraStates if applicable
#' @param PKmodel "IQRmodel" or list of states, parameters and variables
#' @param Output List of outputs to be included in the model script (Default: \code{"OUTPUT1 = "PBlood""})
#' @param outputFile Pathway where to save the resulting model file
#'
#' @return IQRmodel script to be used for estimation/simulation activities
#'
#' @examples
#' # Example 1: Emax model , PKmodel from library
#'  outputFile <- tempfile()
#'  generate_ChemoModel_UpwindFirstOrder(nCpt = 4,
#'                                       Finv = 1,
#'                                       VBlood = 5000, # in mL
#'                                       LiverDuration = 6,   # In Days,
#'                                       LiverGrowth   = list(GR          = "kGRliver",
#'                                                            parameters  = list(kGRliver = list(value = 0.07, notes = "Net parasite growth rate in liver (1/hr)")),
#'                                                            extraStates = NULL),
#'                                       LiverActivity = list(Kill = "EMAXLiver*(Cc+yps)^hLiver/((Cc+yps)^hLiver+EC50Liver^hLiver)",
#'                                                            parameters = list(EMAXLiver = list(value = 0.3, notes = "Maximum killing rate in liver (1/hr)"),
#'                                                                              EC50Liver = list(value = 1, notes = "Blood concentration achieving 50percent of maximum effect (ug/mL) in liver"),
#'                                                                              hLiver    = list(value = 3, notes = "Hill coefficient in liver")),
#'                                                            extraStates = NULL),
#'                                       BloodGrowth   = list(GR          = "kGRblood",
#'                                                            parameters  = list(kGRblood = list(value = 0.07, notes = "Net parasite growth rate in blood (1/hr)")),
#'                                                            extraStates = NULL),
#'                                       BloodActivity = list(Kill = "EMAXBlood*(Cc+yps)^hBlood/((Cc+yps)^hBlood+EC50Blood^hBlood)",
#'                                                            parameters = list(EMAXBlood = list(value = 0.3, notes = "Maximum killing rate in blodd (1/hr)"),
#'                                                                              EC50Blood = list(value = 1, notes = "Concentration achieving 50percent of maximum effect (ug/mL) in blood"),
#'                                                                              hBlood    = list(value = 3, notes = "Hill coefficient in blood")),
#'                                                            extraStates = NULL),
#'
#'                                       PKmodel = IQRmodel(file.path(get_MMVmalariaPath(subdir="inst"),"modelLibrary/PKmodels/model_3cpt_linear_abs1.txt")),
#'                                       Output = list(OUTPUT1 = "PBlood"),
#'                                       outputFile)
#'
#' # Example 2: Effect compartment model in blood, PKmodel manually defined
#'  outputFile <- tempfile()
#'  generate_ChemoModel_UpwindFirstOrder(nCpt = 4,
#'                                       Finv = 1,
#'                                       VBlood = 5000, # in mL
#'                                       LiverDuration = 6,   # In Days,
#'                                       LiverGrowth   = list(GR          = "kGRliver",
#'                                                            parameters  = list(kGRliver = list(value = 0.07, notes = "Net parasite growth rate in liver (1/hr)")),
#'                                                            extraStates = NULL),
#'                                       LiverActivity = list(Kill = "EMAXBlood*(Cc+yps)^hBlood/((Cc+yps)^hBlood+EC50Blood^hBlood)",
#'                                                            parameters = list(EMAXLiver = list(value = 0.3, notes = "Maximum killing rate in liver (1/hr)"),
#'                                                                              EC50Liver = list(value = 1, notes = "Blood concentration achieving 50percent of maximum effect (ug/mL) in liver"),
#'                                                                              hLiver    = list(value = 3, notes = "Hill coefficient in liver")),
#'                                                            extraStates = NULL),
#'                                       BloodGrowth   = list(GR          = "kGRblood",
#'                                                            parameters  = list(kGRblood = list(value = 0.07, notes = "Net parasite growth rate in blood (1/hr)")),
#'                                                            extraStates = NULL),
#'                                       BloodActivity = list(Kill = "EMAXLiver*(Ce+yps)^hLiver/((Ce+yps)^hLiver+EC50Liver^hLiver)",
#'                                                            parameters = list(EMAXBlood = list(value = 0.3, notes = "Maximum killing rate in blodd (1/hr)"),
#'                                                                              EC50Blood = list(value = 1, notes = "Concentration achieving 50percent of maximum effect (ug/mL) in blood"),
#'                                                                              hBlood    = list(value = 3, notes = "Hill coefficient in blood")),
#'                                                            extraStates = list(Ce = list(ODE = "d/dt(Ce) = ke*(Cc - Ce)", IC = "0"))),
#'
#'                                       PKmodel = list(states = list(Ad = list(ODE = "-ka*Ad + Fabs1*INPUT1", IC = "0"),
#'                                                                    Ac = list(ODE = "ka*Ad - Q1/Vc*Ac + Q1/Vp1*Ap1 - Q2/Vc*Ac + Q2/Vp2*Ap2 - CL/Vc*Ac", IC = "0"),
#'                                                                    Ap1 = list(ODE = "Q1/Vc*Ac - Q1/Vp1*Ap1", IC = "0"),
#'                                                                    Ap2 = list(ODE = "Q2/Vc*Ac - Q2/Vp2*Ap2", IC = "0")
#'                                                                    ),
#'                                                      parameters = list(Fabs1 = list(value = 1, notes="Relative bioavailability (-)"),
#'                                                                        ka    = list(value = 2, notes="Absorption rate parameter (1/hour)"),
#'                                                                        CL    = list(value = 3, notes="Apparent clearance (L/hour)"),
#'                                                                        Vc    = list(value = 32, notes="Apparent central volume (L)"),
#'                                                                        Q1    = list(value = 1, notes="Apparent intercompartmental clearance to first peripheral compartment (L/hour)"),
#'                                                                        Vp1   = list(value = 10, notes="Apparent first peripheral volume (L)"),
#'                                                                        Q2    = list(value = 1, notes="Apparent intercompartmental clearance to second peripheral compartment (L/hour)"),
#'                                                                        Vp2   = list(value = 10, notes="Apparent second peripheral volume (L)"),
#'                                                                        Tlag1 = list(value = 0, notes="Absorption lag time (hours)")
#'                                                                        ),
#'                                                     variables  = list(Cc = list(formula = "Ac/Vc"))),
#'                                       Output = list(OUTPUT1 = "PBlood"),
#'                                       outputFile)
#'
#' @export
#' @author Catalina Barcelo (MMV, \email{barceloc@@mmv.org})
#' @family Chemoprotection
generate_ChemoModel_UpwindFirstOrder <- function(nCpt,
                                                 Finv   = 1,
                                                 VBlood = 5000, # in mL
                                                 LiverDuration = 6,   # In Days

                                                 LiverGrowth   = list(GR          = "kGRliver",
                                                                      parameters  = list(kGRliver = list(value = 0.07, notes = "Net parasite growth rate in liver [1/hr]")),
                                                                      extraStates = NULL),


                                                 LiverActivity = list(Kill = "EMAXLiver*(Cc+yps)^hLiver/((Cc+yps)^hLiver+EC50Liver^hLiver)",
                                                                      parameters = list(EMAXLiver = list(value = 0.3, notes = "Maximum killing rate in liver [1/hr]"),
                                                                                        EC50Liver = list(value = 1, notes = "Blood concentration achieving 50percent of maximum effect [ug/mL] in liver"),
                                                                                        hLiver    = list(value = 3, notes = "Hill coefficient in liver [.]")),
                                                                      extraStates = NULL),


                                                 BloodGrowth   = list(GR          = "kGRblood",
                                                                      parameters  = list(kGRblood = list(value = 0.07, notes = "Net parasite growth rate in blood [1/hr]")),
                                                                      extraStates = NULL),


                                                 BloodActivity = list(Kill = "EMAXBlood*(Cc+yps)^hBlood/((Cc+yps)^hBlood+EC50Blood^hBlood)",
                                                                      parameters = list(EMAXBlood = list(value = 0.3, notes = "Maximum killing rate in blodd [1/hr]"),
                                                                                        EC50Blood = list(value = 1, notes = "Concentration achieving 50percent of maximum effect [ug/mL] in blood"),
                                                                                        hBlood    = list(value = 3, notes = "Hill coefficient in blood")),
                                                                      extraStates = NULL),

                                                 PKmodel = IQRmodel(file.path(get_MMVmalariaPath(subdir="inst"),"modelLibrary/PKmodels/model_3cpt_linear_abs1.txt")),

                                                 Output = list(OUTPUT1 = "PBlood"),
                                                 outputFile){

  #---------------------------------------------------#
  # STEP 0: Load IQRmodel if applicable ----
  #---------------------------------------------------#

  if(is.character(PKmodel)){
    PKmodel <- IQRmodel(PKmodel)
  }

  # Initialize INPUT count:
  k <- 1

  #---------------------------------------------------#
  # STEP 1: Initialize blocks ----
  #---------------------------------------------------#

  NameBlock <- c("********** MODEL NAME",
                 paste0("ChemoModel_UpwindFirstOrder_", nCpt, "Cpt"),
                 "",
                 "")

  NotesBlock <- c("********** MODEL NOTES",
                  "",
                  "")

  StatesBlock <- c("********** MODEL STATES",
                   "",
                   "")

  ParametersBlock <- c("********** MODEL PARAMETERS",
                       "",
                       "")

  VariablesBlock <- c("********** MODEL VARIABLES",
                      "",
                      "")

  ReactionsBlock <- c("********** MODEL REACTIONS",
                      "",
                      "")

  FunctionsBlock <- c("********** MODEL FUNCTIONS",
                      "",
                      "")

  EventsBlock <- c("********** MODEL EVENTS",
                   "",
                   "")


  #---------------------------------------------------#
  # STEP 2: Define Model States Block ----
  #---------------------------------------------------#

  #------------------------------------------#
  # 2.a: Differential Equations ----
  #------------------------------------------#

  StatesBlock <- c(StatesBlock,
                   "# Differential Equations:")

  # PK model equations:
  if(!is.null(PKmodel)){
    StatesBlock <- c(StatesBlock,
                     "",
                     "#	PK model:")
    for(states_k in names(PKmodel$states)){
      if(grepl("INPUT", PKmodel[["states"]][[states_k]][["ODE"]])){
        k <- k+1
      }
      StatesBlock <- c(StatesBlock,
                       paste0("d/dt(",states_k, ") = " , PKmodel[["states"]][[states_k]][["ODE"]]))
    }
  }

  # Liver stage equations:
  StatesBlock <- c(StatesBlock,
                   "",
                   "#	Liver Stage:",
                   paste0("d/dt(PLiver1) = (GRLiver - ktr - KillLiver)*PLiver1 + Finv*INPUT", k))
  if(nCpt > 1){
    for (i in 2:nCpt){
      StatesBlock_i <- paste0("d/dt(PLiver", i, ") =  ktr*PLiver", i-1, " + (GRLiver - ktr - KillLiver)*PLiver", i)
      StatesBlock   <- c(StatesBlock,
                         StatesBlock_i)
    }
  }

  # Liver growth stage extraStates equations:
  if(!is.null(LiverGrowth$extraStates)){
    for(states_k in names(LiverGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = " , LiverGrowth[["extraStates"]][[states_k]][["ODE"]]))
    }
  }

  # Liver activity stage extraStates equations:
  if(!is.null(LiverActivity$extraStates)){
    for(states_k in names(LiverActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = " , LiverActivity[["extraStates"]][[states_k]][["ODE"]]))
    }
  }

  # Blood stage equations:
  StatesBlock   <- c(StatesBlock,
                     "",
                     "#	Blood Stage: PBlood is expressed as parasitemia (p/mL)",
                     paste0("d/dt(PBlood) = ktr*PLiver", nCpt, "/VBlood + (GRBlood - KillBlood)*PBlood"))

  # blood growth stage extraStates equations:
  if(!is.null(BloodGrowth$extraStates)){
    for(states_k in names(BloodGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = " , BloodGrowth[["extraStates"]][[states_k]][["ODE"]]))
    }
  }

  # Blood activity stage extraStates equations:
  if(!is.null(BloodActivity$extraStates)){
    for(states_k in names(BloodActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = ", BloodActivity[["extraStates"]][[states_k]][["ODE"]]))
    }
  }


  #------------------------------------------#
  # 2.b: Initial Conditions ----
  #------------------------------------------#

  StatesBlock <- c(StatesBlock,
                   "",
                   "# Initial Conditions:",
                   "")
  # PK model
  if(!is.null(PKmodel)){
    StatesBlock <- c(StatesBlock,
                     "#	PK model:")
    for(states_k in names(PKmodel$states)){
      StatesBlock <- c(StatesBlock,
                       paste0(states_k, "(0) = " , PKmodel[["states"]][[states_k]][["IC"]]))
    }
  }

  # Liver Stage:
  StatesBlock <- c(StatesBlock,
                   "",
                   "#	Liver Stage:",
                   "PLiver1(0) = PLiverBase")
  if(nCpt > 1){
    for (i in 2:nCpt){
      StatesBlock_i <- paste0("PLiver", i, "(0) = 0")
      StatesBlock   <- c(StatesBlock,
                         StatesBlock_i)
    }
  }

  # Liver growth Stage extraStates:
  if(!is.null(LiverGrowth$extraStates)){
    for(states_k in names(LiverGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , LiverGrowth[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  # Liver activity Stage extraStates:
  if(!is.null(LiverActivity$extraStates)){
    for(states_k in names(LiverActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , LiverActivity[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  # Blood Equation:
  StatesBlock   <- c(StatesBlock,
                     "",
                     "#	Blood Stage:",
                     "PBlood(0) = 0")

  # Blood growth stage extraStates:
  if(!is.null(BloodGrowth$extraStates)){
    for(states_k in names(BloodGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , BloodGrowth[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  # Blood activity stage extraStates:
  if(!is.null(BloodActivity$extraStates)){
    for(states_k in names(BloodActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , BloodActivity[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  #---------------------------------------------------#
  # STEP 3: Define Model Parameters Block ----
  #---------------------------------------------------#

  # PK model
  if(!is.null(PKmodel)){
    ParametersBlock <- c(ParametersBlock,
                         "#	PK model:")
    for(para_k in names(PKmodel$parameters)){
      if(grepl("INPUT", para_k)){
        next
      }
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , PKmodel[["parameters"]][[para_k]][["value"]]))
    }
  }

  # Liver infection parameters:
  ParametersBlock <- c(ParametersBlock,
                       "",
                       "# Liver infection:",
                       "PLiverBase = 0      # Baseline parasite count",
                       paste0("LiverDuration = ", LiverDuration, "   # Duration of liver stage [days]"),
                       paste0("Finv = ", Finv, "            # Fraction of parasites invading the liver per infectious bite [.]"),
                       paste0("VBlood = ", VBlood, "       # Blood volume to transform parasite count to concentration [mL]"),
                       "",
                       "yps = 1e-9          # To avoid numerical crashes when computing PL")

  # Liver growth PD parameters:
  if(!is.null(LiverGrowth$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Liver growth stage")
    for(para_k in names(LiverGrowth$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , LiverGrowth[["parameters"]][[para_k]][["value"]], "     # ", LiverGrowth[["parameters"]][[para_k]][["notes"]]))

    }
  }

  # Liver activity PD parameters:
  if(!is.null(LiverActivity$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Liver activity stage")
    for(para_k in names(LiverActivity$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , LiverActivity[["parameters"]][[para_k]][["value"]], "     # ", LiverActivity[["parameters"]][[para_k]][["notes"]]))

    }
  }

  # Blood growth PD parameters:
  if(!is.null(BloodGrowth$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Blood growth stage")
    for(para_k in names(BloodGrowth$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , BloodGrowth[["parameters"]][[para_k]][["value"]], "     # ", BloodGrowth[["parameters"]][[para_k]][["notes"]]))

    }
  }

  # Blood activity PD parameters:
  if(!is.null(BloodActivity$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Blood activity stage")
    for(para_k in names(BloodActivity$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , BloodActivity[["parameters"]][[para_k]][["value"]], "     # ", BloodActivity[["parameters"]][[para_k]][["notes"]]))

    }
  }

  ParametersBlock <- c(ParametersBlock,
                       "")

  #---------------------------------------------------#
  # STEP 4: Define Model Variables Block ----
  #---------------------------------------------------#

  # PK model
  if(!is.null(PKmodel)){
    VariablesBlock <- c(VariablesBlock,
                        "#	PK model:")
    for(vara_k in names(PKmodel$variables)){
      if(grepl("OUTPUT", vara_k)){
        next
      }
      VariablesBlock <- c(VariablesBlock,
                          "",
                          paste0(vara_k, " = " , PKmodel[["variables"]][[vara_k]][["formula"]]))
    }
  }

  # Total parasite level in liver:
  VariablesBlock <- c(VariablesBlock,
                      "",
                      "# Liver compartments:",
                      paste0("nCpt =", nCpt, "     # Number of liver compartments in the model"),
                      paste0("ktr = nCpt/LiverDuration/24       # Transfer rate between liver compartments (1/hour) so the liver cycle takes LiverDuration days"),
                      "",
                      "# Total parasite level in liver:")
  PLiverTot <- "PLiver = PLiver1"
  if(nCpt > 1){
    for (i in 2:nCpt){
      PLiverTot <- paste0(PLiverTot, " + PLiver", i)
    }
  }

  # PK and Drug Effect:
  VariablesBlock <- c(VariablesBlock,
                      PLiverTot,
                      "",
                      "# Parasite growth in liver:",
                      paste0("GRLiver = ", LiverGrowth$GR),
                      "",
                      "# Drug effect in liver:",
                      paste0("KillLiver = ", LiverActivity$Kill),
                      "",
                      "# Parasite growth in blood:",
                      paste0("GRBlood = ", BloodGrowth$GR),
                      "",
                      "# Drug effect in blood:",
                      paste0("KillBlood = ", BloodActivity$Kill),
                      "",
                      "# -------- Output --------")


  for(output_k in names(Output)){
    VariablesBlock <- c(VariablesBlock,
                        paste0(output_k, " = " , Output[[output_k]]))

  }

  VariablesBlock <- c(VariablesBlock,
                      "")

  #---------------------------------------------------#
  # STEP 5: Merge blocks ----
  #---------------------------------------------------#

  # Merge all initialize blocks:
  modelText <- c(NameBlock, NotesBlock, StatesBlock,
                 ParametersBlock, VariablesBlock, ReactionsBlock,
                 FunctionsBlock, EventsBlock
  )

  # Save model text file:
  dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)
  writeLines(modelText, paste0(gsub(".txt", "", outputFile), "_", nCpt, ".txt"))

  # Return modelpath:
  return(paste0(gsub(".txt", "", outputFile), "_", nCpt, ".txt"))
}




#' Generate Chemoprevention Model (Flux Limiter scheme)
#'
#' Generates a chemoprevention model script using a Flux Limiter scheme to characterize the liver stage
#'
#' @param nCpt Integer representing the number of liver compartments to be generated in the model
#' @param Finv Fraction of parasites invading the liver per infectious bite (Default: 1)
#' @param VBlood Blood volume (mL) (Default: 5000)
#' @param LiverDuration Duration of liver stage (days) (Default: 6)
#' @param LiverGrowth List of GR (referring to the PD model, e.g. exponential growth), parameters and extraStates if applicable
#' @param LiverActivity List of Kill (referring to the PD model, e.g. Emax equation), parameters and extraStates if applicable
#' @param BloodGrowth List of GR (referring to the PD model, e.g. exponential growth), parameters and extraStates if applicable
#' @param BloodActivity List of Kill (referring to the PD model, e.g. Emax equation), parameters and extraStates if applicable
#' @param PKmodel "IQRmodel" or list of states, parameters and variables
#' @param Output List of outputs to be included in the model script (Default: \code{"OUTPUT1 = "PBlood""})
#' @param outputFile Pathway where to save the resulting model file
#'
#' @return IQRmodel script to be used for estimation/simulation activities
#'
#' @examples
#' # Example 1: Emax model , PKmodel from library
#'  outputFile <- tempfile()
#'  generate_ChemoModel_FluxLimiter(nCpt = 4,
#'                                  Finv = 1,
#'                                  VBlood = 5000, # in mL
#'                                  LiverDuration = 6,   # In Days,
#'                                  LiverGrowth   = list(GR          = "kGRliver",
#'                                                       parameters  = list(kGRliver = list(value = 0.07, notes = "Net parasite growth rate in liver (1/hr)")),
#'                                                       extraStates = NULL),
#'                                  LiverActivity = list(Kill = "EMAXLiver*(Cc+yps)^hLiver/((Cc+yps)^hLiver+EC50Liver^hLiver)",
#'                                                       parameters = list(EMAXLiver = list(value = 0.3, notes = "Maximum killing rate in liver (1/hr)"),
#'                                                                         EC50Liver = list(value = 1, notes = "Blood concentration achieving 50percent of maximum effect (ug/mL) in liver"),
#'                                                                         hLiver    = list(value = 3, notes = "Hill coefficient in liver")),
#'                                                      extraStates = NULL),
#'                                  BloodGrowth   = list(GR          = "kGRblood",
#'                                                       parameters  = list(kGRblood = list(value = 0.07, notes = "Net parasite growth rate in blood (1/hr)")),
#'                                                       extraStates = NULL),
#'                                  BloodActivity = list(Kill = "EMAXBlood*(Cc+yps)^hBlood/((Cc+yps)^hBlood+EC50Blood^hBlood)",
#'                                                       parameters = list(EMAXBlood = list(value = 0.3, notes = "Maximum killing rate in blodd (1/hr)"),
#'                                                                         EC50Blood = list(value = 1, notes = "Concentration achieving 50percent of maximum effect (ug/mL) in blood"),
#'                                                                         hBlood    = list(value = 3, notes = "Hill coefficient in blood")),
#'                                                       extraStates = NULL),
#'
#'                                  PKmodel = IQRmodel(file.path(get_MMVmalariaPath(subdir="inst"),"modelLibrary/PKmodels/model_3cpt_linear_abs1.txt")),
#'                                  Output = list(OUTPUT1 = "PBlood"),
#'                                  outputFile)
#'
#' # Example 2: Effect compartment model in blood, PKmodel manually defined
#'  outputFile <- tempfile()
#'  generate_ChemoModel_FluxLimiter(nCpt = 4,
#'                                  Finv = 1,
#'                                  VBlood = 5000, # in mL
#'                                  LiverDuration = 6,   # In Days,
#'                                  LiverGrowth   = list(GR          = "kGRliver",
#'                                                       parameters  = list(kGRliver = list(value = 0.07, notes = "Net parasite growth rate in liver (1/hr)")),
#'                                                       extraStates = NULL),
#'                                  LiverActivity = list(Kill = "EMAXBlood*(Cc+yps)^hBlood/((Cc+yps)^hBlood+EC50Blood^hBlood)",
#'                                                       parameters = list(EMAXLiver = list(value = 0.3, notes = "Maximum killing rate in liver (1/hr)"),
#'                                                                         EC50Liver = list(value = 1, notes = "Blood concentration achieving 50percent of maximum effect (ug/mL) in liver"),
#'                                                                         hLiver    = list(value = 3, notes = "Hill coefficient in liver")),
#'                                                       extraStates = NULL),
#'                                  BloodGrowth   = list(GR          = "kGRblood",
#'                                                       parameters  = list(kGRblood = list(value = 0.07, notes = "Net parasite growth rate in blood (1/hr)")),
#'                                                       extraStates = NULL),
#'                                  BloodActivity = list(Kill = "EMAXLiver*(Ce+yps)^hLiver/((Ce+yps)^hLiver+EC50Liver^hLiver)",
#'                                                       parameters = list(EMAXBlood = list(value = 0.3, notes = "Maximum killing rate in blodd (1/hr)"),
#'                                                                         EC50Blood = list(value = 1, notes = "Concentration achieving 50percent of maximum effect (ug/mL) in blood"),
#'                                                                         hBlood    = list(value = 3, notes = "Hill coefficient in blood")),
#'                                                       extraStates = list(Ce = list(ODE = "d/dt(Ce) = ke*(Cc - Ce)", IC = "0"))),
#'
#'                                  PKmodel = list(states = list(Ad = list(ODE = "-ka*Ad + Fabs1*INPUT1", IC = "0"),
#'                                                               Ac = list(ODE = "ka*Ad - Q1/Vc*Ac + Q1/Vp1*Ap1 - Q2/Vc*Ac + Q2/Vp2*Ap2 - CL/Vc*Ac", IC = "0"),
#'                                                               Ap1 = list(ODE = "Q1/Vc*Ac - Q1/Vp1*Ap1", IC = "0"),
#'                                                               Ap2 = list(ODE = "Q2/Vc*Ac - Q2/Vp2*Ap2", IC = "0")
#'                                                               ),
#'                                                 parameters = list(Fabs1 = list(value = 1, notes="Relative bioavailability (-)"),
#'                                                                   ka    = list(value = 2, notes="Absorption rate parameter (1/hour)"),
#'                                                                   CL    = list(value = 3, notes="Apparent clearance (L/hour)"),
#'                                                                   Vc    = list(value = 32, notes="Apparent central volume (L)"),
#'                                                                   Q1    = list(value = 1, notes="Apparent intercompartmental clearance to first peripheral compartment (L/hour)"),
#'                                                                   Vp1   = list(value = 10, notes="Apparent first peripheral volume (L)"),
#'                                                                   Q2    = list(value = 1, notes="Apparent intercompartmental clearance to second peripheral compartment (L/hour)"),
#'                                                                   Vp2   = list(value = 10, notes="Apparent second peripheral volume (L)"),
#'                                                                   Tlag1 = list(value = 0, notes="Absorption lag time (hours)")
#'                                                                   ),
#'                                                variables  = list(Cc = list(formula = "Ac/Vc"))),
#'                                  Output = list(OUTPUT1 = "PBlood"),
#'                                  outputFile)
#'
#' @export
#' @author Catalina Barcelo (MMV, \email{barceloc@@mmv.org})
#' @family Chemoprotection
generate_ChemoModel_FluxLimiter <- function(nCpt,
                                            Finv   = 1,
                                            VBlood = 5000, # in mL
                                            LiverDuration = 6,   # In Days

                                            LiverGrowth   = list(GR          = "kGRliver",
                                                                 parameters  = list(kGRliver = list(value = 0.07, notes = "Net parasite growth rate in liver [1/hr]")),
                                                                 extraStates = NULL),


                                            LiverActivity = list(Kill = "EMAXLiver*(Cc+yps)^hLiver/((Cc+yps)^hLiver+EC50Liver^hLiver)",
                                                                 parameters = list(EMAXLiver = list(value = 0.3, notes = "Maximum killing rate in liver [1/hr]"),
                                                                                   EC50Liver = list(value = 1, notes = "Blood concentration achieving 50percent of maximum effect [ug/mL] in liver"),
                                                                                   hLiver    = list(value = 3, notes = "Hill coefficient in liver")),
                                                                 extraStates = NULL),


                                            BloodGrowth   = list(GR          = "kGRblood",
                                                                 parameters  = list(kGRblood = list(value = 0.07, notes = "Net parasite growth rate in blood [1/hr]")),
                                                                 extraStates = NULL),


                                            BloodActivity = list(Kill = "EMAXBlood*(Cc+yps)^hBlood/((Cc+yps)^hBlood+EC50Blood^hBlood)",
                                                                 parameters = list(EMAXBlood = list(value = 0.3, notes = "Maximum killing rate in blodd [1/hr]"),
                                                                                   EC50Blood = list(value = 1, notes = "Concentration achieving 50percent of maximum effect [ug/mL] in blood"),
                                                                                   hBlood    = list(value = 3, notes = "Hill coefficient in blood")),
                                                                 extraStates = NULL),

                                            PKmodel = IQRmodel(file.path(get_MMVmalariaPath(subdir="inst"),"modelLibrary/PKmodels/model_3cpt_linear_abs1.txt")),

                                            Output = list(OUTPUT1 = "PBlood"),
                                            outputFile){

  #---------------------------------------------------#
  # STEP 0: Load IQRmodel if applicable ----
  #---------------------------------------------------#

  if(is.character(PKmodel)){
    PKmodel <- IQRmodel(PKmodel)
  }

  # Initialize INPUT count:
  k <- 1

  #---------------------------------------------------#
  # STEP 1: Initialize blocks ----
  #---------------------------------------------------#

  NameBlock <- c("********** MODEL NAME",
                 paste0("ChemoModel_FluxLimiter_", nCpt, "Cpt"),
                 "",
                 "")

  NotesBlock <- c("********** MODEL NOTES",
                  "",
                  "")

  StatesBlock <- c("********** MODEL STATES",
                   "",
                   "")

  ParametersBlock <- c("********** MODEL PARAMETERS",
                       "",
                       "")

  VariablesBlock <- c("********** MODEL VARIABLES",
                      "",
                      "")

  ReactionsBlock <- c("********** MODEL REACTIONS",
                      "",
                      "")

  FunctionsBlock <- c("********** MODEL FUNCTIONS",
                      "",
                      "")

  EventsBlock <- c("********** MODEL EVENTS",
                   "",
                   "")


  #---------------------------------------------------#
  # STEP 2: Define Model States Block ----
  #---------------------------------------------------#

  #------------------------------------------#
  # 2.a: Differential Equations ----
  #------------------------------------------#

  StatesBlock <- c(StatesBlock,
                   "# Differential Equations:")

  # PK model equations:
  if(!is.null(PKmodel)){
    StatesBlock <- c(StatesBlock,
                     "",
                     "#	PK model:")
    for(states_k in names(PKmodel$states)){
      if(grepl("INPUT", PKmodel[["states"]][[states_k]][["ODE"]])){
        k <- k+1
      }
      StatesBlock <- c(StatesBlock,
                       paste0("d/dt(",states_k, ") = " , PKmodel[["states"]][[states_k]][["ODE"]]))
    }
  }

  # Liver stage equations:
  # First Liver Stage Compartment:
  #   OSPRE Flux Limiter
  Phi_1 <- "1.5*(PLiver3-PLiver2)*(PLiver3-PLiver1)/((PLiver3-PLiver2)*(PLiver3-PLiver1) + (PLiver2-PLiver1)^2 + yps^2)*1/(1+exp(-(PLiver2-PLiver1)*(PLiver3-PLiver2)/yps^3))"
  #   PLiver1 & PLiver2
  StatesBlock <- c(StatesBlock,
                   "",
                   "#	Liver Stage:",
                   paste0("d/dt(PLiver1) = (GRLiver - ktr - KillLiver)*PLiver1 + Finv*INPUT", k),
                   paste0("d/dt(PLiver2) = ktr*PLiver1 - ktr*(PLiver2 + ", Phi_1, "*(PLiver2-PLiver1)/2) + (GRLiver - KillLiver)*PLiver2"))

  # Loop over liver Stage Compartment:
  for (i in 3:(nCpt-1)){

    # OSPRE Flux Limiter:
    Phi_0 <- paste0("1.5*(PLiver", i  , "-PLiver", i-1, ")*(PLiver", i  , "-PLiver", i-2, ")/((PLiver", i  , "-PLiver", i-1, ")*(PLiver", i  , "-PLiver", i-2, ") + (PLiver", i-1, "-PLiver", i-2, ")^2 + yps^2)*1/(1+exp(-(PLiver", i-1, "-PLiver", i-2, ")*(PLiver", i  , "-PLiver", i-1, ")/yps^3))")
    Phi_1 <- paste0("1.5*(PLiver", i+1, "-PLiver", i  , ")*(PLiver", i+1, "-PLiver", i-1, ")/((PLiver", i+1, "-PLiver", i  , ")*(PLiver", i+1, "-PLiver", i-1, ") + (PLiver", i  , "-PLiver", i-1, ")^2 + yps^2)*1/(1+exp(-(PLiver", i  , "-PLiver", i-1, ")*(PLiver", i+1, "-PLiver", i  , ")/yps^3))")
    # PLiver_i:
    StatesBlock_i <- paste0("d/dt(PLiver", i, ") =",
                            " ktr*(PLiver", i-1," + ", Phi_0, "*(PLiver", i-1, "-PLiver", i-2, ")/2)",
                            " - ktr*(PLiver", i  ," + ", Phi_1, "*(PLiver", i  , "-PLiver", i-1, ")/2)",
                            " + (GRLiver - KillLiver)*PLiver", i)

    # Concatenate Blocks:
    StatesBlock   <- c(StatesBlock,
                       StatesBlock_i)
  }

  # Last Liver Stage Compartment:
  #   OSPRE Flux Limiter
  Phi_0 <- paste0("1.5*(PLiver", nCpt, "-PLiver", nCpt-1, ")*(PLiver", nCpt, "-PLiver", nCpt-2, ")/((PLiver", nCpt, "-PLiver", nCpt-1, ")*(PLiver", nCpt, "-PLiver", nCpt-2, ") + (PLiver", nCpt-1, "-PLiver", nCpt-2, ")^2 + yps^2)*1/(1+exp(-(PLiver", nCpt-1, "-PLiver", nCpt-2, ")*(PLiver", nCpt, "-PLiver", nCpt-1, ")/yps^3))")
  #   PLiver_nCpt
  StatesBlock   <- c(StatesBlock,
                     paste0("d/dt(PLiver", nCpt, ") =",
                            " ktr*(PLiver", nCpt-1 ," + ", Phi_0, "*(PLiver", nCpt-1, "-PLiver", nCpt-2, ")/2)",
                            " - ktr*PLiver", nCpt,
                            " + (GRLiver - KillLiver)*PLiver", nCpt))

  # Liver growth stage extraStates equations:
  if(!is.null(LiverGrowth$extraStates)){
    for(states_k in names(LiverGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = " , LiverGrowth[["extraStates"]][[states_k]][["ODE"]]))
    }
  }


  # Liver activity stage extraStates equations:
  if(!is.null(LiverActivity$extraStates)){
    for(states_k in names(LiverActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = " , LiverActivity[["extraStates"]][[states_k]][["ODE"]]))
    }
  }
  # Blood stage equations:
  StatesBlock   <- c(StatesBlock,
                     "",
                     "#	Blood Stage: PBlood is expressed as parasitemia (p/mL)",
                     paste0("d/dt(PBlood) = ktr*PLiver", nCpt, "/VBlood + (GRBlood - KillBlood)*PBlood"))

  # Blood growth stage extraStates equations:
  if(!is.null(BloodGrowth$extraStates)){
    for(states_k in names(BloodGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = ", BloodGrowth[["extraStates"]][[states_k]][["ODE"]]))
    }
  }

  # Blood activity stage extraStates equations:
  if(!is.null(BloodActivity$extraStates)){
    for(states_k in names(BloodActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0("d/dt(", states_k, ") = ", BloodActivity[["extraStates"]][[states_k]][["ODE"]]))
    }
  }


  #------------------------------------------#
  # 2.b: Initial Conditions ----
  #------------------------------------------#

  StatesBlock <- c(StatesBlock,
                   "",
                   "# Initial Conditions:",
                   "")
  # PK model
  if(!is.null(PKmodel)){
    StatesBlock <- c(StatesBlock,
                     "#	PK model:")
    for(states_k in names(PKmodel$states)){
      StatesBlock <- c(StatesBlock,
                       paste0(states_k, "(0) = " , PKmodel[["states"]][[states_k]][["IC"]]))
    }
  }

  # Liver Stage:
  StatesBlock <- c(StatesBlock,
                   "",
                   "#	Liver Stage:",
                   "PLiver1(0) = PLiverBase")
  if(nCpt > 1){
    for (i in 2:(nCpt)){
      StatesBlock_i <- paste0("PLiver", i, "(0) = 0")
      StatesBlock   <- c(StatesBlock,
                         StatesBlock_i)
    }
  }

  # Liver growth stage extraStates:
  if(!is.null(LiverGrowth$extraStates)){
    for(states_k in names(LiverGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , LiverGrowth[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  # Liver activity stage extraStates:
  if(!is.null(LiverActivity$extraStates)){
    for(states_k in names(LiverActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , LiverActivity[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  # Blood Equation:
  StatesBlock   <- c(StatesBlock,
                     "",
                     "#	Blood Stage:",
                     "PBlood(0) = 0")

  # Blood grwoth stage extraStates:
  if(!is.null(BloodGrowth$extraStates)){
    for(states_k in names(BloodGrowth$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , BloodGrowth[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  # Blood activity stage extraStates:
  if(!is.null(BloodActivity$extraStates)){
    for(states_k in names(BloodActivity$extraStates)){
      StatesBlock <- c(StatesBlock,
                       "",
                       paste0(states_k, "(0) = " , BloodActivity[["extraStates"]][[states_k]][["IC"]]))
    }
  }

  #---------------------------------------------------#
  # STEP 3: Define Model Parameters Block ----
  #---------------------------------------------------#

  # PK model
  if(!is.null(PKmodel)){
    ParametersBlock <- c(ParametersBlock,
                         "#	PK model:")
    for(para_k in names(PKmodel$parameters)){
      if(grepl("INPUT", para_k)){
        next
      }
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , PKmodel[["parameters"]][[para_k]][["value"]]))
    }
  }

  # Liver infection parameters:
  ParametersBlock <- c(ParametersBlock,
                       "",
                       "# Liver infection:",
                       "PLiverBase = 0      # Baseline parasite count",
                       paste0("LiverDuration = ", LiverDuration, "   # Duration of liver stage [days]"),
                       paste0("Finv = ", Finv, "            # Fraction of parasites invading the liver per infectious bite [.]"),
                       paste0("VBlood = ", VBlood, "       # Blood volume to transform parasite count to concentration [mL]"),
                       "",
                       "yps = 1e-9          # To avoid numerical crashes when computing PL")

  # Liver growth PD parameters:
  if(!is.null(LiverGrowth$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Liver growth stage")
    for(para_k in names(LiverGrowth$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , LiverGrowth[["parameters"]][[para_k]][["value"]], "     # ", LiverGrowth[["parameters"]][[para_k]][["notes"]]))

    }
  }

  # Liver activity PD parameters:
  if(!is.null(LiverActivity$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Liver stage activity")
    for(para_k in names(LiverActivity$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , LiverActivity[["parameters"]][[para_k]][["value"]], "     # ", LiverActivity[["parameters"]][[para_k]][["notes"]]))

    }
  }

  # Blood growth PD parameters:
  if(!is.null(BloodGrowth$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Blood activity stage")
    for(para_k in names(BloodGrowth$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , BloodGrowth[["parameters"]][[para_k]][["value"]], "     # ", BloodGrowth[["parameters"]][[para_k]][["notes"]]))

    }
  }

  # Blood activity PD parameters:
  if(!is.null(BloodActivity$parameters)){
    ParametersBlock <- c(ParametersBlock,
                         "",
                         "# Blood activity stage")
    for(para_k in names(BloodActivity$parameters)){
      ParametersBlock <- c(ParametersBlock,
                           paste0(para_k, " = " , BloodActivity[["parameters"]][[para_k]][["value"]], "     # ", BloodActivity[["parameters"]][[para_k]][["notes"]]))

    }
  }

  ParametersBlock <- c(ParametersBlock,
                       "")

  #---------------------------------------------------#
  # STEP 4: Define Model Variables Block ----
  #---------------------------------------------------#

  # PK model
  if(!is.null(PKmodel)){
    VariablesBlock <- c(VariablesBlock,
                        "#	PK model:")
    for(vara_k in names(PKmodel$variables)){
      if(grepl("OUTPUT", vara_k)){
        next
      }
      VariablesBlock <- c(VariablesBlock,
                          "",
                          paste0(vara_k, " = " , PKmodel[["variables"]][[vara_k]][["formula"]]))
    }
  }

  # Total parasite level in liver:
  VariablesBlock <- c(VariablesBlock,
                      "",
                      "# Liver compartments:",
                      paste0("nCpt =", nCpt, "     # Number of liver compartments in the model"),
                      paste0("ktr = nCpt/(LiverDuration*24)       # Transfer rate between liver compartments (1/hour) so the liver cycle takes LiverDuration days"),
                      "",
                      "# Total parasite level in liver:")
  PLiverTot <- "PLiver = PLiver1"
  if(nCpt > 1){
    for (i in 2:nCpt){
      PLiverTot <- paste0(PLiverTot, " + PLiver", i)
    }
  }

  # PK and Drug Effect:
  VariablesBlock <- c(VariablesBlock,
                      PLiverTot,
                      "",
                      "# Parasite growth in liver:",
                      paste0("GRLiver = ", LiverGrowth$GR),
                      "",
                      "# Drug effect in liver:",
                      paste0("KillLiver = ", LiverActivity$Kill),
                      "",
                      "# Parasite growth in blood:",
                      paste0("GRBlood = ", BloodGrowth$GR),
                      "",
                      "# Drug effect in blood:",
                      paste0("KillBlood = ", BloodActivity$Kill),
                      "",
                      "# -------- Output --------")


  for(output_k in names(Output)){
    VariablesBlock <- c(VariablesBlock,
                        paste0(output_k, " = " , Output[[output_k]]))

  }

  VariablesBlock <- c(VariablesBlock,
                      "")

  #---------------------------------------------------#
  # STEP 5: Merge blocks ----
  #---------------------------------------------------#

  # Merge all initialize blocks:
  modelText <- c(NameBlock, NotesBlock, StatesBlock,
                 ParametersBlock, VariablesBlock, ReactionsBlock,
                 FunctionsBlock, EventsBlock
  )

  # Save model text file:
  dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)
  writeLines(modelText, paste0(gsub(".txt", "", outputFile), "_", nCpt, ".txt"))

  # Return modelpath:
  return(paste0(gsub(".txt", "", outputFile), "_", nCpt, ".txt"))
}


#' Evaluate breakthrough event,
#' given an individual time-course of parasitemia measurements and LLOQ parasitemia.
#'
#' Evaluate whether a breakthrough event has been observed and estimate the time of first occurrence.
#' If no event has been observed, either because it did not occur, or if the
#' measurement data for the subject is missing, the time of the last measurement
#' is reported to indicate that the subject is not anymore in the population after
#' that moment in time.
#'
#' A breakthrough event is observed if PL exceeds LLOQ. We report time as the first time above LLOQ.
#' Trivially, a breakthrough event IS NOT observed if PL did not exceed LLOQ. We report time as the
#' time of last measurement.
#'
#' @param dataSim data.table produced by `simulate_virtualtrials()` or `sim_IQRmodel()`, specifically `$simPKPD`.
#' @param LLOQ LLOQ parasitemia in the same scale as in `paraCOL` (i.e. linear or ln - numeric).
#' @param timeCOL column name of `dataSim` containing time records (Default: `TIME`).
#' @param paraCOL column name of `dataSim` containing parasitemia records (Default: `PL`).
#' @param Plog Indicate if parasitemia is in the log or linear scale (Default: `TRUE` which means it is logged).
#' @param FLAGinterpolateTime A logical indicating if the PL measurements should be interpolated within the simulation period.
#' @details
#'
#' @return A data.frame with columns as follows:
#' \itemize{
#' \item USUBJID
#' \item TrialID
#' \item NAME: name of the event evaluated "Breakthrough"
#' \item TIME: Time to first exceed above LLOQ
#' \item VALUE: TRUE or FALSE statement regarding breakthrough event
#' }
#' @examples
#'
#' @export
#' @author Catalina Barcelo (MMV, \email{barceloc@@mmv.org})
#' @family Chemoprotection
evaluate_BreakthroughEvent <- function(dataSim,
                                       LLOQ,
                                       timeCOL       = "TIME",
                                       paraCOL       = "PL",
                                       Plog          = TRUE,
                                       FLAGinterpolateTime = FALSE) {

  #--------------------------------------#
  # STEP 1: Checks ----
  #--------------------------------------#

  # Make sure we have the right column:
  if(!(timeCOL %in% names(dataSim)) || !(paraCOL %in% names(dataSim))){
    stop("'", timeCOL, "', '", paraCOL, " 'are not columns available in 'dataSim': Please Adjust.")
  }

  # Convert to data.table and change column names:
  data.table::setDT(dataSim)
  data.table::setnames(dataSim,
                       c(timeCOL, paraCOL),
                       c("TIME" , "PL"   ),
                       skip_absent = TRUE)
  dataSim <- dataSim[,list(TIME, PL)]

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


  #--------------------------------------#
  # STEP 2: Cleaning ----
  #--------------------------------------#

  # Set '<LLOQ'-infinite parasitemia to LLOQ - 0.0001
  dataSim$PL[dataSim$PL < LLOQ & is.infinite(dataSim$PL)] <- LLOQ - 0.0001

  # Take average at ties in TIME, i.e. time-points for which more than one measurement has been done.
  dataSim[, list(PL = mean(PL, na.rm = TRUE)), keyby = TIME]


  #--------------------------------------#
  # STEP 3: Interpolation parasitemia at point of interest ----
  #--------------------------------------#

  if (length(dataSim$TIME)>=2 & FLAGinterpolateTime) {
    # interpolation at these time-points.
    interpTimes <- seq(min(dataSim$TIME), max(dataSim$TIME), by = 1)

    if(length(interpTimes) > 0L) {
      interpResult <- approx(
        x = dataSim$TIME, y = dataSim$PL, xout = interpTimes, ties = list("ordered", mean), rule = 2)

      # Add interpolation result available data:
      dataSim <- data.table::rbindlist(list(dataSim,
                                            data.frame(TIME    = interpResult$x,
                                                       PL      = interpResult$y,
                                                       stringsAsFactors = FALSE)))

      # Re-order needed:
      dataSim <- dataSim[order(dataSim$TIME),]
    }
  }else{
    warning("evaluate_BreakthroughEvent - The dataset has less than two measurement points:",
            "Impossible to interpolate.")
  }


  #--------------------------------------#
  # STEP 4: Calculate criteria and adjust PL ----
  #--------------------------------------#

  # timepoints of breakthrough:
  idxRecrud <- dataSim$PL > LLOQ
  if (any(idxRecrud)) {
    Tevent <- dataSim$TIME[min(which(idxRecrud))]
    Event  <- TRUE
  } else {
    Tevent <- max(dataSim$TIME)
    Event  <- FALSE
  }


  #--------------------------------------#
  # STEP 5: Prepare output ----
  #--------------------------------------#

  res <- data.frame(NAME = "Breakthrough",
                    TIME = Tevent,
                    VALUE = Event,
                    stringsAsFactors   = FALSE)
  # Output:
  res
}



#' Generate Bite Times Using a Poisson Process
#'
#' This function generates bite event times over a specified observation period using 
#' a Poisson process, based on a mean bite rate (events per day).
#'
#' @param force_of_infection Numeric. Should represent an assumption about the number of infectious bites
#' that occur over a given time period. This can derived, for example, from estimated annual entomological
#' inoculation rates.This number should represent a yearly count : the provided parameter will be 
#' in the form x / 365 /24 to provide a daily rate 
#' @param Tend Numeric. The observation period (in hours) during which bite events will be simulated.
#' If being used as part of a simulation, should likely reflect the simulation's maximum run time 
#'
#' @return A numeric vector of bite times (in hours) within the specified observation period.
#'         Each element represents the cumulative time (in hours) at which a bite event occurs.
#' @details 
#' The function calculates a rate parameter for a Poisson process that represents 
#' the mean number of bite events per hour, derived from a yearly figure of infectious bites (`force_of_infection`). 
#' Using this rate parameter, it generates a cumulative sum of exponentially-distributed
#' inter-arrival times for bite events.
#'
#' If the generated sequence of bite events does not span the entire observation period (`Tend`), 
#' additional events are appended to ensure sufficient coverage, stopping when the maximum 
#' bite event time exceeds the observation period.
#' 
#' The function returns only those events that occur within the observation period.
#'
#' @examples
#' # Generate bite times with an average of 10 bites per year over a 24-hour period
#' generate_bite_time_poisson(force_of_infection = 10, Tend = 24)
#'
#' # Generate bite times with an average of 5 bites per year over a 48-hour period
#' generate_bite_time_poisson(force_of_infection = 5, Tend = 48)
#'
#' @export
#' @author Sam Jones (MMV, \email{joness@@mmv.org})
#' @family Chemoprotection
generate_bite_time_poisson <- function(force_of_infection, Tend) {
  rateParam <- force_of_infection / 365 / 24
  nObs_i <- qpois(1 - 1e-12, lambda = rateParam * Tend)
  
  vec_biteTimes_ijk <- cumsum(rexp(nObs_i, rate = rateParam))
  
  while (force_of_infection > 0 && max(vec_biteTimes_ijk) < Tend) {
    vec_biteTimes_ijk <- c(
      vec_biteTimes_ijk,
      cumsum(rexp(nObs_i, rate = rateParam)) + max(vec_biteTimes_ijk)
    )
  }
  
  # Select only bite events within the observation period Tend
  vec_biteTimes_ijk <- vec_biteTimes_ijk[vec_biteTimes_ijk <= Tend]
  
  return(vec_biteTimes_ijk)
}

#' Add Bite Events to Trial Files
#'
#' This function adds simulated mosquito bite events to existing trial files using a specified
#' bite time generator function. It copies a pre-defined set of trial files from a baseline folder
#' to an output folder, and then updates the trial files with these bite events. The bite times
#' are generated using the provided bite time generator function, which accepts arguments specified
#' in `biteTimeGeneratorArgs`. This allows for flexible and customized bite time generation
#' according to the requirements of your simulation. 
#'
#' @param baselineFolder Character string. The directory path where base simulation output files are
#'        stored and will be updated. This should be the folder containing the baseline trial files
#'        before adding bite events; these baseline files should contain PKPD parameters for all subjects
#' @param outputFolder Character string. The directory path where the updated simulation output files
#'        with added bite events will be stored. The function will copy the contents of `baselineFolder`
#'        to this location and then update the trial files with the new bite events.
#' @param TimePreDose Numeric. The time before the first dose (in hours). If greater than zero,
#'        an initial zero dose event is added at time zero. Including `TimePreDose` permits infectious
#'        bites to occur before the therapy is given.
#' @param bite_ADM Integer. The administration (ADM) code representing the bite input in the
#'        simulation structural model used in the trial files you are updating. Example: if your
#'        structural model has parasites being INPUT3, bite_ADM should = 3 
#' @param bite_AMT Numeric. The amount of malaria parasites introduced per bite event. Default: 640, 
#'        as estimated from sporozoite human challenge studies of DSM265
#' @param nsubjs Integer. The number of subjects in each trial. This is used for indexing and
#'        adjusting subject IDs when updating trial files. This must be equal to nsubjs in `baselineFolder`
#' @param ntrials Integer. The number of trials. This is used for indexing and adjusting trial IDs
#'        when updating trial files. This must be equal to ntrials in `baselineFolder`
#' @param biteTimeGenerator Function. A function used to generate bite times. This function should
#'        accept arguments specified in `biteTimeGeneratorArgs`. The bite time generator function
#'        should be designed to generate bite times according to your simulation requirements.
#' @param biteTimeGeneratorArgs List. A named list of arguments to be passed to the `biteTimeGenerator`
#'        function. These arguments customize the bite time generation process and should match the
#'        parameters expected by your `biteTimeGenerator` function.
#'        
#' @references
#' Cherkaoui-Rbati MH, Andenmatten N, Burgert L, et al. (2023). A pharmacokinetic-pharmacodynamic model 
#' for chemoprotective agents against malaria. \emph{CPT Pharmacometrics Syst Pharmacol}, \bold{12}, 50-61. \doi{10.1002/psp4.12875}
#' 
#' @details 
#' The function operates in two main steps:
#'
#' **Step 1:** It copies the baseline trial files from `baselineFolder` to `outputFolder`. This
#' ensures that PKPD parameters, regressors, and dosing events are retained.
#'
#' **Step 2:** It processes each trial file in `outputFolder` to add bite events. For each subject
#' in each trial, it generates bite times using the provided `biteTimeGenerator` function and the
#' arguments specified in `biteTimeGeneratorArgs`. The bite events are then added to the trial's
#' event table.
#'
#' If `TimePreDose` is greater than zero, a zero dose event is added at time zero to ensure that
#' simulations start at time zero and to allow bites to occur before dosing.
#' 
#' @export
#' @importFrom MMVbase file.copyMMV
#' @author Sam Jones (MMV, \email{joness@@mmv.org})
#' @family Chemoprotection
add_bites_to_trial_files <- function(baselineFolder,
                                            outputFolder,
                                            TimePreDose, 
                                            bite_ADM,
                                            bite_AMT = 640,
                                            nsubjs,
                                            ntrials, 
                                            biteTimeGenerator,
                                            biteTimeGeneratorArgs) {
  
  #--------------------------------------#
  # STEP 1: Copy required trial files ----
  #--------------------------------------#
  # Copy baseline files to the new folder
  file.copyMMV(
    from = baselineFolder,
    to = outputFolder,
    recursive = TRUE,
    overwrite = TRUE
  )
  
  # List trial files
  trialFilenames <- list.files(outputFolder, pattern = "trial.*.rds", full.names = TRUE)
  
  #--------------------------------------#
  # STEP 2: Add bites and adjust trial files ----
  #--------------------------------------#
  # Process each trial file
  for (trialFilename in trialFilenames) {
    trial <- MMVmalaria:::loadTrialFile(trialFilename)
    
    # Get parameter names
    parNames <- setdiff(
      names(trial$eventTable),
      c("ScenID", "ExpID", "DoseID", "TrialID", "ID", "TIME", "ADM", "AMT", "TINF")
    )
    
    if (trial$DoseID == 1) {
      # Dose ID = 1 
      # Sample Bites across patients
      biteTimes_list <- list()
      ID_list <- unique(trial$eventTable$ID)
      
      for (k in seq_along(ID_list)) {
        ID_k <- ID_list[k]
        
        # Generate bites with biteTimeGenerator
        vec_biteTimes_k <- do.call(biteTimeGenerator, biteTimeGeneratorArgs)
        
        n_Bites_k <- length(vec_biteTimes_k)
        
        # Biting events
        biteTimes_k <- data.frame(
          "ScenID"  = rep(trial$ScenID, n_Bites_k),
          "ExpID"   = rep(trial$ExpID, n_Bites_k),
          "DoseID"  = rep(trial$DoseID, n_Bites_k),
          "TrialID" = rep(trial$TrialID, n_Bites_k),
          "ID"      = rep(ID_k, n_Bites_k),
          "TIME"    = vec_biteTimes_k,
          "ADM"     = rep(bite_ADM, n_Bites_k),
          "AMT"     = rep(bite_AMT, n_Bites_k),
          "TINF"    = rep(0.0001, n_Bites_k),
          stringsAsFactors = FALSE
        )
        
        # Get Parameters
        parPKPD_k <- as.data.frame(matrix(NA, nrow = n_Bites_k, ncol = length(parNames)))
        names(parPKPD_k) <- parNames
        
        # Concatenate
        biteTimes_k <- cbind(biteTimes_k, parPKPD_k)
        
        # Add a dose of 0 at time 0 to have all simulations starting at 0
        if (TimePreDose > 0) {
          dose0_k <- data.frame(
            "ScenID"  = trial$ScenID,
            "ExpID"   = trial$ExpID,
            "DoseID"  = trial$DoseID,
            "TrialID" = trial$TrialID,
            "ID"      = ID_k,
            "TIME"    = 0,
            "ADM"     = seq(1, bite_ADM - 1),
            "AMT"     = 0,
            "TINF"    = 0.0001,
            stringsAsFactors = FALSE
          )
          
          # Get PKPD parameters to be added for TIME=0
          parPKPD0_k <- trial$eventTable[
            trial$eventTable$ID == ID_k & trial$eventTable$TIME == TimePreDose,
            parNames
          ]
          
          # Concatenate
          dose0_k <- cbind(dose0_k, parPKPD0_k)
          biteTimes_k <- rbind(dose0_k, biteTimes_k)
        }
        
        # Add to list
        biteTimes_list[[k]] <- biteTimes_k
      }
      
      # Combine all bite times
      biteTimes <- data.table::rbindlist(biteTimes_list)
      
      # Add bite events to the event table
      if (TimePreDose > 0) {
        trial$eventTable[, parNames] <- NA
      }
      trial$eventTable <- rbind(trial$eventTable, biteTimes)
      trial$eventTable <- trial$eventTable[
        order(trial$eventTable$ID, trial$eventTable$TIME, trial$eventTable$ADM),
      ]
      
      # Save updated event table and bite times
      MMVmalaria:::saveToTrialFile(
        outputFolder,
        trial$eventTable,
        listMemberName = "eventTable",
        FLAGreuseStoredTrialSummaryFiles = FALSE
      )
      MMVmalaria:::saveToTrialFile(
        outputFolder,
        biteTimes,
        listMemberName = "biteTimes",
        FLAGreuseStoredTrialSummaryFiles = FALSE
      )
      
    } else if (trial$DoseID > 1) {
      # Open the trial file corresponding to DoseID == 1
      DoseID1 <- 1
      ScenID1 <- 1
      trialFilename_DoseID1 <- file.path(
        dirname(trialFilename),
        paste0("trial_", ScenID1, "_", trial$ExpID, "_", DoseID1, "_", trial$TrialID, ".rds")
      )
      
      if (file.exists(trialFilename_DoseID1)) {
        # Load the DoseID == 1 trial file
        trial_DoseID1 <- MMVmalaria:::loadTrialFile(trialFilename_DoseID1)
        
        # Get the bite events
        biteTimes <- subset(trial_DoseID1$eventTable, ADM == bite_ADM)
        biteTimes$ScenID <- trial$ScenID
        biteTimes$DoseID <- trial$DoseID
        biteTimes$ID <- biteTimes$ID + (trial$ScenID - 1) * nsubjs * ntrials
        
        # Add a dose of 0 at time 0 to have all simulations starting at 0
        if (TimePreDose > 0) {
          dose0 <- subset(
            trial_DoseID1$eventTable,
            TIME == 0 & ADM %in% seq(1, bite_ADM - 1) & AMT == 0
          )
          dose0$ScenID <- trial$ScenID
          dose0$DoseID <- trial$DoseID
          dose0$ID <- dose0$ID + (trial$ScenID - 1) * nsubjs * ntrials
          
          # Concatenate
          biteTimes <- rbind(dose0, biteTimes)
        }
        
        # Add bite events to the event table
        if (TimePreDose > 0) {
          trial$eventTable[, parNames] <- NA
        }
        trial$eventTable <- rbind(trial$eventTable, biteTimes)
        trial$eventTable <- trial$eventTable[
          order(trial$eventTable$ID, trial$eventTable$TIME, trial$eventTable$ADM),
        ]
        
        # Save updated event table
        MMVmalaria:::saveToTrialFile(
          outputFolder,
          trial$eventTable,
          listMemberName = "eventTable",
          FLAGreuseStoredTrialSummaryFiles = FALSE
        )
      } else {
        # Corresponding DoseID == 1 trial file does not exist
        warning(paste("DoseID == 1 trial file not found for", trialFilename))
        next  # Skip to next trial file
      }
    }
  }  # End of loop over trial files
}  # End of add_bites_to_trial_files function