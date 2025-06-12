context("test-predictHumanDoses")

library(IQRtools)
library(MMVmalaria)

projectPath <- system.file("extdata", "test-predictHumanDoses", "MODEL01", "GPF_PKPDcombo.xlsx", package = "MMVmalaria")
modelFile <- system.file("extdata", "test-predictHumanDoses", "MODEL01", "model_PKPDcombo_Bliss.txt", package = "MMVmalaria")

covariates <- MalariaPopulation
# rename some of the columns to names used in models
covariates$PLbase <- covariates$PLBASE_log10_perml
covariates$WT0 <- covariates$WEIGHT_kg

# Fpediatric will be used for dose scaling
covariates$Fpediatric <- sapply(covariates$WT0, function(w) {
  if(w <= 10) {
    0.25
  } else if(w <= 15) {
    0.375
  } else if(w <= 24) {
    0.5
  } else if(w <= 35) {
    0.75
  } else {
    1
  }
})

# Function calculating the time above MIC for COMPOUND 1
evaluate_TimeAboveMICx1 <- function(sim_results, parameters) {
  parameters["EMAX"] <- parameters["EMAXx1"]
  parameters["EC50"] <- parameters["EC50x1"]
  parameters["hill"] <- parameters["hillx1"]
  dfParam <- as.data.frame(as.list(parameters))
  MIC     <- MMVmalaria::getKeysEMAX(x = dfParam)[["MIC"]]
  sim_results$Cc <- sim_results$Ccx1
  MMVmalaria::getTimeAboveMIC(sim_results, MIC = MIC)[["tMIC"]]
}

# Function calculating the time above MIC for COMPOUND 2
evaluate_TimeAboveMICx2 <- function(sim_results, parameters) {
  parameters["EMAX"] <- parameters["EMAXx2"]
  parameters["EC50"] <- parameters["EC50x2"]
  parameters["hill"] <- parameters["hillx2"]
  dfParam <- as.data.frame(as.list(parameters))
  MIC     <- MMVmalaria::getKeysEMAX(x = dfParam)[["MIC"]]
  sim_results$Cc <- sim_results$Ccx2
  MMVmalaria::getTimeAboveMIC(sim_results, MIC = MIC)[["tMIC"]]
}

N_Trial <- 20
N_SubjperTrial <- 10
doseInterval = c(10, 20000)

# tMIC  ----

# . CMP1 ----
dose_tMIC7_CMP1 <- predictDose_TimeAboveMIC(
  fun_EvaluateCriterion = evaluate_TimeAboveMICx1,
  targetTimeAboveMIC = 7 * 24,
  doseInterval = doseInterval,
  modelFile = modelFile,

  # sampling parameters
  projectPath = projectPath,

  N_Trial = N_Trial,
  N_SubjperTrial = N_SubjperTrial,

  covariates = covariates,

  # Specify PLbase and Fpediatric as regressors
  regressorExtra = c("PLbase", "Fpediatric"),

  simtime = seq(0, 600, by = 1),
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24
  )

print(dose_tMIC7_CMP1)

# . CMP2 ----
# need to specify the args_IQRdosing argument
# DOSE is a covariate on Vcx2
dose_tMIC7_CMP2 <- predictDose_TimeAboveMIC(
  fun_EvaluateCriterion = evaluate_TimeAboveMICx2,
  targetTimeAboveMIC = 7 * 24,
  doseInterval = c(10, 10000),
  modelFile = modelFile,

  # sampling parameters
  projectPath = projectPath,

  N_Trial = N_Trial,
  N_SubjperTrial = N_SubjperTrial,

  covariates = covariates,
  DOSEcovariate = c(INPUT2 = "DOSE"),
  Fpediatric = "Fpediatric",

  # Specify PLbase and Fpediatric as regressors
  regressorExtra = c("PLbase", "Fpediatric"),

  simtime = seq(0, 600, by = 1),
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,

  # Set ADM=2 to make it a dosing for CMP2
  args_IQRdosing = list(TIME = 0, ADM = 2, AMT = NA_real_, ADDL = 0, II = 24)
  )
print(dose_tMIC7_CMP2)


# dose_tMIC7_CMP2_2 <- predictDose_TimeAboveMIC(
#   fun_EvaluateCriterion = evaluate_TimeAboveMICx2,
#   targetTimeAboveMIC = 7 * 24,
#   doseInterval = c(10, 10000),
#   modelFile = modelFile,
#
#   # sampling parameters
#   projectPath = projectPath,
#
#   N_Trial = N_Trial,
#   N_SubjperTrial = N_SubjperTrial,
#
#   covariates = covariates,
#   # DOSEcovariate = c(INPUT2 = "DOSE"),
#   # Fpediatric = "Fpediatric",
#
#   # Specify PLbase and Fpediatric as regressors
#   regressorExtra = c("PLbase", "Fpediatric"),
#
#   simtime = seq(0, 600, by = 1),
#   Tk0 = NULL,
#   nbrDoses = 1,
#   timeInterval = 24,
#
#   # Set ADM=2 to make it a dosing for CMP2
#   args_IQRdosing = list(TIME = 0, ADM = 2, AMT = NA_real_, ADDL = 0, II = 24)
#   )
# print(dose_tMIC7_CMP2_2)
#
#
# dose_tMIC7_CMP2_3 <- predictDose_TimeAboveMIC(
#   fun_EvaluateCriterion = evaluate_TimeAboveMICx2,
#   targetTimeAboveMIC = 7 * 24,
#   doseInterval = c(10, 10000),
#   modelFile = modelFile,
#
#   # sampling parameters
#   projectPath = projectPath,
#
#   N_Trial = N_Trial,
#   N_SubjperTrial = N_SubjperTrial,
#
#   covariates = covariates,
#   DOSEcovariate = c(INPUT2 = "DOSE"),
#   # Fpediatric = "Fpediatric",
#
#   # Specify PLbase and Fpediatric as regressors
#   regressorExtra = c("PLbase", "Fpediatric"),
#
#   simtime = seq(0, 600, by = 1),
#   Tk0 = NULL,
#   nbrDoses = 1,
#   timeInterval = 24,
#
#   # Set ADM=2 to make it a dosing for CMP2
#   args_IQRdosing = list(TIME = 0, ADM = 2, AMT = NA_real_, ADDL = 0, II = 24)
#   )
#
# print(dose_tMIC7_CMP2_3)


# PRRtot ----

# . CMP1 ----
dose_PRRtot9_CMP1 <- predictDose_PRRtot(
  targetPRRtot = 9,
  doseInterval = doseInterval,
  modelFile = modelFile,

  # sampling parameters
  projectPath = projectPath,

  N_Trial = N_Trial,
  N_SubjperTrial = N_SubjperTrial,

  covariates = covariates,

  # Specify PLbase and Fpediatric as regressors
  regressorExtra = c("PLbase", "Fpediatric"),

  simtime = seq(0, 600, by = 1),
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24)

print(dose_PRRtot9_CMP1)



# . CMP1 with 100 mg CMP2 ----
dose_PRRtot9_CMP1_100CMP2 <- predictDose_PRRtot(
  targetPRRtot = 9,
  doseInterval = doseInterval,
  modelFile = modelFile,

  # sampling parameters
  projectPath = projectPath,

  N_Trial = N_Trial,
  N_SubjperTrial = N_SubjperTrial,

  covariates = covariates,
  DOSEcovariate = c(INPUT2 = "DOSE"),
  Fpediatric = "Fpediatric",

  # Specify PLbase and Fpediatric as regressors
  regressorExtra = c("PLbase", "Fpediatric"),

  simtime = seq(0, 600, by = 1),
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,

  # Set AMT for ADM=2 to 100 mg
  args_IQRdosing = list(TIME = c(0,0), ADM = c(1,2), AMT = c(NA_real_, 100), ADDL = 0, II = 24)
  )

print(dose_PRRtot9_CMP1_100CMP2)

# . CMP2 ----
dose_PRRtot9_CMP2 <- predictDose_PRRtot(
  targetPRRtot = 9,
  doseInterval = doseInterval,
  modelFile = modelFile,

  # sampling parameters
  projectPath = projectPath,

  N_Trial = N_Trial,
  N_SubjperTrial = N_SubjperTrial,

  covariates = covariates,
  DOSEcovariate = c(INPUT2 = "DOSE"),
  Fpediatric = "Fpediatric",

  # Specify PLbase and Fpediatric as regressors
  regressorExtra = c("PLbase", "Fpediatric"),

  simtime = seq(0, 600, by = 1),
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,

  # Set ADM=2 to make it a dose for CMP2
  args_IQRdosing = list(TIME = 0, ADM = 2, AMT = NA_real_, ADDL = 0, II = 24)
  )

print(dose_PRRtot9_CMP2)

# . CMP2 with 100 mg CMP1 ----
dose_PRRtot9_CMP2_100CMP1 <- predictDose_PRRtot(
  targetPRRtot = 9,
  doseInterval = doseInterval,
  modelFile = modelFile,

  # sampling parameters
  projectPath = projectPath,

  N_Trial = N_Trial,
  N_SubjperTrial = N_SubjperTrial,

  covariates = covariates,
  DOSEcovariate = c(INPUT2 = "DOSE"),
  Fpediatric = "Fpediatric",

  # Specify PLbase and Fpediatric as regressors
  regressorExtra = c("PLbase", "Fpediatric"),

  simtime = seq(0, 600, by = 1),
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,

  # Set ADM=2 to make it a dose for CMP2
  args_IQRdosing = list(TIME = c(0,0), ADM = c(1,2), AMT = c(100,NA_real_), ADDL = 0, II = 24)
  )

print(dose_PRRtot9_CMP2_100CMP1)
