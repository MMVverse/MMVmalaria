% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/18-Simulate_VirtualTrials.R
\name{simulate_VirtualTrials}
\alias{simulate_VirtualTrials}
\title{Simulate virtual clinical trials for a set of experimental setting and dose combinations}
\usage{
simulate_VirtualTrials(
  projectPathList,
  modelFile,
  outputFolder,
  regressorExtra = NULL,
  Doses,
  T_FirstDose = 0,
  nbrDoses = 1,
  timeInterval = 24,
  DOSEcovariate = c(INPUT1 = "DOSELEVEL"),
  N_Trial = 100,
  N_SubjPerTrial = 10,
  IndCovSample = NULL,
  ExpCovSample = NULL,
  FLAG_SAMPLE = 1,
  replace = TRUE,
  setting = list(.paropts = list(.inorder = TRUE, .packages = c("IQRtools",
    "MMVmalaria"))),
  Tk0 = c(INPUT1 = "Tk0"),
  Fpediatric = NULL,
  Nparallel = 1,
  simMode = c("Clinical", "Compare")[1],
  FLAGrunSimulation = TRUE
)
}
\arguments{
\item{projectPathList}{List containing the path of all models to use for sampling. Each element of the list should point to an 'IQRnlmeProject', 'IQRsysProject' and/or a 'GPF' file (see \link[IQRtools:GPF]{IQRtools::GPF}).}

\item{modelFile}{Path to the PKPD model equation file to be used for the simulations.}

\item{outputFolder}{a character string denoting a path to a directory where a .rds file containing the resulting list summarizing each simulated trial is to be saved.}

\item{regressorExtra}{a character vector containing the names of additional regressors to be used in the simulations. Default: NULL.}

\item{Doses}{Numeric vector of the doses to tes,t or a dosing table containing columns "DoseID", "TIME", "ADM",
"AMT", "TINF", such that multiple dosing regimen correspond to different DoseID values. Regimens with time-varying dose
are supported unless, the dose is a covariate of the model. Note that the TINF column is neglected and overriden by
1e-4. If there is a need to specify infusion time for a given INPUT, this can be done by specifying the Tk0 argument.}

\item{T_FirstDose}{Numeric defining time of first dose. Might be diferent than initial time of simulation (Default: 0).}

\item{nbrDoses}{Numeric Vector. If length is 1, all doses in 'Doses' will have the same number of doses, otherwise should be of length 'Doses'.}

\item{timeInterval}{Numeric defining the interval between doses (Note: Sssumed to be the same for all doses)}

\item{DOSEcovariate}{a character vector containing the names of the covariates to be used as dose covariate. Default: \code{c("INPUT1"="DOSELEVEL")}.}

\item{IndCovSample}{Dataframe or Function to generate individual covariate (e.g. WT0, BMI, SEX, etc). Default: NULL.}

\item{ExpCovSample}{Dataframe of experimental covariates to account for (e.g. Formulation, Fasted/Fed, etc). Default: NULL.}

\item{setting}{a named list with elements among the following:
\itemize{
\item FLAGreturnSimPKPD A logical (default: FALSE) indicating if simPKPD should be included in the resulting list.
Note that setting this argument to TRUE is memory intensive and can cause memory crashes.
\item FLAGsaveSimPKPDToFile a logical (default FALSE) indicating if the simPKPD argument should be saved in
\code{outputFolder} in the result .rds file for the corresponding trial. Setting this to TRUE while outputFolder is NULL
will raise in an error.
\item FLAGreuseStoredTrialSummaryFiles a logical (default FALSE) indicating if previously stored .rds trial summary files in \code{outputFolder} should
be reused instead of rerunning each trial simulation.
\item fun_SummarizeTrial a character string denoting the name of a summarizing function to be
called on the simulation results for an experimental setting, a dose, and a virtual trial population. Default:
\code{"MMVmalaria::summarizeTrial_Basic"}. This function should accept at least a data.frame
and a character string argument, namely a data.frame containing the PKPD model simulation result and individual model parameters as
returned from the function \code{\link{IQRtools::sim_IQRmodel}}) and a filename containing the path to the trial .rds file on the disk.
The .rds file can be loaded into a list and used to get the individual covariate values for each subject in the trial. Example code:
\code{trial <- readRDS(trialFileName); simPKPDcovs <- merge(simPKPD, trial$covariates, by = intersect(names(simPKPD), names(trial$covariates)))}.
Any additional arguments of the function can be specified in the list argument
\code{args_SummarizeTrial}. Note that these arguments will have the same values over all trials. The function is expected to return
a named list with at least three data.frame members: summaryPKPD.ByTrial, summaryClinEnd.ByIndiv, summaryClinEnd.ByTrial. In addition,
the function can store result objects into .rds files for post-processing.
See \code{\link{summarizeTrial_Basic}} for an example.
\item args_SummarizeTrial a named list or NULL (default) containing additional arguments
to be passed to \code{fun_SummarizeTrial}. See documentation of \code{\link{summarizeTrial_Basic}}
for possibilities.
\item .paropts a named list specifying arguments passed to \code{\link{foreach::foreach}}. By default only this member is set to
\code{.paropts = list(.inorder = TRUE,
                      .packages = c ("IQRtools","MMVmalaria"))}.
}}

\item{Tk0}{a named character vector with names INPUT1, INPUT2, etc, containing the names of corresponding Tk0 parameters.
Note that the value of the Tk0 parameter for an input will overwrite the value of TINF in the dosing object.
Default \code{c("INPUT1"="Tk0")}. INPUTs in this argument are assumed to have 0-order absorption, i.e. Tk0 is
the time for the dose to be fully absorbed at a constant rate. This is suitable to simulate infusion at constant rate.
The Tk0 parameters should be either present in the GPF file or available as regressors in the regressorExtra argument.}

\item{simMode}{To chose which mode to use for simulations between "Clinical" (New set of individuals sampled for each scenario - Default) and
"Compare" (Same set of individuals used between each scenario adjusted by \code{ExpCovSample} and Dose-regimen as defined in \code{Doses}.).}

\item{FLAGrunSimulation}{Run simulations after PK and PD parameters sampling}
}
\value{
a named list with the following members:
\itemize{
\item summaryPKPD.ByTrial a data.frame resulting from rbinding all summaryPKPD.ByTrial data.frames returned by
\code{setting$fun_SummarizeTrial}.
\item summaryClinEnd.ByIndiv a data.frame resulting from rbinding all summaryPKPD.ByIndiv data.frames returned by
\code{setting$fun_SummarizeTrial}.
\item summaryClinEnd.ByTrial a data.frame resulting from rbinding all summaryClinEnd.ByTrial data.frames returned by
\code{setting$fun_SummarizeTrial}.
\item Covariates a data.frame
\item PKPDparameters a data.frame
\item eventTable a data.frame
\item eventData a data.frame
\item regression
\item abs0inputs
\item abs0Tk0param
\item setting
}
}
\description{
Simulate virtual clinical trials for a set of experimental setting and dose combinations
}
\details{
For shortness, the word 'scenario' is used to indicate a combination of an experimental setting and a dose.
Having two types of covariates, namely, individual and experimental, allows to use the same
individuals for different experimental conditions, and, therefore, to see the difference between
experimental conditions with more clarity.
}
\seealso{
Other Simulations: 
\code{\link{applyCureThreshold}()},
\code{\link{generateTrialSim}()},
\code{\link{generate_GPFfromOldIQRproject}()},
\code{\link{generate_GPFfromParametersTable}()},
\code{\link{getModelParameters_MMVmalariaCSV}()},
\code{\link{getModelParameters_MMVmalariaProject}()},
\code{\link{getModelParameters_MMVmalariaXLS}()},
\code{\link{getPopParameters_MMVmalariaCSV}()},
\code{\link{getPopParameters_MMVmalariaProject}()},
\code{\link{getPopParameters_MMVmalariaXLS}()},
\code{\link{predict_HumanChallenge}()},
\code{\link{sample_CustomParameters}()},
\code{\link{sample_MMVmalariaCSV}()},
\code{\link{sample_MMVmalariaCSV_NEW}()},
\code{\link{sample_MMVmalariaProject}()},
\code{\link{sample_MMVmalariaXLS}()},
\code{\link{simTimeRecrudescenceCombo}()},
\code{\link{simulate_ComboMouse2Human}()},
\code{\link{simulate_PKPDseparate}()}
}
\author{
Mohammed H. Cherkaoui, \email{cherkaouim@mmv.org}, Venelin Mitov \email{venelin.mitov@intiquan.com}
}
\concept{Simulations}
