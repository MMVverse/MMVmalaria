% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/17-PredictHumanDoses.R
\name{predictDose_Generic}
\alias{predictDose_Generic}
\title{Find the minimal dose satisfying a criterion within a given dose interval for a single subject PKPD model simulation}
\usage{
predictDose_Generic(
  fun_EvaluateCriterion,
  targetCriterionValue,
  doseInterval,
  model,
  parameters,
  simtime,
  Tk0 = NULL,
  nbrDoses = 1,
  timeInterval = 24,
  args_IQRdosing = list(TIME = simtime[1], ADM = 1, AMT = NA_real_, ADDL = nbrDoses - 1,
    II = timeInterval),
  addArgs_sim_IQRmodel = NULL,
  addArgs_uniroot = list(tol = 0.1),
  FLAG_BROWSER = FALSE,
  FLAG_RPROF = FALSE,
  FLAG_SAVECALLS_sim_IQRmodel = FALSE
)
}
\arguments{
\item{fun_EvaluateCriterion}{a function object or a character string denoting the name of a function accessible
from the calling (parent) environment. This function should have the following synopsis:
\code{function(sim_results, parameters) { ... }}, where \code{sim_results} is a data.frame resulting from a
previous call to \code{\link{IQRtools::sim_IQRmodel}} and \code{parameters} is a data.frame with a single row
and numeric columns named as the model parameters; the function should return a numeric value. If
\code{fun_EvaluateCriterion} is a character string, the function is searched for via
\code{mget(fun_EvaluateCriterion, mode = "function", ifnotfound = list(NULL), inherits = TRUE)[[fun_EvaluateCriterion]]}.
An error will be thrown if no such function is found.}

\item{targetCriterionValue}{a numeric value to compare the returned values from \code{fun_EvaluateCriterion} against.}

\item{doseInterval}{a numeric vector of length 2 denoting the dosing interval to search within.}

\item{model}{an IQRmodel object denoting the structural model to simulate, i.e. an object returned by a previous call to \code{\link{IQRtools::IQRmodel}}.}

\item{parameters}{Either a named numeric vector denoting the model parameter values to simulate or, in the case of a dose covariate, a named list
with elements as follows:
\itemize{
\item popParamValues: a data.frame as the popParamValues member of the list returned by \code{\link{sampleIndParamValues_MMVmalariaXLS}} (only the 1st row will be used).
\item sampledData: a data.frame as the sampledData member of the list returned by \code{\link{sampleIndParamValues_MMVmalariaXLS}}.
This data.frame should have an integer column 'ID'. Only the 1st row of this data.frame is going to be used.
\item randomEffects: a data.frame as the randomEffects member of the list returned by \code{\link{sampleIndParamValues_MMVmalariaXLS}} (only the 1st row will be used).
\item DOSEcovariate: a named character vector (see documentation of \code{\link{predictDose_GenericTrial}}).
\item Fpediatric: a character string (see documentation of \code{\link{predictDose_GenericTrial}}).
\item regressorExtra: a character vector (see documentation of \code{\link{predictDose_GenericTrial}}).
\item gpf: a GFP object (used to calculate typical individual values based on the covariates in sampledData and to calculate individual values based on
randomEffects and transformation type).
\item Fpediatric: (optional) a character string denoting a column in sampledData to scale adult equivalent doses, before applying covariate formulae.
\item FLAG_SAMPLE: an integer (see documentation of \code{\link{predictDose_GenericTrial}} for possible values)
}}

\item{simtime}{a sorted numeric vector denoting the simulation times (passed to \code{\link{IQRtools::sim_IQRmodel}}).}

\item{Tk0}{A named character vector (default: NULL). Example : c(INPUT1 = "Tk0", INPUT2 = "Tk0_2"), meaning that
the IQRdosing object for each simulation would be adjusted to have TINF = set to the model parameter Tk0 for ADM = 1.}

\item{nbrDoses}{an integer (default 1), denoting the number of doses to simulate}

\item{timeInterval}{a numeric (default 24), denoting the time interval between doses (in the model time units).}

\item{args_IQRdosing}{a named list (default NULL) with arguments to be passed to \code{\link{IQRtools::IQRdosing}}. For flexibility, this
can be used to overwrite some of the otherwise set by default arguments to IQRdosing, such as \code{TIME}, \code{ADM}, \code{ADDL}, \code{II}
(advanced use only).}

\item{addArgs_sim_IQRmodel}{a named list (default NULL) with additional arguments to be passed to \code{\link{IQRtools::sim_IQRmodel}}.}

\item{addArgs_uniroot}{a named list (default \code{list(tol = .1)}) with additional arguments to be passed to \code{\link{stats::uniroot}}.}

\item{FLAG_BROWSER}{logical (default FALSE) indicating if a browser should be started before each simulation (debugging purpose only).}
}
\value{
A (possibly infinite) numeric value. If the criterion value evaluated at \code{min(doseInteval) >= targetCriterionValue}, \code{min(doseInterval)} is returned;
if the criterion value evaluated at \code{max(doseInteval) < targetCriterionValue}, \code{Inf} is returned; otherwise, a value within the interval
\code{[min(doseInterval), max(doseInterval)]} is found using the function \code{\link{stats::uniroot}}. Note that the resulting value is undetermined if \code{criterion
function - targetCriterionValue} has more than 1 root within the interval \code{[min(doseInterval), max(doseInterval)]}.
}
\description{
Find the minimal dose satisfying a criterion within a given dose interval for a single subject PKPD model simulation
}
\details{
In the current implementation, the search for a dose is done using the \code{\link{stats::uniroot}} and in the log10-transformed dose interval.
This implies that the convergence tolerance gets bigger for higher doses.
}
\author{
Venelin Mitov, IntiQuan
}
