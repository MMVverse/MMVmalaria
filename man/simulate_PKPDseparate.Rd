% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10-Simulations.R
\name{simulate_PKPDseparate}
\alias{simulate_PKPDseparate}
\title{simulate_PKPDseparate}
\usage{
simulate_PKPDseparate(
  Doses,
  Ndoses,
  PKmodelFolder,
  PDmodelFolder,
  PKPDmodelFile,
  filename = NULL,
  EC50adj = 1,
  DOSEcovariate = NULL,
  covariates = NULL,
  simtime = 0:800,
  Npop = 10,
  Nind = 10,
  FLAGsamplePK = 1,
  FLAGreturnObject = FALSE,
  LLOQ.PK = NULL,
  LLOQ.PD = NULL,
  cureThreshold = NULL,
  PLbaseSample = NULL,
  setting = NULL
)
}
\arguments{
\item{filename}{Default: \code{NULL}}

\item{EC50adj}{Default: 1}

\item{DOSEcovariate}{Default: \code{NULL}}

\item{covariates}{Default: \code{NULL}}

\item{simtime}{Default: 0:800}

\item{Npop}{Default: 10}

\item{Nind}{Default: 10}

\item{FLAGsamplePK}{Default: 1}

\item{FLAGreturnObject}{Default: FALSE}

\item{LLOQ.PK}{Default: \code{NULL}}

\item{LLOQ.PD}{Default: \code{NULL}}

\item{cureThreshold}{Default: \code{NULL}}

\item{PLbaseSample}{Default: \code{NULL}}

\item{setting}{Default: \code{NULL}}
}
\description{
simulate_PKPDseparate
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
\code{\link{simulate_VirtualTrials}()}
}
\author{
Anne Kümmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
}
\concept{Simulations}
