% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10-Simulations.R
\name{getModelParameters_MMVmalariaXLS}
\alias{getModelParameters_MMVmalariaXLS}
\title{Get Model Parameters from a GPF/XLS file}
\usage{
getModelParameters_MMVmalariaXLS(filename)
}
\arguments{
\item{filename}{A character string with the path to a GPF file.}
}
\value{
A list with the population, IIV, Betas and residual parameters, each formatted into a data.frame.
}
\description{
Load into a list the the population, IIV, Betas and residual parameters, each formatted into a data.frame from
a GPF/XLS file. Used in [getModelParameters_MMVmalariaProject].
}
\seealso{
Other Simulations: 
\code{\link{applyCureThreshold}()},
\code{\link{generateTrialSim}()},
\code{\link{generate_GPFfromOldIQRproject}()},
\code{\link{generate_GPFfromParametersTable}()},
\code{\link{getModelParameters_MMVmalariaCSV}()},
\code{\link{getModelParameters_MMVmalariaProject}()},
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
\code{\link{simulate_PKPDseparate}()},
\code{\link{simulate_VirtualTrials}()}
}
\author{
Mohammed H. Cherkaoui (MMV)
}
\concept{Simulations}
