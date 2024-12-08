% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/19-IQRaux.R
\name{create_GDF_fromDesignSpec}
\alias{create_GDF_fromDesignSpec}
\title{Create an IQRdataGENEARL from a design specification}
\usage{
create_GDF_fromDesignSpec(design)
}
\arguments{
\item{design}{A list which is formatted like the \code{design}-argument of \code{\link[IQRtools:IQRpopEDdb]{IQRtools::IQRpopEDdb()}}}
}
\value{
An \code{\link[IQRtools:IQRdataGENERAL]{IQRtools::IQRdataGENERAL()}}
}
\description{
Populate the data with \code{\link[=populateGDFWithSimValues]{populateGDFWithSimValues()}}
}
\examples{
design <- list(
  time = list(
    group1 = list(
      OUTPUT1 = c(1, 2, 5, 10),
      OUTPUT2 = c(1, 5, 20)
    ),
    group2 = list(
      OUTPUT1 = c(1, 2, 5, 10),
      OUTPUT2 = c(1, 5, 20)
    )
  ),
  dosing = list(
    group1 = IQRdosing(
      TIME = 0, ADM = 1, AMT = 40, TINF = 0
    ),
    group2 = IQRdosing(
      TIME = 0, ADM = 1, AMT = 40, TINF = 0
    )),
  groupsize = list(
    group1 = 12, group2 = 3
  )
)
create_GDF_fromDesignSpec(design)
}
\seealso{
Other IQRdataGENERAL: 
\code{\link{populateGDFWithSimValues_consistent}()}

Other Simulation: 
\code{\link{createSimIQRmodelSpec}()},
\code{\link{getDefaultEventTable}()},
\code{\link{populateGDFWithSimValues_consistent}()},
\code{\link{sim_IQRmodelConsistent}()}
}
\author{
Daniel Lill (daniel.lill@intiquan.com)
}
\concept{IQRdataGENERAL}
\concept{Simulation}
