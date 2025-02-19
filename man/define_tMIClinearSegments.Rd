% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11-FitLines.R
\name{define_tMIClinearSegments}
\alias{define_tMIClinearSegments}
\title{define_tMIClinearSegments
Determine time of MIC based on fitted linear segments}
\usage{
define_tMIClinearSegments(FIT)
}
\arguments{
\item{FIT}{data frame defining linear segments fitted
to clearance and regrowth data of individual
(combined output of \link{fit_LinearModel} and \link{fit_LinearSegments}).
Phase column needs to be defined (see Details)}
}
\value{
data frame with last knot of clearance phase,
first knot of regrowth phase
and interpolated time of MIC defined by Xknot, Yknot
}
\description{
Determines intersection of last segment of clearance phase and line of re-growth
phase as the time of MIC.
}
\details{
Input data frame contains the outputs of \link{fit_LinearModel}
and \link{fit_LinearSegments} and contain coordinates of knots of linear segments.
Phases are differentiated by 'Phase' column with KILL for clearance phase
and RECRUD for re-growth phase.
}
\seealso{
Other Model-free PD assessment by linear segments: 
\code{\link{define_tLLOQ}()},
\code{\link{fit_LinearModel}()},
\code{\link{fit_LinearSegments}()},
\code{\link{get_CLmaxTlag}()},
\code{\link{get_DataClearance}()},
\code{\link{get_DataGrowth}()},
\code{\link{get_Growth}()},
\code{\link{get_MICextrapolPK}()},
\code{\link{get_MICsimPK}()},
\code{\link{get_linepar}()},
\code{\link{plot_LinSegments}()},
\code{\link{plot_tMICextrapolation}()},
\code{\link{plot_tMICextrapolationRange}()},
\code{\link{simulate_indivPK}()}
}
\author{
Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
}
\concept{Model-free PD assessment by linear segments}
