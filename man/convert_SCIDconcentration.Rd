% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/09-KeyParameters.R
\name{convert_SCIDconcentration}
\alias{convert_SCIDconcentration}
\title{convert_SCIDconcentration}
\usage{
convert_SCIDconcentration(
  ConcToConv,
  HuRbp,
  HuPPB,
  MouseRbp,
  MousePPB,
  HuH = 0.45,
  MouseH = 0.45,
  huRBCinSCID = 0.6,
  scidH = 0.8,
  ConvertTo = c("Human Blood", "Human Plasma", "Mouse Blood", "Mouse Plasma", "Free")[1]
)
}
\arguments{
\item{HuH}{Default: 0.45}

\item{MouseH}{Default: 0.45}

\item{huRBCinSCID}{Default: 0.6}

\item{scidH}{Default: 0.8}

\item{ConvertTo}{Default: 'Human'}
}
\description{
convert_SCIDconcentration
}
\seealso{
Other Key Parameters: 
\code{\link{assess_MICMPC90PRR48}()},
\code{\link{estimate_AUC}()},
\code{\link{estimate_LambdaZ}()},
\code{\link{evaluate_APR}()},
\code{\link{evaluate_APR28}()},
\code{\link{formulaPKpar_MacroToMicro}()},
\code{\link{formulaPKpar_MacroToPrimary}()},
\code{\link{formulaPKpar_MicroToMacro}()},
\code{\link{formulaPKpar_MicroToPrimary}()},
\code{\link{formulaPKpar_PrimaryToMacro}()},
\code{\link{formulaPKpar_PrimaryToMicro}()},
\code{\link{formulaPKpar_PrimaryToSecondary}()},
\code{\link{formula_SCIDrbp}()},
\code{\link{getApparentKeys}()},
\code{\link{getCure}()},
\code{\link{getKeysEMAX}()},
\code{\link{getTimeAboveMIC}()},
\code{\link{getTimeAboveMICsim}()},
\code{\link{getTimeAboveMPC90}()},
\code{\link{getTimeAboveMPC90sim}()},
\code{\link{getTimeKRaboveGR}()},
\code{\link{getTimeRecrudescence}()},
\code{\link{get_LagTime}()},
\code{\link{get_NumberOfPKcpt}()},
\code{\link{get_PRRtot}()},
\code{\link{plot_keyPDparametersMMV}()}
}
\author{
Aline Fuchs (MMV), Mohammed H. Cherkaoui (MMV)
}
\concept{Key Parameters}
