% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-DataPreparation.R
\name{import_SCIDpkpdDataCombo_oldTAD}
\alias{import_SCIDpkpdDataCombo_oldTAD}
\title{import_SCIDpkpdDataCombo_oldTAD}
\usage{
import_SCIDpkpdDataCombo_oldTAD(
  dataFile,
  PKsheets = c("Blood levels 1", "Blood levels 2"),
  PKranges = NULL,
  PKmapping,
  PKlloq = c(0.001, 0.001),
  PKlloqIdentifier = c("BLQ", "<LLOQ", "< LLOQ", "BLOQ"),
  PKfactor = c(0.001, 0.001),
  PDsheet = "Parasitemia",
  PDrange = NULL,
  PDmapping,
  PDlloq = 0.01,
  PDlloqIdentifier = c("BLQ", "<0.01"),
  PDname = "Parasitemia Total",
  study = NULL,
  studyID = NULL,
  DayOfFirstDrugAdmin = 1,
  centerNumber = -1,
  centerName = "",
  visitNumber = -1
)
}
\arguments{
\item{PKsheets}{Default: c("Blood levels 1", "Blood levels 2")}

\item{PKranges}{Default: NULL}

\item{PKlloq}{Default: c(0.001, 0.001)}

\item{PKlloqIdentifier}{Default: c("BLQ", "<LLOQ", "< LLOQ", "BLOQ")}

\item{PKfactor}{Default: c(0.001, 0.001)}

\item{PDsheet}{Default: 'Parasitemia'}

\item{PDrange}{Default: NULL}

\item{PDlloq}{Default: 0.01}

\item{PDlloqIdentifier}{Default: c("BLQ", "<0.01")}

\item{PDname}{Default: 'Parasitemia Total'}

\item{study}{Default: NULL}

\item{studyID}{Default: NULL}

\item{DayOfFirstDrugAdmin}{Default: 1}

\item{centerNumber}{Default: -1}

\item{centerName}{Default: ''}

\item{visitNumber}{Default: -1}
}
\description{
import_SCIDpkpdDataCombo_oldTAD
}
\seealso{
Other Data Preparation: 
\code{\link{censorGametocyte}()},
\code{\link{convertGametocytes}()},
\code{\link{convert_Gametocytes}()},
\code{\link{import_SCIDpkpdData}()},
\code{\link{import_SCIDpkpdDataCombo}()},
\code{\link{import_SCIDpkpdDataCombo_GSK}()},
\code{\link{import_SCIDpkpdDataCombo_TAD}()},
\code{\link{import_SCIDpkpdData_GSK}()},
\code{\link{import_SCIDpkpdData_TAD}()},
\code{\link{import_SCIDpkpdData_oldTAD}()},
\code{\link{load_AfricanPediatricMalariaPopulation2to5}()},
\code{\link{load_MalariaPopulation}()},
\code{\link{load_PediatricParameters}()}
}
\author{
Mohammed H. Cherkaoui (MMV)
}
\concept{Data Preparation}
