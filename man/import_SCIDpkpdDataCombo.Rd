% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-DataPreparation.R
\name{import_SCIDpkpdDataCombo}
\alias{import_SCIDpkpdDataCombo}
\title{import_SCIDpkpdDataCombo}
\usage{
import_SCIDpkpdDataCombo(
  dataFile,
  Compound1 = NULL,
  Compound2 = NULL,
  PKsheets = NULL,
  PKranges = NULL,
  PKmapping,
  PKlloq = c(0.001, 0.001),
  PKlloqIdentifier = c("BLQ", "<LLOQ", "< LLOQ", "BLOQ"),
  PKfactor = c(0.001, 0.001),
  PDsheet = NULL,
  PDrange = NULL,
  PDmapping,
  PDlloq = 0.01,
  PDlloqIdentifier = c("BLQ", "<0.01"),
  PDname = "Parasitemia Total",
  study = NULL,
  studyID = NULL,
  DayOfFirstDrugAdmin = NULL,
  centerNumber = -1,
  centerName = "",
  visitNumber = -1,
  intervaldose = NULL,
  DoseTimeSheet = "DosingTimeTable",
  DoseTimeRange = NULL,
  DoseTimeMapping,
  DoseSheet = "DrugTreatmentTable",
  DoseRange = NULL,
  DoseMapping = NULL,
  positiveQC = NULL,
  route = NULL,
  importVersion = "oldTAD"
)
}
\arguments{
\item{Compound1}{Default: NULL}

\item{Compound2}{Default: NULL}

\item{PKsheets}{Default: NULL}

\item{PKranges}{Default: NULL}

\item{PKlloq}{Default: c(0.001, 0.001)}

\item{PKlloqIdentifier}{Default: c("BLQ", "<LLOQ", "< LLOQ", "BLOQ")}

\item{PKfactor}{Default: c(0.001, 0.001)}

\item{PDsheet}{Default: NULL}

\item{PDrange}{Default: NULL}

\item{PDlloq}{Default: 0.01}

\item{PDlloqIdentifier}{Default: c("BLQ", "<0.01")}

\item{PDname}{Default: 'Parasitemia Total'}

\item{study}{Default: NULL}

\item{studyID}{Default: NULL}

\item{DayOfFirstDrugAdmin}{Default: NULL}

\item{centerNumber}{Default: -1}

\item{centerName}{Default: ''}

\item{visitNumber}{Default: -1}

\item{intervaldose}{Default: NULL}

\item{DoseTimeSheet}{Default: 'DosingTimeTable'}

\item{DoseTimeRange}{Default: NULL}

\item{DoseSheet}{Default: 'DrugTreatmentTable'}

\item{DoseRange}{Default: NULL}

\item{DoseMapping}{Default: NULL}

\item{positiveQC}{Default: NULL}

\item{route}{Default: NULL}

\item{importVersion}{Default: 'oldTAD'}
}
\description{
import_SCIDpkpdDataCombo
}
\seealso{
Other Data Preparation: 
\code{\link{censorGametocyte}()},
\code{\link{convertGametocytes}()},
\code{\link{convert_Gametocytes}()},
\code{\link{import_SCIDpkpdData}()},
\code{\link{import_SCIDpkpdDataCombo_GSK}()},
\code{\link{import_SCIDpkpdDataCombo_TAD}()},
\code{\link{import_SCIDpkpdDataCombo_oldTAD}()},
\code{\link{import_SCIDpkpdData_GSK}()},
\code{\link{import_SCIDpkpdData_TAD}()},
\code{\link{import_SCIDpkpdData_oldTAD}()},
\code{\link{load_AfricanPediatricMalariaPopulation2to5}()},
\code{\link{load_MalariaPopulation}()},
\code{\link{load_PediatricParameters}()}
}
\author{
Aline Fuchs (MMV), Anne Kuemmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
}
\concept{Data Preparation}
