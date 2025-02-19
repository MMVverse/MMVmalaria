% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-DataPreparation.R
\name{import_SCIDpkpdDataCombo_TAD}
\alias{import_SCIDpkpdDataCombo_TAD}
\title{import_SCIDpkpdDataCombo_TAD
Import SCID Combo data from TAD}
\usage{
import_SCIDpkpdDataCombo_TAD(
  dataFile,
  Compound1 = NULL,
  Compound2 = NULL,
  PKsheets = "DrugConcentrationData",
  PKranges = NULL,
  PKmapping,
  PKlloq = c(1, 1),
  PKlloqIdentifier = c("BLQ", "<LLOQ", "< LLOQ", "BLOQ"),
  PKfactor = c(0.001, 0.001),
  PDsheet = "ParasitemiaData",
  PDrange = NULL,
  PDmapping,
  PDlloq = 0.01,
  PDlloqIdentifier = c("BLQ", "<0.01", "< 0,01", "< 0.01"),
  PDname = "Parasitemia Total",
  DayOfFirstDrugAdmin = 1,
  centerNumber = -1,
  centerName = "",
  visitNumber = -1,
  intervaldose = NULL,
  route = NULL,
  DoseTimeSheet = "DosingTimeTable",
  DoseTimeRange = NULL,
  DoseTimeMapping,
  DoseSheet = "DrugTreatmentTable",
  DoseRange = NULL,
  DoseMapping,
  positiveQC
)
}
\arguments{
\item{dataFile}{Original TAD excel datafile}

\item{Compound1}{List containining user prefered compound name and MMV compound name. One of the two must be the one used in the originla data provided by TAD (cpd1  and cpd2 by alphabetic order)}

\item{Compound2}{List containining user prefered compound name and MMV compound name. One of the two must be the one used in the originla data provided by TAD (cpd1  and cpd2 by alphabetic order)}

\item{PKsheets}{Name of the single sheet recording PK data}

\item{PKranges}{Cells range of the PKsheet to read (eg A1:C20)}

\item{PKmapping}{Character vector mapping PKsheet headers (left) to desired name in converted data format (right) (eg c("StudyID"= "STUDY", "Concentration" = "VALUE"))}

\item{PKlloq}{Numeric vector of PK LLOQ for compound 1 and compound 2 respectively}

\item{PKlloqIdentifier}{Character vector for PKlloq identifier (eg c("BLQ", "<0.01"))}

\item{PKfactor}{Numeric vector for compound 1 and compound 2 respectively, converting factor for PK data to get desired unit (eg c(1e-3,1e-3) from ng/mL to ug/mL for both compounds) )}

\item{PDsheet}{Name of the sheet recording PD data}

\item{PDmapping}{Character vector mapping PDsheet headers (left) to desired name in converted final dataset (right) (eg c("StudyID"= "STUDY", "Parasitemia" = "VALUE"))}

\item{PDlloq}{Numeric value of PD LLOQ}

\item{PDlloqIdentifier}{Numeric value of PD LLOQ}

\item{PDname}{Character, define the name of the PD variable to be recorded in the converted final dataset}

\item{DayOfFirstDrugAdmin}{Day of first drug administration as recorded in the original dataset}

\item{centerNumber}{Number corresponding to center name, see 'list_Center()'}

\item{centerName}{Center name, see 'list_Center()'}

\item{visitNumber}{Numeric Value, visit Number. Default = -1 when not relevant}

\item{intervaldose}{Character providing dosing interval. only to provide if dosing different of once a day.Used to define TRTNAME (eg "bid")}

\item{route}{Character, route of administration. Only useful if not provided in the original dataset}

\item{DoseTimeSheet}{Name of the sheet recording dosing time of drug administration for both compounds}

\item{DoseTimeRange}{Cells range of the DoseTimeSheet to read (eg A1:C20)}

\item{DoseTimeMapping}{Character vector mapping DoseTimesheet headers (left) to desired name in converted final dataset (right) (eg c("StudyID"="STUDY", "DosingTime"="NT"))}

\item{DoseSheet}{Name of the sheet recording drugs dosing for both compounds}

\item{DoseRange}{Cells range of the DoseSheet to read (eg A1:C20)}

\item{DoseMapping}{Character vector mapping DoseSheet headers (left) to desired name in converted final dataset (right) (eg c("StudyID"="STUDY", "Route" = "ROUTE"))}

\item{positiveQC}{Character vector identifying the name of positive control (must matchupper or lower case as in original dataset)}
}
\description{
Import original SCID combination data from TAD and convert it to a single structured data format
}
\details{
Import original SCID combination data from TAD and convert it to a single structured data format
ready to be imported by IQRTool (IQRdataGENERAL) or pre-ready used for MONOLIX or NONMEM
It includes positive controls data
Uses function [readxl::read_excel]
}
\seealso{
[readxl], [plyr], [dplyr], [reshape2]

Other Data Preparation: 
\code{\link{censorGametocyte}()},
\code{\link{convertGametocytes}()},
\code{\link{convert_Gametocytes}()},
\code{\link{import_SCIDpkpdData}()},
\code{\link{import_SCIDpkpdDataCombo}()},
\code{\link{import_SCIDpkpdDataCombo_GSK}()},
\code{\link{import_SCIDpkpdDataCombo_oldTAD}()},
\code{\link{import_SCIDpkpdData_GSK}()},
\code{\link{import_SCIDpkpdData_TAD}()},
\code{\link{import_SCIDpkpdData_oldTAD}()},
\code{\link{load_AfricanPediatricMalariaPopulation2to5}()},
\code{\link{load_MalariaPopulation}()},
\code{\link{load_PediatricParameters}()}
}
\author{
Aline Fuchs MMV
}
\concept{Data Preparation}
