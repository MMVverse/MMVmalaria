% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01-ActivityMonitoring.R
\name{list_Center}
\alias{list_Center}
\title{list_Center
Loads center list}
\usage{
list_Center()
}
\value{
valid list of center names and numbers
}
\description{
Loads up-to-date center list to perform center name and number verification during data preparation, see \link{check_dataGenetalMMV_Center} and \link{check_dataGeneralMMV}.
This list should be updated for each new center, adding the center name and number.
}
\seealso{
\link{check_dataGenetalMMV_Center}, \link{check_dataGeneralMMV}

Other Data Preparation: 
\code{\link{Check_MissingDatabyNAME}()},
\code{\link{censorGametocyte}()},
\code{\link{check_dataGeneralMMV}()},
\code{\link{check_dataGeneralMMV_Center}()},
\code{\link{convertGametocytes}()},
\code{\link{convert_Gametocytes}()},
\code{\link{convert_catCovToTxt}()},
\code{\link{import_SCIDpkpdData}()},
\code{\link{import_SCIDpkpdDataCombo}()},
\code{\link{import_SCIDpkpdDataCombo_GSK}()},
\code{\link{import_SCIDpkpdDataCombo_TAD}()},
\code{\link{import_SCIDpkpdDataCombo_oldTAD}()},
\code{\link{import_SCIDpkpdData_GSK}()},
\code{\link{import_SCIDpkpdData_TAD}()},
\code{\link{import_SCIDpkpdData_oldTAD}()},
\code{\link{load_AfricanPediatricMalariaPopulation2to5}()},
\code{\link{load_MalariaPopulation}()},
\code{\link{load_PediatricParameters}()},
\code{\link{summarizeReplicates}()},
\code{\link{swapName_MMVnameToName}()},
\code{\link{swapName_NameToMMVname}()},
\code{\link{transform_dataFrame_LongToWide}()},
\code{\link{transform_dataFrame_WideToLong}()}
}
\author{
To be defined
}
\concept{Data Preparation}
