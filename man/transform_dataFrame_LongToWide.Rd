% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-DataPreparation.R
\name{transform_dataFrame_LongToWide}
\alias{transform_dataFrame_LongToWide}
\title{transform_dataFrame_LongToWide}
\usage{
transform_dataFrame_LongToWide(
  dataFrame_Long,
  key,
  value,
  colMaster,
  patternConst = function(x) {
     out <- paste0(".", x)
 }
)
}
\arguments{
\item{patternConst}{Default: function(x) {
   out <- paste0(".", x)
}}
}
\description{
transform_dataFrame_LongToWide
}
\seealso{
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
\code{\link{list_Center}()},
\code{\link{load_AfricanPediatricMalariaPopulation2to5}()},
\code{\link{load_MalariaPopulation}()},
\code{\link{load_PediatricParameters}()},
\code{\link{summarizeReplicates}()},
\code{\link{swapName_MMVnameToName}()},
\code{\link{swapName_NameToMMVname}()},
\code{\link{transform_dataFrame_WideToLong}()}
}
\author{
Mohammed H. Cherkaoui (MMV)
}
\concept{Data Preparation}
