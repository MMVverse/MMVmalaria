% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-DataPreparation.R
\name{convert_Gametocytes}
\alias{convert_Gametocytes}
\title{Convert Gametocytes}
\usage{
convert_Gametocytes(
  data,
  gamName = "Parasitemia female gametocytes",
  GamFactor,
  convUnit = 1,
  oldStdCurve = FALSE,
  oldNewStdFactor = 62
)
}
\arguments{
\item{data}{Data frame or IQRdataGENERAL object with gametocyte data.}

\item{gamName}{Variable name to convert (Default: \code{"Parasitemia female gametocytes"}).}

\item{GamFactor}{Conversion factor from copies/mL to p/mL (p being e.g. female gametocytes).}

\item{convUnit}{Conversion factor from copies/ul to copies/mL (Default: 1).}

\item{oldStdCurve}{TRUE or FALSE indicating whether old standard curve was used (Default: \code{"FALSE"}).}

\item{oldNewStdFactor}{Conversion factor from old to new standard curve (Default: 62).}
}
\description{
Converts gametocyte/parasitemia variable units from laboratory (generally copies/mL) to IQRdataset format units (p/mL).
Examples of gametocyte/parasitemia variables that might need unit conversion are Parasitemia female gametocytes, Parasitemia male gametocytes and Parasitemia Trophozoite.
Conversion factor (GamFactor) from copies/mL to p/mL will depend on the gametocyte/parasitemia variable and the laboratory.
oldStdCurve and oldNewStdFactor might be needed with old data from QIMR laboratory.
}
\examples{
dat <- data.frame(NAME  = c( rep("Parasitemia female gametocytes", 4)),
                  TIME  = c( 100,  120,  200, 240),
                  VALUE = c(9400, 1500, 2100, 3700),
                  UNIT  = c(rep("copies/mL", 4)),
                  LLOQ  = c(rep(1089, 4)),
                  stringsAsFactors = FALSE)

dat <- convert_Gametocytes(data            = dat,
                           gamName         = "Parasitemia female gametocytes",
                           GamFactor       = 279.3,
                           convUnit        = 1,
                           oldStdCurve     = FALSE
)

}
\seealso{
Other Data Preparation: 
\code{\link{censorGametocyte}()},
\code{\link{convertGametocytes}()},
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
\code{\link{load_PediatricParameters}()}
}
\author{
Aline Fuchs (MMV)
}
\concept{Data Preparation}
