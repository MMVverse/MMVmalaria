% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00-GeneralFunctions.R
\name{convert_Unit}
\alias{convert_Unit}
\title{Convert Unit}
\usage{
convert_Unit(
  data,
  UNIT_OUT,
  colUNIT = "UNIT",
  colVALUE = "VALUE",
  conversionFile = file.path(get_MMVmalariaPath(subdir = "inst"),
    "dataLibrary/unitConversion/unitConversion.csv")
)
}
\arguments{
\item{data}{data.frame where to convert value}

\item{UNIT_OUT}{Desired unit for output}

\item{colUNIT}{Name of the column unit (Default: \code{"UNIT"})}

\item{colVALUE}{Name of the column value (Default: \code{"VALUE"})}

\item{conversionFile}{Path and name of conversion file (Default: \code{file.path(get_MMVmalariaPath(subdir = "inst"), "dataLibrary/unitConversion/unitConversion.csv")})}
}
\value{
data.frame identical to \code{data} but with the \code{colVALUE} adjust to the desired output unit as defined in \code{UNIT_OUT}.
}
\description{
Automate the conversion of unit in a data frame.
}
\examples{
dat <- data.frame(TIME  = c( 0, 1, 2, 4, 8),
                  VALUE = c(10, 6, 4, 3, 2),
                  UNIT  = "ng/mL",
                  stringsAsFactors = FALSE)
 dat <- convert_Unit(data     = dat,
                     UNIT_OUT = c("ug/mL"))

}
\seealso{
Other General Functions: 
\code{\link{IQRtableToDataFrame}()},
\code{\link{add_FuncToEnv}()},
\code{\link{adjust_Dose}()},
\code{\link{aux_CommonSubPath}()},
\code{\link{aux_addNamesToPercentiles}()},
\code{\link{aux_constructCIpercentiles}()},
\code{\link{aux_createUSUBJID}()},
\code{\link{aux_formatErrorName}()},
\code{\link{aux_removeEscapeChar}()},
\code{\link{cbindMMV}()},
\code{\link{clopperPearsonMMV}()},
\code{\link{combMMV}()},
\code{\link{create_PKPDsimtime}()},
\code{\link{file.copyMMV}()},
\code{\link{find_MinMMV}()},
\code{\link{getDoselevel}()},
\code{\link{getMKLthreadsMMV}()},
\code{\link{get_IQRdistribution}()},
\code{\link{get_IXGDFtoRemove}()},
\code{\link{get_MMVmalariaPath}()},
\code{\link{libraryMMV}()},
\code{\link{load_RData}()},
\code{\link{logLinTrapzMMV}()},
\code{\link{newtonRaphson.Function}()},
\code{\link{newtonRaphson.Vector}()},
\code{\link{rbindMMV}()},
\code{\link{rectintMMV}()},
\code{\link{reverse_List}()},
\code{\link{saveActivityInfo}()},
\code{\link{saveMMV}()},
\code{\link{setMKLthreadsMMV}()},
\code{\link{simpsonMMV}()},
\code{\link{thisFile}()},
\code{\link{trapzMMV}()}
}
\author{
Mohammed H. Cherkaoui (MMV, \email{cherkaouim@mmv.org})
}
\concept{General Functions}
