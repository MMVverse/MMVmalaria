% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00-GeneralFunctions.R
\name{saveMMV}
\alias{saveMMV}
\title{saveMMV
Save R objects}
\usage{
saveMMV(list, file)
}
\arguments{
\item{list}{A character vector containing the names of objects to be saved.}

\item{file}{The name of the file where the data will be saved.}
}
\description{
\link{saveMMV} writes an external representation of R objects to the specified file.
It uses the function \link[base:save]{base::save}. If 'file' is a path, it will create the directory
if not existent unlike \link[base:save]{base::save}.
}
\seealso{
\link[base:save]{base::save}

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
\code{\link{convert_Unit}()},
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
\code{\link{setMKLthreadsMMV}()},
\code{\link{simpsonMMV}()},
\code{\link{thisFile}()},
\code{\link{trapzMMV}()}
}
\author{
Mohammed H. Cherkaoui (MMV, \email{cherkaouim@mmv.org})
}
\concept{General Functions}
