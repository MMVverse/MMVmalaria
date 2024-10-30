#' Adds functions to an environment
#'
#' Extracts the functions of the environment \code{e2} and add them to the initial
#' environment \code{e1}.
#'
#' @param e1 Initial environment
#' @param e2 Environment whose functions are to be extracted and added to the initial environment
#'
#' @return Environment with the functions of both input environments
#'
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
add_FuncToEnv <- function(e1, e2){
  # Env. Name:
  e1name <- deparse(substitute(e1))
  e2name <- deparse(substitute(e2))

  # Object in each environment:
  listE1 <- ls(e1)
  listE2 <- ls(e2)

  # Add all function in e2 in e1:
  for(v in listE2){
    if ("function" %in% class(e2[[v]])){
      if(v %in% listE1){
        warning("Function '", v, "' is in '", e1name, "' therefore, not added to '", e2name, "'!")
      }else{
        e1[[v]] <- e2[[v]]
      }
    }
  }

  # Output:
  e1
}


#' Adjust Doses
#'
#' When dose predictions are done, it is often the case that the predicted doses have
#' decimal number. As this will never been given to patient, this function take the ceil
#' number of \code{Dose}.
#'
#' @param Dose Numeric vector with doses to adjust
#'
#' @return Numeric vector of same size as \code{Dose}
#'
#' @examples
#' Dose <- c(0.12, 1.2, 87.5, 487)
#' adjust_Dose(Dose = Dose)
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
adjust_Dose <- function(Dose){

  # Check for negative doses:
  if (any(Dose<0)){
    stop("Negative Dose: Please Adjsut")
  }

  # Adjust Dose:
  Dose <- dplyr::case_when((Dose<1)   ~ ceiling(Dose*10)/10,
                           (Dose<10)  ~ ceiling(Dose)      ,
                           (Dose<100) ~ ceiling(Dose/5 )* 5,
                           TRUE       ~ ceiling(Dose/10)*10)

  # Output:
  Dose
}


#' Add names to percentiles vector
#'
#' @param percentiles Vector with percentiles values.
#'
#' @return The vector `percentiles` with names.
#'
#' @author Mohammed H. Cherkaoui (MMV)
#' @family General Functions
aux_addNamesToPercentiles <- function(percentiles){

  # Check that it is between 0 and 100:
  if(any(percentiles>100) || any(percentiles<0)){
    stop("'percentiles' should be a vector containings only values between 0 and 100: Please adjust.")
  }

  # Add Names:
  percentiles        <- sort(percentiles)
  names(percentiles) <- ifelse(percentiles==50,
                               "Median",
                               paste0(percentiles, "th Percentile"))

  # OUtput:
  percentiles
}


#' Construct CI percentiles from CI level
#'
#' @param CIlevel Numeric containing the confidence interval in percent to be estimated (Default: 90)
#'
#' @return The vector \code{CI.percentiles} with names.
#'
#' @author Mohammed H. Cherkaoui (MMV)
#' @family General Functions
aux_constructCIpercentiles <- function(CIlevel){

  # Define percentiles of CI:
  if(CIlevel<0 || CIlevel>100){
    stop("'CIlevel' should be between 0 and 100: Please adjust.")
  }else if(CIlevel==0){
    CI.percentile        <- c(50)
    names(CI.percentile) <- c("CI Median")
  }else{
    CI.percentile        <- c(50 - CIlevel/2, 50         , 50 + CIlevel/2)
    names(CI.percentile) <- c("CI Low"      , "CI Median", "CI High"     )
  }

  # Output:
  CI.percentile
}


#' Common Sub-Path
#'
#' Identify the common sub-path in the vector of paths \code{x}
#'
#' @param x Vector of paths for which the sub-path needs to be detected
#'
#' @return Character of the common sub-path
#'
#' @examples
#' Paths <- c("C:/Users/Default",
#'            "C:/Users/Public")
#' aux_CommonSubPath(Paths)
#'
#' @export
#'
#' @author Anne Kuemmel (IntiQuan), Nathalie Gobeau (MMV, \email{gobeaun@@mmv.org})
#' @family General Functions
aux_CommonSubPath <- function(x) {
  # Sort the vector:
  x <- sort(x)

  # Split the first and last element by path separator:
  d_x <- strsplit(x[c(1,length(x))], "/")

  # Search for the first non common element to get the last matching one:
  der_com <- match(FALSE, do.call("==",d_x)) - 1

  # If there is no matching element, return an empty vector, else return the common part:
  out <- NULL
  if(der_com==0){
    out <- character(0)
  }else{
    out <- paste0(d_x[[1]][1:der_com], collapse = "/")
  }

  # Output:
  out
}


#' Create USUBJID
#'
#' Add column \code{USUBJID} to \code{data} based on columns \code{STUDY}, \code{GROUP} and \code{SUBJECT},
#' when preparing IQR data from SCID dataset
#'
#' @param data Dataset for which to create \code{USUBJID}
#'
#' @return \code{data} with extra column \code{USUBJID}
#'
#' @examples
#' data <- data.frame(STUDY   = "StudyID",
#'                    GROUP   = c("G1","G1","G2","G2"),
#'                    SUBJECT = c("M1","M2", "M1", "M2"),
#'                    stringsAsFactors = FALSE)
#' data$USUBJID <- MMVmalaria:::aux_createUSUBJID(data)
#'
#' @author Aline Fuchs (MMV), Anne Kuemmel (IntiQuan)
#' @family General Functions
aux_createUSUBJID <- function(data) {
  with(data, {

    # Format group with 2 digits
    if (is.character(GROUP)) {
      GROUP <- ifelse(nchar(GROUP) == 1, sprintf("0%s", GROUP), sprintf("%s", GROUP))
    } else {
      GROUP <- sprintf("%02i",GROUP)
    }

    # Format Subject with minimum 2 digits
    if (is.character(SUBJECT)) {
      SUBJECT <- ifelse(nchar(SUBJECT) == 1, sprintf("0%s", SUBJECT), sprintf("%s", SUBJECT))
    } else {
      SUBJECT <- sprintf("%02i",SUBJECT)
    }

    # Create subject ID nbr
    paste(STUDY,GROUP,SUBJECT,sep="_")
  } )
}


#' Format Error Name
#'
#' Converts Error Parameter Names as used in sysFIT to format as used in NLME
#'
#' @param parameterNames Vector of character with the names of the parameters
#'
#' @return \code{parameterNames} with converted format
#'
#' @author Anne Kuemmel (IntiQuan)
#' @family General Functions
aux_formatErrorName <- function(parameterNames) {
  # Converts error parameter names as used in sysfit to format as used in nlme
  idxErr <- grep("OUTPUT", parameterNames)
  for (k in idxErr) {
    parameterNames[k] <- paste0(parameterNames[k], substr(parameterNames[k],7,7))
  }
  parameterNames <- gsub("_sigma_abs", "error_ADD", parameterNames)
  parameterNames <- gsub("_sigma_rel", "error_PROP", parameterNames)
  parameterNames <- gsub("OUTPUT[[:digit:]]", "", parameterNames)

  # Output:
  parameterNames
}


#' Removal of escape characters
#'
#' @param x A character vector.
#' @param escapeChars Escape character to substituted with space (Default: \code{c("\\r\\n", "\\n", "\\r", "\\t")})
#'
#' @return \code{x} where all escape character defined in \code{escapeChars} were replaced by a space.
#'
#' @author Anne Kuemmel (IntiQuan)
#' @family General Functions
aux_removeEscapeChar <- function(x, escapeChars = c("\r\n","\n","\r","\t")) {
  # Removal of escape characters
  # Under windows anyway substituted by white space, but not for linux
  for (ek in escapeChars) {
    x <- gsub(ek," ", x, fixed = TRUE)
  }

  # Output:
  x
}


#' Clopper-Pearson Test
#'
#' Statistical test to estimate the confidence interval of
#' a binomial distribution using the Clopper-Pearson Test
#' (see \href{https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval}{Wikipedia})
#'
#' @param k Number of successful observations of a trial
#' @param n Number of observation
#' @param CI Confidence Interval to estimate (Default: \code{0.9})
#'
#' @return Return a vector with the lower \code{CI_low} and high \code{CI_high} limit of the CI.
#'
#' @examples
#' # Probability of success:
#' p <- 0.8
#' # Number of observations:
#' n <- 100
#' # Number of successful observation:
#' k <- rbinom(1,n,p)
#' # Empirical probability of success:
#' p_Empirical <- k/n
#' # 90% confidence interval of the empirical probability of success:
#' CI <- clopperPearsonMMV(k,n,CI=0.9)
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
clopperPearsonMMV <- function(k,
                              n,
                              CI = 0.90) {

  #-----------------------------------------------------------------------------#
  # STEP 0: Some Checks ----
  #-----------------------------------------------------------------------------#

  # Check if arguments are numerics:
  if(!is.numeric(k) || !is.numeric(n) || !is.numeric(CI)){
    stop("'k', 'n' and/or 'CI' are not numeric: Please adjust.")
  }

  # Check if CI is larger than one:
  if(CI>1){
    CI <- CI/100
    warning("'CI' is larger than 1: It was assumed to be expressed in percent and, therfore, convert to a fraction.")
  }


  #-----------------------------------------------------------------------------#
  # STEP 1: Estimate CI ----
  #-----------------------------------------------------------------------------#

  # Estimate low and high limits using the beta function:
  limits_low <- qbeta(p      = (1-CI)/2,
                      shape1 = k,
                      shape2 = n - k +1)

  limits_high <- qbeta(p      = 1 - (1-CI)/2,
                       shape1 = k + 1,
                       shape2 = n - k)

  # Generate output:
  limits <- c(limits_low, limits_high)
  names(limits) <- c("CI_low", "CI_high")

  # Output:
  limits
}


#' Customizable Collapse Function
#'
#' Collapses a character vector of any length into a character string separated by
#' \code{collapseSymbole}, except for the last element where \code{andSymbole} is used.
#'
#' @param x Character vector of length \code{n}.
#' @param collapseSymbole Character used as a separator for the \code{n-1} first elements in \code{x} (Default: \code{", "})
#' @param andSymbole Character used as a separator for the last element in \code{x} (Default: \code{" & "}).
#' @param messageEmpty Character to print if \code{x} is empty (Default: \code{NULL}).
#'
#' @return A character string where the elements of \code{x} were collapsed.
#'
#' @examples
#' PMXmember <- c("Aline", "Catalina", "Mohammed", "Nathalie")
#' msg       <- paste0("The PMX group is composed of ", collapseMMV(PMXmember), ".")
#' cat(msg)
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions, Editing Functions
collapseMMV <- function(x,
                        collapseSymbole = ", ",
                        andSymbole      = " & ",
                        messageEmpty    = NULL){
  text <- character(0)
  if(length(x)==0){
    warning("'x' is of length 0 in the function 'collapseMMV'")
    text <- paste0(messageEmpty)

  }else if(length(x)==1){
    text <- x[1]

  }else{
    n_x  <- length(x)
    text <- paste0(paste0(x[1:(n_x-1)], collapse = collapseSymbole),
                   andSymbole, x[n_x])
  }

  # Output:
  text
}


#' Convert Unit
#'
#' Automate the conversion of unit in a data frame.
#'
#' @param data data.frame where to convert value
#' @param UNIT_OUT Desired unit for output
#' @param colUNIT Name of the column unit (Default: \code{"UNIT"})
#' @param colVALUE Name of the column value (Default: \code{"VALUE"})
#' @param conversionFile Path and name of conversion file (Default: \code{file.path(get_MMVmalariaPath(subdir = "inst"), "dataLibrary/unitConversion/unitConversion.csv")})
#'
#' @return data.frame identical to \code{data} but with the \code{colVALUE} adjust to the desired output unit as defined in \code{UNIT_OUT}.
#'
#' @examples
#' dat <- data.frame(TIME  = c( 0, 1, 2, 4, 8),
#'                   VALUE = c(10, 6, 4, 3, 2),
#'                   UNIT  = "ng/mL",
#'                   stringsAsFactors = FALSE)
#'  dat <- convert_Unit(data     = dat,
#'                      UNIT_OUT = c("ug/mL"))
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
convert_Unit <- function(data,
                         UNIT_OUT,
                         colUNIT        = "UNIT",
                         colVALUE       = "VALUE",
                         conversionFile = file.path(get_MMVmalariaPath(subdir = "inst"), "dataLibrary/unitConversion/unitConversion.csv")
){

  # Load converison CSV file:
  conversionTable <- IQRtools::IQRloadCSVdata(filename = conversionFile)

  FLAGfactor <- FALSE
  if(is.factor(data[[colUNIT]])){
    FLAGfactor <-  TRUE
    data[[colUNIT]] <-  as.character(data[[colUNIT]])
  }

  # Loop on Unit Out:
  for(UNIT_OUT_k in UNIT_OUT){
    conversionTable_k <- conversionTable[conversionTable[["UNIT_OUT"]] == UNIT_OUT_k,]
    if(nrow(conversionTable_k)==0){
      message("Unit '", UNIT_OUT_k, "' not defined in 'conversionFile': Please adjust your conversion CSV file")

    }else{
      for(UNIT_IN_j in unique(conversionTable_k[["UNIT_IN"]])){
        conversionTable_kj <- conversionTable_k[conversionTable_k[["UNIT_IN"]] == UNIT_IN_j,]
        if(nrow(conversionTable_kj)==0){
          message("Unit '", UNIT_IN_j, "' not defined with unit '",UNIT_OUT_k,"' in 'conversionFile': Please adjust your conversion CSV file")
        }else if(nrow(conversionTable_kj)>1){
          stop("Unit '", UNIT_IN_j, "' is defined multiple times with '",UNIT_OUT_k,"' in 'conversionFile': Please adjust your conversion CSV file")
        }else{
          idx_kj <- data[[colUNIT]] == UNIT_IN_j
          data[[colUNIT]][idx_kj]   <- UNIT_OUT_k

          Formula_kj <- conversionTable_kj[["FORMULA"]]
          Formula_kj <- gsub(" "   , ""                          , Formula_kj, fixed = TRUE)
          Formula_kj <- gsub("X_out="   , ""                     , Formula_kj, fixed = TRUE)
          Formula_kj <- gsub("X_in", "data[[colVALUE]][idx_kj]"  , Formula_kj, fixed = TRUE)
          data[[colVALUE]][idx_kj] <- eval(parse(text=Formula_kj))
        }
      }
    }
  }
  if(FLAGfactor){
    data[[colUNIT]] <-  as.factor(data[[colUNIT]])
  }

  # Output:
  data
}


#' Generate Simulation Times
#'
#' Calculate time points using given time steps between specifed time point. This is useful
#' to optimize time step for ploting. For example with PK profile the time step is often needed
#' to be small during absorption and large during the elimination phase.
#'
#' Calculate timepoints from 0 to Tend separated by regular time steps sepcified in dt. Several time steps can be specified: they will be applied
#' for the respective time intervals specified in the vector Tswitch.
#'
#' @param Tend Final time
#' @param dt Time steps (Defaul: \code{c(  0.25,  1,  2, 6)}).
#' @param Tswitch Time intervals over which time steps are used. It should be a vector of similar length as \code{dt} (Default: \code{c(0,  12, 24, 48)}).
#'
#' @return Numeric vector with the calcultaed time point
#'
#' @examples
#' create_PKPDsimtime(Tend=1,dt=c(0.1,0.5),Tswitch=c(0,0.4))
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
create_PKPDsimtime <- function(Tend,
                               dt      = c(  0.25,  1,  2, 6),
                               Tswitch = c(0,  12, 24, 48)){

  # Control that Tend is larger than first value of Tswitch
  if (Tend<=Tswitch[1]){
    stop("First value in 'Tswitch' is larger than 'Tend'. Please adjust 'Tswitch'.")
  }

  # Control that dt and Tswitch have same length:
  if (length(dt)!=length(Tswitch)){
    stop("'dt' and 'Tswitch' have different length. Please adjust.")
  }

  # Generate simtim:
  #   Loop over dt/Tswitch
  simtime <- c()
  k       <- 1
  while (Tend>Tswitch[k+1] & k<length(dt)){
    simtime_k <- seq(Tswitch[k], Tswitch[k+1], by=dt[k])
    simtime   <- c(simtime, simtime_k)
    k         <- k+1
  }
  #   Add last segment
  simtime_k <- seq(Tswitch[k], Tend, by=dt[k])
  simtime   <- unique(c(simtime, simtime_k, Tend))

  # Output:
  simtime
}




#' File Manipulation
#'
#' [file.copyMMV] works in a similar way to [base::file.copy] but automatically create
#' the destination folder if needed. Copying to existing destination files is skipped
#' unless \code{overwrite = TRUE}. The to argument can specify a single existing directory.
#' If \code{copy.mode = TRUE} file read/write/execute permissions are copied where possible,
#' restricted by \code{'umask'}. (On Windows this applies only to files.) Other security
#' attributes such as ACLs are not copied. On a POSIX filesystem the targets of symbolic
#' links will be copied rather than the links themselves, and hard links are copied separately.
#' Using \code{copy.date = TRUE} may or may not copy the timestamp exactly (for example, fractional
#' seconds may be omitted), but is more likely to do so as from R 3.4.0.
#'
#' @param from Character vectors, containing file names or paths.
#' @param to Character vectors, containing file names or paths.
#' @param overwrite Logical; should existing destination files be overwritten?
#' @param recursive Logical. If to is a directory, should directories in from be copied (and their contents)? (Like cp -R on POSIX OSes.)
#' @param copy.mode Logical: should file permission bits be copied where possible?
#' @param copy.date Logical: should file dates be preserved where possible? See Sys.setFileTime.
#'
#' @md
#'
#' @export
#' @seealso [base::file.copy]
#' @family General Functions
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
file.copyMMV <- function(from,
                         to,
                         overwrite = recursive,
                         recursive = FALSE,
                         copy.mode = TRUE,
                         copy.date = FALSE){

  # Get directory name of output:
  if (is.fileMMV(from)){
    dirName <- dirname(to)
  } else{
    dirName <- to
  }

  # Check if the directory name exists, otherwise create it:
  if (!dir.exists(dirName)){
    dir.create(dirName, recursive = TRUE)
  }

  # Copy file(s):
  if (is.fileMMV(from)){
    file.copy(from = from,
              to   = to,
              overwrite = overwrite,
              recursive = recursive,
              copy.mode = copy.mode,
              copy.date = copy.date)
  } else{
    sapply(as.vector(dir(from)),
           FUN = function(from_k){
             dirName_k <- ifelse(dirname(from_k)==".",
                                 "",
                                 dirname(from_k))
             file.copy(from = file.path(from, from_k),
                       to   = file.path(to  , dirName_k),
                       overwrite = overwrite,
                       recursive = recursive,
                       copy.mode = copy.mode,
                       copy.date = copy.date)
           })
  }

}


#' Find Minimum
#'
#' Finds the minimum of \code{Y(X)} by interpolation with one of the 3 methods: no interpolation, cubic spline or quadratic.
#' Returns the minimum value of \code{Y} (\code{Ymin}) and the value of \code{X} when \code{Y} is minimum (\code{Xmin})
#'
#' @param X numeric vector
#' @param Y numeric vector, same length as X
#' @param dx numeric value, refined interpolation interval in case of the cubic spline interpolation
#' @param Method interpolation method to choose from  \code{"NoInterpolation"}, \code{"CubicSpline"} or \code{"Quadratic"} (Default: \code{"Quadratic"})
#'
#' @return Data frame with \code{Xmin} anf \code{Ymin}: \code{data.frame(Xmin=Xmin,Ymin=Ymin)}
#'
#' @examples
#' X  = c(1,2,3,4,5,6,7)
#' Y  = c(10,4,1,1,6,10,20)
#' Y1 = find_MinMMV(X,Y,Method = "NoInterpolation")
#' Y2 = find_MinMMV(X,Y,Method = "CubicSpline")
#' Y3 = find_MinMMV(X,Y,Method = "Quadratic")
#' plot(X,Y,type="b",ylim=c(0,20))
#' points(Y1,pch=16,col="blue")
#' points(Y2,pch=16,col="green")
#' points(Y3,pch=16,col="red")
#' legend("topleft",legend=c("NoInterpolation","CubicSpline","Quadratic"),pch=rep(16,3),col=c("blue","green","red"))
#'
#' @export
#' @family General Functions
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
find_MinMMV <- function(X,
                        Y,
                        dx = as.numeric((X[2]-X[1])/1000),
                        Method = "Quadratic")
{

  #-----------------------------------------------------------------------------#
  # STEP 1: Check Input Variables ----
  #-----------------------------------------------------------------------------#

  # Check if X and Y have same length:
  if (length(X)!=length(Y)){
    stop("'X' and 'Y' have different lengths. Please adjust the vectors.")
  }

  # Check if Method is valid:
  if (!(tolower(Method) %in% c("nointerpolation", "cubicspline", "quadratic"))){
    stop("'Method' is not a valid method. Please choose between 'NoInterpolation', 'CubicSpline', 'Quadratic', or implement a new method.")
  }


  #-----------------------------------------------------------------------------#
  # STEP 2: Find MIN ----
  #-----------------------------------------------------------------------------#

  # Get index of the minimum value in Y:
  idxMIN <- which.min(Y)

  # If no interpolation:
  if (tolower(Method)=="nointerpolation"){
    Xmin   <- X[idxMIN]
    Ymin   <- Y[idxMIN]


    # Cubic Spline:
  } else if(tolower(Method)=="cubicspline"){

    # If MIN is first or last point:
    if (idxMIN==1 | idxMIN==length(Y)) {
      Xmin   <- X[idxMIN]
      Ymin   <- Y[idxMIN]

      # Otherwise:
    } else {

      # Use 5 or 3 points depending on where idxMIN is:
      if ((idxMIN==length(Y)-1) | (idxMIN==2)){
        # (x;y) points to use for the spline:
        X0 <- X[c(idxMIN-1, idxMIN, idxMIN+1)]
        Y0 <- Y[c(idxMIN-1, idxMIN, idxMIN+1)]

      }else{
        # (x;y) points to use for the spline:
        X0 <- X[c(idxMIN-2, idxMIN-1, idxMIN, idxMIN+1, idxMIN+2)]
        Y0 <- Y[c(idxMIN-2, idxMIN-1, idxMIN, idxMIN+1, idxMIN+2)]
      }

      # X points to use for the spline interpolation:
      Xs <- seq(X[idxMIN-1], X[idxMIN+1], dx)

      # Spline Cubic:
      Ys <- pracma::cubicspline(X0, Y0, Xs)

      # Find minimum:
      idxMINs <- which.min(Ys)
      Xmin   <- Xs[idxMINs]
      Ymin   <- Ys[idxMINs]
    }


    # Quadratic:
  } else if(tolower(Method)=="quadratic"){


    # If MIN is first or last point:
    if (idxMIN==1 | idxMIN==length(Y)) {
      Xmin   <- X[idxMIN]
      Ymin   <- Y[idxMIN]

      # Otherwise:
    } else{
      # FIRST ---
      # Interpolation order 2 to find the "real" minimum:

      # Points for the regression:
      #   Parasitemia before the observed minimum
      X1 <- X[idxMIN-1]
      Y1 <- Y[idxMIN-1]
      #   Parasitemia at observed minimum
      X2 <- X[idxMIN]
      Y2 <- Y[idxMIN]
      #   Parasitemia after the observed minimum
      X3 <- X[idxMIN+1]
      Y3 <- Y[idxMIN+1]

      # Constants of interest to estimate the polynome constants:
      A1 <- Y1/((X1-X2)*(X1-X3))
      A2 <- Y2/((X2-X1)*(X2-X3))
      A3 <- Y3/((X3-X1)*(X3-X2))

      # Constant of the polynom:
      a <- A1 + A2 + A3
      b <- -(A1*(X2+X3) + A2*(X1+X3) + A3*(X1+X2))
      c <- A1*X2*X3 + A2*X1*X3 + A3*X1*X2

      # Minimum Time:
      Xmin <- -b/(2*a)
      Ymin <- a*Xmin^2 + b*Xmin + c
    }
  }


  #Genereate Ouput:
  PointMIN <- data.frame(Xmin = Xmin,
                         Ymin = Ymin,
                         stringsAsFactors = TRUE)

  # Return:
  return(PointMIN)
}


#' Get Activity Path
#'
#' Get activity path to be used in graphs. This allows to keep track of where the
#' figures are located and which script was used to generate them.
#'
#' @param ActivityPath If the user want to force the activity path to plot. (Default: \code{NULL}).
#'
#' @return A character vector with the activity path
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
get_ActivityPath <- function(ActivityPath = NULL) {

  # Define Activity Path for caption:
  if (is.null(ActivityPath)){
    # ActivityPath if on PiNK
    if(file.exists("../tags.RData")){
      load("../tags.RData")
      ActivityPath <- file.path(gsub("/sites/department/ModellingTeam/","",tags$itemPath),
                                tags$itemName)

    # ActivityPath if on S:/M&S
    }else{
      # Get Current Working Directory:
      WorkDir       <- getwd()
      WorkDir_Split <- strsplit(WorkDir, "/", fixed = TRUE)[[1]]
      n_str         <- length(WorkDir_Split)

      if (n_str>5 && WorkDir_Split[n_str-4]!="Projects" && WorkDir_Split[n_str-5]=="Projects_Discovery"){

        # Define Serie:
        SerieName  <- WorkDir_Split[n_str-4]

        # Define Project Name:
        ProjectName <- WorkDir_Split[n_str-2]

        # Define ActivityName if NULL:
        ActivityName  <- WorkDir_Split[n_str-1]

        # Define Activity Path:
        ActivityPath <- file.path(SerieName,"Work", ProjectName,  ActivityName)

      }else if(n_str>4 && WorkDir_Split[n_str-4]=="Projects"){
        # Define Project Name:
        ProjectName <- WorkDir_Split[n_str-3]

        # Define ActivityName if NULL:
        ActivityName  <- WorkDir_Split[n_str-1]

        # Define Activity Path:
        ActivityPath <- file.path(ProjectName, "Work", ActivityName)
      } else{
        warning("The current activity is not in the project folder and 'ActivityPath' is set to NULL, therefore, it is not automatically generated\nYou can manually enter a value for 'ActivityPath'")
        ActivityPath <- "NULL"
      }
    }
  }

  # Output:
  ActivityPath
}



#---------------------#
# DOCUMENTATION ----
#---------------------#



#' Get Distribution
#'
#' Get Distribution of list of parameters from an IQR model
#'
#' @param Distribution Character vector containing the transform formula from an IQR model
#'
#' @return Return a charaacter vector containing the distribution index; i.e. "N" for Normal, "L" for Log-Normal and "G" for Logit
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
get_IQRdistribution <- function(Distribution) {
  Distribution <- gsub("log(psi/(1-psi))","G",Distribution,fixed = TRUE)
  Distribution <- gsub("log(psi)"  ,"L",Distribution,fixed = TRUE)
  Distribution <- gsub("(psi)"     ,"N",Distribution,fixed = TRUE)

  return(Distribution)
}
#' Get Row ID to be removed
#'
#' _IXGDFtoRemove
#'
#' @description
#' @param dataGeneral A general dataset to look into
#' @param obNAME Default: 'log(Parasitemia Total)'
#' @param removeType Default: 'ToLastOb'
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
get_IXGDFtoRemove <- function(dataGeneral,
                              obNAME     = "log(Parasitemia Total)",
                              removeType = "ToLastOb") {

  # Reduce dataset to observation of interest:
  dataRegress <- dataGeneral[dataGeneral$NAME==obNAME,]

  # Do the right operations:
  if (toupper(removeType)=="TOLASTOB"){

    # First detect time of last observation:
    idx_Above <- which(dataRegress$VALUE>dataRegress$LLOQ)
    Tmax      <- max(dataRegress$TIME[idx_Above])

    # Get all IXGDF to discard:
    idx_ToRemove <- sapply(unique(dataRegress$USUBJID),
                           function(ID){
                             # Reduce dataset:
                             dataRegress_ID <- dataRegress[dataRegress$USUBJID==ID,]

                             # Discard observations after Tmax:
                             idx_ToRemove_ID <- which(dataRegress_ID$TIME>Tmax)
                             if (length(idx_ToRemove_ID)>0){
                               out <- dataRegress_ID$IXGDF[idx_ToRemove_ID]
                             }else{
                               out <- NULL
                             }

                             # Output:
                             return(out)
                           })

    # Unlist:
    idx_ToRemove        <- unlist(idx_ToRemove)
    names(idx_ToRemove) <- NULL


  }else if(toupper(removeType)=="REMOVECURE"){

    # Get all USUBJID:
    idx_ToRemove <- sapply(unique(dataRegress$USUBJID),
                           function(ID){
                             # Reduce dataset:
                             dataRegress_ID <- dataRegress[dataRegress$USUBJID==ID,]

                             # Get index of maximum time:
                             idx_Tmax <- which.max(dataRegress_ID$TIME)

                             # Generate output:
                             if (dataRegress_ID$VALUE[idx_Tmax]<=dataRegress_ID$LLOQ[idx_Tmax]){
                               # Get all indexes for which parasitemia was above LLOQ
                               idx_Above <- which(dataRegress_ID$VALUE>dataRegress_ID$LLOQ)

                               # Get the last time for which a parasitemia was observed:
                               Tlast <- max(dataRegress_ID$TIME[idx_Above])

                               # Get IXGDF between Tlast and Tmax:
                               out <- dataRegress_ID$IXGDF[dataRegress_ID$TIME>Tlast]
                             } else{
                               out <- NULL
                             }

                             # Output:
                             return(out)
                           })

    # Unlist:
    idx_ToRemove        <- unlist(idx_ToRemove)
    names(idx_ToRemove) <- NULL


  }else if(toupper(removeType)=="CURE1WEEK"){

    # Get all USUBJID:
    idx_ToRemove <- sapply(unique(dataRegress$USUBJID),
                           function(ID){
                             # Reduce dataset:
                             dataRegress_ID <- dataRegress[dataRegress$USUBJID==ID,]

                             # Get index of maximum time:
                             idx_Tmax <- which.max(dataRegress_ID$TIME)

                             # Generate output:
                             if (dataRegress_ID$VALUE[idx_Tmax]<=dataRegress_ID$LLOQ[idx_Tmax]){
                               # Get all indexes for which parasitemia was above LLOQ
                               idx_Above <- which(dataRegress_ID$VALUE>dataRegress_ID$LLOQ)

                               # Get the last time for which a parasitemia was observed:
                               Tlast <- max(dataRegress_ID$TIME[idx_Above])

                               # Get IXGDF between Tlast and Tmax:
                               out <- dataRegress_ID$IXGDF[dataRegress_ID$TIME>(Tlast+24*7)]
                             } else{
                               out <- NULL
                             }

                             # Output:
                             return(out)
                           })

    # Unlist:
    idx_ToRemove        <- unlist(idx_ToRemove)
    names(idx_ToRemove) <- NULL


  }else{
    idx_ToRemove <- NULL
  }


  # Output:
  return(idx_ToRemove)
}

#' getDoselevel
#'
#' @description
#' @param chr
#' @param N Default: 1
#' @return
#' @export
#' @author Anne Kuemmel (IntiQuan)
#' @family General Functions
getDoselevel <- function(chr, N = 1) {
  out <- strsplit(chr, "+", fixed = TRUE)
  out <- sapply(out, function(x) x[N])
  out <- strsplit(out, "mg", fixed = TRUE)
  out <- sapply(out, function(x) x[1])
  out <- pmax(1e-12, as.numeric(out))
  out
}
#' getMKLthreadsMMV
#' Get the Number of Threads for Math Kernel Library
#'
#' Get the maximum number of threads that can be started by the Math Kernel Library BLAS.
#' This should be less than or equal to the number of processing cores on your computer.
#' It is based on [getMKLthreads] for the Microsoft package 'RevoUtilsMath'. It first
#' checks if 'RevoUtilsMath' is loaded (As it is absent when 'R' is used instead of 'MRO').
#'
#' @return The maximum number of threads that can be started by the Math Kernel Library BLAS.
#'
#' @md
#'
#' @family General Functions
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
getMKLthreadsMMV <- function() {

  # Information about the current session:
  sessInfo <- sessionInfo()

  # Get MKL thread:
  if ("RevoUtilsMath" %in% c(sessInfo$basePkgs, names(sessInfo$otherPkgs))){
    n_Threads <- getMKLthreads()
  }else{
    n_Threads <- 1
  }

  return(n_Threads)
}


#' Convert saved IQR table back to data.frame
#'
#' @param inputFileName Path to text file saved with `IQRtools::IQRoutputTable` to convert back to data.frame
#'
#' @return
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
IQRtableToDataFrame <- function(inputFileName) {

  # This function allows to re-generate a dataFrame, title and footer from a table generated
  # using IQRoutputTable.
  #
  # This can be useful, if a table needs to be changed, for example.

  # Load table:
  rawTable <- readLines(inputFileName)

  # Retrieve the information of each line:
  Title   <- FALSE
  Header  <- TRUE
  xTitle  <- NULL
  Footer  <- FALSE
  xFooter <- NULL
  for (k in 1:length(rawTable)){
    # Check if the first line is a title:
    if (k==1 & !grepl("|", rawTable[k], fixed=TRUE, useBytes=TRUE)){
      xTitle <- stringr::str_trim(rawTable[k])
    }

    # Retrieve Data Frame:
    if (grepl("|", rawTable[k], fixed=TRUE, useBytes=TRUE)){
      #   First Header:
      if (Header){
        # Get columns name:
        ColumnsName <- strsplit(rawTable[k], "|", fixed=TRUE, useBytes=TRUE)[[1]]
        ColumnsName <- stringr::str_trim(ColumnsName)

        # Create dataFrame to return:
        dataFrame <- data.frame(matrix(vector(), nrow=0, ncol=length(ColumnsName)),
                                stringsAsFactors = FALSE)
        colnames(dataFrame) <- ColumnsName

        # The rest of the table is not header:
        Header <- FALSE

        #   Then the rest of the data frame:
      }else{
        # Retrieve next line:
        rawTable_k <- strsplit(rawTable[k], "|", fixed=TRUE, useBytes=TRUE)[[1]]
        rawTable_k <- stringr::str_trim(rawTable_k)
        rawTable_k <- as.data.frame(t(rawTable_k), stringsAsFactors = FALSE)
        colnames(rawTable_k) <- ColumnsName

        # Add it to dataFrame:
        dataFrame <- rbind(dataFrame,rawTable_k)
      }


      # Retrieve footer if there is one:
    } else if(!grepl("----", rawTable[k], fixed=TRUE) & !grepl("====", rawTable[k], fixed=TRUE) & rawTable[k]!="" & k!=1){
      if (!Footer){
        xFooter <- rawTable[k]
        Footer  <- TRUE
      } else{
        xFooter <- paste0(xFooter,"<br>",rawTable[k])
      }
    }
  }

  # Adapt First column name, title and footer if needed:
  xTitle           <- gsub("<TT>   ", "", xTitle          , fixed=TRUE, useBytes=TRUE)
  names(dataFrame) <- gsub("<TH>   ", "", names(dataFrame), fixed=TRUE, useBytes=TRUE)
  dataFrame[,1]    <- gsub("<TR>   ", "", dataFrame[,1]   , fixed=TRUE, useBytes=TRUE)
  dataFrame[,1]    <- gsub("<TR>"   , "", dataFrame[,1]   , fixed=TRUE, useBytes=TRUE)
  xFooter          <- gsub("<TF>   ", "", xFooter         , fixed=TRUE, useBytes=TRUE)

  # Output:
  return(list(dataFrame=dataFrame, xTitle=xTitle, xFooter=xFooter))
}
#' libraryMMV
#'
#' @description
#' @param package
#' @param help (Default: \code{NULL}).
#' @param pos Default: 2
#' @param lib.loc (Default: \code{NULL}).
#' @param logical.return Default: FALSE
#' @param warn.conflicts Default: TRUE
#' @param quietly Default: FALSE
#' @param verbose Default: getOption("verbose")
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
libraryMMV <- function(package, help = NULL, pos = 2, lib.loc = NULL,
                       logical.return = FALSE,
                       warn.conflicts = TRUE, quietly = FALSE,
                       verbose = getOption("verbose")){

  # If not loaded, load it:
  if (!(package %in% (.packages()))){

    # Load with required:
    Loaded <- require(package,
                      character.only = TRUE,
                      lib.loc        = lib.loc,
                      quietly        = quietly,
                      warn.conflicts = warn.conflicts)

    # If it did not load, need to install it:
    if(!Loaded){
      install.packages(package)
      library(package,
              help           = help,
              pos            = pos,
              character.only = TRUE,
              lib.loc        = lib.loc,
              logical.return = logical.return,
              quietly        = quietly,
              warn.conflicts = warn.conflicts,
              verbose        = verbose)
    }
  }
}
#' log-Linear Trapez Integration Method
#'
#' @description
#' @param x Numerical vector such as y=f(x)
#' @param y Numerical vector such as y=f(x)
#' @param logCrieria Criteria to use the log-lin rule (Default: 0.99)
#'
#' @return
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
logLinTrapzMMV <- function(
  x,
  y,
  logCrieria = 0.99){

  # Length of X & Y:
  nX <- length(x)
  nY <- length(y)

  # Check that X & Y have the same length:
  if (nX!=nY) {stop("X & Y have different length. Please Check.")}

  # Check that the length is at least 2:
  if (nX<2) {stop("X & Y should be of length 2 minimum.")}

  # Estimate Integrale:
  #   Test if the logLin Trapeze method should be used
  idx_Test <- ((y[2:nY]<logCrieria*y[1:(nY-1)]) & y[1:(nY-1)]>0 & y[2:nY]>0)
  #   Integral
  if (any(idx_Test)){
    Integral <- sum(ifelse(idx_Test,
                           (y[1:(nY-1)] - y[2:nY])/(log(y[1:(nY-1)])-log(y[2:nY]))*(x[2:nX]-x[1:(nX-1)]),
                           (y[1:(nY-1)] + y[2:nY])/2                              *(x[2:nX]-x[1:(nX-1)])))
  }else{
    Integral <- sum((y[1:(nY-1)] + y[2:nY])*(x[2:nX]-x[1:(nX-1)]))/2
  }

  # Output:
  return(Integral)
}

#' newtonRaphson.Function
#'
#' @description
#' @param funcToEval
#' @param Ysol Default: 0
#' @param x0 (Default: \code{NULL}).
#' @param range Default: c(-10, 10)
#' @param relTol Default: 1e-06
#' @param absTol Default: 1e-09
#' @param n Default: 1000
#' @param iterMax Default: 1000
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
newtonRaphson.Function <- function(
  funcToEval,        # Function to evaluate for the equation
  Ysol    = 0,         # The value such as "newtonRaphson.Function" returns the solution Xsol of Ysol=f(Xsol)
  x0      = NULL,      # Initial value for the algorithm. If not provided, it will analyze the derivative to find the optimum start point.
  range   = c(-10,10), # Numeric vector of two scalars that determine the range of the search
  relTol  = 1e-6,      # Relative tolerance of the solution
  absTol  = 1e-9,      # Absolute tolerance of the solution
  n       = 1000,      # number of interval to estimate the derivative
  iterMax = 1000       # Maximum number of iteration
) {

  # Check that 'range' is properly defined:
  if (!is.numeric(range) && length(unique(range))!=2){
    stop("'range' should be a numeric vector of length 2 in newton.Function.R")
  }

  # Order 'range':
  range <- sort(range)

  # Check if x0 is NULL:
  #   'range' will be used to estimate the optimum x0 to start with using the derivative
  if (is.null(x0)){
    # Calculate optimum x0:
    #   Estimate funcToEval at various points
    X <- seq(range[1], range[2], by=(range[2]-range[1])/n)
    Y <- numeric(length=n+1)
    for (k in 1:(n+1)){
      Y[k] <- funcToEval(X[k])
    }
    #   Estimate the derivative of the function funcToEval at X
    dY    <- numeric(length=n+1)
    dY[1] <- (Y[2]-Y[1])/(X[2]-X[1])
    dY[n] <- (Y[n]-Y[n-1])/(X[n]-X[n-1])
    for (i in 2:(n-1)){
      dY[i] <- (Y[i+1]-Y[i-1])/(X[i+1]-X[i-1])
    }

    # Estimate optimum x0::
    ix = which.max(abs(dY))  # Start where the derivative is maximum
    x0 = X[ix]
  }

  # Initialize Algorithm:
  X0 <- ifelse(x0==0,
               1e3*absTol,
               x0*(1+1e3*relTol)+sign(x0)*1e3*absTol)
  X1 <- x0
  k  <- 0
  dx <- (range[2]-range[1])/n

  # Iterate:
  while((abs(X1-X0)>max(absTol,relTol*max(abs(X0),abs(X1)))) & (k<iterMax)){
    # Update X0:
    X0 <- X1

    # Estimate f0:
    f0 <- funcToEval(X0)

    # Estimate f0':
    f0_minus <- funcToEval(X0-dx)
    f0_plus  <- funcToEval(X0+dx)
    df0      <- (f0_plus-f0_minus)/(2*dx)

    # Estimate X1:
    X1 = X0 - f0/df0

    # Update k:
    k  = k+1
  }

  # Output:
  return(X1)
}
#' newtonRaphson.Vector
#'
#' @description
#' @param X
#' @param Y
#' @param Ysol Default: 0
#' @param x0 Default: FALSE
#' @param relTol Default: 1e-06
#' @param absTol Default: 1e-09
#' @param eps Default: 1e-06
#' @param iterMax Default: 1000
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
newtonRaphson.Vector <- function(
  X,                # Vector X of where the value of the function f are known
  Y,                # Vector Y such as Y=f(X)
  Ysol    = 0,      # The value such as "newtonRaphson.Vector" returns the solution Xsol of Ysol=f(Xsol)
  x0      = FALSE,  # Initial value for the algorithm. If not provided, it will analyze the derivative to find the optimum start point.
  relTol  = 1e-6,   # Relative tolerance of the solution
  absTol  = 1e-9,   # Absolute tolerance of the solution
  eps     = 1e-6,   # Step to use to estimate derivative
  iterMax = 1000    # Maximum number of iteration
) {

  # Length of Y:
  n <- length(Y)

  # Check that X has the same length:
  if (length(X)!=n) {stop("X & Y have different length. Please Check.")}

  # Increasing or Decreasing:
  if (Y[1]<=Y[n]){Increase=TRUE} else{Increase=FALSE}

  # First Check if Ysol is in the range of Y: i.e. min(Y)<=Ysol<=max(Y)
  if (Ysol>=max(Y)){
    warning("Ysol is greater than the maximum value in Y. Returns Xsol such as Ymax=f(Xsol)")
    if (Increase) {
      ix <- which.max(rev(Y))
      ix <- length(Y)-ix+1
    } else {
      ix <- which.max(Y)
    }
    return(X[ix])
  }
  if (Ysol<=min(Y)){
    warning("Ysol is smaller than the minimum value in Y. Returns Xsol such as Ymin=f(Xsol)")
    if (Increase) {
      ix <- which.min(Y)
    } else {
      ix <- which.min(rev(Y))
      ix <- length(Y)-ix+1
    }
    ix <- which.min(Y)
    return(X[ix])
  }

  # Change Y such as the equation becomes g(Xsol)=0, ehere g(x)=f(x)-Ysol
  Ynorm <- Y - Ysol

  # Define X0 & X1
  if (is.numeric(x0)){
    X0 <- ifelse(x0==0,
                 1e3*absTol,
                 x0*(1+1e3*relTol)+sign(x0)*1e3*absTol)
    X1 <- x0
  } else{
    # Estimate the derivative of the function f at X to find the optimum starting point:
    dY    <- numeric(n)
    dY[1] <- (Y[2]-Y[1])/(X[2]-X[1])
    dY[n] <- (Y[n]-Y[n-1])/(X[n]-X[n-1])
    for (i in 2:(n-1)){
      dY[i] <- (Y[i+1]-Y[i-1])/(X[i+1]-X[i-1])
    }

    # Initialize Algorithm:
    ix <- which.max(abs(dY))  # Start where the derivative is maximum
    X0 <- ifelse(X[ix]==0,
                 1e3*absTol,
                 X[ix]*(1+1e3*relTol)+sign(X[ix])*1e3*absTol)
    X1 <- X[ix]
  }

  k <- 0
  while((abs(X1-X0)>max(absTol,relTol*max(abs(X0),abs(X1)))) & (k<iterMax)){
    # Update X0:
    X0 <- X1

    # Estimate f0:
    res <- approx(x=X, y=Ynorm, xout=X0, rule=2)
    f0  <- res$y

    # Estimate f0':
    if ((X0-eps)>=min(X)){
      if ((X0+eps)<=max(X)){
        res      <- approx(x=X, y=Ynorm, xout=X0-eps, rule=2)
        f0_minus <- res$y
        res      <- approx(x=X, y=Ynorm, xout=X0+eps, rule=2)
        f0_plus  <- res$y
        df0      <- (f0_plus-f0_minus)/(2*eps)
      }
      else{
        res      <- approx(x=X, y=Ynorm, xout=X0-eps, rule=2)
        f0_minus <- res$y
        df0      <- (f0-f0_minus)/eps
      }
    }else if ((X0+eps)<=max(X)){
      res      <- approx(x=X, y=Ynorm, xout=X0+eps, rule=2)
      f0_plus  <- res$y
      df0      <- (f0_plus-f0)/(eps)
    }else{
      res      <- approx(x=X, y=Ynorm, xout=min(X), rule=2)
      f0_minus <- res$y
      res      <- approx(x=X, y=Ynorm, xout=max(X), rule=2)
      f0_plus  <- res$y
      df0      <- (f0_plus-f0_minus)/(max(X)-min(X))
    }

    # Estimate X1:
    X1 <- X0 - f0/df0
    k  <- k+1
  }

  return(X1)
}

#' saveMMV
#' Save R objects
#'
#' [saveMMV] writes an external representation of R objects to the specified file.
#' It uses the function [base::save]. If 'file' is a path, it will create the directory
#' if not existent unlike [base::save].
#'
#' @param list A character vector containing the names of objects to be saved.
#' @param file The name of the file where the data will be saved.
#'
#' @md
#'
#' @export
#' @seealso [base::save]
#' @family General Functions
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
saveMMV <- function(list,
                    file){

  # Get directory name:
  dirName <- dirname(file)

  # Check if the directory name exists, otherwise create it:
  if (!dir.exists(dirName)){
    dir.create(dirName, recursive = TRUE)
  }

  # save list of object:
  save(list = list,
       file = file)
}
#' setMKLthreadsMMV
#' Set the Number of Threads for Math Kernel Library
#'
#' Set the maximum number of threads that can be started by the Math Kernel Library BLAS.
#' This should be less than or equal to the number of processing cores on your computer.
#' It is based on [setMKLthreads] from package [RevoUtilsMath], but it first checks if [RevoUtilsMath]
#' is loaded (As it is absent when 'R' is used instead of 'MRO').
#'
#' @param n An integer specifying the maximum number of threads.
#'
#' @md
#'
#' @seealso [RevoUtilsMath::setMKLthreads]
#' @family General Functions
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
setMKLthreadsMMV <- function(n) {

  # Information about the current session:
  sessInfo <- sessionInfo()

  # Set MKL thread:
  if ("RevoUtilsMath" %in% c(sessInfo$basePkgs, names(sessInfo$otherPkgs))){
    setMKLthreads(n)
  }else{
    warning("'RevoUtilsMath' is not installed, therefore, 'setMKLthreads' is not available and the number of used thread is 1 as in old good R.")
  }
}
#' simpsonMMV
#'
#' @description
#' @param x
#' @param y
#' @param tolInt Default: 1e-09
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
simpsonMMV <- function(
  x,               # Vector x of where the value of the function f are known
  y,               # Vector y such as y=f(x)
  tolInt = 1e-9    # Tolerance to use to check if the step size is similar among the interval
) {

  # Length of X & Y:
  nX <- length(x)
  nY <- length(y)

  # Check that X & Y have the same length:
  if (nX!=nY) {stop("X & Y have different length. Please Check.")}

  if (nX<2){
    warning("Length of X & Y less than two: Integral set to 0.")
    Integral <- 0

  }else if(nX==2){
    # Trapez Method:
    Integral <- (x[2] - x[1])*(y[1] + y[2])/2

  }else if(nX==3){
    # Check if the integration step dx have same length:
    if ((abs(x[1] + x[3] -2*x[2])<(x[2]-x[1])*tolInt)){
      # Simpson Method:
      Integral <- (x[3] - x[1])*(y[1] + 4*y[2] + y[3])/6

    }else{
      # Trapez Method:
      Integral <- sum((y[1:(nY-1)] + y[2:nY])*(x[2:nX]-x[1:(nX-1)]))/2
    }

  }else{
    # First Check if all sub-interval of x have similar length:
    dx   <- x[2:nX]-x[1:(nX-1)]
    TEST <- (abs(dx-dx[1])<dx[1]*tolInt)

    # Initialize with first interval using the trapez method:
    Integral <- (x[2]-x[1])*(y[1]+y[2])/2

    if (all(TEST)){
      # Simpson Method simplified for constant dx step:
      Integral <- Integral + sum((x[3:nX] - x[1:(nX-2)])*(y[1:(nX-2)] + 4*y[2:(nX-1)] + y[3:nX]))/6

    }else{
      # Simpson Method by starting with the first point:
      #   Here, as the points are not equally distant, the formula
      #   had to be adapted to account for it.

      # Some constants:
      Cak <- y[1:(nX-2)]/(x[1:(nX-2)] - x[2:(nX-1)])/(x[1:(nX-2)] - x[3:nX]    )
      Cmk <- y[2:(nX-1)]/(x[2:(nX-1)] - x[1:(nX-2)])/(x[2:(nX-1)] - x[3:nX]    )
      Cbk <- y[3:nX]    /(x[3:nX]     - x[1:(nX-2)])/(x[3:nX]     - x[2:(nX-1)])

      # More constants:
      Ak <-   Cak                         + Cmk                         + Cbk
      Bk <- -(Cak*(x[2:(nX-1)] + x[3:nX]) + Cmk*(x[1:(nX-2)] + x[3:nX]) + Cbk*(x[1:(nX-2)] + x[2:(nX-1)]))
      Ck <-   Cak*(x[2:(nX-1)] * x[3:nX]) + Cmk*(x[1:(nX-2)] * x[3:nX]) + Cbk*(x[1:(nX-2)] * x[2:(nX-1)])

      # Integral:
      Integral <- Integral + sum(Ak/3*(x[3:nX]^3 - x[1:(nX-2)]^3) + Bk/2*(x[3:nX]^2 - x[1:(nX-2)]^2) + Ck*(x[3:nX] - x[1:(nX-2)]))
    }

    # Add last interval by the trapez method:
    Integral <- Integral + (x[nX]-x[nX-1])*(y[nX-1]+y[nX])/2

    # Each Interval was counted twice:
    Integral <- Integral/2
  }

  # Output:
  return(Integral)
}
#' thisFile
#'
#' @description
#' @param fullName Default: FALSE
#' @param pathNormalized Default: FALSE
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
thisFile <- function(fullName = FALSE, pathNormalized = FALSE) {

  # use commandArgs:
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle  <- "--file="
  match   <- grep(needle, cmdArgs)

  # Test with cmdArgs
  if (length(match) > 0) {
    # Rscript
    thisFileName <- sub(needle, "", cmdArgs[match])

    # Otherwise use sys.frames
  } else {
    # 'source'd via R console
    #return(normalizePath(sys.frames()[[1]]$ofile))
    thisFileName <- sys.frames()[[1]]$ofile
  }

  # if thisFIleName is NULL it is called from the console:
  if (is.null(thisFileName)){
    thisFileName <- "Console"
  }

  # Keep only fiel name if needed:
  #   In that case the path can not be normilized
  if(!fullName){
    split_fileName <- strsplit(thisFileName, split="/")
    thisFileName   <- split_fileName[[1]][length(split_fileName[[1]])]
    pathNormalized <- FALSE
  }

  # Normalized path if needed
  if (pathNormalized){
    thisFileName <- normalizePath(thisFileName)
  }

  return(thisFileName)
}
#' rectintMMV
#'
#' @description
#' @param x
#' @param y
#' @param FLAGright
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
rectintMMV <- function(x,               # Vector x of where the value of the function f are known
                       y,                # Vector y such as y=f(x)
                       FLAGright = FALSE
) {

  # Length of X & Y:
  nX <- length(x)
  nY <- length(y)

  # Check that X & Y have the same length:
  if (nX!=nY) {stop("X & Y have different length. Please Check.")}

  # Check that the length is at least 2:
  if (nX<2) {stop("X & Y should be of length 2 minimum.")}

  # Estimate Integrale
  #   Right Riemann sum
  if(FLAGright){
    Integral <- sum(y[2:nY]*(x[2:nX]-x[1:(nX-1)]))

  #   Left Riemann sum
  }else{
    Integral <- sum(y[1:(nY-1)]*(x[2:nX]-x[1:(nX-1)]))
  }


  # Output:
  return(Integral)
}

#' trapzMMV
#'
#' @description
#' @param x
#' @param y
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
trapzMMV <- function(x,               # Vector x of where the value of the function f are known
                     y                # Vector y such as y=f(x)
) {

  # Length of X & Y:
  nX <- length(x)
  nY <- length(y)

  # Check that X & Y have the same length:
  if (nX!=nY) {stop("X & Y have different length. Please Check.")}

  # Check that the length is at least 2:
  if (nX<2) {stop("X & Y should be of length 2 minimum.")}

  # Estimate Integrale
  Integral <- sum((y[1:(nY-1)] + y[2:nY])*(x[2:nX]-x[1:(nX-1)]))/2

  # Output:
  return(Integral)
}

#' Combine R Objects by Rows
#'
#' Take a sequence of vector, matrix or data-frame arguments and combine by rows,
#' respectively. The difference with base::rbind is that it will still work
#' if one of the data frame is NULL or empty (i.e. no row).
#'
#' @description
#' @param x A "data.frame"
#' @param y A "data.frame"
#' @param stringsAsFactors logical, passed to as.data.frame
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
rbindMMV <- function(x,
                     y,
                     stringsAsFactors = default.stringsAsFactors()){

  # Bind if both are defined, or just return the non-emtpy one:
  if ((is.null(x) || nrow(x)==0) & !(is.null(y) || nrow(y)==0)){
    out <- y
  }else if (!(is.null(x) || nrow(x)==0) & (is.null(y) || nrow(y)==0)){
    out <- x
  }else if (!(is.null(x) || nrow(x)==0) & !(is.null(y) || nrow(y)==0)){
    out <- rbind(x,
                 y,
                 stringsAsFactors = stringsAsFactors)
  }else{
    out <- NULL
    # stop("Both data frame 'x' and 'y' are empty: At least one should be defined")
  }

  # Output:
  return(out)
}

#' Combine R Objects by Columns
#'
#' Take a sequence of vector, matrix or data-frame arguments and combine by columns,
#' respectively. The difference with base::cbind is that it will still work
#' if one of the data frame is NULL or empty (i.e. no row).
#'
#' @description
#' @param x A "data.frame"
#' @param y A "data.frame"
#' @param stringsAsFactors logical, passed to as.data.frame
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
cbindMMV <- function(x,
                     y,
                     stringsAsFactors = default.stringsAsFactors()){

  # Bind if both are defined, or just return the non-emtpy one:
  if ((is.null(x) || nrow(x)==0) & !(is.null(y) || nrow(y)==0)){
    out <- y
  }else if (!(is.null(x) || nrow(x)==0) & (is.null(y) || nrow(y)==0)){
    out <- x
  }else if (!(is.null(x) || nrow(x)==0) & !(is.null(y) || nrow(y)==0)){
    out <- cbind(x,
                 y,
                 stringsAsFactors = stringsAsFactors)
  }else{
    out <- NULL
    # stop("Both data frame 'x' and 'y' are empty: At least one should be defined")
  }

  # Output:
  return(out)
}

#' combMMV
#'
#' @description
#' @param x
#' @param ...
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
combMMV <- function(x, ...) {
  out <- lapply(seq_along(x),
                function(i){
                  plyr::ldply(list(...), .fun = function(y){
                    plyr::rbind.fill(x[[i]],y[[i]], stringsAsFactors = FALSE)
                  })})

  return(out)
}


#' Retrieve Package Path
#'
#' @description
#' @param subdir
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
get_MMVmalariaPath <- function(subdir = NULL){

  # Depending on how MMVmalaria was installed:
  #   Sources for development:
  if(exists("MMVmalariaPathDev")){
    if(is.null(subdir)){
      out <- file.path(MMVmalariaPathDev)
    }else{
      out <- file.path(MMVmalariaPathDev, subdir)
    }

    #   Installed as a package
  }else if("MMVmalaria" %in% .packages()){
    out <- system.file("",package = "MMVmalaria")

  }else {
    stop("Neither 'MMVmalaria' is installed as a package or 'MMVmalariaPathDev' exist")
  }

  # Output:
  return(out)
}

#' Reverse a List
#'
#' Reverse a list, where all elements of the list have the same length and
#' similar object:
#' If the input is
#'      \code{a = list(list(Obj1 = data.frame(X = c(1,2), Y = c(1,4))),
#'                     list(Obj1 = data.frame(X = c(3,4), Y = c(9,16))))}
#' the output will be
#'      \code{list(Obj1 = list(data.frame(X = c(1,2,3,4), Y = c(1,4,9,16)),
#'                             data.frame(X = c(3,4), Y = c(1,4))))}
#'
#' @param list List to reverse
#' @param simplify Logical or character string; should the result be simplified
#' to a vector, matrix or higher dimensional array if possible? For \code{sapply}
#' it must be named and not abbreviated. The default value, \code{TRUE}, returns
#' a vector or matrix if appropriate, whereas if \code{simplify = "array"} the result
#' may be an array of "rank" (\code{=length(dim(.))}) one higher than the result of
#' \code{FUN(X[[i]])} (Default: \code{FALSE})
#'
#' @return Reversed list
#'
#' @examples
#' a = list(list(Obj1 = data.frame(X = c(1,2), Y = c(1,4))),
#'          list(Obj1 = data.frame(X = c(3,4), Y = c(9,16))))
#' reverse_List(a)
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
reverse_List <- function (list, simplify = FALSE)
{
  if (length(unique(sapply(list, length))) != 1) {
    stop("Not all lists equally long")
  }
  list1 <- list()
  for (i in 1:length(list[[1]])) {
    list1[[i]] <- sapply(list, function(x) x[[i]], simplify = simplify)
  }
  names(list1) <- names(list[[1]])
  return(list1)
}


#' Load RData into a list
#'
#' @param fileName Path to the RData file name.
#'
#' @return List with all object saved in `fileName`
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
load_RData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  variableNames <- ls()[!(ls() %in% c("fileName"))]
  out <- lapply(variableNames,function(x){
    get(x)
  })
  names(out) <- variableNames
  out
}

#' Reduce character to the non-repeating entries
#'
#' Entries which are the same as the preceding one are are replaced by empty character.
#' The inverse of insert_duplicates().
#'
#' @param x character vector.
#' @return character of the same length as x.
#'
#' @export
remove_duplicates <- function(x) {

  x <- as.character(x)

  N <- length(x)
  if (N == 1) return(x)

  is_same_as_preceding <- x[2:N] == x[1:(N-1)]
  x[2:N][is_same_as_preceding] <- ""

  x

}

#' Remove parameters with NA in VALUE column from the GPF object
#' @param gpf a GPF object or filepath to an gpf file.
#' @return a GPF object
#'
#' @export
removeParametersWithNAValue_GPF <- function(gpf) {
  if(is.character(gpf)) {
    gpf <- suppressWarnings(GPF(gpf))
  } else if(!is_GPF(gpf)) {
    stop("removeParametersWithNAValue_GPF: argument gpf should be a gpf object or a file path to a GPF excel file.")
  }
  paramsWithNAValue <- gpf$estimates$PARAMETER[is.na(gpf$estimates$VALUE)]
  paramsToKeep <- setdiff(gpf$estimates$PARAMETER, paramsWithNAValue)
  if(length(paramsToKeep) == 0) {
    stop("removeParametersWithNAValue_GPF: the provided GPF object does not have non-NA values.")
  }
  gpf$estimates <- gpf$estimates[gpf$estimates$PARAMETER %in% paramsToKeep,]
  gpf$uncertainty_correlation <- gpf$uncertainty_correlation[gpf$uncertainty_correlation$PARAMETER %in% paramsToKeep,
                                                             colnames(gpf$uncertainty_correlation) %in% c("PARAMETER", paramsToKeep),
                                                             with = FALSE]
  gpf
}
