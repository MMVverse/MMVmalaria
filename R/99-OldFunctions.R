#' catCovToTxt
#'
#' @description
#' @param data
#' @param cat0ColNames Default: NULL
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV)
#' @family DataPreparation
catCovToTxt <- function(data,
                        cat0ColNames = NULL  # Vector of categorical COV column names for which TXT value is expected
){
  a <- attributes(dataGeneral)
  covCAT <-FALSE

  for(k in seq_along(cat0ColNames)){

    # Select covariate
    cov_k  <- cat0ColNames[k]

    # Check if covariate exist and / or is categorical
    if(cov_k %in% a$catInfo$COLNAME){
      covCAT <- TRUE
    }else if(cov_k %in% a$covInfo$COLNAME){
      stop("'", cov_k, "' is not a continuous covariate, no need to define a text value.")
    }else{
      stop("'", cov_k, "' is not a covariate defined in the dataset.")
    }

    # Identify text matching covariate number
    if (covCAT){
      idx_COL <- which(a$catInfo$COLNAME==cov_k)

      dataD <- data.frame(VALUE    = as.numeric(unlist(strsplit(a$catInfo$VALUES[idx_COL], split = ",", fixed = TRUE))),
                          VALUETXT = strsplit(a$catInfo$VALUETXT[idx_COL], split = ",", fixed = TRUE)[[1]],
                          stringsAsFactors = FALSE)

    }
    # Remove last character N from cov_k to create text column name
    dataGeneral[,gsub('.{1}$',"",cov_k)]   <- dataD$VALUETXT[dataGeneral[,cov_k]]
  }
  # Output
  return(dataGeneral)
}


#' convertGametocytes
#'
#' @description
#' @param data
#' @param gamNameOld Default: 'Parasitemia female gametocytes'
#' @param gamNameNew Default: 'Parasitemia female gametocytes'
#' @param oldStdCurve
#' @param convertu2m
#' @param elutionFactor Default: 1
#' @param femaleGamFactor Default: 590
#' @param oldNewStdFactor Default: 62
#' @return
#' @export
#' @author Aline Fuchs (MMV)
#' @family Data Preparation
convertGametocytes <- function(
  data,                                               # Data frame or IQRdataGENERAL object with gametocyte data
  gamNameOld      = "Parasitemia female gametocytes",
  gamNameNew      = "Parasitemia female gametocytes",
  oldStdCurve,                                        # TRUE or FALSE indicating whether old standard curve was used
  convertu2m,                                         # TRUE or FALSE indicating whether original counts given per uL
  elutionFactor   = 1,                                # Elution factor of RNA sample (need to check which lab)
  femaleGamFactor = 590,                              # Conversion factor from count/mL to p/mL (p being female gametocytes)
  oldNewStdFactor = 62                                # Conversion factor from old to new standard curve
) {
  stop("'convertGametocytes' is deprecated and substituted by 'convert_Gametocytes'. Please switch to the new function.")
}

