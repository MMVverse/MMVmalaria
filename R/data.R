#-----------------------------------------------------------------------------#
# Malaria Dataset ----
#-----------------------------------------------------------------------------#

#' An example data-set of the format used within the MMVmalaria package. 
#' The 10 example individuals in this data-set have randomly sampled values 
#' and are not real malaria patients
#'
#' A dataset containing baseline parasitemia, age, body weight, body height,
#' sex, race, country and region of uncomplicated malaria patients.
#'
#' @format A data frame with 10 rows and 13 variables:
#' \describe{
#'   \item{ID}{Row ID}
#'   \item{USUBJID}{Unique Subject ID}
#'   \item{STUDY}{Study ID}
#'   \item{PBASE_perml}{Baseline Parasitemia, counts//mL}
#'   \item{PLBASE_log10_perml}{log10-transformed baseline parasitemia, log10(counts//mL)}
#'   \item{AGE_years}{Age, years}
#'   \item{WEIGHT_kg}{Body Weight, kg}
#'   \item{HEIGHT_cm}{Body Height, cm}
#'   \item{SEX}{Sex, Male or Female}
#'   \item{RACE}{Race}
#'   \item{COUNTRY}{Country}
#'   \item{REGION}{Region}
#'   \item{SUBJECTTYPE}{Asymptomatic Patient, Symptomatic Patient}
#'   \\
#' }
"MalariaPopulation"


