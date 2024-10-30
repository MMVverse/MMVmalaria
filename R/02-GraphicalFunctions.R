#' Adjust ggplot Theme
#'
#' Adjust the theme of the ggplot object `gr` to various style (e.g. to MMV style).
#'
#' @param gr ggplot object to adjust theme.
#' @param style Style to choose between \code{"MMV"} or \code{"MMVvpc"} (Default: \code{"MMV"}).
#'
#' @return `gr` with adjusted theme
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family Graphical Functions
adjust_ggplotTheme <- function(gr,
                               style = c("MMV","MMVvpc")[1]) {

  # Quick Check:
  if(is.null(style)){
    style <- "original"
  }

  # MMV Style
  if (grepl("MMV", toupper(style))){
    # Adjust theme for Sub-style:
    #   VPC
    if (grepl("VPC", toupper(style))){
      gr <- gr + theme(axis.title    = element_text(face="bold"),
                       plot.title    = element_text(size=16),
                       plot.subtitle = element_text(size=14),
                       legend.text   = element_text(size=10),
                       legend.title  = element_text(size=12),
                       axis.text.x   = element_text(size=8 ),
                       axis.text.y   = element_text(size=8 ),
                       axis.title.x  = element_text(size=10),
                       axis.title.y  = element_text(size=10),
                       strip.text    = element_text(size=10),
                       plot.caption  = element_text(size=8, hjust=0))

      #   Any Other:
    }else{
      gr <- gr + theme(axis.title    = element_text(face="bold"),
                       plot.title    = element_text(size=18),
                       plot.subtitle = element_text(size=16),
                       legend.text   = element_text(size=12),
                       legend.title  = element_text(size=14),
                       axis.text.x   = element_text(size=14),
                       axis.text.y   = element_text(size=14),
                       axis.title.x  = element_text(size=14),
                       axis.title.y  = element_text(size=14),
                       strip.text    = element_text(size=12),
                       plot.caption  = element_text(size=10, hjust=0))
    }

    # Any other style:
  } else{
    gr <- gr
  }

  # Output:
  gr
}


#' ggplot functionality implementing MMV style
#'
#' ggplot function similar to `ggplot2::ggplot` implementing MMV style.
#'
#' @param ...          Typical ggplot input arguments
#' @param style        Style to choose between \code{"MMV"}, \code{"MMVvpc"} or \code{"IQR"} (Default: \code{"MMV"}).
#' @param ActivityPath Path of the current activity (Default: \code{NULL})
#' @param Caption      Caption to add to the plot (Default: \code{NULL})
#'
#' @return A ggplot object with MMV standars.
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family GraphicalFunctions
MMVggplot <- function(...,
                      style = c("MMV","MMVvpc","IQR")[1],
                      ActivityPath = NULL,
                      Caption      = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Use Style of Interest ----
  #-----------------------------------------------------------------------------#

  # MMV Style
  if (grepl("MMV", toupper(style))){
    gr <- IQRtools::IQRggplot(...) +
      scale_color_manual(values = IQRtoolsColors[2:length(IQRtoolsColors)]) +
      scale_fill_manual(values = IQRtoolsColors[2:length(IQRtoolsColors)])

  # IQR Style:
  } else if (toupper(style)=="IQR"){
    gr <- IQRtools::IQRggplot(...)

  # Normal ggplot:
  } else{
    gr <- ggplot(...)
  }


  #-----------------------------------------------------------------------------#
  # STEP 2: Adjust Activity Path ----
  #-----------------------------------------------------------------------------#
  if (is.null(ActivityPath) || ActivityPath!=""){
    ActivityPath <- paste0("Activity: ",get_ActivityPath(ActivityPath))
  }


  #-----------------------------------------------------------------------------#
  # STEP 3: Caption ----
  #-----------------------------------------------------------------------------#

  # Adjust Caption:
  if(is.character(Caption)){
    Caption <- paste0(ActivityPath, "\n", Caption)
  }else{
    Caption <- ActivityPath
  }

  # Add Caption
  gr <- gr + labs(caption = Caption)


  #-----------------------------------------------------------------------------#
  # STEP 4: Adjust theme ----
  #-----------------------------------------------------------------------------#
  gr <- adjust_ggplotTheme(gr    = gr,
                           style = style)


  #-----------------------------------------------------------------------------#
  # STEP 5: Output ----
  #-----------------------------------------------------------------------------#
  return(gr)
}


#' Transform IQR ggplot object to MMV ggplot object
#'
#' When a ggplot object is generated using `IQRtools` functions, this functions
#' allows to adjust the ggplot object to be compliante with MMV standars;
#' It will add the path of the activity and adjust the theme if specified.
#'
#' @param IQRggplot ggplot object generated with `IQRtools` functions.
#' @param style Style to choose between \code{"MMV"} or \code{"MMVvpc"} (Default: \code{"MMV"}).
#' @param ActivityPath Path of the current activity (Default: \code{NULL})
#' @param Caption Caption to add to the plot (Default: \code{NULL})
#'
#' @return A ggplot object compliant with MMV standars.
#'
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family GraphicalFunctions
transform_IQRggplotToMMVggplot <- function(IQRggplot,
                                           style        = c("MMV","MMVvpc")[1],
                                           ActivityPath = NULL,
                                           Caption      = NULL) {

  #-----------------------------------------------------------------------------#
  # STEP 1: Use Style of Interest ----
  #-----------------------------------------------------------------------------#

  # Quick Check:
  if(is.null(style)){
    style <- "original"
  }

  # MMV Style
  if (grepl("MMV", toupper(style))){
    base::suppressMessages(MMVggplot <- IQRggplot +
                           scale_color_manual(values = IQRtoolsColors[2:length(IQRtoolsColors)]))

  # Any other style:
  } else{
    MMVggplot <- IQRggplot
  }


  #-----------------------------------------------------------------------------#
  # STEP 2: Adjust Activity Path ----
  #-----------------------------------------------------------------------------#
  if (is.null(ActivityPath) || ActivityPath!=""){
    ActivityPath <- paste0("Activity: ",get_ActivityPath(ActivityPath))
  }


  #-----------------------------------------------------------------------------#
  # STEP 3: Caption ----
  #-----------------------------------------------------------------------------#

  # Adjust Caption:
  if(is.character(Caption)){
    Caption <- paste0(ActivityPath, "\n", Caption)
  }else{
    Caption <- ActivityPath
  }

  # Add Caption
  MMVggplot <- MMVggplot + labs(caption = Caption)


  #-----------------------------------------------------------------------------#
  # STEP 4: Adjust theme ----
  #-----------------------------------------------------------------------------#
  MMVggplot <- adjust_ggplotTheme(gr    = MMVggplot,
                                  style = style)


  #-----------------------------------------------------------------------------#
  # STEP 5: Output ----
  #-----------------------------------------------------------------------------#
  MMVggplot
}


