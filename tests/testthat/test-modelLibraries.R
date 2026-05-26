test_that("all model library files simulate successfully", {
  
  skip_if_not_installed("IQRtools")
  library(IQRtools)
  
  # Define model libraries. Will need to be adjusted
  # if libraries are added or removed 
  model_libraries <- c(
    "modelLibrary",
    "modelLibrary_tobs0"
  )
  
  model_dirs <- vapply(
    model_libraries,
    function(x) {
      system.file(x, package = "MMVmalaria")
    },
    character(1)
  )
  
  # get data frame of models including which library 
  # ! NOTE - this will pick up any .txt file... 
  # so will get an error (not IQRmodel file... ) if any non-iqrmodel.txt 
  # file is in the model library (i.e., a README?) - but currently the model
  # libraries have exclusively model files in them and it should be good to 
  # keep them this way 
  
  model_files_df <- do.call(
    rbind,
    lapply(names(model_dirs), function(library_name) {
      files <- list.files(
        path = model_dirs[[library_name]],
        pattern = "\\.txt$",
        recursive = TRUE,
        full.names = TRUE
      )
      
      data.frame(
        library = library_name,
        file = files,
        stringsAsFactors = FALSE
      )
    })
  )
  
  # Lapply over all the models.
  # withCallingHandlers will catch + store warnings (and muffle them) then
  # allow the code to continue. 
  # It is wrapped in tryCatch which ensures that the lapply will continue 
  # even if there is an error in one of the models. Having withCallingHandlers 
  # within the tryCatch means that we can catch + store all warnings whereas 
  # tryCatch alone will exit on the first error but does not catch warnings 
  # (and the warnings are useful for use as part of a test suite) 
  #  https://adv-r.hadley.nz/conditions.html#:~:text=The%20handlers%20set%20up%20by,normally%20once%20the%20handler%20returns.

  
  results <- lapply(seq_len(nrow(model_files_df)), function(i) {
    
    library_name <- model_files_df$library[i]
    file <- model_files_df$file[i]
    
    warnings <- new.env(parent = emptyenv())
    warnings$messages <- character(0)
    
    error <- tryCatch(
      {
        withCallingHandlers(
          {
            sim_IQRmodel(IQRmodel(file))
            NA_character_
          },
          warning = function(w) {
            warnings$messages <- c(
              warnings$messages,
              conditionMessage(w)
            )
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) {
        conditionMessage(e)
      }
    )
    
    data.frame(
      library = library_name,
      file = file,
      passed = is.na(error),
      error = error,
      warnings = if (length(warnings$messages) > 0) {
        paste(unique(warnings$messages), collapse = " | ")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  
  # drop = FALSE to make sure failed is a df, even if only one failed 
  # model. 
  
  failed <- results_df[!results_df$passed, , drop = FALSE]
  
  failure_message <- if (nrow(failed) > 0) {
    paste(
      "The following IQR models failed sim_IQRmodel():",
      paste(
        sprintf(
          "\n- Library: %s\n  File: %s\n  Error: %s",
          failed$library,
          failed$file,
          failed$error
        ),
        collapse = "\n"
      ),
      sep = "\n"
    )
  } else {
    sprintf(
      "All %s IQR model(s) across %s model librar(ies) simulated successfully.",
      nrow(results_df),
      length(model_libraries)
    )
  }
  
  # test that handler 
  expect_true(
    all(results_df$passed),
    label = failure_message
  )
})