test_that("all model library files simulate successfully", {
  
  skip_if_not_installed("IQRtools")
  library(IQRtools)
  
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
  
  missing_dirs <- names(model_dirs)[model_dirs == "" | !dir.exists(model_dirs)]
  
  expect_true(
    length(missing_dirs) == 0,
    label = paste(
      "Missing model directories:",
      paste(missing_dirs, collapse = ", ")
    )
  )
  
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
  
  expect_true(
    nrow(model_files_df) > 0,
    label = paste(
      "No .txt model files found in:",
      paste(model_dirs, collapse = ", ")
    )
  )
  
  results <- lapply(seq_len(nrow(model_files_df)), function(i) {
    library_name <- model_files_df$library[i]
    file <- model_files_df$file[i]
    
    warning_messages <- character()
    
    error_message <- tryCatch(
      {
        withCallingHandlers(
          {
            sim_IQRmodel(IQRmodel(file))
            NULL
          },
          warning = function(w) {
            warning_messages <<- c(warning_messages, conditionMessage(w))
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
      passed = is.null(error_message),
      error = if (is.null(error_message)) NA_character_ else error_message,
      warnings = if (length(warning_messages) > 0) {
        paste(unique(warning_messages), collapse = " | ")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  
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
  
  expect_true(
    all(results_df$passed),
    label = failure_message
  )
})