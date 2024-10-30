context("test-generate_GPFFromIQRnlmeProject")

test_that("Generate an .xlsx file from a MONOLIX project directory", {
  # Load model example:
  nlmeProjectPath <- system.file("extdata", "IQRtoolsExampleProject", package = "MMVmalaria")
  fileXLS <- system.file("extdata", "IQRtoolsExampleProject", "PKModelCovariance.xlsx", package = "MMVmalaria")

  # Generation of a GPF object.
  xls <- generate_GPFFromIQRnlmeProject(projectPath = nlmeProjectPath,
                                        filename = fileXLS)

  # Test if it worked
  expect_true(file.exists(fileXLS))
})
