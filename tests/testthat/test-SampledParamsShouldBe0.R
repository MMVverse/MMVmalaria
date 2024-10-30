context("test-SampledParamsShouldBe0")

test_that("Individual parameters with zero IIV should be constant", {

  xls  <- system.file("extdata", "test-SampledParamsShouldBe0", "XLS.xls", package = "MMVmalaria")
  spec <- IQRtools:::specifyParamSampling(xls)
  IQRtools:::getParamEstimates(spec)
  pars <- MMVmalaria::sample_MMVmalariaProject(xls, Nsamples = 10)

  # Expect all to be equal to zero
  expect_true(all(pars$parameterValuesIndividual$dadp == 0))
})
