context("test-keyParameters")

library(MMVmalaria)

# Helper: build a synthetic concentration profile with two bumps above MIC separated by a
# dip below MIC (simulating a trough between two doses), plus a single-plateau (no gap)
# variant for comparison.
make_two_bump_profile <- function(MIC = 1) {
  t <- seq(0, 120, by = 0.05)
  # Two log-normal-shaped bumps peaking well above MIC, separated by a trough that dips
  # below MIC around t = 45.
  bump <- function(t, peakTime, height, width) height * exp(-((t - peakTime)^2) / (2 * width^2))
  Cc <- bump(t, 20, 4, 8) + bump(t, 70, 4, 8) + 0.05
  data.frame(TIME = t, Cc = Cc)
}

make_continuous_profile <- function(MIC = 1) {
  t <- seq(0, 120, by = 0.05)
  bump <- function(t, peakTime, height, width) height * exp(-((t - peakTime)^2) / (2 * width^2))
  # Same two bumps, but close enough together (and tall enough) that the trough between
  # them never drops below MIC - a single continuous excursion.
  Cc <- bump(t, 30, 4, 20) + bump(t, 60, 4, 20) + 0.05
  data.frame(TIME = t, Cc = Cc)
}

test_that("getExcursionsAboveMIC detects a real gap between two bumps", {
  MIC <- 1
  dfGap <- make_two_bump_profile(MIC)

  excursions <- getExcursionsAboveMIC(dfGap, MIC = MIC)
  expect_equal(nrow(excursions), 2)
  expect_true(all(excursions$tEnd > excursions$tOnset))
  # First excursion (around the t=20 bump) ends well before the second (around t=70) begins.
  expect_true(excursions$tEnd[1] < excursions$tOnset[2])
  # Roughly centered where expected (loose tolerance - crossing times depend on bump shape).
  expect_true(excursions$tOnset[1] > 5 && excursions$tOnset[1] < 20)
  expect_true(excursions$tOnset[2] > 55 && excursions$tOnset[2] < 70)
})

test_that("getExcursionsAboveMIC's summed tMIC matches getTimeAboveMIC exactly", {
  MIC <- 1
  for (dfSim in list(make_two_bump_profile(MIC), make_continuous_profile(MIC))) {
    excursions <- getExcursionsAboveMIC(dfSim, MIC = MIC)
    total      <- getTimeAboveMIC(dfSim, MIC = MIC)[["tMIC"]]
    expect_equal(sum(excursions$tMIC), total, tolerance = 1e-6)
  }
})

test_that("getExcursionsAboveMIC returns zero rows when MIC is never exceeded", {
  dfSim <- data.frame(TIME = seq(0, 50, by = 1), Cc = rep(0.01, 51))
  excursions <- getExcursionsAboveMIC(dfSim, MIC = 1)
  expect_equal(nrow(excursions), 0)
  expect_equal(getTimeAboveMIC(dfSim, MIC = 1)[["tMIC"]], 0)
})

test_that("getExcursionsAboveMIC returns a single excursion when always above MIC", {
  dfSim <- data.frame(TIME = seq(0, 50, by = 1), Cc = rep(10, 51))
  excursions <- getExcursionsAboveMIC(dfSim, MIC = 1)
  expect_equal(nrow(excursions), 1)
  expect_equal(excursions$tOnset[1], 0)
  expect_equal(excursions$tEnd[1], 50)
})

# evaluateDoseCriterion_ContinuousTimeAboveMIC ----

# Minimal EMAX parameters giving MIC = 1: with hill = 1, formula_EMAXmodelParsToMIC() computes
# MIC = GR/(EMAX-GR) * EC50, so EC50 = (EMAX-GR)/GR gives MIC = 1 exactly (GR = 0.05, EMAX = 0.5
# -> EC50 = 9).
emax_params_for_MIC1 <- c(GR = 0.05, EMAX = 0.5, EC50 = 9, hill = 1)

test_that("evaluateDoseCriterion_ContinuousTimeAboveMIC only counts the first excursion when there is a gap", {
  dfGap <- make_two_bump_profile(MIC = 1)
  dfGap$PL <- 0  # unused by these criteria but present on real sim_results

  continuousCrit <- evaluateDoseCriterion_ContinuousTimeAboveMIC(dfGap, emax_params_for_MIC1)
  summedCrit      <- evaluateDoseCriterion_TimeAboveMIC(dfGap, emax_params_for_MIC1)

  excursions <- getExcursionsAboveMIC(dfGap, MIC = 1)
  expect_equal(continuousCrit, excursions$tMIC[1], tolerance = 1e-6)
  # The gap means the first excursion alone is meaningfully shorter than the summed total.
  expect_true(continuousCrit < summedCrit - 1)
})

test_that("evaluateDoseCriterion_ContinuousTimeAboveMIC equals the summed criterion when there is no gap", {
  dfContinuous <- make_continuous_profile(MIC = 1)
  dfContinuous$PL <- 0

  continuousCrit <- evaluateDoseCriterion_ContinuousTimeAboveMIC(dfContinuous, emax_params_for_MIC1)
  summedCrit      <- evaluateDoseCriterion_TimeAboveMIC(dfContinuous, emax_params_for_MIC1)

  expect_equal(continuousCrit, summedCrit, tolerance = 1e-6)
})

test_that("evaluateDoseCriterion_ContinuousTimeAboveMIC returns 0 when MIC is never exceeded", {
  dfSim <- data.frame(TIME = seq(0, 50, by = 1), Cc = rep(0.01, 51), PL = 0)
  expect_equal(evaluateDoseCriterion_ContinuousTimeAboveMIC(dfSim, emax_params_for_MIC1), 0)
})
