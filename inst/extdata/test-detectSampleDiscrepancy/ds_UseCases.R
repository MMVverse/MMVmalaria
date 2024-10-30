library(IQRtools)
 
# define original data
obj <- "PKModelCovariance.xlsx"

# case studies to test
objAlt <-
  c(
    None = "PKModelCovariance.xlsx",
    CL_01 = "PKModelCovariance_CL_01.xlsx",
    CL_09 = "PKModelCovariance_CL_09.xlsx",
    CLQ1_IIV_corr_01 = "PKModelCovariance_CLQ1_corr_01.xlsx",
    CLQ1_IIV_corr_09 = "PKModelCovariance_CLQ1_corr_09.xlsx",
    CLWT0_beta_01 = "PKModelCovariance_CLWT0_beta_01.xlsx",
    CLWT0_beta_09 = "PKModelCovariance_CLWT0_beta_09.xlsx",
    CL_RSE_2x = "PKModelCovariance_CL_rse_2x.xlsx",
    CL_RSE_1.1x = "PKModelCovariance_CL_rse_1.1x.xlsx",
    VcCL_UncCorr_01 = "PKModelCovariance_uncCorr_VcCL_01.xlsx",
    VcCL_UncCorr_09 = "PKModelCovariance_uncCorr_VcCL_09.xlsx"
   # CLkabs_UncCorr_01 = "PKModelCovariance_uncCorr_CLkabs_01.xlsx",
  #  CLkabs_UncCorr_09 = "PKModelCovariance_uncCorr_CLkabs_09.xlsx"
  )

# Npop <- c(200)
# Ntests <- c(5)
#objAlt <- c(CLVc_corr_inverted = "PKModelCovariance _CLVc_corr_inverted.xlsx",
# VcQ1_corr_inverted = "PKModelCovariance _VcQ1_corr_inverted.xlsx",
# CLkabs_UncCorr_inverted = "PKModelCovariance_uncCorr_CLkabs_inverted",
# CLQ1_corr_inverted = "PKModelCovariance _CLQ1_corr_inverted.xlsx") # error


Npop <- c(100, 200, 500, 1000)
Ntests <- c(3, 5, 10, 15)

# test all case studies
out <- do.call(rbind, lapply(seq_along(objAlt), function (k) {
  do.call(rbind, lapply(seq_along(Ntests), function (i_test) {
    do.call(rbind, lapply(seq_along(Npop), function (j_pop) {
      testSample <-
        sampleParamFromUncertainty(objAlt[k], Npop = Npop[j_pop] * Ntests[i_test])
      t <-
        detectSampleDiscrepancy(
          obj = obj,
          testSample = testSample,
          Npop = Npop[j_pop],
          Ntests = Ntests[i_test]
        )

      data.frame(
        Change = names(objAlt[k]),
        nSingleDetected = length(t$detectedSingleParams),
        singleDetected = ifelse(
          length(t$detectedSingleParams) > 0,
          paste(t$detectedSingleParams, collapse = ','),NA),
        nPairsDetected = length(t$detectedPairParams),
        Npop = Npop[j_pop],
        Ntests = Ntests[i_test]
      )

    }))
  }))
}))

# output file
write.csv(out,file="UseCases.csv")

levels(out$Change)

Labs <- c("corr(CL,Q1)*0.1","corr(CL,Q1)*0.9","beta_CL(WT0)*0.1","beta_CL(WT0)*0.9","CL*0.1","CL*0.9","CL.RSE*1.1","CL.RSE*2","None",
          "UncCorr_Vc,CL*0.1","UncCorr_Vc,CL*0.9")
levels(out$Change)<-Labs

# Plot
library(ggplot2)
# Define on which page to plot
OutGroups <- rep(c(3,rep(c(1,2),5)),each=16)
out$OutGroup <- OutGroups

for (i in 1:2) {
  DataPlot <- out[out$OutGroup == i | out$OutGroup == 3 , ]
  p <- ggplot(DataPlot) +
    geom_line(
      aes(
        y = DataPlot$nSingleDetected,
        x = DataPlot$Npop,
        color = DataPlot$Change,
       # linetype = DataPlot$Change
      ),
      alpha=0.3
    ) +
    geom_point(aes(
      y = DataPlot$nSingleDetected,
      x = DataPlot$Npop,
      color = DataPlot$Change,
      shape =  DataPlot$Change
    )) +
    facet_wrap(DataPlot$Ntests) +
    labs(x = "Number of populations", y = "single parameters detected", color =
           "Introduced change", title = "Stratified by number of tests", shape = "Introduced change") +
    scale_color_IQRtools()

  pdf(paste0("useCases_singleParameters_", i, ".pdf"))
  print(p)
  dev.off()

}

# Pair parameters
for (i in 1:2) {
  DataPlot <- out[out$OutGroup == i | out$OutGroup == 3 , ]
  p <- ggplot(DataPlot) +
    geom_line(
      aes(
        y = DataPlot$nPairsDetected,
        x = DataPlot$Npop,
        color = DataPlot$Change
        #linetype = DataPlot$Change
      ),
      alpha=0.3
    ) +
    geom_point(aes(
      y = DataPlot$nPairsDetected,
      x = DataPlot$Npop,
      color = DataPlot$Change,
      shape =  DataPlot$Change
    )) +
    facet_wrap(DataPlot$Ntests) +
    labs(x = "Number of populations", y = "Paired parameters detected", color =
           "Introduced change", title = "Stratified by number of tests", shape = "Introduced change") +
    scale_color_IQRtools()

  pdf(paste0("useCases_PairedParameters_", i, ".pdf"))
  print(p)
  dev.off()

}


