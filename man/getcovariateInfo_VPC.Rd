% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07-ModelAssessment.R
\name{getcovariateInfo_VPC}
\alias{getcovariateInfo_VPC}
\title{Get covariate information from all models and the provided data - function to be called by
`CreateDataVPC_IQRnlmeProject_MMV()`}
\usage{
getcovariateInfo_VPC(modelsSample, modelSampleSimulate, dataVPC)
}
\arguments{
\item{modelsSample}{:If sequential model building has been done (e.g. PKPD) model
parameters might need to be sampled from additional IQRnlmeProjects. modelSample
can be a vector with paths to these additional IQRnlmeProjects. It is not allowed that
in these projects and modelSampleSimulate estimated parameter names overlap}

\item{dataVPC}{: Dataset for VPC generation, or called from `getData_IQRnlmeProject( )`}
}
\description{
Get covariate information from all models and the provided data - function to be called by
`CreateDataVPC_IQRnlmeProject_MMV()`
}
