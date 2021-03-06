% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibull_rmst.R
\name{weibull_rmst}
\alias{weibull_rmst}
\title{RMST using Weibull fit}
\usage{
weibull_rmst(num_boots = 1000, dat, tau, alpha, find_pval = FALSE, seed = NULL)
}
\arguments{
\item{num_boots}{Number of bootstrap iterations}

\item{dat}{Data frame of time-to-event data which MUST have the columns
'time', 'arm', and 'status}

\item{tau}{How long of a follow-up to consider, i.e. we integrate the survival
functions from 0 to tau}

\item{alpha}{Confidence interval is given for (alpha/2, 1-alpha/2) percentiles}

\item{find_pval}{Boolean, if TRUE then does bootstrap under the null to find p-value of mean difference and RMST difference}

\item{seed}{For reproducibility}
}
\value{
A list including out_tab (estimate and CI in both arms), trt_rmst,
pbo_rmst, diff_rmst, trt_CI, pbo_CI, diff_CI. Assumes trt coded as arm 1 and
placebo coded as arm 0.
}
\description{
RMST for time-to-event data under parametric Weibull fit for data in each
arm separately. Also can provide CI for RMST estimate and difference in RMST.
}
\examples{
time <- rexp(100)
status <- rbinom(n=100, prob=0.5, size=1)
arm <- c( rep(1, 50), rep(0, 50))
dat <- data.frame(time=time, status=status, arm=arm)
weibull_rmst(dat=dat, tau=1, alpha=0.05, num_boots=200)

}
