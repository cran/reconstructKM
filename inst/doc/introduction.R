## ---- include = FALSE---------------------------------------------------------
library(reconstructKM)
library(survival)
library(survminer)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
  #knitr::opts_chunk$set(fig.width=4, fig.height=3)
)

## ----setup, echo=FALSE--------------------------------------------------------
library(reconstructKM)

## ----out.width = "700px", echo=FALSE------------------------------------------
knitr::include_graphics(("pembro_overall.png"))

## ----out.width = "700px", echo=FALSE------------------------------------------
knitr::include_graphics(("pembro_first_clicks.png"))
knitr::include_graphics(("pembro_last_clicks.png"))

## ----NAR example, eval=TRUE, results='show', echo=TRUE, warning=FALSE---------

# define the NAR
pembro_NAR <- data.frame(time=seq(from=0, to=21, by=3), NAR=c(410, 377, 347, 278,  163, 71, 18, 0))
pbo_NAR <- data.frame(time=seq(from=0, to=21, by=3), NAR=c(206, 183, 149, 104, 59, 25, 8, 0))


## ----reconstruct example, eval=TRUE, results='show', echo=TRUE, warning=FALSE----
# Here I am loading some example data, but you will have to load your own data, for example by modifying
# the commented code below (to use your own file structure and file names).
# setwd("/users/rsun3/desktop")
# pembro_clicks <- read.csv("pembro_clicks.csv")

# load example data
data("pembro_clicks")
data("pembro_NAR")
data("pbo_clicks")
data("pbo_NAR")

# call format_raw_tabs() with the clicks table and NAR table
pembro_aug <- format_raw_tabs(raw_NAR=pembro_NAR,
                                  raw_surv=pembro_clicks) 
pbo_aug <- format_raw_tabs(raw_NAR=pbo_NAR,
                                  raw_surv=pbo_clicks) 

# reconstruct by calling KM_reconstruct()
pembro_recon <- KM_reconstruct(aug_NAR=pembro_aug$aug_NAR, aug_surv=pembro_aug$aug_surv)
pbo_recon <- KM_reconstruct(aug_NAR=pbo_aug$aug_NAR, aug_surv=pbo_aug$aug_surv)

# put the treatment and control arms into one dataset
pembro_IPD <- data.frame(arm=1, time=pembro_recon$IPD_time, status=pembro_recon$IPD_event)
pbo_IPD <- data.frame(arm=0, time=pbo_recon$IPD_time, status=pbo_recon$IPD_event)
allIPD <- rbind(pembro_IPD, pbo_IPD)


## ----Plot Pembro, eval=TRUE, echo=TRUE,  fig.align="center"-------------------
# plot
pembro_KM_fit <- survival::survfit(survival::Surv(time, status) ~ arm, data=allIPD)

pembro_KM <- survminer::ggsurvplot(pembro_KM_fit, data = allIPD, risk.table = TRUE, 
                        palette=c('red', 'blue'),
           legend=c(0.35,0.25), legend.title='',legend.labs=c('Pembrolizumab', 'Placebo'),
           title='Overall Survival',
           ylab='Survival (%)', xlab='Time (Mo)',
           tables.y.text=TRUE,
           tables.y.text.col=FALSE, risk.table.title='Number at Risk', break.time.by=6)
pembro_KM$plot       

## ----Cox Pembro, eval=TRUE, results='show', echo=TRUE, warning=FALSE----------

pembroCox <- coxph(Surv(time, status) ~ arm, data=allIPD)
print_cox_outputs(pembroCox)


## ----RMST Pembro, eval=TRUE, results='show', echo=TRUE, warning=FALSE---------

pembroRMST <- nonparam_rmst(dat = allIPD, tau = 18, alpha = 0.05) 
pembroRMST


## ----RMST weibull Pembro, eval=TRUE, results='show', cache=TRUE, echo=TRUE, warning=FALSE----

weibullFit <- weibull_rmst(num_boots=1000, dat=allIPD, tau=60, alpha=0.05, find_pval=FALSE, seed=NULL)
weibullFit$rmst_df

