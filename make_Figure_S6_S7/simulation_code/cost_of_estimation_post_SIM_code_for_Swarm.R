library(xtable)
library(mvmeta)
library(mvtnorm)
library(Matrix)
library(speedglm)
library(clusterGeneration)

setwd("/data/bocasm/Research/Meta\ analysis/code/simulations/cost_of_estimation_post_SIM")

tryCatch.W.E <- function(expr)
{
    W <- NULL
    w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
         warning = w.handler),
         warning = W)
}

##simulated scenario number
randVar <- as.numeric(commandArgs()[4])
##number of parameters
p <- as.numeric(commandArgs()[5])
##heterogeneity (i.e. Sigma^2/S^2)
het <- as.numeric(commandArgs()[6])
##batch number
batch <- as.numeric(commandArgs()[7])
##whether sample size if the same for all studies or different
sampleType <- commandArgs()[8]
##number of studies
I <- as.numeric(commandArgs()[9])

print(c(randVar, p, het, batch, sampleType, I))

##number of simulations in the batch
nrSims <- 2500
##for the higher values of p, do fewer simulations (will do more batches for them)
if(p == 15)
{
    nrSims <- 1250
}
if(p == 20)
{
    nrSims <- 500
}

##variance between studies (set it equal to het)
##between-study varriance
varBtw <- het
##number of people in each study
if(sampleType == "diff")
{
    if(I == 10)
    {
        nCases <- (5:14)/2*100+25
        nCtrls <- (5:14)/2*100+25
    } else {
        nCases <- (430+(0:19)*60)/2
        nCtrls <- (430+(0:19)*60)/2
    }
} else {
    nCases <- rep(500, I)
    nCtrls <- rep(500, I)
}
n <- nCases+nCtrls

Beta0 <- -0.125*p
Beta0 <- rep(Beta0, I)

OverallMeans <- rep(0.2, p)
OverallMeansMat <- matrix(OverallMeans, nrow=1)

###################################################
#######  Simulate based on this data  #######
###################################################

print(date())

##set the seed in terms of randVar, so that each scenario has different variance-covariance matrices,
##which are preserved through the different batches for that scenario
set.seed(381048+1940*randVar)

##make between-study variance-covariance matrix
Sigma <- cov2cor(genPositiveDefMat(p)$Sigma)*varBtw

##make variance-covariance matrices for all sites
varCovsSitesAll <- list()
for(site in 1:I)
  {
    varCovsSitesAll[[site]] <- 4/mean(n)*solve(cov2cor(genPositiveDefMat(p)$Sigma))
  }

##set a batch-specific seed for the individual-level simulations as well
set.seed(batch*300+313)

##save prevalences for each site
simPrevSite <- matrix(NA, nrow=nrSims, ncol=I)

#save estimates for the study-specific means
StudyEst <- list()
##save estimates for the within-study variances
StudyEstVars <- list()
##for each site, save whether or not there is a warning for the logistic regression
StudyWarns <- matrix("", nrow=nrSims, ncol=I)

##save estimated Sigma
SigmaEst <- list()
##save estimated Sigma from univariate meta-analysis (so only diagonal entries)
SigmaEstUniv <- list()
##save warnings due to mvmeta
MetaWarns <- rep("", nrSims)

print(date())

for(sim in 1:nrSims)
{
  if(sim %% 10 == 0)
  {
    print(sim)
  }

  ##simulate the study-specific means
  studMeans <- rmvnorm(I, mean = OverallMeans,
                       sigma = Sigma, method="svd")

  ##for this simulation, get estimated means and variances
  estMeansSim <- matrix(rep(0, p),
                        nrow = I, ncol = p)

  estVarsSim <- list()

  ##save disease status for everyone, to get overall prevalence
  diseaseAllSites <- c()
  ##get site-specific prevalence
  prevSite <- rep(0, I)
  for(site in 1:I)
  {
      OverallMeansMat <- studMeans[site, , drop=FALSE]

      meansSite <- rep(0, p)

      varCovsSite <- varCovsSitesAll[[site]]

      ##simulate values for the p variables from the multivariate normal for each dataset
      ##simulate 10 times the population size, since will be doing retrospective sampling
      XSimSite <- rmvnorm(n[site]*10, mean = meansSite, sigma = varCovsSite, method="svd")

      ##get estimated probabilities
      logits <- Beta0[site]+.Internal(rowSums(OverallMeansMat[rep(1, nrow(XSimSite)), ] * XSimSite,
                                              n[site]*10, p, TRUE))

      expLogits <- exp(logits)
      probs <- expLogits/(1+expLogits)
      probs <- probs[!is.na(probs)]

      ##get whether the individual is a case or control
      CC <- rbinom(length(probs), 1, probs)

      diseaseAllSites <- c(diseaseAllSites, CC)
      prevSite[site] <- mean(CC)

      ##get cases at this site
      casesSite <- which(CC == 1)
      ctrlsSite <- which(CC == 0)

      ##sample number of cases and controls from the population generated
      casesSite <- sample(casesSite, nCases[site])
      ctrlsSite <- sample(ctrlsSite, nCtrls[site])

      simPrevSite[sim, site] <-
          length(casesSite)/(length(casesSite)+length(ctrlsSite))

      dataSiteSamp <- cbind(c(rep(1, length(casesSite)),
                              rep(0, length(ctrlsSite))),
                            XSimSite[c(casesSite, ctrlsSite), ])

      tmpX <- cbind(rep(1, n[site]), dataSiteSamp[, 2:(p+1)])
      colnames(tmpX) <- 0:p

      ##get estimates, while also catching warnings
      logistSiteCatchWarn <- tryCatch.W.E(speedglm.wfit(y = dataSiteSamp[, 1],
                                                        X = tmpX,
                                                        family = binomial(link="logit")))

      logistSite <- logistSiteCatchWarn$value
      StudyWarns[sim, site] <- paste(logistSiteCatchWarn$warning, collapse="\n")

      estMeansSim[site, ] <- logistSite$coefficients[2:(p+1)] ##coefficients(summary(logistSite))[-1, "Estimate"]
      estVarsSim[[site]] <- summary(logistSite)$cov.scaled[2:(p+1), 2:(p+1)] ##vcov(logistSite)[-1,-1]
  }

  ##save estimated effect sizes
  StudyEst[[sim]] <- estMeansSim
  ##save estimated within-study variance-covariance matrices
  StudyEstVars[[sim]] <- estVarsSim

  ##get estimates of Sigma, while also catching warnings
  SigmaEstCatchWarn <- tryCatch.W.E(mvmeta(StudyEst[[sim]],
                                           StudyEstVars[[sim]],
                                           method="reml")$Psi)

  SigmaEst[[sim]] <- SigmaEstCatchWarn$value
  MetaWarns[sim] <- paste(SigmaEstCatchWarn$warning, collapse="\n")

  ##also get estimates of Sigma from univariate meta-analysis
  SigmaEstUniv.sim <- rep(NA, p)
  for(pp in 1:p)
  {
      SigmaEstUniv.sim[pp] <- mvmeta(StudyEst[[sim]][,pp],
                                     S =
                                     lapply(StudyEstVars[[sim]],
                                            function(x){x[pp, pp, drop=FALSE]}),
                                     method="reml")$Psi
  }

  SigmaEstUniv[[sim]] <- SigmaEstUniv.sim

  ##overall prevalence
  prev <- mean(diseaseAllSites)

  if(sim %% 50 == 0)
  {
      print(prev)
  }
}

print(date())

StudyEstVars[[1]][[1]]
StudyEstVars[[1]][[10]]

save(list=c("StudyEst", "StudyEstVars", "SigmaEst", "SigmaEstUniv", "StudyWarns", "MetaWarns",
     "Sigma", "varCovsSitesAll"),
     file = paste("simResults/cost_of_estimation_post_SIM",
     "randVar", randVar, "p", p, "het", het,
     "batch", batch,
     "sampleType", sampleType,
     "I", I,
     ".RData", sep="_"))

rm(list=ls())

