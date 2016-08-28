library(xtable)
library(mvmeta)
library(mvtnorm)
library(Matrix)
library(speedglm)
library(rjags)

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

##make block diagonal matrix giving correlation, size of block, number of blocks
blockDiag <- function(rho, sizeBlock, nrBlocks)
{
  block <- matrix(rho, sizeBlock, sizeBlock)
  diag(block) <- 1
  blockList <- list()
  for(i in 1:nrBlocks)
  {
    blockList[[i]] <- block
  }
  Sigma <- bdiag(blockList)
  Sigma <- as.matrix(Sigma)
  Sigma
}

##number of studies
I <- as.numeric(commandArgs()[4])
##number of parameters
p <- as.numeric(commandArgs()[5])
##between-study correlation
corrBtw <- as.numeric(commandArgs()[6])
##between-study variance
varBtw <- as.numeric(commandArgs()[7])
##batch number
batch <- as.numeric(commandArgs()[8])
##number of people in the study
nCases <- 500
nCtrls <- 500
n <- nCases+nCtrls
    
##
##I <- 10
##p <- 7
##corrBtw <- 0.5
##varBtw <- 1
##batch <- 10

nrBlocks <- round(p/5)

##have one "incomplete block" if p is not divisible by 5
if(p %% 5 != 0)
  {
    nrBlocks <- nrBlocks + 1
  }

##make between-study var-cov matrix
##Sigma <- matrix(corrBtw, nrow=p, ncol=p)
##diag(Sigma) <- 1
Sigma <- blockDiag(corrBtw, sizeBlock=5, nrBlocks=nrBlocks)[1:p,1:p]
Sigma <- Sigma*varBtw

nrSims <- 2500

##nrSims <- 20

setwd("/data/bocasm/Research/Meta\ analysis/code/simulations/cost_of_estimation_bd")
Beta0 <- -0.125*p
Beta0 <- rep(Beta0, I)

OverallMeans <- rep(0.2, p)
OverallMeansMat <- matrix(OverallMeans, nrow=1)

###################################################
#######  Simulate based on this data  #######
###################################################

print(date())

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
##similarly, save estimated Sigmas from the Bayesian analyses
SigmaEstBayes <- SigmaEstUnivBayes <- list()
##save warnings due to mvmeta
MetaWarns <- rep("", nrSims)
##save estimates of mu from Bayesian approach
OverallEstBayes <- OverallEstUnivBayes <- matrix(NA, nrow=nrSims, ncol=p)

print(date())

#get variance-covariance matrices for all site
varCovsSitesAll <- list()
for(site in 1:I)
{
    ##varCovsSite <- matrix((site-1)/I, p, p)
    ##diag(varCovsSite) <- rep(1, p)
    ##make a block diagonal matrix of size p x p
    varCovsSite <- blockDiag((site-1)/I, sizeBlock=5, nrBlocks=nrBlocks)[1:p,1:p]
    varCovsSitesAll[[site]] <- 4/n * solve(varCovsSite)
}

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
      XSimSite <- rmvnorm(n*10, mean = meansSite, sigma = varCovsSite, method="svd")

      ##get estimated probabilities
      logits <- Beta0[site]+.Internal(rowSums(OverallMeansMat[rep(1, nrow(XSimSite)), ] * XSimSite,
                                              n*10, p, TRUE))

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
      casesSite <- sample(casesSite, nCases)
      ctrlsSite <- sample(ctrlsSite, nCtrls)

      simPrevSite[sim, site] <-
          length(casesSite)/(length(casesSite)+length(ctrlsSite))

      dataSiteSamp <- cbind(c(rep(1, length(casesSite)),
                              rep(0, length(ctrlsSite))),
                            XSimSite[c(casesSite, ctrlsSite), ])

      tmpX <- cbind(rep(1, n), dataSiteSamp[, 2:(p+1)])
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

##now do the Bayesian analysis!
##first get the *precision* matrices and concatenate them!
    Tau <- c()
    for(i in 1:I)
    {
      Tau <- rbind(Tau, solve(StudyEstVars[[sim]][[i]]))
    }
##also do this in the univariate case! (just take the inverse of the variances of the first component only)
    TauUniv <- matrix(NA, ncol=1)
    for(i in 1:I)
    {
      TauUniv[i] <- 1/StudyEstVars[[sim]][[i]][1,1]
    }

    ###starting position for each study in Sigma2
    pos      <- 1+c(0:(I-1))*p
    ###ending position for each study in Sigma2
    pos2     <- pos+p-1
    ###Omega2 is a fixed hyperparameter - it's a prior on the between-study *precision matrix*
    Omega <- diag(1/runif(p,0,2)^2)
    
    ##create model for multivariate meta-analysis
    jagsMultiv <- jags.model('example_MVMA.bug',
                             ##'H:/Simina/metaAnalysis/example6.bug.txt',
                             data = list('p' = p,
                                         'N' = I,
                                         'y'=StudyEst[[sim]],
                                         'pos'=pos,
                                         'sigma'=Tau,
                                         'pos2'=pos2,
                                         'Omega'=Omega),
                             n.chains = 1,
                             n.adapt = 1000,
                             quiet=TRUE)

    jagsMultivRes <- jags.samples(jagsMultiv,
                                  c('intau','mu'),
                                  50000,
                                  thin=5)
##stop()

    SigmaEstBayes[[sim]] <- solve(summary(jagsMultivRes$intau, FUN="median")[1]$stat)
    OverallEstBayes[sim, ] <- summary(jagsMultivRes$mu, FUN="median")[1]$stat

##just save the estimate for the first variance with the univariate model - don't fit it for all p outcomes in this case
   jagsUniv <- jags.model('example_UVMA.bug',
                           data = list('N' = I,
                                       'y'=StudyEst[[sim]][,1],
                                       'pos'=1:I,
                                       'sigma'=TauUniv,
                                       'pos2'=1:I,
                                       'Omega'=Omega[1,1]),
                           n.chains = 1,
                           n.adapt = 1000,
                           quiet=TRUE)

    jagsUnivRes <- jags.samples(jagsUniv,
                                  c('intau','mu'),
                                  50000,
                                  thin=5)

    SigmaEstUnivBayes[[sim]] <- 1/summary(jagsUnivRes$intau, FUN="median")[1]$stat
    OverallEstUnivBayes[sim, 1] <- summary(jagsUnivRes$mu, FUN="median")[1]$stat

  if(sim %% 50 == 0)
  {
      print(prev)
  }
}

print(date())

save(list=c("StudyEst", "StudyEstVars", 
"SigmaEst", "SigmaEstUniv", 
"SigmaEstBayes", "SigmaEstUnivBayes",
"OverallEstBayes", "OverallEstUnivBayes",
"StudyWarns", "MetaWarns"),
     file = paste("simResults/cost_of_estimation",
     "I", I, "p", p, "corrBtw", corrBtw, "varBtw", varBtw,
     "batch", batch,
     ".RData", sep="_"))

##stop()

rm(list=ls())
