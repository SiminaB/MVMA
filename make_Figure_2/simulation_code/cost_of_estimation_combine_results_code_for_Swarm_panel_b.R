library(mvmeta)

setwd("/data/bocasm/Research/Meta analysis/code/simulations/cost_of_estimation_Suppl")

allFiles <- list.files("simResults")

##save the following:
##MSEs for multivariate meta-analysis with known variance-covariance matrices
##MSEs for multivariate meta-analysis with unknown (estimated) variance-covariance matrices
##MSEs for univariate meta-analysis
##variances for the same things...

##save the two MSEs for the different combinations
MSEs <- expand.grid(I = 20,
                    p = 2:10,
                    corrBtw = c(0.5),
                    varBtw = c(1/5, 1, 5))

MSEs$known <- MSEs$unknown <- MSEs$Bayes <- 
MSEs$univKnown <- MSEs$univ <- MSEs$univBayes <- NA

##do this for the variances as well
Vars <- MSEs

##do this for one row at a time
rowNr <- as.numeric(commandArgs()[4])

##

I <- MSEs[rowNr, 1]
p <- MSEs[rowNr, 2]
corrBtw <- MSEs[rowNr, 3]
varBtw <- MSEs[rowNr, 4]

OverallMeans <- rep(0.2, p)

##get all the files with this I and p
filesOfType <- allFiles[grep(paste("I", I, "p", p,"corrBtw", corrBtw, "varBtw", varBtw, sep="_"), allFiles)]
print(paste("I=",I," p=",p," corrBtw=",corrBtw," varBtw=",varBtw," number of files=",
            length(filesOfType), sep=""))

##get Sigma for these values
Sigma <- matrix(corrBtw, nrow=p, ncol=p)
diag(Sigma) <- 1
Sigma <- Sigma*varBtw

##save all study-specific effect size estimates and variance estimates
StudyEstAll <- StudyEstVarsAll <- SigmaEstAll <- SigmaEstUnivAll <- list()
StudyWarnsAll <- MetaWarnsAll <- c()

##added on August 7, 2016: also save Bayesian estimates
OverallEstBayesAll <- OverallEstUnivBayesAll <- c()

for(file in 1:length(filesOfType))
{
    load(paste("simResults", filesOfType[file], sep="/"))

    nrSims <- length(StudyEst)

    ##print(file)
    ##print(range(sapply(StudyEst, length)))

    StudyEstAll[((file-1)*nrSims+1):(file*nrSims)] <- StudyEst
    StudyEstVarsAll[((file-1)*nrSims+1):(file*nrSims)] <- StudyEstVars
    SigmaEstAll[((file-1)*nrSims+1):(file*nrSims)] <- SigmaEst
    SigmaEstUnivAll[((file-1)*nrSims+1):(file*nrSims)] <- SigmaEstUniv

    ##print(range(sapply(StudyEstAll, length)))
    ##print(length(StudyEstAll))

    StudyWarnsAll <- rbind(StudyWarnsAll, StudyWarns)
    MetaWarnsAll <- c(MetaWarnsAll, MetaWarns)

    OverallEstBayesAll <- rbind(OverallEstBayesAll, OverallEstBayes)
    OverallEstUnivBayesAll <- rbind(OverallEstUnivBayesAll, OverallEstUnivBayes)
}
nrSimsTot <- length(StudyEstAll)

##save empirical variance estimates for within-study variance of effect size estimates
StudyEmpVars <- list()
for(i in 1:I)
{
    StudyEmpVars[[i]] <- var(t(sapply(StudyEstAll, function(x){x[i,]})))
}

StudyEmpVarsInv <- lapply(StudyEmpVars, solve)

##transpose the study estimates
StudyEstAll <- lapply(StudyEstAll, t)

##save meta-analytic estimates for known Sigma and for unknown (estimated) Sigma and for univariate analysis
##(assume that the within-study correlations are estimated, either way)
OverallEstSigmaKnown <- OverallEstSigmaUnknown <- OverallEstSigmaUnivKnown <- OverallEstSigmaUniv <-
    matrix(NA, nrow=nrSimsTot, ncol=p)

for(sim in 1:nrSimsTot)
{
    if(sim %% 10000 == 0)
    {
        print(sim)
    }

    StudyEstVars.sim <- StudyEstVarsAll[[sim]]

    StudyEst.sim <- StudyEstAll[[sim]]

    SigmaEst.sim <- SigmaEstAll[[sim]]
    SigmaEstUniv.sim <- SigmaEstUnivAll[[sim]]

    Y1 <- StudyEst.sim[,1,drop=FALSE]

    WiUnknownInv <- lapply(StudyEstVars.sim, function(x,s){solve(x+s)},
                           s=SigmaEst.sim)
    WiKnownInv <- lapply(StudyEstVars.sim, function(x,s){solve(x+s)},
                         s=Sigma)

    Term1u <- WiUnknownInv[[1]]
    Term2u <- WiUnknownInv[[1]] %*% Y1

    Term1k <- WiKnownInv[[1]]
    Term2k <- WiKnownInv[[1]] %*% Y1

    for(i in 2:I)
    {
        Yi <- StudyEst.sim[,i,drop=FALSE]

        Term1u <- Term1u + WiUnknownInv[[i]]
        Term2u <- Term2u + WiUnknownInv[[i]] %*% Yi

        Term1k <- Term1k + WiKnownInv[[i]]
        Term2k <- Term2k + WiKnownInv[[i]] %*% Yi
    }

    OverallEstSigmaUnknown[sim, ] <- ##(solve(tX %*% solve(SiEst.sim) %*% X) %*% tX %*% solve(SiEst.sim) %*% Y)[,1]
        (solve(Term1u) %*% Term2u)[,1]
    OverallEstSigmaKnown[sim, ] <- ##(solve.Term1k %*% tX %*% solve(Si) %*% Y)[,1]
        (solve(Term1k) %*% Term2k)[,1]

    for(pp in 1:p)
    {
        ##get the estimated within-study variances for parameter pp for all studies
        StudyEstVars.sim.pp <- sapply(StudyEstVars.sim,
                                      function(x){x[pp, pp]})

        OverallEstSigmaUniv[sim, pp] <-
            sum(StudyEst.sim[pp,]*1/(SigmaEstUniv.sim[pp]+StudyEstVars.sim.pp))/
                sum(1/(SigmaEstUniv.sim[pp]+StudyEstVars.sim.pp))

        OverallEstSigmaUnivKnown[sim, pp] <-
            sum(StudyEst.sim[pp,]*1/(Sigma[pp,pp]+StudyEstVars.sim.pp))/
                sum(1/(Sigma[pp,pp]+StudyEstVars.sim.pp))

    }
}

MSEs[rowNr, c("known", "unknown", "Bayes", "univKnown", "univ", "univBayes")] <-
    c(mean((OverallEstSigmaKnown[,1]-OverallMeans[1])^2),
      mean((OverallEstSigmaUnknown[,1]-OverallMeans[1])^2),
      mean((OverallEstBayesAll[,1]-OverallMeans[1])^2),
      mean((OverallEstSigmaUnivKnown[,1]-OverallMeans[1])^2),
      mean((OverallEstSigmaUniv[,1]-OverallMeans[1])^2),
      mean((OverallEstUnivBayesAll[,1]-OverallMeans[1])^2))

Vars[rowNr, c("known", "unknown", "Bayes", "univKnown", "univ", "univBayes")] <-
   c(var(OverallEstSigmaKnown[,1]),
     var(OverallEstSigmaUnknown[,1]),
     var(OverallEstBayesAll[,1]),
     var(OverallEstSigmaUnivKnown[,1]),
     var(OverallEstSigmaUniv[,1]),
     var(OverallEstUnivBayesAll[,1]))

###################################

head(OverallEstSigmaUnknown)
head(OverallEstSigmaKnown)
head(OverallEstSigmaUniv)

MSEs[rowNr,]
Vars[rowNr,]

##save results
save(list=c("MSEs", "Vars",
     "OverallEstSigmaUnknown", "OverallEstSigmaKnown",
     "OverallEstSigmaUniv"),
     file=paste("simResultsComb/combine_cost_of_estimation_Suppl_Bayes",rowNr,".RData", sep=""))

rm(list=ls())
