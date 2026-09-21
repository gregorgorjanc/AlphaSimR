# fmt: skip file

#' @title Mean genetic values
#'
#' @description Returns the mean genetic values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanG(pop)
#'
#' @export
meanG = function(pop){
  colMeans(pop@gv)
}

#' @title Mean phenotypic values
#'
#' @description Returns the mean phenotypic values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanP(pop)
#'
#' @export
meanP = function(pop){
  colMeans(pop@pheno)
}

#' @title Mean estimated breeding values
#'
#' @description Returns the mean estimated breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
#' meanEBV(pop)
#'
#' @export
meanEBV = function(pop){
  colMeans(pop@ebv)
}

#' @title Total genetic variance
#'
#' @description Returns total genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#'
#' @export
varG = function(pop){
  G = popVar(pop@gv)
  rownames(G) = colnames(G) = colnames(pop@gv)
  return(G)
}

#' @title Phenotypic variance
#'
#' @description Returns phenotypic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varP(pop)
#'
#' @export
varP = function(pop){
  P = popVar(pop@pheno)
  rownames(P) = colnames(P) = colnames(pop@pheno)
  return(P)
}

#' @title Variance of estimated breeding values
#'
#' @description Returns variance of estimated breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
#' varA(pop)
#' varEBV(pop)
#'
#' @export
varEBV = function(pop){
  ebv = popVar(pop@ebv)
  rownames(ebv) = colnames(ebv) = colnames(pop@ebv)
  return(ebv)
}

#' @title Calculate quantitative genetic quantities
#'
#' @description
#' Calculates quantitative genetic quantities and their variances
#' for an object of \code{\link{Pop-class}}
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @return
#' \describe{
#' \item{varG}{an nTrait by nTrait matrix of total genetic variances}
#' \item{varA}{an nTrait by nTrait matrix of additive genetic variances}
#' \item{varD}{an nTrait by nTrait matrix of dominance genetic variances}
#' \item{varAA}{an nTrait by nTrait matrix of additive-by-additive genetic variances}
#' \item{varN}{an nTrait by nTrait matrix of non-additive genetic variances}
#' \item{genicVarG}{an nTrait vector of total genic variances}
#' \item{genicVarA}{an nTrait vector of additive genic variances}
#' \item{genicVarD}{an nTrait vector of dominance genic variances}
#' \item{genicVarAA}{an nTrait vector of additive-by-additive genic variances}
#' \item{genicVarN}{an nTrait vector of non-additive genic variances}
#' \item{covG_HW}{an nTrait vector of total genicTODO covariances due to non-random mating}
#' \item{covA_HW}{an nTrait vector of additive TODO covariances due to non-random mating}
#' \item{covD_HW}{an nTrait vector of dominance TODO covariances due to non-random mating}
#' \item{covAA_HW}{an nTrait vector of additive-by-additive TODO covariances due to non-random mating}
#' \item{covN_HW}{an nTrait vector of non-additive TODO covariances due to non-random mating}
#' \item{covG_L}{an nTrait vector of total genic TODO covariances due to linkage disequilibrium}
#' \item{covA_L}{an nTrait vector of additive TODOcovariances due to linkage disequilibrium}
#' \item{covD_L}{an nTrait vector of dominance TODO covariances due to linkage disequilibrium}
#' \item{covAA_L}{an nTrait vector of additive-by-additive TODO covariances due to linkage disequilibrium}
#' \item{covAD_L}{an nTrait vector of additive by dominance TODO covariances due to linkage disequilibrium}
#' \item{covAAA_L}{an nTrait vector of additive by additive-by-additive TODO covariances due to linkage disequilibrium}
#' \item{covDAA_L}{an nTrait vector of dominance by additive-by-additive TODO covariances due to linkage disequilibrium}
#' \item{covAN_L}{an nTrait vector of additive by non-additive TODO covariances due to linkage disequilibrium}
#' \item{mu}{an nTrait vector of trait means}
#' \item{mu_HW}{an nTrait vector of expected trait means under random mating}
#' \item{gv}{a matrix of genetic values with dimensions nInd by nTraits}
#' \item{bv}{a matrix of breeding values with dimensions nInd by nTraits}
#' \item{dd}{a matrix of dominance deviations with dimensions nInd by nTraits}
#' \item{aa}{a matrix of additive-by-additive epistatic deviations with dimensions nInd by nTraits}
#' \item{nd}{a matrix of non-additive deviations with dimensions nInd by nTraits}
#' \item{gv_mu}{an nTrait TODO vector of genetic value means with dimensions nInd by nTraits TODO}
#' \item{gv_a}{a matrix of additive  genetic values with dimensions nInd by nTraits}
#' \item{gv_d}{a matrix of dominance TODO genetic values with dimensions nInd by nTraits}
#' \item{gv_aa}{a matrix of additive-by-additive TODO genetic values with dimensions nInd by nTraits}
#' \item{gv_n}{a matrix of non-additive TODO genetic values with dimensions nInd by nTraits}
#' \item{alpha}{a list of average allele substitution effects with length nTraits}
#' \item{alpha_HW}{a list of average allele substitution effects at Hardy-Weinberg equilibrium with length nTraits}
#' }
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genParam(pop, simParam=SP)
#'
#' @export
genParam = function(pop,simParam=NULL,nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  nInd = nInd(pop)
  nTraits = simParam$nTraits
  traitNames = simParam$traitNames

  # Blank nInd x nTrait matrices
  gv = matrix(NA_real_, nrow=nInd, ncol=nTraits)
  colnames(gv) = traitNames
  bv = dd = aa = nd = gv_a = gv_d = gv_aa = gv_n = gv

  # Blank nTrait vectors
  genicVarA = rep(NA_real_, nTraits)
  names(genicVarA) = traitNames
  genicVarD = genicVarAA = genicVarN =
    covG_HW = covA_HW = covD_HW = covAA_HW = covN_HW =
    covAD_L = covAAA_L = covDAA_L = covAN_L =
    mu = mu_HW = gv_mu = genicVarA

  # Average effect of an allele substitution
  alpha = vector("list", length=nTraits)
  names(alpha) = traitNames
  alpha_HW = alpha

  #Loop through trait calculations
  for(i in seq_len(nTraits)){
    trait = simParam$traits[[i]]
    tmp = calcGenParam(trait,pop,nThreads)
    genicVarA[i] = tmp$genicVarA2
    genicVarN[i] = 0
    covA_HW[i] = tmp$genicVarA-tmp$genicVarA2
    covN_HW[i] = 0
    gv[,i] = tmp$gv
    bv[,i] = tmp$bv
    nd[,i] = rep(0,pop@nInd)
    mu[i] = tmp$mu
    mu_HW[i] = tmp$mu_HWE
    gv_mu[i] = tmp$gv_mu
    gv_a[,i] = tmp$gv_a
    gv_n[,i] = rep(0,pop@nInd)
    if(.hasSlot(trait,"domEff")){
      genicVarD[i] = tmp$genicVarD2
      genicVarN[i] = genicVarN[i] + genicVarD[i]
      covD_HW[i] = tmp$genicVarD-tmp$genicVarD2
      covN_HW[i] = covN_HW[i] + covD_HW[i]
      dd[,i] = tmp$dd
      nd[,i] = nd[,i] + dd[,i]
      gv_d[,i] = tmp$gv_d
      gv_n[,i] = gv_n[,i] + gv_d[,i]
    }else{
      genicVarD[i] = 0
      covD_HW[i] = 0
      dd[,i] = rep(0,pop@nInd)
      gv_d[,i] = rep(0,pop@nInd)
    }
    if(.hasSlot(trait,"epiEff")){
      genicVarAA[i] = tmp$genicVarAA2
      genicVarN[i] = genicVarN[i] + genicVarAA[i]
      covAA_HW[i] = tmp$genicVarAA-tmp$genicVarAA2
      covN_HW[i] = covN_HW[i] + covAA_HW[i]
      aa[,i] = tmp$aa
      nd[,i] = nd[,i] + aa[,i]
      gv_aa[,i] = tmp$gv_aa
      gv_n[,i] = gv_n[,i] + gv_aa[,i]
    }else{
      genicVarAA[i] = 0
      covAA_HW[i] = 0
      aa[,i] = rep(0,pop@nInd)
      gv_aa[,i] = rep(0,pop@nInd)
    }
    if(nInd==1){
      covAD_L[i] = 0
      covAAA_L[i] = 0
      covDAA_L[i] = 0
      covAN_L[i] = 0
    } else {
      covAD_L[i] = popVar(cbind(bv[,i],dd[,i]))[1,2]
      covAAA_L[i] = popVar(cbind(bv[,i],aa[,i]))[1,2]
      covDAA_L[i] = popVar(cbind(dd[,i],aa[,i]))[1,2]
      covAN_L[i] = popVar(cbind(bv[,i],nd[,i]))[1,2]
    }
    alpha[[i]] = tmp$alpha
    alpha_HW[[i]] = tmp$alpha_HW
  }

  varG = popVar(gv)
  rownames(varG) = colnames(varG) = traitNames

  varA = popVar(bv)
  rownames(varA) = colnames(varA) = traitNames

  varD = popVar(dd)
  rownames(varD) = colnames(varD) = traitNames

  varAA = popVar(aa)
  rownames(varAA) = colnames(varAA) = traitNames

  varN = popVar(nd)
  rownames(varN) = colnames(varN) = traitNames

  genicVarG = genicVarA + genicVarD + genicVarAA
  covG_HW = covA_HW + covD_HW + covAA_HW

  output = list(varG=varG,
                varA=varA,
                varD=varD,
                varAA=varAA,
                varN=varN,
                genicVarG=genicVarG,
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarAA=genicVarAA,
                genicVarN=genicVarN,
                covG_HW=covG_HW,
                covA_HW=covA_HW,
                covD_HW=covD_HW,
                covAA_HW=covAA_HW,
                covN_HW=covN_HW,
                covG_L=diag(varG)-genicVarG-covG_HW,
                covA_L=diag(varA)-genicVarA-covA_HW,
                covD_L=diag(varD)-genicVarD-covD_HW,
                covAA_L=diag(varAA)-genicVarAA-covAA_HW,
                covN_L=diag(varN)-genicVarN-covN_HW,
                covAD_L=covAD_L,
                covAAA_L=covAAA_L,
                covDAA_L=covDAA_L,
                covAN_L=covAN_L,
                mu=mu,
                mu_HW=mu_HW,
                gv=gv,
                bv=bv,
                dd=dd,
                aa=aa,
                nd=nd,
                gv_mu=gv_mu,
                gv_a=gv_a,
                gv_d=gv_d,
                gv_aa=gv_aa,
                gv_n=gv_n,
                alpha=alpha,
                alpha_HW=alpha_HW)
  return(output)
}

#' @title Additive genetic variance
#'
#' @description Returns additive genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varA(pop, simParam=SP)
#'
#' @export
varA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varA
}

#' @title Dominance genetic variance
#'
#' @description Returns dominance genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varD(pop, simParam=SP)
#'
#' @export
varD = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varD
}

#' @title Additive-by-additive epistatic genetic variance
#'
#' @description Returns additive-by-additive epistatic genetic
#' variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varAA(pop, simParam=SP)
#'
#' @export
varAA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varAA
}

#' @title Non-additive genetic variance
#'
#' @description Returns non-additive genetic variance for all traits
#'   (includes dominance and epistatic genetic variance)
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop, simParam=SP)
#' varA(pop, simParam=SP)
#' varN(pop, simParam=SP)
#'
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop, simParam=SP)
#' varA(pop, simParam=SP)
#' varD(pop, simParam=SP)
#' varAA(pop, simParam=SP)
#' varN(pop, simParam=SP)
#'
#' @export
varN = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varN
}

#' @title Breeding value
#'
#' @description Returns breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' bv(pop, simParam=SP)
#'
#' @export
bv = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$bv
}

#' @title Dominance deviations
#'
#' @description Returns dominance deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' dd(pop, simParam=SP)
#'
#' @export
dd = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$dd
}

#' @title Additive-by-additive epistatic deviations
#'
#' @description Returns additive-by-additive epistatic
#' deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, rel=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' aa(pop, simParam=SP)
#'
#' @export
aa = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$aa
}

#' @title Non-additive deviations
#'
#' @description Returns non-additive deviations for all traits
#'   (includes dominance and epistatic deviations)
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' cbind(dd(pop, simParam=SP),
#'       aa(pop, simParam=SP),
#'       nd(pop, simParam=SP))
#'
#' @export
nd = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$nd
}

#' @title Total genic variance
#'
#' @description Returns total genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#'
#' @export
genicVarG = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarG
}

#' @title Additive genic variance
#'
#' @description Returns additive genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#'
#' @export
genicVarA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarA
}

#' @title Dominance genic variance
#'
#' @description Returns dominance genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#'
#' @export
genicVarD = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarD
}

#' @title Additive-by-additive genic variance
#'
#' @description Returns additive-by-additive epistatic
#' genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#'
#' @export
genicVarAA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarAA
}

#' @title Non-additive genic variance
#'
#' @description Returns non-additive genic variance for all traits
#'   (includes dominance and epistatic genic variance)
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#' genicVarN(pop, simParam=SP)
#'
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#' genicVarN(pop, simParam=SP)
#'
#' @export
genicVarN = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarN
}

#' @title Genetic value
#'
#' @description A wrapper for accessing the gv slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @details See vignette TODO for background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' gv(pop)
#'
#' @export
gv = function(pop){
  pop@gv
}

#' @title Phenotype
#'
#' @description A wrapper for accessing the pheno slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pheno(pop)
#'
#' @export
pheno = function(pop){
  pop@pheno
}

#' @title Estimated breeding value
#'
#' @description A wrapper for accessing the ebv slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
#' ebv(pop)
#'
#' @export
ebv = function(pop){
  pop@ebv
}

#' @title Calculate parent average
#'
#' @param pop \code{\link{Pop-class}} with individuals whose parent average
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"}, \code{"\link{bv}"},
#'   \code{"\link{ebv}"}, or \code{"\link{pheno}"}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return a matrix of parent averages with dimensions nInd by nTraits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#' parentAverage(pop2, parents = pop)
#' parentAverage(pop2, mothers = pop, fathers = pop)
#'
#' @export
parentAverage = function(pop, parents = NULL, mothers = NULL, fathers = NULL,
                         use = "gv", simParam = NULL, nThreads=NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if (!is.null(parents)) {
    matchMothers = match(x = pop@mother, table = parents@id)
    matchFathers = match(x = pop@father, table = parents@id)
  } else {
    if (is.null(mothers) | is.null(fathers)) {
      stop("must provide either 'parents' or both 'mothers' and 'fathers'!")
    }
    matchMothers = match(x = pop@mother, table = mothers@id)
    matchFathers = match(x = pop@father, table = fathers@id)
  }
  if (anyNA(matchMothers)) {
    stop("some parents/mothers not found!")
  }
  if (anyNA(matchFathers)) {
    stop("some parents/fathers not found!")
  }
  if (use %in% c("gv", "ebv", "pheno")) {
    if (!is.null(parents)) {
      ret = 0.5 * (slot(object = parents, name = use)[matchMothers, , drop = FALSE] +
                   slot(object = parents, name = use)[matchFathers, , drop = FALSE])
    } else {
      ret = 0.5 * (slot(object = mothers, name = use)[matchMothers, , drop = FALSE] +
                   slot(object = fathers, name = use)[matchFathers, , drop = FALSE])
    }
  } else if (use == "bv") {
    if (!is.null(parents)) {
      ret = 0.5 * (bv(parents, simParam = simParam,
                      nThreads=nThreads)[matchMothers, , drop = FALSE] +
                   bv(parents, simParam = simParam,
                      nThreads=nThreads)[matchFathers, , drop = FALSE])
    } else {
      ret = 0.5 * (bv(mothers, simParam = simParam,
                      nThreads=nThreads)[matchMothers, , drop = FALSE] +
                   bv(fathers, simParam = simParam,
                      nThreads=nThreads)[matchFathers, , drop = FALSE])
    }
  } else {
    stop("use must be one of 'gv', 'bv', 'ebv', or 'pheno'!")
  }
  return(ret)
}

#' @title Calculate Mendelian sampling
#'
#' @param pop \code{\link{Pop-class}} with individuals whose parent average
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"}, \code{"\link{bv}"},
#'   \code{"\link{ebv}"}, or \code{"\link{pheno}"}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return a matrix of Mendelian samplings with dimensions nInd by nTraits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#' mendelianSampling(pop2, parents = pop)
#' mendelianSampling(pop2, mothers = pop, fathers = pop)
#'
#' @export
mendelianSampling = function(pop, parents = NULL, mothers = NULL, fathers = NULL,
                             use = "gv", simParam = NULL, nThreads=NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  pa = parentAverage(pop = pop, parents = parents, mothers = mothers, fathers = fathers,
                     use = use, simParam = simParam, nThreads=nThreads)
  if (use %in% c("gv", "ebv", "pheno")) {
    ret = slot(object = pop, name = use) - pa
  } else if (use == "bv") {
    ret = bv(pop, simParam = simParam, nThreads=nThreads) - pa
  } else {
    stop("use must be one of 'gv', 'bv', 'ebv', or 'pheno'!")
  }
  return(ret)
}

#' @title Number of individuals
#'
#' @description A wrapper for accessing the nInd slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' nInd(pop)
#'
#' @export
nInd = function(pop){
  pop@nInd
}
