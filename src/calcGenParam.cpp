#include "alphasimr.h"

// Calculates genetic parameters for traits with epistasis
Rcpp::List calcGenParamE(const Rcpp::S4& trait, 
                         const Rcpp::S4& pop,
                         int nThreads){
  //Information from pop
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  //Information from trait
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a = trait.slot("addEff");
  arma::mat E;
  E = Rcpp::as<arma::mat>(trait.slot("epiEff"));
  E.col(0) -= 1; //R to C++
  E.col(1) -= 1; //R to C++
  arma::vec d;
  double intercept = trait.slot("intercept");
  arma::mat bvMat(nInd,nThreads,arma::fill::zeros); // "Breeding value"
  arma::mat aaMat(nInd,nThreads,arma::fill::zeros); // Epistatic deviations
  arma::mat gv_t; // Total genetic value
  arma::mat gv_a(nInd,nThreads,arma::fill::zeros); // Genetic value due to a
  arma::mat gv_aa(nInd,nThreads,arma::fill::zeros); // Genetic value due to aa
  arma::vec genicA(nThreads,arma::fill::zeros); // No LD
  arma::vec genicA2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicD(nThreads,arma::fill::zeros); // No LD
  arma::vec genicD2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicAA(nThreads,arma::fill::zeros); // No LD
  arma::vec genicAA2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec mu(nThreads,arma::fill::zeros); // Observed mean
  arma::vec eMu(nThreads,arma::fill::zeros); // Expected mean with HWE
  arma::vec alpha(a.n_elem);
  arma::vec alphaHW(a.n_elem);
  
  arma::mat ddMat, gv_d;
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
    ddMat.set_size(nInd,nThreads);
    ddMat.zeros();
    gv_d.set_size(nInd,nThreads);
    gv_d.zeros();
  }
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<E.n_rows; ++i){
    double gvMu1, gvMu2, gvEMu1, gvEMu2, 
    genoMu1, genoMu2, p1, p2, q1, q2, dK, 
    gvMu, gvEMu, gvNoLDMu;
    
    arma::uword tid; //Thread ID
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    
    //Observed frequencies
    arma::mat freq(ploidy+1,ploidy+1,arma::fill::zeros);
    for(arma::uword j=0; j<nInd; ++j){
      freq(genoMat(j,E(i,0)),genoMat(j,E(i,1))) += 1;
    }
    freq = freq/accu(freq);
    arma::vec freq1 = sum(freq,1);
    arma::vec freq2 = sum(freq,0).t();
    
    genoMu1 = accu(freq1%x);
    p1 = genoMu1/dP;
    q1 = 1-p1;
    
    genoMu2 = accu(freq2%x);
    p2 = genoMu2/dP;
    q2 = 1-p2;
    
    // Expected frequencies
    arma::vec freqE1(ploidy+1), freqE2(ploidy+1);
    for(arma::uword k=0; k<(ploidy+1); ++k){
      dK = double(k);
      freqE1(k) = choose(dP,dK)*std::pow(p1,dK)*std::pow(q1,dP-dK);
      freqE2(k) = choose(dP,dK)*std::pow(p2,dK)*std::pow(q2,dP-dK);
    }
    
    // Frequencies with no LD
    arma::mat freqNoLD(ploidy+1,ploidy+1);
    arma::mat freqNoLDE(ploidy+1,ploidy+1);
    for(arma::uword j=0; j<(ploidy+1); ++j){
      for(arma::uword k=0; k<(ploidy+1); ++k){
        freqNoLDE(j,k) = freqE1(j)*freqE2(k);
        freqNoLD(j,k) = freq1(j)*freq2(k);
      }
    }
    
    //Marginal values (individual loci)
    //Additive effects
    arma::vec aEff1 = xa*a(E(i,0));
    arma::vec aEff2 = xa*a(E(i,1));
    //Additive-by-additive effects
    arma::mat aaEff = xa*xa.t()*E(i,2);
    //Dominance effects
    arma::vec dEff1, dEff2;
    //Genetic value
    arma::vec gv1, gv2, gvE1, gvE2;
    if(hasD){
      dEff1 = xd*d(E(i,0));
      dEff2 = xd*d(E(i,1));
      gv1 = aEff1+dEff1;
      gv2 = aEff2+dEff2;
      gvE1 = gv1;
      gvE2 = gv2;
      for(arma::uword j=0; j<(ploidy+1); ++j){
        gv1(j) += accu(freq2%(aEff2+dEff2+E(i,2)*xa(j)*xa));
        gv2(j) += accu(freq1%(aEff1+dEff1+E(i,2)*xa(j)*xa));
        gvE1(j) += accu(freqE2%(aEff2+dEff2+E(i,2)*xa(j)*xa));
        gvE2(j) += accu(freqE1%(aEff1+dEff1+E(i,2)*xa(j)*xa));
      }
    }else{
      gv1 = aEff1;
      gv2 = aEff2;
      gvE1 = gv1;
      gvE2 = gv2;
      for(arma::uword j=0; j<(ploidy+1); ++j){
        gv1(j) += accu(freq2%(aEff2+E(i,2)*xa(j)*xa));
        gv2(j) += accu(freq1%(aEff1+E(i,2)*xa(j)*xa));
        gvE1(j) += accu(freqE2%(aEff2+E(i,2)*xa(j)*xa));
        gvE2(j) += accu(freqE1%(aEff1+E(i,2)*xa(j)*xa));
      }
    }
    
    gvMu1 = accu(freq1%gv1);
    gvMu2 = accu(freq2%gv2);
    gvEMu1 = accu(freqE1%gvE1);
    gvEMu2 = accu(freqE2%gvE2);
    
    alpha(E(i,0)) = accu(freq1%(gv1-gvMu1)%(x-genoMu1))/
      accu(freq1%(x-genoMu1)%(x-genoMu1));
    alphaHW(E(i,0)) = accu(freqE1%(gvE1-gvEMu1)%(x-genoMu1))/
      accu(freqE1%(x-genoMu1)%(x-genoMu1));
    alpha(E(i,1)) = accu(freq2%(gv2-gvMu2)%(x-genoMu2))/
      accu(freq2%(x-genoMu2)%(x-genoMu2));
    alphaHW(E(i,1)) = accu(freqE2%(gvE2-gvEMu2)%(x-genoMu2))/
      accu(freqE2%(x-genoMu2)%(x-genoMu2));
    
    //Check for divide by zero
    if(!std::isfinite(alpha(E(i,0)))) alpha(E(i,0))=0;
    if(!std::isfinite(alphaHW(E(i,0)))) alphaHW(E(i,0))=0;
    if(!std::isfinite(alpha(E(i,1)))) alpha(E(i,1))=0;
    if(!std::isfinite(alphaHW(E(i,1)))) alphaHW(E(i,0))=0;
    
    //Breeding values
    arma::vec bv1, bv2, bvE1, bvE2;
    bv1 = (x-genoMu1)*alpha(E(i,0)); //Breeding values
    bvE1 = (x-genoMu1)*alphaHW(E(i,0)); //Random mating breeding value
    bv2 = (x-genoMu2)*alpha(E(i,1)); //Breeding values
    bvE2 = (x-genoMu2)*alphaHW(E(i,0)); //Random mating breeding value
    genicA(tid) += accu(freq1%bv1%bv1);
    genicA2(tid) += accu(freqE1%bvE1%bvE1);
    genicA(tid) += accu(freq2%bv2%bv2);
    genicA2(tid) += accu(freqE2%bvE2%bvE2);
    //Dominance deviation
    arma::vec dd1, dd2, ddE1, ddE2;
    if(hasD){
      dd1 = gv1-bv1-gvMu1; //Dominance deviations (lack of fit)
      ddE1 = gvE1-bvE1-gvEMu1; //Random mating dominance deviation
      dd2 = gv2-bv2-gvMu2; //Dominance deviations (lack of fit)
      ddE2 = gvE2-bvE2-gvEMu2; //Random mating dominance deviation
      genicD(tid) += accu(freq1%dd1%dd1);
      genicD2(tid) += accu(freqE1%ddE1%ddE1);
      genicD(tid) += accu(freq2%dd2%dd2);
      genicD2(tid) += accu(freqE2%ddE2%ddE2);
    }
    
    //Joint values (both loci)
    //Genetic value matrix
    arma::mat GV(ploidy+1,ploidy+1);
    //Breeding value matrix
    arma::mat BV(ploidy+1,ploidy+1);
    arma::mat BVE(ploidy+1,ploidy+1);
    //Dominance deviation matrix
    arma::mat DD(ploidy+1,ploidy+1);
    arma::mat DDE(ploidy+1,ploidy+1);
    //Epistasis matrix (lack of fit)
    arma::mat AA(ploidy+1,ploidy+1);
    arma::mat AANoLD(ploidy+1,ploidy+1);
    arma::mat AAE(ploidy+1,ploidy+1);
    for(arma::uword j=0; j<(ploidy+1); ++j){
      for(arma::uword k=0; k<(ploidy+1); ++k){
        BV(j,k) = bv1(j)+bv2(k);
        BVE(j,k) = bvE1(j)+bvE2(k);
        if(hasD){
          GV(j,k) = xa(j)*a(E(i,0)) + xa(k)*a(E(i,1)) +
            xd(j)*d(E(i,0)) + xd(k)*d(E(i,1)) + 
            xa(j)*xa(k)*E(i,2);
          DD(j,k) = dd1(j)+dd2(k);
          DDE(j,k) = ddE1(j)+ddE2(k);
        }else{
          GV(j,k) = xa(j)*a(E(i,0)) + xa(k)*a(E(i,1)) +
            xa(j)*xa(k)*E(i,2);
        }
      }
    }
    gvMu = accu(freq%GV);
    gvNoLDMu = accu(freqNoLD%GV);
    gvEMu = accu(freqNoLDE%GV);
    mu(tid) += gvMu;
    eMu(tid) += gvEMu;
    if(hasD){
      AA = GV-BV-DD-gvMu;
      AANoLD = GV-BV-DD-gvNoLDMu;
      AAE = GV-BVE-DDE-gvEMu;
    }else{
      AA = GV-BV-gvMu;
      AANoLD = GV-BV-gvNoLDMu;
      AAE = GV-BVE-gvEMu;
    }
    genicAA(tid) += accu(freqNoLD%AANoLD%AANoLD);
    genicAA2(tid) += accu(freqNoLDE%AAE%AAE);
    
    //Fill in individual effects
    for(arma::uword j=0; j<nInd; ++j){
      bvMat(j,tid) += BV(genoMat(j,E(i,0)),genoMat(j,E(i,1)));
      aaMat(j,tid) += AA(genoMat(j,E(i,0)),genoMat(j,E(i,1)));
      gv_a(j,tid) += aEff1(genoMat(j,E(i,0)))+aEff2(genoMat(j,E(i,1)));
      gv_aa(j,tid) += aaEff(genoMat(j,E(i,0)),genoMat(j,E(i,1)));
      if(hasD){
        ddMat(j,tid) += DD(genoMat(j,E(i,0)),genoMat(j,E(i,1)));
        gv_d(j,tid) += dEff1(genoMat(j,E(i,0)))+dEff2(genoMat(j,E(i,1)));
      }
    }
    
  }
  
  if(hasD){
    gv_t = gv_a + gv_d + gv_aa;
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_t,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("dd")=sum(ddMat,1),
                              Rcpp::Named("aa")=sum(aaMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarD")=accu(genicD),
                              Rcpp::Named("genicVarAA")=accu(genicAA),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarD2")=accu(genicD2),
                              Rcpp::Named("genicVarAA2")=accu(genicAA2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_d")=sum(gv_d,1),
                              Rcpp::Named("gv_aa")=sum(gv_aa,1),
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW);
  }else{
    gv_t = gv_a + gv_aa;
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_t,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("aa")=sum(aaMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarAA")=accu(genicAA),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarAA2")=accu(genicAA2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_aa")=sum(gv_aa,1),
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW);
  }
}

// Calculates genetic parameters for traits with imprinting
Rcpp::List calcGenParamS(const Rcpp::S4& trait, 
                         const Rcpp::S4& pop,
                         int nThreads){
  //Information from pop
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  //Information from trait
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a = trait.slot("addEff");
  arma::vec d;
  arma::vec s;
  // TODO: Expand to polyploids
  arma::vec x(ploidy+2); // Genotype dossage
  x(0) = 0;
  x(1) = 1;
  x(2) = 1;
  x(3) = 2;
  arma::vec xa = (x-dP/2.0)*(2.0/dP); // -1, 0, 0, 1 for diploids
  arma::vec xaE = (x-dP/2.0)*(2.0/dP); // -1, 0, 0, 1 for diploids
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP); // 0, 1, 1, 0 for diploids
  arma::vec xdE = x%(dP-x)*(2.0/dP)*(2.0/dP); // 0, 1, 1, 0 for diploids
  // TODO expand to polyploids
  arma::vec xs = xd; // 0, -1, 1, 0 for diploids
  xs(1) = -xs(1);
  arma::vec xsE = xd; // 0, -1, 1, 0 for diploids
  xsE(1) = -xsE(1);
  double intercept = trait.slot("intercept");
  arma::mat bvMat(nInd,nThreads,arma::fill::zeros); // "Breeding value"
  arma::mat bvMatM(nInd,nThreads,arma::fill::zeros); // "Breeding value (maternal)"
  arma::mat bvMatP(nInd,nThreads,arma::fill::zeros); // "Breeding value (paternal)"
  arma::mat gv_t; // Total genetic value
  arma::mat gv_a(nInd,nThreads,arma::fill::zeros); // Genetic value due to a
  arma::vec genicA(nThreads,arma::fill::zeros); // No LD
  arma::vec genicAM(nThreads,arma::fill::zeros); // No LD (maternal)
  arma::vec genicAP(nThreads,arma::fill::zeros); // No LD (paternal)
  arma::vec genicA2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicAM2(nThreads,arma::fill::zeros); // No LD and HWE (maternal)
  arma::vec genicAP2(nThreads,arma::fill::zeros); // No LD and HWE (paternal)
  arma::vec genicD(nThreads,arma::fill::zeros); // No LD
  arma::vec genicD2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicS(nThreads,arma::fill::zeros); // No LD (genic imprinting devation variance is the same between sexes)
  arma::vec genicS2(nThreads,arma::fill::zeros); // No LD and HWE (genic imprinting devation variance is the same between sexes)
  arma::vec mu(nThreads,arma::fill::zeros); // Observed mean
  arma::vec eMu(nThreads,arma::fill::zeros); // Expected mean with HWE
  arma::mat ddMat, gv_d; // Dominance deviation and genetic value due to d
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
    ddMat.set_size(nInd,nThreads);
    ddMat.zeros();
    gv_d.set_size(nInd,nThreads);
    gv_d.zeros();
  }
  arma::mat sdMat, gv_s; // Imprinting deviation and genetic value due to s
  s = Rcpp::as<arma::vec>(trait.slot("impEff"));
  sdMat.set_size(nInd,nThreads);
  sdMat.zeros();
  gv_s.set_size(nInd,nThreads);
  gv_s.zeros();
  
  arma::vec alpha(a.n_elem);
  arma::vec alphaHW(a.n_elem);
  
  arma::vec alphaM(a.n_elem);
  arma::vec alphaMHW(a.n_elem);
  
  arma::vec alphaP(a.n_elem);
  arma::vec alphaPHW(a.n_elem);
  
  arma::vec beta(a.n_elem);
  arma::vec betaHW(a.n_elem);
  
  arma::vec gamma(a.n_elem);
  arma::vec gammaHW(a.n_elem);
  
  arma::vec m(a.n_elem);
  arma::vec mE(a.n_elem);
  
  arma::vec m_a(a.n_elem);
  arma::vec m_aE(a.n_elem);
  
  arma::vec m_d(a.n_elem);
  arma::vec m_dE(a.n_elem);
  
  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> genoMatM = getMaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> genoMatP = getPaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<a.n_elem; ++i){
    
    arma::uword tid; //Thread ID
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    
    // TODO expand to polyploids
    arma::vec freq(ploidy+2,arma::fill::zeros), freqE(ploidy+2); // Genotype frequencies, observed and HWE
    arma::vec aEff(ploidy+2), dEff(ploidy+2), sEff(ploidy+2), eff(ploidy+2); // Genetic values, additive, dominance and imprinting
    arma::vec bv(ploidy+2), dd(ploidy+2), gv(ploidy+2), gvE(ploidy+2); // Statistical values, additive and dominance
    arma::vec bvM(ploidy+2), bvP(ploidy+2); //Additive values with imprinting
    arma::vec sdM(ploidy+2), sdP(ploidy+2), sd(ploidy+2),sdE(ploidy+2); // Imprinting deviations
    arma::vec gvM(ploidy+2), gvP(ploidy+2); //Statistical values with imprinting
    arma::vec bvE(ploidy+2), ddE(ploidy+2); //Expected for random mating
    arma::vec bvME(ploidy+2), bvPE(ploidy+2),
    sdME(ploidy+2), sdPE(ploidy+2);
    double gvMu, gvEMu, genoMu, p, q, dK,index1, index2, index;
    double muA, muD, muEA, muED, muS, muES;
    
    arma::vec xa_i = xa; //I use different local variables, not shared
    arma::vec xd_i = xd;
    arma::vec xs_i = xs;
    arma::vec xaE_i = xaE; 
    arma::vec xdE_i = xdE;
    arma::vec xsE_i = xsE;
    
    // Compute genotype frequencies
    for(arma::uword j=0; j<nInd; ++j){
      index1 = genoMatM(j,i);
      index2 = genoMatP(j,i);
      index2 = index2*2;
      index = index1+index2;
      freq(index) += 1;
    }
    freq = freq/accu(freq);
    genoMu = accu(freq%x);
    p = genoMu/dP;
    q = 1-p;
    
    // Expected genotype frequencies
    // TODO Generalize this method
    freqE.zeros();
    freqE(0) = q*q;
    freqE(1) = q*p;
    freqE(2) = p*q;
    freqE(3) = p*p;
    
    // Set genetic values
    aEff = xa_i*a(i);
    sEff = xs_i*s(i);
    if(hasD){
      dEff = xd_i*d(i);
      gv = aEff+dEff+sEff; // -a, d-i, d+i, a for diploids
    }else{
      gv = aEff+sEff; // -a, -i, +i, a for diploids
    }
    
    // Mean genetic values
    gvMu = accu(freq%gv);
    gvEMu =  accu(freqE%gv);
    gv = gv-gvMu;
    gvE = gv-gvEMu;
    
    mu(tid) += gvMu;
    eMu(tid) += gvEMu;
    
    muA = accu(freq%xa_i); // Do I need more muA? We are substracting from xa and it is not more useful
    muEA = accu(freqE%xa_i);
    
    xa_i = xa_i - muA;
    xaE_i = xaE_i - muEA;
    
    // Average effect
    alpha(i) = accu(freq%gv%xa_i)/
      accu(freq%xa_i%xa_i);
    alphaHW(i) = accu(freqE%gvE%xaE_i)/
      accu(freq%xaE_i%xaE_i);
    
    // Check for divide by zero
    if(!std::isfinite(alpha(i))) alpha(i)=0;
    if(!std::isfinite(alphaHW(i))) alphaHW(i)=0;
    
    // Set additive genic variances
    bv = xa_i*alpha(i); //Breeding values
    bvE = xaE_i*alphaHW(i); //Random mating breeding value
    genicA(tid) += accu(freq%bv%bv);
    genicA2(tid) += accu(freq%bvE%bvE);
    
    // Set dominance genic variances
    if(hasD){
      muD = accu(freq%xd_i); // sum(freq*x_d) for centering dominance vector
      muED = accu(freqE%xd_i); // sum(freq*x_d) for centering dominance vector
      
      xd_i = xd_i - muD;
      xdE_i = xdE_i - muED;
      
      m(i) = accu(freq%xd_i%xa)/ // Regression coefficient from regressing x_d on x_a
        accu(freq%xa_i%xa_i);
      mE(i) = accu(freqE%xdE_i%xaE_i)/ // Regression coefficient from regressing x_d on x_a
        accu(freqE%xaE_i%xaE_i);
      
      if(!std::isfinite(m(i))) m(i)=0;
      if(!std::isfinite(mE(i))) mE(i)=0;
      
      xd_i = xd_i - xa_i*m(i); // centering by x_a*m
      xdE_i = xdE_i - xaE_i*mE(i); // centering by x_aE*mE
      
      muS = accu(freq%xs_i); // sum(freq*x_i) for centering imprinting vector
      muES = accu(freqE%xs_i); // sum(freq*x_i) for centering imprinting vector
      
      xs_i = xs_i - muS;
      xsE_i = xs_i - muES;
      
      m_a(i) = accu(freq%xs_i%xa_i)/ // Regression coefficient from regressing x_i on x_a
        accu(freq%xa_i%xa_i);
      m_aE(i) = accu(freqE%xsE_i%xaE_i)/ // Regression coefficient from regressing x_i on x_a
        accu(freqE%xaE_i%xaE_i);

      if(!std::isfinite(m_a(i))) m_a(i)=0;
      if(!std::isfinite(m_aE(i))) m_aE(i)=0;
      
      m_d = accu(freq%xs_i%xd_i)/ // Regression coefficient from regressing x_i on x_d
        accu(freq%xd_i%xd_i);
      m_dE = accu(freqE%xsE_i%xdE_i)/ // Regression coefficient from regressing x_i on x_d
        accu(freqE%xdE_i%xdE_i);
      
      if(!std::isfinite(m_d(i))) m_d(i)=0;
      if(!std::isfinite(m_dE(i))) m_dE(i)=0;
      
      xs_i = xs_i - xa_i*m_a - xd_i*m_d(i);
      xsE_i = xsE_i - xaE_i*m_aE - xdE_i*m_dE(i);
      
      beta(i) = accu(freq%gv%xd_i) / // Calculate beta
        accu(freq%xd_i%xd_i);
      betaHW(i) = accu(freq%gvE%xdE_i) / // Calculate betaHW
        accu(freqE%xdE_i%xdE_i);
      
      if(!std::isfinite(beta(i))) beta(i)=0;
      if(!std::isfinite(betaHW(i))) betaHW(i)=0;
      
      dd = xd_i*beta(i); // Dominance deviation values
      ddE = xdE_i*betaHW(i); // Random mating dominance deviation values
      genicD(tid) += accu(freq%dd%dd);
      genicD2(tid) += accu(freqE%ddE%ddE);
      
      gamma(i) = accu(freq%gv%xs_i) / // Calculate gammaE
        accu(freq%xs_i%xs_i);
      gammaHW(i) = accu(freqE%gvE%xsE_i) / // Calculate gammaE
        accu(freqE%xsE_i%xsE_i);
      
      if(!std::isfinite(gamma(i))) gamma(i)=0;
      if(!std::isfinite(gammaHW(i))) gammaHW(i)=0;
      
      sd = xs_i*gamma(i); // Silencing (imprinting) deviation values
      sdE = xsE_i*gammaHW(i); // Random mating silencing (imprinting) deviation values
      
    } else {
    
    sd = gv-bv; // Silencing (imprinting) deviations (lack of fit)
    sdE = gvE-bvE; // Random mating silencng (imprinting) deviation  
    
    }
    
    genicS(tid) += accu(freq%sd%sd);
    genicS2(tid) += accu(freqE%sdE%sdE);
    
    alphaM(i)  = alpha(i)  - s(i);
    alphaP(i)  = alpha(i)  + s(i);
    alphaMHW(i) = alphaHW(i) - s(i);
    alphaPHW(i) = alphaHW(i) + s(i);
    
    bvM = xa_i*alphaM(i); // Breeding values (maternal)
    bvP = xa_i*alphaP(i); // Breeding values (paternal)
    bvME = xaE_i*alphaMHW(i); // Random mating breeding value (maternal)
    bvPE = xaE_i*alphaPHW(i); // Random mating breeding value (paternal)
    genicAM(tid) += accu(freq%bvM%bvM);
    genicAP(tid) += accu(freq%bvP%bvP);
    genicAM2(tid) += accu(freqE%bvME%bvME);
    genicAP2(tid) += accu(freqE%bvPE%bvPE);
    
    // Set values for individuals
    for(arma::uword j=0; j<nInd; ++j){
      index1 = genoMatM(j,i);
      index2 = genoMatP(j,i);
      index2 = index2*2;
      // TODO expand to polyploids!
      index = index1+index2;
      
      gv_a(j,tid) += aEff(index);
      bvMat(j,tid) += bv(index);
      bvMatM(j,tid) += bvM(index);
      bvMatP(j,tid) += bvP(index);
      
      if(hasD){
        gv_d(j,tid) += dEff(index);
        ddMat(j,tid) += dd(index);
      }
      
      gv_s(j,tid) += sEff(index);
      sdMat(j,tid) += sd(index);
    }
  }
  if(hasD){
    gv_t = gv_a + gv_d + gv_s;
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_t,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("bvM")=sum(bvMatM,1),
                              Rcpp::Named("bvP")=sum(bvMatP,1),
                              Rcpp::Named("dd")=sum(ddMat,1),
                              Rcpp::Named("id")=sum(sdMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarD")=accu(genicD),
                              Rcpp::Named("genicVarI")=accu(genicS),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarD2")=accu(genicD2),
                              Rcpp::Named("genicVarI2")=accu(genicS2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_d")=sum(gv_d,1),
                              Rcpp::Named("gv_i")=sum(gv_s,1),
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW,
                              Rcpp::Named("alphaM")=alphaM,
                              Rcpp::Named("alpha_MHW")=alphaMHW,
                              Rcpp::Named("alphaP")=alphaP,
                              Rcpp::Named("alpha_PHW")=alphaPHW);
    
  }else{
    gv_t = gv_a + gv_s;
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_t,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("bvM")=sum(bvMatM,1),
                              Rcpp::Named("bvP")=sum(bvMatP,1),
                              Rcpp::Named("id")=sum(sdMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarI")=accu(genicS),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarI2")=accu(genicS2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_i")=sum(gv_s,1),
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW,
                              Rcpp::Named("alphaM")=alphaM,
                              Rcpp::Named("alpha_MHW")=alphaMHW,
                              Rcpp::Named("alphaP")=alphaP,
                              Rcpp::Named("alpha_PHW")=alphaPHW);
  }
}


// Calculates breeding values, dominance deviations and genic
// variances. Additive and dominance genetic variances are calculated
// from breeding values and dominance deviations. 
// [[Rcpp::export]]
Rcpp::List calcGenParam(const Rcpp::S4& trait, 
                        const Rcpp::S4& pop,
                        int nThreads){
  if(trait.hasSlot("epiEff")){
    return calcGenParamE(trait, pop, nThreads);
  }
  if(trait.hasSlot("impEff")){
    return calcGenParamS(trait, pop, nThreads);
  }
  //Information from pop
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  //Information from trait
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a = trait.slot("addEff");
  arma::vec d;
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP); // -1, 0, 1 for diploids
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP); // 0, 1, 0 for diploids
  double intercept = trait.slot("intercept");
  arma::mat bvMat(nInd,nThreads,arma::fill::zeros); // "Breeding value"
  arma::mat gv_t; // Total genetic value
  arma::mat gv_a(nInd,nThreads,arma::fill::zeros); // Genetic value due to a
  arma::vec genicA(nThreads,arma::fill::zeros); // No LD
  arma::vec genicA2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicD(nThreads,arma::fill::zeros); // No LD
  arma::vec genicD2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec mu(nThreads,arma::fill::zeros); // Observed mean
  arma::vec eMu(nThreads,arma::fill::zeros); // Expected mean with HWE
  arma::mat ddMat, gv_d; // Dominance deviation and genetic value due to d
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
    ddMat.set_size(nInd,nThreads);
    ddMat.zeros();
    gv_d.set_size(nInd,nThreads);
    gv_d.zeros();
  }
  arma::vec alpha(a.n_elem);
  arma::vec alphaHW(a.n_elem);
  
  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<a.n_elem; ++i){
    
    arma::uword tid; //Thread ID
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    
    arma::vec freq(ploidy+1,arma::fill::zeros), freqE(ploidy+1); // Genotype frequencies, observed and HWE
    arma::vec aEff(ploidy+1), dEff(ploidy+1), eff(ploidy+1); // Genetic values, additive and dominance
    arma::vec bv(ploidy+1), dd(ploidy+1), gv(ploidy+1); // Statistical values, additive and dominance
    arma::vec bvE(ploidy+1), ddE(ploidy+1); //Expected for random mating
    double gvMu, gvEMu, genoMu, p, q, dK;
    
    // Compute genotype frequencies
    for(arma::uword j=0; j<nInd; ++j){
      freq(genoMat(j,i)) += 1;
    }
    freq = freq/accu(freq);
    genoMu = accu(freq%x);
    p = genoMu/dP;
    q = 1-p;
    
    // Expected genotype frequencies
    freqE.zeros();
    for(arma::uword k=0; k<(ploidy+1); ++k){
      dK = double(k);
      freqE(k) = choose(dP,dK)*std::pow(p,dK)*std::pow(q,dP-dK);
    }
    
    // Set genetic values
    aEff = xa*a(i);
    if(hasD){
      dEff = xd*d(i);
      gv = aEff+dEff; // -a, d, a for diploids
    }else{
      gv = aEff; // -a, 0, a for diploids
    }
    
    // Mean genetic values
    gvMu = accu(freq%gv);
    gvEMu =  accu(freqE%gv);
    mu(tid) += gvMu;
    eMu(tid) += gvEMu;
    
    // Average effect
    alpha(i) = accu(freq%(gv-gvMu)%(x-genoMu))/
      accu(freq%(x-genoMu)%(x-genoMu));
    alphaHW(i) = accu(freqE%(gv-gvEMu)%(x-genoMu))/
      accu(freqE%(x-genoMu)%(x-genoMu)); 
    
    // Check for division by zero
    if(!std::isfinite(alpha(i))) alpha(i)=0;
    if(!std::isfinite(alphaHW(i))) alphaHW(i)=0;
    
    // Set additive genic variances
    bv = (x-genoMu)*alpha(i); //Breeding values
    bvE = (x-genoMu)*alphaHW(i); //Random mating breeding value
    genicA(tid) += accu(freq%bv%bv);
    genicA2(tid) += accu(freqE%bvE%bvE);
    
    // Set dominance genic variances
    if(hasD){
      dd = gv-bv-gvMu; //Dominance deviations (lack of fit)
      ddE = gv-bvE-gvEMu; //Random mating dominance deviation
      genicD(tid) += accu(freq%dd%dd);
      genicD2(tid) += accu(freqE%ddE%ddE);
    }
    
    // Set values for individuals
    for(arma::uword j=0; j<nInd; ++j){
      gv_a(j,tid) += aEff(genoMat(j,i));
      bvMat(j,tid) += bv(genoMat(j,i));
      if(hasD){
        gv_d(j,tid) += dEff(genoMat(j,i));
        ddMat(j,tid) += dd(genoMat(j,i));
      }
    }
  }
  if(hasD){
    gv_t = gv_a + gv_d;
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_t,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("dd")=sum(ddMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarD")=accu(genicD),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarD2")=accu(genicD2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_d")=sum(gv_d,1),
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW);
  }else{
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_a,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW);
  }
}
