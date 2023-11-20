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
    alpha1, alpha2, alphaE1, alphaE2,
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

    alpha1 = accu(freq1%(gv1-gvMu1)%(x-genoMu1))/
      accu(freq1%(x-genoMu1)%(x-genoMu1));
    alphaE1 = accu(freqE1%(gvE1-gvEMu1)%(x-genoMu1))/
      accu(freqE1%(x-genoMu1)%(x-genoMu1));
    alpha2 = accu(freq2%(gv2-gvMu2)%(x-genoMu2))/
      accu(freq2%(x-genoMu2)%(x-genoMu2));
    alphaE2 = accu(freqE2%(gvE2-gvEMu2)%(x-genoMu2))/
      accu(freqE2%(x-genoMu2)%(x-genoMu2));

    //Check for divide by zero
    if(!std::isfinite(alpha1)) alpha1=0;
    if(!std::isfinite(alphaE1)) alphaE1=0;
    if(!std::isfinite(alpha2)) alpha2=0;
    if(!std::isfinite(alphaE2)) alphaE2=0;

    //Breeding values
    arma::vec bv1, bv2, bvE1, bvE2;
    bv1 = (x-genoMu1)*alpha1; //Breeding values
    bvE1 = (x-genoMu1)*alphaE1; //Random mating breeding value
    bv2 = (x-genoMu2)*alpha2; //Breeding values
    bvE2 = (x-genoMu2)*alphaE2; //Random mating breeding value
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
                              Rcpp::Named("gv_mu")=intercept);
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
                              Rcpp::Named("gv_mu")=intercept);
  }
}

// Calculates breeding values, dominance and imprinting deviations and genic
// variances. Additive, dominance, and imprinting genetic variances are calculated
// from breeding values, dominance, and imprinting deviations.
// [[Rcpp::export]]
Rcpp::List calcGenParam(const Rcpp::S4& trait,
                        const Rcpp::S4& pop,
                        int nThreads){
  if(trait.hasSlot("epiEff")){
    return calcGenParamE(trait, pop, nThreads);
  }

  //Information from pop
  bool hasD = trait.hasSlot("domEff");
  bool hasS = trait.hasSlot("impEff"); // Imprinting (=silencing) but using s since i is an iterator
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  //Information from trait
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a = trait.slot("addEff");
  arma::vec d;
  arma::vec s;
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP); // -1, 0, 1 for diploids
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP); // 0, 1, 0 for diploids
  // TODO expand to polyploids
  arma::vec xsM = xd; // 0,  1, 0 for diploids
  arma::vec xsP = xd; // 0, -1, 0 for diploids
  xsP(1) = -xsP(1);
  double intercept = trait.slot("intercept");
  arma::mat bvMat(nInd,nThreads,arma::fill::zeros); // "Breeding value"
  arma::mat bvMMat(nInd,nThreads,arma::fill::zeros); // "Breeding value" (maternal)
  arma::mat bvPMat(nInd,nThreads,arma::fill::zeros); // "Breeding value" (paternal)
  arma::mat gv_t; // Total genetic value
  arma::mat gv_a(nInd,nThreads,arma::fill::zeros); // Genetic value due to a
  arma::vec genicA(nThreads,arma::fill::zeros); // No LD
  arma::vec genicA2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicAM(nThreads,arma::fill::zeros); // No LD (maternal)
  arma::vec genicAP(nThreads,arma::fill::zeros); // No LD (paternal)
  arma::vec genicAM2(nThreads,arma::fill::zeros); // No LD and HWE (maternal)
  arma::vec genicAP2(nThreads,arma::fill::zeros); // No LD and HWE (paternal)
  arma::vec genicD(nThreads,arma::fill::zeros); // No LD
  arma::vec genicD2(nThreads,arma::fill::zeros); // No LD and HWE
  arma::vec genicS(nThreads,arma::fill::zeros); // No LD TODO: genic imprinting devation variance is the same between sexes
  arma::vec genicS2(nThreads,arma::fill::zeros); // No LD and HWE TODO: genic imprinting devation variance is the same between sexes
  arma::vec genicSM(nThreads,arma::fill::zeros); // No LD (maternal) TODO: genic imprinting devation variance is the same between sexes
  arma::vec genicSP(nThreads,arma::fill::zeros); // No LD (paternal) TODO: genic imprinting devation variance is the same between sexes
  arma::vec genicSM2(nThreads,arma::fill::zeros); // No LD and HWE (maternal) TODO: genic imprinting devation variance is the same between sexes
  arma::vec genicSP2(nThreads,arma::fill::zeros); // No LD and HWE (paternal) TODO: genic imprinting devation variance is the same between sexes
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
  arma::mat sdMMat, sdPMat, gv_s, gv_sM, gv_sP; // Imprinting deviation and genetic value due to s
  if(hasS){
    s = Rcpp::as<arma::vec>(trait.slot("impEff"));
    sdMMat.set_size(nInd,nThreads);
    sdPMat.set_size(nInd,nThreads);
    sdMMat.zeros();
    sdPMat.zeros();
    gv_s.set_size(nInd,nThreads);
    gv_sM.set_size(nInd,nThreads);
    gv_sP.set_size(nInd,nThreads);
    gv_s.zeros();
    gv_sM.zeros();
    gv_sP.zeros();
  }

  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")),
                                             lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> genoMatM;
  if(hasS){
    genoMatM = getMaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")),
                               lociPerChr, lociLoc, nThreads);
  }

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

    arma::vec tmp(ploidy+1,arma::fill::zeros), freq(ploidy+1,arma::fill::zeros), freqE(ploidy+1); // Genotype frequencies, observed and HWE
    arma::vec freqHetM(ploidy+1,arma::fill::zeros), freqHetME(ploidy+1); // Frequencies of phased heterozygotes, observed and HWE
    arma::vec aEff(ploidy+1), dEff(ploidy+1); // Genetic values, additive and dominance (imprinting below)
    arma::vec bv(ploidy+1), dd(ploidy+1), gv(ploidy+1); // Statistical values, additive and dominance (imprinting below)
    arma::vec bvE(ploidy+1), ddE(ploidy+1); // Expected for random mating
    double gvMu, gvEMu, genoMu, p, q, dK, alpha, alphaE;
    arma::vec sEffM(ploidy+1), sEffP(ploidy+1);
    arma::vec bvM(ploidy+1), bvP(ploidy+1),
              sdM(ploidy+1), sdP(ploidy+1),
              gvM(ploidy+1), gvP(ploidy+1);
    arma::vec bvME(ploidy+1), bvPE(ploidy+1),
              sdME(ploidy+1), sdPE(ploidy+1);
    double alphaM, alphaP, alphaME, alphaPE;

    // Compute genotype frequencies
    for(arma::uword j=0; j<nInd; ++j){
      freq(genoMat(j,i)) += 1;
    }
    tmp = freq; // when tmp = accu(tmp) there werea problem due size is not = 1. I replace it to be the same object pre-processed
    freq = freq/accu(tmp);
    genoMu = accu(freq%x);
    p = genoMu/dP;
    q = 1-p;

    // Compute frequencies of phased heterozygotes
    if(hasS){
      for(arma::uword j=0; j<nInd; ++j){
        if(0 < genoMat(j,i) && genoMat(j,i) < ploidy){
          freqHetM(genoMatM(j,i)) += 1;
          // freqHetM(0) is freq of 0-1 (mat-pat) heterozygote
          // freqHetM(1) is freq of 1-0 (mat-pat) heterozygote
        }
      }
      freqHetM = freqHetM/accu(tmp); // there was also a problem due the freqHetM size. I did it equal to freq
    }

    // Expected genotype frequencies
    freqE.zeros();
    for(arma::uword k=0; k<(ploidy+1); ++k){
      dK = double(k);
      freqE(k) = choose(dP,dK)*std::pow(p,dK)*std::pow(q,dP-dK);
    }

    // Expected frequencies of maternal heterozygotes
    freqHetME.zeros();
    for(arma::uword k=1; k<(ploidy); ++k){
      dK = double(k);
      freqHetME(k) = std::pow(p,dK)*std::pow(q,dP-dK);
    }

    // Set genetic values
    aEff = xa*a(i);
    gv = aEff; // -a, 0, a for diploids
    if(hasD){
      dEff = xd*d(i);
      gv = gv+dEff; // -a, d, a for diploids
    }
    if(hasS){
      sEffM = xsM*s(i);
      sEffP = xsP*s(i);
      gvM = gv+sEffM; // -a, d+i, a for diploids
      gvP = gv+sEffP; // -a, d-i, a for diploids
    }

    // Mean genetic values
    gvMu = accu(freq%gv);
    gvEMu =  accu(freqE%gv);
    if(hasS){
      // the above is gvMu = freq(0)*gv(0) + freq(1)*gv(1) + freq(2)*gv(2) for diploids
      // TODO: change freq(1)-freqHetM(1) to freqHetM(0), but is freqHetME(0) correct?
      gvMu = freq(0)*gv(0) + (freq(1)-freqHetM(1))*gvP(1) + freqHetM(1)*gvM(1) + freq(2)*gv(2);
      gvEMu = freqE(0)*gv(0) + (freqE(1)-freqHetME(1))*gvP(1) + freqHetME(1)*gvM(1) + freqE(2)*gv(2);
    }
    mu(tid) += gvMu;
    eMu(tid) += gvEMu;

    // Average effect
    alpha = accu(freq%(gv-gvMu)%(x-genoMu))/
      accu(freq%(x-genoMu)%(x-genoMu));
    alphaE = accu(freqE%(gv-gvEMu)%(x-genoMu))/
      accu(freqE%(x-genoMu)%(x-genoMu));
    if(hasS){
      // TODO: we need to check these calculations
      alphaM = accu(freq%(gvM-gvMu)%(x-genoMu))/
        accu(freq%(x-genoMu)%(x-genoMu));
      alphaP = accu(freq%(gvP-gvMu)%(x-genoMu))/
        accu(freq%(x-genoMu)%(x-genoMu));
      alphaME = accu(freqE%(gvM-gvEMu)%(x-genoMu))/
        accu(freqE%(x-genoMu)%(x-genoMu));
      alphaPE = accu(freqE%(gvP-gvEMu)%(x-genoMu))/
        accu(freqE%(x-genoMu)%(x-genoMu));
    }

    // Check for divide by zero
    if(!std::isfinite(alpha)) alpha=0;
    if(!std::isfinite(alphaE)) alphaE=0;
    if(hasS){
      if(!std::isfinite(alphaM)) alphaM=0;
      if(!std::isfinite(alphaP)) alphaP=0;
      if(!std::isfinite(alphaME)) alphaME=0;
      if(!std::isfinite(alphaPE)) alphaPE=0;
    }

    // Set additive genic variances
    bv = (x-genoMu)*alpha; //Breeding values
    bvE = (x-genoMu)*alphaE; //Random mating breeding value
    genicA(tid) += accu(freq%bv%bv);
    genicA2(tid) += accu(freqE%bvE%bvE);
    if(hasS){
      bvM = (x-genoMu)*alphaM; //Breeding values (maternal)
      bvP = (x-genoMu)*alphaP; //Breeding values (paternal)
      bvME = (x-genoMu)*alphaME; //Random mating breeding value (maternal)
      bvPE = (x-genoMu)*alphaPE; //Random mating breeding value (paternal)
      genicAM(tid) += accu(freq%bvM%bvM);
      genicAP(tid) += accu(freq%bvP%bvP);
      genicAM2(tid) += accu(freqE%bvME%bvME);
      genicAP2(tid) += accu(freqE%bvPE%bvPE);
    }

    // Set dominance genic variances
    if(hasD){
      dd = gv-bv-gvMu; //Dominance deviations (lack of fit)
      ddE = gv-bvE-gvEMu; //Random mating dominance deviation
      genicD(tid) += accu(freq%dd%dd);
      genicD2(tid) += accu(freqE%ddE%ddE);
    }

    // Set imprinting genic variances
    if(hasS){
      // TODO: extend for dominance! Check also elsewhere
      sdM = gvM-bvM-gvMu; //Imprinting deviations (lack of fit) (maternal)
      sdP = gvP-bvP-gvMu; //Imprinting deviations (lack of fit) (paternal)
      sdME = gvM-bvME-gvEMu; //Random mating imprinting deviation (maternal)
      sdPE = gvP-bvPE-gvEMu; //Random mating imprinting deviation (paternal)
      //genicS(tid) += accu(freq%sd%sd);
      genicSM(tid) += accu(freq%sdM%sdM);
      genicSP(tid) += accu(freq%sdP%sdP);
      //genicS2(tid) += accu(freqE%sdE%sdE);
      genicSM2(tid) += accu(freqE%sdME%sdME);
      genicSP2(tid) += accu(freqE%sdPE%sdPE);
    }

    // Set values for individuals
    for(arma::uword j=0; j<nInd; ++j){
      gv_a(j,tid) += aEff(genoMat(j,i));
      bvMat(j,tid) += bv(genoMat(j,i));
      if(hasS){
        bvMMat(j,tid) += bvM(genoMat(j,i));
        bvPMat(j,tid) += bvP(genoMat(j,i));
      }
      if(hasD){
        gv_d(j,tid) += dEff(genoMat(j,i));
        ddMat(j,tid) += dd(genoMat(j,i));
      }
      if(hasS){
        gv_sM(j,tid) += sEffM(genoMat(j,i));
        gv_sP(j,tid) += sEffP(genoMat(j,i));
        gv_s(j,tid) += gv_sM(j,tid) * (1 - genoMatM(j,i)) +
                       gv_sP(j,tid) *      genoMatM(j,i);
        // no need for sdMat since it would be zero (by definition)
        sdMMat(j,tid) += sdM(genoMat(j,i));
        sdPMat(j,tid) += sdP(genoMat(j,i));
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
                              Rcpp::Named("gv_mu")=intercept);
  }else if(hasS){
    gv_t = gv_a + gv_s;  //gv_d causes problem due it is dim(0x0), I also romoved intercept.
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_t,1), //1  //I had a problem with the substitute function
                              Rcpp::Named("bv")=sum(bvMat,1),     // then, I have used the old one
                              Rcpp::Named("bvM")=sum(bvMMat,1),   // I have not included the dominance effects and intercept
                              Rcpp::Named("bvP")=sum(bvPMat,1),
                              Rcpp::Named("idM")=sum(sdMMat,1),
                              Rcpp::Named("idP")=sum(sdPMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarAM")=accu(genicAM),
                              Rcpp::Named("genicVarAP")=accu(genicAP),
                              Rcpp::Named("genicVarIM")=accu(genicSM),  //10
                              Rcpp::Named("genicVarIP")=accu(genicSP),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarAM2")=accu(genicAM2),
                              Rcpp::Named("genicVarAP2")=accu(genicAP2),
                              Rcpp::Named("genicVarIM2")=accu(genicSM2),
                              Rcpp::Named("genicVarIP2")=accu(genicSP2),
                              Rcpp::Named("mu")=accu(mu),
                              Rcpp::Named("mu_HWE")=accu(eMu),
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_i")=sum(gv_s,1)); //20
                    
    
    //Rcpp::Named("genicVarI")=accu(genicS),
    //Rcpp::Named("genicVarI2")=accu(genicS2),
    //Rcpp::Named("gv_iM")=sum(gv_sM,1),
    //Rcpp::Named("gv_iP")=sum(gv_sP,1),
  }else{
    return Rcpp::List::create(Rcpp::Named("gv")=sum(gv_a,1)+intercept,
                              Rcpp::Named("bv")=sum(bvMat,1),
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("gv_a")=sum(gv_a,1),
                              Rcpp::Named("gv_mu")=intercept);
  }
}
