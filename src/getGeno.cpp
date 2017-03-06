// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/*
 * Retrieves allele dossages from requested loci
 * geno: raw genotype data
 * nInd: number of individuals in population
 * nChr: number of chromosomes
 * ploidy: ploidy of species
 * lociPerChr: number of loci per chromosome
 * lociLoc: physical position of loci
 */
// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(Rcpp::List& geno, int nInd, int nChr, int ploidy,
                                 arma::ivec& lociPerChr, arma::uvec& lociLoc){
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr), arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  for(int i=0; i<nChr; ++i){
    loc1 = loc2+1;
    loc2 += lociPerChr[i];
    Rcpp::List chrGeno = geno[i];
    arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2))-1; //R to C++
    for(int j=0; j<ploidy; ++j)
      output.cols(loc1,loc2) +=
        Rcpp::as<arma::Mat<unsigned char> >(chrGeno[j]).cols(chrLociLoc);
  }
  return output;
}