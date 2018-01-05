// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// writeASGenotypes
void writeASGenotypes(const arma::Cube<unsigned char>& g, const arma::field<arma::uvec>& locations, const arma::uvec& allLocations, const arma::vec& snpchips, const std::vector<std::string>& names, const char missing, const std::string fname);
RcppExport SEXP _AlphaSimR_writeASGenotypes(SEXP gSEXP, SEXP locationsSEXP, SEXP allLocationsSEXP, SEXP snpchipsSEXP, SEXP namesSEXP, SEXP missingSEXP, SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Cube<unsigned char>& >::type g(gSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type allLocations(allLocationsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type snpchips(snpchipsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type names(namesSEXP);
    Rcpp::traits::input_parameter< const char >::type missing(missingSEXP);
    Rcpp::traits::input_parameter< const std::string >::type fname(fnameSEXP);
    writeASGenotypes(g, locations, allLocations, snpchips, names, missing, fname);
    return R_NilValue;
END_RCPP
}
// writeASHaplotypes
void writeASHaplotypes(const arma::Cube<unsigned char>& g, const arma::field<arma::uvec>& locations, const arma::uvec& allLocations, const arma::vec& snpchips, const std::vector<std::string>& names, const char missing, const std::string fname);
RcppExport SEXP _AlphaSimR_writeASHaplotypes(SEXP gSEXP, SEXP locationsSEXP, SEXP allLocationsSEXP, SEXP snpchipsSEXP, SEXP namesSEXP, SEXP missingSEXP, SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Cube<unsigned char>& >::type g(gSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type allLocations(allLocationsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type snpchips(snpchipsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type names(namesSEXP);
    Rcpp::traits::input_parameter< const char >::type missing(missingSEXP);
    Rcpp::traits::input_parameter< const std::string >::type fname(fnameSEXP);
    writeASHaplotypes(g, locations, allLocations, snpchips, names, missing, fname);
    return R_NilValue;
END_RCPP
}
// gebvRR
arma::mat gebvRR(const Rcpp::S4& RRsol, const Rcpp::S4& pop);
RcppExport SEXP _AlphaSimR_gebvRR(SEXP RRsolSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type RRsol(RRsolSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(gebvRR(RRsol, pop));
    return rcpp_result_gen;
END_RCPP
}
// gegvRRD
arma::mat gegvRRD(const Rcpp::S4& RRsol, const Rcpp::S4& pop);
RcppExport SEXP _AlphaSimR_gegvRRD(SEXP RRsolSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type RRsol(RRsolSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(gegvRRD(RRsol, pop));
    return rcpp_result_gen;
END_RCPP
}
// gebvRRD
arma::mat gebvRRD(const Rcpp::S4& RRsol, const Rcpp::S4& pop);
RcppExport SEXP _AlphaSimR_gebvRRD(SEXP RRsolSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type RRsol(RRsolSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(gebvRRD(RRsol, pop));
    return rcpp_result_gen;
END_RCPP
}
// gebvGCA
arma::mat gebvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, bool female, bool isSCAsol);
RcppExport SEXP _AlphaSimR_gebvGCA(SEXP solSEXP, SEXP popSEXP, SEXP femaleSEXP, SEXP isSCAsolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type sol(solSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< bool >::type female(femaleSEXP);
    Rcpp::traits::input_parameter< bool >::type isSCAsol(isSCAsolSEXP);
    rcpp_result_gen = Rcpp::wrap(gebvGCA(sol, pop, female, isSCAsol));
    return rcpp_result_gen;
END_RCPP
}
// gebvSCA
arma::mat gebvSCA(const Rcpp::S4& sol, const Rcpp::S4& pop, bool isSCAsol);
RcppExport SEXP _AlphaSimR_gebvSCA(SEXP solSEXP, SEXP popSEXP, SEXP isSCAsolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type sol(solSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< bool >::type isSCAsol(isSCAsolSEXP);
    rcpp_result_gen = Rcpp::wrap(gebvSCA(sol, pop, isSCAsol));
    return rcpp_result_gen;
END_RCPP
}
// getGeno
arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc);
RcppExport SEXP _AlphaSimR_getGeno(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    rcpp_result_gen = Rcpp::wrap(getGeno(geno, lociPerChr, lociLoc));
    return rcpp_result_gen;
END_RCPP
}
// getDomGeno
arma::imat getDomGeno(const arma::Mat<unsigned char>& geno);
RcppExport SEXP _AlphaSimR_getDomGeno(SEXP genoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<unsigned char>& >::type geno(genoSEXP);
    rcpp_result_gen = Rcpp::wrap(getDomGeno(geno));
    return rcpp_result_gen;
END_RCPP
}
// getHaplo
arma::Mat<unsigned char> getHaplo(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc);
RcppExport SEXP _AlphaSimR_getHaplo(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    rcpp_result_gen = Rcpp::wrap(getHaplo(geno, lociPerChr, lociLoc));
    return rcpp_result_gen;
END_RCPP
}
// getOneHaplo
arma::Mat<unsigned char> getOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc, int haplo);
RcppExport SEXP _AlphaSimR_getOneHaplo(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP, SEXP haploSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    Rcpp::traits::input_parameter< int >::type haplo(haploSEXP);
    rcpp_result_gen = Rcpp::wrap(getOneHaplo(geno, lociPerChr, lociLoc, haplo));
    return rcpp_result_gen;
END_RCPP
}
// writeGeno
void writeGeno(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc, Rcpp::String filePath);
RcppExport SEXP _AlphaSimR_writeGeno(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP, SEXP filePathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type filePath(filePathSEXP);
    writeGeno(geno, lociPerChr, lociLoc, filePath);
    return R_NilValue;
END_RCPP
}
// writeOneHaplo
void writeOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc, int haplo, Rcpp::String filePath);
RcppExport SEXP _AlphaSimR_writeOneHaplo(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP, SEXP haploSEXP, SEXP filePathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    Rcpp::traits::input_parameter< int >::type haplo(haploSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type filePath(filePathSEXP);
    writeOneHaplo(geno, lociPerChr, lociLoc, haplo, filePath);
    return R_NilValue;
END_RCPP
}
// getGv
arma::field<arma::vec> getGv(const Rcpp::S4& trait, const Rcpp::S4& pop);
RcppExport SEXP _AlphaSimR_getGv(SEXP traitSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(getGv(trait, pop));
    return rcpp_result_gen;
END_RCPP
}
// calcGenParam
Rcpp::List calcGenParam(const Rcpp::S4& trait, const Rcpp::S4& pop);
RcppExport SEXP _AlphaSimR_calcGenParam(SEXP traitSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(calcGenParam(trait, pop));
    return rcpp_result_gen;
END_RCPP
}
// getHybridGv
arma::field<arma::vec> getHybridGv(const Rcpp::S4& trait, const Rcpp::S4& motherGeno, arma::uvec& mother, const Rcpp::S4& fatherGeno, arma::uvec& father);
RcppExport SEXP _AlphaSimR_getHybridGv(SEXP traitSEXP, SEXP motherGenoSEXP, SEXP motherSEXP, SEXP fatherGenoSEXP, SEXP fatherSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type motherGeno(motherGenoSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type fatherGeno(fatherGenoSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type father(fatherSEXP);
    rcpp_result_gen = Rcpp::wrap(getHybridGv(trait, motherGeno, mother, fatherGeno, father));
    return rcpp_result_gen;
END_RCPP
}
// cross2
Rcpp::List cross2(const arma::field<arma::Cube<unsigned char> >& motherGeno, arma::uvec mother, const arma::field<arma::Cube<unsigned char> >& fatherGeno, arma::uvec father, const arma::field<arma::vec>& genMaps, double recombRatio, bool trackRec);
RcppExport SEXP _AlphaSimR_cross2(SEXP motherGenoSEXP, SEXP motherSEXP, SEXP fatherGenoSEXP, SEXP fatherSEXP, SEXP genMapsSEXP, SEXP recombRatioSEXP, SEXP trackRecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type motherGeno(motherGenoSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type fatherGeno(fatherGenoSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type father(fatherSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::vec>& >::type genMaps(genMapsSEXP);
    Rcpp::traits::input_parameter< double >::type recombRatio(recombRatioSEXP);
    Rcpp::traits::input_parameter< bool >::type trackRec(trackRecSEXP);
    rcpp_result_gen = Rcpp::wrap(cross2(motherGeno, mother, fatherGeno, father, genMaps, recombRatio, trackRec));
    return rcpp_result_gen;
END_RCPP
}
// createDH2
Rcpp::List createDH2(const arma::field<arma::Cube<unsigned char> >& geno, int nDH, const arma::field<arma::vec>& genMaps, double recombRatio, bool useFemale, bool trackRec);
RcppExport SEXP _AlphaSimR_createDH2(SEXP genoSEXP, SEXP nDHSEXP, SEXP genMapsSEXP, SEXP recombRatioSEXP, SEXP useFemaleSEXP, SEXP trackRecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< int >::type nDH(nDHSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::vec>& >::type genMaps(genMapsSEXP);
    Rcpp::traits::input_parameter< double >::type recombRatio(recombRatioSEXP);
    Rcpp::traits::input_parameter< bool >::type useFemale(useFemaleSEXP);
    Rcpp::traits::input_parameter< bool >::type trackRec(trackRecSEXP);
    rcpp_result_gen = Rcpp::wrap(createDH2(geno, nDH, genMaps, recombRatio, useFemale, trackRec));
    return rcpp_result_gen;
END_RCPP
}
// popVar
arma::mat popVar(const arma::mat& X);
RcppExport SEXP _AlphaSimR_popVar(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(popVar(X));
    return rcpp_result_gen;
END_RCPP
}
// mergeGeno
arma::field<arma::Cube<unsigned char> > mergeGeno(const arma::field<arma::Cube<unsigned char> >& x, const arma::field<arma::Cube<unsigned char> >& y);
RcppExport SEXP _AlphaSimR_mergeGeno(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(mergeGeno(x, y));
    return rcpp_result_gen;
END_RCPP
}
// calcChrFreq
arma::vec calcChrFreq(const arma::Cube<unsigned char>& geno);
RcppExport SEXP _AlphaSimR_calcChrFreq(SEXP genoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Cube<unsigned char>& >::type geno(genoSEXP);
    rcpp_result_gen = Rcpp::wrap(calcChrFreq(geno));
    return rcpp_result_gen;
END_RCPP
}
// convToImat
arma::imat convToImat(const arma::Mat<unsigned char>& X);
RcppExport SEXP _AlphaSimR_convToImat(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<unsigned char>& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(convToImat(X));
    return rcpp_result_gen;
END_RCPP
}
// sampAllComb
arma::Mat<arma::uword> sampAllComb(arma::uword nLevel1, arma::uword nLevel2, arma::uword n);
RcppExport SEXP _AlphaSimR_sampAllComb(SEXP nLevel1SEXP, SEXP nLevel2SEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type nLevel1(nLevel1SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nLevel2(nLevel2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sampAllComb(nLevel1, nLevel2, n));
    return rcpp_result_gen;
END_RCPP
}
// sampHalfDialComb
arma::Mat<arma::uword> sampHalfDialComb(arma::uword nLevel, arma::uword n);
RcppExport SEXP _AlphaSimR_sampHalfDialComb(SEXP nLevelSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type nLevel(nLevelSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sampHalfDialComb(nLevel, n));
    return rcpp_result_gen;
END_RCPP
}
// calcCoef
arma::mat calcCoef(arma::mat& X, arma::mat& Y);
RcppExport SEXP _AlphaSimR_calcCoef(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(calcCoef(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// readMat
arma::mat readMat(std::string fileName, int rows, int cols, char sep, int skipRows, int skipCols);
RcppExport SEXP _AlphaSimR_readMat(SEXP fileNameSEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP sepSEXP, SEXP skipRowsSEXP, SEXP skipColsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< char >::type sep(sepSEXP);
    Rcpp::traits::input_parameter< int >::type skipRows(skipRowsSEXP);
    Rcpp::traits::input_parameter< int >::type skipCols(skipColsSEXP);
    rcpp_result_gen = Rcpp::wrap(readMat(fileName, rows, cols, sep, skipRows, skipCols));
    return rcpp_result_gen;
END_RCPP
}
// solveUVM
Rcpp::List solveUVM(const arma::mat& y, const arma::mat& X, const arma::mat& Z, const arma::mat& K);
RcppExport SEXP _AlphaSimR_solveUVM(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(solveUVM(y, X, Z, K));
    return rcpp_result_gen;
END_RCPP
}
// solveAniModel
Rcpp::List solveAniModel(const arma::mat& y, const arma::mat& X, const arma::mat& K);
RcppExport SEXP _AlphaSimR_solveAniModel(SEXP ySEXP, SEXP XSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(solveAniModel(y, X, K));
    return rcpp_result_gen;
END_RCPP
}
// solveRRBLUP
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X, const arma::mat& M);
RcppExport SEXP _AlphaSimR_solveRRBLUP(SEXP ySEXP, SEXP XSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(solveRRBLUP(y, X, M));
    return rcpp_result_gen;
END_RCPP
}
// solveMVM
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X, const arma::mat& Z, const arma::mat& K, double tol, int maxIter);
RcppExport SEXP _AlphaSimR_solveMVM(SEXP YSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP KSEXP, SEXP tolSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(solveMVM(Y, X, Z, K, tol, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// solveRRBLUPMV
Rcpp::List solveRRBLUPMV(const arma::mat& Y, const arma::mat& X, const arma::mat& M, double tol, int maxIter);
RcppExport SEXP _AlphaSimR_solveRRBLUPMV(SEXP YSEXP, SEXP XSEXP, SEXP MSEXP, SEXP tolSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(solveRRBLUPMV(Y, X, M, tol, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// solveMKM
Rcpp::List solveMKM(arma::mat& y, arma::mat& X, arma::field<arma::mat>& Zlist, arma::field<arma::mat>& Klist, int maxIter);
RcppExport SEXP _AlphaSimR_solveMKM(SEXP ySEXP, SEXP XSEXP, SEXP ZlistSEXP, SEXP KlistSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type Zlist(ZlistSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type Klist(KlistSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(solveMKM(y, X, Zlist, Klist, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// solveRRBLUPMK
Rcpp::List solveRRBLUPMK(arma::mat& y, arma::mat& X, arma::field<arma::mat>& Mlist, int maxIter);
RcppExport SEXP _AlphaSimR_solveRRBLUPMK(SEXP ySEXP, SEXP XSEXP, SEXP MlistSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type Mlist(MlistSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(solveRRBLUPMK(y, X, Mlist, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x, arma::vec reps, std::string genoTrain, int nMarker, int skip);
RcppExport SEXP _AlphaSimR_callRRBLUP(SEXP ySEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoTrainSEXP, SEXP nMarkerSEXP, SEXP skipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoTrain(genoTrainSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP(y, x, reps, genoTrain, nMarker, skip));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP_D
Rcpp::List callRRBLUP_D(arma::mat y, arma::uvec x, arma::vec reps, std::string genoTrain, int nMarker, int skip);
RcppExport SEXP _AlphaSimR_callRRBLUP_D(SEXP ySEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoTrainSEXP, SEXP nMarkerSEXP, SEXP skipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoTrain(genoTrainSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP_D(y, x, reps, genoTrain, nMarker, skip));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP_MV
Rcpp::List callRRBLUP_MV(arma::mat Y, arma::uvec x, arma::vec reps, std::string genoTrain, int nMarker, int skip, int maxIter);
RcppExport SEXP _AlphaSimR_callRRBLUP_MV(SEXP YSEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoTrainSEXP, SEXP nMarkerSEXP, SEXP skipSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoTrain(genoTrainSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP_MV(Y, x, reps, genoTrain, nMarker, skip, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP_GCA
Rcpp::List callRRBLUP_GCA(arma::mat y, arma::uvec x, arma::vec reps, std::string genoFemale, std::string genoMale, int nMarker, int skip, int maxIter);
RcppExport SEXP _AlphaSimR_callRRBLUP_GCA(SEXP ySEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoFemaleSEXP, SEXP genoMaleSEXP, SEXP nMarkerSEXP, SEXP skipSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoFemale(genoFemaleSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoMale(genoMaleSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP_GCA(y, x, reps, genoFemale, genoMale, nMarker, skip, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP_SCA
Rcpp::List callRRBLUP_SCA(arma::mat y, arma::uvec x, arma::vec reps, std::string genoFemale, std::string genoMale, int nMarker, int skip, int maxIter);
RcppExport SEXP _AlphaSimR_callRRBLUP_SCA(SEXP ySEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoFemaleSEXP, SEXP genoMaleSEXP, SEXP nMarkerSEXP, SEXP skipSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoFemale(genoFemaleSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoMale(genoMaleSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP_SCA(y, x, reps, genoFemale, genoMale, nMarker, skip, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// calcG
arma::mat calcG(arma::mat X);
RcppExport SEXP _AlphaSimR_calcG(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(calcG(X));
    return rcpp_result_gen;
END_RCPP
}
// calcD
arma::mat calcD(arma::mat X);
RcppExport SEXP _AlphaSimR_calcD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(calcD(X));
    return rcpp_result_gen;
END_RCPP
}
// calcGIbs
arma::mat calcGIbs(arma::mat X);
RcppExport SEXP _AlphaSimR_calcGIbs(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(calcGIbs(X));
    return rcpp_result_gen;
END_RCPP
}
// fastDist
arma::mat fastDist(const arma::mat& X);
RcppExport SEXP _AlphaSimR_fastDist(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastDist(X));
    return rcpp_result_gen;
END_RCPP
}
// fastPairDist
arma::mat fastPairDist(const arma::mat& X, const arma::mat& Y);
RcppExport SEXP _AlphaSimR_fastPairDist(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(fastPairDist(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// gaussKernel
arma::mat gaussKernel(arma::mat& D, double theta);
RcppExport SEXP _AlphaSimR_gaussKernel(SEXP DSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussKernel(D, theta));
    return rcpp_result_gen;
END_RCPP
}
// packHaplo
arma::Cube<unsigned char> packHaplo(arma::Mat<unsigned char>& haplo, arma::uword ploidy, bool inbred);
RcppExport SEXP _AlphaSimR_packHaplo(SEXP haploSEXP, SEXP ploidySEXP, SEXP inbredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned char>& >::type haplo(haploSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< bool >::type inbred(inbredSEXP);
    rcpp_result_gen = Rcpp::wrap(packHaplo(haplo, ploidy, inbred));
    return rcpp_result_gen;
END_RCPP
}
// MaCS
Rcpp::List MaCS(Rcpp::String args, arma::uword maxSites);
RcppExport SEXP _AlphaSimR_MaCS(SEXP argsSEXP, SEXP maxSitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type args(argsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxSites(maxSitesSEXP);
    rcpp_result_gen = Rcpp::wrap(MaCS(args, maxSites));
    return rcpp_result_gen;
END_RCPP
}
// tuneTraitA
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno, arma::vec& addEff, double varG);
RcppExport SEXP _AlphaSimR_tuneTraitA(SEXP genoSEXP, SEXP addEffSEXP, SEXP varGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned char>& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type addEff(addEffSEXP);
    Rcpp::traits::input_parameter< double >::type varG(varGSEXP);
    rcpp_result_gen = Rcpp::wrap(tuneTraitA(geno, addEff, varG));
    return rcpp_result_gen;
END_RCPP
}
// tuneTraitAD
Rcpp::List tuneTraitAD(arma::Mat<unsigned char>& geno, arma::vec& addEff, arma::vec& domEff, double varG, bool useVarA);
RcppExport SEXP _AlphaSimR_tuneTraitAD(SEXP genoSEXP, SEXP addEffSEXP, SEXP domEffSEXP, SEXP varGSEXP, SEXP useVarASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned char>& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type addEff(addEffSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type domEff(domEffSEXP);
    Rcpp::traits::input_parameter< double >::type varG(varGSEXP);
    Rcpp::traits::input_parameter< bool >::type useVarA(useVarASEXP);
    rcpp_result_gen = Rcpp::wrap(tuneTraitAD(geno, addEff, domEff, varG, useVarA));
    return rcpp_result_gen;
END_RCPP
}
