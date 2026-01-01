#include "alphasimr.h"

// This is just an example - we will replace it later with more appropriate
// functions working with tree sequences.

// Note that the below call for Rcpp depends and plugins on top of export
// clashes with roxygen's export - roxygen does not export when depends and
// plugins are present:(

//' @title Get number of individuals in tree sequence
//' @description Get number of individuals in tree sequence
//' @param ts an external pointer to a \code{tsk_treeseq_t} object.
//' @return integer number of individuals
//' @examples
//' ts_file <- system.file("examples", "test.trees", package = "tskitr")
//' ts <- tskitr::ts_load(ts_file) # slendr also has ts_load()!
//' tskitr::ts_num_individuals(ts)
//' ts_num_individuals(ts)
//' ts_num_individuals2(ts)
//' AlphaSimR::ts_num_individuals2(ts)
//' AlphaSimR:::ts_num_individuals2(ts)
//'
//' @export
// [[Rcpp::depends(tskitr)]]
// [[Rcpp::plugins(tskitr)]]
// [[Rcpp::export]]
int ts_num_individuals2(SEXP ts) {
   int n;
   Rcpp::XPtr<tsk_treeseq_t> xptr(ts);
   n = (int) tsk_treeseq_get_num_individuals(xptr);
   return n;
}


