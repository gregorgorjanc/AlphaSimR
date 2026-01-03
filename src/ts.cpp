#include "alphasimr.h"
// [[Rcpp::depends(tskitr)]]
// [[Rcpp::plugins(tskitr)]]

// This is just an example - we will replace it later with more appropriate
// functions working with tree sequences.

//' @title Get number of individuals in tree sequence
//' @description Get number of individuals in tree sequence
//' @param ts an external pointer to a \code{tsk_treeseq_t} object.
//' @return integer number of individuals
//' @examples
//' ts_file <- system.file("examples", "test.trees", package = "tskitr")
//' ts <- tskitr::ts_load(ts_file) # slendr also has ts_load()!
//' tskitr::ts_num_individuals(ts)
//' AlphaSimR::ts_num_individuals2(ts)
//' @export
// [[Rcpp::export]]
int ts_num_individuals2(SEXP ts) {
   int n;
   Rcpp::XPtr<tsk_treeseq_t> xptr(ts);
   n = (int) tsk_treeseq_get_num_individuals(xptr);
   return n;
}


