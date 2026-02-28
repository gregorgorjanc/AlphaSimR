#include "alphasimr.h"
#include <RcppTskit.hpp>
// [[Rcpp::depends(RcppTskit)]]
// [[Rcpp::plugins(RcppTskit)]]

// TODO: This is just an example - we will replace it later with more
// appropriate
//       functions working with tree sequences. For example to obtain haplotypes
//       using https://tskit.dev/tskit/docs/stable/c-api.html#decoding-genotypes
//
//       See also ts: Create haplotypes from a tree sequence for downstream
//       AlphaSimR work (but subsample sites!) #4
//       https://github.com/HighlanderLab/AlphaSimR/issues/4
//' @title Summarise `tskit` table collection
//' @param tc an external pointer to a \code{tsk_table_collection_t} object.
//' @return A list.
//' @examples
//' ts_file <- system.file("examples", "test.trees", package = "RcppTskit")
//' tc <- RcppTskit:::tc_ptr_load(ts_file)
//' RcppTskit:::tc_ptr_summary(tc)
//' AlphaSimR:::tc_ptr_summary2(tc) # a bit simpler at this stage ...
//' tc <- RcppTskit::TableCollection$new(ts_file)
//' RcppTskit:::tc_ptr_summary(tc$pointer)
//' AlphaSimR:::tc_ptr_summary2(tc$pointer) # a bit simpler at this stage ...
//' @export
// [[Rcpp::export]]
Rcpp::List tc_ptr_summary2(const SEXP tc) {
  RcppTskit_table_collection_xptr tc_xptr(tc);
  const tsk_table_collection_t *tables = tc_xptr;
  return Rcpp::List::create(
      Rcpp::_["num_provenances"] = tables->provenances.num_rows,
      Rcpp::_["num_populations"] = tables->populations.num_rows,
      Rcpp::_["num_migrations"] = tables->migrations.num_rows,
      Rcpp::_["num_individuals"] = tables->individuals.num_rows,
      Rcpp::_["num_nodes"] = tables->nodes.num_rows,
      Rcpp::_["num_edges"] = tables->edges.num_rows,
      Rcpp::_["num_sites"] = tables->sites.num_rows,
      Rcpp::_["num_mutations"] = tables->mutations.num_rows,
      Rcpp::_["sequence_length"] = tables->sequence_length);
}
