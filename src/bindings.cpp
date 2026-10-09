// R bindings of the numeric core in src/core (shared with the Python package lisaclust).
#include <Rcpp.h>

#include <cmath>
#include <vector>

#include "lisaclust/core.hpp"

namespace {
std::vector<lisaclust::Ring> toRings(Rcpp::List rings) {
  std::vector<lisaclust::Ring> win(rings.size());
  for (R_xlen_t i = 0; i < rings.size(); ++i) {
    Rcpp::List ring = rings[i];
    Rcpp::NumericVector rx = ring["x"], ry = ring["y"];
    for (R_xlen_t j = 0; j < rx.size(); ++j) win[i].push_back({rx[j], ry[j]});
  }
  return win;
}
}  // namespace

// [[Rcpp::export(.discWindowArea)]]
Rcpp::NumericVector discWindowArea(Rcpp::NumericVector x, Rcpp::NumericVector y, double r, int npoly,
                                   Rcpp::List rings) {
  Rcpp::NumericVector out(x.size());
  lisaclust::discAreas(x.begin(), y.begin(), x.size(), r, npoly, toRings(rings), out.begin());
  return out;
}

// [[Rcpp::export(.localCurves)]]
Rcpp::List localCurvesR(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::IntegerVector type, int nTypes,
                        Rcpp::NumericVector Rs, Rcpp::NumericVector labelVal, Rcpp::NumericVector wt,
                        Rcpp::NumericVector lam, Rcpp::NumericMatrix edge, bool Lfunction,
                        bool includeSelf = true) {
  const int n = x.size();
  std::vector<int> t0(n);
  for (int i = 0; i < n; ++i) t0[i] = type[i] - 1;
  lisaclust::LocalCurves res = lisaclust::localCurves(
      x.begin(), y.begin(), t0.data(), n, nTypes, std::vector<double>(Rs.begin(), Rs.end()),
      std::vector<double>(labelVal.begin(), labelVal.end()), wt.begin(),
      std::vector<double>(lam.begin(), lam.end()), edge.begin(), Lfunction, includeSelf);
  Rcpp::NumericVector value(res.value.begin(), res.value.end());
  for (double& v : value) if (std::isnan(v)) v = NA_REAL;
  value.attr("dim") = Rcpp::IntegerVector::create(res.n, res.K, res.nb);
  return Rcpp::List::create(
      Rcpp::Named("value") = value,
      Rcpp::Named("cell") = Rcpp::LogicalVector(res.cellPresent.begin(), res.cellPresent.end()),
      Rcpp::Named("bin") = Rcpp::LogicalVector(res.binPresent.begin(), res.binPresent.end()),
      Rcpp::Named("type") = Rcpp::LogicalVector(res.typePresent.begin(), res.typePresent.end()));
}

// [[Rcpp::export(.nearestLabels)]]
Rcpp::IntegerVector nearestLabelsR(Rcpp::NumericVector tx, Rcpp::NumericVector ty, Rcpp::IntegerVector label,
                                   Rcpp::NumericVector qx, Rcpp::NumericVector qy) {
  std::vector<int> l0(label.size());
  for (R_xlen_t i = 0; i < label.size(); ++i) l0[i] = label[i] - 1;
  std::vector<int> res = lisaclust::nearestLabels(tx.begin(), ty.begin(), l0.data(), tx.size(), qx.begin(),
                                                  qy.begin(), qx.size());
  Rcpp::IntegerVector out(res.size());
  for (std::size_t i = 0; i < res.size(); ++i) out[i] = res[i] < 0 ? NA_INTEGER : res[i] + 1;
  return out;
}
