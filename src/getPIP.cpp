#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
int find_closest(const int seed, const int scan, const IntegerVector scans, const NumericVector mzs, const IntegerVector clu, const double tol){
  int u = -1, d = -1;
  bool ub = true, db = true;
  for (int i=0; i<scans.size(); ++i){
    if (ub && scans[seed-i] == scan && clu[seed-i] == 0) {
      u = seed-i;
      ub = false;
    }
    if (db && scans[seed+i] == scan && clu[seed+i] == 0) {
      d = seed+i;
      db = false;
    }
    if (mzs[seed-i] < mzs[seed]-tol) {
      ub = false;
    }
    if (mzs[seed+i] > mzs[seed]+tol) {
      db = false;
    }
    if (!ub && !db){
      break;
    }
  }
  
  if (u > 0 && d > 0){
    if (mzs[seed]-mzs[u] < mzs[d]-mzs[seed]) {
      return u;
    } else {
      return d;
    }
  } else {
    return std::max(u, d);
  }
}

// [[Rcpp::export]]
List get_trace(const int seed, const IntegerVector scans, const NumericVector mzs, IntegerVector clu, double tol){
  List output;
  IntegerVector res;
  res.push_back(seed);
  int l,r;
  int lg = 0, rg = 0;
  bool lb = true, rb = true; 
  for (int i=1; i<scans.size(); i++){
    if (lb) {
      l = find_closest(seed, scans[seed]-i, scans, mzs, clu, tol);
      if (l >= 0) {
        res.push_front(l);
        clu[l] = 1;
        lg = 0;
      } else {
        lg++;
      }
    }
    if (rb) {
      r = find_closest(seed, scans[seed]+i, scans, mzs, clu, tol);
      if (r >= 0) {
        res.push_back(r);
        clu[r] = 1;
        rg = 0;
      } else {
        rg ++;
      }
    }
    if (lg > 2) {lb = false;}
    if (rg > 2) {rb = false;}
    if (!lb && !rb) {break;}
  }
  
  output["res"] = res;
  output["clu"] = clu;
  return output;
}

// [[Rcpp::export]]
IntegerVector getROI(int seed, IntegerVector scans, NumericVector mzs, NumericVector ints, double mztol, int max_width) {
  int ref = seed-1;
  int end = mzs.size();
  double refMz = mzs[ref];
  int refScan = scans[ref];
  
  // locate roi
  int leftScan = refScan - 0.5 * max_width;
  int rightScan = refScan + 0.5 * max_width;
  double upMz = refMz - mztol;
  double downMz = refMz + mztol;
  
  IntegerVector roi = Rcpp::IntegerVector::create(ref);
  for (int j=ref-1; j>0; j--){
    if (mzs[j]<upMz){
      break;
    }
    if (scans[j]>=leftScan && scans[j]<=rightScan){
      roi.push_back(j);
    }
  }
  for (int k=ref+1; k<end; k++){
    if (mzs[k]>downMz){
      break;
    }
    if (scans[k]>=leftScan && scans[k]<=rightScan){
      roi.push_back(k);
    }
  }
  
  return roi;
}