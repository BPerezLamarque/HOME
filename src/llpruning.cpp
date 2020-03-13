#include <Rcpp.h>
                   
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vmult(NumericVector m1, NumericMatrix m2){
  NumericVector out(m1.size());
  NumericVector cm2;
  for (int j = 0; j < m2.ncol(); ++j) {
    cm2 = m2(_,j);
    out(j) = std::inner_product(m1.begin(), m1.end(), cm2.begin(), 0.);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector vvmult(NumericVector m1, NumericVector m2){
  NumericVector out(4);
  out(0)=m1(0)*m2(0);
  out(1)=m1(1)*m2(1);
  out(2)=m1(2)*m2(2);
  out(3)=m1(3)*m2(3);
  return out;
}

// [[Rcpp::export]]
double vvvmult(NumericVector m1, NumericVector m2){
  double out=0;
  out=m1(0)*m2(0)+m1(1)*m2(1)+m1(2)*m2(2)+m1(3)*m2(3);
  return out;
}

// [[Rcpp::export]]
double llpruning(NumericVector nodes, NumericVector edges, NumericVector el, NumericMatrix L, NumericVector eig_val, NumericMatrix eig_vect, NumericMatrix ivp, NumericVector propinv) {
  
  int nnode = nodes.size();
  int index1 = 0;
  int index2 = 0;
  int v1 = 0;
  int v2 = 0;
  double t1 = 0;
  double t2 = 0;
  for(int i = 0; i < nnode; ++i){
    index1 = 2*i;
    index2 = 2*i+1;
    v1 = edges(index1);
    v2 = edges(index2);
    t1 = el(index1);
    t2 = el(index2);
    L(nodes(i),Rcpp::_)=vvmult(vmult(vvmult(vmult(L(v1,Rcpp::_), eig_vect),exp(t1*eig_val)),ivp), vmult(vvmult(vmult(L(v2,Rcpp::_), eig_vect),exp(t2*eig_val)),ivp));
  }
  t1=log(vvvmult(propinv,L(nodes(nnode-1),Rcpp::_)));
  return t1;
}
