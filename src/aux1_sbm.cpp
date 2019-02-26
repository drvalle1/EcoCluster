#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' Summarize data for SBM according to the different location and species groups
//' 
//' This function determines how many observations are in each location group 
//' and species group for two cases: 
//' - when dat(i,j)=1: these results are stored in matrix nql1
//' - when dat(i,j)=0: these results are stored in matrix nql0
//' 
//' @param z integer vector with cluster assignments of locations
//' @param w integer vector with cluster assignments of species
//' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and contains the presence-absence data
//' @param ngrloc number of location groups (KL)
//' @param ngrspp number of species groups (KS)
//' @return this function returns a list with two KL x KS matrices: nql1 and nql0
//' @export
//'
// [[Rcpp::export]]
Rcpp::List getql(IntegerVector z, IntegerVector w, IntegerMatrix dat,
                 int ngrloc, int ngrspp) {
  IntegerMatrix nql1(ngrloc,ngrspp);
  IntegerMatrix nql0(ngrloc,ngrspp);
  
  for(int i=0; i<dat.nrow();i++){
    for (int j=0; j<dat.ncol(); j++){
      if (dat(i,j)==1) {
        nql1(z[i],w[j])=nql1(z[i],w[j])+1;
      }
      if (dat(i,j)==0){
        nql0(z[i],w[j])=nql0(z[i],w[j])+1;
      }
    }
  }
  
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("nql1") = nql1,
                                          Rcpp::Named("nql0") = nql0);
  return(resTemp);
}

//' Convert the vector v into the vector theta/phi using the stick-breaking formula
//' 
//' This function converts the input vector v into another vector of probabilities 
//' that sum to one based on the stick-breaking equation
//' 
//' @param v vector of probabilities (these do not have to sum to 1)
//' @return this function returns a vector res with probabilities that sum to one
//' @export
//' 
// [[Rcpp::export]]
NumericVector convertSBtoNormal(NumericVector v) {
  NumericVector res(v.size());
  
  res[0]=v[0];
  double prod=1-v[0];
  for(int j=1; j<v.size();j++){
    res(j)=v[j]*prod;    
    prod=prod*(1-v[j]);
  }
  
  return (res);
}
