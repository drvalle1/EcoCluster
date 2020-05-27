#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' Summarize the data for mixture model
//' 
//' This function summarizes the data by calculating the number of observations across all locations 
//' assigned to a given group for which the species was present (stored in res1) or absent (stored in res0)
//' 
//' @param dat this matrix has L rows (locations) and S columns (species)
//'            and contains the presence-absence data (i.e., number of times 
//'            a given species was observed at a given location)
//' @param n.minus.y matrix with L rows (e.g., locations) and S columns (e.g., species),
//'                  calculated as matrix(nl,L,S)-dat           
//' @param z vector of size L with cluster assignment for each location
//' @param nspp number of species (S)
//' @param nloc number of locations (L)
//' @param ngroup maximum number of groups (K)
//' @return this function returns a list containing two matrices of size K x S: res1 and res0
//' @export

//' This function calculates ncs1 and ncs0
// [[Rcpp::export]]
Rcpp::List ncs(IntegerMatrix dat, IntegerMatrix nminusy, IntegerVector z, 
               int nspp, int nloc, int ngroup) {
    
    IntegerMatrix res1(ngroup,nspp);
    IntegerMatrix res0(ngroup,nspp);
    
    for(int i=0; i<nloc; i++){
        for(int j=0; j<nspp; j++){
            res1(z(i),j)=res1(z(i),j)+dat(i,j);
            res0(z(i),j)=res0(z(i),j)+nminusy(i,j);
        }
    }
    
    Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("ncs1") = res1,
                                            Rcpp::Named("ncs0") = res0);
    return(resTemp);
}

