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
//' This function summarizes the data by calculating the number of locations for each 
//' location group and species for which dat(i,j)=1 (stored in res1) and dat(i,j)=0 (stored in res0)
//' 
//' @param dat matrix with L rows (e.g., locations) and S columns (e.g., species),
//'        containing the presence-absence data
//' @param z vector of size L with cluster assignment for each location
//' @param nspp number of species (S)
//' @param nloc number of locations (L)
//' @param ngroup maximum number of groups (K)
//' @return this function returns a list containing two matrices of size K x S: res1 and res0
//' @export
//'
// [[Rcpp::export]]
Rcpp::List ncs(IntegerMatrix dat, IntegerVector z,
               int nspp, int nloc, int ngroup) {
    
    IntegerMatrix res1(ngroup,nspp);
    IntegerMatrix res0(ngroup,nspp);
    
    for(int i=0; i<nloc; i++){
        for(int j=0; j<nspp; j++){
            if(dat(i,j)==1) res1(z(i),j)=res1(z(i),j)+1;
            if(dat(i,j)==0) res0(z(i),j)=res0(z(i),j)+1;
        }
    }
    
    Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("ncs1") = res1,
                                            Rcpp::Named("ncs0") = res0);
    return(resTemp);
}
