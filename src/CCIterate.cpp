#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include <sstream>
#include <stdio.h>

using Rcpp::as;
using Rcpp::clone;
using Rcpp::Function;
using Rcpp::IntegerMatrix;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::runif;

using std::max;
using std::min;

using std::stringstream;


// [[Rcpp::export]]
List CCIterate(int nIter, NumericMatrix old, NumericMatrix old_ll,
               bool single, int burn, IntegerMatrix N, 
               NumericMatrix n)
{
    List chain;
    
    for (int I = 2; I <= nIter; ++I)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                //get left hand
                double lh1, lh2, lh3, rh1, rh2, rh3;
                lh1 = (j == 0) ? 0.0 : old(i, j-1);
                lh2 = (i == 0) ? 0.0 : old(i - 1, j);
                rh1 = (j == 2) ? 1.0 : old(i,j+1);
                rh2 = (i == 2) ? 1.0 : old(i + 1, j);
                lh3 = 0;
                rh3 = 1;

                if (!single)
                {
                    bool test1 = old(1,0) < old(0,1);
                    bool test2 = old(2,1) < old(1,2);
                    if (test1 && test2)
                    {
                        if (i == 0 && j == 2)
                            lh3 = old(2,0);
                        else if (i == 2 && j == 0)
                            rh3 = old(0,2);
                    }
                    else if (!test1 && !test2) 
                    {
                        if (i == 2 && j == 0)
                            lh3 = old(0,2);
                        else if (i == 0 && j == 2)
                            rh3 = old(2,0);
                    }
                }
                double lh = max(max(lh1, lh2), lh3);

                double rh = (rh3 > lh) ? min(min(rh1, rh2), rh3) : min(rh1, rh2);
                if (rh < lh)
                    rh = 1.0;

                double draw = runif(1, lh, rh)(0);

                double ar = 2;
                
                /*************************************
                        Original like function
                **************************************
                like<-function(theta,N,n) 
                {
                    #choose(N,n)->x1 #this part cancels!
                    1->x1
                    n*log(theta)->x2
                    (N-n)*log(1-theta)->x3
                    sum(x2+x3)
                }
                */
                double new_ll = n(i,j)*log(draw) + (N(i,j) - n(i,j))*log(1.0 - draw);
                

                if (old(i,j) != 1.0 && old(i,j) != 0.0)
                    ar = exp(new_ll - old_ll(i,j));
                if (ar > runif(1)(0))
                {
                    old(i,j) = draw;
                    old_ll(i,j) = new_ll;
                }
            }
        }

        if (I > burn && I % 4 == 0)
        {
            stringstream hoop;
            hoop << I;
            chain[hoop.str()] = clone(old);
        }
    }
    return chain;
}

