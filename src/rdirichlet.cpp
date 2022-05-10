#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <sys/time.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
NumericMatrix rdirichlet_cpp(NumericMatrix alpha) {
	int n = alpha.nrow();
	int m = alpha.ncol();
	NumericMatrix rand_mat(n,m);
	for(int i = 0;i < m;i++) {
		double* results = new double[n];
		const gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		struct timeval tv; // Seed generation based on time
		gettimeofday(&tv,0);
		unsigned long mySeed = tv.tv_sec + tv.tv_usec;
		T = gsl_rng_default; // Generator setup
		r = gsl_rng_alloc (T);
		gsl_rng_set(r, mySeed);
		gsl_ran_dirichlet(r, n, alpha(_,i).begin(), results);
		gsl_rng_free(r);
		for(int j = 0;j < n;j++) {
			rand_mat(j,i) = results[j];
		}
		delete[] results;
	}
	return(rand_mat);
}

