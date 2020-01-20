#include <R.h>
#include <Rinternals.h>

void bandwidth_mise_cdf(double (*x), int (*n), double (*bw));

void bandwidth_mise_gaussian(double (*x), int (*n),double (*bw));

void gaussian_kcdf(double (*x), int (*n), double (*xgrid), int (*m),double (*bw), double (*Fhat));

void bandwidth_mise_pdf_gaussian(double (*x), int (*n), double (*bw));

void gaussian_kpdf(double *x, int (*n), double (*xgrid), int (*m), double (*bw), double *fhat);

void gamma_kcdf(double (*x), int (*n), double (*xgrid), int (*m),double (*bw), double (*Fhat));

void bandwidth_mise_gamma_pdf(double *x, int (*n), double (*bw));

void gamma_kpdf_inv(double (*x), int (*n), double (*xgrid), int (*m),double (*bw), double (*fhat));

void gamma_kpdf(double (*x), int (*n), double (*xgrid), int (*m),double (*bw), double (*fhat));

void RIG_kcdf(double (*x), int (*n), double(*xgrid),int (*m),double (*bw), double (*Fhat));

void runAdjRand(int (*n), int (*K1), int (*K2), int *id1, int *id2,
	double (*Rand),
	double (*aRand), double (*F));

void gamma_kcdf_inv(double (*x), int (*n),  double(*xgrid),int (*m),double (*bw), double (*Fhat));
