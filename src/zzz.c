#include <R.h>
#include <R_ext/Rdynload.h>

#include "zzz.h"

static const R_CMethodDef cMethods[] = {
	{"bandwidth_mise_cdf", (DL_FUNC) &bandwidth_mise_cdf, 3},
	{"bandwidth_mise_gaussian", (DL_FUNC) &bandwidth_mise_gaussian, 3},
	{"gaussian_kcdf", (DL_FUNC) &gaussian_kcdf, 6},
	{"bandwidth_mise_pdf_gaussian", (DL_FUNC) &bandwidth_mise_pdf_gaussian, 3},
	{"gaussian_kpdf", (DL_FUNC) &gaussian_kpdf, 6},
	{"gamma_kcdf", (DL_FUNC) &gamma_kcdf, 6},
	{"bandwidth_mise_gamma_pdf", (DL_FUNC) &bandwidth_mise_gamma_pdf,3},
	{"gamma_kpdf_inv", (DL_FUNC) &gamma_kpdf_inv, 6},
	{"gamma_kpdf", (DL_FUNC) &gamma_kpdf, 6},
	{"RIG_kcdf", (DL_FUNC) &RIG_kcdf, 6},
	{"runAdjRand", (DL_FUNC) &runAdjRand, 8},
	{"gamma_kcdf_inv", (DL_FUNC) &gamma_kcdf_inv, 6},
	/* Finish R_CMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_SynClustR(DllInfo *info){
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_SynClustR(). */
