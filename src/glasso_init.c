// Automatically generated, editing not advised.
#ifndef R_glasso_H
#define R_glasso_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("glasso", String)
#else
#define _(String) (String)
#endif

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(glasso)(
	int *nn,
	double *sss,
	double *rrho,
	int *ia,
	int *is,
	int *itr,
	int *ipen,
	double *thr,
	int *maxit,
	double *www,
	double *wwwi,
	int *nniter,
	double *ddel,
	int *jerr
);


static R_NativePrimitiveArgType glasso_t[] = {
	INTSXP,
	REALSXP,
	REALSXP,
	INTSXP,
	INTSXP,
	INTSXP,
	INTSXP,
	REALSXP,
	INTSXP,
	REALSXP,
	REALSXP,
	INTSXP,
	REALSXP,
	INTSXP
};
void F77_SUB(glassopath)(
	double *beta,
	double *what,
	int *jerrs,
	double *rholist,
	int *nrho,
	int *n,
	double *ss,
	double *rho,
	int *ia,
	int *itr,
	int *ipen,
	double *thr,
	int *maxit,
	double *ww,
	double *wwi,
	int *niter,
	double *del,
	int *jerr
);


static R_NativePrimitiveArgType glassopath_t[] = {
	REALSXP,
	REALSXP,
	INTSXP,
	REALSXP,
	INTSXP,
	INTSXP,
	REALSXP,
	REALSXP,
	INTSXP,
	INTSXP,
	INTSXP,
	REALSXP,
	INTSXP,
	REALSXP,
	REALSXP,
	INTSXP,
	REALSXP,
	INTSXP
};
static R_FortranMethodDef fMethods[] = {
	FDEF(glasso),
	FDEF(glassopath),
	{NULL, NULL, 0}
};

void R_init_glasso(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

#endif
