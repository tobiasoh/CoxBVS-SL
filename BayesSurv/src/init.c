#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BayesSurv_calJpost_helper_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _BayesSurv_matProdVec(void *, void *);
extern SEXP _BayesSurv_settingInterval_cpp(void *, void *, void *, void *);
extern SEXP _BayesSurv_sumMatProdVec(void *, void *);
extern SEXP _BayesSurv_updateBH_cpp(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _BayesSurv_updateBH_list_cpp(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _BayesSurv_updateRP_genomic_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_BayesSurv_calJpost_helper_cpp",  (DL_FUNC) &_BayesSurv_calJpost_helper_cpp,   9},
    {"_BayesSurv_matProdVec",           (DL_FUNC) &_BayesSurv_matProdVec,            2},
    {"_BayesSurv_settingInterval_cpp",  (DL_FUNC) &_BayesSurv_settingInterval_cpp,   4},
    {"_BayesSurv_sumMatProdVec",        (DL_FUNC) &_BayesSurv_sumMatProdVec,         2},
    {"_BayesSurv_updateBH_cpp",         (DL_FUNC) &_BayesSurv_updateBH_cpp,          7},
    {"_BayesSurv_updateBH_list_cpp",    (DL_FUNC) &_BayesSurv_updateBH_list_cpp,     7},
    {"_BayesSurv_updateRP_genomic_cpp", (DL_FUNC) &_BayesSurv_updateRP_genomic_cpp, 11},
    {NULL, NULL, 0}
};

void R_init_BayesSurv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
