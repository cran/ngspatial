#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ngspatial_bmse(SEXP);
extern SEXP ngspatial_randWalk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ngspatial_randWalkTrain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ngspatial_rautologistic_(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ngspatial_bmse",           (DL_FUNC) &ngspatial_bmse,            1},
    {"ngspatial_randWalk",       (DL_FUNC) &ngspatial_randWalk,       11},
    {"ngspatial_randWalkTrain",  (DL_FUNC) &ngspatial_randWalkTrain,   9},
    {"ngspatial_rautologistic_", (DL_FUNC) &ngspatial_rautologistic_,  3},
    {NULL, NULL, 0}
};

void R_init_ngspatial(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
