

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.

JPS 2021
This file was created by using this command within top level of ibdreg:
which is necessary for any package with c code:
R> tools::package_native_routine_registration_skeleton(".")

All code below (and above) comments was the output of that result

Then add the .registration=TRUE string as second argument to useDynLib() in NAMESPACE

*/

/* .C calls */
extern void lsConstrain(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sim_mark_prop_free_mem();
extern void sim_mark_prop_gen_ped(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void simulate_marker_propagation(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"lsConstrain",                 (DL_FUNC) &lsConstrain,                 17},
    {"sim_mark_prop_free_mem",      (DL_FUNC) &sim_mark_prop_free_mem,       0},
    {"sim_mark_prop_gen_ped",       (DL_FUNC) &sim_mark_prop_gen_ped,       11},
    {"simulate_marker_propagation", (DL_FUNC) &simulate_marker_propagation,  4},
    {NULL, NULL, 0}
};

void R_init_ibdreg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
