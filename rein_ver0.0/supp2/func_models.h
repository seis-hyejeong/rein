#ifndef _func_models_h_
#define _func_models_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "parameters.h"

/* initialize */
void initialize_model(MODEL *);

/* replace model */
/* void replace_model(MODEL m1, MODEL m2)
*/
void replace_model(MODEL, MODEL *);

/* prototype of the functions */
/* 
void first_model(PARAM param, PRIOR prior, MODEL m)
*/
void first_model(PARAM, PRIOR, MODEL *);

/*
int model_accept_tf(gsl_rng *r, LIKELI likeli0, LIKELI likeli1)
*/
int model_accept_tf(gsl_rng *, LIKELI, LIKELI);

/* 
int PT_model_accept_tf(gsl_rng *r, LIKELI likeli0, LIKELI likeli1, double invtemp, PARATEMP parallelt)
*/ 
int PT_model_accept_tf(gsl_rng *, LIKELI, LIKELI, double, PARATEMP);

/*
int newmodel_ocean(gsl_rng * r, PRIOR prior, int nunknown, int nlayers_3, double model0[], double *model1, int iflag)
*/
int newmodel_ocean(gsl_rng *, PRIOR, int, int, double [], double *, int);

/* 
int newmodel_ocean_all(gsl_rng * r, PRIOR prior, int nlayers_3, double model0[], double *model1, int iflag)
*/
int newmodel_ocean_all(gsl_rng *, PRIOR, int, double [], double *, int);

/* supplementary of model change -- for the density 
int newmodel_ocean_supp(int index, PARAM param, PRIOR priorparam, double *model1)
*/
int newmodel_ocean_supp(int, PARAM, PRIOR, double *);

/*
void modelvec2inputvecs_ocean(int nlayer, double model[], float *mh, float *mvp, float *mvs, float *mrho) 
*/
void modelvec2inputvecs_ocean(int, double [], float *, float *, float *, float *);

/*
void modelvec2inputvecs_general(PARAM param, MODEL *m)
*/
void modelvec2inputvecs_general(PARAM, MODEL *);

/*
int model_prior_probability(int indx, double model[], double lowlim[], double highlim[]) 
*/
int model_prior_probability_ind(int, double [], double [], double []);

/*
int model_prior_probability_all(PRIOR prior, double model[])
*/
int model_prior_probability_all(PRIOR, double []);

/*
int model_prior_monotonic_increase(int nlayer, float vec[]) 
*/
int model_prior_monotonic_increase(int, float []);

/* check whether the Vp is always larger than Vs 
int model_vp_and_vs(int nlayer, float vp[], float vs[])
*/
int model_vp_and_vs(int, float [], float []);

/* check prior probability monotonic increase and Vp & Vs simultaneously
*/
int model_check_prior(MODEL, PARAM);

/* making new model and checking the prior probability all at once 
 int model_getnew_check_prior(gsl_rng *r, PRIOR prior, PARAM param, MODEL m0, MODEL *m1)
*/
int model_getnew_check_prior(gsl_rng *, PRIOR, PARAM, MODEL, MODEL *);

/*
float vp_empirical_rho_solid(float vp) 
*/
float vp_empirical_rho_solid(float);

/*
void get_densityrange_from_vp(double lowvp, double highvp, double *lowrho, double *highrho) 
*/
void get_densityrange_from_vp(double, double, double *, double *);

/* 
int get_minmax_range(PARAM param, PRIOR prior, OUTPARAM o_param)
*/
int get_minmax_range(PARAM, PRIOR, OUTPARAM *);

/*
int grid_arrray_get_length(double intv, double min, double max) 
*/
int grid_arrray_get_length(double, double, double);

/*
void grid_array(double intv, int ngrid, double min, double *array) 
*/
void grid_array(double, int, double, double *);

/*
void models_2_addcount(int nlayer, float mh[], float mvar[], int nh, double hgrid[], double dh, int nvar, double vargrid[], double dvar, float *count)
*/
void models_2_addcount(int, float [], float [], int, double [], double, int, double [], double, float *);

/*
void synth_2_addcount(int npts, float synth[], int namp, double ampgrid[], double damp, float *count)
*/
void synth_2_addcount(int, float[], int, double [], double, float *);

/*
void count_2_ppd_norm(int nh, int nvar, float *count)
*/
void count_2_ppd_norm(int, int, float *);

/*
void count_2_ppd(int nh, int nvar, float *count, int nmodels) 
*/
void count_2_ppd(int, int, float *, int);

/*
void get_median_struct(int nh, int nvar, double vargrid[], float ppd[], double *model_med) 
*/
void get_median_struct(int, int, double [], float [], double *);

/*
void get_maximumdensity_struct(int nh, int nvar, double vargrid[], float ppd[],  double *model_maxppd)
*/
void get_maximumdensity_struct(int, int, double [], float [],  double *);

#endif
