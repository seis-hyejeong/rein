#ifndef _func_likelihood_h_
#define _func_likelihood_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

#include "parameters.h"

/*
 void model_likelihood_log_sigmasquare(int ndata, float sigma[], float *inv_sigma)
*/
void model_likelihood_log_sigmasquare(int, float [], float *);
/*
double model_likelihood_log_term2(int ndata, float synth[], float obs[], float inv_sigma[]) 
*/
double model_likelihood_log_term2(int, float [], float [], float []);

/*
double model_likelihood_log_covmat(int ndata, float synth[], float obs[], gsl_matrix *L, gsl_vector *vtmp, gsl_vector *ytmp)
*/
double model_likelihood_log_covmat_ldlt(int, float [], float [], gsl_matrix *, gsl_vector *, gsl_vector *);

/*
double model_likelihood_log_covmat2(int ndata, float synth[], float obs[], gsl_matrix *L, gsl_permutation *per, gsl_vector *vec, gsl_vector *yvec)
*/
double model_likelihood_log_covmat2(int, float [], float [], gsl_matrix *, gsl_permutation *, gsl_vector *, gsl_vector *);

/*
void model_likelihood_log_alpha(LIKELI p)
*/
void model_likelihood_log_alpha(LIKELI *);

/* Write the likelihood
int print_likelihood(int ith_cal, FILE *filestr, PARATEMP parallelt)
*/
int print_likelihood(int, FILE *, PARATEMP);

/*
int print_likelihood_T(int ith_cal, FILE *filestr, PARATEMP parallelt, double printT)
*/
int print_likelihood_T(int, FILE *, PARATEMP, double);

/*
int swap_temperature(PARATEMP *parallelt,  gsl_rng *r)
*/
int swap_temperature(PARATEMP *,  gsl_rng *);

/*
void synthetics_likelihood_ver3(LIKELI *p, SYNTHPAR s_param, float s_racf[], float indv_racf[], float o_racf[], float o_racf_indv[], float inv_racf[], gsl_matrix *mat_racf, gsl_vector *vec_racf, INDATA i_racf, float s_zacf[], float indv_zacf[], float o_zacf[], float o_zacf_indv[], float inv_zacf[], gsl_matrix *mat_zacf, gsl_vector *vec_zacf, INDATA i_zacf, float s_pacf[], float indv_pacf[], float o_pacf[], float o_zacf_indv[], float inv_pacf[], gsl_matrix *mat_pacf, gsl_vector *vec_pacf, INDATA i_pacf, float s_rwave[], float o_rwave[], float inv_rwave[], gsl_matrix *mat_zacf, gsl_vector *vec_rwave, INDATA i_rwave);
*/
void synthetics_likelihood_ver3_ldlt(LIKELI *, SYNTHPAR, float [], float [], float [], float [], float [], gsl_matrix *, gsl_vector *, gsl_vector *, INDATA, float [], float [], float [], float [], float [], gsl_matrix *, gsl_vector *, gsl_vector *, INDATA, float [], float [], float [], float [], float [], gsl_matrix *, gsl_vector *, gsl_vector *, INDATA, float [], float [], float [], gsl_matrix *, gsl_vector *, gsl_vector *, INDATA);

/* 
void synthetics_likelihood_ver3_lu(LIKELI *p, SYNTHPAR s_param, float s_racf[], float indv_racf[], float o_racf[], float o_racf_indv[], float inv_racf[], gsl_matrix *mat_racf, gsl_permutation *per_racf, gsl_vector *vec_racf, gsl_vector *y_racf, INDATA i_racf, float s_zacf[], float indv_zacf[], float o_zacf[], float o_zacf_indv[], float inv_zacf[], gsl_matrix *mat_zacf, gsl_permutation *per_zacf, gsl_vector *vec_zacf, gsl_vector *y_zacf, INDATA i_zacf, float s_pacf[], float indv_pacf[], float o_pacf[], float o_pacf_indv[], float inv_pacf[], gsl_matrix *mat_pacf, gsl_permutation *per_pacf, gsl_vector *vec_pacf, gsl_vector *y_pacf, INDATA i_pacf, float s_rwave[], float o_rwave[], float inv_rwave[], gsl_matrix *mat_rwave, gsl_permutation *per_rwave, gsl_vector *vec_rwave, gsl_vector *y_rwave, INDATA i_rwave)
*/
void synthetics_likelihood_ver3_lu(LIKELI *, SYNTHPAR, float [], float [], float [], float [], float [], gsl_matrix *, gsl_permutation *, gsl_vector *, gsl_vector *, INDATA, float [], float [], float [], float [], float [], gsl_matrix *, gsl_permutation *, gsl_vector *, gsl_vector *, INDATA, float [], float [], float [], float [], float [], gsl_matrix *, gsl_permutation *, gsl_vector *, gsl_vector *, INDATA, float [], float [], float [], gsl_matrix *, gsl_permutation *, gsl_vector *, gsl_vector *, INDATA);

/*
void synthetics_likelihood_norm(LIKELI p, SYNTHPAR s_param, float s_racf[], float o_racf[], float inv_racf[], INDATA i_racf, float s_zacf[], float o_zacf[], float inv_zacf[], INDATA i_zacf, double s_rwave[], float o_rwave[], float inv_rwave[], INDATA i_rwave)
*/
void synthetics_likelihood_norm(LIKELI *, SYNTHPAR, float [], float [], float [], INDATA, float [], float [], float [], INDATA, float [], float [], float [], INDATA);

#endif
