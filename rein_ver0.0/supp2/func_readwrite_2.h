#ifndef _func_readwrite_2_h_
#define _func_readwrite_2_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <ctype.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>

#include "sac_hk.h"
#include "parameters.h"
#include "func_synthetic_iso3.h"
#include "func_models.h"

/* prototype of the functions */
/* void initialize_param(PARAM *param)*/
void initialize_param(PARAM *);

/* void initialize_prior(PRIOR priorparam) */
void initialize_prior(PRIOR *);

/*
void initialize_orgdata(DATAORG *orgdata)
*/
void initialize_orgdata(DATAORG *);

/* void initialize_outparam (OUTPARAM outparam) */
void initialize_outparam(OUTPARAM *);

/*
void initialize_indata(INDATA inparam)
*/
void initialize_indata(INDATA *);

/* 
void initialize_paralleltempering(PARATEMP *);
*/
void initialize_paralleltempering(PARATEMP *);

/* 
int take_parameters_ver2(char infilename[strlng], char inmodel[strlng], PARAM *parameter, DATAORG *orgdata, PARATEMP *parallelt, INDATA *racfdata, INDATA *zacfdata, INDATA *pacfdata, INDATA *rwavedata, OUTPARAM *outparam)
*/
int take_parameters_ver2(char [strlng], char [strlng], PARAM *, DATAORG *, PARATEMP *, INDATA *, INDATA *, INDATA *, INDATA *, OUTPARAM *);

/* initialize temperatures
int get_parallel_temperatures(PARATEMP *parallelt, gsl_rng *r)
*/
int get_parallel_temperatures(PARATEMP *, gsl_rng *);

/*
int take_data_input_ver2(int my_id, INDATA *data, DATAORG orgdata, char outroot[strlng], float **stack, float **error, float **indvtrace, gsl_matrix **L, gsl_vector **A)
*/
int take_data_input_ver3(int, INDATA *, DATAORG, char [strlng], float **, float **, float **);
//int take_data_input_ver2(int, INDATA *, DATAORG, char [strlng], float **, float **, float **, gsl_matrix *, gsl_vector *);

/*
int waterlevel_error(int my_id, INDATA *data, float *error)
*/
int waterlevel_error(int, INDATA *, float *);

/*
int generate_covariance_LU(int my_id, char outroot[strlng], INDATA data, float *indvtrace, float *stack, gsl_matrix *L, gsl_permutation *permute);
*/
int generate_covariance_LU(int, char [strlng], INDATA, float *, float *, gsl_matrix *, gsl_permutation *);

/*
int generate_covariance_ldlt(int my_id, char outroot[strlng], INDATA data, float *indvtrace, float *stack, gsl_matrix *L)
*/
int generate_covariance_ldlt(int, char [strlng], INDATA, float *, float *, gsl_matrix *);

/*
int taper_covariance(INDATA data, gsl_matrix *L)
*/
int taper_covariance(INDATA, gsl_matrix *);

/*
double taper_covariance_model(INDATA data, gsl_matrix *L)
*/
double taper_covariance_model(INDATA, gsl_matrix *);

/*
int organize_rayp_snrs(DATAORG orgdata, INDATA *data, char outroot[strlng], int *prep_inclbin, int *prep_ind0, int *prep_bincount, int *prep_fileindex) 
*/
int organize_rayp_snrs(DATAORG, INDATA *, char [strlng], int *, int *, int *, int *);

/*
void calculate_snr_factors(INDATA *data, SYNTHPAR synthpar)
*/
void calculate_snr_factors(INDATA *, SYNTHPAR);

/* when there's velocity consideration 
void calculate_snr_factors_vel(INDATA *data, SYNTHPAR synthpar, double *wsquare)
*/
void calculate_snr_factors_vel(INDATA *, SYNTHPAR, double *);

/*
int initialize_synth_ver2(SYNTHPAR *synthpar, INDATA racfdata, INDATA zacfdata, INDATA pacfdata, INDATA rwavedata)
*/
int initialize_synth_ver2(SYNTHPAR *, INDATA, INDATA, INDATA, INDATA);

/*
int read_inputvelmod(PARAM param, PRIOR prior, char velmodname[strlng], MODEL m)
*/
int read_inputvelmod(PARAM, PRIOR, char [strlng], MODEL *);

/*
int read_inputpriormod_ver2(PARAM param, PRIOR *prior, MODEL *m)
*/
int read_inputpriormod_ver2(PARAM, PRIOR *, MODEL *);

/*
void change_array_to_printable(OUTPARAM *o_param)
*/
void change_array_to_printable(OUTPARAM *);

/*
void change_array_to_printable_single(int npts, double *array, double delta)
*/
void change_array_to_printable_single(int, double *, double);

/*
int print_ppdmodels2(int nh, double harray[], double vp[], double vs[], double rho[], double kappa[], char filename[])
*/
int print_ppdmodels2(int, double [], double [], double [], double [], double [], char []);
/*
void print_array_variables2(OUTPARAM o_param)
*/
void print_array_variables2(OUTPARAM);

/* print the probability density function 
void print_pdf_variable(OUTPARAM o_param, char outfilename[], int nvar, double *vararray, float *pdf_var)
*/
void print_pdf_variable(OUTPARAM, char [], int, double *, float *);

/* 
void print_pdf_synthetics(char outfilename[], int ntime, double t0, double dt, float intrace[], int namp, double amparr[], float pdf[]);
*/
void print_pdf_synthetics(char [], int, double, double, float [], int, double [], float []);

#endif
