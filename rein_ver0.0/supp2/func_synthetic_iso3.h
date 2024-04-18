#ifndef _func_synthetic_iso3_h_
#define _func_synthetic_iso3_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

#include "sac.h"
#include "chris_only.h"
#include "calculation_w.h"
#include "complex_pack_2.h"
#include "matrix_pack.h"
#include "acf_functions.h"
#include "fft_pack.h"
#include "parameters.h"

#define bandwidth 0.
#define attenuation 0.
#define order 2
#define pass 2

/*
void calculate_stack(int ntrace, int npts, float xindv[], float *xstack)
*/
void calculate_stack(int, int, float [], float *);

/*
void calculate_stack_stddev(int ntrace, int npts, float xindv [], float *xstack, float *stddev)
*/ 
void calculate_stack_stddev(int, int, float [], float *, float *);

/* 
int calculate_covariance(float waterlevel, int ntrace, int npts, float xindv[], float xstack[], gsl_matrix *cov, gsl_matrix *L );
int calculate_covariance(float waterlevel, int ntrace, int npts, float xindv[], float xstack[], gsl_matrix *cov, gsl_matrix *L, gsl_permutation *permute)
*/
int calculate_covariance(int, int, float [], float [], gsl_matrix *, gsl_matrix *, gsl_permutation *);

/*
void stack_synths(int ntraces, int npts, double *xsynth); 
*/
void stack_synths(int, int, float *);

/*
void average_synths(int ntraces, int npts, double *stacked_xsynth)   
*/
void average_synths(int, int, float *);

/*
void Ntimes_trace(int ntraces, int npts, double *stacked_obs)
*/
void Ntimes_trace(int, int, float *);

/* generating synthetic noise 
int synthetic_noises_freq(gsl_rng *r, int nerr, int nfft, double sigma, gsl_complex *xerr, double deltat, double gaussf, char savename[strlng])
*/
int synthetic_noises_freq(gsl_rng *, int, int, double, gsl_complex *, double, double, char []);

/*
Applying gaussian filter to the synthetic trace in the frequency domain.
double gaussian_lp(double awidth, int nfft, double deltat, gsl_complex *xarr1)
*/
double gaussian_lp(double, int, double, gsl_complex *);
/* calculating the square sume of frequency spectra
double sigma_complex_abs(int npts, gsl_complex *xtrace_c)
*/ 
double sigma_complex_abs(int, gsl_complex *);

/* add the square sum 
double sigma_complex_abs2(int npts, gsl_complex *xtrace_c)
*/
double sigma_complex_abs2(int, gsl_complex *);

/* adding certain real number constant throughout the complex array
void add_complex_double_series(int npts, gsl_complex *xtrace_c, double add)
*/
void add_complex_double_series(int, gsl_complex *, double);

/*
int calculate_synthetics_noise_ver6(PARAM inparam, float *s_racf, float *indv_racf, INDATA racf, float *s_zacf, float *indv_zacf, INDATA zacf, float *s_pacf, float *indv_pacf, INDATA pacf, float *s_rstack, float *s_zstack, INDATA rstack, MODEL m, SYNTHPAR synth, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw, gsl_complex *shift, double *w2array)
*/
int calculate_synthetics_noise_ver6(PARAM, float *, float *, INDATA, float *, float *, INDATA, float *, float *, INDATA, float *, float *, INDATA, MODEL, SYNTHPAR, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *, gsl_complex *, double *);

/*
void set_christoffel_iso_water(int nlayer, float mod_h[], float mod_rho[], float mod_p[], float mod_s[], double *chris)
*/
void set_christoffel_iso_water(int, float [], float [], float [], float [], double *);

/*
void set_elastic_anisotropic(float rho, float alpha, float beta, float B, float C, float E, float azi, float tilt, double *o_tensor)
*/
void set_elastic_anisotropic(float, float, float, float, float, float, float, float, double *);

/*
void set_elastic_isotropic(float rho, float alpha, float beta, double *o_tensor)
*/
void set_elastic_isotropic(float, float, float, double *);

/*
int get_nfft(double time, double dt);
*/
int get_nfft(double, double);

/*
void set_parameters(int nfft, double dt, double *t0_rf, double *t0_acf, double *dw)
*/
void set_parameters(int, double, double *, double *, double *);

/*
void get_zerotindex(double t0_rf, double t0_acf, double dt, int *ind0_rf, int *ind0_acf);
*/
void get_zerotindex(double, double, double, int *, int *);

/*
double convert_gauss(double gauss)
*/
double convert_gauss(double);

/*
int filter_signal(int npts, double dt, double low, double high, double *xtrace)
*/
int filter_signal(int, double, double, double, float *);

/* get max 
double get_abs_max(int npts, double *trace)
*/
double get_abs_max(int, double *);

/*
float normalize_max(int npts, float *trace, int ind0, int windownpts)
*/
float normalize_max(int, float *, int, int);

/* normalize to the number that you want to scale -- not one, but to a number given
float normalize_max_2_new(int npts, float *trace, int ind0, double newmax, int windownpts)
*/
float normalize_max_2_new(int, float *, int, double, int);

/*
int normalize_RZ_ver2(int npts, double *ztrace, double *rtrace, int ind0)
*/
int normalize_RZ_ver2(int, float *, float *, int);

/*
int normalize_single(int npts, double *trace, int ind0)
*/
int normalize_single(int, float *, int);

/*
int make_synthetic_3comp_iso_ver2(int iinput, int nray, double i_rayp[maxray], int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, double dt, double t0, double gausst, float *xtrace_r, float *xtrace_z, int inorm, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw)
*/
int make_synthetic_3comp_iso_ver2(int, int, double [], int, float [], float [], double [], int, float, int, double, double, double, double, float *, float *, int, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *);


/*
int make_synthetic_Rcomp(int iinput, int nray, double i_rayp[maxray], int nlayer, float m_thick[], float m_rho [], double m_chris[], int iwater, float water_alpha, int nfft, double dw, double dt, double t0, double gausst, float *xtrace_r, int inorm)
*/
int make_synthetic_Rcomp_iso(int, int, double [maxray], int, float [], float [], double [], int, float, int, double, double, double, double, float *, int, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *);

/*
int make_synthetic_ACF_3comp(bool ivelocity, int iinput, int nray, double i_rayp[maxray], int nlayer, float m_thick[], float m_rho [], double m_chris[], int iwater, float water_alpha, int nfft, double dw, double dt, double gaussf, double wsmooth, float *xtrace_r, float *xtrace_t, float *xtrace_z);
*/
int make_synthetic_ACF_3comp_iso(bool, int, int, double [maxray], int, float [], float [], double [], int, float, int, double, double, double, double, float *, float *, float *, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *, gsl_complex *);


/*
int make_synthetic_ACF_R_iso_noise_ver6(int iinput, INDATA racf, int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, float *xtrace_r, float *indv_r, int synthind_0, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw, gsl_complex *shift, double *w2array)
*/
int make_synthetic_ACF_R_iso_noise_ver6(int, INDATA, int, float [], float [], double [], int, float, int, double, float *, float *, int, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *, gsl_complex *, double *);

/*
int make_synthetic_ACF_Z_iso_noise_ver6(int iinput, INDATA zacf, int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, float *xtrace_z, float *indv_z, int synthind_0, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw,  gsl_complex *shift, double *w2array)
*/
int make_synthetic_ACF_Z_iso_noise_ver6(int, INDATA, int, float [], float [], double [], int, float, int, double, float *, float *, int, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *, gsl_complex *, double *);

/* same for the PACF */
int make_synthetic_ACF_DPG_iso_noise_ver6(int, INDATA, int, float [], float [], double [], int, float, int, double, float *, float *, int, gsl_matrix_complex **, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_matrix_complex *, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *, gsl_complex *, gsl_complex *, double *);

#endif
