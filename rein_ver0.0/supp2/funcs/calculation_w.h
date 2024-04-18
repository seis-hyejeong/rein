#ifndef _calculation_w_h_
#define _calculation_w_h_

#include <gsl/gsl_poly.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "parameters-2.h"
#include "matrix_pack.h"
#include "fft_pack.h"

void mod_elastic(int, double [], double [], double *, double *);
int qsols2(double , double [], double, double *);
/* isotropic */
int qsols2_iso(double, double [], double, double *);

void oscillvec(double [], double [], double [], double [], double [], double *, int);
/*void E_matrix(double C[], double rayp, double qz[], double dvec[], gsl_matrix_complex *Emat) */
void E_matrix(double [], double, double [], double [], gsl_matrix_complex *);
/* void Einv_w_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *Einvmat) */
void Einv_w_matrix(gsl_matrix_complex *, gsl_matrix_complex *);
/* void Einv_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *Einvmat) */
void Einv_matrix(gsl_matrix_complex *, gsl_matrix_complex *);
/* void A_w_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *E_w_invmat, gsl_complex [], gsl_matrix_complex *A, gsl_matrix_complex *Dw_w) */
void A_w_matrix(gsl_matrix_complex *, gsl_matrix_complex *, gsl_complex [], gsl_matrix_complex *, gsl_matrix_complex *);
/* void A_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *Einvmat, gsl_complex expw[], gsl_matrix_complex *A, gsl_matrix_complex *Dw) */
void A_matrix(gsl_matrix_complex *, gsl_matrix_complex *, gsl_complex [], gsl_matrix_complex *, gsl_matrix_complex *);
/* void cal_expw(int nlayer, double qzs[], float dms[], double dw, int nfft, gsl_complex *expw) */
void cal_expw(int, double [], float [], double, int, gsl_complex *);
/* double cal_poly(int order, double *coeffs, double x) */
double cal_poly(int, double *, double);
/* void Jw_matrix(int nlayer, gsl_matrix_complex *invEn, gsl_matrix_complex **A, gsl_matrix_complex *J, gsl_matrix_complex *Jtmp) */
void Jw_matrix(int, gsl_matrix_complex *, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *); 
/* void Jw_matrix2(int nlayer, gsl_matrix_complex **A, gsl_matrix_complex *J, gsl_matrix_complex *Jtmp) */
void Jw_matrix2(int, gsl_matrix_complex **, gsl_matrix_complex *, gsl_matrix_complex *);
/* void get_response(int iinci, gsl_matrix_complex *Jwmat, gsl_complex *surf_u, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p) */
void get_response(int, gsl_matrix_complex *, gsl_complex *, gsl_matrix_complex*, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *);
/* void get_response_w(int iinci, gsl_matrix_complex *Jwmat, gsl_complex *surf_u, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p) */
void get_response_w(int, gsl_matrix_complex *, gsl_complex *, gsl_matrix_complex*, gsl_vector_complex *, gsl_vector_complex *, gsl_permutation *);

void conv_source(gsl_complex *, gsl_complex *, gsl_complex *, gsl_complex [], int);
void source(int, double, double, gsl_complex *);
void output_invrot(int, double [], double [], double []);
double p_tt(int, double [], double []);

#endif
