#ifndef _matrix_pack_h_
#define _matrix_pack_h_

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include "complex_pack_2.h"
#include "parameters-2.h"

void tensor_set_all(double *, double);
void tensor_set(double *, const int, const int, const int, const int, const double);
void matrix_set(int, double *, int, int, double);
double tensor_get(const double *, const int, const int, const int, const int);
double matrix_get(int, double *, int, int);
void tensor_2_mat(double *, double *);
int conv_ind(int, int, int, int, int *, int *);
int conv_ind_op(int, int, int *, int *, int *, int *);
void ext_array(int, double *, double *);
int conv(double *, int, double *, int, double *);
void sum_v(int, double *, double *);
void amul_v(int, double *, double);

void mat_vec_mul(double *, double *, int, int, double *);
void cmat_vec_mul(gsl_matrix_complex *, gsl_complex [], int, int, gsl_complex *);
void mat_mat_mul(double [], double [], int, int, int, double *);
void cmat_mat_mul(gsl_complex [], gsl_complex [], int, int, int, gsl_complex *);
/*
void cmat_mat_mul2(gsl_complex A[], gsl_complex B[], int m, int n, int k, gsl_complex *C)
*/
void cmat_mat_mul2(gsl_complex [], gsl_complex [], int, int, int, gsl_complex *);
void mat_rot(double *, double *, double *);
void mat_rot_inv(double *, double *, double *);
double vdot(double [], double [], int);
double vabs(double [], int);
double vabs2(double [],int);
void vcross(double [], double [], double *);
void mat_inv_solve(int, double [], double [], double *);
void set_elastic(double, double, double, double, double, double, double, double, double *);
void get_aaxis_rotmat(double, double, double *);
void get_2d_rotmat(double, double *);
double get_rothor_angle(double qvec[]);
void get_vel_perturb(double, double, double, double, double, double, double *);
void rot_tensor(double *, double *, double *);

#endif 
