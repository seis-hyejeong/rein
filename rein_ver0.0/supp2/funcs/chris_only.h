#ifndef _chris_only_h_
#define _chris_only_h_

#include "parameters-2.h"

void get_struct_grotmat(double, double, double *);
int cal_oscillvec_eigen(double const [], double const, double const, double *, double *);
int cal_oscillvec_iso(double const [], double const, double const, double *, double *);
void F_matrix(double const [], double const [], double *);
int get_3dvec(double const, double *);
int get_vel(int const, double const [], double const [], double const, double *, double *);
void change_units(double, int, double *);

#endif
