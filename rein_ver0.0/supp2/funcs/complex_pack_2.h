#ifndef _complex_pack_2_h_
#define _complex_pack_2_h_

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

gsl_complex cadd(gsl_complex, gsl_complex);
gsl_complex csub(gsl_complex, gsl_complex);
gsl_complex cmul(gsl_complex, gsl_complex);
gsl_complex cdiv(gsl_complex, gsl_complex);
gsl_complex amul(double, gsl_complex);
gsl_complex adiv(double, gsl_complex);
double c_abs(gsl_complex);

#endif
