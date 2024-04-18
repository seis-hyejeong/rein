#include <math.h>
#include "complex_pack_2.h"

/* simplify some complex number operations...! *
 * it is using gsl complex routine, 
 * but only to simplify the written expression... */

/* simple arithmatic */
gsl_complex cadd(a, b)
gsl_complex a, b;
{
  return gsl_complex_add(a, b);
}

gsl_complex csub(a, b) /* c = a-b */
gsl_complex a, b;
{
  return gsl_complex_sub(a, b);
}

gsl_complex cmul(a, b)
gsl_complex a, b;
{
  return gsl_complex_mul(a, b);
}

gsl_complex cdiv(a, b) /* c = a/b */
gsl_complex a, b;
{
  return gsl_complex_div(a, b);
}

gsl_complex amul(a, x)
double a; 
gsl_complex x;
{
  return gsl_complex_mul_real(x, a);
}

gsl_complex adiv(a, x)
double a;
gsl_complex x;
{
  return gsl_complex_div_real(x, a);
}

double c_abs(a)
gsl_complex a;
{
  return gsl_complex_abs(a);
}
