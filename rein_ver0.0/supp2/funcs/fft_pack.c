#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft_pack.h"

/* create base trace */
/* calculate number of nfft and in the main part, we should create zero series of nfft length. */
int trace_base(double dt, double t_length)
{
  /* time = t_pre + t_length + t_pre (additional tail) */
  int nlength, nfft;

  nfft = 1;
  nlength = (int) (t_length)/dt;
  while( nfft < nlength){ nfft *=2; }
  
  return nfft;
}

/* create waveform in frequency domain */
void trace_fft(double *xtrace, int nfft, gsl_complex *xfft)
{
  int i;
  double * xtrace_double = (double *)alloca(sizeof(double)*nfft*2);

  for(i=0; i < nfft; i++)
  {
    *(xtrace_double +2*i+0) = *(xtrace+i);
    *(xtrace_double +2*i+1) =0.;
  }

  gsl_fft_complex_radix2_forward(xtrace_double, 1, nfft);

  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(xfft+i, *(xtrace_double+2*i+0), *(xtrace_double+2*i+1));
  }

  return;
}

/* do inverse fourier transform: make into time series */
/* output: xtrace */
void trace_inv_fft(gsl_complex *xfft, int nfft, double *xtrace)
{
  double *xtrace_double;
  int i;
  xtrace_double = (double *)alloca(sizeof(double)*2*nfft);

  for(i=0; i < nfft; i++)
  {
    *(xtrace_double+2*i+0) = GSL_REAL(*(xfft+i));
    *(xtrace_double+2*i+1) = GSL_IMAG(*(xfft+i));
  }
  
  gsl_fft_complex_radix2_inverse(xtrace_double, 1, nfft);
  
  for(i=0; i < nfft; i++)
  {
    *(xtrace+i) = *(xtrace_double+2*i);
  }

  return;
}

/* multiply amplitude to waveform */
/* overwriting the trace with multiplying "amp" */
void trace_mul_amp(gsl_complex *xfft, double amp, int nfft)
{
  int i;
  for(i=0; i < nfft; i++)
  {
    *(xfft+i) = gsl_complex_mul_real(*(xfft+i), amp);
  }

  return;
}

/* give time shift to waveform */
/* it is shifting to the right
 * and multiply the given amplitude coeff. (if not, set amp=1
 * we overwrite the sequence...
 */
void trace_shift_amp(double tshift, double dt, int nfft, double amp, gsl_complex *xfft)
{
  int i;
  double phase = 2*M_PI/(nfft*dt)*tshift*(-1.); /* df*(-1)*tshift where df = 2*M_PI/(nfft*dt) */
  double dc = cos(phase); double ds = sin(phase);
  double cos_0, sin_0, cos_1, sin_1;
  gsl_complex tmpcomp;

  cos_0 = 1; sin_0 = 0; /* i = 0 */
  for(i=0; i < nfft/2; i++)
  {
    GSL_SET_COMPLEX(&tmpcomp, cos_0, sin_0);
    cos_1 = cos_0*dc - sin_0*ds;
    sin_1 = sin_0*dc + cos_0*ds;

    *(xfft+i) = gsl_complex_mul(*(xfft+i), tmpcomp);
    *(xfft+i) = gsl_complex_mul_real(*(xfft+i), amp);
    cos_0 = cos_1; sin_0 = sin_1;  
  }

  cos_0 = cos(M_PI*tshift/dt); sin_0 = sin(M_PI*tshift/dt); /* nfft/2*df*tshift */
  for(i= nfft/2; i < nfft; i++)
  {
    GSL_SET_COMPLEX(&tmpcomp, cos_0, sin_0);
    cos_1 = cos_0*dc - sin_0*ds;
    sin_1 = sin_0*dc + cos_0*ds;

    *(xfft+i) = gsl_complex_mul(*(xfft+i), tmpcomp);
    *(xfft+i) = gsl_complex_mul_real(*(xfft+i), amp);
    cos_0 = cos_1; sin_0 = sin_1;
  }

  return;
}

/* prepare the differentiation factor in advance */
void trace_differential(double dt, int nfft, gsl_complex *diff)
{
  int i;
  double tmpdouble, dw;
  dw = (M_PI*2/dt)/nfft;

  for(i=0; i < nfft; i++)
  {
    tmpdouble = 0;
    if( i < nfft/2 +1){ tmpdouble = dw*i; }
    else{ tmpdouble = (i-nfft)*dw; }
    GSL_SET_COMPLEX( diff + i, 0, tmpdouble);
  }

  return; 
}

/* \omega^2 instead of i\omega -- helpful for ACFs */
void trace_differential_double(double dt, int nfft, double *diff)
{
  int i;
  double tmpdouble, dw;
  dw = (M_PI*2/dt)/nfft;

  for(i=0; i < nfft; i++)
  {
    tmpdouble = 0;
    if( i < nfft/2 +1){ tmpdouble = dw*i; }
    else{ tmpdouble = (i-nfft)*dw; }
    diff[i] = tmpdouble*tmpdouble; 
  }

  return;
}

/* prepare the trace shift factor in advance */
void trace_shift_factor(double tshift, double dt, double amp, int nfft, gsl_complex *shift)
{
  int i;
  double phase = 2*M_PI/(nfft*dt)*tshift*(-1.); /* df*(-1)*tshift where df = 2*M_PI/(nfft*dt) */
  double dc = cos(phase); double ds = sin(phase);
  double cos_0, sin_0, cos_1, sin_1;

  cos_0 = 1; sin_0 = 0; /* i = 0 */
  for(i=0; i < nfft/2; i++)
  {
    GSL_SET_COMPLEX(shift+i, cos_0, sin_0);
    *(shift+i) = gsl_complex_mul_real(*(shift+i), amp);

    cos_1 = cos_0*dc - sin_0*ds;
    sin_1 = sin_0*dc + cos_0*ds;   
    cos_0 = cos_1; sin_0 = sin_1;
  }

  cos_0 = cos(M_PI*tshift/dt); sin_0 = sin(M_PI*tshift/dt); /* nfft/2*df*tshift */
  for(i= nfft/2; i < nfft; i++)
  {
    GSL_SET_COMPLEX(shift+i, cos_0, sin_0);
    *(shift+i) = gsl_complex_mul_real(*(shift+i), amp);

    cos_1 = cos_0*dc - sin_0*ds;
    sin_1 = sin_0*dc + cos_0*ds;

    cos_0 = cos_1; sin_0 = sin_1;
  }

  return;
}

void trace_shift_with_factor(int nfft, gsl_complex *factor, gsl_complex *xfft)
{
  int i;

  for(i=0; i < nfft/2; i++)
  {
    *(xfft+i) = gsl_complex_mul(*(xfft+i), factor[i]);
  }

  for(i= nfft/2; i < nfft; i++)
  {
    *(xfft+i) = gsl_complex_mul(*(xfft+i), factor[i]);
  }

  return;
}


/* add two signals in frequency domain */
void trace_2add(int nfft, gsl_complex *x1, gsl_complex *x2)
{ // x2 = x1+x2
  int i;
  gsl_complex xtmp;

  for(i=0; i < nfft; i++)
  {
    xtmp = gsl_complex_add(*(x2+i), *(x1+i));
    *(x2+i) = xtmp;
  }

  return;
}

/* initialize trace */
void trace_initialize(int npts, gsl_complex *xfft)
{
  int i;
  for(i=0; i < npts; i++)
  {
    GSL_SET_COMPLEX(xfft+i, 0, 0);
  }

  return;
}

/* copy trace */
void trace_copy(int npts, gsl_complex *xorig, gsl_complex *xcopy)
{
  int i;
  for(i=0; i <npts; i++)
  {
    *(xcopy+i) = *(xorig+i);
  }
  return;
}
