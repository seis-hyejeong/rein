#ifndef _acf_functions_h_
#define _acf_functions_h_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "complex_pack_2.h"
#include "math.h"

/* double gaussian(double awidth, int nfft, double deltat, gsl_complex *xarr1, gsl_complex *xarr2, gsl_complex *xarr3) */
double gaussian(double, int, double, gsl_complex *, gsl_complex *, gsl_complex *);
/* double gaussian_1comp(double awidth, int nfft, double deltat, gsl_complex *xarr1) */
double gaussian_1comp(double, int, double, gsl_complex *);

/* void autocorr(int nfft, gsl_complex *xarr) */
double autocorr(int, gsl_complex *);

/* considering differentitaion in the autocorrelation 
double autocorr_ver2(int nfft, gsl_complex *xarr, bool diff_tf, double *wsquare)
*/
double autocorr_ver2(int, gsl_complex *, bool, double *);

/* the noise added, no whitening ACF 
-- we already have the autocorrelation calculated, but just need to add the different noises at a time.
void autocorr_noise_byadding(int nfft, double square_sum, double snr_factor, gsl_complex *xautocorr, gsl_complex *newautocorr)
*/
void autocorr_noise_byadding(int, double, double, gsl_complex *, gsl_complex *);


/* void autocorr_noise(int nfft, double snr_factor, gsl_complex *xarr) */
void autocorr_noise(int, double, gsl_complex *);
/* void autocorr_smooth_p2017(int nfft, double deltat, double smooth, gsl_complex *xarr) */
void autocorr_smooth_p2017(int, double, double, gsl_complex *);

/*
void autocorr_smooth_t2019_proto(int nfft, double deltat, double smooth, gsl_complex xautocorr[], double *smoothing)
*/
void autocorr_smooth_t2019_proto(int, double, double, gsl_complex [], double *);

/*
void autocorr_smooth_t2019_proto_ver2(int nfft, double deltat, double smooth, gsl_complex xautocorr[], double *smoothing)
*/
void autocorr_smooth_t2019_proto_ver2(int, double, double, gsl_complex [], double *);

/*
void autocorr_smooth_t2019_final_noise(int nfft, gsl_complex xautocorr[], double smoothing[], double square_sum, double snr_factor, gsl_complex *newautocorr)
*/
void autocorr_smooth_t2019_final_noise(int, gsl_complex [], double [], double, double, gsl_complex *);

/* new version for differentiation */
/*
In case of no smoothing
void autocorr_noise_byadding_ver3(int nfft, double square_sum, double snr_factor, gsl_complex *xautocorr, gsl_complex *newautocorr, bool diff_tf, double *wsquare)
*/
void autocorr_noise_byadding_ver3(int, double, double, gsl_complex *, gsl_complex *, bool, double *);

/* when there is smoothing considered
void autocorr_smooth_t2019_final_noise_ver3(int nfft, gsl_complex xautocorr[], double smoothing[], double square_sum, double snr_factor, gsl_complex *newautocorr, bool diff_tf, double *wsquare)
*/
void autocorr_smooth_t2019_final_noise_ver3(int, gsl_complex [], double [], double, double, gsl_complex *, bool, double *);

/*
void autocorr_smooth_t2019_final_noise_ver2(int nfft, gsl_complex xautocorr[], double smoothing[], double square_sum, double snr_factor, gsl_complex *newautocorr)
 */
void autocorr_smooth_t2019_final_noise_ver2(int, gsl_complex [], double [], double, double, gsl_complex *);

/* void ratio_spectral(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio) */
void ratio_spectral(int, gsl_complex *, gsl_complex *, gsl_complex *);
/* void ratio_simplediv(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio) */
void ratio_simplediv(int, gsl_complex *, gsl_complex *, gsl_complex *);

/* void hanning(int npts, float *window) */
void hanning(int, float *);
/* float area_func(int npts, float func[]) */
float area_func(int, float []);
/* void autocorr_smooth_func(int nfft, int nint_half, float *window, gsl_complex *xarr) */
void autocorr_smooth_func(int, int, float *, gsl_complex *);

#endif
