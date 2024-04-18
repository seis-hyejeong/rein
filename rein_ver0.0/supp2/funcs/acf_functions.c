#include "acf_functions.h" 

double gaussian(double awidth, int nfft, double deltat, gsl_complex *xarr1, gsl_complex *xarr2, gsl_complex *xarr3)
{
  int i;
  double df, gaussamp, normalizer = 0;

  df = 2*M_PI/(nfft*deltat);
  for(i=0; i < nfft/2; i++)
  {
    gaussamp = exp(-pow(df*i,2)/(4*awidth*awidth)); 
    xarr1[i] = amul(gaussamp, xarr1[i]);
    xarr2[i] = amul(gaussamp, xarr2[i]);
    xarr3[i] = amul(gaussamp, xarr3[i]);
    normalizer += gaussamp;
  }
  for(i=nfft/2; i < nfft; i++)
  {
    gaussamp = exp(-pow(df*(i-nfft),2)/(4*awidth*awidth));
    xarr1[i] = amul(gaussamp, xarr1[i]);
    xarr2[i] = amul(gaussamp, xarr2[i]);
    xarr3[i] = amul(gaussamp, xarr3[i]);
    normalizer += gaussamp;
  }

  /* invert the normalizer */
  normalizer = nfft/normalizer;

  /* multiply the normalizer */
//  for(i=0; i < nfft; i++)
//  {
//    xarr1[i] = amul(normalizer, xarr1[i]);
//    xarr2[i] = amul(normalizer, xarr2[i]);
//    xarr3[i] = amul(normalizer, xarr3[i]);
//  }

  /* return the normalizing value */
  normalizer = 1./normalizer;

  return normalizer;
}

double gaussian_1comp(double awidth, int nfft, double deltat, gsl_complex *xarr1)
{
  int i, ind;
  double df, df2, gaussamp, normalizer = 0;

  df = 2*M_PI/(nfft*deltat);
  df2 = df*df;
  for(i=0; i < nfft/2; i++)
  {
    gaussamp = exp(-df2*i*i/(4*awidth*awidth));
    xarr1[i] = amul(gaussamp, xarr1[i]);
    normalizer += gaussamp;
  }
  for(i=nfft/2; i < nfft; i++)
  {
    ind = i - nfft;
    gaussamp = exp(-df2*ind*ind/(4*awidth*awidth));
    xarr1[i] = amul(gaussamp, xarr1[i]);
    normalizer += gaussamp;
  }

  /* invert the normalizer normalizer = nfft/normalizer; */

  /* return the normalizing value */
  normalizer = normalizer/nfft;

  return normalizer;
}

double autocorr(int nfft, gsl_complex *xarr)
{
  int i;
/*  gsl_complex xtmp, xconj; */
  double size, square_sum;

  square_sum = 0.;
  for(i=0; i < nfft; i++)
  {
/*    xconj = gsl_complex_conjugate(xarr[i]);
    xtmp = gsl_complex_mul(xconj, xarr[i]);
    xarr[i] = xtmp; 
*/
    size = gsl_complex_abs2(xarr[i]);    
    GSL_SET_COMPLEX(xarr+i, size, 0.);
    square_sum += size;
  }

  return square_sum;
}

/* version 2 autocorrelation -- can consider the case of differentiation */
double autocorr_ver2(int nfft, gsl_complex *xarr, bool diff_tf, double *wsquare)
{
  int i;
  double size, square_sum;

  square_sum = 0.;
  for(i=0; i < nfft; i++)
  {
    size = gsl_complex_abs2(xarr[i]);
    if( diff_tf ){ size *= wsquare[i]; }

    square_sum += size;
    GSL_SET_COMPLEX(xarr+i, size, 0.);
  }

  return square_sum;
}


/* make the noise added autocorrelation from the autocorr calculated */
void autocorr_noise_byadding(int nfft, double square_sum, double snr_factor, gsl_complex *xautocorr, gsl_complex *newautocorr)
{
  int i;
  double addfactor = 0.;

  if(snr_factor > 0.) /* if not, we are not using it --> noise free ! */
  {
    addfactor = square_sum*snr_factor;
    for(i=0; i < nfft; i++)
    {
      newautocorr[i] = gsl_complex_add_real(xautocorr[i], addfactor);
    }
  }
  else
  { /* if no noise, return the original ! */
    for(i=0; i < nfft; i++)
    {
      newautocorr[i] = xautocorr[i];
    }
  }
 
  return;
}

/* the snr factor is:
 snr_factor = 1./(npts*(snr*snr-1));
*/
void autocorr_noise(int nfft, double snr_factor, gsl_complex *xarr)
{
  int i;
  double size, square_sum;

  square_sum = 0.;
  if(snr_factor > 0.) /* if snr_factor < 0, we are not using it */
  {
    for(i= 0; i < nfft; i++)
    {
      square_sum += gsl_complex_abs2(xarr[i]);
    }
    square_sum *= snr_factor;
  }

  /* the actual acf */
  for(i=0; i < nfft; i++)
  {
    size = gsl_complex_abs2(xarr[i]);
    size += square_sum;
    GSL_SET_COMPLEX(xarr+i, size, 0.);
  }

  return;
}


/* autocorrelation with smoothing -- denominator is smoothed. */
/* the definition of pham  and tkalcic 2017 */
void autocorr_smooth_p2017(int nfft, double deltat, double smooth, gsl_complex *xarr)
{
  int i, j, nint, nelements, nfft_2;
  gsl_complex xtmp, xconj;
  double df, xxtmp, xxtmp2, *smoothing;

  smoothing = calloc(nfft, sizeof(double));
  df = 1./(nfft*deltat);
  nint = floor(smooth/df/2);
  nfft_2 = nfft/2;

  /* get smoothing weights */
  xxtmp = 0; nelements = nint + 1;
  for(j=0; j <= nint; j++)
  {
    xxtmp += gsl_complex_abs(xarr[j]);
  }
  xxtmp2 = xxtmp/nelements;
  smoothing[0] = xxtmp2*xxtmp2;

  for(i=1; i < nint+1; i++)
  {
    xxtmp += gsl_complex_abs(xarr[nint+i]);
    nelements++;
    xxtmp2 = xxtmp/nelements;  
    smoothing[i] = xxtmp2*xxtmp2;
  }

  nelements = 2*nint+1; /* constant for a while */
  for(i=nint+1; i< nfft_2-nint+1; i++)
  {
    xxtmp += gsl_complex_abs(xarr[nint+i]);
    xxtmp -= gsl_complex_abs(xarr[i-nint-1]);
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2*xxtmp2; 
  }

  for(i = nfft_2-nint+1; i <= nfft_2; i++) /* decreasing */
  {
    xxtmp -= gsl_complex_abs(xarr[i-nint-1]);
    nelements--; 
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2*xxtmp2;
  }

  /* negative frequency */
  xxtmp = 0; nelements = nint + 2;
  for(j= nfft_2; j <= nint+ nfft_2+1; j++)
  {
    xxtmp += gsl_complex_abs(xarr[j]);
  }
  xxtmp2 = xxtmp/nelements;
  smoothing[nfft_2+1] = xxtmp2*xxtmp2;

  for(i= nfft_2+2; i < nfft_2+ nint; i++)
  {
    xxtmp += gsl_complex_abs(xarr[nint+i]);
    nelements++; 
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2*xxtmp2;
  }

  /* constant for a while */
  nelements = 2*nint +1 ;
  for(i = nfft_2+ nint; i < nfft - nint; i++)
  {
    xxtmp += gsl_complex_abs(xarr[i+nint]);
    xxtmp -= gsl_complex_abs(xarr[i-nint-1]);

    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2*xxtmp2;
  }

  for(i = nfft-nint; i < nfft; i++) /* decreasing */
  {
    xxtmp -= gsl_complex_abs(xarr[i-nint-1]);
    nelements--;
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2*xxtmp2;
  }
 
  for(i=0; i < nfft; i++)
  {
    xconj = gsl_complex_conjugate(xarr[i]);
    xtmp = gsl_complex_mul(xconj, xarr[i]);
    /* now smoothing!! */
    xarr[i] = adiv(smoothing[i], xtmp);
  }

  free(smoothing);
  return;
}

/* the ACF with definition of Tauzin 2019 -- all abs changes to abs2 -- calculate denominator only */
void autocorr_smooth_t2019_proto(int nfft, double deltat, double smooth, gsl_complex xautocorr[], double *smoothing)
{
  int i, j, nint, nelements, nfft_2;
  double df, xxtmp, xxtmp2;

  df = 1./(nfft*deltat);
  nint = floor(smooth/df/2);
  nfft_2 = nfft/2;

  /* get smoothing weights */ 
  xxtmp = 0; nelements = nint + 1;
  for(j=0; j <= nint; j++)
  {
    xxtmp += GSL_REAL(xautocorr[j]); //gsl_complex_abs2(xarr[j]);
  }
  xxtmp2 = xxtmp/nelements;
  smoothing[0] = xxtmp2;

  for(i=1; i < nint+1; i++)
  {
    xxtmp += GSL_REAL(xautocorr[nint+i]); // gsl_complex_abs2(xarr[nint+i]);
    nelements++;
    xxtmp2 = xxtmp/nelements;  
    smoothing[i] = xxtmp2;
  }
  
  nelements = 2*nint+1; /* constant for a while */
  for(i=nint+1; i< nfft_2-nint+1; i++)
  { 
    xxtmp += GSL_REAL(xautocorr[nint+i]); // gsl_complex_abs2(xarr[nint+i]);
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); // gsl_complex_abs2(xarr[i-nint-1]);
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2;
  }
  
  for(i = nfft_2-nint+1; i <= nfft_2; i++) /* decreasing */
  { 
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); //  gsl_complex_abs2(xarr[i-nint-1]);
    nelements--; 
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2;
  }
  
  /* negative frequency */
  xxtmp = 0; nelements = nint + 2;
  for(j= nfft_2; j <= nint+ nfft_2+1; j++)
  { 
    xxtmp += GSL_REAL(xautocorr[j]); // gsl_complex_abs2(xarr[j]);
  }
  xxtmp2 = xxtmp/nelements;
  smoothing[nfft_2+1] = xxtmp2;
  
  for(i= nfft_2+2; i < nfft_2+ nint; i++)
  { 
    xxtmp += GSL_REAL(xautocorr[nint+i]); // gsl_complex_abs2(xarr[nint+i]);
    nelements++; 
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2;
  }
  
  /* constant for a while */
  nelements = 2*nint +1 ; 
  for(i = nfft_2+ nint; i < nfft - nint; i++)
  { 
    xxtmp += GSL_REAL(xautocorr[i+nint]); // gsl_complex_abs2(xarr[i+nint]);
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); // gsl_complex_abs2(xarr[i-nint-1]);

    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2;
  }
  
  for(i = nfft-nint; i < nfft; i++) /* decreasing */
  { 
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); // gsl_complex_abs2(xarr[i-nint-1]);
    nelements--;
    xxtmp2 = xxtmp/nelements;
    smoothing[i] = xxtmp2;
  }

  return;
}

/* type 2 of the smoothing method */
void autocorr_smooth_t2019_proto_ver2(int nfft, double deltat, double smooth, gsl_complex xautocorr[], double *smoothing)
{
  int i, j, nint, nelements, nfft_2;
  double df, xxtmp; 

  df = 1./(nfft*deltat);
  nint = floor(smooth/df/2);
  nfft_2 = nfft/2;

  i = 0;
  nelements= nint + i + 1;
  xxtmp = 0;
  for(j=0; j <= i+nint; j++){ xxtmp += GSL_REAL(xautocorr[j]); } 
  smoothing[i] = xxtmp/nelements;

  for(i=1; i <= nint; i++)
  { 
    nelements++;
    xxtmp += GSL_REAL(xautocorr[i+nint]); 
    smoothing[i] = xxtmp/nelements;
  }

  nelements = 2*nint +1;
  for(i=nint+1; i <= nfft_2-nint; i++)
  {
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); 
    xxtmp += GSL_REAL(xautocorr[i+nint]); 
    smoothing[i] = xxtmp/nelements;
  }

  for(i = nfft_2-nint+1; i <= nfft_2; i++)
  {
    nelements--;
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); 
    smoothing[i] = xxtmp/nelements;
  }

  /* negative frequency */
  i = nfft_2 + 1;
  nelements = nint +i - nfft_2 +1;
  xxtmp = 0;
  for(j = nfft_2; j <= nint+i; j++)
  {
    xxtmp += GSL_REAL(xautocorr[j]); 
  }
  smoothing[i] = xxtmp/nelements;

  for(i= nfft_2+2; i < nfft_2+ nint+1; i++)
  {
    nelements++;
    xxtmp += GSL_REAL(xautocorr[nint+i]);
    smoothing[i] = xxtmp/nelements;
  }

  nelements = 2*nint + 1;
  for(i = nfft_2+ nint+1; i < nfft - nint; i++)
  {
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); 
    xxtmp += GSL_REAL(xautocorr[i+nint]); 
    smoothing[i] = xxtmp/nelements;
  }

  for(i = nfft-nint; i < nfft; i++)
  {
    nelements--;
    xxtmp -= GSL_REAL(xautocorr[i-nint-1]); 
    smoothing[i] = xxtmp/nelements;
  }

  return;
}

/* step 2 for the smoothed ACF -- and considering the noise. 
 * definition following Tauzin 2019 
*/
void autocorr_smooth_t2019_final_noise_ver2(int nfft, gsl_complex xautocorr[], double smoothing[], double square_sum, double snr_factor, gsl_complex *newautocorr)
{
  int i;
  double tmpdouble, add_factor = 0.;
  gsl_complex tmpcomp;

  if(snr_factor > 0){ add_factor= square_sum*snr_factor; } /* means we consider noise */
  for(i = 0; i < nfft; i++)
  {
    tmpcomp = gsl_complex_add_real(xautocorr[i], add_factor);  /* for the numerator */
    tmpdouble = smoothing[i] + add_factor; /* for the denominator */
    newautocorr[i] = gsl_complex_div_real(tmpcomp, tmpdouble);
  }

  return;
}

/* add noise considering possibility of existing differentiation */
void autocorr_smooth_t2019_final_noise_ver3(int nfft, gsl_complex xautocorr[], double smoothing[], double square_sum, double snr_factor, gsl_complex *newautocorr, bool diff_tf, double *wsquare)
{
  int i; 
  double tmpdouble, add_factor0, add_factor = 0.;
  gsl_complex tmpcomp;

  add_factor0 = 0;
  if(snr_factor > 0){ add_factor0 = square_sum*snr_factor; } /* means we consider noise */
  for(i = 0; i < nfft; i++)
  {
    add_factor = add_factor0;
    if(diff_tf){ add_factor = add_factor0*wsquare[i]; }
    
    tmpcomp = gsl_complex_add_real(xautocorr[i], add_factor);  /* for the numerator */
    tmpdouble = smoothing[i] + add_factor; /* for the denominator */
    newautocorr[i] = gsl_complex_div_real(tmpcomp, tmpdouble);
  }

  return;
}

/* with the whitening -- noise adding */
void autocorr_noise_byadding_ver3(int nfft, double square_sum, double snr_factor, gsl_complex *xautocorr, gsl_complex *newautocorr, bool diff_tf, double *wsquare)
{
  int i;
  double addfactor = 0.;

  if(snr_factor > 0.) /* if not, we are not using it --> noise free ! */
  {
    addfactor = square_sum*snr_factor;
    for(i=0; i < nfft; i++)
    {
      if(diff_tf){ newautocorr[i] = gsl_complex_add_real(xautocorr[i], addfactor*wsquare[i]); }
      else{ newautocorr[i] = gsl_complex_add_real(xautocorr[i], addfactor); }
    }
  }
  else
  { /* if no noise, return the original ! */
    for(i=0; i < nfft; i++)
    {
      newautocorr[i] = xautocorr[i];
    }
  }

  return;
}

void ratio_spectral(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio)
{
  int i;
  gsl_complex tmpconj, finalnumer, finaldenom;

  for(i=0; i < nfft; i++)
  {
    tmpconj = gsl_complex_conjugate(denom[i]);
    finalnumer = gsl_complex_mul(numer[i], tmpconj); /* H* Z' */
    finaldenom = gsl_complex_mul(denom[i], tmpconj); /* Z* Z' */
 
    ratio[i] = cdiv(finalnumer, finaldenom);
  }

  return;
}

void ratio_simplediv(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio)
{ /* well, the ratio will be double (bc double divided by double), but for the sake of uniform format, we do complex */
  int i;
  double tmpcomp, tmpnumer, tmpdenom;

  for(i=0; i < nfft; i++)
  {
    tmpnumer = gsl_complex_abs(numer[i]);
    tmpdenom = gsl_complex_abs(denom[i]);
    tmpcomp = tmpnumer/tmpdenom;
    GSL_SET_COMPLEX(ratio+i, tmpcomp, 0);
  }

  return;
}

/* the hanning window with a given length */
void hanning(int npts, float *window)
{
  int i, tmpnpts;
  float x = 0;

  tmpnpts = npts -1;
  for(i=0; i < npts; i++)
  {
    x = (float) i/tmpnpts;
    *(window+i) = 0.5*(1-cos(M_PI*2*x)); 
    // fprintf(stderr, "%d: %.5f\n", i, window[i]);
  }
 
  return;
}

/* calculate area of the function */
float area_func(int npts, float func[])
{
  float area = 0;
  int i = 0;
  for(i = 0; i < npts; i++)
  {
    area += func[i];
  }

  return area;
}

/* smoothing using certain time window */
/* The length of the function is 2*nint_half+1 and nint_half = floor(smooth/def/2) where df = 1/(nfft*deltat) */
void autocorr_smooth_func(int nfft, int nint_half, float *window, gsl_complex *xarr)
{
  int i, j;
  gsl_complex xtmp, xconj;
  double xxtmp, *smoothing;
  float area = 0;

  smoothing = calloc(nfft, sizeof(double)); 
  for(i=0; i < nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(nint_half+i+1, window + nint_half -i);
    //fprintf(stderr, "area: %.4f\n", area);
    for(j=0; j <=i+nint_half; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[nint_half -i+j] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i=nint_half; i< nfft/2-nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(2*nint_half +1, window);
    for(j = i-nint_half; j <= i+nint_half; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j-i+nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft/2-nint_half; i <= nfft/2; i++)
  {
    xxtmp = 0;
    area = area_func(nint_half + (nfft/2-i)+1, window + nint_half - (nfft/2-i));
    for(j=i - nint_half; j <= nfft/2; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j-i + nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i= nfft/2+1; i < nfft/2+ nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(nint_half +i - nfft/2 +1, window + nint_half - i + nfft/2);
    for(j = nfft/2; j <= nint_half+i; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[ j - i + nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft/2+ nint_half; i < nfft - nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(2*nint_half + 1, window);
    for(j= i - nint_half; j <= i+nint_half; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j - i + nint_half ] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft-nint_half; i < nfft; i++)
  {
    xxtmp = 0;
    area = area_func(nfft -i + nint_half, window + nint_half +1 + i - nfft);
    for(j = i - nint_half; j < nfft; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j - i + nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }

  for(i=0; i < nfft; i++)
  {
    xconj = gsl_complex_conjugate(xarr[i]);
    xtmp = gsl_complex_mul(xconj, xarr[i]);
    /* now smoothing!! */
    xarr[i] = adiv(smoothing[i], xtmp);
  }

  free(smoothing);
  return;
}

