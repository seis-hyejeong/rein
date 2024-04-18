#include "func_synthetic_iso3.h"

/* pre state the working space... */
gsl_complex WScomp[maxnfft];
double WSdouble0[maxnfft], WSdouble[maxnfft*2];
float WSfloat[maxnfft];

/* stack the data */
void calculate_stack(int ntrace, int npts, float xindv[], float *xstack)
{
  int i, j;
  double inverse = 1./ntrace;

  for(i = 0; i < npts; i++)
  {
    xstack[i] = xindv[i]; /* initialize */
  }

  for(i=1; i < ntrace; i++)
  {
    for(j=0; j < npts; j++)
    {
      xstack[j] += xindv[ i*npts + j];
    }
  }

  for(j = 0; j < npts; j++)
  {
    xstack[j] *= inverse;
  }

  return;
}

/* calculate average and the stddev of the data */
void calculate_stack_stddev(int ntrace, int npts, float xindv[], float *xstack, float *stddev)
{
  int i, j;
  double inverse = 1./ntrace;

  for(i = 0; i < npts; i++)
  {
    xstack[i] = xindv[i]; /* initialize */
    stddev[i] = 0.;
  }

  for(i=1; i < ntrace; i++)
  {
    for(j=0; j < npts; j++)
    {
      xstack[j] += xindv[ i*npts + j];
    }
  }

  for(j = 0; j < npts; j++)
  {
    xstack[j] *= inverse; /* get the average */
  }

  /* calculate stddev */
  inverse = 1./(ntrace - 1);
  for(j=0; j < npts; j++) /* calculate standard deviation from in data */
  {
    *(stddev+j) =0.;
    for(i=0; i < ntrace; i++)
    {
      *(stddev+j) += (xindv[i*npts+j]- *(xstack+j))*(xindv[i*npts+j] - *(xstack+j));
    }
    *(stddev+j) = (float) sqrt(*(stddev+j)*inverse);
  }

  return;
}

/* calculate covariance of the data -- cholesky decomposition 
 We use LDLT decomposition for practical reason (semi positive definite case)
 the matrix L must be the inverse of the decomposed L.
 The decomposed L is lower triangular. The inverse of it is also lower triangular.
*/
int calculate_covariance(int ntrace, int npts, float xindv[], float xstack[], gsl_matrix *cov, gsl_matrix *L, gsl_permutation *permute)
{
  size_t ii, jj;
  int s;
  double tmpdouble;
  gsl_matrix *diff;

  /* allocate the matrix space */
  diff = gsl_matrix_calloc(npts, ntrace);
  if(diff ==NULL){ fprintf(stderr, "error in allocation!!!!!!!!!\n"); return -1; }

  /* now copy the read data into the matrix calculating the deviation */
  for(ii=0; ii < ntrace; ii++)
  {
    for(jj=0; jj < npts; jj++)
    {
      tmpdouble = xindv[npts*ii + jj] - xstack[jj];
      gsl_matrix_set(diff, jj, ii, tmpdouble);
    }
  }

  /* calculate the covariance by matrix multiple */
  tmpdouble = 1./ntrace;
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, tmpdouble, diff, diff, 0., cov);

  /* copy the matrix to get LDLT decomp */
  gsl_matrix_memcpy(L, cov);

  gsl_matrix_free(diff);
  gsl_linalg_LU_decomp(L, permute, &s);
  return 0; 
}

/* stack the function */
void stack_synths(int ntraces, int npts, float *xsynth)
{ /* we overwrite the array with the synthetic stack */

  int i, j;
  double inverse;

  for(i=1; i < ntraces; i++)
  {
    for(j=0; j < npts; j++)
    {
      xsynth[j] += xsynth[i*npts + j];
    }
  }

  inverse = 1./ntraces;
  if(ntraces> 1)
  {
    for(i=0; i < npts; i++){ xsynth[i] *= inverse; }
  }

  return;
}

/* if you decided to make output of stacked synthetics, let me divide it here.*/
void average_synths(int ntraces, int npts, float *stacked_xsynth)
{
  int i;
  double inverse = 0;
  inverse = 1./ntraces;

  for(i = 0; i < npts; i++){ *(stacked_xsynth+i) *= inverse; }

  return;
}

/* we can multiply to the trace to avoid repeating divisions */
void Ntimes_trace(int ntraces, int npts, float *stacked_obs)
{ 
  int i; 

  for(i = 0; i < npts; i++){ *(stacked_obs+i) *= ntraces; }

  return;
}

/* get synthetic error -- in frequency 
WS: working space of length nfft*2
xerr: gsl complex trace of the synthetic error -- filtered if gaussf > 0
*/
int synthetic_noises_freq(gsl_rng *r, int nerr, int nfft, double sigma, gsl_complex *xerr, double deltat, double gaussf, char savename[strlng])
{
  int i, j;
  char outname[strlng+10];
  FILE *filestr; 
  if(nerr < 1){ return -1; }

  for(i=0; i < nerr; i++)
  {
    for(j=0; j < nfft; j++)
    {
      WSdouble[j] = gsl_ran_gaussian(r, sigma);
    }

    /* forward fourier transform */
    trace_fft(WSdouble, nfft, xerr+i*nfft);
    if(gaussf > 0){ gaussian_lp(gaussf, nfft, deltat, xerr+i*nfft); trace_inv_fft(xerr+i*nfft, nfft, WSdouble); }

    /* work for writing */
    sprintf(outname, "%s.%03d.txt", savename, i);
    filestr = fopen(outname, "w");
    for(j=0; j < nfft; j++)
    {
      fprintf(filestr, "%.5f %.8f\n", j*deltat, WSdouble[j]);
    } 
    fclose(filestr);
  }

  return 0;
}

/* apply gaussian low pass to the noise */
double gaussian_lp(double awidth, int nfft, double deltat, gsl_complex *xarr1)
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

  normalizer = nfft/normalizer;

  /* invert the normalizer normalizer = nfft/normalizer; */
  for(i=0; i < nfft; i++)
  {
    xarr1[i] = amul(normalizer, xarr1[i]);
  }

  /* return the normalizing value */
  normalizer = 1./normalizer;

  return normalizer;
}

/* for the error consideration version 2: we have to get sum of the frequency trace */
double sigma_complex_abs(int npts, gsl_complex *xtrace_c)
{
  int i;
  double sum = 0.;
  for(i=0; i < npts; i++)
  {
    sum += gsl_complex_abs(xtrace_c[i]);
  }

  return sum;
}

double sigma_complex_abs2(int npts, gsl_complex *xtrace_c)
{
  int i;
  double sum = 0.;
  for(i=0; i < npts; i++)
  {
    sum += gsl_complex_abs2(xtrace_c[i]);
  }

  return sum;
}

/* add a constant double throughout the array */
void add_complex_double_series(int npts, gsl_complex *xtrace_c, double add)
{
  int i;
  
  for(i=0; i < npts; i++)
  {
    xtrace_c[i].dat[0] += add;
  }

  return;
}

/* calculate synthetics as we want */
int calculate_synthetics_noise_ver6(PARAM inparam, float *s_racf, float *indv_racf, INDATA racf, float *s_zacf, float *indv_zacf, INDATA zacf, float *s_pacf, float *indv_pacf, INDATA pacf, float *s_rstack, float *s_zstack, INDATA rstack, MODEL m, SYNTHPAR synth, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw, gsl_complex *shift, double *w2array)
{
  int errn;
  double mchris[81*inparam.wnl];

  /* calculate only if it's true */
  set_christoffel_iso_water( inparam.wnl, m.h, m.rho, m.vp, m.vs, mchris);
  
  errn = 0;
  if(racf.tf)
  {
    errn = make_synthetic_ACF_R_iso_noise_ver6(inputR, racf, inparam.wnl, m.h, m.rho, mchris, inparam.iwater, m.vp[0], synth.nfft, synth.dw, s_racf, indv_racf, synth.racfind[0], tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, expw, shift, w2array);
    if(errn !=0){ return 2; }
    normalize_single(synth.nfft, s_racf, synth.racfind0); 
  }
  if(zacf.tf)
  {
    errn = make_synthetic_ACF_Z_iso_noise_ver6(inputZ, zacf, inparam.wnl, m.h, m.rho, mchris, inparam.iwater, m.vp[0], synth.nfft, synth.dw, s_zacf, indv_zacf, synth.zacfind[0], tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, expw, shift, w2array);
    if(errn != 0 ){ return 3; }
    normalize_single(synth.nfft, s_zacf, synth.zacfind0);
  }
  if(pacf.tf) /* the dpg version */
  {
    errn = make_synthetic_ACF_DPG_iso_noise_ver6(inputZ, pacf, inparam.wnl, m.h, m.rho, mchris, inparam.iwater, m.vp[0], synth.nfft, synth.dw, s_pacf, indv_pacf, synth.pacfind[0], tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, expw, shift, w2array);
    if( errn != 0 ){ return 4; }
    normalize_single(synth.nfft, s_pacf, synth.pacfind0);
  }  
  if(rstack.tf)
  {
    errn = make_synthetic_3comp_iso_ver2(inputZ, rstack.nfiles, rstack.rayp, inparam.wnl, m.h, m.rho, mchris, inparam.iwater, m.vp[0], synth.nfft, synth.dw, rstack.delta, synth.rwaveind0*rstack.delta, rstack.gauss, s_rstack, s_zstack, 1, tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, expw);
    /* errn = make_synthetic_Rcomp_iso(inputZ, rstack.nfiles, rstack.bazs, rstack.rayp, inparam.wnl, m.h, m.rho, mchris, 1, m.vp[0], synth.nfft, synth.dw, rstack.delta, synth.rwaveind0*rstack.delta, rstack.gauss, s_rstack, 1, tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, expw); */
    if(errn !=0){ return 1; }
    /* filter additionally, if there is */
    filter_signal(synth.nfft, rstack.delta, rstack.ffilter[0], rstack.ffilter[1], s_rstack);
    filter_signal(synth.nfft, rstack.delta, rstack.ffilter[0], rstack.ffilter[1], s_zstack);
    normalize_RZ_ver2(synth.nfft, s_zstack, s_rstack, synth.rwaveind0);
    normalize_max_2_new(synth.nfft, s_rstack, synth.rwaveind[0], rstack.w, synth.rwaveind[1]- synth.rwaveind[0]);
  }

  /* now return !!*/
  return 0;
}

/* get christoffel matrix with the model */
void set_christoffel_iso_water(int nlayer, float mod_h[], float mod_rho[], float mod_p[], float mod_s[], double *chris)
{
  int i;

  for(i=0; i < nlayer; i++)
  {
      set_elastic_isotropic(mod_rho[i], mod_p[i], mod_s[i], chris+81*i);
  }

  return;
}

/* single layer christoffel -- anisotropic */
void set_elastic_anisotropic(float rho, float alpha, float beta, float B, float C, float E, float azi, float tilt, double *o_tensor)
{
  double c1111, c2222, c3333, c1122, c1133, c1313, c1212;
  double A, D, rotM[9], tensor[81];
  int i, j, k, l;
  double tmp2;

  tensor_set_all(tensor, 0); // initialize with value 0.

  A=0;
  D=0;

  c1111 = (1+A-B+C)*rho*alpha*alpha;
  c2222 = c1111;
  c3333 = (1+A+B+C)*rho*alpha*alpha;
  c1122 = (1+A-B+C)*rho*alpha*alpha - 2*(1+D-E)*rho*beta*beta;
  c1133 = (1+A-3*C)*rho*alpha*alpha - 2*(1+D+E)*rho*beta*beta;
  c1313 = (1+A+B+C)*rho*beta*beta;
  c1212 = (c1111-c1122)/2;

  tensor_set(tensor, 1,1,1,1, c1111);
  tensor_set(tensor, 2,2,2,2, c2222);
  tensor_set(tensor, 3,3,3,3, c3333);
  tensor_set(tensor, 1,1,2,2, c1122);
  tensor_set(tensor, 1,1,3,3, c1133);
  tensor_set(tensor, 2,2,3,3, c1133);
  tensor_set(tensor, 1,3,1,3, c1313);
  tensor_set(tensor, 2,3,2,3, c1313);
  tensor_set(tensor, 1,2,1,2, c1212);

  for(i=1; i<=3; i++)
  {
    for(j=1; j <=3; j++)
    {
      for(k=1; k <=3; k++)
      {
        for(l=1; l <=3; l++)
        {
          if( tensor_get(tensor, i, j, k, l)!=0)
          {
            tmp2 = tensor_get(tensor, i, j, k, l);
            tensor_set(tensor, i, j, l, k, tmp2);
            tensor_set(tensor, j, i, k, l, tmp2);
            tensor_set(tensor, j, i, l, k, tmp2);
            tensor_set(tensor, k, l, i, j, tmp2);
          }
        }
      }
    }
  }
  get_aaxis_rotmat(azi, tilt, rotM);
  rot_tensor(rotM, tensor, o_tensor);

  return;
}

/* single layer christoffel -- isotropic */
void set_elastic_isotropic(float rho, float alpha, float beta, double *tensor)
{
  double c1111, c2222, c3333, c1122, c1133, c1313, c1212;
  double A, D;
  int i, j, k, l;
  double tmp2;

  tensor_set_all(tensor, 0); // initialize with value 0.

  A=0;
  D=0;

  c1111 = rho*alpha*alpha;
  c2222 = c1111;
  c3333 = (1)*rho*alpha*alpha;
  c1122 = (1)*rho*alpha*alpha - 2*(1)*rho*beta*beta;
  c1133 = (1)*rho*alpha*alpha - 2*(1)*rho*beta*beta;
  c1313 = (1)*rho*beta*beta;
  c1212 = (c1111-c1122)/2;

  tensor_set(tensor, 1,1,1,1, c1111);
  tensor_set(tensor, 2,2,2,2, c2222);
  tensor_set(tensor, 3,3,3,3, c3333);
  tensor_set(tensor, 1,1,2,2, c1122);
  tensor_set(tensor, 1,1,3,3, c1133);
  tensor_set(tensor, 2,2,3,3, c1133);
  tensor_set(tensor, 1,3,1,3, c1313);
  tensor_set(tensor, 2,3,2,3, c1313);
  tensor_set(tensor, 1,2,1,2, c1212);

  for(i=1; i<=3; i++)
  {
    for(j=1; j <=3; j++)
    {
      for(k=1; k <=3; k++)
      {
        for(l=1; l <=3; l++)
        {
          if( tensor_get(tensor, i, j, k, l)!=0)
          {
            tmp2 = tensor_get(tensor, i, j, k, l);
            tensor_set(tensor, i, j, l, k, tmp2);
            tensor_set(tensor, j, i, k, l, tmp2);
            tensor_set(tensor, j, i, l, k, tmp2);
            tensor_set(tensor, k, l, i, j, tmp2);
          }
        }
      }
    }
  }

  return;
}

/* get relaxed number of nfft */
int get_nfft(double time, double dt)
{
  int nfft = 1;
  while( dt*nfft < time )
  {
    nfft *= 2;
  } 
  
  return nfft;
}

/* get index of zero time of synthetics and the time information */
void set_parameters(int nfft, double dt, double *t0_rf, double *t0_acf, double *dw)
{
  double tlength = nfft*dt;
  *(t0_rf) = tlength*0.2; /* put 20 % of length to the pre-signal part */
  *(t0_acf) = 0.5*tlength;
  *(dw) = M_PI*2/tlength;

  return;
}

void get_zerotindex(double t0_rf, double t0_acf, double dt, int *ind0_rf, int *ind0_acf)
{
  *(ind0_rf) = round(t0_rf/dt); 
  *(ind0_acf) = round(t0_acf/dt);
  
  return;
}

double convert_gauss(double gauss)
{
  double convertgauss = 1;
  /* the product of two parameters should be 0.5 */

  convertgauss = 0.5/gauss;

  return convertgauss;
}

/* sac filtering the synthetic traces: using subroutine provided by sac */
int filter_signal(int npts, double delta, double low, double high, float *xtrace)
{
  if(high < 0 && low >0)
  {
    xapiir(xtrace, npts, SAC_BUTTERWORTH, bandwidth, attenuation, order, SAC_HIGHPASS, low, high, delta, pass );
  }
  else if(low < 0 && high >0)
  {
    xapiir(xtrace, npts, SAC_BUTTERWORTH, bandwidth, attenuation, order, SAC_LOWPASS, low, high, delta, pass );
  }
  else if(high > 0 && low > 0)
  {
    xapiir(xtrace, npts, SAC_BUTTERWORTH, bandwidth, attenuation, order, SAC_BANDPASS, low, high, delta, pass );
  }
  else{ return -1; }

  return 0;
}

/* get max */
double get_abs_max(int npts, double *trace)
{
  double max = 0;
  int i = 0;
  while(i < npts)
  {
    if(fabs(trace[i]) > max){ max = fabs(trace[i]); }
    i++;
  }

  return max;
}

/* normalize with maximum in the time window of the trace */
float normalize_max(int npts, float *trace, int ind0, int windownpts)
{
  float max = 0;
  int j, i = ind0;
  j = ind0 + windownpts;
  while(i < j)
  {
    if(trace[i] > max){ max = trace[i]; }
    i++;
  }

  if(max < 1.E-6){ return -1; } /* could not get a larger positive value */

  max = 1./max;
  for(i=0; i < npts; i++)
  {
    trace[i] *= max;
  }

  return max;
}

/* normalize to a new number you want -- to downscale the R WAVE data */
float normalize_max_2_new(int npts, float *trace, int ind0, double newmax, int windownpts)
{
  float max = 0;
  int j, i = ind0;
  j = ind0 + windownpts;
  while(i < j)
  {
    if(trace[i] > max){ max = trace[i]; }
    i++;
  }

  if(max < 1.E-6){ return -1; } /* could not get a larger positive value */

  max = newmax/max;
  for(i=0; i < npts; i++)
  {
    trace[i] *= max;
  }

  return max;
}

/* after filter, normalize -- relative to Z at certain time*/
int normalize_RZ(int npts, float *ztrace, float *rtrace, int ind0)
{
  int i;
  double normalizer = 1;
  if(fabs(ztrace[ind0]) < 1.E-4){ fprintf(stderr, "problem in normalizing 2 traces\n"); return -1; }
 
  normalizer = 1./ztrace[ind0];

  for(i=0; i < npts; i++)
  { 
    *(ztrace+i) = ztrace[i]*normalizer;
    *(rtrace+i) = rtrace[i]*normalizer;
  }

  return 0;
}

/* after filter, normalize -- relative to Z at maximum time. and shifting needed. */
int normalize_RZ_ver2(int npts, float *ztrace, float *rtrace, int ind0)
{
  int idiff, i, maxind;
  double normalizer = -1;

  /* but note that it is good to use high frequency */
  maxind = 0;
  for(i = 0; i < npts; i++)
  {
    if(ztrace[i] > normalizer){ normalizer = ztrace[i]; maxind = i; } 
  }

  /* make maxind to locate at ind0 */
//  if( fabs(maxind - ind0) > 1 ){  fprintf(stderr, "there's time shift of the normalizing waveforms by Z\n"); }
  if( fabs(maxind - ind0) > npts*0.5){ fprintf(stderr, "problem in normalizing 2 traces (2)\n"); return -1; } 
  /* too much of shifting, this means something not right */
  if(fabs(ztrace[maxind]) < 1.E-6){ fprintf(stderr, "problem in normalizing 2 traces\n"); return -1; }

  normalizer = 1./ztrace[maxind];
  idiff = maxind - ind0;

  for(i=0; i < npts; i++)
  {
    *(ztrace+i) = ztrace[i]*normalizer;
    *(rtrace+i) = rtrace[i]*normalizer;
  }
  
  if( idiff == 0){ }
  else if( idiff > 0 ) /* when the arrival is apparently later */ 
  { /*  I should make the traces come earlier */    
    for(i = 0; i < npts - idiff; i++)
    {
      ztrace[i] = ztrace[i + idiff];
      rtrace[i] = rtrace[i + idiff];
    }
    for(i = npts - idiff; i < npts; i++)
    {
      ztrace[i] = 0.;
      rtrace[i] = 0.;
    } 
  }
  else if( idiff < 0 ) /* when the arrival is apparently earlier */
  {
    for(i = npts -1; i >= -idiff; i--)
    {
      ztrace[i] = ztrace[i + idiff];
      rtrace[i] = rtrace[i + idiff]; 
    }
    for(i = 0; i < -idiff; i++)
    {
      ztrace[i] = 0.;
      rtrace[i] = 0.;
    } 
  }
  else
  {
    fprintf(stderr, "the algorithm of shifting time is wrong\n"); return -1; 
  }

  return 0;
}

/* after filter, normalize -- relative to amplitude at certain time */
int normalize_single(int npts, float *trace, int ind0)
{
  int i;
  double normalizer = 1;
  if(fabs(trace[ind0]) < 1.E-4){ fprintf(stderr, "problem in normalizing 2 traces\n"); return -1; }

  normalizer = 1./trace[ind0];

  for(i=0; i < npts; i++)
  { 
    *(trace+i) = trace[i]*normalizer;
  }

  return 0;
}

/* calculate synthetics of waveform */
int make_synthetic_3comp_iso_ver2(int iinput, int nray, double i_rayp[maxray], int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, double dt, double t0, double gausst, float *xtrace_r, float *xtrace_z, int inorm, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw)
{
  int errn, i, j, k, ii, jj, ind0;
  double t00, t_first=0;
  double q3sol[6*nlayer], dvec[3*6], q0tmp[3],  *dxfft0;
  double norm = 0;
  gsl_complex tmpcomp, xfft0[nfft], x1_response[ nfft], x3_response[ nfft], x1_response_stack[ nfft], x3_response_stack[nfft], surftmp[6], surftmp2[6];
/*  gsl_complex x2_response[ nfft], x2_response_stack[nfft]; */
  int nreceiverlayer = 1;
  errn = 0;
  ind0 = (int) t0/dt;
  /* make input waveform */
  dxfft0 = (double *)calloc(sizeof(double *), nfft);

  /* if gausst > 0 */ 
  if(gausst > 0)
  {
    for(i=0; i < nfft; i++)
    {
      *(dxfft0+i) = exp( (-1.)*pow(dt*(i-nfft/2), 2)/(4*gausst*gausst));
      x1_response_stack[i] = GSL_COMPLEX_ZERO;
      x3_response_stack[i] = GSL_COMPLEX_ZERO;
    }
    trace_fft(dxfft0, nfft, xfft0);
    trace_shift_amp(-(nfft/2)*dt, dt, nfft, 1., xfft0); /* make it as if that peak is at t=0 */

//    if(iverbose==1){ fprintf(stderr, "%d number of points for fft and trace base is made.\n", nfft); }
  }
  else
  {
    for(i=0; i < nfft; i++)
    {
      GSL_SET_COMPLEX(xfft0+i, 1., 0.);
      x1_response_stack[i] = GSL_COMPLEX_ZERO; 
//      x2_response_stack[i] = GSL_COMPLEX_ZERO;
      x3_response_stack[i] = GSL_COMPLEX_ZERO;
    }
  }

  /* onto main calculation */
  for(i=0; i < nray; i++)
  {
    errn=0;
    /* set horizontal slownesses (into 2D) */
    *(q0tmp+0) = *(i_rayp+i);
    *(q0tmp+1) = 0.;

    t_first=0.;

    /* copy xfft0 into three components */
    trace_copy(nfft, xfft0, x1_response);
    trace_copy(nfft, xfft0, x3_response);

    /* first calculate the ones that are independent to frequency: except for the water layer */
    if(iwater==1)
    { 
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;} 
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); 
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha; 
     
      /* get the Ematrix 6 x 6 */
      E_matrix(m_chris, *(q0tmp+0), q3sol, dvec, tmpE[0]);
      /* extract to 2 x 2 matrix and get Einv */
      gsl_matrix_complex_set(tmpE_w, 0, 0, gsl_matrix_complex_get(tmpE[0], 2, 0));
      gsl_matrix_complex_set(tmpE_w, 1, 0, gsl_matrix_complex_get(tmpE[0], 5, 0)); /* for n = 1 */
      gsl_matrix_complex_set(tmpE_w, 0, 1, gsl_matrix_complex_get(tmpE[0], 2, 3));
      gsl_matrix_complex_set(tmpE_w, 1, 1, gsl_matrix_complex_get(tmpE[0], 5, 3)); /* for n = 4 */
      Einv_w_matrix(tmpE_w, tmpinvE_w);
    } 

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2_iso(*(q0tmp+0), m_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ // fprintf(stderr, "complex root coming out!!! Rwave (3)\n"); 
//	fprintf(stderr, "synth-3 %d layer density %f ray parameter %f\n", j, m_rho[j], i_rayp[i]);
        return -1; // goto M_00001;
      }
      errn = cal_oscillvec_iso(m_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); return -2; }
      if(j!=nlayer-1 && j >= nreceiverlayer){      t_first += fabs(*(q3sol+6*j+iinput-1)**(m_thick+j)); } 
      E_matrix(m_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE[j]);
      Einv_matrix(tmpE[j], tmpinvE[j]);
    }
    
    /* calculate expw */
    cal_expw(nlayer-1, q3sol, m_thick, dw, nfft, expw);
 
    /* get matrix E for n-th layer (j=nlayer-1) */
    /* now frequency dependent part.. 
     * calculations will be quite large due to calculation with respect to
     * frequency vale... */
    for(k=0; k < nfft; k++)
    {
      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE[j], tmpinvE[j], expw + 6*nfft*j + 6*k, A[j], mat66); 
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE[0], tmpinvE_w, expw + 6*k, A_w, mat22); /* j = 0 */
      }
 
      /* now calculate matrix J */
      Jw_matrix(nlayer-iwater,  tmpinvE[nlayer-1], &A[iwater], Jmatrix, mat66);
      
      /* response at the top of the solid layer */
      if(iwater==0){  get_response(iinput, Jmatrix, surftmp, J11, x, in, p); }
      else /* if iwater ==1 */
      {
        /* Jmatrix_w is a matrix to solve (u11, u12, u13, tau33, u03) */
        /* u1i are displacement at seafloor, tau33 is stress at seafloor, f(i) are coeff at the bottom of the model, u03 is at surface of water */
        /* first, initialize the Jmatrix_w */
        gsl_matrix_complex_set_zero(Jmatrix_w);
 
        /* second, make the matrix for iwater: first */
        for(ii=0; ii < 3; ii++) 
        {
          for(jj=0; jj<3; jj++){ gsl_matrix_complex_set(Jmatrix_w, ii, jj, gsl_matrix_complex_get(Jmatrix, ii, jj)); }
          gsl_matrix_complex_set(Jmatrix_w, ii, 3, gsl_matrix_complex_get(Jmatrix, ii, 5));
        }
      
        GSL_SET_COMPLEX(&tmpcomp,-1, 0);
        gsl_matrix_complex_set(Jmatrix_w, 3, 2, tmpcomp); /* at (4, 3) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 3, tmpcomp); /* at (5, 4) */
        gsl_matrix_complex_set(Jmatrix_w, 3, 4, gsl_matrix_complex_get(A_w, 0,0)); /* put A11 to (4, 5) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 4, gsl_matrix_complex_get(A_w, 1,0)); /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp, J11, x, in, p);
      }

      /* now get response at the receiver if the receiver is located deeper than the layer */
      /* first, we need to get Jw2 */
      if(iwater < nreceiverlayer)
      { 
        Jw_matrix2(nreceiverlayer-iwater, &A[iwater], Jmatrix, mat66);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2); 
        /* copy back to the surftmp */
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }

      x1_response[k] = gsl_complex_mul(*(x1_response+k), surftmp[0]);
      x3_response[k] = gsl_complex_mul(*(x3_response+k), surftmp[2]);
    }

    t00 = t0 - t_first;
    trace_shift_amp(t00, dt, nfft, 1, x1_response);
    trace_shift_amp(t00, dt, nfft, 1, x3_response);

    if(inorm ==1) /* this must be well done!!! */
    {
      trace_inv_fft(x3_response, nfft, dxfft0);
      norm = -1;

      for(k=0; k < nfft; k++)
      {
        if( norm < -dxfft0[k] ){ norm = -dxfft0[k]; } /* originally: norm = -1./dxfft0[ind0]; */
      } /* our assumption: the timing of the peak will be same for a same structure but different ray parameter */
      norm = 1./norm;
      for(k=0; k < nfft; k++)
      {
        x1_response[k] = gsl_complex_mul_real(x1_response[k], norm);
        x3_response[k] = gsl_complex_mul_real(x3_response[k], norm);
      }
    }

    /* add it to the stacked frequency trace */
    for(k=0; k < nfft; k++)
    {
      *(x1_response_stack+k) = gsl_complex_add(x1_response[k], *(x1_response_stack+k));
      *(x3_response_stack+k) = gsl_complex_add(x3_response[k], *(x3_response_stack+k));
    }
    
  }
 
  /* inverse fourier transform */
  trace_inv_fft(x1_response_stack, nfft, dxfft0);  
  for(k=0; k< nfft; k++){ *(xtrace_r+k) = (float) *(dxfft0+k); } 
  trace_inv_fft(x3_response_stack, nfft, dxfft0);
  for(k=0; k< nfft; k++){ *(xtrace_z+k) = (float) (-1)**(dxfft0+k); }
 
  free(dxfft0);
  return 0;
}

/* make only R. Z is only for normalizing */
int make_synthetic_Rcomp_iso(int iinput, int nray, double i_rayp[maxray], int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, double dt, double t0, double gausst, float *xtrace_r, int inorm, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw)
{
  int errn, i, j, k, ii, jj, ind0;
  double t00, t_first=0;
  double q3sol[6*nlayer], dvec[3*6], q0tmp[3],  *dxfft0;
  double norm = 0;
  gsl_complex tmpcomp, xfft0[nfft], x1_response[ nfft], x3_response[ nfft], x1_response_stack[nfft], surftmp[6], surftmp2[6];
  int nreceiverlayer = iwater;

  errn = 0;
  ind0 = (int) t0/dt;
  /* make input waveform */
  dxfft0 = (double *)calloc(sizeof(double *), nfft);

  /* if gausst > 0 */
  if(gausst > 0)
  {
    for(i=0; i < nfft; i++)
    {
      *(dxfft0+i) = exp( (-1.)*pow(dt*(i-nfft/2), 2)/(4*gausst*gausst));
      x1_response_stack[i] = GSL_COMPLEX_ZERO;
    }
    trace_fft(dxfft0, nfft, xfft0);
    trace_shift_amp(-(nfft/2)*dt, dt, nfft, 1., xfft0); /* make it as if that peak is at t=0 */

    //if(iverbose==1){ fprintf(stderr, "%d number of points for fft and trace base is made.\n", nfft); }
  }
  else
  {
    for(i=0; i < nfft; i++)
    {
      GSL_SET_COMPLEX(xfft0+i, 1., 0.);
      x1_response_stack[i] = GSL_COMPLEX_ZERO; 
    }
  }

  /* onto main calculation */
  for(i=0; i < nray; i++)
  {
    errn=0;
    /* set horizontal slownesses (into 2D) */
    *(q0tmp+0) = *(i_rayp+i);
    *(q0tmp+1) = 0.;

    t_first=0.;

    /* copy xfft0 into three components */
    trace_copy(nfft, xfft0, x1_response);
    trace_copy(nfft, xfft0, x3_response);

    /* first calculate the ones that are independent to frequency: except for the water layer */
    if(iwater==1)
    {
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;}
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] );
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha;

      /* get the Ematrix 6 x 6 */
      E_matrix(m_chris, *(q0tmp+0), q3sol, dvec, tmpE[0]);
      /* extract to 2 x 2 matrix and get Einv */
      gsl_matrix_complex_set(tmpE_w, 0, 0, gsl_matrix_complex_get(tmpE[0], 2, 0));
      gsl_matrix_complex_set(tmpE_w, 1, 0, gsl_matrix_complex_get(tmpE[0], 5, 0)); /* for n = 1 */
      gsl_matrix_complex_set(tmpE_w, 0, 1, gsl_matrix_complex_get(tmpE[0], 2, 3));
      gsl_matrix_complex_set(tmpE_w, 1, 1, gsl_matrix_complex_get(tmpE[0], 5, 3)); /* for n = 4 */
      Einv_w_matrix(tmpE_w, tmpinvE_w);
    }

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2_iso(*(q0tmp+0), m_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ // fprintf(stderr, "complex root coming out!!! Rwave (1)\n"); 
//        fprintf(stderr, "synth-R %d layer density %f ray parameter %f\n", j, m_rho[j], i_rayp[i]);
        return -1;
        goto M_00002;
      }
      errn = cal_oscillvec_iso(m_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); return -2; }
      if(j!=nlayer-1 && j >= nreceiverlayer){      t_first += fabs(*(q3sol+6*j+iinput-1)**(m_thick+j)); } 
      E_matrix(m_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE[j]);
      Einv_matrix(tmpE[j], tmpinvE[j]);
    }

    /* calculate expw */
    cal_expw(nlayer-1, q3sol, m_thick, dw, nfft, expw);

    for(k=0; k < nfft; k++)
    {
      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE[j], tmpinvE[j], expw + 6*nfft*j + 6*k, A[j], mat66);
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE[0], tmpinvE_w, expw + 6*k, A_w, mat22); /* j = 0 */
      }

      Jw_matrix(nlayer-iwater, tmpinvE[nlayer-1], &A[iwater], Jmatrix, mat66);

      /* response at the top of the solid layer */
      if(iwater==0){  get_response(iinput, Jmatrix, surftmp, J11, x, in, p); }
      else /* if iwater ==1 */
      {
        gsl_matrix_complex_set_zero(Jmatrix_w); 

        /* second, make the matrix for iwater: first */
        for(ii=0; ii < 3; ii++)
        {
          for(jj=0; jj<3; jj++){ gsl_matrix_complex_set(Jmatrix_w, ii, jj, gsl_matrix_complex_get(Jmatrix, ii, jj)); }
          gsl_matrix_complex_set(Jmatrix_w, ii, 3, gsl_matrix_complex_get(Jmatrix, ii, 5));
        }

        GSL_SET_COMPLEX(&tmpcomp,-1, 0);
        gsl_matrix_complex_set(Jmatrix_w, 3, 2, tmpcomp); /* at (4, 3) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 3, tmpcomp); /* at (5, 4) */
        gsl_matrix_complex_set(Jmatrix_w, 3, 4, gsl_matrix_complex_get(A_w, 0,0)); /* put A11 to (4, 5) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 4, gsl_matrix_complex_get(A_w, 1,0)); /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp, J11, x, in, p);
      }

      /* now get response at the receiver if the receiver is located deeper than the layer */
      /* first, we need to get Jw2 */
      if(iwater < nreceiverlayer)
      {
        Jw_matrix2(nreceiverlayer-iwater, &A[iwater], Jmatrix, mat66);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2);
        /* copy back to the surftmp */
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }

      x1_response[k] = gsl_complex_mul(*(x1_response+k), surftmp[0]);
      x3_response[k] = gsl_complex_mul(*(x3_response+k), surftmp[2]);
    }

    t00 = t0 - t_first;
    trace_shift_amp(t00, dt, nfft, 1, x1_response);
    trace_shift_amp(t00, dt, nfft, 1, x3_response);

    if(inorm ==1)
    {
      trace_inv_fft(x3_response, nfft, dxfft0);
      /* norm = -1./dxfft0[ind0]; */
      norm = -1;
      for(k=0; k < nfft; k++)
      {
        if( -dxfft0[k] > norm ){ norm = -dxfft0[k]; }
      }
      norm = 1./norm; 
      for(k=0; k < nfft; k++)
      {
        x1_response[k] = gsl_complex_mul_real(x1_response[k], norm);
      }
    }

    /* add it to the stacked frequency trace */
    for(k=0; k < nfft; k++)
    {
      *(x1_response_stack+k) = gsl_complex_add(x1_response[k], *(x1_response_stack+k));
    }

    M_00002: ;
  }

  /* inverse fourier transform */
  trace_inv_fft(x1_response_stack, nfft, dxfft0);
  for(k=0; k< nfft; k++){ *(xtrace_r+k) = (float) *(dxfft0+k); }

  free(dxfft0);
  return 0;
}

int make_synthetic_ACF_3comp_iso(bool ivelocity, int iinput, int nray, double i_rayp[maxray], int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, double dt, double gaussf, double wsmooth, float *xtrace_r, float *xtrace_t, float *xtrace_z, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw, gsl_complex *iwarray)
{
  int nreceiverlayer = 1;
  int errn, i, j, k, ii, jj, ind0;
  double t0, t_first=0;
  double q3sol[6*nlayer], dvec[3*6], q0tmp[3],  *dxfft0;
  double norm = 0;
  gsl_complex tmpcomp, xfft0[nfft], x1_response[ nfft], x2_response[ nfft], x3_response[ nfft], surftmp[6], surftmp2[6];

  t0 = (nfft*dt)/2; 
  errn = 0;
  ind0 = (int) t0/dt;
  /* make input waveform */
  dxfft0 = (double *)calloc(sizeof(double *), nfft);

  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(xfft0+i, 1., 0.);
    /* initialize xtrace */
    xtrace_r[i] = 0.; xtrace_t[i] = 0.; xtrace_z[i] = 0.;
  }

  /* onto main calculation */
  for(i=0; i < nray; i++)
  {
    errn=0;
    *(q0tmp+0) = *(i_rayp+i);
    *(q0tmp+1) = 0.;

    t_first=0.;

    trace_copy(nfft, xfft0, x1_response);
    trace_copy(nfft, xfft0, x2_response);
    trace_copy(nfft, xfft0, x3_response);

    if(iwater==1)
    {
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;}
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] );
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha;

      E_matrix(m_chris, *(q0tmp+0), q3sol, dvec, tmpE[0]);
      gsl_matrix_complex_set(tmpE_w, 0, 0, gsl_matrix_complex_get(tmpE[0], 2, 0)); 
      gsl_matrix_complex_set(tmpE_w, 1, 0, gsl_matrix_complex_get(tmpE[0], 5, 0)); /* for n = 1 */
      gsl_matrix_complex_set(tmpE_w, 0, 1, gsl_matrix_complex_get(tmpE[0], 2, 3)); 
      gsl_matrix_complex_set(tmpE_w, 1, 1, gsl_matrix_complex_get(tmpE[0], 5, 3)); /* for n = 4 */
      Einv_w_matrix(tmpE_w, tmpinvE_w);
    }

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2_iso(*(q0tmp+0), m_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ // fprintf(stderr, "complex root coming out!!! ACF-3\n"); 
//        fprintf(stderr, "ACF-3 %d layer density %f ray parameter %f\n", j, m_rho[j], i_rayp[i]);
        return -1;
//        goto ACF_00001;
      }
      errn = cal_oscillvec_iso(m_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); return -2; /*goto ACF_00001; */}
      E_matrix(m_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE[j]);
      Einv_matrix(tmpE[j], tmpinvE[j]);
    }

    /* calculate expw */
    cal_expw(nlayer-1, q3sol, m_thick, dw, nfft, expw);

    for(k=0; k < nfft; k++)
    {
      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE[j], tmpinvE[j], &expw[ 6*nfft*j + 6*k], A[j], mat66);
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE[0], tmpinvE_w, &expw[6*k], A_w, mat22); /* j = 0 */
      }

      Jw_matrix(nlayer-iwater, tmpinvE[nlayer-1], &A[iwater], Jmatrix, mat66);

      if(iwater==0){  get_response(iinput, Jmatrix, surftmp, J11, x, in, p); }
      else /* if iwater ==1 */
      {
        gsl_matrix_complex_set_zero(Jmatrix_w);

        for(ii=0; ii < 3; ii++)
        { 
          for(jj=0; jj<3; jj++){ gsl_matrix_complex_set(Jmatrix_w, ii, jj, gsl_matrix_complex_get(Jmatrix, ii, jj)); }
          gsl_matrix_complex_set(Jmatrix_w, ii, 3, gsl_matrix_complex_get(Jmatrix, ii, 5));
        }

        GSL_SET_COMPLEX(&tmpcomp,-1, 0);
        gsl_matrix_complex_set(Jmatrix_w, 3, 2, tmpcomp); /* at (4, 3) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 3, tmpcomp); /* at (5, 4) */
        gsl_matrix_complex_set(Jmatrix_w, 3, 4, gsl_matrix_complex_get(A_w, 0,0)); /* put A11 to (4, 5) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 4, gsl_matrix_complex_get(A_w, 1,0)); /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp, J11, x, in, p);
      }

      if(iwater < nreceiverlayer)
      {
        Jw_matrix2(nreceiverlayer-iwater, &A[iwater], Jmatrix, mat66);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2);
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }

      /* now make the velocity */
      if(ivelocity)
      { 
        surftmp[0] = gsl_complex_mul( surftmp[0], iwarray[k]); 
        surftmp[1] = gsl_complex_mul( surftmp[1], iwarray[k]);
        surftmp[2] = gsl_complex_mul( surftmp[2], iwarray[k]);
      }
      x1_response[k] = gsl_complex_mul(*(x1_response+k), surftmp[0]);
      x2_response[k] = gsl_complex_mul(*(x2_response+k), surftmp[1]);
      x3_response[k] = gsl_complex_mul(*(x3_response+k), surftmp[2]);
    }

    if(wsmooth > 0)
    {
      autocorr_smooth_p2017(nfft, dt, wsmooth, x1_response);
      autocorr_smooth_p2017(nfft, dt, wsmooth, x2_response);
      autocorr_smooth_p2017(nfft, dt, wsmooth, x3_response);
    }
    else /* when there is no smoothing */
    {
      autocorr(nfft, x1_response);
      autocorr(nfft, x2_response);
      autocorr(nfft, x3_response);
    }

    if(gaussf >0)
    {
      norm = gaussian(gaussf, nfft, dt, x1_response, x2_response, x3_response);
    }

    trace_shift_amp(t0, dt, nfft, 1, x1_response);
    trace_shift_amp(t0, dt, nfft, 1, x2_response);
    trace_shift_amp(t0, dt, nfft, 1, x3_response);

    /* inverse fourier transform & add with normalization */
    trace_inv_fft(x1_response, nfft, dxfft0);
    norm = 1./dxfft0[nfft/2];
    for(k=0; k< nfft; k++){ *(xtrace_r+k) += (float) *(dxfft0+k)*norm; }

    trace_inv_fft(x2_response, nfft, dxfft0);
    norm = 1./dxfft0[nfft/2];
    for(k=0; k< nfft; k++){ *(xtrace_t+k) += (float) *(dxfft0+k)*norm; }

    trace_inv_fft(x3_response, nfft, dxfft0);
    norm = 1./dxfft0[nfft/2];
    for(k=0; k< nfft; k++){ *(xtrace_z+k) += (float) *(dxfft0+k)*norm; }

//    ACF_00001: ;
  }

  free(dxfft0);

  return 0;
}

/* make R-ACF */
/* autocorrelation function */
int make_synthetic_ACF_R_iso_noise_ver6(int iinput, INDATA racf, int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, float *xtrace_r, float *indv_r, int synthind_0, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw, gsl_complex *shift, double *w2array)
{
  int tmpint, tmp_nsnrs, nreceiverlayer = iwater;
  int errn, i, j, k, ii, jj, ind0;
  double t0, t_first=0;
  double q3sol[6*nlayer], dvec[3*6], q0tmp[3];
  double tmpdouble, norm = 0;
  gsl_complex tmpcomp, xfft0[nfft], x1_response[ nfft], surftmp[6], surftmp2[6];

  t0 = (nfft*racf.delta)/2;
  errn = 0;
  ind0 = (int) t0/racf.delta;

  for(i=0; i < nfft; i++)
  { 
    GSL_SET_COMPLEX(xfft0+i, 1., 0.);
    /* initialize xtrace */
    xtrace_r[i] = 0.;
  }

  /* onto main calculation */
  tmp_nsnrs = 0;
  for(i=0; i < racf.nrayp; i++)
  {
    errn=0;
    *(q0tmp+0) = racf.rayp[i];
    *(q0tmp+1) = 0;

    t_first=0.;

    trace_copy(nfft, xfft0, x1_response);

    if(iwater==1)
    {
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;}
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] );
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha;

      E_matrix(m_chris, *(q0tmp+0), q3sol, dvec, tmpE[0]);
      gsl_matrix_complex_set(tmpE_w, 0, 0, gsl_matrix_complex_get(tmpE[0], 2, 0)); 
      gsl_matrix_complex_set(tmpE_w, 1, 0, gsl_matrix_complex_get(tmpE[0], 5, 0)); /* for n = 1 */
      gsl_matrix_complex_set(tmpE_w, 0, 1, gsl_matrix_complex_get(tmpE[0], 2, 3)); 
      gsl_matrix_complex_set(tmpE_w, 1, 1, gsl_matrix_complex_get(tmpE[0], 5, 3)); /* for n = 4 */
      Einv_w_matrix(tmpE_w, tmpinvE_w);
    }

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2_iso(*(q0tmp+0), m_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ // fprintf(stderr, "complex root coming out!!! ACF-R\n"); 
//        fprintf(stderr, "ACF-R %d layer density %f ray parameter %f\n", j, m_rho[j], racf.rayp[i]);
        return -1;
        /* goto ACF_00002; */
      }
      errn = cal_oscillvec_iso(m_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); return -2; /*goto ACF_00002; */}
      E_matrix(m_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE[j]);
      Einv_matrix(tmpE[j], tmpinvE[j]);
    }

    /* calculate expw */
    cal_expw(nlayer-1, q3sol, m_thick, dw, nfft, expw);

    for(k=0; k < nfft; k++)
    {
      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE[j], tmpinvE[j], expw + 6*nfft*j + 6*k, A[j], mat66);
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE[0], tmpinvE_w, expw + 6*k, A_w, mat22); /* j = 0 */
      }
      Jw_matrix(nlayer-iwater,  tmpinvE[nlayer-1], &A[iwater], Jmatrix, mat66);

      if(iwater==0){  get_response(iinput, Jmatrix, surftmp, J11, x, in,  p); }
      else /* if iwater ==1 */
      {
        gsl_matrix_complex_set_zero(Jmatrix_w);

        for(ii=0; ii < 3; ii++)
        {
          for(jj=0; jj<3; jj++){ gsl_matrix_complex_set(Jmatrix_w, ii, jj, gsl_matrix_complex_get(Jmatrix, ii, jj)); }
          gsl_matrix_complex_set(Jmatrix_w, ii, 3, gsl_matrix_complex_get(Jmatrix, ii, 5));
        }

        GSL_SET_COMPLEX(&tmpcomp,-1, 0);
        gsl_matrix_complex_set(Jmatrix_w, 3, 2, tmpcomp); /* at (4, 3) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 3, tmpcomp); /* at (5, 4) */
        gsl_matrix_complex_set(Jmatrix_w, 3, 4, gsl_matrix_complex_get(A_w, 0,0)); /* put A11 to (4, 5) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 4, gsl_matrix_complex_get(A_w, 1,0)); /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp, J11, x, in, p);
      }

      if(iwater < nreceiverlayer)
      {
        Jw_matrix2(nreceiverlayer-iwater, &A[iwater], Jmatrix, mat66);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2);
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }
      /* now make the velocity if necessaru */
//      if(racf.ivelocity){ surftmp[0] = gsl_complex_mul( surftmp[0], iwarray[k]); }
      x1_response[k] = gsl_complex_mul(*(x1_response+k), surftmp[0]);
    }

    /* calculate autocorrelation */
    tmpdouble = autocorr_ver2(nfft, x1_response, racf.ivelocity, w2array);
    if(racf.w > 0)
    {
      autocorr_smooth_t2019_proto_ver2(nfft, racf.delta, racf.w, x1_response, WSdouble0); 
    } /* if so, prepare the whitening denominator */

    for( ii = 0; ii < racf.nsnrs[i]; ii++)
    {
      if( racf.w > 0 ){ autocorr_smooth_t2019_final_noise_ver2(nfft, x1_response, WSdouble0, tmpdouble, racf.snrs_factor[tmp_nsnrs+ii], WScomp); }
      else{ autocorr_noise_byadding(nfft, tmpdouble, racf.snrs_factor[tmp_nsnrs+ii], x1_response, WScomp); }
      /*if( racf.w > 0 ){ autocorr_smooth_t2019_final_noise_ver3(nfft, x1_response, WSdouble0, tmpdouble, racf.snrs_factor[tmp_nsnrs+ii], WScomp, racf.ivelocity, w2array); }
      else{ autocorr_noise_byadding_ver3(nfft, tmpdouble, racf.snrs_factor[tmp_nsnrs+ii], x1_response, WScomp, racf.ivelocity, w2array); } */

      /*  post-processing */
      if(racf.gauss >0){ norm = gaussian_1comp(racf.gauss, nfft, racf.delta, WScomp); }
      trace_shift_with_factor(nfft, shift, WScomp);
      trace_inv_fft(WScomp, nfft, WSdouble);
      for(jj=0; jj < nfft; jj++){ WSfloat[jj] = (float) WSdouble[jj]; }
      filter_signal(nfft, racf.delta, racf.ffilter[0], racf.ffilter[1], WSfloat);

      norm = 1./WSfloat[nfft/2];
      /* stack or not?! */
      for(jj=0; jj < nfft; jj++){ *(xtrace_r + jj) += WSfloat[jj]*norm; }
      if(!racf.stack_tf)
      { 
        tmpint = (tmp_nsnrs + ii)*racf.npts;
        for(jj=0; jj < racf.npts; jj++){ 
          *(indv_r + tmpint + jj) = WSfloat[jj + synthind_0]*norm; 
        }  
      }
    }

//    ACF_00002: ;
    tmp_nsnrs += racf.nsnrs[i];
  }

  return 0;
}

/* for ACF Z component */
int make_synthetic_ACF_Z_iso_noise_ver6(int iinput, INDATA zacf, int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, float *xtrace_z,  float *indv_z, int synthind_0, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw,  gsl_complex *shift, double *w2array)
{
  int tmpint, tmp_nsnrs, nreceiverlayer = iwater;
  int errn, i, j, k, ii, jj, ind0;
  double tmpdouble, t0, t_first=0;
  double q3sol[6*nlayer], dvec[3*6], q0tmp[3];
  double norm = 0;
  gsl_complex tmpcomp, xfft0[nfft], x3_response[ nfft], surftmp[6], surftmp2[6];

  t0 = (nfft*zacf.delta)/2;
  errn = 0;
  ind0 = (int) t0/zacf.delta;

  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(xfft0+i, 1., 0.);
    /* initialize xtrace */
    xtrace_z[i] = 0.;
  }

  /* onto main calculation */
  tmp_nsnrs = 0;
  for(i=0; i < zacf.nrayp; i++)
  {
    errn=0;
    *(q0tmp+0) = zacf.rayp[i];
    *(q0tmp+1) = 0.;

    t_first=0.;

    trace_copy(nfft, xfft0, x3_response);

    if(iwater==1)
    {
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;}
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] );
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha;

      E_matrix(m_chris, *(q0tmp+0), q3sol, dvec, tmpE[0]);
      gsl_matrix_complex_set(tmpE_w, 0, 0, gsl_matrix_complex_get(tmpE[0], 2, 0));
      gsl_matrix_complex_set(tmpE_w, 1, 0, gsl_matrix_complex_get(tmpE[0], 5, 0)); /* for n = 1 */
      gsl_matrix_complex_set(tmpE_w, 0, 1, gsl_matrix_complex_get(tmpE[0], 2, 3));
      gsl_matrix_complex_set(tmpE_w, 1, 1, gsl_matrix_complex_get(tmpE[0], 5, 3)); /* for n = 4 */

      Einv_w_matrix(tmpE_w, tmpinvE_w);
    }

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2_iso(*(q0tmp+0), m_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ // fprintf(stderr, "complex root coming out!!! ACF-Z\n"); 
//        fprintf(stderr, "ACF-Z %d layer density %f ray parameter %f\n", j, m_rho[j], zacf.rayp[i]);
        return -1;
        /* goto ACF_00003; */
      }
      errn = cal_oscillvec_iso(m_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); return -2; /*goto ACF_00003;*/ }
      E_matrix(m_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE[j]);
      Einv_matrix(tmpE[j], tmpinvE[j]);
    }

    /* calculate expw */
    cal_expw(nlayer-1, q3sol, m_thick, dw, nfft, expw);

    for(k=0; k < nfft; k++)
    {
      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE[j], tmpinvE[j], expw + 6*nfft*j + 6*k, A[j], mat66);
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE[0], tmpinvE_w, expw + 6*k, A_w, mat22); /* j = 0 */
      }

      Jw_matrix(nlayer-iwater,  tmpinvE[nlayer-1], &A[iwater], Jmatrix, mat66);

      if(iwater==0){  get_response(iinput, Jmatrix, surftmp, J11, x, in, p); }
      else /* if iwater ==1 */
      {
        gsl_matrix_complex_set_zero(Jmatrix_w);

        for(ii=0; ii < 3; ii++)
        {
          for(jj=0; jj<3; jj++){ gsl_matrix_complex_set(Jmatrix_w, ii, jj, gsl_matrix_complex_get(Jmatrix, ii, jj)); }
          gsl_matrix_complex_set(Jmatrix_w, ii, 3, gsl_matrix_complex_get(Jmatrix, ii, 5));
        }

        GSL_SET_COMPLEX(&tmpcomp,-1, 0);
        gsl_matrix_complex_set(Jmatrix_w, 3, 2, tmpcomp); /* at (4, 3) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 3, tmpcomp); /* at (5, 4) */
        gsl_matrix_complex_set(Jmatrix_w, 3, 4, gsl_matrix_complex_get(A_w, 0,0)); /* put A11 to (4, 5) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 4, gsl_matrix_complex_get(A_w, 1,0)); /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp, J11, x, in, p);
      }

      if(iwater < nreceiverlayer)
      {
        Jw_matrix2(nreceiverlayer-iwater, &A[iwater], Jmatrix, mat66);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2);
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }
      
      /* now make the velocity */
      //if(zacf.ivelocity){ surftmp[2] = gsl_complex_mul( surftmp[2], iwarray[k]); }
      x3_response[k] = gsl_complex_mul(*(x3_response+k), surftmp[2]);
    }


    /* calculate autocorrelation */
    tmpdouble = autocorr_ver2(nfft, x3_response, zacf.ivelocity, w2array);
    if(zacf.w > 0)
    {
      autocorr_smooth_t2019_proto_ver2(nfft, zacf.delta, zacf.w, x3_response, WSdouble0);
    } /* if so, prepare the whitening denominator */

    for( ii = 0; ii < zacf.nsnrs[i]; ii++)
    {
/*      if( zacf.w > 0 ){ autocorr_smooth_t2019_final_noise_ver3(nfft, x3_response, WSdouble0, tmpdouble, zacf.snrs_factor[tmp_nsnrs+ii], WScomp, zacf.ivelocity, w2array); }
      else{ autocorr_noise_byadding_ver3(nfft, tmpdouble, zacf.snrs_factor[tmp_nsnrs+ii], x3_response, WScomp, zacf.ivelocity, w2array); } */
      if( zacf.w > 0 ){ autocorr_smooth_t2019_final_noise_ver2(nfft, x3_response, WSdouble0, tmpdouble, zacf.snrs_factor[tmp_nsnrs+ii], WScomp); }
      else{ autocorr_noise_byadding(nfft, tmpdouble, zacf.snrs_factor[tmp_nsnrs+ii], x3_response, WScomp); }


      /*  post-processing */
      if(zacf.gauss >0){ norm = gaussian_1comp(zacf.gauss, nfft, zacf.delta, WScomp); }
      trace_shift_with_factor(nfft, shift, WScomp);
      trace_inv_fft(WScomp, nfft, WSdouble);
      for(jj=0; jj < nfft; jj++){ WSfloat[jj] = (float) WSdouble[jj]; }
      filter_signal(nfft, zacf.delta, zacf.ffilter[0], zacf.ffilter[1], WSfloat);

      norm = 1./WSfloat[nfft/2];
      /* stack or not?! */
      for(jj=0; jj < nfft; jj++){ *(xtrace_z + jj) += WSfloat[jj]*norm; }
      if(!zacf.stack_tf)
      { 
        tmpint = (tmp_nsnrs + ii)*zacf.npts;
        for(jj=0; jj < zacf.npts; jj++){ 
          *(indv_z + tmpint + jj) = WSfloat[jj + synthind_0]*norm; 
	}
      }
    }

//    ACF_00003: ;
    tmp_nsnrs += zacf.nsnrs[i]; 
  }

  return 0;
}

/* for ACF DPG component : iwater ==1 only. but to unify the form, I give iwater */ 
int make_synthetic_ACF_DPG_iso_noise_ver6(int iinput, INDATA pacf, int nlayer, float m_thick[], float m_rho[], double m_chris[], int iwater, float water_alpha, int nfft, double dw, float *xtrace_p,  float *indv_p, int synthind_0, gsl_matrix_complex **tmpE, gsl_matrix_complex **tmpinvE, gsl_matrix_complex *Jmatrix, gsl_matrix_complex *Jmatrix_w, gsl_matrix_complex **A, gsl_matrix_complex *tmpE_w, gsl_matrix_complex *tmpinvE_w, gsl_matrix_complex *A_w, gsl_matrix_complex *mat22, gsl_matrix_complex *mat66, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p, gsl_complex *expw,  gsl_complex *shift, double *w2array)
{
  int tmpint, tmp_nsnrs, nreceiverlayer = 1;
  int errn, i, j, k, ii, jj, ind0;
  double tmpdouble, t0, t_first=0;
  double q3sol[6*nlayer], dvec[3*6], q0tmp[3];
  double norm = 0;
  gsl_complex tmpcomp, xfft0[nfft], xp_response[ nfft], surftmp[6], surftmp2[6];

  t0 = (nfft*pacf.delta)/2;
  errn = 0;
  ind0 = (int) t0/pacf.delta;

  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(xfft0+i, 1., 0.);
    /* initialize xtrace */
    xtrace_p[i] = 0.;
  }

  /* onto main calculation */
  tmp_nsnrs = 0;
  for(i=0; i < pacf.nrayp; i++)
  {
    errn=0;
    *(q0tmp+0) = pacf.rayp[i];
    *(q0tmp+1) = 0.;

    t_first=0.;

    trace_copy(nfft, xfft0, xp_response);

    if(iwater==1)
    {
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;}
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] );
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha;

      E_matrix(m_chris, *(q0tmp+0), q3sol, dvec, tmpE[0]);
      gsl_matrix_complex_set(tmpE_w, 0, 0, gsl_matrix_complex_get(tmpE[0], 2, 0));
      gsl_matrix_complex_set(tmpE_w, 1, 0, gsl_matrix_complex_get(tmpE[0], 5, 0)); /* for n = 1 */
      gsl_matrix_complex_set(tmpE_w, 0, 1, gsl_matrix_complex_get(tmpE[0], 2, 3));
      gsl_matrix_complex_set(tmpE_w, 1, 1, gsl_matrix_complex_get(tmpE[0], 5, 3)); /* for n = 4 */

      Einv_w_matrix(tmpE_w, tmpinvE_w);
    }

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2_iso(*(q0tmp+0), m_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ // fprintf(stderr, "complex root coming out!!! ACF-P\n");
//        fprintf(stderr, "ACF-P %d layer density %f ray parameter %f\n", j, m_rho[j], pacf.rayp[i]);
        return -1;
        /* goto ACF_00004; */
      }
      errn = cal_oscillvec_iso(m_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); return -2; /*goto ACF_00003;*/ }
      E_matrix(m_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE[j]);
      Einv_matrix(tmpE[j], tmpinvE[j]);
    }

    /* calculate expw */
    cal_expw(nlayer-1, q3sol, m_thick, dw, nfft, expw);

    for(k=0; k < nfft; k++)
    {
      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE[j], tmpinvE[j], expw + 6*nfft*j + 6*k, A[j], mat66);
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE[0], tmpinvE_w, expw + 6*k, A_w, mat22); /* j = 0 */
      }

      Jw_matrix(nlayer-iwater, tmpinvE[nlayer-1], &A[iwater], Jmatrix, mat66);

      if(iwater==0){  get_response(iinput, Jmatrix, surftmp, J11, x, in, p); }
      else /* if iwater ==1 */
      {
        gsl_matrix_complex_set_zero(Jmatrix_w);

        for(ii=0; ii < 3; ii++)
        {
          for(jj=0; jj<3; jj++){ gsl_matrix_complex_set(Jmatrix_w, ii, jj, gsl_matrix_complex_get(Jmatrix, ii, jj)); }
          gsl_matrix_complex_set(Jmatrix_w, ii, 3, gsl_matrix_complex_get(Jmatrix, ii, 5));
        }

        GSL_SET_COMPLEX(&tmpcomp,-1, 0);
        gsl_matrix_complex_set(Jmatrix_w, 3, 2, tmpcomp); /* at (4, 3) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 3, tmpcomp); /* at (5, 4) */
        gsl_matrix_complex_set(Jmatrix_w, 3, 4, gsl_matrix_complex_get(A_w, 0,0)); /* put A11 to (4, 5) */
        gsl_matrix_complex_set(Jmatrix_w, 4, 4, gsl_matrix_complex_get(A_w, 1,0)); /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp, J11, x, in, p);
      }

      if(iwater < nreceiverlayer)
      {
        Jw_matrix2(nreceiverlayer-iwater, &A[iwater], Jmatrix, mat66);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2);
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }

      /* now make the dpg */
//      surftmp[5] = gsl_complex_mul( surftmp[5], iwarray[k]); 
      xp_response[k] = gsl_complex_mul(*(xp_response+k), surftmp[5]);
    }

    /* calculate autocorrelation -- DPG-ACF */
    tmpdouble = autocorr_ver2(nfft, xp_response, pacf.ivelocity, w2array); /* always differentiation considered */
    if(pacf.w > 0)
    {
      autocorr_smooth_t2019_proto_ver2(nfft, pacf.delta, pacf.w, xp_response, WSdouble0);
    } /* if so, prepare the whitening denominator */

    for( ii = 0; ii < pacf.nsnrs[i]; ii++)
    {
/*      if( pacf.w > 0 ){ autocorr_smooth_t2019_final_noise_ver3(nfft, xp_response, WSdouble0, tmpdouble, pacf.snrs_factor[tmp_nsnrs+ii], WScomp, true, w2array); }
      else{ autocorr_noise_byadding_ver3(nfft, tmpdouble, pacf.snrs_factor[tmp_nsnrs+ii], xp_response, WScomp, true, w2array); } */
      if( pacf.w > 0 ){ autocorr_smooth_t2019_final_noise_ver2(nfft, xp_response, WSdouble0, tmpdouble, pacf.snrs_factor[tmp_nsnrs+ii], WScomp); }
      else{ autocorr_noise_byadding(nfft, tmpdouble, pacf.snrs_factor[tmp_nsnrs+ii], xp_response, WScomp); }

      /*  post-processing */
      if(pacf.gauss >0){ norm = gaussian_1comp(pacf.gauss, nfft, pacf.delta, WScomp); }
      trace_shift_with_factor(nfft, shift, WScomp);
      trace_inv_fft(WScomp, nfft, WSdouble);
      for(jj=0; jj < nfft; jj++){ WSfloat[jj] = (float) WSdouble[jj]; }
      filter_signal(nfft, pacf.delta, pacf.ffilter[0], pacf.ffilter[1], WSfloat);

      norm = 1./WSfloat[nfft/2];
      /* stack or not?! */
      for(jj=0; jj < nfft; jj++){ *(xtrace_p + jj) += WSfloat[jj]*norm; }
      if(!pacf.stack_tf)
      {
        tmpint = (tmp_nsnrs + ii)*pacf.npts;
        for(jj=0; jj < pacf.npts; jj++){
          *(indv_p + tmpint + jj) = WSfloat[jj + synthind_0]*norm;
        }
      }
    }

//    ACF_00004: ;
    tmp_nsnrs += pacf.nsnrs[i];
  }

  return 0;
}

