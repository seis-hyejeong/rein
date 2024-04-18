#include "func_likelihood-2.h"

/* calculating the likelihood of the model
   for each dataset separately. Merge it later. 
log(P(d|m)) = SUM( -log(sigma_i) ) - SUM( 1/(2*sigma_i^2)*(obs_i - syn_i)^2 )
                 cancels                    remains
log(P(d|m')) - log(P(d|m))
= SUM( -1/(2*sigma_i^2)*( (obs_i - syn'_i)^2 - (obs_i - syn_i)^2 ) )
*/

/*
 second term: -1/(2*sigma_i^2)
*/
void model_likelihood_log_sigmasquare(int ndata, float sigma[], float *inv_sigma)
{
  int i;
  if(ndata==1){ return; }
 
  for(i=0; i < ndata; i++)
  {
    *(inv_sigma+i) = -0.5/(sigma[i]*sigma[i]);
  }

  return;
}

/* SUM( -1/(2*sigma_i^2)*(obs_i - syn'_i)^2 */
double model_likelihood_log_term2(int ndata, float synth[], float obs[], float inv_sigma[])
{
  int i;
  double tmpdiff, likelihood = 0;
  for(i=0; i < ndata; i++)
  {
    tmpdiff = obs[i] - synth[i];
    likelihood += (inv_sigma[i]*tmpdiff*tmpdiff); 
  }

  return likelihood;
}

/* calculate the likelihood with the covariance matrix
 */ 
double model_likelihood_log_covmat_ldlt(int ndata, float synth[], float obs[], gsl_matrix *L, gsl_vector *Vec, gsl_vector *Yvec)
{
  size_t ii, jj;
  double tmpdiff, likelihood = 0;

  /* now calculate the multiple -- the L is in lower triangular */
  /* now get the likelihood */
  likelihood = 0;
  for(ii=0; ii < ndata; ii++)
  {
    tmpdiff = obs[ii] - synth[ii];
    for(jj = 0; jj < ii; jj++)
    { 
      tmpdiff -= ( gsl_matrix_get(L, ii, jj)*gsl_vector_get(Yvec, jj) ); 
    }
    tmpdiff *=  gsl_vector_get(Vec, ii); 
    gsl_vector_set(Yvec, ii, tmpdiff);

    likelihood -=  tmpdiff*tmpdiff*0.5; 
  }

  if( isnan(likelihood) ){ return 1.; }  
  return likelihood; 
}


double model_likelihood_log_covmat2(int ndata, float synth[], float obs[], gsl_matrix *L, gsl_permutation *per, gsl_vector *vec, gsl_vector *yvec)
{
  size_t i;
  double tmpdiff, likelihood = 0;

  /* get the difference vector */
  for(i = 0; i < ndata; i++)
  {
    tmpdiff = obs[i] - synth[i];
    gsl_vector_set( vec, i, tmpdiff); 
  }
 
  /* now calculate the multiple -- the L is in lower triangular */
  i = gsl_linalg_LU_solve(L, per, vec, yvec);
  if( i != 0 ){ return 1.; }

  /* now get the likelihood */
  likelihood = 0;
  for(i=0; i < ndata; i++)
  {
    likelihood += gsl_vector_get(vec, i)*gsl_vector_get(yvec, i);
  }
  likelihood *= (-0.5);

  if( isnan(likelihood) ){ return 1.; }
  return likelihood; 
}

/* 
 C = P^(-1)LU
 C^(-1) = U^(-1) L^(-1) P
 solve v^(t) U^(-1) and L^(-1) (P v)
*/
double model_likelihood_log_covmat_lu(int ndata, float synth[], float obs[], gsl_matrix *L, gsl_permutation *per, gsl_vector *vec, gsl_vector *yvec)
{
  size_t ii, jj;
  double tmpdiff, likelihood = 0;

  /* get the difference vector */
  for(ii = 0; ii < ndata; ii++)
  {
    tmpdiff = obs[ii] - synth[ii];
    gsl_vector_set( vec, ii, tmpdiff);
    gsl_vector_set( yvec, ii, tmpdiff);
  }

  gsl_permute_vector(per, vec);
  /* solve the lower matrix with permutation applied vector : L has diagonal elements of unity*/
  for(ii=0; ii < ndata; ii++)
  {
    tmpdiff = gsl_vector_get(vec, ii);
 
    for(jj = 0; jj < ii; jj++)
    {
      tmpdiff -= ( gsl_matrix_get(L, ii, jj)*gsl_vector_get(vec, jj) );
    }
    gsl_vector_set(vec, ii, tmpdiff); /* replace */
  }

  /* now solve the upper matrix with v^(t): diagonal not unity */
  for(ii=0; ii < ndata; ii++)
  {
    tmpdiff = obs[ii] - synth[ii];
    for(jj = 0; jj < ii; jj++)
    {
      tmpdiff -= ( gsl_matrix_get(L, jj, ii)*gsl_vector_get(yvec, jj) );
    }
    tmpdiff *= gsl_matrix_get(L, ii, ii);

    gsl_vector_set(yvec, ii, tmpdiff);
  }
 
  /* calculate the multiple and multiple by -0.5 */
  tmpdiff = 0;
  for(ii=0; ii < ndata; ii++)
  {
    tmpdiff += ( gsl_vector_get(yvec, ii)*gsl_vector_get(vec, ii) );
  }
  likelihood = tmpdiff*(-0.5);

  if( isnan(likelihood) ){ return 1; }

  return likelihood;
}


/* final likelihood joining R-ACF, Z-ACF, and RF (or waveform stack) */
void model_likelihood_log_alpha(LIKELI *p)
{
  double log_alpha = 0;
  
  log_alpha =  (*p).racf + (*p).zacf + (*p).pacf + (*p).rwave;
  (*p).pall = (float) log_alpha;
 
  return; 
}

int print_likelihood(int ith_cal, FILE *filestr, PARATEMP parallelt)
{
  int i, count;
 
  fwrite(&ith_cal, sizeof(int), 1, filestr);
  count = 0; i = 0;
  while(i < parallelt.nprocs && count < parallelt.ncool)
  {
    if(parallelt.Ts[i] == 1.){ 
      fwrite(&parallelt.likelihood[i], sizeof(float), 1, filestr); 
      count++;
    }
    i++;
  }

  if( count != parallelt.ncool){ fprintf(stderr, "count %d and ncool %d\n", count, parallelt.ncool); return -1; }
  else{ return 0; }
}

int print_likelihood_T(int ith_cal, FILE *filestr, PARATEMP parallelt, double printT)
{
  int i, count;
   
  fwrite(&ith_cal, sizeof(int), 1, filestr);
 
  count = 0; i = 0;
  while(i < parallelt.nprocs && count < 1)
  {
    if(parallelt.Ts[i] == printT)
    {
      fwrite(&parallelt.likelihood[i], sizeof(float), 1, filestr);
      count++;
    }
    i++;
  }

  if( count != 1){ return -1; }
  else{ return 0; }
}

/* the temperature swapping function */
int swap_temperature(PARATEMP *parallelt, gsl_rng *r)
{
  int nchange, i, p, q;
  double log_alpha, r1r2, invTp, invTq, Tp, Tq, logPp, logPq;

  nchange = 0;
  for(i=0; i < (*parallelt).nprocs; i++)
  {
    p = gsl_rng_get(r) %(*parallelt).nprocs;
    q = gsl_rng_get(r) %(*parallelt).nprocs;
    while(q == p){ q = gsl_rng_get(r) %(*parallelt).nprocs; }
    logPp = (*parallelt).likelihood[p];
    logPq = (*parallelt).likelihood[q];
    
    Tp = (*parallelt).Ts[p];
    Tq = (*parallelt).Ts[q];
   
    invTp = (*parallelt).invTs[p];
    invTq = (*parallelt).invTs[q];

    r1r2 = (invTp - invTq)*(logPq - logPp);
    log_alpha = log(gsl_rng_uniform(r));
    
    if( log_alpha < r1r2 )
    { /* we change! */
      nchange++;
      (*parallelt).Ts[p] = Tq;
      (*parallelt).Ts[q] = Tp;
      (*parallelt).invTs[p] = invTq;
      (*parallelt).invTs[q] = invTp;
    }  
  }

  return nchange;
}
//  synthetics_likelihood_ver3(&likelihood0, s_param,
//  synth_racf, indv_racf, inracfs, inracf_indv, inracfe_inv, Lracf, pracf, Vracf, Yracf, in_racf,
//  synth_zacf, indv_zacf, inzacfs, inzacf_indv, inzacfe_inv, Lzacf, pzacf, Vzacf, Yzacf, in_zacf,
//  synth_pacf, indv_pacf, inpacfs, inpacf_indv, inpacfe_inv, Lpacf, ppacf, Vpacf, Ypacf, in_pacf,
//  synth_rwave, inrwaves, inrwavee_inv, Lrwave, prwave, Vrwave, Yrwave, in_rwave);

/* calculate likelihood considering non-stack case and with DPG */
void synthetics_likelihood_ver3_ldlt(LIKELI *p, SYNTHPAR s_param, float s_racf[], float indv_racf[], float o_racf[], float o_racf_indv[], float inv_racf[], gsl_matrix *mat_racf, gsl_vector *vec_racf, gsl_vector *y_racf, INDATA i_racf, float s_zacf[], float indv_zacf[], float o_zacf[], float o_zacf_indv[], float inv_zacf[], gsl_matrix *mat_zacf, gsl_vector *vec_zacf, gsl_vector *y_zacf, INDATA i_zacf, float s_pacf[], float indv_pacf[], float o_pacf[], float o_pacf_indv[], float inv_pacf[], gsl_matrix *mat_pacf, gsl_vector *vec_pacf, gsl_vector *y_pacf, INDATA i_pacf, float s_rwave[], float o_rwave[], float inv_rwave[], gsl_matrix *mat_rwave, gsl_vector *vec_rwave, gsl_vector *y_rwave, INDATA i_rwave)
{ 
  int i;
  double logP = 0;

  /* initialize */
  (*p).racf = 0.; (*p).zacf = 0.; (*p).pacf = 0.; (*p).rwave = 0.;
  /* R - ACF */
  if(i_racf.tf)
  {
    logP = 0;
    if(i_racf.stack_tf && i_racf.covar_tf) /* covariance */
    {
/*double model_likelihood_log_covmat_ldlt(int ndata, float synth[], float obs[], gsl_matrix *L, gsl_vector *Vec, gsl_vector *Yvec) */
      logP = model_likelihood_log_covmat_ldlt(i_racf.npts, s_racf + s_param.racfind[0], o_racf + i_racf.tind0, mat_racf, vec_racf, y_racf); 
    }
    else if(i_racf.stack_tf) /* no covariance, but using the stack */
    { 
      logP = model_likelihood_log_term2(i_racf.npts, s_racf + s_param.racfind[0], o_racf + i_racf.tind0, inv_racf);  /* stacked -- compare only with the stacked */
    }
    else
    {
      for(i=0; i < i_racf.nfiles; i++){ logP += model_likelihood_log_term2(i_racf.npts, indv_racf + i_racf.startind[i], o_racf_indv + i_racf.startind[i], inv_racf); }
    }
    (*p).racf = logP;
  }
  /* Z - ACF */
  if(i_zacf.tf)
  {
    logP = 0;
    if(i_zacf.stack_tf && i_zacf.covar_tf)
    {
      logP = model_likelihood_log_covmat_ldlt(i_zacf.npts, s_zacf + s_param.zacfind[0], o_zacf + i_zacf.tind0, mat_zacf, vec_zacf, y_zacf); 
    }
    else if(i_zacf.stack_tf)
    { 
      logP = model_likelihood_log_term2(i_zacf.npts, s_zacf + s_param.zacfind[0], o_zacf + i_zacf.tind0, inv_zacf); 
    }
    else
    {
      for(i=0; i < i_zacf.nfiles; i++){ logP += model_likelihood_log_term2(i_zacf.npts, indv_zacf + i_zacf.startind[i], o_zacf_indv + i_zacf.startind[i], inv_zacf); }
    }
    (*p).zacf = logP;
  }
  /* P - ACF */
  if(i_pacf.tf)
  {
    logP = 0;
    if( i_pacf.stack_tf && i_pacf.covar_tf )
    {
      logP = model_likelihood_log_covmat_ldlt( i_pacf.npts, s_pacf + s_param.pacfind[0], o_pacf + i_pacf.tind0, mat_pacf, vec_pacf, y_pacf );
    }
    else if(i_pacf.stack_tf)
    { 
      logP = model_likelihood_log_term2(i_pacf.npts, s_pacf + s_param.pacfind[0], o_pacf + i_pacf.tind0, inv_pacf); 
    }
    else
    {
      for(i=0; i < i_pacf.nfiles; i++){ logP += model_likelihood_log_term2(i_pacf.npts, indv_pacf + i_pacf.startind[i], o_pacf_indv + i_pacf.startind[i], inv_pacf); }
    }
    (*p).pacf = logP;
  }
  /* R stack */
  if(i_rwave.tf)
  {
    logP = 0;
    if(!i_rwave.stack_tf){ fprintf(stderr, "error in the algorithm of r-stack taking.\n"); exit(-1); }
    
    if( i_rwave.covar_tf ) /* can take the covariance ? */
    {
      logP = model_likelihood_log_covmat_ldlt( i_rwave.npts, s_rwave + s_param.rwaveind[0], o_rwave + i_rwave.tind0, mat_rwave, vec_rwave, y_rwave);
    }
    else /* not taking the covariance */
    { 
      logP = model_likelihood_log_term2(i_rwave.npts, s_rwave + s_param.rwaveind[0], o_rwave + i_rwave.tind0, inv_rwave); 
    }
    (*p).rwave = logP; 
  }

  /* merge all likelihood components */  
  model_likelihood_log_alpha(p);
 
  return;
} 


/* calculate likelihood with LU decomposition matrix */
void synthetics_likelihood_ver3_lu(LIKELI *p, SYNTHPAR s_param, float s_racf[], float indv_racf[], float o_racf[], float o_racf_indv[], float inv_racf[], gsl_matrix *mat_racf, gsl_permutation *per_racf, gsl_vector *vec_racf, gsl_vector *y_racf, INDATA i_racf, float s_zacf[], float indv_zacf[], float o_zacf[], float o_zacf_indv[], float inv_zacf[], gsl_matrix *mat_zacf, gsl_permutation *per_zacf, gsl_vector *vec_zacf, gsl_vector *y_zacf, INDATA i_zacf, float s_pacf[], float indv_pacf[], float o_pacf[], float o_pacf_indv[], float inv_pacf[], gsl_matrix *mat_pacf, gsl_permutation *per_pacf, gsl_vector *vec_pacf, gsl_vector *y_pacf, INDATA i_pacf, float s_rwave[], float o_rwave[], float inv_rwave[], gsl_matrix *mat_rwave, gsl_permutation *per_rwave, gsl_vector *vec_rwave, gsl_vector *y_rwave, INDATA i_rwave)
{
  int i;
  double logP = 0;

  /* initialize */
  (*p).racf = 0.; (*p).zacf = 0.; (*p).pacf = 0.; (*p).rwave = 0.;
  /* R - ACF */
  if(i_racf.tf)
  {
    logP = 0;
    if(i_racf.stack_tf && i_racf.covar_tf) /* covariance */
    {
      logP = model_likelihood_log_covmat_lu(i_racf.npts, s_racf + s_param.racfind[0], o_racf + i_racf.tind0, mat_racf, per_racf, vec_racf, y_racf);
    }
    else if(i_racf.stack_tf) /* no covariance, but using the stack */
    {
      logP = model_likelihood_log_term2(i_racf.npts, s_racf + s_param.racfind[0], o_racf + i_racf.tind0, inv_racf);  /* stacked -- compare only with the stacked */
    }
    else
    {
      for(i=0; i < i_racf.nfiles; i++){ logP += model_likelihood_log_term2(i_racf.npts, indv_racf + i_racf.startind[i], o_racf_indv + i_racf.startind[i], inv_racf); }
    }
    (*p).racf = logP;
  }
  /* Z - ACF */
  if(i_zacf.tf)
  {
    logP = 0;
    if(i_zacf.stack_tf && i_zacf.covar_tf)
    {
      logP = model_likelihood_log_covmat_lu(i_zacf.npts, s_zacf + s_param.zacfind[0], o_zacf + i_zacf.tind0, mat_zacf, per_zacf, vec_zacf, y_zacf);
    }
    else if(i_zacf.stack_tf)
    {
      logP = model_likelihood_log_term2(i_zacf.npts, s_zacf + s_param.zacfind[0], o_zacf + i_zacf.tind0, inv_zacf);
    }
    else
    {
      for(i=0; i < i_zacf.nfiles; i++){ logP += model_likelihood_log_term2(i_zacf.npts, indv_zacf + i_zacf.startind[i], o_zacf_indv + i_zacf.startind[i], inv_zacf); }
    }
    (*p).zacf = logP;
  }
  /* P - ACF */
  if(i_pacf.tf)
  {
    logP = 0;
    if( i_pacf.stack_tf && i_pacf.covar_tf )
    {
      logP = model_likelihood_log_covmat_lu( i_pacf.npts, s_pacf + s_param.pacfind[0], o_pacf + i_pacf.tind0, mat_pacf, per_pacf, vec_pacf, y_pacf );
    }
    else if(i_pacf.stack_tf)
    {
      logP = model_likelihood_log_term2(i_pacf.npts, s_pacf + s_param.pacfind[0], o_pacf + i_pacf.tind0, inv_pacf);
    }
    else
    {
      for(i=0; i < i_pacf.nfiles; i++){ logP += model_likelihood_log_term2(i_pacf.npts, indv_pacf + i_pacf.startind[i], o_pacf_indv + i_pacf.startind[i], inv_pacf); }
    }
    (*p).pacf = logP;
  }
  /* R stack */
  if(i_rwave.tf)
  {
    logP = 0;
    if(!i_rwave.stack_tf){ fprintf(stderr, "error in the algorithm of r-stack taking.\n"); exit(-1); }

    if( i_rwave.covar_tf ) /* can take the covariance ? */
    {
      logP = model_likelihood_log_covmat_lu( i_rwave.npts, s_rwave + s_param.rwaveind[0], o_rwave + i_rwave.tind0, mat_rwave, per_rwave, vec_rwave, y_rwave);
    }
    else /* not taking the covariance */
    {
      logP = model_likelihood_log_term2(i_rwave.npts, s_rwave + s_param.rwaveind[0], o_rwave + i_rwave.tind0, inv_rwave);
    }
    (*p).rwave = logP;
  }

  /* merge all likelihood components */
  model_likelihood_log_alpha(p);

  return;
}


/* calculate likelihood of all */
void synthetics_likelihood(LIKELI *p, SYNTHPAR s_param, float s_racf[], float o_racf[], float inv_racf[], INDATA i_racf, float s_zacf[], float o_zacf[], float inv_zacf[], INDATA i_zacf, float s_rwave[], float o_rwave[], float inv_rwave[], INDATA i_rwave)
{
  if(i_racf.tf)
  {
    (*p).racf = model_likelihood_log_term2(i_racf.npts, s_racf + s_param.racfind[0], o_racf + i_racf.tind0, inv_racf);
  }else{ (*p).racf = 0; }
  if(i_zacf.tf)
  {
    (*p).zacf = model_likelihood_log_term2(i_zacf.npts, s_zacf + s_param.zacfind[0], o_zacf + i_zacf.tind0, inv_zacf);
  }else{ (*p).zacf = 0; }
  if(i_rwave.tf)
  {
    (*p).rwave = model_likelihood_log_term2(i_rwave.npts, s_rwave + s_param.rwaveind[0], o_rwave + i_rwave.tind0, inv_rwave);
  }else{ (*p).rwave = 0; }
  model_likelihood_log_alpha(p);

  return; 
} 

/* the likelihood with normalization */
void synthetics_likelihood_norm(LIKELI *p, SYNTHPAR s_param, float s_racf[], float o_racf[], float inv_racf[], INDATA i_racf, float s_zacf[], float o_zacf[], float inv_zacf[], INDATA i_zacf, float s_rwave[], float o_rwave[], float inv_rwave[], INDATA i_rwave)
{
  if(i_racf.tf)
  {
    (*p).racf = model_likelihood_log_term2(i_racf.npts, s_racf + s_param.racfind[0], o_racf + i_racf.tind0, inv_racf);
    (*p).racf *= (1./i_racf.npts);
  }else{ (*p).racf = 0; }
  if(i_zacf.tf)
  {
    (*p).zacf = model_likelihood_log_term2(i_zacf.npts, s_zacf + s_param.zacfind[0], o_zacf + i_zacf.tind0, inv_zacf);
    (*p).zacf *= (1./i_zacf.npts);
  }else{ (*p).zacf = 0; }
  if(i_rwave.tf)
  {
    (*p).rwave = model_likelihood_log_term2(i_rwave.npts, s_rwave + s_param.rwaveind[0], o_rwave + i_rwave.tind0, inv_rwave);
    (*p).rwave *= (1./i_rwave.npts);
  }else{ (*p).rwave = 0; }
  model_likelihood_log_alpha(p);

  return;
}

