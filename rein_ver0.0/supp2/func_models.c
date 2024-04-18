#include "func_models.h"

/* initialize model */
void initialize_model(MODEL *m)
{
  int i = 0;
  (*m).nl = 0;
  for(i=0; i < n_max; i++)
  {
    (*m).modvec[i] = 0.;
    (*m).h[i] = 0.;
    (*m).vp[i] = 0.;
    (*m).vs[i] = 0.;
    (*m).rho[i]= 0.;
  }
  for(i= n_max; i < 4*n_max; i++){ (*m).modvec[i] = 0.; }

  return;
}

/* replace the model2 with model1
 */
void replace_model(MODEL m1, MODEL *m2)
{
  int i;

  (*m2).nl = m1.nl;
  for(i=0; i < m1.nl; i++) 
  {
    (*m2).modvec[i] = m1.modvec[i];
    (*m2).h[i] = m1.h[i];
    (*m2).vp[i] = m1.vp[i];
    (*m2).vs[i] = m1.vs[i]; 
    (*m2).rho[i] = m1.rho[i];
  }
  i = m1.nl;  /* this is to get full Vp and Vs */
  (*m2).vp[i] = m1.vp[i];
  (*m2).vs[i] = m1.vs[i];
  (*m2).rho[i] = m1.rho[i];
  for(i = m1.nl; i < 4*m1.nl; i++){ (*m2).modvec[i] = m1.modvec[i]; }

  return;
}

/* initialize the model:
 * the MODEL should have values read already 
 */
void first_model(PARAM param, PRIOR prior, MODEL *m)
{
  int i, tmpind; // ilayer;
  double tmploc;

  const gsl_rng_type *T;
  time_t t;  
  srand((unsigned) time(&t));
  int s =rand();
  gsl_rng *r;
  gsl_rng_env_setup();
  T= gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, s);
  
  /* first layer is water layer */
  first_model_0001: ;
  for(i=0; i < prior.nsolve; i++)
  {
    tmpind = prior.mindex[i];
    tmploc = gsl_rng_uniform(r);
    (*m).modvec[tmpind] = (prior.highlim[tmpind] - prior.lowlim[tmpind])*tmploc + prior.lowlim[tmpind];
  }

  /* get rho from vp if we are not solving for density and no density value is given fixed*/
  for(i= param.nl; i < 2*param.nl; i++)
  {
    if(!prior.rho_tf[i-param.nl] && prior.highlim[i+ 2*param.nl] < 0){ (*m).modvec[i+ 2*param.nl] = vp_empirical_rho_solid( (*m).modvec[i] ); } 
    if( prior.p2s[i-param.nl] > 0 ){ (*m).modvec[i +param.nl] = (*m).modvec[i]*prior.p2s[i-param.nl]; } 
  }

  /* change it into vp, vs, rho, h */
  (*m).nl = param.wnl;
  modelvec2inputvecs_general(param, m);

  tmpind = 0; /* for the first model, use monotonically increasing one */
  if(param.increase){
    tmpind = model_prior_monotonic_increase(param.wnl, (*m).vp);
    tmpind += model_prior_monotonic_increase(param.wnl, (*m).vs);
    tmpind += model_prior_monotonic_increase(param.wnl, (*m).rho);
  }
  tmpind += model_vp_and_vs(param.wnl, (*m).vp, (*m).vs);

  if(tmpind !=0){ goto first_model_0001; }

  gsl_rng_free(r);

  return;   
}

/* accept the model or not from the log likelihood 
The random number r should be [0, 1] random from uniform distribution 
we compare after applying log. 
accept: return 0
reject: return 1  
*/
int model_accept_tf(gsl_rng *r, LIKELI likeli0, LIKELI likeli1)
{
  double log_alpha, tf_alpha, tf_r = 0; /* variable name for true-false "r" */
  tf_r = gsl_rng_uniform(r);
  tf_r = log(tf_r);
  log_alpha = likeli1.pall - likeli0.pall; 

//  fprintf(stdout, "r %f alpha %f\n", tf_r, log_alpha);
  tf_alpha = 0;
  if( tf_alpha > log_alpha){ tf_alpha = log_alpha; }
  else{ return 0; }
 
  if( tf_alpha > tf_r){ return 0; }
  else{ return 1; }
}

/* accept the model or not -- optimized form for the parallel tempering 
The random number r should be drawn from [0, 1) uniform distribution.
We compare after applying logarithm.
The likelihood should be considered with temperature.
accept: return 0
reject: return 1
*/
int PT_model_accept_tf(gsl_rng *r, LIKELI likeli0, LIKELI likeli1, double invtemp, PARATEMP parallelt)
{ /* because we do this by each process, invtemp (single) input is enough */
  double log_alpha, tf_r; /* variable for the true-false r */
  tf_r = gsl_rng_uniform(r);
  tf_r = log(tf_r);

  log_alpha = likeli1.pall - likeli0.pall;
  /* adjust by temperature */
  if( parallelt.tf ){ log_alpha *= invtemp; }
  
  /* decide accept or reject */
  if( log_alpha > tf_r ){ return 0; } /* accept */
  else{ return 1;}

}

/* making new model from old model 
 * model has values in the order of H Vp Vs Rho -- each nlayer number of elements
 * first index (index 0) goes to the water thickness. 
 * nlayers_3: three times nlayer. 
 * nunknown: number of unknowns 
 * sigvec: vector with length of nlayers_n
 * mindex: vector with length of index.
 * the perturbation follows normal distribution with sigma sigx. 
 */
int newmodel_ocean(gsl_rng * r, PRIOR prior, int nunknown, int nlayers_3, double model0[], double *model1, int iflag)
{
  int i, index;
  double randx = 0;
  //time_t t;

  /* copy the model 0 */
  int nlayer = nlayers_3/3;
  for(i=0; i < nlayers_3 + nlayer; i++){ *(model1+i) = model0[i]; }

  // srand((unsigned) time(&t)); /* maybe this is not needed here? */
  index = gsl_rng_get(r) % nunknown; //rand() % nunknown;
  // fprintf(stderr, "%d --  %d\n", nunknown, index); 

  if(iflag == 0){ randx = gsl_ran_gaussian(r, prior.sigma[prior.mindex[index]]); }
  else{ randx = gsl_ran_gaussian(r, prior.sigma[prior.mindex[index]]*magicmul); }
  *(model1 + prior.mindex[index]) =  model0[prior.mindex[index]] + randx;

  if(prior.mindex[index] >=nlayer && prior.mindex[index] < nlayer*2 && !(prior.rho_tf[prior.mindex[index]-nlayer]) && prior.highlim[prior.mindex[index] + 2*nlayer] < 0)
  {  /* change rho */
     *(model1 + prior.mindex[index] + 2*nlayer) = vp_empirical_rho_solid(model1[prior.mindex[index]]);
     if(vp_empirical_rho_solid(model1[prior.mindex[index]]) <0){ return -1; }
  }
  if( prior.mindex[index] >=nlayer && prior.mindex[index] < nlayer*2 && prior.p2s[ prior.mindex[index]-nlayer ]  > 0)
  { /* change vs */
    *(model1 + prior.mindex[index] + nlayer) = model1[prior.mindex[index]]*prior.p2s[ prior.mindex[index]-nlayer ];
  }

  return prior.mindex[index];
}

/* change all changeable parameters at a time */
int newmodel_ocean_all(gsl_rng * r, PRIOR prior, int nlayers_3, double model0[], double *model1, int iflag)
{
  int i;
  double randx = 0;

  /* copy the model 0 */
  int nlayer = nlayers_3/3;
  for(i=0; i < nlayers_3 + nlayer; i++){ *(model1+i) = model0[i]; }

  for(i=0; i < prior.nsolve; i++)
  { 
    if(iflag == 0){ randx = gsl_ran_gaussian(r, prior.sigma[prior.mindex[i]]); }
    else{ randx = gsl_ran_gaussian(r, prior.sigma[prior.mindex[i]]*magicmul); }
    *(model1 + prior.mindex[i]) =  model0[prior.mindex[i]] + randx;

    if(prior.mindex[i] >=nlayer && prior.mindex[i] < nlayer*2 && !(prior.rho_tf[prior.mindex[i]-nlayer]) && prior.highlim[prior.mindex[i] + 2*nlayer] < 0)
    {  /* change rho */
       *(model1 + prior.mindex[i] + 2*nlayer) = vp_empirical_rho_solid(model1[prior.mindex[i]]);
       if(vp_empirical_rho_solid(model1[prior.mindex[i]]) <0){ return -1; }
    }
    if( prior.mindex[i] >=nlayer && prior.mindex[i] < nlayer*2 && prior.p2s[ prior.mindex[i]-nlayer ]  > 0)
    { /* change vs */
      *(model1 + prior.mindex[i] + nlayer) = model1[prior.mindex[i]]*prior.p2s[ prior.mindex[i]-nlayer ];
    }

  }

  return 0;
}

/* somewhat do additional modification */
int newmodel_ocean_supp(int index, PARAM param, PRIOR priorparam, double *model1)
{  
  int ilayer = 0;
  if(index >= param.nl && index < 2*param.nl)
  {
    ilayer = index - param.nl;
    if(!priorparam.s_tf[ilayer]) /* S is not being searched -- so we give Vp/Vs = 4*/
    { model1[ilayer + 2*param.nl] = model1[index]*0.25; }  
  }
  else if(index >= 2*param.nl && index < 3*param.nl)
  {
    ilayer = index - 2*param.nl;
    if(!priorparam.p_tf[ilayer])
    { model1[ilayer + param.nl] = model1[index]*4; /* P is not being searched -- so we give with Vp/Vs = 4 */
      if(!priorparam.rho_tf[ilayer] && priorparam.highlim[ilayer + 3*param.nl] < 0){ model1[ilayer + 3*param.nl] = vp_empirical_rho_solid(model1[ilayer + param.nl]); } /* modify the corresponding rho */
    }
  }

  return 0;
}

/* convert model vector into input model vectors -- case of ocean & general */
void modelvec2inputvecs_ocean(int nlayer, double model[], float *mh, float *mvp, float *mvs, float *mrho)
{
  int i, indbase = 0;

  /* thickness of layers */
  for(i=0; i < nlayer; i++){ *(mh + i) = model[i]; }
  *(mh+nlayer) = 100.; /* half space */
  /* vp */
  *(mvp + 0) = watervp; /* water */
  indbase = nlayer;
  for(i=0; i < nlayer; i++){ *(mvp + i+1) = model[i+indbase]; }
  /* vs */
  *(mvs + 0) = 0; /* water */
  indbase += nlayer;
  for(i=0; i < nlayer; i++){ *(mvs + i+1) = model[i+indbase]; }
  /* rho */
  *(mrho + 0) = waterrho; /* water density */
  indbase += nlayer;
  for(i=0; i < nlayer; i++)
  {
    *(mrho + i +1) = model[i + indbase]; 
  }

  return;
} 

void modelvec2inputvecs_general(PARAM param, MODEL *m)
{
  int i, nlay = 0;
  
  nlay = (*m).nl - param.iwater;
  /* can do all at a time */
  for(i = 0; i < param.iwater; i++)
  {
    (*m).h[i] = (*m).modvec[i];
    (*m).vp[i] = watervp;
    (*m).vs[i] = 0.;
    (*m).rho[i] = waterrho;
  }
  for(i = param.iwater; i < (*m).nl; i++)
  { // (*m).modvec
    (*m).h[i] = (*m).modvec[i]; 
    (*m).vp[i] = (*m).modvec[i - param.iwater + nlay];
    (*m).vs[i] = (*m).modvec[i - param.iwater + nlay*2];
    (*m).rho[i] = (*m).modvec[i - param.iwater + nlay*3];
  }
  /* halfspace h */
  (*m).h[ (*m).nl-1 ] = 100;

  return;
}

/* model yes or no -- check the availability based on the prior conditions 
 * Return 0: if true (accept)
 * Return 1: if false
 * prior conditions 
 * -- within a range that we set.
 */
int model_prior_probability_ind(int index, double model[], double lowlim[], double highlim[])
{
  int tf = 0;
  if(index < 0){ tf = 1; }
  else if( model[index] < lowlim[index] || model[index] > highlim[index]){ tf = 1; }

  return tf;
}

int model_prior_probability_all(PRIOR prior, double model[])
{
  int i, tmpind, tf;
  tf = 0;

  for(i=0; i < prior.nsolve; i++)
  {
    tmpind = prior.mindex[i];
    if( model[tmpind] < prior.lowlim[tmpind] || model[tmpind] > prior.highlim[tmpind] ){ tf = 1; return tf;}
  }

  return tf;
}

/* model prior probablity -- monotonic increase by depth: should do for vp and vs each */
int model_prior_monotonic_increase(int nlayer, float vec[])
{
  int i, tf = 0;
  i = 0;

  while(i < nlayer-1 && tf == 0 )
  {
    if(vec[i] > vec[i+1]){ tf = 1; }
    i++;
  }

  return tf;
}

int model_vp_and_vs(int nlayer, float vp[], float vs[])
{
  int i, tf;
  tf = 0;
  i = 0;
  while(i < nlayer && tf == 0) 
  {
    if(vp[i] <= vs[i]*kmin0){ tf = 1; }
    i++;
  }

  return tf;
}

/* do the prior monotonic increase and the Vp Vs check at a time */
int model_check_prior(MODEL m, PARAM param)
{
  int i, tf;
  tf = 0;
 
  i = 0;
  while(i < param.wnl-1 && tf == 0)
  {
    if(param.increase)
    {
      if( m.vp[i] > m.vp[i+1] ){ tf = 1; }
      else if( m.vs[i] > m.vs[i+1] ){ tf = 1; }
      else if( m.rho[i] > m.rho[i+1] ){ tf = 1;}
    }
    if( m.vp[i] <= m.vs[i]*kmin0){ tf = 1; }
    i++;
  }

  return tf;
}

/* merge getting new model and checking the prior prbability */
int model_getnew_check_prior(gsl_rng *r, PRIOR prior, PARAM param, MODEL m0, MODEL *m1)
{
  int i, tf, tftmp;

  tf = 0; tftmp = 0;
  if(param.changeall){
    tf = newmodel_ocean_all(r, prior, param.nl3, m0.modvec, (*m1).modvec, tftmp);
    tf += model_prior_probability_all(prior, (*m1).modvec);
    modelvec2inputvecs_general(param, m1);

    if(tf == 0 && model_check_prior( *m1, param) != 0){
      while( tf != 0 || model_check_prior( *m1, param) != 0){
        tf = newmodel_ocean_all(r, prior, param.nl3, m0.modvec, (*m1).modvec, tftmp);
        tf += model_prior_probability_all(prior, (*m1).modvec);
        modelvec2inputvecs_general(param, m1);
        tftmp++;
      }
    }
  }else{
    i = newmodel_ocean(r, prior, prior.nsolve, param.nl3, m0.modvec, (*m1).modvec, tftmp);
    tf = model_prior_probability_ind(i, (*m1).modvec, prior.lowlim, prior.highlim);
    modelvec2inputvecs_general(param, m1);

    if(tf == 0 && model_check_prior( *m1, param) != 0){
      while( tf != 0 || model_check_prior( *m1, param) != 0){
        i = newmodel_ocean(r, prior, prior.nsolve, param.nl3, m0.modvec, (*m1).modvec, tftmp);
        tf = model_prior_probability_ind(i, (*m1).modvec, prior.lowlim, prior.highlim);
        modelvec2inputvecs_general(param, m1);
        tftmp++;
      }
    }
  }

//  if(tftmp != 0){ fprintf(stderr, "hey %d added \n", tftmp); }
  return tf;
}

/* get density from the vp: we only cover from sediments to crust */
/* float vp_empirical_rho_solid(float vp)
{
  double rho = 0;
  double vptmp = vp*0.001;

  if( vptmp <= 1.50){ rho = 1.030; }
  else if( 1.50 < vptmp && vptmp <= 1.53){ rho = 14.8*vptmp - 21.014; }
  else if( 1.53 < vptmp && vptmp <= 2.00 ){ rho = 1.135*vptmp - 0.190; }
  else if( 2.00 < vptmp && vptmp <= 2.80 ){ rho = 0.917 + 0.744*vptmp - 0.08*vptmp*vptmp; }
  else if( 2.80 < vptmp ){ rho = 1.6612*vptmp - 0.4721*vptmp*vptmp + 0.0671*vptmp*vptmp*vptmp - 0.0043*vptmp*vptmp*vptmp*vptmp + 0.000106*vptmp*vptmp*vptmp*vptmp*vptmp; } 
  else{ return -1;  }
  
  rho *= 1000;
  return rho;
}
*/

/* the mud, clay, silt */
/* float vp_empirical_rho_solid(float vp)
{
  double rho = 0;
  double vptmp = vp*0.001;

  if( vptmp <= 1.50){ rho = 1.030; }
  else if( 1.50 < vptmp && vptmp <= 1.53){ rho = 14.8*vptmp - 21.014; }
  else if( 1.53 < vptmp && vptmp <= 2.00 ){ rho = 1.135*vptmp - 0.190; }
  else if( 2.00 < vptmp && vptmp <= 3.50 ){ rho = 0.917 + 0.744*vptmp - 0.08*vptmp*vptmp; }
  else if( 3.50 < vptmp && vptmp <= 7.00 ){ rho = 2.4372 + 0.0761*vptmp; }
  else if( 7.00 < vptmp ){ rho = 1.6612*vptmp - 0.4721*vptmp*vptmp + 0.0671*vptmp*vptmp*vptmp - 0.0043*vptmp*vptmp*vptmp*vptmp + 0.000106*vptmp*vptmp*vptmp*vptmp*vptmp; } 
  else{ return -1;  }

  rho *= 1000;
  return rho;
}
*/

/* the calcareous sediments */
float vp_empirical_rho_solid(float vp)
{
  double rho = 0;
  double vptmp = vp*0.001;

  if( vptmp <= 1.50){ rho = 1.030; }
  else if( 1.50 < vptmp && vptmp <= 1.53){ rho = 14.8*vptmp - 21.014; }
  else if( 1.53 < vptmp && vptmp <= 1.90 ){ rho = 1.135*vptmp - 0.190; }
  else if( 1.90 < vptmp && vptmp <= 3.50 ){ rho = 2.351 - 7.497*pow(vptmp, -4.656); }
  else if( 3.50 < vptmp && vptmp <= 7.00 ){ rho = 2.4372 + 0.0761*vptmp; }
  else if( 7.00 < vptmp ){ rho = 1.6612*vptmp - 0.4721*vptmp*vptmp + 0.0671*vptmp*vptmp*vptmp - 0.0043*vptmp*vptmp*vptmp*vptmp + 0.000106*vptmp*vptmp*vptmp*vptmp*vptmp; } 
  else{ return -1;  }

  rho *= 1000;
  return rho;
}




/* get density range from the prior range set for Vp */
void get_densityrange_from_vp(double lowvp, double highvp, double *lowrho, double *highrho)
{
  *lowrho = vp_empirical_rho_solid(lowvp);
  *highrho = vp_empirical_rho_solid(highvp);

  return;
}

/* get minimum and maximum range of the outputs */
int get_minmax_range(PARAM param, PRIOR prior, OUTPARAM *o_param)
{
  int i;
  double min, max;

  /* for vp */
  min = 1500; max = -100;
  for(i=0; i < param.nl; i++)
  {
    if( min > prior.lowlim[i+param.nl] ){ min = prior.lowlim[i+param.nl]; }
    if( max < prior.highlim[i+param.nl] ){ max = prior.highlim[i+param.nl]; }
  }
  if(min - (*o_param).dvp > 0){ min -= 2*(*o_param).dvp; }
  max += 4*(*o_param).dvp;
  (*o_param).vparr[0] = min;
  (*o_param).vparr[max_length-1] = max;
  (*o_param).nvp = ceil((max-min)/(*o_param).dvp);
  if( (*o_param).nvp > max_length){ fprintf(stderr, "Vp output not good. Too much elements top consider. %d\n", (*o_param).nvp); return -1; }

  /* get rho */
  (*o_param).rhoarr[max_length-1] = vp_empirical_rho_solid(max); 
  if( param.iwater == 0){ (*o_param).rhoarr[0] = vp_empirical_rho_solid(min); }
  else{ (*o_param).rhoarr[0] = waterrho - 10; } /* if there's water */

  /* addition if we are solving for rho */
  for(i=0; i < param.nl; i++)
  {
    if( prior.lowlim[i + param.nl + param.nl + param.nl] > 0 && prior.lowlim[i + param.nl + param.nl + param.nl] < (*o_param).rhoarr[0])
    {
      fprintf(stderr, "am I here?: %d\n", i + param.nl + param.nl + param.nl);
      (*o_param).rhoarr[0] = prior.lowlim[i + param.nl + param.nl + param.nl];
    }
    if( prior.highlim[i + param.nl + param.nl + param.nl] > 0 && prior.highlim[i + param.nl + param.nl + param.nl] > (*o_param).rhoarr[max_length-1])
    {
      (*o_param).rhoarr[max_length-1] = prior.highlim[i + param.nl + param.nl + param.nl];
    }
  }
  (*o_param).drho = ((*o_param).rhoarr[max_length-1]- (*o_param).rhoarr[0])/((*o_param).nvp-2);
  (*o_param).nrho = (*o_param).nvp;
  
  /* for vs */
  min = 10000000; max = -100;
  for(i=0; i < param.nl; i++)
  { 
    if( min > prior.lowlim[i+param.nl+param.nl] ){ min = prior.lowlim[i+param.nl+param.nl]; }
    if( max < prior.highlim[i+param.nl+param.nl] ){ max = prior.highlim[i+param.nl+param.nl]; }
  }
  (*o_param).vsarr[0] = min;
  if( min - (*o_param).dvs > 0){ (*o_param).vsarr[0] = min - 2*(*o_param).dvs;}
  max += (*o_param).dvs*4;
  (*o_param).vsarr[max_length-1] = max;
  (*o_param).nvs = ceil((max-min)/(*o_param).dvs);
  if( (*o_param).nvs > max_length){ fprintf(stderr, "Vs output not good. Too many elements to consider.: %d \n", (*o_param).nvs); return -1; }

  /* for h: from h0 to the maximum possible summation of thicknesses. */
  max = 0;
  for(i=0; i < param.wnl-1; i++)
  {
    max += prior.highlim[i];    
  }
  max += ((*o_param).dh + 100);
  (*o_param).harr[0] = (*o_param).h0;
  (*o_param).harr[max_length-1] = max;
  (*o_param).nh = ceil((max-(*o_param).h0)/(*o_param).dh );
  if( (*o_param).nh > max_length){ fprintf(stderr, "H output not good. Too many elements to consider. %d\n", (*o_param).nh); return -1; }

  /* now for Vp/Vs = kappa */
  (*o_param).karr[0] = kmin;
  (*o_param).karr[max_length-1] = kmax;
  (*o_param).nk = ceil( (kmax - kmin)/dk );
  if( (*o_param).nk > max_length){ fprintf(stderr, "Vp/Vs output not good. Too many elements to consider. %d\n", (*o_param).nh); return -1; }

  /* check whether everything is within the maximum possible array length: max_length */
  if( (*o_param).nh > max_length || (*o_param).nvp > max_length || (*o_param).nvs > max_length || (*o_param).nrho > max_length || (*o_param).nk > max_length )
  {
    fprintf(stderr, "some problem in the output format. too long. \n\
    should be less than %d but \n\
    nH: %d   \n\
    nVp: %d  \n\
    nVs: %d  \n\
    nrho: %d \n\
    nkappa: %d \n\
    terminating.\n", max_length, (*o_param).nh, (*o_param).nvp, (*o_param).nvs, (*o_param).nrho, (*o_param).nk);
    return -1; 
  }
 
  return 0;
}

/* prepare to print as probablity -- make an array */
int grid_arrray_get_length(double intv, double min, double max)
{
  int narray = 0;
  narray = ceil((max-min)/intv)+1;

  return narray; /* this should be used to allocate a memory space */
}

/* make array with given minimum and maximum range */
void grid_array(double intv, int ngrid, double min, double *array)
{
  int i = 0;
  
  for(i=0; i < ngrid; i++)
  {
    array[i] = min + i*intv;
  }

  return; 
}

/* organize the accepted sets of model into probability 
 * The count array index is defined by ...
 * nvar*i + j where 0<= i < nh and 0<=j<nvar 
 * make sure that the counts are pre-defined by 0.
 */
void models_2_addcount(int nlayer, float mh[], float mvar[], int nh, double hgrid[], double dh, int nvar, double vargrid[], double dvar, float *count)
{
  int ind0, ind1, varind, i, k;
  double h0, h1;
 
  h1 = 0; h0 = 0;
  for(k=0; k < nlayer; k++)
  {
    h0 = h1; 
    h1 += mh[k];    /* modify the range h1 value */
    if(k==0){ h0 = hgrid[0]; }
    else if(k==nlayer-1){ h1 = hgrid[nh-1]; } 

    if(h1 < hgrid[0] || hgrid[nh-1] <= h0 || vargrid[nvar-1] <= mvar[k] || mvar[k] < vargrid[0] ){ /* we do nothing. */  }
    else
    { 
      ind0 = floor((h0-hgrid[0])/dh); 
      ind1 = floor((h1-hgrid[0])/dh);
      if(ind1 >= nh){ ind1 = nh-1; }

      /* find the variable's index - because it's constant throughout the depth */
      if( !isnan(mvar[k]) && mvar[k] >= vargrid[0] && mvar[k] < vargrid[nvar-1])
      {
        varind = floor((mvar[k] - vargrid[0])/dvar);

        /* now fill in the counts */
        for(i=ind0; i < ind1; i++){ *(count + i*nvar + varind) += 1;}
      }
    }   
  }

  return; 
}

/* now add count to the time - synthetics amplitude array */
void synth_2_addcount(int npts, float synth[], int namp, double ampgrid[], double damp, float *count)
{
  int ampind, i;
  
  for(i=0; i < npts; i++)
  { 
    if( ampgrid[namp-1] <= synth[i] || synth[i] < ampgrid[0] ){ /* we do nothing. */  }
    else
    {
      /* find the amplitude's index - because it's constant throughout the depth */
      ampind = floor((synth[i] - ampgrid[0])/damp); 

      /* now fill in the count */ 
      *(count + i*namp + ampind) += 1;
    }
  }

  return;
}
/* now change the count array into the probability array */
void count_2_ppd_norm(int nh, int nvar, float *count)
{
  int i, j, tmpind;
  double inv, allcount;

  for(i=0; i < nh; i++)
  {
    tmpind = i*nvar;
    allcount = 0;
    for(j=0; j < nvar; j++){ allcount += count[tmpind + j]; }
    if(allcount >0){
      inv = 1./allcount;
      for(j=0; j < nvar; j++){ *(count + tmpind + j) *= inv; }
    }
  }

  return;
}

void count_2_ppd(int nh, int nvar, float *count, int nmodels)
{
  int i, j, tmpind;
  double inv;
  inv = 1./nmodels;

  for(i=0; i < nh; i++)
  {
    tmpind = i*nvar;
    for(j=0; j < nvar; j++){ *(count + tmpind + j) *= inv; }
  }

  return;
}

/* getting median of the probability density function */
void get_median_struct(int nh, int nvar, double vargrid[], float ppd[], double *model_med)
{
  int i, j;
  double accumul = 0;

  for(i=0; i < nh; i++)
  {
    accumul = 0;
    j = 0;
    while(j < nvar && accumul<0.5)
    {
      accumul += ppd[i*nvar + j];
      j++;
    }
    if(j == nvar && accumul < 0.5){ model_med[i] = NAN; }
    else{ model_med[i] = vargrid[j-1];}
  }

  return;
}

/* getting the maximum probability density structure */
void get_maximumdensity_struct(int nh, int nvar, double vargrid[], float ppd[],  double *model_maxppd)
{
  int i, j, maxj;
  double maxppd = 0;
 
  for(i=0; i < nh; i++)
  {
    maxj = -1;
    maxppd = 0;
    for(j=0; j < nvar; j++)
    {
      if( ppd[i*nvar + j] > maxppd){ maxj = j; maxppd = ppd[i*nvar + j];}
    }
    if(maxj >=0 && maxppd > 0){    model_maxppd[i] = vargrid[maxj];  }
    else{ model_maxppd[i] = NAN;}
  }

  return; 
}
