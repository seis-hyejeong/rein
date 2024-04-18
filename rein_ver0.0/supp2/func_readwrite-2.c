#include "func_readwrite_2.h"
#include <gsl/gsl_statistics_double.h>

void initialize_param(PARAM *param){
  (*param).nl = 0;
  (*param).wnl = 0;
  (*param).nl3 = 0;
  (*param).iwater = 1;
  (*param).increase = true;
  (*param).niter = 0;
  (*param).nburn = 0;
  (*param).nsave = 1;
  (*param).changeall = true;
  sprintf( (*param).priorfile, "f");

  return;
}

void initialize_paralleltempering(PARATEMP *parallelt){
  int i;

  (*parallelt).tf = false;
  (*parallelt).nprocs = 0;
  (*parallelt).ncool = 0;
  (*parallelt).rcool = -1.;
  (*parallelt).thigh = 1.;
 
  for(i = 0; i < maxprocs; i++)
  {
    (*parallelt).Ts[i] = 1.;
    (*parallelt).invTs[i] = 1.;
  }

  return;
}

void initialize_prior(PRIOR *priorparam){
  int i; 
  (*priorparam).nsolve = 0;
  for(i=0; i < n_max; i++)
  {
    (*priorparam).rho_tf[i] = false;
    (*priorparam).p_tf[i] = false;
    (*priorparam).s_tf[i] = false;
    (*priorparam).h_tf[i] = false;
    (*priorparam).p2s[i] = -1;
  }

  for(i=0; i < n_max*4; i++)
  {
    (*priorparam).mindex[i] = -1;
    (*priorparam).highlim[i] = -1.;
    (*priorparam).lowlim[i] = -1.;
    (*priorparam).sigma[i] = -1.; 
  }

  return;
}

void initialize_indata(INDATA *inparam)
{
  int i;
  (*inparam).ivelocity = false;
  (*inparam).tf = false;
  (*inparam).stack_tf = true;
  (*inparam).covar_tf = false;
  sprintf((*inparam).name, " ");
  sprintf((*inparam).stack, " ");
  sprintf((*inparam).error, " ");
  (*inparam).time[0] = 0.; (*inparam).time[1] = 0.;
  (*inparam).delta = 0.;
  (*inparam).tind0 = 0; 
  (*inparam).npts = 1;
  (*inparam).nrayp = 0;
  (*inparam).nfiles = 0;
  (*inparam).ffilter[0] = 0.;
  (*inparam).ffilter[1] = 0.;
  (*inparam).gauss = -1;
  (*inparam).w = -1;
  (*inparam).frac = -1.;
  (*inparam).waterlvl = -1.;
  (*inparam).nbin = 0;
  for(i=0; i < n_max; i++){
    (*inparam).nsnrs[i] = 0; 
    (*inparam).rayp[i] = 0.;
    (*inparam).snrs[i] = 0.;
    (*inparam).snrs_factor[i] = -1;
    (*inparam).startind[i] = 0;
  }

  return; 
}

void initialize_orgdata(DATAORG *orgdata)
{
  int i;
  (*orgdata).nbin = 1;
  (*orgdata).dp = (maxrayp - minrayp);

  for(i=0; i < n_max; i++){
    (*orgdata).arrp[i] = 0;
  }

  return;
}

void initialize_outparam(OUTPARAM *outparam){
  int i;

  (*outparam).dvs = 0.;
  (*outparam).dvp = 0.;
  (*outparam).drho = 0.;
  (*outparam).dh = 0.;
  (*outparam).h0 = 0.;
  (*outparam).nvs = 0;
  (*outparam).nvp = 0;
  (*outparam).nrho = 0;
  (*outparam).nh = 0;
  (*outparam).nk = 0;

  sprintf((*outparam).outname, " ");
  
  for(i=0; i < max_length; i++)
  {
    (*outparam).vsarr[i] = 0.;
    (*outparam).vparr[i] = 0.;
    (*outparam).rhoarr[i] = 0.;
    (*outparam).harr[i] = 0.;
    (*outparam).karr[i] = 0.;
  }  

  return;
}


/* take inputs */
int take_parameters_ver2(char infilename[strlng], char inmodel[strlng], PARAM *parameter, DATAORG *orgdata, PARATEMP *parallelt, INDATA *racfdata, INDATA *zacfdata, INDATA *pacfdata, INDATA *rwavedata, OUTPARAM *outparam)
{
  int i, j = 0;
  FILE *tmpfile;
  char aline[strlng], strtmp[strlng], content[strlng];

  if(access(infilename, F_OK) !=0 ){ fprintf(stderr, "not existing file!! terminating.\n"); return -1; }
  sprintf(inmodel, "f");
  /* initialize (*rwavedata).w to 1 */
  (*rwavedata).w = 1.;

  tmpfile = fopen(infilename, "r");
  while(fgets(aline, strlng, tmpfile)!=NULL)
  {
    j = 0;
    while(aline[j]==' ' && j < strlen(aline) ){ j++; } /* remove white spaces before the text */
    if( aline[j] == '\n' ){ } /* this was an empty line */
    else if(aline[j]== '#' ){} /* this is commented out line */
    else if( sscanf(&aline[j], "%s = %s", content, strtmp) == 2 )
    {
      // fprintf(stderr, "content %s and strtmp %s\n", content, strtmp);
      if(strcmp(content, "racf_in") ==0){
        if( strcmp(strtmp, "F")==0 || strcmp(strtmp, "f") == 0){ (*racfdata).tf = false;}
        else{ (*racfdata).tf = true; sscanf(strtmp, "%s", (*racfdata).name); }
      }
      else if(strcmp(content, "racf_vel") ==0){
        if( strcmp(strtmp, "T")==0 || strcmp(strtmp, "t") == 0){ (*racfdata).ivelocity = true; }
        else{ (*racfdata).ivelocity = false;  }
      }
      else if(strcmp(content, "zacf_in") ==0){
        if( strcmp(strtmp, "F")==0 || strcmp(strtmp, "f") == 0){ (*zacfdata).tf = false; }
        else{ (*zacfdata).tf = true; sscanf(strtmp, "%s", (*zacfdata).name);  }
      }
      else if(strcmp(content, "zacf_vel") ==0){
        if( strcmp(strtmp, "T")==0 || strcmp(strtmp, "t") == 0){ (*zacfdata).ivelocity = true; }
        else{ (*zacfdata).ivelocity = false;  }
      }
      else if(strcmp(content, "pacf_in") ==0){
        if( strcmp(strtmp, "F")==0 || strcmp(strtmp, "f") == 0){ (*pacfdata).tf = false; }
        else{ (*pacfdata).tf = true; sscanf(strtmp, "%s", (*pacfdata).name); (*pacfdata).ivelocity = true; }
      }
      else if(strcmp(content, "pacf_trac") ==0){
        if( strcmp(strtmp, "f")==0 || strcmp(strtmp, "F") == 0){ (*pacfdata).ivelocity = true; }
        else{ (*pacfdata).ivelocity = false;  }
      }
      else if(strcmp(content, "rstack_in") ==0){
        if( strcmp(strtmp, "F")==0 || strcmp(strtmp, "f") == 0){ (*rwavedata).tf = false; }
        else{ (*rwavedata).tf = true; sscanf(strtmp, "%s", (*rwavedata).name); }
      }
      else if(strcmp(content, "pintv") == 0){ 
        sscanf(strtmp, "%lf", &(*orgdata).dp); 
        (*orgdata).arrp[0] = minrayp; 
        i=1;
        while( i < n_max &&  (*orgdata).arrp[i-1] + (*orgdata).dp <= maxrayp)
        { 
          (*orgdata).arrp[i] = (*orgdata).arrp[i-1] + (*orgdata).dp; i++; 
        }
        (*orgdata).arrp[ i ] = (*orgdata).arrp[i-1]; /* for the last element, put it into the same number as the last eelement */
        (*orgdata).nbin = i;
      }
      else if(strcmp(content, "changeall") ==0){
        if( strcmp(strtmp, "F")==0 || strcmp(strtmp, "f") == 0){ (*parameter).changeall = false; }
        else if( strcmp(strtmp, "T")==0 || strcmp(strtmp, "t") == 0){ (*parameter).changeall = true; }
        else{ fprintf(stderr, "changeall not recognized. we keep the default value: change all at once\n"); }
      }
      else if(strcmp(content, "initial_vmod_in")==0){ strcpy(inmodel, strtmp); }
      else if(strcmp(content, "prior_vmod_in")==0){ strcpy( (*parameter).priorfile, strtmp); }  
      else if(strcmp(content, "water") == 0){ 
        sscanf(strtmp, "%d", &(*parameter).iwater); 
        if( (*parameter).iwater > 1){ fprintf(stderr, "water included or excluded not well done\n"); return -1; }
      } 
      else if(strcmp(content, "increase") == 0){
        if( strcmp(strtmp, "F")==0 || strcmp(strtmp, "f") == 0){ (*parameter).increase = false; }
        else if( strcmp(strtmp, "T")==0 || strcmp(strtmp, "t") == 0){ (*parameter).increase = true; }
        else{ fprintf(stderr, "monotonic increase set up not recognized. We use default value: increasing by depth\n"); }
      }
      else if(strcmp(content, "n_layer")==0){ sscanf(strtmp, "%d", &(*parameter).nl); (*parameter).nl3 = (*parameter).nl*3; }
      else if(strcmp(content, "n_burn_in")==0){ sscanf(strtmp, "%d", &(*parameter).nburn); }
      else if(strcmp(content, "n_iter")==0){ sscanf(strtmp, "%d", &(*parameter).niter); }
      else if(strcmp(content, "n_saving")==0){ sscanf(strtmp, "%d", &(*parameter).nsave); }
      /* temperature setting */
      else if(strcmp(content, "r_cool")==0){ sscanf(strtmp, "%lf", &(*parallelt).rcool);  
        if((*parallelt).rcool > 0 && (*parallelt).rcool <= 1.)
        { /* if we got the r_cool --> means we will do parallel tempering */
          (*parallelt).ncool = ceil((*parallelt).rcool*(*parallelt).nprocs); 
          if( (*parallelt).ncool > (*parallelt).nprocs){ (*parallelt).ncool = (*parallelt).nprocs; }
          else if( (*parallelt).ncool == 0){ (*parallelt).ncool = 1; }
          (*parallelt).tf = true; 
        }
      }
      else if(strcmp(content, "T_high")==0 || strcmp(content, "t_high")==0){ sscanf(strtmp, "%lf", &(*parallelt).thigh); } 
      else if(strcmp(content, "dvs")==0){ sscanf(strtmp, "%lf", &(*outparam).dvs); }
      else if(strcmp(content, "dvp")==0){ sscanf(strtmp, "%lf", &(*outparam).dvp); }
      else if(strcmp(content, "dh")==0){ sscanf(strtmp, "%lf", &(*outparam).dh); }
      else if(strcmp(content, "hprint0")==0){ sscanf(strtmp, "%lf", &(*outparam).h0); }
      else if(strcmp(content, "outrootname")==0){ strcpy((*outparam).outname, strtmp); }
      else{ fprintf(stderr, "not within the input parameters: %s\n", content); return -1; }
    }
    else{ fprintf(stderr, "input not well taken?: %s \n terminating...\n", infilename); return -1; }
    
  }  
  if(tmpfile!=NULL){ fclose(tmpfile); }

  /* set the default DATAORG if the choices were not given */ 
  if( (*orgdata).nbin == 1)
  {
    (*orgdata).arrp[0] = minrayp;
    (*orgdata).dp = defaultdp;
    i=1;
    while( i < n_max &&  (*orgdata).arrp[i-1] + defaultdp <= maxrayp)
    { 
      (*orgdata).arrp[i] = (*orgdata).arrp[i-1] + defaultdp; 
      i++;
    }
    (*orgdata).nbin = i;
  }

  /* set the wnl after all is read */
  (*parameter).wnl = (*parameter).nl + (*parameter).iwater;

  /* set density range */
//  if((*parameter).vp[0] > 0 && (*parameter).vp[1] > 0)
//  {
//    get_densityrange_from_vp((*parameter).vp[0], (*parameter).vp[1], &(*parameter).rho[0], &(*parameter).rho[1]);
//  }
 
  /* if the prior_vmod_in is false and initial vmod in is false, all parameter prior information should be given */
//  if( strcmp(inmodel, "f") ==0 && strcmp((*parameter).priorfile,"f")==0)
//  {
//    if((*parameter).sigwaterh < 0 || (*parameter).sigvs < 0 || (*parameter).sigh < 0 || (*parameter).sigvp < 0){ fprintf(stderr, "prior parameter set-up contradictory. \n"); return -1; }
//  }

  /* check approriateness of set up */
  if( (*parameter).niter ==0 || (*parameter).nsave ==0 || (*parameter).nsave > (*parameter).niter )
  {
    fprintf(stderr, "the input not taken for the iteration set-ups. : %s\n", infilename); return -1;
  }
 
  return 0;  
}

/* initialize temperature */
int get_parallel_temperatures(PARATEMP *parallelt, gsl_rng *r)
{
  double rr, logthigh;
  int i;

  logthigh = log( (*parallelt).thigh); 
  for(i=0; i < (*parallelt).ncool; i++)
  { 
    (*parallelt).Ts[i] = 1.; (*parallelt).invTs[i] = 1.; 
    (*parallelt).likelihood[i] = 0.;
  }
  for(i = (*parallelt).ncool; i < (*parallelt).nprocs-1; i++)
  { 
    rr = gsl_rng_uniform(r);
    rr *= logthigh;
    (*parallelt).Ts[i] = exp(rr); 
    (*parallelt).invTs[i] = 1./(*parallelt).Ts[i];
    (*parallelt).likelihood[i] = 0.;
  } 

  /* the highest temperature should follow thigh */
  i = (*parallelt).nprocs-1;
  (*parallelt).Ts[i] = (*parallelt).thigh;
  (*parallelt).invTs[i] = 1./(*parallelt).Ts[i];
  (*parallelt).likelihood[i] = 0.;

  return 0;
}

/* taking the input data -- ver 3 made in August 17th, 2022 */
int take_data_input_ver3(int my_id, INDATA *data, DATAORG orgdata, char outroot[strlng], float **stack, float **error, float **indvtrace)
{
  int i, j, iline, iname, tmpint, itmpfloat;
  char aline[strlng], listnames[n_max][strlng];
  int prep_ind0[n_max], prep_bincount[n_max], prep_fileindex[n_max], prep_inclbin[n_max];
  double tmpsnrs[n_max], tmpdouble;
  float tmpfloat, *tmptrace;
  SACHEAD hdr_s, hdr_e;
  FILE *tmpfile2;
   
  if(!(*data).tf){ return 0; } /*because it's false we don't use this */

  if(access((*data).name, F_OK) !=0 ){ fprintf(stderr, "%s: not existing file!! terminating.\n", (*data).name); return -1; }
  tmpfile2 = fopen((*data).name, "r");

  iline = 0; tmpfloat = -1; itmpfloat = -1;
  iname = 0; /* this is when the name of the files are not given individually */
  while(iline < 9 && fgets(aline, strlng, tmpfile2)!=NULL)
  {
    j = 0;
    while(aline[j]==' ' && j < strlen(aline) ){ j++; } /* remove white spaces before the text */

    if( aline[j] == '\n' ){ } /* this was an empty line */
    else if(aline[j]=='#' ){ } /* this is commented out line */
    else{
      if(iline ==0){ sscanf(&aline[j], "%d", &(*data).nrayp); }
      else if(iline==1){ sscanf(&aline[j], "%s", (*data).stack); }
      else if(iline==2){ 
        if(sscanf(&aline[j], "%s %f", (*data).error, &tmpfloat)==2){ itmpfloat = 1; (*data).waterlvl = tmpfloat; }
        else{ sscanf(&aline[j], "%s", (*data).error); }
 
        if( strcmp((*data).error, "C") == 0 || strcmp( (*data).error, "c") == 0){ (*data).covar_tf = true; }
      }
      else if(iline==3){
        if(aline[j] == 'f' || aline[j] == 'F'){ (*data).stack_tf = false; } /* we do not stack -- take individual input files */
        else if(aline[j] == 't' || aline[j] == 'T'){ (*data).stack_tf = true; } /* we stack and compare with the stack only */
        else{ fprintf(stderr, "error taking true/false of the stacking.\n"); return -1; }
      }
      else if(iline==4){ 
        if( sscanf( &aline[j], "%d %f", &(*data).nsnrs[0], &(*data).frac) ==2 ){   
          /* check whether the file is agreeing or not -- the nsnrs[0] must be 1 and the nrayp should not be 1*/
          if( (*data).nsnrs[0] != 1 || (*data).nrayp == 1 ){ fprintf(stderr, "error taking the nsnr, nrayp and fraction %s\n", (*data).name); return -1; }   
        }
        else{ sscanf(&aline[j], "%d", &(*data).nsnrs[0]); }  /* nsnr should be 1 if no noise is considered */
      }
      else if(iline==5){ sscanf(&aline[j], "%lf", &(*data).gauss); }
      else if(iline==6){ 
        if(aline[j]=='f' || aline[j] =='F'){ (*data).ffilter[0] = -1; (*data).ffilter[1] = -1; }
        else{
          sscanf(&aline[j], "%lf,%lf", &(*data).ffilter[0], &(*data).ffilter[1]);
        }
//        fprintf(stderr, "filter: %f %f\n", (*data).ffilter[0], (*data).ffilter[1]);
      }
      else if(iline==7){ sscanf(&aline[j], "%lf", &(*data).w); }
      else if(iline==8){ sscanf(&aline[j], "%lf,%lf", &(*data).time[0], &(*data).time[1]); }
      iline++;
    }
    /* end of reading one line */
  }

  /* check the consistency of the rayp and snr input */
  if( (*data).nrayp < 1){ if(my_id == 0){ fprintf(stderr, "wrong number of ray parameter. should be larger than 0. \n"); }return -1; }
  else if( (*data).nrayp > 1 && (*data).nsnrs[0] > 1){ fprintf(stderr, "wrong combination of n_rayp %d and n_snr %d\n", (*data).nrayp, (*data).nsnrs[0]); return -1;  }
  
  if( (*data).nrayp == 1 && (*data).nsnrs[0] == 1){
    if( (*data).stack_tf == false ){ (*data).stack_tf = true; } 
  }
 
  (*data).nfiles = (*data).nrayp*(*data).nsnrs[0];
  if( (*data).nfiles > n_max ){ if(my_id == 0){ fprintf(stderr, "too many file input (reading data).\n"); } return -1; }
 
  i = 0;
  while(fgets(aline, strlng, tmpfile2)!=NULL )
  {
    j = 0;
    while(aline[j]==' ' && j < strlen(aline) ){ j++; } /* remove white spaces before the text */

    if( aline[j] == '\n' ){ } /* this was an empty line */
    else if(aline[j]=='#' ){ } /* this is commented out line */
    else{
      if((*data).nrayp == 1){ 
        if(i == 0){ sscanf(&aline[j], "%lf", &(*data).rayp[i]); } /* take the ray parameter */
        /* take the n snrs */
        else if(i < (*data).nsnrs[0] +1){ 
          if( !(*data).covar_tf && (*data).stack_tf  ){ sscanf(&aline[j], "%lf", &(*data).snrs[i-1]); } /* only need to take the snr  -- snr should be negative if we are not considering noise. */
          else{ sscanf(&aline[j], "%lf %s", &(*data).snrs[i-1], listnames[i-1]); iname = 1; } /* take snr and the file name */
        }
      }else{
        if( !(*data).covar_tf && (*data).frac < 0 && (*data).stack_tf ){ sscanf(&aline[j], "%lf %lf", &(*data).rayp[i], &(*data).snrs[i]); } /* take only the rayp and the snr */
        else{ sscanf(&aline[j], "%lf %lf %s", &(*data).rayp[i], &(*data).snrs[i], listnames[i]); iname = 1; } /* take rayp, snr, and the file name */
      } 
      i++;
    }
  }

  /* close the file if it's done reading */
  if( tmpfile2 != NULL ){ fclose(tmpfile2); }

  if( (*data).covar_tf && (!(*data).stack_tf || iname != 1) ){ fprintf(stderr, "error taking the data input. When \"c\" is taken, %d we should be using the stacks for error calculation\n", iname); return -1;  }
  if( (*data).frac > 0 && (!(*data).stack_tf || iname != 1) ){ fprintf(stderr, "it is contradictory to make synthetics partially considered and the stack is wrong. %d\n", iname); return -1; }
  if((*data).nrayp == 1){ i--; }

  /* if the numbers taken do not equal the nfiles, that's an error */
  if(i != (*data).nfiles ){ if(my_id ==0){ fprintf(stderr, "the ray information and SNR information are not sufficiently given!!: %d and %d\n", i, (*data).nfiles); }  return -1; }

  if(iname == 1 && my_id == 0)
  {
    sprintf(aline, "%s.files.lst", outroot);
    tmpfile2 = fopen(aline, "w");
    for(i=0; i < (*data).nfiles; i++){ fprintf(tmpfile2, "ind %d file %s\n", i, listnames[i]); }
    fclose(tmpfile2);
  }

  /* check the appropriateness of the given snrs */
  for(i=0; i < (*data).nfiles; i++)
  {
    if( (*data).snrs[i] > 0 && (*data).snrs[i] < 1){ if(my_id == 0){ fprintf(stderr, "erroneous SNR. It should be larger than 1 to be considered.\n"); } return -1; }
  }  

  /* if the frac > 0, we should organize the rayp and the snrs */
  if( (*data).frac > 0 )
  {
//    fprintf(stderr, "Fraction is %f\n", (*data).frac);
    for( i = 0; i < n_max; i++){ 
      /* initialize mutliple arrays */
      prep_ind0[i] = 0;
      prep_bincount[i] = 0;
      prep_fileindex[i] = 0;
      prep_inclbin[i] = 0;
    }
    /* work on organizing */   
    j = organize_rayp_snrs(orgdata, data, outroot, prep_inclbin, prep_ind0, prep_bincount, prep_fileindex);
    if( j != 0){ fprintf(stderr, "error while organizing the information. terminate\n"); return -1; }

    /* print the rayp and the bin information */
    sprintf(aline, "%s.raypbin", outroot);
    tmpfile2 = fopen(aline, "w");
    /* organize the rayp and nsnr info */ 
    for(i=0; i < (*data).nbin; i++) /* nbin is number of bins considered */
    {
      tmpdouble = 0.5*(orgdata.arrp[ prep_inclbin[i] ] + orgdata.arrp[ prep_inclbin[i] + 1] );
      fprintf(tmpfile2, "rayp0 %.5e rayp1 %.5e center_rayp %.6e nfiles %d\n", orgdata.arrp[ prep_inclbin[i] ], orgdata.arrp[ prep_inclbin[i] + 1], tmpdouble, prep_bincount[i] );
      for(j=0; j < prep_bincount[i]; j++)
      { 
        tmpint = prep_ind0[i] + j ;
        fprintf(tmpfile2, "%d %.5e %.5e\n", j, (*data).rayp[ prep_fileindex[tmpint] ], (*data).snrs[ prep_fileindex[tmpint] ] ); 
        tmpsnrs[ tmpint ] = (*data).snrs[ prep_fileindex[tmpint] ];
      }
    }
 
    fprintf(tmpfile2, "fin");  
    fclose(tmpfile2);

    /* initialize the rayp and snrs */
    for(i=0; i < n_max; i++)
    { (*data).nsnrs[i] = 0;
      (*data).rayp[i] = 0.;
      (*data).snrs[i] = 0.;
    }

    /* replace the rayp and nsnr to the data array */
    (*data).nfiles = prep_bincount[ (*data).nbin - 1 ] + prep_ind0[ (*data).nbin - 1 ];
    (*data).nrayp = (*data).nbin;
    for(i = 0; i < (*data).nbin; i++)
    {
      (*data).rayp[i] = 0.5*( orgdata.arrp[ prep_inclbin[i] ] + orgdata.arrp[ prep_inclbin[i] + 1] ); /* center value of the bin */
      (*data).nsnrs[i] = prep_bincount[i]; 
    }
 
    for(i=0; i < (*data).nfiles; i++)
    {
      (*data).snrs[i] = tmpsnrs[i]; 
    } /* all copied */
  }
  else if( (*data).nrayp > 1 && (*data).nsnrs[0] == 1 )
  { /* fix in the nsnrs -- nsnr per ray parameter */
    for(i=0; i < (*data).nrayp; i++)
    {
      (*data).nsnrs[i] = 1; 
    }
  }

  /* now open the stack and error files. get ready for retun */
  /* check existance and read */
  if(itmpfloat < 0 && (*data).frac > 0){ tmpfloat = (float) atof((*data).error); (*data).waterlvl = tmpfloat; }
  /* step 1: stack */
  if( (*data).frac < 0 ){
    if( access((*data).stack, F_OK) !=0 ){ if(my_id == 0){ fprintf(stderr, "not existing input data %s file!! terminating.\n", (*data).stack); } return -1; }
    else{ *stack = rsac((*data).stack, &hdr_s); }
  }
  else{ *stack = rsac( listnames[0], &hdr_s); } /* temporary hdr_s */

  /* step 2: error file */
  if( !(*data).covar_tf && (*data).frac < 0 ){ /* not covariance, and frac < 0 */
    if(access((*data).error, F_OK) ==0 ){ *error = rsac((*data).error, &hdr_e); }
    else if( isdigit((*data).error[0]) !=0 ){ /* if constant value */
      *error = rsac((*data).stack, &hdr_e); 
      tmpfloat = (float) atof((*data).error);
      for(i=0; i < hdr_e.npts; i++){ (*error)[i] = tmpfloat; }
      if(my_id == 0){ fprintf(stderr, "taking the error as constant number %.5e\n", tmpfloat); }

      if( tmpfloat > 0 )
      {
        for(i=0; i < hdr_e.npts; i++)
        {
          if( (*error)[i] < tmpfloat ){ (*error)[i] = tmpfloat; }
        }
      }

    } 
    else{ if(my_id ==0){ fprintf(stderr, "not existing input error %s file!! terminating.\n", (*data).error); } return -1; }
  }
  else{ /* else -- only malloc*/
    *error = rsac( listnames[0], &hdr_e);
    for( i = 0; i < hdr_e.npts; i++){ (*error)[i] = 0.; } 
  }

  (*data).delta = hdr_s.delta;

  /* check major component equality */
  if( hdr_s.delta == hdr_e.delta && hdr_s.b == hdr_e.b && hdr_s.npts == hdr_e.npts){ }
  else{ if(my_id == 0){ fprintf(stderr, "filename %s erroneous input %s and %s\n", (*data).name, (*data).stack, (*data).error); } return -1; }

  /* get index related to given time */
  (*data).tind0 = floor(((*data).time[0] - hdr_s.b)/hdr_s.delta);
  (*data).npts = ceil(((*data).time[1] - (*data).time[0])/hdr_s.delta) +1;

  /* now get ready with the individual input taking */
  if( (*data).frac < 0 && ((*data).stack_tf == (*data).covar_tf) ) /* if == stack F cov F, stack T cov T */
  {
    /* allocate the memory space */
    *indvtrace = malloc( (*data).npts*sizeof(float)*(*data).nfiles );
    if( *indvtrace == NULL){ fprintf(stderr, "%03d there is error allocating memory during %s\n", my_id, (*data).name); return -1; }

    for(i=0; i < (*data).nfiles; i++)
    {
      tmptrace = rsac( listnames[i], &hdr_e); 
      if( hdr_s.delta == hdr_e.delta && hdr_s.b == hdr_e.b && hdr_s.npts == hdr_e.npts && tmptrace !=NULL ){ }
      else{ if(my_id == 0){ fprintf(stderr, "erroroneous input %s and %s --> %s\n", (*data).stack, (*data).error, listnames[i]); } return -1; }
      (*data).startind[i] = i*(*data).npts; 
      /* copy the trace */
      for(j=0; j < (*data).npts; j++){ (*indvtrace)[ i*(*data).npts + j ] = tmptrace[(*data).tind0 +j]; }
    }

    /* if covariance is true, it is important to calculate the error term */
    if( (*data).covar_tf ) /* it sort of overlap (waste of energy to calculate the stack again */
    {
      calculate_stack_stddev((*data).nfiles, (*data).npts, (*indvtrace), tmptrace, (*error) + (*data).tind0 );
    }
  }
  else if( (*data).frac > 0 ) /* when frac > 0, take only the files that we need.  */
  {  /* read only the necessary and get the only necessary */
    *indvtrace = malloc( (*data).npts*sizeof(float)*(*data).nfiles );
    
    for(i=0; i < (*data).nfiles; i++)
    {
      tmptrace = rsac( listnames[ prep_fileindex[i] ], &hdr_e);
      if( hdr_s.delta == hdr_e.delta && hdr_s.b == hdr_e.b && hdr_s.npts == hdr_e.npts && tmptrace !=NULL ){ }
      else{ if(my_id == 0){ fprintf(stderr, "erroroneous input %s file SAC --> %s\n", (*data).name, listnames[ prep_fileindex[i] ]); } return -1; }

      /* copy the trace to indvtrace */
      for(j=0; j < (*data).npts; j++){ (*indvtrace)[ i*(*data).npts + j] = tmptrace[(*data).tind0 +j]; } 
    }
    /* calculate new stack and the new standev */
    calculate_stack_stddev((*data).nfiles, (*data).npts, (*indvtrace), (*stack) + (*data).tind0, (*error) + (*data).tind0);
     
    //calculate_stack_stddev((*data).nfiles, (*data).npts, (*indvtrace), tmptrace, (*error) + (*data).tind0);
  }
  else{ *indvtrace = malloc( 1*sizeof(float) ); } /* do not return empty, but at least one datapoint. */

  /* water level to the error */
  if( (*data).waterlvl > 0 )
  {
    if(my_id==0){ fprintf(stderr, "using water level\n"); }
    for(i=0; i < (*data).npts; i++)
    {
      if((*error)[i + (*data).tind0 ] < (*data).waterlvl ){ (*error)[i + (*data).tind0 ] = (*data).waterlvl; }
    }
  }

  return 0;
}

int waterlevel_error(int my_id, INDATA *data, float *error)
{
  int i, j;
  float maxe = -1000;

  for(i=0; i < (*data).npts; i++)
  {
    if( error[i] > maxe ){ maxe = error[i]; }
  }
  maxe *= waterlevelfrac;

  j = 0;
  for(i=0; i < (*data).npts; i++)
  {
    if( error[i] < maxe ){ error[i] = maxe; j++; }
  }
  if( j!=0 ){ (*data).waterlvl = maxe; } 

  return j; 
}

/* with LU decomposition */
int generate_covariance_LU(int my_id, char outroot[strlng], INDATA data, float *indvtrace, float *stack, gsl_matrix *L, gsl_permutation *permute)
{
  size_t ii, jj;
  int i, j, s;
  float tmpfloat;
  double tmpdouble;
  gsl_matrix *diff;
  char aline[strlng];
  FILE *tmpfile2;

  /* calculate covariance -- start with allocating memories */
  /* first allocate memory */
  diff = gsl_matrix_calloc(data.npts, data.nfiles);
  if(diff ==NULL){ fprintf(stderr, "error in allocation!!!!!!!!!\n"); return 1; }
  
  /* calculate the covariance */
  /* now copy the read data into the matrix calculating the deviation */
  for(i=0; i < data.nfiles; i++)
  { 
    for(j=0; j < data.npts; j++)
    {
      tmpdouble = indvtrace[data.npts*i + j] - stack[ data.tind0 + j];
      gsl_matrix_set(diff, j, i, tmpdouble);
    }
  }
  
  /* calculate the covariance by matrix multiple */
  tmpdouble = 1./data.nfiles;
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, tmpdouble, diff, diff, 0., L);
  
  /* print the covariance information */
  if(my_id==0)
  { 
    sprintf(aline, "%s.covariance", outroot);
    tmpfile2 = fopen(aline, "wb");
    fwrite(&data.npts, sizeof(int), 1, tmpfile2);

    for(i=0; i < data.npts; i++)
    {
      tmpfloat = data.time[0] + data.delta*i;
      fwrite(&tmpfloat, sizeof(float), 1, tmpfile2);
    }

    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        tmpfloat = (float) gsl_matrix_get(L, ii, jj);
        fwrite(&tmpfloat, sizeof(float), 1, tmpfile2);
      }
    }
    fclose(tmpfile2);

    sprintf(aline, "%s.stddev", outroot);
    tmpfile2 = fopen(aline, "w");
    fprintf(tmpfile2, "indexx stack stddev_square\n");
    for(i=0; i < data.npts; i++)
    {
      tmpfloat = (float) gsl_matrix_get(L, i, i);
      fprintf(tmpfile2, "%d %.10e %.10e\n",  data.tind0 +i, stack[data.tind0+i], tmpfloat);
    }
    fclose(tmpfile2);

    sprintf(aline, "%s.cov", outroot);
    tmpfile2 = fopen(aline, "w");
    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        fprintf(tmpfile2, "%.10e ", gsl_matrix_get(L, ii, jj));
      }
      fprintf(tmpfile2, "\n");
    }
    fclose(tmpfile2);
  }
  
  gsl_matrix_free(diff);
  
  /* water level the values */
  jj = 0;
  if(data.waterlvl > 0)
  {
    tmpdouble = data.waterlvl;
    tmpdouble = tmpdouble*tmpdouble;
    for(ii=0; ii < data.npts; ii++)
    {
      if( gsl_matrix_get(L, ii, ii) < tmpdouble){ gsl_matrix_set(L, ii, ii, tmpdouble); jj++; }
    }
    if( jj != 0 && my_id == 0){ fprintf(stderr, "water level applied to covariance of %s\n", data.name);}
  }

  /* taper the covariance matrixt */
  tmpdouble = taper_covariance_model(data, L);
  if( i !=0 && my_id == 0){
    fprintf(stderr, "covariance matrix tapered to zero in off-diagonal %s, expA: %.5e\n", data.name, tmpdouble);

    sprintf(aline, "%s.cov-taper", outroot);
    tmpfile2 = fopen(aline, "w");
    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        fprintf(tmpfile2, "%.10e ", gsl_matrix_get(L, ii, jj));
      }
      fprintf(tmpfile2, "\n");
    }
    fclose(tmpfile2);

  }

  /* now LU decomposition */
  i = gsl_linalg_LU_decomp(L, permute, &s);
  if( i != 0){ if(my_id == 0){ fprintf(stderr, "error: singular matrix of %s\n", data.name); } return i; }

  if(my_id ==0)
  {
    sprintf(aline, "%s.LUpermutation", outroot);
    tmpfile2 = fopen(aline, "w");
    gsl_permutation_fprintf(tmpfile2, permute, " %u");
    fclose(tmpfile2);

    sprintf(aline, "%s.L", outroot);
    tmpfile2 = fopen(aline, "w");
    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        fprintf(tmpfile2, "%.10e ", gsl_matrix_get(L, ii, jj));
      }
      fprintf(tmpfile2, "\n");
    }
    fclose(tmpfile2);
  }

  /* change the diagonal -- give inverse */
  for(ii = 0; ii < data.npts; ii++)
  {
    tmpdouble = gsl_matrix_get(L, ii, ii);
    if(fabs(tmpdouble) > Ulimit){ tmpdouble = 1./tmpdouble; }
    else{ tmpdouble = 0; }
    gsl_matrix_set(L, ii, ii, tmpdouble);
  }

  return 0;
}

/* with LDLT decomposition */
int generate_covariance_ldlt(int my_id, char outroot[strlng], INDATA data, float *indvtrace, float *stack, gsl_matrix *L)
{
  size_t ii, jj;
  int i, j;
  float tmpfloat;
  double tmpdouble,diag[data.npts];
  gsl_matrix *diff;
  char aline[strlng];
  FILE *tmpfile2;

  /* calculate covariance -- start with allocating memories */
  /* first allocate memory */
  diff = gsl_matrix_calloc(data.npts, data.nfiles);
  if(diff ==NULL){ fprintf(stderr, "error in allocation!!!!!!!!!\n"); return 1; }

  /* calculate the covariance */
  /* now copy the read data into the matrix calculating the deviation */
  for(i=0; i < data.nfiles; i++)
  {
    for(j=0; j < data.npts; j++)
    {
      tmpdouble = indvtrace[data.npts*i + j] - stack[ data.tind0 + j];
      gsl_matrix_set(diff, j, i, tmpdouble);
    }
  }

  /* calculate the covariance by matrix multiple */
  tmpdouble = 1./data.nfiles;
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, tmpdouble, diff, diff, 0., L);

  /* print the covariance information */
  if(my_id==0)
  {
    sprintf(aline, "%s.covariance", outroot);
    tmpfile2 = fopen(aline, "wb");
    fwrite(&data.npts, sizeof(int), 1, tmpfile2);

    for(i=0; i < data.npts; i++)
    {
      tmpfloat = data.time[0] + data.delta*i;
      fwrite(&tmpfloat, sizeof(float), 1, tmpfile2);
    }
 
    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        tmpfloat = (float) gsl_matrix_get(L, ii, jj);
        fwrite(&tmpfloat, sizeof(float), 1, tmpfile2);
      }
    }
    fclose(tmpfile2);

    sprintf(aline, "%s.stddev", outroot);
    tmpfile2 = fopen(aline, "w");
    fprintf(tmpfile2, "indexx stack stddev_square\n");
    for(i=0; i < data.npts; i++)
    {
      tmpfloat = (float) gsl_matrix_get(L, i, i);
      fprintf(tmpfile2, "%d %.10e %.10e\n",  data.tind0 +i, stack[data.tind0+i], tmpfloat);
    }
    fclose(tmpfile2);

    sprintf(aline, "%s.cov", outroot);
    tmpfile2 = fopen(aline, "w");
    for(ii=0; ii < data.npts; ii++)
    {
      diag[ii] = gsl_matrix_get(L, ii, ii);
      for(jj=0; jj < data.npts; jj++)
      {
        fprintf(tmpfile2, "%.10e ", gsl_matrix_get(L, ii, jj));
      }
      fprintf(tmpfile2, "\n");
    }
    fclose(tmpfile2);
  }

  gsl_matrix_free(diff);

  /* water level the values */
  jj = 0;
  if(data.waterlvl > 0)
  {
    tmpdouble = data.waterlvl;
    tmpdouble = tmpdouble*tmpdouble;
    for(ii=0; ii < data.npts; ii++)
    {
      if( gsl_matrix_get(L, ii, ii) < tmpdouble){ gsl_matrix_set(L, ii, ii, tmpdouble); jj++; }
    }
    if( jj != 0 && my_id == 0){ fprintf(stderr, "water level applied to covariance of %s\n", data.name);}
  }

  /* taper the covariance matrixt */
  tmpdouble = taper_covariance_model(data, L);
  if( i !=0 && my_id == 0){
    fprintf(stderr, "covariance matrix tapered to zero in off-diagonal %s, expA: %.5e\n", data.name, tmpdouble);

    sprintf(aline, "%s.cov-taper", outroot);
    tmpfile2 = fopen(aline, "w");
    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        fprintf(tmpfile2, "%.10e ", gsl_matrix_get(L, ii, jj));
      }
      fprintf(tmpfile2, "\n");
    }
    fclose(tmpfile2);

  }

  /* now LU decomposition */
  i = gsl_linalg_cholesky_decomp1(L);
  if( i != 0){ if(my_id == 0){ fprintf(stderr, "error: singular matrix of %s\n", data.name); } return i; }
//  else{ fprintf(stderr, "error message is %d\n", i); }

  if(my_id ==0)
  {
    sprintf(aline, "%s.L", outroot);
    tmpfile2 = fopen(aline, "w");
    for(ii=0; ii < data.npts; ii++)
    {
      for(jj=0; jj < data.npts; jj++)
      {
        fprintf(tmpfile2, "%.10e ", gsl_matrix_get(L, ii, jj));
      }
      fprintf(tmpfile2, "\n");
    }
    fclose(tmpfile2);
  }

  return 0; 
}

/* taper the covariance matrix */
int taper_covariance(INDATA data, gsl_matrix *L)
{
  int nn = L->size1;
  size_t ii, jj, width_gauss, width_freq, width_final;
  double tmpdouble;

  width_gauss = 0;
  if( data.gauss > 0 )
  { 
    tmpdouble = sqrt(-log(gausslim))*data.gauss/M_PI; 
    width_gauss = ceil( 1/tmpdouble/data.delta ); 
  }
  width_freq = 0;
  if( data.ffilter[1] > 0 ){ width_freq = ceil( (1./data.ffilter[1])/data.delta ); }
 
  
  width_final = width_gauss;
  if( width_final < width_freq){ width_final = width_freq; }

  if(width_final < 2){ return 0; } /* do nothing */
  else
  {
    width_final = ceil(sqrt(2)*width_final);
    for(ii=0; ii < nn; ii++)
    {
      if( ii > width_final){ for(jj=0; jj < ii - width_final; jj++){ gsl_matrix_set(L, jj, ii, 0.); } /* taper to zero */ }
      if( ii + width_final < nn){ for(jj = ii + width_final; jj < nn; jj++){ gsl_matrix_set(L, jj, ii, 0.);  } /* taper to zero */ }
    }
  }
  
  return width_final;
}

double taper_covariance_model(INDATA data, gsl_matrix *L)
{
  int nn = L->size1;
  int ii, jj;
  double a_gauss, a_freq, a_final;
  double diag[nn], tmpmax, tmpdouble, tmpdouble2;

  a_gauss = -1; 
  if( data.gauss > 0 )
  {
    a_gauss = data.gauss; //*data.gauss;
  }
  a_freq = -1;
  if( data.ffilter[1] > 0 ){ a_freq = M_PI*data.ffilter[1]/sqrt(-log(0.5)); }


  a_final = 0;
  if( a_gauss > 0){ a_final = a_gauss; }
  if( a_freq > 0 && a_final > a_freq){ a_final = a_freq; }

  for(ii = 0; ii < nn; ii++){ diag[ii] = gsl_matrix_get(L, ii, ii); }
  tmpmax = gsl_stats_mean(diag, 1, nn);
  // gsl_matrix_max(L);
  if(a_final < 3*data.delta){ return 0; } /* do nothing -- there's risk of error!  */
  else
  {
    for(ii=0; ii < nn; ii++)
    {
      for(jj=ii; jj < nn; jj++)
      {
        tmpdouble = fabs( ii-jj)/sqrt(2)*data.delta;
          tmpdouble2 = tmpmax*exp(-tmpdouble*a_final*a_final);
          gsl_matrix_set(L, ii, jj, tmpdouble2);
          gsl_matrix_set(L, jj, ii, tmpdouble2);
      } 
    }

  }
 
  return a_final;
}


/* organize the in data -- rayp and snrs */
int organize_rayp_snrs(DATAORG orgdata, INDATA *data, char outroot[strlng], int *prep_inclbin, int *prep_ind0, int *prep_bincount, int *prep_fileindex) /* necessary arrays needed */
{
  int i, j, tmpint, tmpinclude, nfileinclude;
  int countind[n_max], indx[n_max], inclindx[n_max], prep_tmpcount[n_max];
  size_t sortind[orgdata.nbin];
  double tmpdouble, fraction[n_max]; 
  FILE *filestr;
  char outfilename[strlng];

  /* first, initialize the countind */
  for(i=0; i < n_max; i++){ fraction[i] = 0.; countind[i] = 0.; indx[i] = 0; inclindx[i] = 0; }

  /* make the index array and the array of median */
  for(i = 0; i < (*data).nfiles; i++)
  {
    tmpint = floor(( (*data).rayp[i] - orgdata.arrp[0] )/orgdata.dp);
    if( tmpint >= orgdata.nbin ){ tmpint = orgdata.nbin - 1; }
    else if(tmpint < 0){ tmpint = 0; fprintf(stderr, "here?!\n");  }
    indx[ i ] = tmpint; 
    countind[tmpint] +=1; 
  }

  /* now sort the countind */
  for(i = 0; i < orgdata.nbin; i++){ fraction[i] = (float)countind[i]/(*data).nfiles; }
  gsl_sort_index(sortind, fraction, 1, orgdata.nbin); 

  /* now get the number of bins to include */
  nfileinclude = 0;

  if( (*data).frac < 1 )
  { 
    tmpdouble = 0; i = orgdata.nbin - 1;
    while( i >= 0 && tmpdouble < (*data).frac )
    {
      tmpint = sortind[i];
      tmpdouble += fraction[ tmpint ];
      inclindx[ tmpint ] = 1; /* we include the index */
      nfileinclude += countind[tmpint]; /* number of rays included */
      i--; 
    }
    tmpinclude = orgdata.nbin-1 - i ;
    // fprintf(stderr, "including %d bins\n", tmpinclude); 
  }
  /* when the nbin is given instead of fraction */
  else
  {
    tmpinclude = (int) (*data).frac; /* transfer the fraction into number of including bins */
    tmpdouble = 0; /* calculate the fraction */
    for(i=0; i < tmpinclude; i++)
    {
      tmpint = sortind[orgdata.nbin - 1 - i];
      inclindx[ tmpint ] = 1; /* we include the index */
      tmpdouble += fraction[ tmpint ];
      nfileinclude += countind[tmpint]; /* number of rays inlucded */
    }
    (*data).frac = tmpdouble; /* give back the actual fraction */
  }
 
  (*data).nbin = tmpinclude;
 
  /* print the bin information */
  sprintf(outfilename, "%s.raypinfo", outroot);
  filestr = fopen(outfilename, "w");
  fprintf(filestr, "nfiles %d nbin %zu frac %.5e Nincluded %d includedNfiles %d\n", (*data).nfiles, orgdata.nbin, (*data).frac, tmpinclude, nfileinclude);

  /* print each bin and preparation at the same time */
  fprintf(filestr, "rayp0 rayp1 rayp_center count fraction\n"); 

  tmpint = 0; j = 0;
  for(i = 0; i < orgdata.nbin; i++)
  {
    tmpdouble = (orgdata.arrp[i] + orgdata.arrp[i+1])*0.5;
    if( inclindx[i] == 1 ) /* if it is an included bin -- print with YES */
    { 
      fprintf(filestr, "%.5e %.5e %.6e %03d %.5e YES\n", orgdata.arrp[i], orgdata.arrp[i+1], tmpdouble, countind[i], fraction[i]  );  
      if( j > orgdata.nbin ){ fprintf(stderr, "error in organizing multiple rayp into bins. terminating\n"); return -1; }
      prep_inclbin[j] = i;
      prep_bincount[j] = countind[i];
      prep_tmpcount[j] = 0;
      prep_ind0[j] = tmpint;
      tmpint += countind[i];
      j++; 
    }
    else /* if it is not included bin -- print with NO */
    {
      fprintf(filestr, "%.5e %.5e %.6e %03d %.5e NO\n", orgdata.arrp[i], orgdata.arrp[i+1], tmpdouble, countind[i], fraction[i]  );  
    }
  }  
  fclose(filestr);

  /* now get ready to organize the information */
  tmpint = 0; j = 0;
  for(i = 0; i < (*data).nfiles; i++) 
  {
    j = 0;
    while(j < tmpinclude && indx[i] != prep_inclbin[j]){ j++; }
    // indx[ i ] = tmpint; 
    if( j < tmpinclude && indx[i] == prep_inclbin[j] )
    { /* if the ray is included, give the ray index to the array */
      prep_fileindex[ prep_ind0[j] + prep_tmpcount[j] ] = i;
      prep_tmpcount[j] += 1;
    }
  }

  /* check whether all files are found */
  for(i= 0; i < tmpinclude; i++)
  {
    if(prep_tmpcount[i] != prep_bincount[i] ){ fprintf(stderr, "error in organizing rays. total numbers not fit. terminating\n"); return -1; }
  }

  return 0;  
}

/* get the snr factors corresponding to the synthetic parameters */
void calculate_snr_factors(INDATA *data, SYNTHPAR synthpar)
{
  int i;
  
  for(i=0; i < (*data).nfiles; i++)
  {
    if( (*data).snrs[i] > 1. ){ (*data).snrs_factor[i] = 1./(synthpar.nfft*((*data).snrs[i]*(*data).snrs[i] -1)); }
    /* default is negative -- so we do not need to repeat giving negative for the else condition */
  }

  return;
}

/* snr factor when we are considering the velocity */
void calculate_snr_factors_vel(INDATA *data, SYNTHPAR synthpar, double *wsquare)
{
  int i;
  double tmpdouble;
  tmpdouble = 0;
  for(i=0; i < synthpar.nfft; i++)
  {
    tmpdouble += wsquare[i];
  }

  for(i=0; i < (*data).nfiles; i++)
  {
    if( (*data).snrs[i] > 1. ){ (*data).snrs_factor[i] = 1./(tmpdouble*((*data).snrs[i]*(*data).snrs[i] -1)); }
    /* default is negative -- so we do not need to repeat giving negative for the else condition */
  }

  return;
}

/* initialize synthetic parameter */
int initialize_synth_ver2(SYNTHPAR *synthpar, INDATA racfdata, INDATA zacfdata, INDATA pacfdata, INDATA rwavedata)
{
  int nfft = 0;
  double delta, t0_rf, t0_acf;

  /* check equality of the deltas */
  delta = 0;
  if(racfdata.tf){ delta = racfdata.delta; }
  else if(zacfdata.tf){ delta = zacfdata.delta; }
  else if(rwavedata.tf){ delta = rwavedata.delta; }
  else if(pacfdata.tf){ delta = pacfdata.delta;  }
  else{ fprintf(stderr, "error!!!\n"); }

  nfft = get_nfft(TLENGTH, delta);
  if( racfdata.tf){
    if(fabs(racfdata.delta - delta) > 1.E-4){ fprintf(stderr, "R-ACF error % and %f\n", racfdata.delta, delta); return -1; }
  }
  if( zacfdata.tf){ 
    if(fabs(zacfdata.delta - delta) > 1.E-4){ fprintf(stderr, "Z-ACF error %f and %f\n", zacfdata.delta, delta); return -1; }
  }
  if( pacfdata.tf){ 
    if(fabs(pacfdata.delta - delta) > 1.E-4){ fprintf(stderr, "P-ACF error %f and %f\n", pacfdata.delta, delta); return -1; }
  }
  if( rwavedata.tf ){
    if(fabs(rwavedata.delta -delta) > 1.E-4){ fprintf(stderr, "R STACK error %f and %f\n", rwavedata.delta, delta); return -1; } 
  }
  set_parameters(nfft, delta, &t0_rf, &t0_acf, &(*synthpar).dw);
  get_zerotindex(t0_rf, t0_acf, delta, &(*synthpar).rwaveind0, &(*synthpar).racfind0);
  (*synthpar).zacfind0 = (*synthpar).racfind0;
  (*synthpar).pacfind0 = (*synthpar).racfind0;

  (*synthpar).nfft = nfft;
  if(racfdata.tf)
  {
    (*synthpar).racfind[0] = (t0_acf + racfdata.time[0])/racfdata.delta;
    (*synthpar).racfind[1] = (*synthpar).racfind[0] + racfdata.npts;
  }else{ (*synthpar).racfind[0] = 0; (*synthpar).racfind[1] = 0;}
  
  if(zacfdata.tf)
  {
    (*synthpar).zacfind[0] = (t0_acf + zacfdata.time[0])/zacfdata.delta;
    (*synthpar).zacfind[1] = (*synthpar).zacfind[0] + zacfdata.npts;
  }else{ (*synthpar).zacfind[0] = 0; (*synthpar).zacfind[1] =  0;}
 
  if(pacfdata.tf)
  { 
    (*synthpar).pacfind[0] = (t0_acf + pacfdata.time[0])/pacfdata.delta;
    (*synthpar).pacfind[1] = (*synthpar).pacfind[0] + pacfdata.npts;
  }else{ (*synthpar).pacfind[0] = 0; (*synthpar).pacfind[1] =  0;}  
 
  if(rwavedata.tf)
  {
    (*synthpar).rwaveind[0] = (t0_rf + rwavedata.time[0])/rwavedata.delta;
    (*synthpar).rwaveind[1] = (*synthpar).rwaveind[0] + rwavedata.npts;
  }else{ (*synthpar).rwaveind[0] = 0; (*synthpar).rwaveind[1] = 0; }

  return 0;
}

/* read manually given initial velocity model: get model vector. */
int read_inputvelmod(PARAM param, PRIOR prior, char velmodname[strlng], MODEL *m)
{
  int i, iwater, ntmp = 0;
  char aline[strlng];
  FILE *tmpfile;
 
  iwater = param.iwater;
  (*m).nl = param.nl; 
  /* if false */
  if(strlen(velmodname) ==1 && ( velmodname[0] == 'f' || velmodname[0] == 'F'))
  {
    first_model(param, prior, m);
    return 0;
  }
  else
  {
    /* if not */
    tmpfile = fopen(velmodname, "r");
    while(fgets(aline, strlng, tmpfile)!=NULL && ntmp < param.nl+1)
    {
      if(sscanf(aline, "%f %f %f %f", &(*m).h[ntmp], &(*m).vp[ntmp], &(*m).vs[ntmp], &(*m).rho[ntmp]) == 4)
      {
        ntmp++;
      }else if(sscanf(aline, "%f %f %f", &(*m).h[ntmp], &(*m).vp[ntmp], &(*m).vs[ntmp]) == 3)
      {
        (*m).rho[ntmp] = vp_empirical_rho_solid((*m).vp[ntmp]); 
        ntmp++;
      }
    }
    if(ntmp != param.wnl){ fprintf(stderr, "the initial model has some problem.\n"); return -1;}
    /* convert it into model vector */
    for(i=0; i < param.nl; i++)
    {
      (*m).modvec[i] = (*m).h[i];
      (*m).modvec[i+param.nl] = (*m).vp[i+ iwater];
      (*m).modvec[i+param.nl+param.nl] = (*m).vs[i + iwater];
      (*m).modvec[i+param.nl*3] = (*m).rho[i + iwater];
    }
    fclose(tmpfile);
  }  

  
  return 0;
}

/* we can set the prior model information by a text file. */
int read_inputpriormod_ver2(PARAM param, PRIOR *prior, MODEL *m)
{
  int j, iwater, indadd, iline;
  FILE *tmpfile;
  char aline[strlng];
  double tmp1, tmp2, tmp3;
  char *token, delim[] = " ";

  iwater = param.iwater;
  /* if false, take information from the PARAM */
  if(strcmp(param.priorfile,"f") ==0 )
  {
    fprintf(stderr, "prior file is false? wrong. terminating\n"); return -1;
  }
  else
  {
    iline = 0;
    /* check existence of file */
    if(access(param.priorfile, F_OK) !=0 ){ fprintf(stderr, "not existing prior model file!! terminating.\n"); return -1; } 
    /* open the file and start reading */ 
    tmpfile = fopen(param.priorfile, "r");
    indadd = 0;
    while(fgets(aline, strlng, tmpfile)!=NULL)
    {
      j = 0;
      while(aline[j]==' ' && j < strlen(aline) ){ j++; } /* remove white spaces before the text */
      if( aline[j] == '\n' ){ } /* this was an empty line */
      else if(aline[j]=='#' ){} /* this is commented out line */
      else if(iline < iwater)
      { 
        if(sscanf(&aline[j], "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
        {  
          (*prior).h_tf[0] = true;
          (*prior).highlim[0] = tmp2;
          (*prior).lowlim[0] = tmp1;
          (*prior).sigma[0] = tmp3;
          (*prior).mindex[indadd] = 0;
          indadd++;
        }
        else
        { sscanf(&aline[j], "%lf", &tmp1);
          (*m).modvec[0] = tmp1;
          (*prior).highlim[0]=tmp1;
          (*prior).lowlim[0] = tmp1;
        }
        iline++;
      }
      else if(iline < param.wnl - 1)
      { 
        /* now use strtok */
        token = strtok(&aline[j], delim);
        if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
        {
          (*prior).h_tf[iline] = true;
          (*prior).highlim[iline] = tmp2;
          (*prior).lowlim[iline] = tmp1;
          (*prior).sigma[iline] = tmp3;
          (*prior).mindex[indadd] = iline;
          indadd++;
        }
        else
        { sscanf(token, "%lf", &tmp1);
          (*m).modvec[iline] = tmp1;
          (*prior).highlim[iline] = tmp1;
          (*prior).lowlim[iline] = tmp1;
        }

        /* for Vp */
        token = strtok(NULL, delim);
        if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
        {
          (*prior).p_tf[iline - iwater] = true;
          (*prior).highlim[iline - iwater + param.nl] = tmp2;
          (*prior).lowlim[iline - iwater + param.nl] = tmp1;
          (*prior).sigma[iline - iwater + param.nl] = tmp3;
          (*prior).mindex[indadd] = iline - iwater + param.nl;
          indadd++;
        }
        else
        { sscanf(token, "%lf", &tmp1);
          (*m).modvec[iline - iwater + param.nl] = tmp1;
          (*m).modvec[iline - iwater + 3*param.nl] = vp_empirical_rho_solid(tmp1); 
          (*prior).highlim[iline - iwater + param.nl] = tmp1;
          (*prior).lowlim[iline - iwater + param.nl] = tmp1;
        }

        /* for Vs */
        token = strtok(NULL, delim);
        if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
        {
          (*prior).s_tf[iline - iwater] = true;
          (*prior).highlim[iline - iwater + param.nl+ param.nl] = tmp2;
          (*prior).lowlim[iline - iwater + param.nl + param.nl] = tmp1;
          (*prior).sigma[iline - iwater + param.nl + param.nl] = tmp3;
          (*prior).mindex[indadd] = iline - iwater + param.nl + param.nl;
          indadd++;
        } 
        else
        { sscanf(token, "%lf", &tmp1);
          if( tmp1 >= 50)
          {
            (*m).modvec[iline - iwater + param.nl + param.nl] = tmp1;
            (*prior).highlim[iline - iwater + param.nl+ param.nl] = tmp1;
            (*prior).lowlim[iline - iwater + param.nl + param.nl] = tmp1;
          }
          else
          { 
            (*prior).p2s[iline-1] = 1./tmp1; 
            (*prior).highlim[iline - iwater + param.nl+ param.nl] = (*prior).highlim[iline - iwater + param.nl]/tmp1;
            (*prior).lowlim[iline - iwater + param.nl + param.nl] = (*prior).lowlim[iline - iwater + param.nl]/tmp1;
          }
        }
        /* when I give density, token is not null */
        token = strtok(NULL, delim);
        if(token !=NULL)
        { 
          if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
          {  
            (*prior).rho_tf[iline- iwater] = true; /* true is 1 */ 
            (*prior).highlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp2;
            (*prior).lowlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp1;
            (*prior).sigma[iline - iwater + param.nl+ param.nl +param.nl] = tmp3;
            (*prior).mindex[indadd] = iline - iwater + param.nl + param.nl + param.nl;
            indadd++;
          }
          else
          { 
            if(sscanf(token, "%lf", &tmp1)==1)
            {  
              (*m).modvec[iline - iwater + param.nl + param.nl + param.nl] = tmp1;
              (*prior).highlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp1;
              (*prior).lowlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp1;
            }
          }
        }
        /* if density not given, we don't solve. */ 
        iline++;
      }
      else /* the last row */
      { /* only Vp and Vs and possibly Rho */
        /* now use strtok */
        token = strtok(&aline[j], delim);
        if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
        {
          (*prior).p_tf[iline - iwater] = true;
          (*prior).highlim[iline - iwater + param.nl] = tmp2;
          (*prior).lowlim[iline - iwater + param.nl] = tmp1;
          (*prior).sigma[iline - iwater + param.nl] = tmp3;
          (*prior).mindex[indadd] = iline - iwater + param.nl;
          indadd++;
        }
        else
        { sscanf(token, "%lf", &tmp1);
          (*m).modvec[iline - iwater + param.nl] = tmp1;
          (*m).modvec[iline - iwater + 3*param.nl] = vp_empirical_rho_solid(tmp1);      
          (*prior).highlim[iline - iwater + param.nl] = tmp1;
          (*prior).lowlim[iline - iwater + param.nl] = tmp1;
        }

        token = strtok(NULL, delim);
        if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
        {
          (*prior).s_tf[iline- iwater] = true;
          (*prior).highlim[iline - iwater + param.nl+ param.nl] = tmp2;
          (*prior).lowlim[iline - iwater + param.nl + param.nl] = tmp1;
          (*prior).sigma[iline - iwater + param.nl + param.nl] = tmp3;
          (*prior).mindex[indadd] = iline - iwater + param.nl + param.nl;
          indadd++;
        }
        else 
        { sscanf(token, "%lf", &tmp1);
          if(tmp1 >=50)
          {
            (*m).modvec[iline- iwater + param.nl + param.nl] = tmp1;
            (*prior).highlim[iline - iwater + param.nl+ param.nl] = tmp1;
            (*prior).lowlim[iline - iwater + param.nl + param.nl] = tmp1;
          }
          else
          {
            (*prior).p2s[iline-iwater] = 1./tmp1;
            (*prior).highlim[iline - iwater + param.nl+ param.nl] = (*prior).highlim[iline - iwater + param.nl]/tmp1;
            (*prior).lowlim[iline - iwater + param.nl + param.nl] = (*prior).lowlim[iline - iwater + param.nl]/tmp1;
          }
        }

        token = strtok(NULL, delim);
        if(token !=NULL) /* density taking possible */
        {
          if(sscanf(token, "%lf,%lf,%lf", &tmp1, &tmp2, &tmp3)==3)
          {
            (*prior).rho_tf[iline- iwater] = true; /* true is 1 */
            (*prior).highlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp2;
            (*prior).lowlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp1;
            (*prior).sigma[iline - iwater + param.nl+ param.nl +param.nl] = tmp3;
            (*prior).mindex[indadd] = iline - iwater + param.nl + param.nl + param.nl;
            indadd++;
          }
          else
          { 
            if(sscanf(token, "%lf", &tmp1)==1)
            {                  
              (*m).modvec[iline - iwater + param.nl + param.nl + param.nl] = tmp1;
              (*prior).highlim[iline -iwater + param.nl+ param.nl +param.nl] = tmp1;
              (*prior).lowlim[iline - iwater + param.nl+ param.nl +param.nl] = tmp1;
            }
          }
        }

        iline++;
      }
    }
    if(tmpfile!=NULL){ fclose(tmpfile); }
    /*final check whether sufficient information is given */
    if(iline!=param.nl+1 ){ return -1; } /* fail */
  }

  /* indadd equals the number of unknowns */
  (*prior).nsolve = indadd; 

  return 0;
}
/* format of the prior file
water_h min,water_h max,sigma
h min,hmax,sigma vpmin,vpmax,vpsig vsmin,vsmax,vssig
...

If you don't want to use it, only give a single number 
*/

/* for printing the models from the ppd -- e.g., median, maximum 1D array */
int print_ppdmodels2(int nh, double harray[], double vp[], double vs[], double rho[], double kappa[],  char filename[])
{
  int i;
  FILE *tmpfile;
  tmpfile = fopen(filename, "w");
  
  fprintf(tmpfile, "H[m] Vp[m/s] Vs[m/s] rho[kg/m^3] Vp/Vs[]\n");
  for(i=0; i < nh-1; i++)
  {
    fprintf(tmpfile, "%.5f %.5f %.5f %.5f %.5f\n", harray[i], vp[i], vs[i], rho[i], kappa[i]);
  }
  i = nh-1;
  fprintf(tmpfile, "%.5f %.5f %.5f %.5f %.5f", harray[i], vp[i], vs[i], rho[i], kappa[i]);
  fclose(tmpfile);
 
  return 0; 
}

/* change the array in the form of printing */
void change_array_to_printable(OUTPARAM *o_param)
{
  int i;
  double tmpdelta;

  /* for h */
  tmpdelta = (*o_param).dh*0.5;
  for(i=0; i < (*o_param).nh; i++)
  {
    (*o_param).harr[i] += tmpdelta;
  }
 
  /* for vp */
  tmpdelta = (*o_param).dvp*0.5;
  for(i=0; i < (*o_param).nvp; i++)
  {
    (*o_param).vparr[i] += tmpdelta;
  }

  /* for vs */
  tmpdelta = (*o_param).dvs*0.5;
  for(i=0; i < (*o_param).nvs; i++)
  {
    (*o_param).vsarr[i] += tmpdelta; 
  }

  /* for rho */
  tmpdelta = (*o_param).drho*0.5;
  for(i=0; i < (*o_param).nrho; i++)
  {
    (*o_param).rhoarr[i] += tmpdelta;
  }  
 
  return;
}

void change_array_to_printable_single(int npts, double *array, double delta)
{
  int i;
  double tmpdelta = delta*0.5;
  for(i=0; i < npts; i++)
  {
    *(array+i) = array[i] + tmpdelta;
  }

  return;
}

/* print array variables!! */
void print_array_variables2(OUTPARAM o_param)
{
  int i;
  FILE *tmpfile;
  char outfilename[strlng];

  /* print the array of variables */
  sprintf(outfilename, "%s.h.lst", o_param.outname);
  tmpfile = fopen(outfilename, "w");
  for(i=0; i < o_param.nh-1; i++)
  {
    fprintf(tmpfile, "%.5f\n", o_param.harr[i]);
  }
  fclose(tmpfile);

  sprintf(outfilename, "%s.Vp.lst", o_param.outname); /* Vp */
  tmpfile = fopen(outfilename, "w");
  for(i=0; i < o_param.nvp-1; i++)
  {
    fprintf(tmpfile, "%.5f\n", o_param.vparr[i]);
  }
  fclose(tmpfile);

  sprintf(outfilename, "%s.Vs.lst", o_param.outname); /* Vs */
  tmpfile = fopen(outfilename, "w");
  for(i=0; i < o_param.nvs-1; i++)
  {
    fprintf(tmpfile, "%.5f\n", o_param.vsarr[i]);
  }
  fclose(tmpfile);

  sprintf(outfilename, "%s.rho.lst", o_param.outname); /* rho */
  tmpfile = fopen(outfilename, "w");
  for(i=0; i < o_param.nrho-1; i++)
  {
    fprintf(tmpfile, "%.5f\n", o_param.rhoarr[i]);
  }
  fclose(tmpfile);

  sprintf(outfilename, "%s.kappa.lst", o_param.outname); /* Vp/Vs == kappa */
  tmpfile = fopen(outfilename, "w");
  for(i=0; i < o_param.nk-1; i++)
  {
    fprintf(tmpfile, "%.5f\n", o_param.karr[i]); 
  }
  fclose(tmpfile);
  return;
}

/* print the probability densityfunction of certain variable */
void print_pdf_variable(OUTPARAM o_param, char outfilename[], int nvar, double *vararray, float *pdf_var)
{
  int i, j;
  FILE *tmpfile;

  tmpfile = fopen(outfilename, "w");
  for(i=0; i < o_param.nh-1; i++)
  {
    for(j = 0; j < nvar-1; j++)
    {
      fprintf(tmpfile, "%.5f,%.5f,%.4e\n", o_param.harr[i], vararray[j], pdf_var[i*nvar + j]);
    }
  }
  fclose(tmpfile);
  
  return; 
}

/* print the probability density function of the synthetics -- into a binary file */
void print_pdf_synthetics(char outfilename[], int ntime, double t0, double dt, float intrace[], int namp, double amparr[], float pdf[]) 
{
  int i, j;
  FILE *tmpfile;
  float tmpfloat;

  tmpfile = fopen(outfilename, "wb");
  
  fwrite(&ntime, sizeof(int), 1, tmpfile); /* write time npts */
  for(i=0; i < ntime; i++)
  {
    tmpfloat = t0 + dt*i;
    fwrite(&tmpfloat, sizeof(float), 1, tmpfile);
  }  
  fwrite(intrace, sizeof(float), ntime, tmpfile);
 
  /* the amplitude information */
  j = namp -1;
  fwrite(&j, sizeof(int), 1, tmpfile);
  for(i=0; i < j; i++)
  {
    tmpfloat = (float) amparr[i];
    fwrite(&tmpfloat, sizeof(float), 1, tmpfile);  
  }
  /* now write the amplitude pdf */
  for(i=0; i < ntime; i++)
  { // i*namp + j --> can write j = 0 ~ (namp-1) at once;
    fwrite(&pdf[i*namp], sizeof(float), j, tmpfile);
  } 

  fclose(tmpfile); 

  return;
}

