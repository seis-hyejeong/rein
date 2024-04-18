#include <mpi.h>
#include <unistd.h>
#include <time.h>
     
#include "func_likelihood-2.h"
#include "func_models.h"
#include "func_readwrite_2.h"
#include "func_synthetic_iso3.h"
#include "parameters.h"
#include "parameters-2.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include "sac.h"
#define scale 3.

/* pre-allocated arrays */
gsl_complex comp_expw[maxlay*maxnfft*6], shift[maxnfft];
double w2array[maxnfft];
int MPI_ith_cal[maxprocs], MPI_i_accept[maxprocs], MPI_i_accept_burn[maxprocs], MPI_ireject[maxprocs], MPI_ireject_cal[maxprocs];
int MPI_T1_nsaved[maxprocs], MPI_Th_nsaved[maxprocs], MPI_i_accept1[maxprocs], MPI_i_accept_burn1[maxprocs];
float mkappa[n_max];

/* prototype of function */
int errmsg(char **);
void helpmsg(char **);
void helpmsg_param();
void helpmsg_input();
void helpmsg_model();
void helpmsg_output();
void verbosecheck(PARAM, PRIOR, INDATA, INDATA, INDATA, INDATA, SYNTHPAR, OUTPARAM);
void verbosecheck_file(FILE *, PARAM, PRIOR, INDATA, INDATA, INDATA, INDATA, SYNTHPAR, OUTPARAM);
void print_iprior(FILE *, PARAM, PRIOR);
 
float *inracfs, *inracfe, *inracfe_inv, *inracf_indv, *inzacfs, *inzacfe, *inzacfe_inv, *inzacf_indv, *inpacfs, *inpacfe, *inpacfe_inv, *inpacf_indv, *inrwaves, *inrwavee, *inrwavee_inv, *inrwave_indv;
float *indv_racf, *indv_zacf, *indv_pacf, *synth_racf, *synth_zacf, *synth_pacf, *synth_rwave, *synth_zwave;
clock_t start, end;

/* main function */
int main(int argc, char **argv)
{
  gsl_set_error_handler_off();
  MPI_Status stat;
  size_t ii;
  int fileexists[maxprocs];
  int num_procs, my_id, namelen, additer, tmpiter;
  char processor_name[strlng];
  int i, j, errn, errn2, ith_cal, i_accept, i_accept1, i_accept_burn, i_accept_burn1, ireject, ireject_cal, ith_cal_ver, ith_cal_save, nswap, nsaved_T1, nsaved_Th, iverbose, verbal_pern, tmpint;
  char infilename[strlng], inmodel[strlng], outfilename[strlng];
  double tmpdouble, *parray, *sarray, *rhoarray, *karray, *amparray;
  float tmpfloat, *pdf_vp, *pdf_vs, *pdf_rho, *pdf_k, *pdf_rwave, *pdf_zacf, *pdf_racf, *pdf_pacf, *tmpppd;
  gsl_matrix_complex *tmpE[maxlay], *tmpinvE[maxlay], *Jmatrix, *Jmatrix_w, *A[maxlay], *tmpE_w, *tmpinvE_w, *A_w, *mat22, *mat66, *J11;
  gsl_vector_complex *x, *in;
  gsl_permutation *p;
  gsl_matrix *Lracf, *Lzacf, *Lpacf, *Lrwave;
  gsl_vector *Vracf, *Vzacf, *Vpacf, *Vrwave, *Yracf, *Yzacf, *Ypacf, *Yrwave;

  PARAM i_param;
  PRIOR i_prior;
  INDATA in_racf, in_zacf, in_pacf, in_rwave;
  DATAORG org_data;
  SYNTHPAR s_param;
  OUTPARAM o_param;
  PARATEMP i_ptemp;

  LIKELI likelihood0, likelihood1;
  MODEL m0, m1;
  FILE *filelog, *filelogt, *filemod, *filemod1, *filemodh, *fileP, *fileP1, *filePh, *tmpfile; /* to save logs: filemod1 and fileP1 is for saving only the temperatture 1 */
  fileP1 = NULL;
  filePh = NULL; filemodh = NULL;

  start = clock(); 

  /* set up for the MPI */
  my_id = 0; num_procs = 0; nsaved_T1 = 0; nsaved_Th = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Get_processor_name(processor_name, &namelen);
  if(maxprocs < num_procs){ exit(-1); }
  /* print the first and last chain's processor locations */
  if( my_id == 0 ){ fprintf(stderr, "First %d process AT %s\n", my_id, processor_name); }
  else if( my_id == num_procs -1 ){ fprintf(stderr, "Last %d process AT %s\n", my_id, processor_name); }

  /* intialize */
  for(i=0; i < maxlay; i++)
  {
   tmpE[i] = gsl_matrix_complex_calloc(6,6);
   tmpinvE[i] = gsl_matrix_complex_calloc(6,6);
   A[i] = gsl_matrix_complex_calloc(6,6);
  }
  Jmatrix = gsl_matrix_complex_calloc(6,6);
  Jmatrix_w = gsl_matrix_complex_calloc(5,5);
  tmpE_w = gsl_matrix_complex_calloc(2,2);
  tmpinvE_w = gsl_matrix_complex_calloc(2,2);
  A_w = gsl_matrix_complex_calloc(2,2);
  mat22 = gsl_matrix_complex_calloc(2,2); /* allocated space to make things easy */
  mat66 = gsl_matrix_complex_calloc(6,6);
  J11 = gsl_matrix_complex_calloc(5, 5);
  x = gsl_vector_complex_calloc(5);
  in = gsl_vector_complex_calloc(5);
  p = gsl_permutation_alloc(5);

  /* start running */
  errn = 0; errn2 = 0; 
  iverbose = 0; verbal_pern = -100;
  /* initialize */
  initialize_param(&i_param);
  initialize_paralleltempering(&i_ptemp);
  s_param = synthparam_null;
  initialize_prior(&i_prior);
  initialize_indata(&in_racf);
  initialize_indata(&in_zacf);
  initialize_indata(&in_pacf);
  initialize_indata(&in_rwave);
  initialize_outparam(&o_param);
  likelihood0 = likelihood_null; 
  likelihood1 = likelihood_null;
  initialize_model(&m0);
  initialize_model(&m1);
  initialize_orgdata(&org_data);

  /* initialize the passing parameters before possibility of any error */
  ith_cal = 0; i_accept = 0; i_accept_burn = 0; ireject = 0; ireject_cal = 0;
  i_accept1 = 0; i_accept_burn1 = 0; 

  /* now take input */
  if(argc < 2){ if(my_id==0){ errmsg(argv); } goto MAIN_0002; }
  else if(argc==2 )
  {
    if( strcmp(argv[1], "h") == 0 || strcmp(argv[1], "-h") == 0){ if(my_id==0){ helpmsg(argv); } goto MAIN_0002; } 
    else if( strcmp(argv[1], "p") == 0 || strcmp(argv[1], "-p") == 0){ if(my_id==0){ helpmsg_param(); } goto MAIN_0002; }
    else if( strcmp(argv[1], "i") == 0 || strcmp(argv[1], "-i") == 0){ if(my_id==0){ helpmsg_input(); } goto MAIN_0002; }
    else if( strcmp(argv[1], "m") == 0 || strcmp(argv[1], "-m") == 0){ if(my_id==0){ helpmsg_model(); } goto MAIN_0002; }
    else if( strcmp(argv[1], "o") == 0 || strcmp(argv[1], "-o") == 0){ if(my_id==0){ helpmsg_output(); } goto MAIN_0002; }
  }
  /* if the verbose option is given */
  else if(argc==3)
  {
    iverbose = 1;
    if( 1==sscanf(&argv[2][2], "%d", &verbal_pern) )
    { 
      if(my_id==0){ fprintf(stderr, "printing every %d steps...\n", verbal_pern); }
    }
    else{ if(my_id==0){ fprintf(stderr, "numerical not recognized. error.\n"); } goto MAIN_0002;  }

    /* print the information of the amplitude saving array */
    if(my_id==0)
    { fprintf(stderr, "the synthetics' probability density function will be also saved. \n\
From %.4f to %.4f Interval %.4f Ngrid %d\n", amps0 + 0.5*damps, amps0 + (namps-0.5)*damps, damps, (namps-1));
    }
  }
  else
  {
    if(my_id==0){ fprintf(stderr, "too many inputs. Maybe an inappropriate input\n"); } goto MAIN_0002;
  }

  /* for random number generation */
  const gsl_rng_type *T;
  sleep(my_id+2); /* to make sure the gsl random is differently made */

  time_t t;
  srand((unsigned) time(&t));
  int s =rand();
  gsl_rng *r;
  gsl_rng_env_setup();
  T= gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, s);

  /* take input and start running */
  sscanf(argv[1], "%s", infilename);

  /* read the infile name and the initial model */
  i_ptemp.nprocs = num_procs;
  errn = take_parameters_ver2(infilename, inmodel, &i_param, &org_data, &i_ptemp, &in_racf, &in_zacf, &in_pacf, &in_rwave, &o_param);
  if(errn !=0){ fprintf(stderr, "terminating 1\n"); goto MAIN_0002; }

  /* print whether we have veloity ACF or displacement ACF */
  if(my_id == 0 && (in_racf.tf || in_zacf.tf) )
  {
    if(in_racf.ivelocity || in_zacf.ivelocity ){ fprintf(stderr, "we are working with velocity autocorrelation function\n"); }
    else{ fprintf(stderr, "we are working with displacement autocorrelation function\n");  }
  }
  if(my_id == 0 && in_pacf.tf) 
  {
    if(in_pacf.ivelocity){ fprintf(stderr, "we are working with pressure autocorrelation function\n"); }
    else{ fprintf(stderr, "we are working with traction autocorrelation function\n"); }
  }

  /* log file set */
  sprintf(outfilename, "%s.%03d.inversion.log", o_param.outname, my_id);
  filelog = fopen(outfilename, "w");
  sprintf(outfilename, "%s.%03d.T.log", o_param.outname, my_id);
  filelogt = fopen(outfilename, "wb");

  if(in_racf.tf )
  {
    if(in_racf.ivelocity){ fprintf(filelog, "we are working with R velocity autocorrelation function\n"); }
    else{ fprintf(filelog, "we are working with R displacement autocorrelation function\n");  }
  }
  if(in_zacf.tf)
  {
    if(in_zacf.ivelocity){ fprintf(filelog, "we are working with Z velocity autocorrelation function\n"); }
    else{ fprintf(filelog, "we are working with Z displacement autocorrelation function\n");  }
  }
  if(in_pacf.tf)
  {
    if(in_pacf.ivelocity){ fprintf(filelog, "we are working with Pressure autocorrelation function\n"); }
    else{ fprintf(filelog, "we are working with Traction autocorrelation function\n");  }
  }
 
  if(my_id==0){ fprintf(stderr, "total %d processes\n", i_ptemp.nprocs); }
  fprintf(filelog, "process %d/%d AT %s\n", my_id, i_ptemp.nprocs, processor_name);
  fprintf(filelog, "generator type: %s\n", gsl_rng_name (r));
  fprintf(filelog, "process %d seed = %d randnumber0 %f \n", my_id, s, gsl_rng_uniform (r));

  /* take the input data information */
  sprintf(outfilename, "%s.racf", o_param.outname);
  errn = take_data_input_ver3(my_id, &in_racf, org_data, outfilename, &inracfs, &inracfe, &inracf_indv);
  if( in_racf.covar_tf )
  {
    Lracf = gsl_matrix_calloc( in_racf.npts, in_racf.npts);
    if(Lracf == NULL ){ errn += 1; }
    i = generate_covariance_ldlt(my_id, outfilename, in_racf, inracf_indv, inracfs, Lracf);
    if( i !=0 )
    { in_racf.covar_tf = false; 
      if(my_id == 0)
      {
        fprintf(stderr, "ERROR: not using the covariance matrix for R-ACF %s\n", in_racf.name);
        sprintf(outfilename, "%s.racf.L", o_param.outname);
        remove(outfilename);
      }
      fprintf(filelog, "ERROR: not using the covariance matrix for R-ACF %s\n", in_racf.name);
      if( in_racf.waterlvl < 0)
      {
        if(waterlevel_error(my_id, &in_racf, inracfe + in_racf.tind0)!=0)
        {
          if(my_id== 0){ fprintf(stderr, "ERROR: covariance canceled, water-leveled for R-ACF %s (%.5e) \n", in_racf.name, in_racf.waterlvl); }
          fprintf(filelog, "ERROR: covariance canceled, water-leveled for R-ACF %s (%.5e) \n", in_racf.name, in_racf.waterlvl);
        }
      }
    }
  }else{ Lracf = gsl_matrix_calloc( 1, 1);   }

  sprintf(outfilename, "%s.zacf", o_param.outname);
  errn += take_data_input_ver3(my_id, &in_zacf, org_data, outfilename, &inzacfs, &inzacfe, &inzacf_indv);
  if( in_zacf.covar_tf )
  {
    Lzacf = gsl_matrix_calloc( in_zacf.npts, in_zacf.npts);
    if(Lzacf == NULL ){ errn += 1; }
    i = generate_covariance_ldlt(my_id, outfilename, in_zacf, inzacf_indv, inzacfs, Lzacf);
    if( i !=0 )
    { in_zacf.covar_tf = false; 
      if(my_id == 0)
      {
        fprintf(stderr, "ERROR: not using the covariance matrix for Z-ACF %s\n", in_zacf.name);
        sprintf(outfilename, "%s.zacf.L", o_param.outname);
        remove(outfilename);
      }
      fprintf(filelog, "ERROR: not using the covariance matrix for Z-ACF %s\n", in_zacf.name);
      if( in_zacf.waterlvl < 0)
      {
        if(waterlevel_error(my_id, &in_zacf, inzacfe + in_zacf.tind0)!=0)
        {
          if(my_id== 0){ fprintf(stderr, "ERROR: covariance canceled, water-leveled for Z-ACF %s (%.5e) \n", in_zacf.name, in_zacf.waterlvl); }
          fprintf(filelog, "ERROR: covariance canceled, water-leveled for Z-ACF %s (%.5e) \n", in_zacf.name, in_zacf.waterlvl);
        }
      }
    }
  }else{ Lzacf = gsl_matrix_calloc( 1, 1);    }

  sprintf(outfilename, "%s.pacf", o_param.outname);
  errn += take_data_input_ver3(my_id, &in_pacf, org_data, outfilename, &inpacfs, &inpacfe, &inpacf_indv);
  /* pacf */ 
  if( in_pacf.covar_tf )
  {
    Lpacf = gsl_matrix_calloc( in_pacf.npts, in_pacf.npts);
    if(Lpacf == NULL ){ errn += 1; }
    i = generate_covariance_ldlt(my_id, outfilename, in_pacf, inpacf_indv, inpacfs, Lpacf);
    if( i !=0 )
    { in_pacf.covar_tf = false; 
      if(my_id == 0)
      {
        fprintf(stderr, "ERROR: not using the covariance matrix for P-ACF %s\n", in_pacf.name);
        sprintf(outfilename, "%s.pacf.L", o_param.outname);
        remove(outfilename);
      }
      fprintf(filelog, "ERROR: not using the covariance matrix for P-ACF %s\n", in_pacf.name);
      if( in_pacf.waterlvl < 0)
      {
        if(waterlevel_error(my_id, &in_pacf, inpacfe + in_pacf.tind0)!=0)
        {
          if(my_id== 0){ fprintf(stderr, "ERROR: covariance canceled, water-leveled for P-ACF %s (%.5e) \n", in_pacf.name, in_pacf.waterlvl); }
          fprintf(filelog, "ERROR: covariance canceled, water-leveled for P-ACF %s (%.5e) \n", in_pacf.name, in_pacf.waterlvl);
        }
      }
    }

  }else{ Lpacf = gsl_matrix_calloc( 1, 1);  }

  sprintf(outfilename, "%s.rwave", o_param.outname);
  errn += take_data_input_ver3(my_id, &in_rwave, org_data, outfilename, &inrwaves, &inrwavee, &inrwave_indv);
  if( in_rwave.covar_tf )
  {
    Lrwave = gsl_matrix_calloc( in_rwave.npts, in_rwave.npts);
    if(Lrwave == NULL ){ errn += 1; }
    i = generate_covariance_ldlt(my_id, outfilename, in_rwave, inrwave_indv, inrwaves, Lrwave);
    if( i !=0 )
    { in_rwave.covar_tf = false; 
      if(my_id == 0)
      {
        fprintf(stderr, "ERROR: not using the covariance matrix for R-wave %s\n", in_rwave.name);
        sprintf(outfilename, "%s.rwave.L", o_param.outname);
        remove(outfilename);
      }
      fprintf(filelog, "ERROR: not using the covariance matrix for R-wave %s\n", in_rwave.name);
      if( in_rwave.waterlvl < 0)
      {
        if(waterlevel_error(my_id, &in_rwave, inrwavee + in_rwave.tind0)!=0)
        {
          if(my_id== 0){ fprintf(stderr, "ERROR: covariance canceled, water-leveled for R-wave %s (%.5e) \n", in_rwave.name, in_rwave.waterlvl); }
          fprintf(filelog, "ERROR: covariance canceled, water-leveled for R-wave %s (%.5e) \n", in_rwave.name, in_rwave.waterlvl);
        }
      }
    }
  }else{ Lrwave = gsl_matrix_calloc( 1, 1); }

  if(errn !=0){ if(my_id ==0){ fprintf(stderr, "terminating 2\n"); }  goto MAIN_0002; }

  /* now work on vectors */
  if( in_racf.covar_tf )
  {
    Vracf = gsl_vector_calloc(in_racf.npts); Yracf = gsl_vector_calloc(in_racf.npts);
    for(ii = 0; ii < in_racf.npts; ii++){ tmpdouble = 1./gsl_matrix_get(Lracf, ii, ii);  gsl_vector_set(Vracf, ii, tmpdouble); }
  }else
  {
    Vracf = gsl_vector_calloc(1); Yracf = gsl_vector_calloc(1);
  }
  if( in_pacf.covar_tf )
  {
    Vpacf = gsl_vector_calloc(in_pacf.npts); Ypacf = gsl_vector_calloc(in_pacf.npts);
    for(ii = 0; ii < in_pacf.npts; ii++){ tmpdouble = 1./gsl_matrix_get(Lpacf, ii, ii); gsl_vector_set(Vpacf, ii, tmpdouble); }
  }else
  {
    Vpacf = gsl_vector_calloc(1); Ypacf = gsl_vector_calloc(1);
  }
  if( in_zacf.covar_tf )
  {
    Vzacf = gsl_vector_calloc(in_zacf.npts); Yzacf = gsl_vector_calloc(in_zacf.npts);
    for(ii = 0; ii < in_zacf.npts; ii++){ tmpdouble = 1./gsl_matrix_get(Lzacf, ii, ii); gsl_vector_set(Vzacf, ii, tmpdouble); }
  }else
  {
    Vzacf = gsl_vector_calloc(1); Yzacf = gsl_vector_calloc(1);
  }
  if( in_rwave.covar_tf )
  {
    Vrwave = gsl_vector_calloc(in_rwave.npts); Yrwave = gsl_vector_calloc(in_rwave.npts);
    for(ii = 0; ii < in_rwave.npts; ii++){ tmpdouble = 1./gsl_matrix_get(Lrwave, ii, ii); gsl_vector_set(Vrwave, ii, tmpdouble); }
  }else
  {
    Vrwave = gsl_vector_calloc(1); Yrwave = gsl_vector_calloc(1);
  }

  /* print warning for rwave error */
  if(in_rwave.tf && in_rwave.nsnrs[0] > 1){
    if( my_id == 0)
    {
      fprintf(stderr, "Hey, the code does not consider error for the rstack. given Nsnr %d is invalid.\n", in_rwave.nsnrs[0]);
      for(i = 0; i < in_rwave.nfiles; i++){ fprintf(stderr, "invalid snr: %.5f\n", in_rwave.snrs[i]);  }
    }
    in_rwave.nsnrs[0] = 1; 
  }

  /* initialize synthetic information */
  errn = initialize_synth_ver2(&s_param, in_racf, in_zacf, in_pacf, in_rwave);
  if(errn !=0){ fprintf(stderr, "terminating 3\n"); goto MAIN_0002; }
  if( s_param.nfft > maxnfft){ fprintf(stderr, "too many nfft points. terminating\n"); goto MAIN_0002; }

  /* set the snr_factor in prior -- only for the R-ACF and the Z-ACF. */
  calculate_snr_factors(&in_racf, s_param);
  calculate_snr_factors(&in_zacf, s_param);
  calculate_snr_factors(&in_pacf, s_param);

  /*
  if(my_id ==0 && iverbose ==1)
  {
    fprintf(stderr, "the zacf nfiles %d\n", in_zacf.nfiles); 
    for(i=0; i < in_zacf.nfiles; i++){ fprintf(stderr, "%d th file %f snr %f scale\n", i, in_zacf.snrs[i], in_zacf.snrs_factor[i]); }
  }
  */
  /* normalize the in_rwave */
  if(in_rwave.tf)
  {
    tmpfloat = normalize_max_2_new(in_rwave.npts, inrwaves+in_rwave.tind0, 0, in_rwave.w, in_rwave.npts);
    for(i=0; i < in_rwave.npts; i++)
    {
      inrwavee[i + in_rwave.tind0] *= tmpfloat;
    }
  }

  /* set up a prior probability */
  errn = read_inputpriormod_ver2(i_param, &i_prior, &m0);
  if(errn !=0){ fprintf(stderr, "terminating 4\n"); goto MAIN_0002; }

  if(i_param.changeall){ fprintf(filelog, "change all parameters at a time\n"); }
  else{ fprintf(filelog, "change one parameter at a time\n"); }
  if(my_id == 0)
  { 
    if(i_param.changeall){ fprintf(stderr, "change all parameters at a time\n"); }
    else{ fprintf(stderr, "change one parameter at a time: %d parameters\n", i_prior.nsolve); }
  }

  /* make first velocity model */
  errn = read_inputvelmod(i_param, i_prior, inmodel, &m0);
  if(errn !=0){ fprintf(stderr, "error making initial velocity model.terminate 5\n"); goto MAIN_0002; }

  /* check all inputs if it's verbose */
  if(iverbose ==1 && my_id ==0)
  {  
    /* print the initial model */
    fprintf(stderr, "Thickness\tP-velocity\tS-velocity\tDensity\n");
    for(i=0; i < i_param.wnl; i++)
    {
      fprintf(stderr, "%.4f\t%.4f\t%.4f\t%.4f\n", m0.h[i], m0.vp[i], m0.vs[i], m0.rho[i]);
    }
  }
  /* save in the log too */
  fprintf(filelog, "Thickness\tP-velocity\tS-velocity\tDensity\n");
  for(i=0; i < i_param.wnl; i++)
  {
    fprintf(filelog, "%.4f\t%.4f\t%.4f\t%.4f\n", m0.h[i], m0.vp[i], m0.vs[i], m0.rho[i]);
  }

  /* now allocate memories for running */
  tmpint = 1;
  if(in_racf.tf){ 
    tmpint = s_param.nfft; 
  }
  synth_racf = calloc(tmpint, sizeof(float));

  tmpint = 1;
  if(in_zacf.tf){ 
    tmpint = s_param.nfft; 
  }
  synth_zacf = calloc(tmpint, sizeof(float));

  tmpint = 1;
  if(in_pacf.tf){
    tmpint = s_param.nfft;
  }
  synth_pacf = calloc(tmpint, sizeof(float));

  tmpint = 1;
  if(in_rwave.tf){ tmpint = s_param.nfft; }
  synth_rwave = calloc(tmpint, sizeof(float));
  synth_zwave = calloc(tmpint, sizeof(float));

  if(synth_racf == NULL || synth_zacf == NULL || synth_pacf == NULL || synth_rwave == NULL || synth_zwave == NULL ){ fprintf(stderr, "error in memory allocation 6\n"); goto MAIN_0002; }
 
  /* allocate memories for un-stacking */
  tmpint = 1;
  if(!in_racf.stack_tf){
    tmpint = in_racf.npts*in_racf.nfiles;
  }
  indv_racf = calloc(tmpint, sizeof(float));
  
  tmpint = 1;
  if(!in_zacf.stack_tf){
    tmpint = in_zacf.npts*in_zacf.nfiles;
  }
  indv_zacf = calloc(tmpint, sizeof(float));
 
  tmpint = 1;
  if( !in_pacf.stack_tf){
    tmpint = in_pacf.npts*in_pacf.nfiles;
  }
  indv_pacf = calloc(tmpint, sizeof(float));

  if(indv_racf == NULL || indv_zacf == NULL || indv_pacf == NULL){ fprintf(stderr, "error in memory allocation -- individual 6-2\n"); goto MAIN_0002; }

  /* get ready for shifting factor for acfs  and differentiation -- predefined*/
  
  /* prepare the time shift and the differentiation and/or the pressure */
  if( in_zacf.ivelocity || in_racf.ivelocity || in_pacf.ivelocity ){
    tmpdouble = 1.;
    if(in_racf.tf){ tmpdouble = in_racf.delta; }
    else if(in_zacf.tf){ tmpdouble = in_zacf.delta; }
    else if(in_pacf.tf){ tmpdouble = in_pacf.delta; }
    else
    { /* this is a contradictory situation */
      if(my_id == 0){ fprintf(stderr, "contradictory set up of differential setting. not good, but no termination\n"); } 
    }

    if(my_id == 0){ fprintf(stderr, "gotten t = %f (npts = %d) for the w2array\n", tmpdouble, s_param.nfft); }
    trace_differential_double(tmpdouble, s_param.nfft, w2array);  
  }

  if(in_racf.tf || in_zacf.tf || in_pacf.tf ){
    /* void trace_shift_factor(double tshift, double dt, double amp, int nfft, gsl_complex *shift) */
    if(in_racf.tf){ trace_shift_factor(s_param.nfft*in_racf.delta*0.5, in_racf.delta, 1., s_param.nfft, shift); }
    else if(in_zacf.tf){ trace_shift_factor(s_param.nfft*in_zacf.delta*0.5, in_zacf.delta, 1., s_param.nfft, shift); }
    else{ trace_shift_factor(s_param.nfft*in_pacf.delta*0.5, in_pacf.delta, 1., s_param.nfft, shift); }
  }

  /* initialize the output grids */ 
  errn = get_minmax_range(i_param, i_prior, &o_param);
  if(errn !=0){ fprintf(stderr, "error making output grid format. terminate 7\n"); exit(-1); }

  grid_array(o_param.dvp, o_param.nvp, o_param.vparr[0], o_param.vparr);
  grid_array(o_param.dvs, o_param.nvs, o_param.vsarr[0], o_param.vsarr);
  grid_array(o_param.dh, o_param.nh, o_param.harr[0], o_param.harr);
  grid_array(o_param.drho, o_param.nrho, o_param.rhoarr[0], o_param.rhoarr);
  grid_array(dk, o_param.nk, o_param.karr[0], o_param.karr);

  pdf_vp = calloc(o_param.nvp*o_param.nh, sizeof(float));
  pdf_vs = calloc(o_param.nvs*o_param.nh, sizeof(float));
  pdf_rho = calloc(o_param.nrho*o_param.nh, sizeof(float));
  pdf_k = calloc(o_param.nk*o_param.nh, sizeof(float));
  parray = calloc(o_param.nh, sizeof(double));
  sarray = calloc(o_param.nh, sizeof(double));
  rhoarray = calloc(o_param.nh, sizeof(double)); 
  karray = calloc(o_param.nh, sizeof(double));
 
  pdf_rwave = calloc(namps*in_rwave.npts, sizeof(float));
  pdf_zacf = calloc(namps*in_zacf.npts, sizeof(float));
  pdf_racf = calloc(namps*in_racf.npts, sizeof(float));
  pdf_pacf = calloc(namps*in_pacf.npts, sizeof(float));
  amparray = calloc(namps, sizeof(double));

  if(tmpint < o_param.nvp*o_param.nh){ tmpint = o_param.nvp*o_param.nh; }
  if(tmpint < o_param.nvs*o_param.nh){ tmpint = o_param.nvs*o_param.nh; }
  if(tmpint < o_param.nrho*o_param.nh){ tmpint = o_param.nrho*o_param.nh; }
  if(in_rwave.tf && tmpint < namps*in_rwave.npts){ tmpint = namps*in_rwave.npts; }
  if(in_zacf.tf && tmpint < namps*in_zacf.npts){ tmpint = namps*in_zacf.npts; }
  if(in_racf.tf && tmpint < namps*in_racf.npts){ tmpint = namps*in_racf.npts; }
  if(in_pacf.tf && tmpint < namps*in_pacf.npts){ tmpint = namps*in_pacf.npts; }
  tmpppd = calloc(4*i_param.wnl, sizeof(float));

  if(pdf_rwave == NULL || pdf_zacf == NULL || pdf_racf == NULL || pdf_pacf == NULL || amparray == NULL || tmpppd == NULL){ fprintf(stderr, "the synthetics amlitude saving failed. 7-2\n"); goto MAIN_0001;  }
  grid_array(damps, namps, amps0, amparray);
 
  if(iverbose ==1 && my_id ==0)
  { /* check all inputs */
    verbosecheck(i_param, i_prior, in_rwave, in_racf, in_zacf, in_pacf, s_param, o_param);
  }
  verbosecheck_file(filelog, i_param, i_prior, in_rwave, in_racf, in_zacf, in_pacf, s_param, o_param);
 
  if(pdf_vp == NULL || pdf_vs == NULL || pdf_rho == NULL || pdf_k == NULL ){ fprintf(stderr, "error in memory allocation 8\n"); goto MAIN_0001; }
  if(parray == NULL || sarray == NULL || rhoarray == NULL || karray == NULL ){ fprintf(stderr, "error in memory allocation 9\n"); goto MAIN_0001; }
  for(i=0; i < o_param.nvp*o_param.nh; i++){ pdf_vp[i] = 0.; }
  for(i=0; i < o_param.nvs*o_param.nh; i++){ pdf_vs[i] = 0.; }
  for(i=0; i < o_param.nrho*o_param.nh; i++){ pdf_rho[i] = 0.; }
  for(i=0; i < o_param.nk*o_param.nh; i++){ pdf_k[i] = 0.; }

  for(i=0; i < namps*in_rwave.npts; i++){ pdf_rwave[i] = 0; }
  for(i=0; i < namps*in_zacf.npts; i++){ pdf_zacf[i] = 0; } 
  for(i=0; i < namps*in_racf.npts; i++){ pdf_racf[i] = 0; }
  for(i=0; i < namps*in_pacf.npts; i++){ pdf_pacf[i] = 0; }

  /* get ready with the likelihood calculation -- assuming random gaussian noise -- ver. 0 2021/11/26 */
  if( in_racf.tf ){ inracfe_inv = calloc(in_racf.npts, sizeof(float)); }else{  inracfe_inv = calloc(1, sizeof(float)); }
  if( in_zacf.tf ){ inzacfe_inv = calloc(in_zacf.npts, sizeof(float)); }else{  inzacfe_inv = calloc(1, sizeof(float)); }
  if( in_pacf.tf ){ inpacfe_inv = calloc(in_pacf.npts, sizeof(float)); }else{  inpacfe_inv = calloc(1, sizeof(float)); }
  if( in_rwave.tf ){ inrwavee_inv = calloc(in_rwave.npts, sizeof(float)); }else{  inrwavee_inv = calloc(1, sizeof(float)); }
  if(inracfe_inv == NULL || inzacfe_inv == NULL || inpacfe_inv == NULL || inrwavee_inv == NULL){ fprintf(stderr, "error in memory allocation 10\n"); goto MAIN_0001; } 

  if( !in_racf.covar_tf ){ model_likelihood_log_sigmasquare(in_racf.npts, inracfe + in_racf.tind0, inracfe_inv); }
  if( !in_zacf.covar_tf ){ model_likelihood_log_sigmasquare(in_zacf.npts, inzacfe + in_zacf.tind0, inzacfe_inv); }
  if( !in_pacf.covar_tf ){ model_likelihood_log_sigmasquare(in_pacf.npts, inpacfe + in_pacf.tind0, inpacfe_inv); }
  if( !in_rwave.covar_tf ){ model_likelihood_log_sigmasquare(in_rwave.npts, inrwavee + in_rwave.tind0, inrwavee_inv); } 

  /* get ready for printing files */
  sprintf(outfilename, "%s.%03d.likelihood.txt", o_param.outname, my_id);
  fileP = fopen(outfilename, "wb");
  if(my_id == 0 && i_ptemp.tf)
  {
    sprintf(outfilename, "%s.T1_likelihood.txt", o_param.outname);
    fileP1 = fopen(outfilename, "wb");

    fwrite( &i_ptemp.ncool, sizeof(int), 1, fileP1);
    fwrite( &i_param.nburn, sizeof(int), 1, fileP1);
    tmpint = i_param.niter - i_param.nburn;
    fwrite( &tmpint, sizeof(int), 1, fileP1);  

    sprintf(outfilename, "%s.Th_likelihood.txt", o_param.outname);
    filePh = fopen(outfilename, "wb");

    tmpint = 1;
    fwrite( &tmpint, sizeof(int), 1, filePh);
    fwrite( &i_param.nburn, sizeof(int), 1, filePh);
    tmpint = i_param.niter - i_param.nburn;
    fwrite( &tmpint, sizeof(int), 1, filePh);
  }

  sprintf(outfilename, "%s.prior.txt", o_param.outname);
  tmpfile = fopen(outfilename, "w");
  print_iprior(tmpfile, i_param, i_prior);
  fclose(tmpfile); 

  sprintf(outfilename, "%s.%03d.models.txt", o_param.outname, my_id);
  filemod = fopen(outfilename, "wb");
  tmpint = (i_param.niter - i_param.nburn)/i_param.nsave;
  fwrite(&tmpint, sizeof(int), 1, filemod);
  fwrite(&i_param.wnl, sizeof(int), 1, filemod);

  /* preparing to print when T = 1 */
  sprintf(outfilename, "%s.%03d.models_T1.txt", o_param.outname, my_id);  
  filemod1 = fopen(outfilename, "wb");
  if(i_ptemp.tf)
  {
    sprintf(outfilename, "%s.%03d.models_Th.txt", o_param.outname, my_id);
    filemodh = fopen(outfilename, "wb");
  }

  /* set up for the Parallel Tempering */
  if(my_id==0 && i_ptemp.tf){
    if(i_ptemp.nprocs == 1){ i_ptemp.thigh = 1;}
    get_parallel_temperatures(&i_ptemp, r);
    if(iverbose ==1){ fprintf(stderr, "we are parallel tempering.\n"); }

    /* print the array of temperature */
    sprintf(outfilename, "%s.Tarray.txt", o_param.outname);
    tmpfile = fopen(outfilename, "w");
    for(i= 0; i < i_ptemp.nprocs; i++)
    {
      /* print the array of temperature */
      fprintf(tmpfile, "%.10f\n", i_ptemp.Ts[i]);
    }
    fclose(tmpfile);

    /* prepare to print the number of changes at each step */
    sprintf(outfilename, "%s.N_Tchange.txt", o_param.outname);
    tmpfile = fopen(outfilename, "w");
    fprintf(tmpfile, "Nstep N_swapped\n");
  }

  /* start the iteration with initialization */
  additer = 0;
  
  /* convert gaussian parameter for r stack */
  in_rwave.gauss = convert_gauss(in_rwave.gauss);
  /* get the first synthetics and likelihood */
  m0.nl = i_param.wnl; m1.nl = i_param.wnl;

  /* for the first case, give many chance as possible when there's calculation error */ 
  tmpiter = 0;
  while( calculate_synthetics_noise_ver6(i_param, synth_racf, indv_racf, in_racf, synth_zacf, indv_zacf, in_zacf, synth_pacf, indv_pacf, in_pacf, synth_rwave, synth_zwave, in_rwave, m0, s_param, tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, comp_expw, shift, w2array)!=0 )
  { /* get new model and repeat */
    while( read_inputvelmod(i_param, i_prior, inmodel, &m0) != 0){ }
    tmpiter++;
  }
  additer += tmpiter;

  if(tmpiter!=0){ 
    fprintf(filelog, "synth error code: %d. \n", errn2);
    fprintf(stderr, "%03d synth error! %d\n", my_id, tmpiter);
    for(i=0; i < i_param.wnl; i++)
      {
        fprintf(filelog, "%.4f\t%.4f\t%.4f\t%.4f\n", m0.h[i], m0.vp[i], m0.vs[i], m0.rho[i]);
      }
  }
  synthetics_likelihood_ver3_ldlt(&likelihood0, s_param, 
synth_racf, indv_racf, inracfs, inracf_indv, inracfe_inv, Lracf, Vracf, Yracf, in_racf, 
synth_zacf, indv_zacf, inzacfs, inzacf_indv, inzacfe_inv, Lzacf, Vzacf, Yzacf, in_zacf, 
synth_pacf, indv_pacf, inpacfs, inpacf_indv, inpacfe_inv, Lpacf, Vpacf, Ypacf, in_pacf,
synth_rwave, inrwaves, inrwavee_inv, Lrwave, Vrwave, Yrwave, in_rwave);    

  /* print the starting */
  fprintf(filelog, "initial log likelihood %.6f\n", likelihood0.pall);
  fprintf(stderr, "%03d initial log likelihood %.6f\n", my_id, likelihood0.pall);
  fwrite(&likelihood0.pall, sizeof(float), 1, fileP);

  /* send the likelihood to the i_ptemp -- first step */
  if(i_ptemp.tf ){ 
    if( my_id == 0 ){ i_ptemp.likelihood[0] = likelihood0.pall; }
    else{ MPI_Send(&likelihood0.pall, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD); }
  }

  if(my_id==0 && i_ptemp.tf)
  {
    for(i=1; i < i_ptemp.nprocs; i++)
    {
      MPI_Recv(&i_ptemp.likelihood[i], 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &stat);
    }
    /* write the likelihood of T = 1 only */ 
    errn = print_likelihood(ith_cal, fileP1, i_ptemp);
    if(errn != 0){ fprintf(stderr, "error in the algorithm of parallel tempering: saving likelihoods\n"); exit(-1); }

    errn = print_likelihood_T(ith_cal, filePh, i_ptemp, i_ptemp.thigh);
    if(errn != 0){ fprintf(stderr, "error in the algorithm of parallel tempering: saving likelihoods (high temp) \n"); exit(-1); }
  }

  /* set up for the Parallel Tempering */
  if(my_id==0 && i_ptemp.tf){
    for(i= 1; i < i_ptemp.nprocs; i++)
    {
      MPI_Send(&i_ptemp.Ts[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      MPI_Send(&i_ptemp.invTs[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
    }
  }

  /* now start the iteration */
  ith_cal_ver = 0; ith_cal_save = 0; 
  while(ith_cal < i_param.nburn)
  {
    if(i_ptemp.tf)
    {
      if(my_id == 0){ s_param.T = i_ptemp.Ts[0]; s_param.invT = i_ptemp.invTs[0]; }
      else{
        /* now initialize the temperature now */
        MPI_Recv(&s_param.T, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
        MPI_Recv(&s_param.invT, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
      }
    }

    /* print the temperature information to the log file */
    fwrite(&s_param.T, sizeof(double), 1, filelogt);

    /* initialize errn */
    errn = 0;
    /* get new model check the model conditions too */
    errn = model_getnew_check_prior(r, i_prior, i_param, m0, &m1);

    if(errn == 0)
    { /* model accepts its condition, so continue to synthetics. */
      /* If there's an error, give two more chances per step? */
      errn2 = calculate_synthetics_noise_ver6(i_param, synth_racf, indv_racf, in_racf, synth_zacf, indv_zacf, in_zacf, synth_pacf, indv_pacf, in_pacf, synth_rwave, synth_zwave, in_rwave, m1, s_param, tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, comp_expw, shift, w2array);

      synthetics_likelihood_ver3_ldlt(&likelihood1, s_param,
synth_racf, indv_racf, inracfs, inracf_indv, inracfe_inv, Lracf, Vracf, Yracf, in_racf,
synth_zacf, indv_zacf, inzacfs, inzacf_indv, inzacfe_inv, Lzacf, Vzacf, Yzacf, in_zacf,
synth_pacf, indv_pacf, inpacfs, inpacf_indv, inpacfe_inv, Lpacf, Vpacf, Ypacf, in_pacf,
synth_rwave, inrwaves, inrwavee_inv, Lrwave, Vrwave, Yrwave, in_rwave);

       
      if(errn2 != 0)
      { 
        fprintf(filelog, "synth error: %d step Not accepted\n", ith_cal);
        // if( iverbose == 1){ fprintf(stderr, "%03d synth error at %d\n", my_id, ith_cal); }
        
        for(i=0; i < i_param.wnl; i++)
        {
          fprintf(filelog, "%.4f\t%.4f\t%.4f\t%.4f\n", m1.h[i], m1.vp[i], m1.vs[i], m1.rho[i]);
        }
        ireject_cal++;
      }
      else
      {
        /* now compare likelihoods */
        errn =  PT_model_accept_tf(r, likelihood0, likelihood1, s_param.invT, i_ptemp);
        if(errn == 0)
        { /* accept the model and replace */
          likelihood0.pall = likelihood1.pall;
          replace_model(m1, &m0);
          i_accept_burn++;
          if(s_param.T == 1.){ i_accept_burn1++; }
        }
      }
    }else{ ireject++; }
  
    ith_cal++;
    ith_cal_ver++;
    fwrite(&likelihood0.pall, sizeof(float), 1, fileP); 
    if(likelihood1.pacf > 0 || likelihood1.zacf > 0 || likelihood1.racf > 0)
    {
      fprintf(stderr, "racf %f zacf %f pacf %f rwave %f total %f myid %d at %d\n", likelihood1.racf, likelihood1.zacf, likelihood1.pacf, likelihood1.rwave, likelihood1.pall, my_id, ith_cal);
      exit(-1);
    }
 
    if(ith_cal_ver == verbal_pern)
    { 
      ith_cal_ver = 0;
      fprintf(stderr, "%03d burn_in %d likelihood %.6f (T %.6f)\n", my_id, ith_cal, likelihood0.pall, s_param.T);  
      /* save model */
      fprintf(filelog, "Thickness\tP-velocity\tS-velocity\tDensity\n");
      for(i=0; i < i_param.wnl; i++)
      {
        fprintf(filelog, "%.4f\t%.4f\t%.4f\t%.4f\n", m0.h[i], m0.vp[i], m0.vs[i], m0.rho[i]);
      }
    }

    /* send the likelihood to the i_ptemp */
    if(i_ptemp.tf){ 
      if(my_id == 0){ i_ptemp.likelihood[0] = likelihood0.pall; }
      else{ MPI_Send(&likelihood0.pall, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);  }
    }

    /* swap temperature */
    if(my_id==0 && i_ptemp.tf)
    {
      for(i=1; i < i_ptemp.nprocs; i++)
      {
        MPI_Recv(&i_ptemp.likelihood[i], 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &stat);
      }
      /* write the likelihood of T = 1 only */
      errn = print_likelihood(ith_cal, fileP1, i_ptemp);
      if(errn != 0){ fprintf(stderr, "error in the algorithm of parallel tempering: saving likelihoods\n"); exit(-1); }

      errn = print_likelihood_T(ith_cal, filePh, i_ptemp, i_ptemp.thigh);
      if(errn != 0){ fprintf(stderr, "error in the algorithm of parallel tempering: saving likelihoods (high temp) \n"); exit(-1); }

      /* swap */
      nswap = swap_temperature(&i_ptemp, r);
      fprintf(tmpfile, "%d %d\n", ith_cal, nswap);

      for(i= 1; i < i_ptemp.nprocs; i++)
      {
        MPI_Send(&i_ptemp.Ts[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        MPI_Send(&i_ptemp.invTs[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      }
    }

  }

  ith_cal_save = 0; 
  while(ith_cal < i_param.niter)
  {
    if(i_ptemp.tf)
    {
      if(my_id == 0){ s_param.T = i_ptemp.Ts[0]; s_param.invT = i_ptemp.invTs[0]; }
      else{
        /* now initialize the temperature now */
        MPI_Recv(&s_param.T, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
        MPI_Recv(&s_param.invT, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
      }
    }

    /* print the temperature information to the log file */
    fwrite(&s_param.T, sizeof(double), 1, filelogt);

    /* initialize errn */
    errn = 0;

    /* here, we have to save each m0 at the end */
    /* get new model and check additional model conditions */
    errn = model_getnew_check_prior(r, i_prior, i_param, m0, &m1);
    
    if(errn == 0)
    { /* true, so continue to synthetics. */
      /* If there's an error, give two more chances per step? */
      errn2 = calculate_synthetics_noise_ver6(i_param, synth_racf, indv_racf, in_racf, synth_zacf, indv_zacf, in_zacf, synth_pacf, indv_pacf, in_pacf, synth_rwave, synth_zwave, in_rwave, m1, s_param, tmpE, tmpinvE, Jmatrix, Jmatrix_w, A, tmpE_w, tmpinvE_w, A_w, mat22, mat66, J11, x, in, p, comp_expw, shift, w2array);

      synthetics_likelihood_ver3_ldlt(&likelihood1, s_param,
synth_racf, indv_racf, inracfs, inracf_indv, inracfe_inv, Lracf, Vracf, Yracf, in_racf,
synth_zacf, indv_zacf, inzacfs, inzacf_indv, inzacfe_inv, Lzacf, Vzacf, Yzacf, in_zacf,
synth_pacf, indv_pacf, inpacfs, inpacf_indv, inpacfe_inv, Lpacf, Vpacf, Ypacf, in_pacf,
synth_rwave, inrwaves, inrwavee_inv, Lrwave, Vrwave, Yrwave, in_rwave);

      if(errn2 != 0)
      {  
        fprintf(filelog, "synth error: %d Not accepted\n", ith_cal);
        // if( iverbose == 1){ fprintf(stderr, "%03d synth error: %d\n", my_id, ith_cal); }

        for(i=0; i < i_param.wnl; i++)
        {
          fprintf(filelog, "%.4f\t%.4f\t%.4f\t%.4f\n", m1.h[i], m1.vp[i], m1.vs[i], m1.rho[i]);
        }
        ireject_cal++;
      }
      else
      {
        /* now compare likelihoods */
        errn =  PT_model_accept_tf(r, likelihood0, likelihood1, s_param.invT, i_ptemp);
        if(errn == 0)
        { /* accept the model and replace */
          likelihood0.pall = likelihood1.pall;
          replace_model(m1, &m0);
          i_accept++;
          if(s_param.T == 1.){ i_accept1++; }
        }
      }
    }else{ ireject++; }
 
    /* save the m0 into the grid -- location changed*/

    ith_cal++;
    ith_cal_ver++;
    ith_cal_save++;

    fwrite(&likelihood0.pall, sizeof(float), 1, fileP);
    if(likelihood1.pacf > 0 || likelihood1.zacf > 0 || likelihood1.racf > 0)
    {
      fprintf(stderr, "racf %f zacf %f pacf %f rwave %f total %f myid %d at %d\n", likelihood1.racf, likelihood1.zacf, likelihood1.pacf, likelihood1.rwave, likelihood1.pall, my_id, ith_cal);
      exit(-1);
    }

    if(ith_cal_ver == verbal_pern)
    {
      ith_cal_ver = 0;
      fprintf(stderr, "%03d %d likelihood %.6f (T %.6f)\n", my_id, ith_cal, likelihood0.pall, s_param.T);

      /* save model */
      fprintf(filelog, "Thickness\tP-velocity\tS-velocity\tDensity\n");
      for(i=0; i < i_param.wnl; i++)
      {
        fprintf(filelog, "%.4f\t%.4f\t%.4f\t%.4f\n", m0.h[i], m0.vp[i], m0.vs[i], m0.rho[i]);
      }
    }
    /* print the m0 into the binary file */
    if(ith_cal_save== i_param.nsave)
    {
      ith_cal_save = 0;
      tmpfloat = (float) s_param.T;
      fwrite(&tmpfloat, sizeof(float), 1, filemod); 
      fwrite(m0.h, sizeof(float), i_param.wnl, filemod);
      fwrite(m0.vp, sizeof(float), i_param.wnl, filemod);
      fwrite(m0.vs, sizeof(float), i_param.wnl, filemod);
      fwrite(m0.rho, sizeof(float), i_param.wnl, filemod);
     
      if(s_param.T == 1.)
      { 
        nsaved_T1++;
        /* save the m0 into the grid if and only if T == 1 */
        fwrite(m0.h, sizeof(float), i_param.wnl, filemod1);
        fwrite(m0.vp, sizeof(float), i_param.wnl, filemod1);
        fwrite(m0.vs, sizeof(float), i_param.wnl, filemod1);
        fwrite(m0.rho, sizeof(float), i_param.wnl, filemod1);

        for(i=0; i < i_param.wnl; i++){ 
          if(m0.vs[i] == 0){mkappa[i] = NAN;}
          else{ mkappa[i] = m0.vp[i]/m0.vs[i]; }
        }
        /* save Vp */ models_2_addcount(i_param.wnl, m0.h, m0.vp, o_param.nh, o_param.harr, o_param.dh, o_param.nvp, o_param.vparr, o_param.dvp, pdf_vp);
        /* save Vs */ models_2_addcount(i_param.wnl, m0.h, m0.vs, o_param.nh, o_param.harr, o_param.dh, o_param.nvs, o_param.vsarr, o_param.dvs, pdf_vs);
        /* save rho */ models_2_addcount(i_param.wnl, m0.h, m0.rho, o_param.nh, o_param.harr, o_param.dh, o_param.nrho, o_param.rhoarr, o_param.drho, pdf_rho);
        /* save kappa */ models_2_addcount(i_param.wnl, m0.h, mkappa, o_param.nh, o_param.harr, o_param.dh, o_param.nk, o_param.karr, dk, pdf_k);
        /* save Rstack */ if(in_rwave.tf){ synth_2_addcount(in_rwave.npts, synth_rwave + s_param.rwaveind[0], namps, amparray, damps, pdf_rwave); }
        /* save Z-ACF */ if(in_zacf.tf){ synth_2_addcount(in_zacf.npts, synth_zacf + s_param.zacfind[0], namps, amparray, damps, pdf_zacf); }
        /* save R-ACF */ if(in_racf.tf){ synth_2_addcount(in_racf.npts, synth_racf + s_param.racfind[0], namps, amparray, damps, pdf_racf); }
        /* save P-ACF */ if(in_pacf.tf){ synth_2_addcount(in_pacf.npts, synth_pacf + s_param.pacfind[0], namps, amparray, damps, pdf_pacf); }
      }
      else if(i_ptemp.tf && s_param.T == i_ptemp.thigh)
      {
        nsaved_Th++;
        fwrite(m0.h, sizeof(float), i_param.wnl, filemodh);
        fwrite(m0.vp, sizeof(float), i_param.wnl, filemodh);
        fwrite(m0.vs, sizeof(float), i_param.wnl, filemodh);
        fwrite(m0.rho, sizeof(float), i_param.wnl, filemodh);
      }

    } 

    /* send the likelihood to the i_ptemp */
    if(i_ptemp.tf){
      if(my_id == 0){ i_ptemp.likelihood[0] = likelihood0.pall; }
      else{ MPI_Send(&likelihood0.pall, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD); }
    }

    /* swap temperature */
    if(my_id==0 && i_ptemp.tf){
      for(i=1; i < i_ptemp.nprocs; i++)
      { 
        MPI_Recv(&i_ptemp.likelihood[i], 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &stat);
      }
      /* write the likelihood of T = 1 only */
      errn = print_likelihood(ith_cal, fileP1, i_ptemp);
      if(errn != 0){ fprintf(stderr, "error in the algorithm of parallel tempering: saving likelihoods\n"); exit(-1); }

      errn = print_likelihood_T(ith_cal, filePh, i_ptemp, i_ptemp.thigh);
      if(errn != 0){ fprintf(stderr, "error in the algorithm of parallel tempering: saving likelihoods (high temp) \n"); exit(-1); }

      nswap = swap_temperature(&i_ptemp, r);
      fprintf(tmpfile, "%d %d\n", ith_cal, nswap);

      for(i= 1; i < i_ptemp.nprocs; i++)
      {
        MPI_Send(&i_ptemp.Ts[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        MPI_Send(&i_ptemp.invTs[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      }
    }

    /* to the next step */
  }

  if(i_ptemp.tf)
  {
    /* last step to receive the given mpi information and write temperature */
    if(my_id == 0){ s_param.T = i_ptemp.Ts[0]; s_param.invT = i_ptemp.invTs[0]; }
    else{
      MPI_Recv(&s_param.T, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
      MPI_Recv(&s_param.invT, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
    }
  }

  fwrite(&s_param.T, sizeof(double), 1, filelogt);
 
  if(filemod1!=NULL){ fclose(filemod1); }
  if(filemodh!=NULL){ fclose(filemodh); }
  /* the file saving models of temperature 1 and hottest temperature only */

  /* close the file of fileP1 and tmpfile -- in the my_id == 0 */
  if(my_id == 0 && i_ptemp.tf){ 
    if(fileP1 !=NULL){ fclose(fileP1); }
    if(tmpfile !=NULL){ fclose(tmpfile); } 
    if(filePh != NULL){ fclose(filePh); }
  }

  /* now close the model count */ 
  if(filemod!=NULL){ fclose(filemod); }
  if(filemod!=NULL){ fclose(fileP); }
  if(filemod!=NULL){ fclose(filelogt); }

  /* close the logfile */
  fprintf(filelog, "accepted %d burn %d/%d out of %d. Prior reject %d. Calculation reject %d. Additional Iteration due to Cal reject %d. Done. \n", i_accept, i_accept_burn,i_param.nburn, i_param.niter-i_param.nburn, ireject, ireject_cal, additer);
  if(iverbose ==1)
  {
    fprintf(stderr, "%03d accepted %d burn %d/%d out of %d. Prior reject %d. Calculation reject %d. Additional Iteration due to Cal reject %d. Done. \n", my_id, i_accept, i_accept_burn, i_param.nburn, i_param.niter-i_param.nburn, ireject, ireject_cal, additer);
  }
  if(filelog!=NULL){ fclose(filelog); }

  /* send information to the 0-th process */
  MAIN_0001: ;

  /* save ad files */
  sprintf(outfilename, "%s.%03d.info", o_param.outname, my_id);
  filelog = fopen(outfilename, "wb");
  fwrite( &nsaved_T1, sizeof(int), 1, filelog);
  fwrite( &nsaved_Th, sizeof(int), 1, filelog);

  fwrite( &i_accept1, sizeof(int), 1, filelog);
  fwrite( &i_accept_burn1, sizeof(int), 1, filelog);

  fwrite( &ith_cal, sizeof(int), 1, filelog);
  fwrite( &i_accept, sizeof(int), 1, filelog);
  fwrite( &i_accept_burn, sizeof(int), 1, filelog);
  fwrite( &ireject, sizeof(int), 1, filelog);
  fwrite( &ireject_cal, sizeof(int), 1, filelog);
  
  fclose(filelog);

  /* write the files ! */
  if(ith_cal > 0)
  {
    /* write Vp */
    sprintf(outfilename, "%s.%03d.Vp.ppd", o_param.outname, my_id);
    filelog = fopen(outfilename, "wb");
    fwrite(pdf_vp, sizeof(float), o_param.nvp*o_param.nh, filelog);
    fclose(filelog); 

    /* write Vs */
    sprintf(outfilename, "%s.%03d.Vs.ppd", o_param.outname, my_id);   
    filelog = fopen(outfilename, "wb");
    fwrite(pdf_vs, sizeof(float), o_param.nvs*o_param.nh, filelog);
    fclose(filelog);

    /* write rho */
    sprintf(outfilename, "%s.%03d.rho.ppd", o_param.outname, my_id);
    filelog = fopen(outfilename, "wb");
    fwrite(pdf_rho, sizeof(float), o_param.nrho*o_param.nh, filelog);
    fclose(filelog);

    /* write kappa */
    sprintf(outfilename, "%s.%03d.kappa.ppd", o_param.outname, my_id);
    filelog = fopen(outfilename, "wb");
    fwrite(pdf_k, sizeof(float), o_param.nk*o_param.nh, filelog);
    fclose(filelog);    
 
    /* write the synthetics pdf */
    if(in_rwave.tf){  
      sprintf(outfilename, "%s.%03d.rwave.ppd", o_param.outname, my_id);
      filelog = fopen(outfilename, "wb");
      fwrite(pdf_rwave, sizeof(float), namps*in_rwave.npts, filelog); 
      fclose(filelog);
    }
    if(in_racf.tf){  
      sprintf(outfilename, "%s.%03d.racf.ppd", o_param.outname, my_id);
      filelog = fopen(outfilename, "wb");
      fwrite(pdf_racf, sizeof(float), namps*in_racf.npts, filelog); 
      fclose(filelog); 
    }
    if(in_zacf.tf){
      sprintf(outfilename, "%s.%03d.zacf.ppd", o_param.outname, my_id);
      filelog = fopen(outfilename, "wb");
      fwrite(pdf_zacf, sizeof(float), namps*in_zacf.npts, filelog);
      fclose(filelog);
    }
    if(in_pacf.tf){
      sprintf(outfilename, "%s.%03d.pacf.ppd", o_param.outname, my_id);
      filelog = fopen(outfilename, "wb");
      fwrite(pdf_pacf, sizeof(float), namps*in_pacf.npts, filelog);
      fclose(filelog);
    }
//    fprintf(stderr, "info sent from %d to 0\n", my_id);
  }

  //MPI_Barrier(MPI_COMM_WORLD);

  /* gather all information */
  if(my_id==0)
  { 
    errn = num_procs; tmpint = 0;
    for(i=0; i < num_procs; i++){ fileexists[i] = 1; }
    while( errn > 0)
    {
      for(i=0; i < num_procs; i++){ 
        sprintf(outfilename, "%s.%03d.info", o_param.outname, i);
        if( fileexists[i] == 1 && access(outfilename, F_OK) == 0){ fileexists[i] = 0; errn--;  }
      }
    }
    fprintf(stderr, "existing all files\n");

    for(i=0; i < num_procs; i++)
    { 
      fprintf(stderr, "received information %03d/%03d\n", i+1, num_procs);
      /* save ad files */
      sprintf(outfilename, "%s.%03d.info", o_param.outname, i);
      filelog = fopen(outfilename, "rb");

      fread( MPI_T1_nsaved+i, sizeof(int), 1, filelog);
      fread( MPI_Th_nsaved+i, sizeof(int), 1, filelog);

      fread( MPI_i_accept1+i, sizeof(int), 1, filelog);
      fread( MPI_i_accept_burn1+i, sizeof(int), 1, filelog);

      fread( MPI_ith_cal+i, sizeof(int), 1, filelog);
      fread( MPI_i_accept+i, sizeof(int), 1, filelog);
      fread( MPI_i_accept_burn+i, sizeof(int), 1, filelog);
      fread( MPI_ireject+i, sizeof(int), 1, filelog);
      fread( MPI_ireject_cal+i, sizeof(int), 1, filelog);

      fclose(filelog);
      /* remove the file */
      remove(outfilename);
 
      tmpint += MPI_T1_nsaved[i];
      errn += MPI_Th_nsaved[i];
    }
    nsaved_T1 = tmpint;
    nsaved_Th = errn; 
    fprintf(stderr, "saving %d models for cold(T=1) and %d models for hot(T_high)\n", nsaved_T1, nsaved_Th);
 
    /* initialize the ppd grids to zero -- initialize */
    for(j=0; j < o_param.nvp*o_param.nh; j++){ pdf_vp[j] = 0; }
    for(j=0; j < o_param.nvs*o_param.nh; j++){ pdf_vs[j] = 0; }
    for(j=0; j < o_param.nrho*o_param.nh; j++){ pdf_rho[j] = 0; }
    for(j=0; j < o_param.nk*o_param.nh; j++){ pdf_k[j] = 0; }
    for(j=0; j < namps*in_rwave.npts; j++){ pdf_rwave[j] = 0; }
    for(j=0; j < namps*in_zacf.npts; j++){ pdf_zacf[j] = 0; }
    for(j=0; j < namps*in_racf.npts; j++){ pdf_racf[j] = 0; }
    for(j=0; j < namps*in_pacf.npts; j++){ pdf_pacf[j] = 0; }
 
    /* merged file */
    sprintf(outfilename, "%s.likelihood.txt", o_param.outname);
    filelog = fopen(outfilename, "wb");
    
    sprintf(outfilename, "%s.models_T1.txt", o_param.outname);
    tmpfile = fopen(outfilename, "wb");
    fwrite(&nsaved_T1, sizeof(int), 1, tmpfile); /* number of models written */
    fwrite(&i_param.wnl, sizeof(int), 1, tmpfile); /* number of layers in one set */

    sprintf(outfilename, "%s.models_Th.txt", o_param.outname);
    filemodh = fopen(outfilename, "wb");
    fwrite(&nsaved_Th, sizeof(int), 1, filemodh); /* number of models written */
    fwrite(&i_param.wnl, sizeof(int), 1, filemodh); /* number of layers in one set */

    sprintf(outfilename, "%s.T.log", o_param.outname);
    filelogt = fopen(outfilename, "wb");
    fwrite(&num_procs, sizeof(int), 1, filelogt);
    fwrite(&i_param.niter, sizeof(int), 1, filelogt);

    /* prepare to print the accepted at T = 1 */
    fwrite(&i_ptemp.ncool, sizeof(int), 1, filelog);
    tmpint = i_ptemp.ncool*(i_param.niter - i_param.nburn);
    fwrite(&tmpint, sizeof(int), 1, filelog);
    tmpint = i_ptemp.ncool*i_param.nburn;
    fwrite(&tmpint, sizeof(int), 1, filelog);
    tmpint = 0;
    for(i=0; i < num_procs; i++){ tmpint += MPI_i_accept1[i]; }
    fwrite(&tmpint, sizeof(int), 1, filelog);
    tmpint = 0;
    for(i=0; i < num_procs; i++){ tmpint += MPI_i_accept_burn1[i]; }
    fwrite(&tmpint, sizeof(int), 1, filelog);

    /* get ready to print each */
    fwrite(&num_procs, sizeof(int), 1, filelog); /* number of chains */
    tmpint = i_param.niter - i_param.nburn;
    fwrite(&tmpint, sizeof(int), 1, filelog); /* total iteration - nburn_in */
    fwrite(&i_param.nburn, sizeof(int), 1, filelog);

    for(i=0; i < num_procs; i++)
    {
      /* print necessary information for the likelihood file */
      fwrite(&i, sizeof(int), 1, filelog); /* ID */
      fwrite(MPI_i_accept+i, sizeof(int), 1, filelog);
      fwrite(MPI_i_accept_burn+i, sizeof(int), 1, filelog);
      fwrite(MPI_ireject+i, sizeof(int), 1, filelog);
      fwrite(MPI_ireject_cal+i, sizeof(int), 1, filelog);

      /* write likelihood of the chain */
      sprintf(outfilename, "%s.%03d.likelihood.txt", o_param.outname, i);
      fileP = fopen(outfilename,"rb");

      sprintf(outfilename, "%s.%03d.T.log", o_param.outname, i);
      filemod = fopen(outfilename, "rb");

      /* read the likelihood file, print it */
      for(j = 0; j < i_param.niter; j++)
      {
        fread(&tmpfloat, sizeof(float), 1, fileP);
        fwrite(&tmpfloat, sizeof(float), 1, filelog);

        fread(&tmpdouble, sizeof(double), 1, filemod);
        tmpfloat = (float) tmpdouble;
        fwrite(&tmpfloat, sizeof(float), 1, filelogt);
      }
      fclose(fileP); fclose(filemod);
      remove(outfilename);
      sprintf(outfilename, "%s.%03d.likelihood.txt", o_param.outname, i);
      remove(outfilename);
       
      /* merge the probability density -- read the saved files. do not forget to remove it in the end */ 
      /* read Vp */
      sprintf(outfilename, "%s.%03d.Vp.ppd", o_param.outname, i);
      if( access(outfilename, F_OK) == 0)
      {
        filemod = fopen(outfilename, "rb");
        for(j=0; j < o_param.nvp*o_param.nh; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_vp[j] += tmpfloat; }
        fclose(filemod);
        /* now remove the file since we read it */ remove(outfilename); 
      }

      /* read Vs */
      sprintf(outfilename, "%s.%03d.Vs.ppd", o_param.outname, i);
      if( access(outfilename, F_OK) == 0)
      {
        filemod = fopen(outfilename, "rb");
        for(j = 0; j < o_param.nvs*o_param.nh; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_vs[j] += tmpfloat; }
        fclose(filemod);
        /* remove */ remove(outfilename);
      }

      /* read rho */
      sprintf(outfilename, "%s.%03d.rho.ppd", o_param.outname, i);
      if( access(outfilename, F_OK) == 0)
      {
        filemod = fopen(outfilename, "rb");
        for(j = 0; j < o_param.nrho*o_param.nh; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_rho[j] += tmpfloat; }
        fclose(filemod);
        /* remove */ remove(outfilename);
      }

      /* read kappa = Vp/Vs */
      sprintf(outfilename, "%s.%03d.kappa.ppd", o_param.outname, i);
      if( access(outfilename, F_OK) == 0)
      {
        filemod = fopen(outfilename, "rb");
        for(j = 0; j < o_param.nk*o_param.nh; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_k[j] += tmpfloat; }
        fclose(filemod);
        /* remove */ remove(outfilename);
      }
 
      /* read the synthetics pdf */
      if(in_rwave.tf)
      {
        sprintf(outfilename, "%s.%03d.rwave.ppd", o_param.outname, i);
        if( access(outfilename, F_OK) == 0)
        {
          filemod = fopen(outfilename, "rb");
          for(j = 0; j < namps*in_rwave.npts; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_rwave[j] += tmpfloat; }
          fclose(filemod);
          remove(outfilename);
        }
      }

      if(in_zacf.tf)
      {
        sprintf(outfilename, "%s.%03d.zacf.ppd", o_param.outname, i);
        if( access(outfilename, F_OK) == 0)
        {
          filemod = fopen(outfilename, "rb");
          for(j = 0; j < namps*in_zacf.npts; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_zacf[j] += tmpfloat; }
          fclose(filemod);
          remove(outfilename);
        }
      }

      if(in_racf.tf)
      {
        sprintf(outfilename, "%s.%03d.racf.ppd", o_param.outname, i);
        if( access(outfilename, F_OK) == 0)
        {
          filemod = fopen(outfilename, "rb");
          for(j = 0; j < namps*in_racf.npts; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_racf[j] += tmpfloat; }
          fclose(filemod);
          remove(outfilename);
        }
      }

      if(in_pacf.tf)
      {
        sprintf(outfilename, "%s.%03d.pacf.ppd", o_param.outname, i);
        if( access(outfilename, F_OK) == 0)
        {
          filemod = fopen(outfilename, "rb");
          for(j = 0; j < namps*in_pacf.npts; j++){ fread(&tmpfloat, sizeof(float), 1, filemod); pdf_pacf[j] += tmpfloat; }
          fclose(filemod);
          remove(outfilename);
        }
      }      

      /* now read model files and write into one */
      sprintf(outfilename, "%s.%03d.models_T1.txt", o_param.outname, i);
      filemod = fopen(outfilename, "rb");
      
      for(j=0; j < MPI_T1_nsaved[i]; j++)
      {
        fread(tmpppd, sizeof(float), i_param.wnl*4, filemod);
        fwrite(tmpppd, sizeof(float), i_param.wnl*4, tmpfile);
      }
      fclose(filemod);
      remove(outfilename);

      sprintf(outfilename, "%s.%03d.models_Th.txt", o_param.outname, i);
      filemod = fopen(outfilename, "rb");
      for(j=0; j < MPI_Th_nsaved[i]; j++)
      {
        fread(tmpppd, sizeof(float), i_param.wnl*4, filemod);
        fwrite(tmpppd, sizeof(float), i_param.wnl*4, filemodh);
      }
      fclose(filemod);
      remove(outfilename); 

      fprintf(stderr, "done merging for process %03d\n", i+1); 
    }
    if(tmpfile!=NULL){ fclose(tmpfile); }
    if(tmpfile!=NULL){ fclose(filelog); }

    /* now get the files ready */
    /* modify the arrays for printing: x_i + 0.5dx */
    change_array_to_printable(&o_param);
    change_array_to_printable_single(namps, amparray, damps);

    /* change to pdf */
    /* convert the count into the  probability density function */
    count_2_ppd_norm(o_param.nh, o_param.nvp, pdf_vp);
    count_2_ppd_norm(o_param.nh, o_param.nvs, pdf_vs);
    count_2_ppd_norm(o_param.nh, o_param.nrho, pdf_rho);
    count_2_ppd_norm(o_param.nh, o_param.nk, pdf_k);
    count_2_ppd(in_rwave.npts, namps, pdf_rwave, nsaved_T1);
    count_2_ppd(in_zacf.npts, namps, pdf_zacf, nsaved_T1);
    count_2_ppd(in_racf.npts, namps, pdf_racf, nsaved_T1);
    count_2_ppd(in_pacf.npts, namps, pdf_pacf, nsaved_T1);

    /* get median structure and then print immediately */
    get_median_struct(o_param.nh, o_param.nvp, o_param.vparr, pdf_vp, parray);
    get_median_struct(o_param.nh, o_param.nvs, o_param.vsarr, pdf_vs, sarray);
    get_median_struct(o_param.nh, o_param.nrho, o_param.rhoarr, pdf_rho, rhoarray);
    get_median_struct(o_param.nh, o_param.nk, o_param.karr, pdf_k, karray);
  
    sprintf(outfilename, "%s.med.txt", o_param.outname);
    print_ppdmodels2(o_param.nh, o_param.harr, parray, sarray, rhoarray, karray, outfilename);
    fprintf(filelog, "wrote median structure %s\n", outfilename);

    /* get highest pdf structure */
    get_maximumdensity_struct(o_param.nh, o_param.nvp, o_param.vparr, pdf_vp, parray);
    get_maximumdensity_struct(o_param.nh, o_param.nvs, o_param.vsarr, pdf_vs, sarray);
    get_maximumdensity_struct(o_param.nh, o_param.nrho, o_param.rhoarr, pdf_rho, rhoarray);
    get_maximumdensity_struct(o_param.nh, o_param.nk, o_param.karr, pdf_k, karray);

    sprintf(outfilename, "%s.max.txt", o_param.outname);
    print_ppdmodels2(o_param.nh, o_param.harr, parray, sarray, rhoarray, karray, outfilename);
    fprintf(filelog, "wrote maximum ppd structure %s\n", outfilename);

    /* print array variables */
    print_array_variables2(o_param);

    /* print the pdf individually */
    sprintf(outfilename, "%s.Vp.ppd", o_param.outname);
    print_pdf_variable(o_param, outfilename, o_param.nvp, o_param.vparr, pdf_vp);
    fprintf(filelog, "wrote Vp probability density function %s\n", outfilename);

    sprintf(outfilename, "%s.Vs.ppd", o_param.outname);
    print_pdf_variable(o_param, outfilename, o_param.nvs, o_param.vsarr, pdf_vs);
    fprintf(filelog, "wrote Vs probability density function %s\n", outfilename);

    sprintf(outfilename, "%s.rho.ppd", o_param.outname);
    print_pdf_variable(o_param, outfilename, o_param.nrho, o_param.rhoarr, pdf_rho);
    fprintf(filelog, "wrote rho probability density function %s\n", outfilename);

    sprintf(outfilename, "%s.kappa.ppd", o_param.outname);
    print_pdf_variable(o_param, outfilename, o_param.nk, o_param.karr, pdf_k);
    fprintf(filelog, "wrote kappa probability density function %s\n", outfilename);
 
    /* print the pdf in a binary */
    if(in_rwave.tf)
    {
      sprintf(outfilename, "%s.rwave.ppd", o_param.outname);
      print_pdf_synthetics(outfilename, in_rwave.npts, in_rwave.time[0], in_rwave.delta, inrwaves + in_rwave.tind0, namps, amparray, pdf_rwave);
    }
    if(in_racf.tf)
    {
      sprintf(outfilename, "%s.racf.ppd", o_param.outname);
      print_pdf_synthetics(outfilename, in_racf.npts, in_racf.time[0], in_racf.delta, inracfs + in_racf.tind0, namps, amparray, pdf_racf);
    }
    if(in_zacf.tf)
    {
      sprintf(outfilename, "%s.zacf.ppd", o_param.outname);
      print_pdf_synthetics(outfilename, in_zacf.npts, in_zacf.time[0], in_zacf.delta, inzacfs + in_zacf.tind0, namps, amparray, pdf_zacf);
    }
    if(in_pacf.tf)
    {
      sprintf(outfilename, "%s.pacf.ppd", o_param.outname);
      print_pdf_synthetics(outfilename, in_pacf.npts, in_pacf.time[0], in_pacf.delta, inpacfs + in_pacf.tind0, namps, amparray, pdf_pacf);
    }
    fprintf(stderr, "merging all done at process 0 \n"); 

    for(i=0; i < num_procs; i++)
    {
      /* At the last step if the file is not erased */ 
      sprintf(outfilename, "%s.%03d.likelihood.txt", o_param.outname, i);
      if(access(outfilename, F_OK)==0){ remove(outfilename); }
      if(!i_ptemp.tf)
      {
        sprintf(outfilename, "%s.%03d.models.txt", o_param.outname, i);
        if(access(outfilename, F_OK)==0){ remove(outfilename); }
      }
    }

  }

 
  for(i=0; i < maxlay; i++)
  {
   gsl_matrix_complex_free(tmpE[i]);
   gsl_matrix_complex_free(tmpinvE[i]);
   gsl_matrix_complex_free(A[i]);
  }
  gsl_matrix_complex_free(Jmatrix);
  gsl_matrix_complex_free(Jmatrix_w);
  gsl_matrix_complex_free(tmpE_w);
  gsl_matrix_complex_free(tmpinvE_w);
  gsl_matrix_complex_free(A_w);
  gsl_matrix_complex_free(mat22); /* allocated space to make things easy */
  gsl_matrix_complex_free(mat66);
  gsl_matrix_complex_free(J11);
  gsl_vector_complex_free(x);
  gsl_vector_complex_free(in);
  gsl_permutation_free(p);

  MAIN_0002: ;
  end = clock();
  tmpdouble = (end-start)/ CLOCKS_PER_SEC; 
  if(my_id==0)
  {  
    fprintf(stdout, "Time: %.1f seconds\n", tmpdouble);
    if(iverbose == 1){ fprintf(stderr, "Time: %.1f seconds\n", tmpdouble); }
  }
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

/* subroutines */
int errmsg(char **input)
{
  fprintf(stderr, "\n Transmitted Reflected Inversion									\n\
  Getting structure using R-ACF, Z-ACF, and stacked R stacked waveform.							\n\n\
  how to run:														\n\
  %s (input text file) [-v(per N iteration)] 										\n\n\
  For help message, type %s h  OR %s -h  										\n\
  Ver 0.0: 2021 11 21 -- finished writing 2021 11 27 -- finished compiling 2021 11 28 					\n\
           Help messages added 2021 11 29										\n\
	   Debugging done through 2021 12 02 - 05									\n\
           Some output writing functions tidied up 2021 12 06								\n\
  Ver 1.5: 2021 12 24 -- changed the whole code to gsl matrix forms and reduced the time asmp				\n\
  Ver 2.0: 2022 01 12 -- changed the code to take errors differently for ACFs 						\n\
  Ver 2.1: 2022 01 14 -- changed the code to make new models changing to whole parameters				\n\
  Ver 3.0: 2022 01 17 -- changed the code to parallellize 								\n\
  Ver 4.0: 2022 01 24 -- changed the code to consider noise effect of ACF from SNR					\n\
  Ver 4.1: 2022 01 28 -- proper treatment of noise (for whitening) and add option to change one parameter at a time	\n\
  Ver 5.0: 2022 03 07 -- changed the code greatly. takes parallel tempering and considers the case of un-stacked. 	\n\
  Ver 5.1: 2022 03 09 -- changed the code to print the models of highest temperature					\n\
  Ver 5.2: 2022 03 16-18 -- changed the code (minor) to merge the pdf by files, not MPI_Send and MPI_Recv 		\n\
  Ver 6.0: 2022 06 29 - 07 22 -- changed the code (major) 								\n\
           Can use the DPG ACFs 											\n\
           Implement the peak shifting effect in the R synthetics							\n\
           Can use velocity ACFs and displacement ACFs									\n\
           Takes more realistic consideration of the multiple rayparameter & snr combination				\n\
           Can use the covariance matrix to calculate likelihood at each step.						\n\
  Ver 6.1: 2022 08 07 -- minor changes regarding handling the velocity, displacement ACFs 				\n\
  Ver 6.2: 2022 08 08 -- minor corrections in printing pacf and vp/vs							\n\
  Ver 7.0: 2022 08 17 -- changed the code to take the LU, LDLT decomposition && have safety points for matrix           \n\
  Ver 7.1: 2022 08 18 -- changed the code to properly take the noise level consideration in velocity, Pressure ACFs     \n\
  Ver 7.2: 2022 08 18 -- two versions made (LU version and LDLT version)                                                \n\
  Ver 7.3: 2022 08 19 -- modification of the velocity and white noise case (white noise in velocity) 			\n\
  HyeJeong Kim\n\n", input[0], input[0], input[0]);

  return 0;
}

void helpmsg(char **input)
{
  fprintf(stderr, "\n Overall explanation of this Transmitted Reflected Inversion\n\
  We take stacked trace with stack error to invert for shallow structures under the ocean.\n\
  For further explanation, type as following.   \n\
  %s p  :: for parameter input text information \n\
  %s i  :: for input data text information      \n\
  %s m  :: for input model text information     \n\
  %s o  :: for the explanation of output files  \n\
  Maximum number of layers is %d !!!            \n\
  2021 11 29 \n\
  Last Updated 2022 08 18 HjK\n\n", input[0], input[0], input[0], input[0], maxlay);  

  return;
}

void helpmsg_output( )
{
  fprintf(stderr, "\n The format of the output of this code 		\n\
  1. probability density of R-stack, Z-ACF, R-ACF, P-ACF    		\n\
  2. probability density of Vp, Vs, Rho by depth 	    		\n\
  3. Temperature array file 				    		\n\
  4. Model file of T == 1						\n\
  5. Model file of T == T_high						\n\
  6. Likelihood and Temperature of each chain throughout interation	\n\
  7. Likelihood of T=1 values and T=T_high values 			\n\
  Written 2022 07 22 HyeJeong Kim\n\n");

  return;
}


void helpmsg_param()
{
  fprintf(stderr, "\n The look of input text file (parameter file):\n\
  To comment, use # mark at the start of the line.\n\
  The format of the parameter file is (content_name) = (the input)\n\
  ##############################						    \n\
  ##### the parameter file #####						    \n\
  ##############################						    \n\
  # input data file names							    \n\
  racf_in =            # Radial ACF upon S incidence				    \n\
  zacf_in =            # Vertical ACF upon P incidence				    \n\
  pacf_in =            # DPG ACF upon P incidence 				    \n\
  rstack_in =          # stacked  radial component after normalization by Z	    \n\
  racf_vel =           # T or t, if you want to use the velocity RACF               \n\
  zacf_vel =           # T or t, if you wish to use the velocity ZACF               \n\
  pacf_trac = 	       # T to use the traction ACF 				    \n\
  # the information of rayparameter and the signal to noise ratio handling	    \n\
  pintv = 	       # the interval of ray parameters to consider 		    \n\
  # specify model information (water, increase)		 			    \n\
  water =              # T for existing (default) F for not.                        \n\
  increase =           # T for monotonic increase (default) F for not               \n\
  # initial model information. F if initial model is not given.			    \n\
  initial_vmod_in =    # should include water layer if water is not 0		    \n\
  # the number of solid layers to solve for.					    \n\
  n_layer =            # nlayer excluding the water layer 			    \n\
  # how to give new candidates --- all at once or one at once (default all)	    \n\
  changeall = 	       # to change one at once f, t to change all at once           \n\
  # iteration number information						    \n\
  n_burn_in =          # do not save until				   	    \n\
  n_iter =             # total number of iterations				    \n\
  n_saving =           # the saving interval. 					    \n\
  # The PARALLEL TEMPERING 							    \n\
  T_high = 	       # the maximum temperature for the temperature array 	    \n\
  r_cool = 	       # the ratio of T = 1 chains out of all processes             \n\
  # to specify the number of noise traces when there exist whitening in ACF         \n\
  # If we are solving for all Vp, Vs, and H for n_layers.		 	    \n\
  # We can set prior probability with a text file 				    \n\
  # (specifying for each layer & paremter)					    \n\
  # The file gives min,max,stddev if we are solving. 				    \n\
  # Or give single value of fixed param value.  				    \n\
  prior_vmod_in =      # the file that has prior model constraints 		    \n\
  # output format: grid sizes for printing.					    \n\
  # the range of printing value comes from minimum value to maximum search value.   \n\
  dvs =                # the grid size for printing				    \n\
  dvp =										    \n\
  dh =										    \n\
  hprint0 = # for H, since there's water layer, it's good to specify lower limit.   \n\
  # output text file root names.                                                    \n\
  outrootname =                                                                     \n\
  \n\n\
  Fin!\n\n");

  return;
}

void helpmsg_input()
{
  fprintf(stderr, "\n The look of input text file (input DATA file):		    \n\
  To comment, use # mark at the start of the line.				    \n\
  We take the stacked file and the error of the stack.				    \n\
  The inputs should be in the order. Different from the parameter file.             \n\
  If it is not in order, something unexpected and unpleasant will happen.	    \n\
 										    \n\
  **The ERROR (standard deviation) file can be given in three ways		    \n\
  1) (file name) 								    \n\
  2) (file name) (water level of error) 					    \n\
  3) (constant value)								    \n\
  4) \"C\" or \"c\" [(waterlevel)] for covariance calculation.			    \n\
     For this to be achieved, there are other conditions to meet.		    \n\
     The files should be all given.						    \n\
     The files given should not be filtered. 					    \n\
     The stacked file should be provided. Stacking option should be true.	    \n\
										    \n\
  ** The NSNR line can be given in three ways					    \n\
  1) number of snr per ray parameter (if ray parameter N > 1, it should be 1 only)  \n\
     ** it should be 1 if no noise is considered. 				    \n\
  2) number of snr per rayp (=1), fraction (0-1) of combination to include.	    \n\
  3) number of snr per rayp(=1), number of highest percentage bins to include.	    \n\
  Note1) The set-up of the bins are given in the parameter file. 		    \n\
  Note2) For the fraction apporach to be taken, the stack should be true. 	    \n\
										    \n\
  ###############################                                                   \n\
  ##### the input data file #####                                                   \n\
  ###############################                                                   \n\
  (number of ray parameters in the stacked) 					    \n\
  stacked file name								    \n\
  ERROR file name (This will be C to use covariance)				    \n\
  (valid for ACF only) Likelihood with individual (F) or stack (T)                  \n\
  NSNR: Number of SNRs per ray parameter -- 1 for one rayp multiple snr with 1 rayp \n\
  a_gauss									    \n\
  post-deconvolution butterworth filter freq0,freq1 (types: F (for false). 	    \n\
  for other cases, put low cut and high cut frequencies.)			    \n\
  w smoothing factor								    \n\
  tstart,tend 									    \n\
  (The information of rays in the following lines as the following:		    \n\
  Case 1: a single ray parameter:						    \n\
    [ray parameter]								    \n\
    [snr 1] ([file name ACF -- when stacking is false or covariance true])	    \n\
    [snr 2] ([file name ACF])			 				    \n\
    ...										    \n\
    [snr nsnr] ([file name ACF])			 			    \n\
  \n\n\
  Case 2: multiple ray parameter: (single snr is only possible)			    \n\
    [rayp 1] [snr 1] ([file name ACF -- when stacking is false or covariance true]) \n\
    ...										    \n\
    [rayp nray] [snr nray] ([file name ACF])			 		    \n\
  \n\n\
  OVERALL, the possible cases considered taking the data input are:		    \n\
  Fin!\n\n");

  return;
}

void helpmsg_model()
{
  fprintf(stderr, "\n The format of input model file & input prior model file in text:    	\n\
  1) The model file (initial model) 						        \n\
  H_layer1 Vp_layer1 Vs_layer1  rho_layer1						\n\
  H_layer2 Vp_layer2 Vs_layer2  rho_layer2						\n\
  ...											\n\
  H_layerN Vp_layerN Vs_layerN 								\n\n\
  Where the layer1 is usually the water layer! 						\n\
  \n\n\
  2) The model prior information file 							\n\
  : We only solve when the full triplet (MIN,MAX,SIGMA) is given 			\n\
  [H_water info]									\n\
  [H_layer 1] [Vp layer 1 info] [Vs layer 1 info] [Rho layer 1 info]			\n\
  ...											\n\
  [Vp layer N info] [Vs layer N info] [Rho layer N info]				\n\n\
  ** the density Rho **									\n\
  ** The density can be either fixed to constant, 					\n\
  ** determined by empirical relation with Vp						\n\
  ** Or made to freely search.								\n\n\
  ** the Vs **										\n\
  ** the Vs information can be fixed to constant Vp/Vs ratio.           		\n\
  ** Range of Vp/Vs cannot be given, but a constant. 					\n\
  ** Give a single number smaller than 100.						\n\
  For example, the following input file for nlayer = 3 will mean			\n\
  5.0,5.5,0.01										\n\
  1.5           1.6,2.0,0.1  0.04,0.4,0.05 						\n\
  0.1,0.5,0.05  2.0,2.5,0.1  0.50,1.2,0.10 						\n\
  5.3           2.0,3.0,0.1 								\n\
  We solve for water depths.								\n\
  We solve Vp and Vs for layer 1 (top layer) with h = 1.5				\n\
  We solve for H, Vp, Vs for layer 2 (mid layer) 					\n\
  We solve for Vs for layer 3 (bottom layer) with Vp fixed 5.3				\n\
  And we only use the empirical relation between Vp and rho				\n\
  NOTE THAT UNITS MUST BE M/S AND M!!							\n\
  \n\n\
  Fin\n\n");
   
  return;
}

void verbosecheck(PARAM iparam, PRIOR iprior, INDATA irwave, INDATA iracf, INDATA izacf, INDATA ipacf, SYNTHPAR sparam, OUTPARAM oparam)
{
  int i, tmpind;
  /* print everything!! */
  fprintf(stderr, "solve for %d solid layers during %d iterations %d unsaved. saving per %d steps\n\n", iparam.nl, iparam.niter, iparam.nburn, iparam.nsave);
  if(irwave.tf)
  {
    fprintf(stderr, "We solve for\nR-stack %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", irwave.stack, irwave.nfiles, irwave.time[0], irwave.time[1]);
  }
  if(iracf.tf)
  {
    fprintf(stderr, "We solve for\nR-acf %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", iracf.stack, iracf.nfiles, iracf.time[0], iracf.time[1]);
  }
  if(izacf.tf)
  {
    fprintf(stderr, "We solve for\nZ-acf %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", izacf.stack, izacf.nfiles, izacf.time[0], izacf.time[1]);
  }
  if(ipacf.tf)
  {
    fprintf(stderr, "We solve for\nP-acf %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", ipacf.stack, ipacf.nfiles, ipacf.time[0], ipacf.time[1]);
  }
  fprintf(stderr, "for the synthetics, we use %d point length \n\n", sparam.nfft);

  fprintf(stderr, "we are solving for the following\n");
  for(i=0; i < iprior.nsolve; i++)
  {
    tmpind = iprior.mindex[i];
    if(tmpind < iparam.nl){ fprintf(stderr, "H :: "); }
    else if(tmpind < iparam.nl*2 && tmpind >= iparam.nl){ fprintf(stderr, "Vp :: "); }
    else if(tmpind < iparam.nl*3 && tmpind >= iparam.nl*2){ fprintf(stderr, "Vs :: "); }
    else{ fprintf(stderr, "Rho :: "); }
    fprintf(stderr, "%d-th parameter min %.4f max %.4f random walk %.4f\n", tmpind, iprior.lowlim[tmpind], iprior.highlim[tmpind], iprior.sigma[tmpind]);
  }
  fprintf(stderr, "\n");

  for(i=0; i < iparam.nl; i++)
  {
    if(iprior.p2s[i] > 0){ fprintf(stderr, "%d-th layer Vp/Vs fixed as %.4f\n", i+1, 1./iprior.p2s[i]); }
  }

  fprintf(stderr, "Then, we save the probability density function and median, maximum pdf structures with %s name\n", oparam.outname);
  fprintf(stderr, "\n\
  dH: %.4f (%.4f - %.4f)\n\
  dVp: %.4f (%.4f - %.4f)\n\
  dVs: %.4f (%.4f - %.4f)\n\
  drho %.4f (%.4f - %.4f)\n\n", oparam.dh, oparam.harr[0], oparam.harr[oparam.nh-1], oparam.dvp, oparam.vparr[0], oparam.vparr[oparam.nvp-1], 
  oparam.dvs, oparam.vsarr[0], oparam.vsarr[oparam.nvs-1], oparam.drho, oparam.rhoarr[0], oparam.rhoarr[oparam.nrho-1]);

  return;
}

void verbosecheck_file(FILE *filestr, PARAM iparam, PRIOR iprior, INDATA irwave, INDATA iracf, INDATA izacf, INDATA ipacf, SYNTHPAR sparam, OUTPARAM oparam)
{
  int i, tmpind;
  /* print everything!! */
  fprintf(filestr, "solve for %d solid layers during %d iterations %d unsaved. Saving per %d steps.\n\n", iparam.nl, iparam.niter, iparam.nburn, iparam.nsave);
  if(irwave.tf)
  {
    fprintf(filestr, "We solve for\nR-stack %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", irwave.stack, irwave.nfiles, irwave.time[0], irwave.time[1]);
  }
  if(iracf.tf)
  {
    fprintf(filestr, "We solve for\nR-acf %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", iracf.stack, iracf.nfiles, iracf.time[0], iracf.time[1]);
  }
  if(izacf.tf)
  {
    fprintf(filestr, "We solve for\nZ-acf %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", izacf.stack, izacf.nfiles, izacf.time[0], izacf.time[1]);
  }
  if(ipacf.tf)
  {
    fprintf(filestr, "We solve for\nP-acf %s\nwhich is stack of %d files. t = %.2f - %.2f \n\n", ipacf.stack, ipacf.nfiles, ipacf.time[0], ipacf.time[1]);
  }
  fprintf(filestr, "for the synthetics, we use %d point length \n\n", sparam.nfft);

  fprintf(filestr, "we are solving for the following\n");
  for(i=0; i < iprior.nsolve; i++)
  {
    tmpind = iprior.mindex[i];
    if(tmpind < iparam.nl){ fprintf(filestr, "H :: "); }
    else if(tmpind < iparam.nl*2 && tmpind >= iparam.nl){ fprintf(filestr, "Vp :: "); }
    else if(tmpind < iparam.nl*3 && tmpind >= iparam.nl*2){ fprintf(filestr, "Vs :: "); }
    else{ fprintf(filestr, "Rho :: "); }
    fprintf(filestr, "%d-th parameter min %.4f max %.4f random walk %.4f\n", tmpind, iprior.lowlim[tmpind], iprior.highlim[tmpind], iprior.sigma[tmpind]);
  }
  fprintf(filestr, "\n");
  for(i=0; i < iparam.nl; i++)
  {
    if(iprior.p2s[i] > 0){ fprintf(filestr, "%d-th layer Vp/Vs fixed as %.4f\n", i+1, 1./iprior.p2s[i]); }
  }

  fprintf(filestr, "Then, we save the probability density function and median, maximum pdf structures with %s name\n", oparam.outname);
  fprintf(filestr, "\n\
  dH: %.4f (%.4f - %.4f)\n\
  dVp: %.4f (%.4f - %.4f)\n\
  dVs: %.4f (%.4f - %.4f)\n\
  drho %.4f (%.4f - %.4f)\n\n", oparam.dh, oparam.harr[0], oparam.harr[oparam.nh-1], oparam.dvp, oparam.vparr[0], oparam.vparr[oparam.nvp-1],
  oparam.dvs, oparam.vsarr[0], oparam.vsarr[oparam.nvs-1], oparam.drho, oparam.rhoarr[0], oparam.rhoarr[oparam.nrho-1]);

  return;
}

void print_iprior(FILE *filestr, PARAM param, PRIOR iprior)
{
  int i = 0;
  fprintf(filestr, "%.4f,%.4f 1500,1500 0,0 1030,1030\n", iprior.lowlim[i], iprior.highlim[i]); 
  for(i=1; i < param.nl; i++)
  {
    fprintf(filestr, "%.4f,%.4f %.4f,%.4f %.4f,%.4f ", iprior.lowlim[i], iprior.highlim[i], iprior.lowlim[param.nl+i-1], iprior.highlim[param.nl + i-1], iprior.lowlim[2*param.nl+i-1], iprior.highlim[2*param.nl + i-1]);
    if(iprior.rho_tf[i-1]){ fprintf(filestr, "%.4f,%.4f\n", iprior.lowlim[3*param.nl+i-1], iprior.highlim[3*param.nl + i-1]); }
    else if(!iprior.rho_tf[i-1] && iprior.lowlim[3*param.nl+i-1]<0){ fprintf(filestr, "%.4f,%.4f\n", vp_empirical_rho_solid(iprior.lowlim[param.nl+i-1]), vp_empirical_rho_solid(iprior.highlim[param.nl + i-1])); }
    else{ fprintf(filestr, "%.4f,%.4f\n", iprior.lowlim[3*param.nl+i-1], iprior.highlim[3*param.nl + i-1]); }
  }

  fprintf(filestr, "%.4f,%.4f %.4f,%.4f %.4f,%.4f ", 100., 100., iprior.lowlim[param.nl+i-1], iprior.highlim[param.nl + i-1], iprior.lowlim[2*param.nl+i-1], iprior.highlim[2*param.nl + i-1]);
  if(iprior.rho_tf[i-1]){ fprintf(filestr, "%.4f,%.4f\n", iprior.lowlim[3*param.nl+i-1], iprior.highlim[3*param.nl + i-1]); }
  else if(!iprior.rho_tf[i-1] && iprior.lowlim[3*param.nl+i-1]<0){ fprintf(filestr, "%.4f,%.4f\n", vp_empirical_rho_solid(iprior.lowlim[param.nl+i-1]), vp_empirical_rho_solid(iprior.highlim[param.nl + i-1])); }
    else{ fprintf(filestr, "%.4f,%.4f\n", iprior.lowlim[3*param.nl+i-1], iprior.highlim[3*param.nl + i-1]); }

  return;
}
