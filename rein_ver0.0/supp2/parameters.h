#ifndef _parameters_h_
#define _parameters_h_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_complex.h>

#define nbuffer 1
#define strlng 300
#define n_max 200
#define maxprocs 480
#define max_length 2010
#define TLENGTH 100
#define inputR 2
#define inputZ 1
#define damps 0.005
#define amps0 -0.8
#define namps 401 
#define maxnfft 4096
#define statk 1.5
#define kmin0 1.1
#define kmin 1.5
#define kmax 10
#define dk 0.05
#define watervp 1500.
#define waterrho 1030.
#define minrayp 3.0E-5
#define maxrayp 16.0E-5
#define defaultdp 5.0E-6
#define Ulimit -10000000000000
#define gausslim 0.1
#define magicmul 10
#define waterlevelfrac 0.1

/* define new data type -- the parameter set up */
typedef struct parameter {
  int nl;
  int wnl;
  int nl3;
  int iwater;
  bool increase;
  int niter;
  int nburn;
  int nsave;
  bool changeall;
  char priorfile[strlng];
} PARAM;

typedef struct parallelt {
  bool tf;
  int nprocs;
  int ncool;
  double rcool;
  double thigh;
  double Ts[maxprocs]; 
  double invTs[maxprocs];
  float likelihood[maxprocs];
} PARATEMP;

typedef struct priorparam {
  int nsolve;
  bool rho_tf[n_max];
  bool p_tf[n_max];
  bool s_tf[n_max];
  bool h_tf[n_max];
  double p2s[n_max];
  int mindex[n_max*4];
  double highlim[n_max*4];
  double lowlim[n_max*4];
  double sigma[n_max*4];
} PRIOR;

typedef struct in_data {
  bool ivelocity;
  bool tf;
  bool stack_tf;
  bool covar_tf;
  char name[strlng];
  char stack[strlng];
  char error[strlng];
  double time[2];
  double delta;
  int tind0;
  int npts;
  int nrayp;
  int nsnrs[n_max];
  int nfiles;
  double rayp[n_max];
  double snrs[n_max]; 
  double snrs_factor[n_max];
  int startind[n_max];
  double gauss;
  double ffilter[2];
  double w;
  float frac;
  float waterlvl;
  int nbin;
} INDATA;

typedef struct org_data {
  size_t nbin;
  double dp;
  double arrp[n_max];
} DATAORG;

typedef struct synth_param {
  int nfft;
  int racfind[2];
  int racfind0;
  int zacfind[2];
  int zacfind0;
  int pacfind[2];
  int pacfind0;
  int rwaveind[2];
  int rwaveind0;
  double dw;
  double T;
  double invT;
} SYNTHPAR;

typedef struct out_param {
  double dvs;
  double dvp;
  double drho;
  double dh;
  double h0; /* minimum depth to start printing */
  int nvs;
  int nvp;
  int nrho;
  int nh;
  int nk;
  double vsarr[max_length];
  double vparr[max_length];
  double rhoarr[max_length];
  double harr[max_length];
  double karr[max_length];
  char outname[strlng];
  int nsaved_chain;
} OUTPARAM;

/* define the model variable */
typedef struct model {
  int nl;
  double modvec[n_max*4];
  float h[n_max];
  float vp[n_max];
  float vs[n_max];
  float rho[n_max];
} MODEL;

typedef struct likelihood {
  double racf;
  double zacf;
  double pacf;
  double rwave;
  float pall;
} LIKELI;

static LIKELI likelihood_null = {
  0., 0., 0., 0., 0.
};

/* null values of the PARAMETERs */
static SYNTHPAR synthparam_null = {
  0,
  {0., 0.},
  0,
  {0., 0.},
  0,
  {0., 0.,},
  0,
  {0., 0.,},
  0,
  0.,
  1.,
  1.
};


#endif
