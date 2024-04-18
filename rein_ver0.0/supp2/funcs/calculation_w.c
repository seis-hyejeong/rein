#include <math.h>

#include "calculation_w.h"
#include "matrix_pack.h"
 
/* This basically includes making matrices necessary for calculation */
/* notations follow Crampin 1977 */

/* prepare all christoffel tensors */
void mod_elastic(int nlayer, double qvec[], double i_tensor[], double *rotM, double *o_tensor)
{
  int i;
  double hangle=0;

  hangle = get_rothor_angle(qvec);
  get_2d_rotmat(hangle, rotM);

  for(i=0; i < nlayer; i++)
  {
    rot_tensor(rotM, i_tensor+i*tintv, o_tensor+i*tintv);
  }
  
  return;
}

/* get phase velocities for given ray paramter */
int qsols2(double rayp, double C[], double rho, double *q3sol)
{
  int errn =0; /* if return value errn is 1, it means this step is not succesful */
  double vq[3][2], F[9][3], poly[9], qsoltmp[6*2], tmp[3], tmp1[7], tmp2[7], det_tmp[3][7];
  int i, j, k, m, n;
  double tmpC=0;
  /* parameter space for solving polynomials */
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (7);

  for(i=0; i <3; i++){ vq[i][0]=0; vq[i][1]=0; }
  vq[2][1] = 1;  vq[0][0] = rayp;

  /* obtain matrix F (equation (5) or Keith and Crampin 1977 a) */
  for(j=0; j<3; j++)
  {
    for(m=0; m <3; m++)
    { /* first, initialize the terms */
      F[j*3+m][0] =0; F[j*3+m][1] =0; F[j*3+m][2] =0;
      tmp[0]=0; tmp[1]=0; tmp[2]=0;
      
       for(k=0; k<3; k++)
       {
         for(n=0; n<3; n++)
         {
           tmpC = tensor_get(C,j+1, k+1, m+1, n+1);
           conv(vq[k], 2, vq[n], 2, tmp);
           amul_v(3, tmp, tmpC);
           sum_v(3, tmp, F[j*3+m]);
         }
       }
       if(j==m){ F[j*3+m][0] -= rho; }
    }
  }
  
  /* sub11 */
  conv(F[(2-1)*3+(2-1)], 3, F[(3-1)*3+(3-1)], 3, tmp1);
  conv(F[(2-1)*3+(3-1)], 3, F[(3-1)*3+(2-1)], 3, tmp2); amul_v(5, tmp2, -1);
  sum_v(5, tmp1, tmp2); // tmp2 = tmp1+tmp2;
  conv(F[(1-1)*3+(1-1)], 3, tmp2, 5, det_tmp[0]); /* x11*sub11 */

  /* sub12 */
  conv(F[(2-1)*3+(1-1)], 3, F[(3-1)*3+(3-1)], 3, tmp1); 
  conv(F[(2-1)*3+(3-1)], 3, F[(3-1)*3+(1-1)], 3, tmp2); amul_v(5, tmp2, -1);
  sum_v(5, tmp1, tmp2); // tmp2 = tmp1+tmp2;
  conv(F[(1-1)*3+(2-1)], 3, tmp2, 5, det_tmp[1]);

  /* sub13 */
  conv(F[(2-1)*3+(1-1)], 3, F[(3-1)*3+(2-1)], 3, tmp1); 
  conv(F[(2-1)*3+(2-1)], 3, F[(3-1)*3+(1-1)], 3, tmp2); amul_v(5, tmp2, -1);
  sum_v(5, tmp1, tmp2); 
  conv(F[(1-1)*3+(3-1)], 3, tmp2, 5, det_tmp[2]);

  /* det(F) = sub11 - sub12 + sub13 */
  amul_v(7, det_tmp[1], -1);
  sum_v(7, det_tmp[0], det_tmp[1]);
  sum_v(7, det_tmp[1], det_tmp[2]);
  for(i=0; i < 7; i++){ poly[i] = det_tmp[2][i]; }

  /* solve the polynomial */
  if(gsl_poly_complex_solve( poly, 7, w, qsoltmp) == GSL_EFAILED){ return -1; }; /* qsoltmp is complex solutions */

  gsl_poly_complex_workspace_free(w);

  /* handle obtained solutions */
  for(i=0; i <6; i++)
  {
    if(qsoltmp[2*i+1]>1.E-4 || fabs( qsoltmp[2*i+0]) <1.E-5){ errn =1;}
  }
  if(errn==1){ goto L_0001;}

  gsl_sort(qsoltmp, 2, 6); /* gives it in ascending order (increasing order) */

  /* first: P wave  --> goes to index 1, 4 */
  *(q3sol+0) = qsoltmp[2*2];
  *(q3sol+3) = qsoltmp[3*2];

  /* second: S1 wave --> goes to index 2, 5 */
  *(q3sol+1) = qsoltmp[1*2];
  *(q3sol+4) = qsoltmp[4*2];

  /* third: S2 wave  --> goes to index 3, 6 */
  *(q3sol+2) = qsoltmp[0*2];
  *(q3sol+5) = qsoltmp[5*2];

  L_0001: ;
  return errn;
}

/* isotropic case -- save time using the isotroopic velocities */
int qsols2_iso(double rayp, double C[], double rho, double *q3sol)
{
  int errn =0; /* if return value errn is 1, it means this step is not succesful */
  double ivp, ivs, vp_imped, vs_imped;
  double raypsq =0;

  vp_imped = tensor_get(C, 1, 1, 1, 1);
  vs_imped = tensor_get(C, 1, 3, 1, 3);
  ivp = rho/vp_imped;
  ivs = rho/vs_imped;

  raypsq = rayp*rayp;

  /* first: P wave  --> goes to index 1 (negative), 4 (positive) */
  if( ivp - raypsq  > 0 )
  { 
    *(q3sol +0) = -sqrt(ivp - raypsq);
    *(q3sol +3) = -q3sol[0];
  }
  else{ errn += 1; }

  /* second: S1 wave --> goes to index 2, 5 */
  /* third: S2 wave  --> goes to index 3, 6 */
  if( ivs - raypsq  > 0 )
  { 
    *(q3sol +1) = -sqrt(ivs - raypsq);
    *(q3sol +2) = q3sol[1];
    *(q3sol +4) = -q3sol[1];
    *(q3sol +5) = -q3sol[1];
  }
  else{ errn +=1; }
 
  return errn;
}

/* for propagator part of matrices */ 

/* matrix E in Keith and Crampin 1977 */
/* output is matrix of complex numbers for convenience of overall calculation */
void E_matrix(double C[], double rayp, double qz[], double dvec[], gsl_matrix_complex *Emat)
{ 
  int j, n, m;
  double c, tmp=0;
  gsl_complex tmpcomp;
  c = 1./rayp; /* the expression in Keith & Crampin 1977 use phase velocity! */

  for(n=0; n <6; n++)
  {
    for(j=0; j<3; j++){  GSL_SET_COMPLEX(&tmpcomp, dvec[3*n+j], 0.); 
      gsl_matrix_complex_set(Emat, j, n, tmpcomp);
    }
  }

  for(n=0; n <6; n++)
  {
    for(j=0; j<3; j++)
    {  
      tmp=0;
      for(m=0; m <3; m++)
      {
        tmp += dvec[3*n+m]*(tensor_get(C, j+1, 3, m+1, 1)+ tensor_get(C, j+1, 3, m+1, 3)*c*qz[n]);
      }
      GSL_SET_COMPLEX(&tmpcomp, tmp, 0.); 
      gsl_matrix_complex_set( Emat, j+3, n, tmpcomp);
    }
  }

  return;
}

/* the water Einv_matrix (difference is the dimension of the matrix */
void Einv_w_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *Einvmat)
{
  int i, j, s;
  gsl_complex tmpcomp;
 
  gsl_permutation * p = gsl_permutation_alloc(2);
  gsl_matrix * Ereal = gsl_matrix_calloc(2,2);
  gsl_matrix * invE = gsl_matrix_calloc(2,2);

  /* make it as real matrix to invert it ! */
  for(i=0; i < 2; i++)
  {
    for(j=0; j < 2; j++)
    {
      gsl_matrix_set(Ereal, i, j, GSL_REAL( gsl_matrix_complex_get(Emat, i, j) ));
    }
  }

  gsl_linalg_LU_decomp(Ereal, p, &s);
  gsl_linalg_LU_invert(Ereal, p, invE);

  for(i=0; i < 2; i++)
  {
    for(j=0; j < 2; j++)
    {
      GSL_SET_COMPLEX(&tmpcomp, gsl_matrix_get(invE, i, j), 0); 
      gsl_matrix_complex_set(Einvmat, i, j, tmpcomp);
    } 
  }

  gsl_matrix_free(Ereal); gsl_matrix_free(invE);
  gsl_permutation_free(p);

  return;
}
 
void Einv_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *Einvmat)
{    
  int i, j, s;
  gsl_complex tmpcomp;

  gsl_permutation * p = gsl_permutation_alloc(6);
  gsl_matrix * Ereal = gsl_matrix_calloc(6,6);
  gsl_matrix * invE = gsl_matrix_calloc(6,6);

  /* make it as real matrix to invert it! */
  for(i=0; i< 6; i++)
  {    
    for(j=0; j<6; j++)
    {
      gsl_matrix_set(Ereal, i, j, GSL_REAL( gsl_matrix_complex_get(Emat, i, j) ));
    }
  }

  gsl_linalg_LU_decomp(Ereal, p, &s);
  gsl_linalg_LU_invert(Ereal, p, invE);

  for(i=0; i <6; i++)
  {
    for(j=0; j<6; j++)
    {
      GSL_SET_COMPLEX(&tmpcomp, gsl_matrix_get(invE, i, j), 0);
      gsl_matrix_complex_set(Einvmat, i, j, tmpcomp);
    }
  }

  gsl_matrix_free(Ereal); gsl_matrix_free(invE);
  gsl_permutation_free(p);

  return;
}

/* get A matrix in water */
void A_w_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *E_w_invmat, gsl_complex expw[], gsl_matrix_complex *A, gsl_matrix_complex *Dw_w)
{
  int i, j, k, n;

  int irow[2] = {2, 5};
  int icol[2] = {0, 3};

  for(i=0; i<2; i++) /* col index */
  {
    n = icol[i];
    for(j=0; j <2; j++) /* row index */
    {
      k = irow[j];
      gsl_matrix_complex_set(Dw_w, j, i, gsl_complex_mul(expw[n], gsl_matrix_complex_get(Emat, k, n)) ); /* D = expmat*E */ 
    }
  }

  /* A = D*invE */
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, Dw_w, E_w_invmat, GSL_COMPLEX_ZERO, A);

  return;
}

/* A = D*(invE) */ /* follow exactly Keith and Crampin 1970 */
void A_matrix(gsl_matrix_complex *Emat, gsl_matrix_complex *Einvmat, gsl_complex expw[], gsl_matrix_complex *A, gsl_matrix_complex *Dw)
{
  int j, n;

  for(n=0; n<6; n++)
  {
    for(j=0; j<6; j++)
    {
      gsl_matrix_complex_set( Dw, j, n, gsl_complex_mul(expw[n], gsl_matrix_complex_get( Emat, j, n)) );
      /* invD = iExpmat*invE */
    }
  }
  
  /* A = D*invE */
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, Dw, Einvmat, GSL_COMPLEX_ZERO, A);
 
//  gsl_matrix_complex_free(Dw);

  return;
}

/* calculating the expw in advance --> the length of expw is 6 x nlayer x nfft */
void cal_expw(int nlayer, double qzs[], float dms[], double dw, int nfft, gsl_complex *expw)
{
  int i, j, k, iind, layerintv, kind, freqintv, nfft_2;
  double dm;
  double dtheta, dc, ds, cos_0, sin_0, cos_1, sin_1;
 
  nfft_2 = nfft/2;
  layerintv = 6*nfft;
  freqintv = 6;

  iind = 0; kind = 0;
  for(i=0; i < nlayer; i++) /* (6*nfft)*i + (6*k) + j */
  { dm = dms[i];
    //fprintf(stderr, "dm: %f\n", dm);
    for(j=0; j < 6; j++)
    {
      kind = 0;

      dtheta = (-1.)*dw*dm*qzs[i*6+j];
      dc = cos(dtheta); ds = sin(dtheta);
      cos_0 = 1.; sin_0 = 0.;
      for(k=0; k < nfft_2; k++)
      {
        GSL_SET_COMPLEX( &expw[ iind + kind + j],  cos_0, sin_0);
        cos_1 = cos_0*dc - sin_0*ds;
        sin_1 = sin_0*dc + cos_0*ds;
        
        cos_0 = cos_1; sin_0 = sin_1;
        kind += freqintv; 
      }
      cos_0 = cos( dw*nfft_2*dm*qzs[i*6+j] ); sin_0 = sin( dw*nfft_2*dm*qzs[i*6+j] );
      for(k=nfft_2; k < nfft; k++)
      {
        GSL_SET_COMPLEX( &expw[ iind + kind + j],  cos_0, sin_0);
        cos_1 = cos_0*dc - sin_0*ds;
        sin_1 = sin_0*dc + cos_0*ds;

        cos_0 = cos_1; sin_0 = sin_1; 
        kind += freqintv;
      }    
    }
    iind+= layerintv;
  }

  return;
}


double cal_poly(int order, double *coeffs, double x)
{ // it is highest order. not number of terms (for 4th order polynomial it is 4)
    int i;
    double xx, output=0;
    xx=1;

    for(i=0; i <=order; i++)
    {
       output += *(coeffs+i)*xx;
        xx = xx*x;
    }
    return output;
}

/***** MAIN propagator part!!!!!******************/
/* fn = inv(En)*A(n-1)*A(n-2)*...A1*v0 */
/* fn = J*v0
 * v0 = inv(J)*fn --> we will get matrix J
 */ 
/* fn = (f1 f2 f3 f4 f5 f6) and v0 = (u01, u02, u03, 0, 0, 0) */
void Jw_matrix(int nlayer, gsl_matrix_complex *invEn, gsl_matrix_complex **A, gsl_matrix_complex *J, gsl_matrix_complex *Jtmp)
{
    int i;
 
    /* initialize matrix J */
    gsl_matrix_complex_set_identity(J);
    
    for(i=0; i < nlayer-1; i++) /* add at left */
    {
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, A[i], J, GSL_COMPLEX_ZERO, Jtmp);  
      /* replace J by Jtmp */
      gsl_matrix_complex_memcpy(J, Jtmp);
    }
    
    /*multiply invEn at left */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, invEn, J, GSL_COMPLEX_ZERO, Jtmp);
    gsl_matrix_complex_memcpy(J, Jtmp);
    
    return;
}

void Jw_matrix2(int nlayer, gsl_matrix_complex **A, gsl_matrix_complex *J, gsl_matrix_complex *Jtmp)
{
  int i;

  /* initialize matrix J */
  gsl_matrix_complex_set_identity(J);

  for(i=0; i < nlayer; i++) /* add at left, we will multiply all nlayers */
  {
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, A[i], J, GSL_COMPLEX_ZERO, Jtmp);
    /* replace J by Jtmp */
    gsl_matrix_complex_memcpy(J, Jtmp);
  }

  /* we do not need to add invE at the left*/

  return;
}

/* calculate for selected input */
/* we can have either P input or SV input */
/* solve matrix equation...! */
/* iinci =1 -> P wave incidence, iinci =2 -> S wave */
void get_response(int iinci, gsl_matrix_complex *Jwmat, gsl_complex *surf_u, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p)
{ 
  int i,j, s;
  gsl_complex tmp;

  /* in = J11*(x) */
  
  /* initialize gsl complex array "in" */
  gsl_vector_complex_set_zero(in);
  GSL_SET_COMPLEX(&tmp, 1., 0.);
  gsl_vector_complex_set(in, iinci-1, tmp);

  /* extract part of matrix */
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      gsl_matrix_complex_set(J11, i, j, gsl_matrix_complex_get(Jwmat, i, j));
    }
  } 
//  *J11 = gsl_matrix_complex_submatrix(Jwmat, 0, 0, 3, 3).matrix;
  gsl_linalg_complex_LU_decomp(J11, p, &s);
  gsl_linalg_complex_LU_solve(J11, p, in, x);

  for(i=0; i<3; i++)
  {
    *(surf_u+i) = gsl_vector_complex_get(x, i);
    //fprintf(stderr, "%f %f\n", GSL_REAL(tmp), GSL_IMAG(tmp));
//    if( fabs(GSL_IMAG(tmp)) > 1E-4){ fprintf(stderr, "it should not be like this...! %f \n", GSL_IMAG(tmp)); exit(-1);}
  } 
  for(i=3; i < 6; i++)
  {
    *(surf_u+i) = GSL_COMPLEX_ZERO;
  }

  return;
}

/* get response in case of existing water layer */
void get_response_w(int iinci, gsl_matrix_complex *Jwmat, gsl_complex *surf_u, gsl_matrix_complex *J11, gsl_vector_complex *x, gsl_vector_complex *in, gsl_permutation *p)
{
  int i,j, s;
  gsl_complex tmp;

  /* in = J11*(x) */

  /* initialize gsl complex array "in" */
  gsl_vector_complex_set_zero(in);
  GSL_SET_COMPLEX(&tmp, 1., 0.);
  gsl_vector_complex_set(in, iinci-1, tmp);

  /* copy the matrix to the one we will be using */
  for(i=0; i<5; i++)
  {
    for(j=0; j<5; j++)
    {
      gsl_matrix_complex_set(J11, i, j, gsl_matrix_complex_get(Jwmat, i, j));
    }
  }

//  *J11 = gsl_matrix_complex_submatrix(Jwmat, 0, 0, 5, 5).matrix;

  gsl_linalg_complex_LU_decomp(J11, p, &s);
  gsl_linalg_complex_LU_solve(J11, p, in, x);
//  for(i=0; i < 5; i++){ tmp = gsl_vector_complex_get(x, i); fprintf(stderr, "%d: %f  + %f i\n",i, GSL_REAL(tmp), GSL_IMAG(tmp)); }
  for(i=0; i<3; i++)
  {
    *(surf_u+i) = gsl_vector_complex_get(x, i);
  }
  for(i=3; i < 5; i++)
  {
    *(surf_u+i) = GSL_COMPLEX_ZERO;
  }
  *(surf_u + 5) = gsl_vector_complex_get(x, 3); /* the tau33 */

  return;
}


/* get convolution of source signal and response to get output */
/* u1(R), u2(T), u3(Z) is all organized by frequency and signal is organized also in same frequency */
void conv_source(gsl_complex *u1, gsl_complex *u2, gsl_complex *u3, gsl_complex source[], int nfft)
{
  int i;
  for(i=0; i < nfft; i++)
  {
    *(u1+i) = gsl_complex_mul(*(u1+i), source[i]); 
    *(u2+i) = gsl_complex_mul(*(u2+i), source[i]); 
    *(u3+i) = gsl_complex_mul(*(u3+i), source[i]); 
  }

  return;
}

void source(int nfft, double dt, double gaussf, gsl_complex *source)
{
  int i;
  double x[nfft];
  for(i=0; i < nfft; i++)
  {
      *(x+i) = exp( -pow(dt*(i-nfft/2),2)/(4*gaussf*gaussf));
  }
  trace_fft(x, nfft, source);
  trace_shift_amp(-(nfft/2)*dt, dt, nfft, 1, source);

  return;
}

/* rotate back to NEZ */
void output_invrot(int npts, double rotM [], double u1[], double u2[])
{
  int i;
  double tmpur[3] = {0, 0, 0};
  double tmpu[3] = {0, 0, 0};
  
  for(i=0; i < npts; i++)
  {
    tmpur[0] = u1[i]; tmpur[1] = u2[i];
    mat_rot_inv(rotM, tmpur, tmpu);
    u1[i] = tmpu[0]; u2[i] = tmpu[1]; 
  }
  
  return;
}

/* calculate P arrival time */
double p_tt(int nlayer, double qz_p[], double thickness[])
{
  int i;
  double tt=0;

  for(i=0; i < nlayer; i++)
  {
    tt += qz_p[i]*thickness[i];
  }

  return tt;
}
