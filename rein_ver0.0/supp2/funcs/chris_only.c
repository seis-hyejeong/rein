#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gsl/gsl_eigen.h"

#include "matrix_pack.h"
#include "chris_only.h"

/* axis system rotation clockwise horizontally (: azimuth of x1 axis occurs looked from above) 
 * --> angle "lambda" 
 * axis system rotation clockwise w.r.t x2 axis( dipping of x1 occurs) --> angle "eta"
 *
 * for hexagonal symmetry, azi is azi of x3 and tilt is tilt of axis from x3 axis */

/* receiving matrices and angles */
/* input format for velocity structure
 * (nlayer)
 * (thick) (rho)    (isoflag)    (azi)    (tilt)    (discon_str)    (discon_dip)
 * c11	c12	c13	c14	c15	c16
 * c21	c22	c23	c24	c25	c26
 * c31	c32	c33	c34	c35	c36
 * c41	c42	c43	c44	c45	c46
 * c51	c52	c53	c54	c55	c56
 * c61	c62	c63	c64	c65	c66
 */

/* with received rotation angles, make rotation matrix */
void get_struct_grotmat(double trend, double tilt, double *M)
{
  /* inclination measured from horizontal. trend is calculated from CW from North */
  /* eta: trend, direction of hexagonal symmetric axis from North in Clockwise direction
   * lambda: dip, angle measured from vertical downward axis (x3)
   */
  double eta, lambda;
  lambda = trend*deg2rad;
  eta = tilt*deg2rad;

  matrix_set(3, M, 1, 1, cos(eta)*cos(lambda));
  matrix_set(3, M, 1, 2, cos(eta)*sin(lambda));
  matrix_set(3, M, 1, 3, (-1)*sin(eta));
  matrix_set(3, M, 2, 1, (-1)*sin(lambda));
  matrix_set(3, M, 2, 2, cos(lambda));
  matrix_set(3, M, 2, 3, 0);
  matrix_set(3, M, 3, 1, sin(eta)*cos(lambda));
  matrix_set(3, M, 3, 2, sin(eta)*sin(lambda));
  matrix_set(3, M, 3, 3, cos(eta));

  return;
}

/* get oscillation vector for that wave from propagating slowness */
/* for isotropic case --> which we worry degeneracy */

/* it is getting oscillation vectors from eigen vector approach */
int cal_oscillvec_eigen(double const C[], double const rho, double const rayp, double *qzsol, double *oscillvec)
{
  /* qsol is arranged in ascending order of absolute values and negative first, than positive 
   * output of this function will be rearranged as qP-qSV-qSH in order of negative first, positive after 
   * Note that input C(christoffel tensor) should be rotated to coordinate of qsolutions.
   * Eigen value for matrix F(j,m) = C(j,k,m,n)*q(k)*q(n) should be same as "rho" and the corresponding eigen vector will have will have oscillation vector. 
   */

  int i, j, k, errn=0; /* the variable "ind" helps ordering the values -> P - SV - SH - P- SV - SH */
  double Ftmp[9], qtmp[3] = {rayp, 0., 0.};
  double new_q3sol[6], normalos[3], tmp_oscill[18];
  gsl_matrix *F = gsl_matrix_calloc(3,3);

  int indcount;
  gsl_complex evaltmp, evec_tmp[3*2];

  /* space for solving eigenvalue problem */
  gsl_eigen_nonsymmv_workspace *ws = gsl_eigen_nonsymmv_alloc(3);
  gsl_vector_complex *e_vals = gsl_vector_complex_alloc(3);
  gsl_matrix_complex *e_vecs = gsl_matrix_complex_alloc(3,3);
  
  for(i=0; i < 6; i++)
  {
    qtmp[2] = *(qzsol+i);
    F_matrix(C, qtmp, Ftmp);

    for(j=0; j<3; j++)
    {
      for(k=0; k<3; k++){ gsl_matrix_set(F, j, k, *(Ftmp+j*3+k)); }
    }    /* now we have the matrix */

    gsl_eigen_nonsymmv(F, e_vals, e_vecs, ws);
    
    indcount=0;
    for(j=0; j<3; j++)
    {
      evaltmp = gsl_vector_complex_get(e_vals, j);
//      fprintf(stderr, "eigen value: %f\n", GSL_REAL(evaltmp));
//      fprintf(stderr, "rho is %f\n", rho);
      if( fabs((GSL_REAL(evaltmp)-rho)/rho) < 1E-5 && GSL_IMAG(evaltmp) < 1E-5)
      {
        evec_tmp[3*indcount+0] = gsl_matrix_complex_get(e_vecs, 0, j);
        evec_tmp[3*indcount+1] = gsl_matrix_complex_get(e_vecs, 1, j);
        evec_tmp[3*indcount+2] = gsl_matrix_complex_get(e_vecs, 2, j);
//        fprintf(stderr, "%d vec: %f %f %f\n", indcount, GSL_REAL(evec_tmp[3*indcount+0]), GSL_REAL(evec_tmp[3*indcount+1]), GSL_REAL(evec_tmp[3*indcount+2]));
        indcount++;
      } 
    }

    if(indcount==1)
    { 
      for(j=0; j<3; j++)
      {
        if( GSL_IMAG(evec_tmp[j]) < 1E-4){ *(tmp_oscill+3*i+j) = GSL_REAL(evec_tmp[j]); }
      }
    }

    else if(indcount==2)
    {
      //fprintf(stderr, "rayp: %f and qz: %f\n", rayp, *(qzsol+i));
      /* first check whether the solution of 6-th polynomial also shows similarity */      
      if( fabs((*(qzsol+i) - *(qzsol+i+1))/ (*qzsol+i)) > 1E-4)
      { 
        fprintf(stderr, "s1 s2 not same but two eigen vec????? %f %f\n", *(qzsol+i), *(qzsol+i+1)); 
        errn=1; goto L_0001; 
      }
//      if(i<3)
//      {
      *(tmp_oscill+3*(i+1)+0)=0; *(tmp_oscill+3*(i+1)+1)=1; *(tmp_oscill+3*(i+1)+2)=0;
//      }
//      else
//      { *(tmp_oscill+3*(i+1)+0)=0; *(tmp_oscill+3*(i+1)+1)=-1; *(tmp_oscill+3*(i+1)+2)=0;}
      
      for(j=0; j<3; j++)
      {
        evaltmp = gsl_vector_complex_get(e_vals, j);
        if( (GSL_REAL(evaltmp)-rho)/rho > 1E-3)
        {  
          normalos[0] = GSL_REAL(gsl_matrix_complex_get(e_vecs, 0, j));
          normalos[1] = GSL_REAL(gsl_matrix_complex_get(e_vecs, 1, j));
          normalos[2] = GSL_REAL(gsl_matrix_complex_get(e_vecs, 2, j));
        }
      }
      /* get the rest one */
      vcross(normalos, tmp_oscill+3*(i+1), tmp_oscill+3*i);
      //fprintf(stderr, "%f %f %f sh %f %f %f\n", *(tmp_oscill+3*i+0), *(tmp_oscill+3*i+1), *(tmp_oscill+3*i+2), *(tmp_oscill+3*i+3), *(tmp_oscill+3*i+4), *(tmp_oscill+3*i+5));
      i++; // we increase additional i .
    }

    else
    { /* case of P and S1, S2 waves having same phase velocity or some complex numbers causing problem */
      fprintf(stderr, "there is error in your algorithm..!!! Please fix this.\n"); 
      errn=1; goto L_0001; 
    }
  }

  /* change order into P-SV-SH */
  /* order of P does not need to be changed */
  new_q3sol[0] = *(qzsol+0);
  new_q3sol[3] = *(qzsol+3);
//  fprintf(stderr, "1: %f 2: %f\n", *(qzsol+0),*(qzsol+3));
  for(j=0; j<3; j++)
  {
    *(oscillvec+0*3+j) = *(tmp_oscill+0*3+j);
    *(oscillvec+3*3+j) = *(tmp_oscill+3*3+j);
  }

  /* now go for upgoing SV, SH */
  if( fabs(*(tmp_oscill+1*3+1)) > fabs(*(tmp_oscill+2*3+1)) )
  {
    new_q3sol[1] = *(qzsol+2);
    new_q3sol[2] = *(qzsol+1);
    for(j=0; j<3; j++)
    {
      *(oscillvec+1*3+j) = *(tmp_oscill+2*3+j);
      *(oscillvec+2*3+j) = *(tmp_oscill+1*3+j);
    }
  }
  else
  {
    new_q3sol[1] = *(qzsol+1);
    new_q3sol[2] = *(qzsol+2);
    for(j=0; j<3; j++)
    {
      *(oscillvec+1*3+j) = *(tmp_oscill+1*3+j);
      *(oscillvec+2*3+j) = *(tmp_oscill+2*3+j);
    }
  }

  /* go for downgoing SV and SH */
  if( fabs(*(tmp_oscill+4*3+1)) > fabs(*(tmp_oscill+5*3+1)))
  {
    new_q3sol[4] = *(qzsol+5);
    new_q3sol[5] = *(qzsol+4);
    for(j=0; j<3; j++)
    {
      *(oscillvec+4*3+j) = *(tmp_oscill+5*3+j);
      *(oscillvec+5*3+j) = *(tmp_oscill+4*3+j);
    }
  }
  else
  {
    new_q3sol[4] = *(qzsol+4);
    new_q3sol[5] = *(qzsol+5);
    for(j=0; j<3; j++)
    {
      *(oscillvec+4*3+j) = *(tmp_oscill+4*3+j);
      *(oscillvec+5*3+j) = *(tmp_oscill+5*3+j);
    }
  }

  for(i=0; i<6; i++)
  {
    if(fabs(*(oscillvec+i*3+0)) > 1E-4)
    {
      double sgn=0;
      sgn = (*(oscillvec+i*3+0)/fabs(*(oscillvec+i*3+0)));
//      fprintf(stderr, "%d's %f sign is %f\n", i+1, *(oscillvec+i*3+0), sgn);
      for(j=0; j<3; j++){ *(oscillvec+i*3+j) *= sgn; }
    }
  }

  /* copy new_q3sol to qzsol memory location */
  for(j=0; j<6; j++){ *(qzsol+j) = new_q3sol[j]; }
//  for(i=0; i <6; i++)
//  {
//    if(i==0 || i==3)
//    {
//      fprintf(stderr, "oscillvec of %d %f: %f %f %f\n", i+1,  qzsol[i], *(oscillvec+3*i+0),  *(oscillvec+3*i+1),  *(oscillvec+3*i+2));
//    }
//  }


  L_0001: ;
  gsl_vector_complex_free(e_vals);
  gsl_matrix_complex_free(e_vecs);
  gsl_eigen_nonsymmv_free(ws);

  return errn;
}

int cal_oscillvec_iso(double const C[], double const rho, double const rayp, double *qzsol, double *oscillvec)
{
  int errn =0; /* if return value errn is 1, it means this step is not succesful */
  double vp, vs, vp_imped, vs_imped;
  double qparallel[3];

  vp_imped = tensor_get(C, 1, 1, 1, 1);
  vs_imped = tensor_get(C, 1, 3, 1, 3);
  vp = 1./sqrt(rho/vp_imped);
  vs = 1./sqrt(rho/vs_imped);

  /* SH is purely horizontal */
  *(oscillvec+3*(3-1)+0) = 0; *(oscillvec+3*(3-1)+1) = 1; *(oscillvec+3*(3-1)+2) = 0;
  *(oscillvec+3*(6-1)+0) = 0; *(oscillvec+3*(6-1)+1) = 1; *(oscillvec+3*(6-1)+2) = 0;

  /* P wave will be parallel */
  *(oscillvec+3*(1-1)+0) = rayp*vp; *(oscillvec+3*(1-1)+1)=0; *(oscillvec+3*(1-1)+2) = *(qzsol+0)*vp;
  *(oscillvec+3*(4-1)+0) = rayp*vp; *(oscillvec+3*(4-1)+1)=0; *(oscillvec+3*(4-1)+2) = *(qzsol+3)*vp;

  /* SV will be perpendicular to both SH oscillation and SV propagation direction */
  qparallel[0] = rayp*vs; qparallel[1]=0; qparallel[2] = *(qzsol+1)*vs;
  vcross(qparallel, oscillvec+3*(3-1), oscillvec+3*(2-1)); /* SV = SVprop x SH */
  qparallel[0] = rayp*vs; qparallel[1]=0; qparallel[2] = *(qzsol+4)*vs;
  vcross(qparallel, oscillvec+3*(6-1), oscillvec+3*(5-1));

/*  for(i=0; i <6; i++)
  {
    fprintf(stderr, "oscillvec of %d %f: %f %f %f\n", i+1,  qzsol[i], *(oscillvec+3*i+0),  *(oscillvec+3*i+1),  *(oscillvec+3*i+2));
  }
*/
  return errn;
}

/* from given Christoffel elastic tensor and propagating slowness, calculate the 3 by 3 matrix */
void F_matrix(double const C[], double const qvec[], double *F)
{
  int j, k, m, n;
  double tmp;

  for(j=1; j<=3; j++)
  {
    for(m=1; m<=3; m++)
    {
      matrix_set(3, F, j, m, 0.);
      tmp =0.;

      for(k=1; k<=3; k++)
      {
        for(n=1; n<=3; n++)
        {
          tmp += tensor_get(C, j, k, m, n)*qvec[k-1]*qvec[n-1];
        }
      }

      matrix_set(3, F, j, m, tmp);
    }
  }

  return;
}

/* calculating isotropic part velocities? */

/* get propagating vector (unit vector) */
int get_3dvec(double const aintv, double *pvec)
{
  int i;
  double theta, phi=0;
  double intv = aintv*deg2rad;
  
  /* pvec = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
   * theta: 0 ~ 180, phi: 0~360
   */

  i=0; /* i is basically counting variable */
  while(phi <= 2*M_PI)
  {
    theta=0.;
    while(theta < M_PI)
    {
      *(pvec + i*3+0) = sin(theta)*cos(phi);
      *(pvec + i*3+1) = sin(theta)*sin(phi);
      *(pvec + i*3+2) = cos(theta);
      theta += intv;
      i++;
    }
    phi += intv;
  }

  return i;
}

/* solving eigen value-vector problem for given propagating direction */
/* calculate phase velocity surface in 3D */
int get_vel(int const ncount, double const pvec[], double const C[], double const rho, double *vel, double *oscills)
{
  int i, j, k, m, n;
  int errn=0;
  double tmp, ptmp[3];
  gsl_matrix * D = gsl_matrix_calloc(3,3);
  gsl_complex eval_tmp, evec_tmp;

  /* space for solving eigen problems */
  gsl_eigen_nonsymmv_workspace * ws = gsl_eigen_nonsymmv_alloc(3);
  gsl_vector_complex * e_vals = gsl_vector_complex_alloc(3);
  gsl_matrix_complex * e_vecs = gsl_matrix_complex_alloc(3,3);

  for(i=0; i < ncount; i++)
  {
    ptmp[0] = *(pvec+3*i+0); ptmp[1] = *(pvec+3*i+1); ptmp[2] = *(pvec+3*i+2);

    for(j=1; j<=3; j++)
    {
      for(m=1; m<=3; m++)
      {
        tmp=0.;
        for(k=1; k<=3; k++)
        {
          for(n=1; n<=3; n++)
          {
            tmp += tensor_get(C, j, k, m, n)*ptmp[k-1]*ptmp[n-1];
          }
        }
        tmp = tmp/rho;
        gsl_matrix_set(D, j-1, m-1, tmp);
      }
    }
   
    /* we have the matrix now */
    gsl_eigen_nonsymmv(D, e_vals, e_vecs, ws);
    gsl_eigen_nonsymmv_sort(e_vals, e_vecs, GSL_EIGEN_SORT_ABS_ASC); /*뒤로갈수록 증가하는 순서 */

    for(j=0; j<3; j++)
    {
      eval_tmp = gsl_vector_complex_get(e_vals, j);
      if( fabs(GSL_IMAG(eval_tmp)) < 1E-5){ *(vel+3*i+j) = sqrt(GSL_REAL(eval_tmp)); }
      else{ errn = 1; goto L_00001; }
      for(k=0; k<3; k++)
      {
        evec_tmp = gsl_matrix_complex_get(e_vecs, k, j);
        if( fabs(GSL_IMAG(evec_tmp)) < 1E-5){ *(oscills+i*9+j*3+k)= GSL_REAL(evec_tmp); }
        else{errn=1; goto L_00001; }
      }
    }
    
  }

  gsl_matrix_free(D);
  gsl_matrix_complex_free(e_vecs);
  gsl_vector_complex_free(e_vals);
  gsl_eigen_nonsymmv_free(ws);

  L_00001: ;
  if(errn==1){  fprintf(stderr, "hey, the solutions are not real...! maybe this does not have sol?\n"); }
  return errn;
}

/* for convenience for unit transform */
void change_units(double scale, int nelements, double *physical)
{
  int i;
  double tmp=0.;

  for(i=0; i < nelements; i++)
  {
    tmp = *(physical+i);
    *(physical+i) = tmp*scale;
  }

  return;
}
