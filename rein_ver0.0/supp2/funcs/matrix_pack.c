#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix_pack.h"

/* some of tensor operators */
/* initializing tensors */
void tensor_set_all(double *T, double x)
{
  int i;
  for(i=0; i < tintv; i++)
  {
    *(T +i) = x;
  }

  return;
}

void tensor_set(double *T, const int i, const int j, const int k, const int l, const double x)
{
  *(T+27*(i-1)+9*(j-1)+3*(k-1)+(l-1))=x;
  return;
}

double tensor_get(const double *T, const int i, const int j, const int k, const int l)
{
  return *(T+27*(i-1)+9*(j-1)+3*(k-1)+(l-1));
}

/* only limited for square matrix */
void matrix_set(int size, double *M, int i, int j, double x)
{
  int d= size;
  *(M+d*(i-1)+(j-1)) = x;

  return;
}

double matrix_get(int size, double *M, int i, int j)
{
  double x;
  int d = size;
  x = *(M+d*(i-1)+(j-1));

  return x;
}

/* change elastic tensor to 6 x 6 matrix */
/* 11 -> 1
 *  * 22 -> 2
 *   * 33 -> 3
 *    * 32 -> 4
 *     * 31 -> 5
 *      * 21 -> 6
 *       * M should have at least 6 by 6 (36 elements) space
 *        */
void tensor_2_mat(double *T, double *M)
{
  int i, j;
  int ii, jj, kk, ll;
  double tmp=0;

  // initialize matrix
  for(i=0; i < 6; i ++)
  {
    for(j=0; j < 6; j++)
    {
      *(M+i*6+j) = 0;
    }
  }

  // get numbers from the tensor.
  for(i=0; i < 6; i++)
  {
    for(j=0; j < 6; j++)
    {
      conv_ind_op(i+1, j+1, &ii, &jj, &kk, &ll);
      tmp = tensor_get(T, ii, jj, kk, ll);
      *(M+i*6+j) = tmp;
      *(M+j*6+i) = tmp;
    }
  }

  return;
}

int conv_ind(int i, int j, int k, int l, int *a, int *b)
{
  int errn=0;

  if(i==1 && j==1){*a=1;}
  else if(i==2 && j==2){*a=2;}
  else if(i==3 && j==3){*a=3;}
  else if((i==3 && j==2) || (i==2 && j==3)){*a=4;}
  else if((i==3 && j==1) || (i==1 && j==3)){*a=5;}
  else if((i==2 && j==1) || (i==1 && j==2)){*a=6;}
  else{ fprintf(stderr, "error in your index i: %d j: %d\n", i, j); errn=1;}

  if(k==1 && l==1){*b=1;}
  else if(k==2 && l==2){*b=2;}
  else if(k==3 && l==3){*b=3;}
  else if((k==3 && l==2) || (k==2 && l==3)){*b=4;}
  else if((k==3 && l==1) || (k==1 && l==3)){*b=5;}
  else if((k==2 && l==1) || (k==1 && l==2)){*b=6;}
  else{ fprintf(stderr, "error in your index k: %d l: %d\n", k, l); errn=1;}

  return errn;
}

/* convert index in opposite way, from 6 by 6 matrix to 3x3x3x3 tensor */
int conv_ind_op(int a, int b, int *i, int *j, int *k, int *l)
{
  int errn=0;

  if(a==1){ *i=1; *j=1; }
  else if(a==2){ *i=2; *j=2; }
  else if(a==3){ *i=3; *j=3; }
  else if(a==4){ *i=3; *j=2; }
  else if(a==5){ *i=3; *j=1; }
  else if(a==6){ *i=2; *j=1; }
  else{ fprintf(stderr, "a: %d is out of range\n", a); errn=1;}

  if(b==1){ *k=1; *l=1; }
  else if(b==2){ *k=2; *l=2; }
  else if(b==3){ *k=3; *l=3; }
  else if(b==4){ *k=3; *l=2; }
  else if(b==5){ *k=3; *l=1; }
  else if(b==6){ *k=2; *l=1; }
  else{ fprintf(stderr, "b: %d is out of range\n", b); errn=1;}

  return errn;
}

/* copy double array : extraction occur for nelements in order */
void ext_array(int nelements, double *in, double *out)
{
    int i=0;
    for(i=0; i < nelements; i++)
    {
        *(out+i) = *(in+i);
    }

    return;
}

/* simple matrix operator collection */
/* w(k) = u(j)*v(k-j) */
int conv(double *u, int l1, double *v, int l2, double *w)
{
  int i, j;
  int sind, eind;

  for(i=0; i < l1+l2-1; i++)
  {
    *(w+i)=0;
    sind=0; eind=i;

    if(eind>=l2){eind=l2-1;}
    if(i>=l1){ sind = i+1-l1; }

    for(j=sind; j <= eind; j++)
    {
      *(w+i) += *(u+i-j)**(v+j);
    }
  }

  return l1+l2-1;
}

/* sum: y = x+y */
void sum_v(int dimension, double *x, double *y)
{ // only sum up for same dimension vectors
  int i;

  for(i=0; i < dimension; i++)
  {
    *(y+i) = *(x+i) + *(y+i);
  }

  return;
}

/* constant multiplier x = ax */
void amul_v(int dimension, double *x, double coeff)
{
  int i;
  for(i=0; i < dimension; i++)
  {
    *(x+i) = *(x+i)*coeff;
  }

  return;
}

/* matrix and vector multiplication
 * y = Ax */
void mat_vec_mul(double *A, double *x, int m, int n, double *y)
{
  int i, j;
  double tmp;

  for(i=0; i < m; i++)
  {
    tmp=0;
    for(j=0; j < n; j++)
    {
      tmp += *(A+i*n+j)**(x+j);
    }
    *(y+i) = tmp;
  }

  return;
}

/* complex number version */
void cmat_vec_mul(gsl_matrix_complex *A, gsl_complex x[], int m, int n, gsl_complex *y)
{
  int i, j;
  gsl_complex tmp, tmp2; //, tmp3;

  for(i=0; i < m; i++)
  {
    tmp = GSL_COMPLEX_ZERO;
    for(j=0; j < n; j++)
    {
      tmp2 = gsl_complex_mul( gsl_matrix_complex_get(A, i, j), *(x+j));
      tmp = gsl_complex_add(tmp, tmp2);
    }
    *(y+i) = tmp;
  }

  return;
}

/* matrix multiplication 
 * C = A*B 
 * A is m by n, B is n by k, C is m by k
 */
void mat_mat_mul(double A[], double B[], int m, int n, int k, double *C)
{ // matrix sizes are (m by n) and (n by k)
  int x, y, z;
  double tmp;

  for(x=0; x < m; x++)
  {
    for(y=0; y < k; y++)
    {
      tmp=0;
      for(z =0; z < n; z++)
      {
        tmp += *(A+x*n+z)**(B+z*k+y);
      }
      *(C+x*k+y) = tmp;
    }
  }

  return;
}

/* this is special case for 3 by 3 matrix and 3-dimensional vector */
void mat_rot(double *R, double *x, double *xr)
{ // xr = Rx;
  int i, j;

  for(i=0; i < 3; i++)
  {
    *(xr+i)=0;
    for(j=0; j < 3; j++)
    {
      *(xr+i) += *(R+i*3+j)**(x+j);
    }
  }

  return;
}

/* inverse rotation refer to
 * Q = Hq is original rotation if we know q and H
 * then,
 * when Q and H is known, operator to get q is called inverse rotation
 * H should be square matrix
 */
void mat_rot_inv(double *R, double *xr, double *x)
{ // x = inv(R)*xr  (in other words, xr = Rx with unknown x.
  int i, j, s;
  gsl_matrix *matR = gsl_matrix_calloc(3,3);
  gsl_vector *vecxr = gsl_vector_calloc(3);
  gsl_vector *vecx = gsl_vector_alloc(3);
  gsl_permutation *p = gsl_permutation_alloc(3);

  for(i=0; i <3; i++)
  {
    gsl_vector_set(vecxr, i, *(xr+i));
    for(j=0; j < 3; j++)
    {
      gsl_matrix_set(matR, i, j, *(R+3*i+j));
    }
  }
  gsl_linalg_LU_decomp(matR, p, &s);
  gsl_linalg_LU_solve(matR, p, vecxr, vecx);

  for(i=0; i <3; i++)
  {
    *(x+i) = gsl_vector_get(vecx, i);
  }

  gsl_vector_free(vecx); gsl_vector_free(vecxr);
  gsl_matrix_free(matR);
  gsl_permutation_free(p);

  return;
}

/* some vector operations */
double vdot(double x[], double y[], int dimension)
{
  int i;
  double result=0;

  for(i=0; i < dimension; i++)
  {
    result += *(x+i)**(y+i);
  }

  return result;
}

double vabs(double x[], int dimension)
{
  int i;
  double result=0;

  for(i=0; i < dimension; i++)
  {
    result += x[i]*x[i];
  }

  return sqrt(result);
}

double vabs2(double x[],int dimension)
{
  int i;
  double result =0;
  for(i=0; i < dimension; i++)
  {
    result += x[i]*x[i];
  }

  return result;
}

void vcross(double x[], double y[], double *vec)
{ /* only possible for 3-dimensional vector */
  
  *(vec+0) = x[1]*y[2] - x[2]*y[1];
  *(vec+1) = x[2]*y[0] - x[0]*y[2];
  *(vec+2) = x[0]*y[1] - x[1]*y[0];

  return;
}

/* inverse matrix multiplication refer to
 * Coeff = inv(L)A is where we know L and A.
 * then, it is equivalent to
 * A = L*Coeff (A, L known)
 * Then we can solve for "Coeff" as linear equation
 */
void mat_inv_solve(int dimension, double L[], double A[], double *x)
{ // x = inv(L)*A  (in other words, A = Lx with unknown x.
    int i, j, s;

    gsl_matrix *matL = gsl_matrix_calloc(dimension, dimension);
    gsl_vector *vecA = gsl_vector_calloc(dimension);
    gsl_vector *vecx = gsl_vector_alloc(dimension);
    gsl_permutation *p = gsl_permutation_alloc(dimension);

    for(i=0; i <dimension; i++)
    {
        gsl_vector_set(vecA, i, *(A+i));
        for(j=0; j < dimension; j++)
        {
            gsl_matrix_set(matL, i, j, *(L+dimension*i+j));
        }
    }
    gsl_linalg_LU_decomp(matL, p, &s);
    gsl_linalg_LU_solve(matL, p, vecA, vecx);

    for(i=0; i < dimension; i++)
    {
        *(x+i) = gsl_vector_get(vecx, i);
    }

    gsl_vector_free(vecx); gsl_vector_free(vecA);
    gsl_matrix_free(matL);
    gsl_permutation_free(p);

    return;
}

/* this is to define Christoffel's matrix 
 * from given inputs as A, B, C, D, E, 
 * also this includes packages that rotates elastic tensors from its orientation
 * azi, tilt of symmetric axis must be given 
 * structurely rotated tensor is output of this function.
 */

/* one layer christoffel. for several layers, you have to iteratively define. *
 * experessions follow form of Nagaya et al 2008
 */
void set_elastic(double rho, double alpha, double beta, double B, double C, double E, double azi, double tilt, double *o_tensor)
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

void get_aaxis_rotmat(double trend, double incli, double *M)
{
  /* inclination measured from horizontal. trend is calculated from CW from North */
  /* eta: trend, direction of hexagonal symmetric axis from North in Clockwise direction
   * lambda: dip, angle measured from vertical downward axis (x3)
   */
  double eta, lambda;
  lambda = trend*deg2rad;
  eta = (90-incli)*deg2rad;

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

void get_2d_rotmat(double angle, double *M)
{
  /* we give angle that makes x2 component of slowness into 0 for newly defined coord.
   * it is denoted as angle epsilon in Oda 2012
   * rotation is clockwise from x1 axis, and rotation axis is x3
   */

  double epsilon = angle*deg2rad;

  matrix_set(3, M, 1, 1, cos(epsilon)); /* rotation in clock wise direction */
  matrix_set(3, M, 1, 2, sin(epsilon));
  matrix_set(3, M, 1, 3, 0);
  matrix_set(3, M, 2, 1, (-1)*sin(epsilon));
  matrix_set(3, M, 2, 2, cos(epsilon));
  matrix_set(3, M, 2, 3, 0);
  matrix_set(3, M, 3, 1, 0);
  matrix_set(3, M, 3, 2, 0);
  matrix_set(3, M, 3, 3, 1);

  return;
}

double get_rothor_angle(double qvec[])
{
  double angle=0;
  angle = atan2(qvec[1], qvec[0]);

  if( angle < 0 ){ angle += M_PI*2; }
  angle = angle * rad2deg;

  return angle;
}

void get_vel_perturb(double rho, double alpha, double beta, double B, double C, double E, double *vel_p)
{ // given in order of Vp, Vsv, Vsh
  /* output: vel_p 
   * *(vel_p+0): perturbation in % for square of P wave
   * *(vel_p+1): perturbation in % for square of S wave
   */
  double angle, vp, vs, maxp, minp, maxs, mins;
  maxp=-1000; minp=1000; maxs=-1000; mins=1000;
  *(vel_p+0) = 0;
  *(vel_p+1) = 0;

  angle =0;
  while(angle <= M_PI_2)
  {
    vp = alpha*(B*cos(2*angle)+C*cos(4*angle));
    vs = beta*E*cos(2*angle);

    if(maxp < vp){ maxp = vp;}
    if(vp < minp){ minp = vp;}
    if(maxs < vs){ maxs = vs;}
    if(vs < mins){ mins = vs;}
    angle += 0.05*rad2deg;
  }

  *(vel_p+0) = 0.5*(maxp-minp);
  *(vel_p+1) = 0.5*(maxs-mins);

  return;
}

void rot_tensor(double *mrot, double *T, double *rotT)
{
  int i, j, k, l, ii, jj, kk, ll;
  double tmp=0;

  for(i=1; i<=3; i++)
  {
    for(j=1; j <=3; j++)
    {
      for(k=1; k <=3; k++)
      {
        for(l=1; l <=3; l++)
        {
          tmp=0;
          for(ii=1; ii <=3; ii++)
          {
          for(jj=1; jj<=3; jj++)
          {
          for(kk=1; kk<=3; kk++)
          {
          for(ll=1; ll<=3; ll++)
          {
            tmp += matrix_get(3, mrot, ii, i)*matrix_get(3, mrot, jj, j)*matrix_get(3, mrot, kk, k)*matrix_get(3, mrot, ll, l)*tensor_get(T, ii,jj,kk,ll);
          }
          }
          }
          }
          tensor_set(rotT, i, j, k, l, tmp);
        }
      }
    }
  }

  return;
}
