#include "auto_f2c.h"
/*==========================================================================
 *  Code generated by the Maple script autocode.mpl 
 *
 *  Dependent variables: [phi[1], phi[2], phi[3], phi[4], phi[5], phi[6], phi[7]]
 *  Parameters:          [omega]
 *  The vector field given in Maple was
 *     f[1] = omega^2*(7*sin(phi[1])+6*sin(phi[2])+5*sin(phi[3])+4*sin(phi[4])+3*sin(phi[5])+2*sin(phi[6])+sin(phi[7]))-7*tan(phi[1])
 *     f[2] = omega^2*(6*sin(phi[1])+6*sin(phi[2])+5*sin(phi[3])+4*sin(phi[4])+3*sin(phi[5])+2*sin(phi[6])+sin(phi[7]))-6*tan(phi[2])
 *     f[3] = omega^2*(5*sin(phi[1])+5*sin(phi[2])+5*sin(phi[3])+4*sin(phi[4])+3*sin(phi[5])+2*sin(phi[6])+sin(phi[7]))-5*tan(phi[3])
 *     f[4] = omega^2*(4*sin(phi[1])+4*sin(phi[2])+4*sin(phi[3])+4*sin(phi[4])+3*sin(phi[5])+2*sin(phi[6])+sin(phi[7]))-4*tan(phi[4])
 *     f[5] = omega^2*(3*sin(phi[1])+3*sin(phi[2])+3*sin(phi[3])+3*sin(phi[4])+3*sin(phi[5])+2*sin(phi[6])+sin(phi[7]))-3*tan(phi[5])
 *     f[6] = omega^2*(2*sin(phi[1])+2*sin(phi[2])+2*sin(phi[3])+2*sin(phi[4])+2*sin(phi[5])+2*sin(phi[6])+sin(phi[7]))-2*tan(phi[6])
 *     f[7] = omega^2*(sin(phi[1])+sin(phi[2])+sin(phi[3])+sin(phi[4])+sin(phi[5])+sin(phi[6])+sin(phi[7]))-tan(phi[7])
 *  Remember that Maple starts its indices at 1, while C starts at 0.
 *==========================================================================
 */

/*--------------------------------------------------------------------------*/
/* FUNC  Defines the vector field.                                          */
/*--------------------------------------------------------------------------*/

int func(integer ndim, const doublereal *u, const integer *icp, 
         const doublereal *par, integer ijac, 
         doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
      integer dfdu_dim1, dfdp_dim1;
      doublereal jacob[7][7];
      doublereal dfdpar[7][1];

      doublereal phi[7];
      doublereal omega;

      dfdu_dim1 = ndim;
      dfdp_dim1 = ndim;

      omega = par[0];
      phi[0] = u[0];
      phi[1] = u[1];
      phi[2] = u[2];
      phi[3] = u[3];
      phi[4] = u[4];
      phi[5] = u[5];
      phi[6] = u[6];

      f[0] = omega*omega*(7.0*sin(phi[0])+6.0*sin(phi[1])+5.0*sin(phi[2])+4.0*
sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]))-7.0*tan(phi[0]);
      f[1] = omega*omega*(6.0*sin(phi[0])+6.0*sin(phi[1])+5.0*sin(phi[2])+4.0*
sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]))-6.0*tan(phi[1]);
      f[2] = omega*omega*(5.0*sin(phi[0])+5.0*sin(phi[1])+5.0*sin(phi[2])+4.0*
sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]))-5.0*tan(phi[2]);
      f[3] = omega*omega*(4.0*sin(phi[0])+4.0*sin(phi[1])+4.0*sin(phi[2])+4.0*
sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]))-4.0*tan(phi[3]);
      f[4] = omega*omega*(3.0*sin(phi[0])+3.0*sin(phi[1])+3.0*sin(phi[2])+3.0*
sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]))-3.0*tan(phi[4]);
      f[5] = omega*omega*(2.0*sin(phi[0])+2.0*sin(phi[1])+2.0*sin(phi[2])+2.0*
sin(phi[3])+2.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]))-2.0*tan(phi[5]);
      f[6] = omega*omega*(sin(phi[0])+sin(phi[1])+sin(phi[2])+sin(phi[3])+sin(
phi[4])+sin(phi[5])+sin(phi[6]))-tan(phi[6]);

      if (ijac == 0) {
          return 0;
      }

      /* Compute the Jacobian matrix. */

      jacob[0][0] = 7.0*omega*omega*cos(phi[0])-7.0-7.0*pow(tan(phi[0]),2.0);
      jacob[0][1] = 6.0*omega*omega*cos(phi[1]);
      jacob[0][2] = 5.0*omega*omega*cos(phi[2]);
      jacob[0][3] = 4.0*omega*omega*cos(phi[3]);
      jacob[0][4] = 3.0*omega*omega*cos(phi[4]);
      jacob[0][5] = 2.0*omega*omega*cos(phi[5]);
      jacob[0][6] = omega*omega*cos(phi[6]);
      jacob[1][0] = 6.0*omega*omega*cos(phi[0]);
      jacob[1][1] = 6.0*omega*omega*cos(phi[1])-6.0-6.0*pow(tan(phi[1]),2.0);
      jacob[1][2] = 5.0*omega*omega*cos(phi[2]);
      jacob[1][3] = 4.0*omega*omega*cos(phi[3]);
      jacob[1][4] = 3.0*omega*omega*cos(phi[4]);
      jacob[1][5] = 2.0*omega*omega*cos(phi[5]);
      jacob[1][6] = omega*omega*cos(phi[6]);
      jacob[2][0] = 5.0*omega*omega*cos(phi[0]);
      jacob[2][1] = 5.0*omega*omega*cos(phi[1]);
      jacob[2][2] = 5.0*omega*omega*cos(phi[2])-5.0-5.0*pow(tan(phi[2]),2.0);
      jacob[2][3] = 4.0*omega*omega*cos(phi[3]);
      jacob[2][4] = 3.0*omega*omega*cos(phi[4]);
      jacob[2][5] = 2.0*omega*omega*cos(phi[5]);
      jacob[2][6] = omega*omega*cos(phi[6]);
      jacob[3][0] = 4.0*omega*omega*cos(phi[0]);
      jacob[3][1] = 4.0*omega*omega*cos(phi[1]);
      jacob[3][2] = 4.0*omega*omega*cos(phi[2]);
      jacob[3][3] = 4.0*omega*omega*cos(phi[3])-4.0-4.0*pow(tan(phi[3]),2.0);
      jacob[3][4] = 3.0*omega*omega*cos(phi[4]);
      jacob[3][5] = 2.0*omega*omega*cos(phi[5]);
      jacob[3][6] = omega*omega*cos(phi[6]);
      jacob[4][0] = 3.0*omega*omega*cos(phi[0]);
      jacob[4][1] = 3.0*omega*omega*cos(phi[1]);
      jacob[4][2] = 3.0*omega*omega*cos(phi[2]);
      jacob[4][3] = 3.0*omega*omega*cos(phi[3]);
      jacob[4][4] = 3.0*omega*omega*cos(phi[4])-3.0-3.0*pow(tan(phi[4]),2.0);
      jacob[4][5] = 2.0*omega*omega*cos(phi[5]);
      jacob[4][6] = omega*omega*cos(phi[6]);
      jacob[5][0] = 2.0*omega*omega*cos(phi[0]);
      jacob[5][1] = 2.0*omega*omega*cos(phi[1]);
      jacob[5][2] = 2.0*omega*omega*cos(phi[2]);
      jacob[5][3] = 2.0*omega*omega*cos(phi[3]);
      jacob[5][4] = 2.0*omega*omega*cos(phi[4]);
      jacob[5][5] = 2.0*omega*omega*cos(phi[5])-2.0-2.0*pow(tan(phi[5]),2.0);
      jacob[5][6] = omega*omega*cos(phi[6]);
      jacob[6][0] = omega*omega*cos(phi[0]);
      jacob[6][1] = omega*omega*cos(phi[1]);
      jacob[6][2] = omega*omega*cos(phi[2]);
      jacob[6][3] = omega*omega*cos(phi[3]);
      jacob[6][4] = omega*omega*cos(phi[4]);
      jacob[6][5] = omega*omega*cos(phi[5]);
      jacob[6][6] = omega*omega*cos(phi[6])-1.0-pow(tan(phi[6]),2.0);
      ARRAY2D(dfdu,0,0) = jacob[0][0];
      ARRAY2D(dfdu,0,1) = jacob[0][1];
      ARRAY2D(dfdu,0,2) = jacob[0][2];
      ARRAY2D(dfdu,0,3) = jacob[0][3];
      ARRAY2D(dfdu,0,4) = jacob[0][4];
      ARRAY2D(dfdu,0,5) = jacob[0][5];
      ARRAY2D(dfdu,0,6) = jacob[0][6];
      ARRAY2D(dfdu,1,0) = jacob[1][0];
      ARRAY2D(dfdu,1,1) = jacob[1][1];
      ARRAY2D(dfdu,1,2) = jacob[1][2];
      ARRAY2D(dfdu,1,3) = jacob[1][3];
      ARRAY2D(dfdu,1,4) = jacob[1][4];
      ARRAY2D(dfdu,1,5) = jacob[1][5];
      ARRAY2D(dfdu,1,6) = jacob[1][6];
      ARRAY2D(dfdu,2,0) = jacob[2][0];
      ARRAY2D(dfdu,2,1) = jacob[2][1];
      ARRAY2D(dfdu,2,2) = jacob[2][2];
      ARRAY2D(dfdu,2,3) = jacob[2][3];
      ARRAY2D(dfdu,2,4) = jacob[2][4];
      ARRAY2D(dfdu,2,5) = jacob[2][5];
      ARRAY2D(dfdu,2,6) = jacob[2][6];
      ARRAY2D(dfdu,3,0) = jacob[3][0];
      ARRAY2D(dfdu,3,1) = jacob[3][1];
      ARRAY2D(dfdu,3,2) = jacob[3][2];
      ARRAY2D(dfdu,3,3) = jacob[3][3];
      ARRAY2D(dfdu,3,4) = jacob[3][4];
      ARRAY2D(dfdu,3,5) = jacob[3][5];
      ARRAY2D(dfdu,3,6) = jacob[3][6];
      ARRAY2D(dfdu,4,0) = jacob[4][0];
      ARRAY2D(dfdu,4,1) = jacob[4][1];
      ARRAY2D(dfdu,4,2) = jacob[4][2];
      ARRAY2D(dfdu,4,3) = jacob[4][3];
      ARRAY2D(dfdu,4,4) = jacob[4][4];
      ARRAY2D(dfdu,4,5) = jacob[4][5];
      ARRAY2D(dfdu,4,6) = jacob[4][6];
      ARRAY2D(dfdu,5,0) = jacob[5][0];
      ARRAY2D(dfdu,5,1) = jacob[5][1];
      ARRAY2D(dfdu,5,2) = jacob[5][2];
      ARRAY2D(dfdu,5,3) = jacob[5][3];
      ARRAY2D(dfdu,5,4) = jacob[5][4];
      ARRAY2D(dfdu,5,5) = jacob[5][5];
      ARRAY2D(dfdu,5,6) = jacob[5][6];
      ARRAY2D(dfdu,6,0) = jacob[6][0];
      ARRAY2D(dfdu,6,1) = jacob[6][1];
      ARRAY2D(dfdu,6,2) = jacob[6][2];
      ARRAY2D(dfdu,6,3) = jacob[6][3];
      ARRAY2D(dfdu,6,4) = jacob[6][4];
      ARRAY2D(dfdu,6,5) = jacob[6][5];
      ARRAY2D(dfdu,6,6) = jacob[6][6];

      if (ijac == 1) {
          return 0;
      }

      /* Compute the derivative of the vector field  */
      /* with respect to the parameters.             */

      dfdpar[0][0] = 2.0*omega*(7.0*sin(phi[0])+6.0*sin(phi[1])+5.0*sin(phi[2])
+4.0*sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]));
      dfdpar[1][0] = 2.0*omega*(6.0*sin(phi[0])+6.0*sin(phi[1])+5.0*sin(phi[2])
+4.0*sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]));
      dfdpar[2][0] = 2.0*omega*(5.0*sin(phi[0])+5.0*sin(phi[1])+5.0*sin(phi[2])
+4.0*sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]));
      dfdpar[3][0] = 2.0*omega*(4.0*sin(phi[0])+4.0*sin(phi[1])+4.0*sin(phi[2])
+4.0*sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]));
      dfdpar[4][0] = 2.0*omega*(3.0*sin(phi[0])+3.0*sin(phi[1])+3.0*sin(phi[2])
+3.0*sin(phi[3])+3.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]));
      dfdpar[5][0] = 2.0*omega*(2.0*sin(phi[0])+2.0*sin(phi[1])+2.0*sin(phi[2])
+2.0*sin(phi[3])+2.0*sin(phi[4])+2.0*sin(phi[5])+sin(phi[6]));
      dfdpar[6][0] = 2.0*omega*(sin(phi[0])+sin(phi[1])+sin(phi[2])+sin(phi[3])
+sin(phi[4])+sin(phi[5])+sin(phi[6]));
      ARRAY2D(dfdp,0,0) = dfdpar[0][0];
      ARRAY2D(dfdp,1,0) = dfdpar[1][0];
      ARRAY2D(dfdp,2,0) = dfdpar[2][0];
      ARRAY2D(dfdp,3,0) = dfdpar[3][0];
      ARRAY2D(dfdp,4,0) = dfdpar[4][0];
      ARRAY2D(dfdp,5,0) = dfdpar[5][0];
      ARRAY2D(dfdp,6,0) = dfdpar[6][0];

      return 0;
}

/*--------------------------------------------------------------------------*/
/* STPNT  Gives a starting pointing                                         */
/*--------------------------------------------------------------------------*/

int stpnt(integer ndim, doublereal t,
          doublereal *u, doublereal *par)
{
      doublereal phi[7];
      doublereal omega;

      /* Change these values from 0.0 to the correct values! */

      /* Parameter values of the starting point */
      omega = 0.0;

      /* The starting point. */
      phi[0] = 0.0;
      phi[1] = 0.0;
      phi[2] = 0.0;
      phi[3] = 0.0;
      phi[4] = 0.0;
      phi[5] = 0.0;
      phi[6] = 0.0;

      par[0] = omega;
      u[0] = phi[0];
      u[1] = phi[1];
      u[2] = phi[2];
      u[3] = phi[3];
      u[4] = phi[4];
      u[5] = phi[5];
      u[6] = phi[6];
      return 0;
}

/*--------------------------------------------------------------------------*/
/* BCND  Defines the boundary conditions                                    */
/*--------------------------------------------------------------------------*/

int bcnd(integer ndim, const doublereal *par, const integer *icp,
         integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
         doublereal *fb, doublereal *dbc)
{
      return 0;
}

/*--------------------------------------------------------------------------*/
/* ICND  Defines the integral conditions                                    */
/*--------------------------------------------------------------------------*/

int icnd(integer ndim, const doublereal *par, const integer *icp,
         integer nint, const doublereal *u, const doublereal *uold,
         const doublereal *udot, const doublereal *upold, integer ijac,
         doublereal *fi, doublereal *dint)
{
      return 0;
}

/*--------------------------------------------------------------------------*/
/* FOPT                                                                     */
/*--------------------------------------------------------------------------*/

int fopt(integer ndim, const doublereal *u, const integer *icp,
         const doublereal *par, integer ijac,
         doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
      return 0;
}

/*--------------------------------------------------------------------------*/
/* PVLS                                                                     */
/*--------------------------------------------------------------------------*/

int pvls(integer ndim, const doublereal *u,
         doublereal *par)
{
      return 0;
}

/*--------------------------------------------------------------------------*/
