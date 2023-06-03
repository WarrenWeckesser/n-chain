#include "auto_f2c.h"
/*==========================================================================
 *  Dependent variables: [omega, phi[2], phi[3], phi[4], phi[5], phi[6], phi[7]]
 *  Parameters:          [phi[1]]
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

      omega = u[0];
      phi[0] = par[0];
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
      omega = 4.40405811295248 ;
      /*
       *  0.43936735946172
       *  1.01324473615173
       *  1.60245959229890
       *  2.21367411434621
       *  2.86044637155862
       *  3.56849832447737
       *  4.40405811295248
       */
      /* The starting point. */
      phi[0] = 0.000001;
      phi[1] = 0.0;
      phi[2] = 0.0;
      phi[3] = 0.0;
      phi[4] = 0.0;
      phi[5] = 0.0;
      phi[6] = 0.0;

      u[0]   = omega;
      par[0] = phi[0];
      u[1]   = phi[1];
      u[2]   = phi[2];
      u[3]   = phi[3];
      u[4]   = phi[4];
      u[5]   = phi[5];
      u[6]   = phi[6];
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
      par[1] = u[0];
      return 0;
}

/*--------------------------------------------------------------------------*/
