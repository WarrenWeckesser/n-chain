#include "auto_f2c.h"
/*==========================================================================
 *  Code generated by the Maple script autocode.mpl 
 *
 *  Dependent variables: [phi[1], phi[2], phi[3], phi[4], phi[5], phi[6], phi[7], phi[8], phi[9], phi[10], phi[11], phi[12], phi[13], phi[14], phi[15]]
 *  Parameters:          [omega]
 *  The vector field given in Maple was
 *     f[1] = omega^2*(15*sin(phi[1])+14*sin(phi[2])+13*sin(phi[3])+12*sin(phi[4])+11*sin(phi[5])+10*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-15*tan(phi[1])
 *     f[2] = omega^2*(14*sin(phi[1])+14*sin(phi[2])+13*sin(phi[3])+12*sin(phi[4])+11*sin(phi[5])+10*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-14*tan(phi[2])
 *     f[3] = omega^2*(13*sin(phi[1])+13*sin(phi[2])+13*sin(phi[3])+12*sin(phi[4])+11*sin(phi[5])+10*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-13*tan(phi[3])
 *     f[4] = omega^2*(12*sin(phi[1])+12*sin(phi[2])+12*sin(phi[3])+12*sin(phi[4])+11*sin(phi[5])+10*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-12*tan(phi[4])
 *     f[5] = omega^2*(11*sin(phi[1])+11*sin(phi[2])+11*sin(phi[3])+11*sin(phi[4])+11*sin(phi[5])+10*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-11*tan(phi[5])
 *     f[6] = omega^2*(10*sin(phi[1])+10*sin(phi[2])+10*sin(phi[3])+10*sin(phi[4])+10*sin(phi[5])+10*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-10*tan(phi[6])
 *     f[7] = omega^2*(9*sin(phi[1])+9*sin(phi[2])+9*sin(phi[3])+9*sin(phi[4])+9*sin(phi[5])+9*sin(phi[6])+9*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-9*tan(phi[7])
 *     f[8] = omega^2*(8*sin(phi[1])+8*sin(phi[2])+8*sin(phi[3])+8*sin(phi[4])+8*sin(phi[5])+8*sin(phi[6])+8*sin(phi[7])+8*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-8*tan(phi[8])
 *     f[9] = omega^2*(7*sin(phi[1])+7*sin(phi[2])+7*sin(phi[3])+7*sin(phi[4])+7*sin(phi[5])+7*sin(phi[6])+7*sin(phi[7])+7*sin(phi[8])+7*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-7*tan(phi[9])
 *     f[10] = omega^2*(6*sin(phi[1])+6*sin(phi[2])+6*sin(phi[3])+6*sin(phi[4])+6*sin(phi[5])+6*sin(phi[6])+6*sin(phi[7])+6*sin(phi[8])+6*sin(phi[9])+6*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-6*tan(phi[10])
 *     f[11] = omega^2*(5*sin(phi[1])+5*sin(phi[2])+5*sin(phi[3])+5*sin(phi[4])+5*sin(phi[5])+5*sin(phi[6])+5*sin(phi[7])+5*sin(phi[8])+5*sin(phi[9])+5*sin(phi[10])+5*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-5*tan(phi[11])
 *     f[12] = omega^2*(4*sin(phi[1])+4*sin(phi[2])+4*sin(phi[3])+4*sin(phi[4])+4*sin(phi[5])+4*sin(phi[6])+4*sin(phi[7])+4*sin(phi[8])+4*sin(phi[9])+4*sin(phi[10])+4*sin(phi[11])+4*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-4*tan(phi[12])
 *     f[13] = omega^2*(3*sin(phi[1])+3*sin(phi[2])+3*sin(phi[3])+3*sin(phi[4])+3*sin(phi[5])+3*sin(phi[6])+3*sin(phi[7])+3*sin(phi[8])+3*sin(phi[9])+3*sin(phi[10])+3*sin(phi[11])+3*sin(phi[12])+3*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-3*tan(phi[13])
 *     f[14] = omega^2*(2*sin(phi[1])+2*sin(phi[2])+2*sin(phi[3])+2*sin(phi[4])+2*sin(phi[5])+2*sin(phi[6])+2*sin(phi[7])+2*sin(phi[8])+2*sin(phi[9])+2*sin(phi[10])+2*sin(phi[11])+2*sin(phi[12])+2*sin(phi[13])+2*sin(phi[14])+sin(phi[15]))-2*tan(phi[14])
 *     f[15] = omega^2*(sin(phi[1])+sin(phi[2])+sin(phi[3])+sin(phi[4])+sin(phi[5])+sin(phi[6])+sin(phi[7])+sin(phi[8])+sin(phi[9])+sin(phi[10])+sin(phi[11])+sin(phi[12])+sin(phi[13])+sin(phi[14])+sin(phi[15]))-tan(phi[15])
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
      doublereal jacob[15][15];
      doublereal dfdpar[15][1];

      doublereal phi[15];
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
      phi[7] = u[7];
      phi[8] = u[8];
      phi[9] = u[9];
      phi[10] = u[10];
      phi[11] = u[11];
      phi[12] = u[12];
      phi[13] = u[13];
      phi[14] = u[14];

      f[0] = omega*omega*(15.0*sin(phi[0])+14.0*sin(phi[1])+13.0*sin(phi[2])+
12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi
[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]))-15.0*tan(phi[0]);
      f[1] = omega*omega*(14.0*sin(phi[0])+14.0*sin(phi[1])+13.0*sin(phi[2])+
12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi
[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]))-14.0*tan(phi[1]);
      f[2] = omega*omega*(13.0*sin(phi[0])+13.0*sin(phi[1])+13.0*sin(phi[2])+
12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi
[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]))-13.0*tan(phi[2]);
      f[3] = omega*omega*(12.0*sin(phi[0])+12.0*sin(phi[1])+12.0*sin(phi[2])+
12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi
[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]))-12.0*tan(phi[3]);
      f[4] = omega*omega*(11.0*sin(phi[0])+11.0*sin(phi[1])+11.0*sin(phi[2])+
11.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi
[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]))-11.0*tan(phi[4]);
      f[5] = omega*omega*(10.0*sin(phi[0])+10.0*sin(phi[1])+10.0*sin(phi[2])+
10.0*sin(phi[3])+10.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi
[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]))-10.0*tan(phi[5]);
      f[6] = omega*omega*(9.0*sin(phi[0])+9.0*sin(phi[1])+9.0*sin(phi[2])+9.0*
sin(phi[3])+9.0*sin(phi[4])+9.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi[7])+7.0
*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-9.0*tan(phi[6]);
      f[7] = omega*omega*(8.0*sin(phi[0])+8.0*sin(phi[1])+8.0*sin(phi[2])+8.0*
sin(phi[3])+8.0*sin(phi[4])+8.0*sin(phi[5])+8.0*sin(phi[6])+8.0*sin(phi[7])+7.0
*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-8.0*tan(phi[7]);
      f[8] = omega*omega*(7.0*sin(phi[0])+7.0*sin(phi[1])+7.0*sin(phi[2])+7.0*
sin(phi[3])+7.0*sin(phi[4])+7.0*sin(phi[5])+7.0*sin(phi[6])+7.0*sin(phi[7])+7.0
*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-7.0*tan(phi[8]);
      f[9] = omega*omega*(6.0*sin(phi[0])+6.0*sin(phi[1])+6.0*sin(phi[2])+6.0*
sin(phi[3])+6.0*sin(phi[4])+6.0*sin(phi[5])+6.0*sin(phi[6])+6.0*sin(phi[7])+6.0
*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-6.0*tan(phi[9]);
      f[10] = omega*omega*(5.0*sin(phi[0])+5.0*sin(phi[1])+5.0*sin(phi[2])+5.0*
sin(phi[3])+5.0*sin(phi[4])+5.0*sin(phi[5])+5.0*sin(phi[6])+5.0*sin(phi[7])+5.0
*sin(phi[8])+5.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-5.0*tan(phi[10]);
      f[11] = omega*omega*(4.0*sin(phi[0])+4.0*sin(phi[1])+4.0*sin(phi[2])+4.0*
sin(phi[3])+4.0*sin(phi[4])+4.0*sin(phi[5])+4.0*sin(phi[6])+4.0*sin(phi[7])+4.0
*sin(phi[8])+4.0*sin(phi[9])+4.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-4.0*tan(phi[11]);
      f[12] = omega*omega*(3.0*sin(phi[0])+3.0*sin(phi[1])+3.0*sin(phi[2])+3.0*
sin(phi[3])+3.0*sin(phi[4])+3.0*sin(phi[5])+3.0*sin(phi[6])+3.0*sin(phi[7])+3.0
*sin(phi[8])+3.0*sin(phi[9])+3.0*sin(phi[10])+3.0*sin(phi[11])+3.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-3.0*tan(phi[12]);
      f[13] = omega*omega*(2.0*sin(phi[0])+2.0*sin(phi[1])+2.0*sin(phi[2])+2.0*
sin(phi[3])+2.0*sin(phi[4])+2.0*sin(phi[5])+2.0*sin(phi[6])+2.0*sin(phi[7])+2.0
*sin(phi[8])+2.0*sin(phi[9])+2.0*sin(phi[10])+2.0*sin(phi[11])+2.0*sin(phi[12])
+2.0*sin(phi[13])+sin(phi[14]))-2.0*tan(phi[13]);
      f[14] = omega*omega*(sin(phi[0])+sin(phi[1])+sin(phi[2])+sin(phi[3])+sin(
phi[4])+sin(phi[5])+sin(phi[6])+sin(phi[7])+sin(phi[8])+sin(phi[9])+sin(phi[10]
)+sin(phi[11])+sin(phi[12])+sin(phi[13])+sin(phi[14]))-tan(phi[14]);

      if (ijac == 0) {
          return 0;
      }

      /* Compute the Jacobian matrix. */

      jacob[0][0] = 15.0*omega*omega*cos(phi[0])-15.0-15.0*pow(tan(phi[0]),2.0)
;
      jacob[0][1] = 14.0*omega*omega*cos(phi[1]);
      jacob[0][2] = 13.0*omega*omega*cos(phi[2]);
      jacob[0][3] = 12.0*omega*omega*cos(phi[3]);
      jacob[0][4] = 11.0*omega*omega*cos(phi[4]);
      jacob[0][5] = 10.0*omega*omega*cos(phi[5]);
      jacob[0][6] = 9.0*omega*omega*cos(phi[6]);
      jacob[0][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[0][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[0][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[0][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[0][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[0][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[0][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[0][14] = omega*omega*cos(phi[14]);
      jacob[1][0] = 14.0*omega*omega*cos(phi[0]);
      jacob[1][1] = 14.0*omega*omega*cos(phi[1])-14.0-14.0*pow(tan(phi[1]),2.0)
;
      jacob[1][2] = 13.0*omega*omega*cos(phi[2]);
      jacob[1][3] = 12.0*omega*omega*cos(phi[3]);
      jacob[1][4] = 11.0*omega*omega*cos(phi[4]);
      jacob[1][5] = 10.0*omega*omega*cos(phi[5]);
      jacob[1][6] = 9.0*omega*omega*cos(phi[6]);
      jacob[1][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[1][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[1][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[1][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[1][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[1][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[1][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[1][14] = omega*omega*cos(phi[14]);
      jacob[2][0] = 13.0*omega*omega*cos(phi[0]);
      jacob[2][1] = 13.0*omega*omega*cos(phi[1]);
      jacob[2][2] = 13.0*omega*omega*cos(phi[2])-13.0-13.0*pow(tan(phi[2]),2.0)
;
      jacob[2][3] = 12.0*omega*omega*cos(phi[3]);
      jacob[2][4] = 11.0*omega*omega*cos(phi[4]);
      jacob[2][5] = 10.0*omega*omega*cos(phi[5]);
      jacob[2][6] = 9.0*omega*omega*cos(phi[6]);
      jacob[2][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[2][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[2][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[2][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[2][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[2][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[2][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[2][14] = omega*omega*cos(phi[14]);
      jacob[3][0] = 12.0*omega*omega*cos(phi[0]);
      jacob[3][1] = 12.0*omega*omega*cos(phi[1]);
      jacob[3][2] = 12.0*omega*omega*cos(phi[2]);
      jacob[3][3] = 12.0*omega*omega*cos(phi[3])-12.0-12.0*pow(tan(phi[3]),2.0)
;
      jacob[3][4] = 11.0*omega*omega*cos(phi[4]);
      jacob[3][5] = 10.0*omega*omega*cos(phi[5]);
      jacob[3][6] = 9.0*omega*omega*cos(phi[6]);
      jacob[3][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[3][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[3][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[3][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[3][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[3][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[3][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[3][14] = omega*omega*cos(phi[14]);
      jacob[4][0] = 11.0*omega*omega*cos(phi[0]);
      jacob[4][1] = 11.0*omega*omega*cos(phi[1]);
      jacob[4][2] = 11.0*omega*omega*cos(phi[2]);
      jacob[4][3] = 11.0*omega*omega*cos(phi[3]);
      jacob[4][4] = 11.0*omega*omega*cos(phi[4])-11.0-11.0*pow(tan(phi[4]),2.0)
;
      jacob[4][5] = 10.0*omega*omega*cos(phi[5]);
      jacob[4][6] = 9.0*omega*omega*cos(phi[6]);
      jacob[4][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[4][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[4][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[4][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[4][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[4][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[4][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[4][14] = omega*omega*cos(phi[14]);
      jacob[5][0] = 10.0*omega*omega*cos(phi[0]);
      jacob[5][1] = 10.0*omega*omega*cos(phi[1]);
      jacob[5][2] = 10.0*omega*omega*cos(phi[2]);
      jacob[5][3] = 10.0*omega*omega*cos(phi[3]);
      jacob[5][4] = 10.0*omega*omega*cos(phi[4]);
      jacob[5][5] = 10.0*omega*omega*cos(phi[5])-10.0-10.0*pow(tan(phi[5]),2.0)
;
      jacob[5][6] = 9.0*omega*omega*cos(phi[6]);
      jacob[5][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[5][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[5][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[5][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[5][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[5][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[5][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[5][14] = omega*omega*cos(phi[14]);
      jacob[6][0] = 9.0*omega*omega*cos(phi[0]);
      jacob[6][1] = 9.0*omega*omega*cos(phi[1]);
      jacob[6][2] = 9.0*omega*omega*cos(phi[2]);
      jacob[6][3] = 9.0*omega*omega*cos(phi[3]);
      jacob[6][4] = 9.0*omega*omega*cos(phi[4]);
      jacob[6][5] = 9.0*omega*omega*cos(phi[5]);
      jacob[6][6] = 9.0*omega*omega*cos(phi[6])-9.0-9.0*pow(tan(phi[6]),2.0);
      jacob[6][7] = 8.0*omega*omega*cos(phi[7]);
      jacob[6][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[6][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[6][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[6][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[6][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[6][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[6][14] = omega*omega*cos(phi[14]);
      jacob[7][0] = 8.0*omega*omega*cos(phi[0]);
      jacob[7][1] = 8.0*omega*omega*cos(phi[1]);
      jacob[7][2] = 8.0*omega*omega*cos(phi[2]);
      jacob[7][3] = 8.0*omega*omega*cos(phi[3]);
      jacob[7][4] = 8.0*omega*omega*cos(phi[4]);
      jacob[7][5] = 8.0*omega*omega*cos(phi[5]);
      jacob[7][6] = 8.0*omega*omega*cos(phi[6]);
      jacob[7][7] = 8.0*omega*omega*cos(phi[7])-8.0-8.0*pow(tan(phi[7]),2.0);
      jacob[7][8] = 7.0*omega*omega*cos(phi[8]);
      jacob[7][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[7][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[7][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[7][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[7][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[7][14] = omega*omega*cos(phi[14]);
      jacob[8][0] = 7.0*omega*omega*cos(phi[0]);
      jacob[8][1] = 7.0*omega*omega*cos(phi[1]);
      jacob[8][2] = 7.0*omega*omega*cos(phi[2]);
      jacob[8][3] = 7.0*omega*omega*cos(phi[3]);
      jacob[8][4] = 7.0*omega*omega*cos(phi[4]);
      jacob[8][5] = 7.0*omega*omega*cos(phi[5]);
      jacob[8][6] = 7.0*omega*omega*cos(phi[6]);
      jacob[8][7] = 7.0*omega*omega*cos(phi[7]);
      jacob[8][8] = 7.0*omega*omega*cos(phi[8])-7.0-7.0*pow(tan(phi[8]),2.0);
      jacob[8][9] = 6.0*omega*omega*cos(phi[9]);
      jacob[8][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[8][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[8][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[8][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[8][14] = omega*omega*cos(phi[14]);
      jacob[9][0] = 6.0*omega*omega*cos(phi[0]);
      jacob[9][1] = 6.0*omega*omega*cos(phi[1]);
      jacob[9][2] = 6.0*omega*omega*cos(phi[2]);
      jacob[9][3] = 6.0*omega*omega*cos(phi[3]);
      jacob[9][4] = 6.0*omega*omega*cos(phi[4]);
      jacob[9][5] = 6.0*omega*omega*cos(phi[5]);
      jacob[9][6] = 6.0*omega*omega*cos(phi[6]);
      jacob[9][7] = 6.0*omega*omega*cos(phi[7]);
      jacob[9][8] = 6.0*omega*omega*cos(phi[8]);
      jacob[9][9] = 6.0*omega*omega*cos(phi[9])-6.0-6.0*pow(tan(phi[9]),2.0);
      jacob[9][10] = 5.0*omega*omega*cos(phi[10]);
      jacob[9][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[9][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[9][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[9][14] = omega*omega*cos(phi[14]);
      jacob[10][0] = 5.0*omega*omega*cos(phi[0]);
      jacob[10][1] = 5.0*omega*omega*cos(phi[1]);
      jacob[10][2] = 5.0*omega*omega*cos(phi[2]);
      jacob[10][3] = 5.0*omega*omega*cos(phi[3]);
      jacob[10][4] = 5.0*omega*omega*cos(phi[4]);
      jacob[10][5] = 5.0*omega*omega*cos(phi[5]);
      jacob[10][6] = 5.0*omega*omega*cos(phi[6]);
      jacob[10][7] = 5.0*omega*omega*cos(phi[7]);
      jacob[10][8] = 5.0*omega*omega*cos(phi[8]);
      jacob[10][9] = 5.0*omega*omega*cos(phi[9]);
      jacob[10][10] = 5.0*omega*omega*cos(phi[10])-5.0-5.0*pow(tan(phi[10]),2.0
);
      jacob[10][11] = 4.0*omega*omega*cos(phi[11]);
      jacob[10][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[10][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[10][14] = omega*omega*cos(phi[14]);
      jacob[11][0] = 4.0*omega*omega*cos(phi[0]);
      jacob[11][1] = 4.0*omega*omega*cos(phi[1]);
      jacob[11][2] = 4.0*omega*omega*cos(phi[2]);
      jacob[11][3] = 4.0*omega*omega*cos(phi[3]);
      jacob[11][4] = 4.0*omega*omega*cos(phi[4]);
      jacob[11][5] = 4.0*omega*omega*cos(phi[5]);
      jacob[11][6] = 4.0*omega*omega*cos(phi[6]);
      jacob[11][7] = 4.0*omega*omega*cos(phi[7]);
      jacob[11][8] = 4.0*omega*omega*cos(phi[8]);
      jacob[11][9] = 4.0*omega*omega*cos(phi[9]);
      jacob[11][10] = 4.0*omega*omega*cos(phi[10]);
      jacob[11][11] = 4.0*omega*omega*cos(phi[11])-4.0-4.0*pow(tan(phi[11]),2.0
);
      jacob[11][12] = 3.0*omega*omega*cos(phi[12]);
      jacob[11][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[11][14] = omega*omega*cos(phi[14]);
      jacob[12][0] = 3.0*omega*omega*cos(phi[0]);
      jacob[12][1] = 3.0*omega*omega*cos(phi[1]);
      jacob[12][2] = 3.0*omega*omega*cos(phi[2]);
      jacob[12][3] = 3.0*omega*omega*cos(phi[3]);
      jacob[12][4] = 3.0*omega*omega*cos(phi[4]);
      jacob[12][5] = 3.0*omega*omega*cos(phi[5]);
      jacob[12][6] = 3.0*omega*omega*cos(phi[6]);
      jacob[12][7] = 3.0*omega*omega*cos(phi[7]);
      jacob[12][8] = 3.0*omega*omega*cos(phi[8]);
      jacob[12][9] = 3.0*omega*omega*cos(phi[9]);
      jacob[12][10] = 3.0*omega*omega*cos(phi[10]);
      jacob[12][11] = 3.0*omega*omega*cos(phi[11]);
      jacob[12][12] = 3.0*omega*omega*cos(phi[12])-3.0-3.0*pow(tan(phi[12]),2.0
);
      jacob[12][13] = 2.0*omega*omega*cos(phi[13]);
      jacob[12][14] = omega*omega*cos(phi[14]);
      jacob[13][0] = 2.0*omega*omega*cos(phi[0]);
      jacob[13][1] = 2.0*omega*omega*cos(phi[1]);
      jacob[13][2] = 2.0*omega*omega*cos(phi[2]);
      jacob[13][3] = 2.0*omega*omega*cos(phi[3]);
      jacob[13][4] = 2.0*omega*omega*cos(phi[4]);
      jacob[13][5] = 2.0*omega*omega*cos(phi[5]);
      jacob[13][6] = 2.0*omega*omega*cos(phi[6]);
      jacob[13][7] = 2.0*omega*omega*cos(phi[7]);
      jacob[13][8] = 2.0*omega*omega*cos(phi[8]);
      jacob[13][9] = 2.0*omega*omega*cos(phi[9]);
      jacob[13][10] = 2.0*omega*omega*cos(phi[10]);
      jacob[13][11] = 2.0*omega*omega*cos(phi[11]);
      jacob[13][12] = 2.0*omega*omega*cos(phi[12]);
      jacob[13][13] = 2.0*omega*omega*cos(phi[13])-2.0-2.0*pow(tan(phi[13]),2.0
);
      jacob[13][14] = omega*omega*cos(phi[14]);
      jacob[14][0] = omega*omega*cos(phi[0]);
      jacob[14][1] = omega*omega*cos(phi[1]);
      jacob[14][2] = omega*omega*cos(phi[2]);
      jacob[14][3] = omega*omega*cos(phi[3]);
      jacob[14][4] = omega*omega*cos(phi[4]);
      jacob[14][5] = omega*omega*cos(phi[5]);
      jacob[14][6] = omega*omega*cos(phi[6]);
      jacob[14][7] = omega*omega*cos(phi[7]);
      jacob[14][8] = omega*omega*cos(phi[8]);
      jacob[14][9] = omega*omega*cos(phi[9]);
      jacob[14][10] = omega*omega*cos(phi[10]);
      jacob[14][11] = omega*omega*cos(phi[11]);
      jacob[14][12] = omega*omega*cos(phi[12]);
      jacob[14][13] = omega*omega*cos(phi[13]);
      jacob[14][14] = omega*omega*cos(phi[14])-1.0-pow(tan(phi[14]),2.0);
      ARRAY2D(dfdu,0,0) = jacob[0][0];
      ARRAY2D(dfdu,0,1) = jacob[0][1];
      ARRAY2D(dfdu,0,2) = jacob[0][2];
      ARRAY2D(dfdu,0,3) = jacob[0][3];
      ARRAY2D(dfdu,0,4) = jacob[0][4];
      ARRAY2D(dfdu,0,5) = jacob[0][5];
      ARRAY2D(dfdu,0,6) = jacob[0][6];
      ARRAY2D(dfdu,0,7) = jacob[0][7];
      ARRAY2D(dfdu,0,8) = jacob[0][8];
      ARRAY2D(dfdu,0,9) = jacob[0][9];
      ARRAY2D(dfdu,0,10) = jacob[0][10];
      ARRAY2D(dfdu,0,11) = jacob[0][11];
      ARRAY2D(dfdu,0,12) = jacob[0][12];
      ARRAY2D(dfdu,0,13) = jacob[0][13];
      ARRAY2D(dfdu,0,14) = jacob[0][14];
      ARRAY2D(dfdu,1,0) = jacob[1][0];
      ARRAY2D(dfdu,1,1) = jacob[1][1];
      ARRAY2D(dfdu,1,2) = jacob[1][2];
      ARRAY2D(dfdu,1,3) = jacob[1][3];
      ARRAY2D(dfdu,1,4) = jacob[1][4];
      ARRAY2D(dfdu,1,5) = jacob[1][5];
      ARRAY2D(dfdu,1,6) = jacob[1][6];
      ARRAY2D(dfdu,1,7) = jacob[1][7];
      ARRAY2D(dfdu,1,8) = jacob[1][8];
      ARRAY2D(dfdu,1,9) = jacob[1][9];
      ARRAY2D(dfdu,1,10) = jacob[1][10];
      ARRAY2D(dfdu,1,11) = jacob[1][11];
      ARRAY2D(dfdu,1,12) = jacob[1][12];
      ARRAY2D(dfdu,1,13) = jacob[1][13];
      ARRAY2D(dfdu,1,14) = jacob[1][14];
      ARRAY2D(dfdu,2,0) = jacob[2][0];
      ARRAY2D(dfdu,2,1) = jacob[2][1];
      ARRAY2D(dfdu,2,2) = jacob[2][2];
      ARRAY2D(dfdu,2,3) = jacob[2][3];
      ARRAY2D(dfdu,2,4) = jacob[2][4];
      ARRAY2D(dfdu,2,5) = jacob[2][5];
      ARRAY2D(dfdu,2,6) = jacob[2][6];
      ARRAY2D(dfdu,2,7) = jacob[2][7];
      ARRAY2D(dfdu,2,8) = jacob[2][8];
      ARRAY2D(dfdu,2,9) = jacob[2][9];
      ARRAY2D(dfdu,2,10) = jacob[2][10];
      ARRAY2D(dfdu,2,11) = jacob[2][11];
      ARRAY2D(dfdu,2,12) = jacob[2][12];
      ARRAY2D(dfdu,2,13) = jacob[2][13];
      ARRAY2D(dfdu,2,14) = jacob[2][14];
      ARRAY2D(dfdu,3,0) = jacob[3][0];
      ARRAY2D(dfdu,3,1) = jacob[3][1];
      ARRAY2D(dfdu,3,2) = jacob[3][2];
      ARRAY2D(dfdu,3,3) = jacob[3][3];
      ARRAY2D(dfdu,3,4) = jacob[3][4];
      ARRAY2D(dfdu,3,5) = jacob[3][5];
      ARRAY2D(dfdu,3,6) = jacob[3][6];
      ARRAY2D(dfdu,3,7) = jacob[3][7];
      ARRAY2D(dfdu,3,8) = jacob[3][8];
      ARRAY2D(dfdu,3,9) = jacob[3][9];
      ARRAY2D(dfdu,3,10) = jacob[3][10];
      ARRAY2D(dfdu,3,11) = jacob[3][11];
      ARRAY2D(dfdu,3,12) = jacob[3][12];
      ARRAY2D(dfdu,3,13) = jacob[3][13];
      ARRAY2D(dfdu,3,14) = jacob[3][14];
      ARRAY2D(dfdu,4,0) = jacob[4][0];
      ARRAY2D(dfdu,4,1) = jacob[4][1];
      ARRAY2D(dfdu,4,2) = jacob[4][2];
      ARRAY2D(dfdu,4,3) = jacob[4][3];
      ARRAY2D(dfdu,4,4) = jacob[4][4];
      ARRAY2D(dfdu,4,5) = jacob[4][5];
      ARRAY2D(dfdu,4,6) = jacob[4][6];
      ARRAY2D(dfdu,4,7) = jacob[4][7];
      ARRAY2D(dfdu,4,8) = jacob[4][8];
      ARRAY2D(dfdu,4,9) = jacob[4][9];
      ARRAY2D(dfdu,4,10) = jacob[4][10];
      ARRAY2D(dfdu,4,11) = jacob[4][11];
      ARRAY2D(dfdu,4,12) = jacob[4][12];
      ARRAY2D(dfdu,4,13) = jacob[4][13];
      ARRAY2D(dfdu,4,14) = jacob[4][14];
      ARRAY2D(dfdu,5,0) = jacob[5][0];
      ARRAY2D(dfdu,5,1) = jacob[5][1];
      ARRAY2D(dfdu,5,2) = jacob[5][2];
      ARRAY2D(dfdu,5,3) = jacob[5][3];
      ARRAY2D(dfdu,5,4) = jacob[5][4];
      ARRAY2D(dfdu,5,5) = jacob[5][5];
      ARRAY2D(dfdu,5,6) = jacob[5][6];
      ARRAY2D(dfdu,5,7) = jacob[5][7];
      ARRAY2D(dfdu,5,8) = jacob[5][8];
      ARRAY2D(dfdu,5,9) = jacob[5][9];
      ARRAY2D(dfdu,5,10) = jacob[5][10];
      ARRAY2D(dfdu,5,11) = jacob[5][11];
      ARRAY2D(dfdu,5,12) = jacob[5][12];
      ARRAY2D(dfdu,5,13) = jacob[5][13];
      ARRAY2D(dfdu,5,14) = jacob[5][14];
      ARRAY2D(dfdu,6,0) = jacob[6][0];
      ARRAY2D(dfdu,6,1) = jacob[6][1];
      ARRAY2D(dfdu,6,2) = jacob[6][2];
      ARRAY2D(dfdu,6,3) = jacob[6][3];
      ARRAY2D(dfdu,6,4) = jacob[6][4];
      ARRAY2D(dfdu,6,5) = jacob[6][5];
      ARRAY2D(dfdu,6,6) = jacob[6][6];
      ARRAY2D(dfdu,6,7) = jacob[6][7];
      ARRAY2D(dfdu,6,8) = jacob[6][8];
      ARRAY2D(dfdu,6,9) = jacob[6][9];
      ARRAY2D(dfdu,6,10) = jacob[6][10];
      ARRAY2D(dfdu,6,11) = jacob[6][11];
      ARRAY2D(dfdu,6,12) = jacob[6][12];
      ARRAY2D(dfdu,6,13) = jacob[6][13];
      ARRAY2D(dfdu,6,14) = jacob[6][14];
      ARRAY2D(dfdu,7,0) = jacob[7][0];
      ARRAY2D(dfdu,7,1) = jacob[7][1];
      ARRAY2D(dfdu,7,2) = jacob[7][2];
      ARRAY2D(dfdu,7,3) = jacob[7][3];
      ARRAY2D(dfdu,7,4) = jacob[7][4];
      ARRAY2D(dfdu,7,5) = jacob[7][5];
      ARRAY2D(dfdu,7,6) = jacob[7][6];
      ARRAY2D(dfdu,7,7) = jacob[7][7];
      ARRAY2D(dfdu,7,8) = jacob[7][8];
      ARRAY2D(dfdu,7,9) = jacob[7][9];
      ARRAY2D(dfdu,7,10) = jacob[7][10];
      ARRAY2D(dfdu,7,11) = jacob[7][11];
      ARRAY2D(dfdu,7,12) = jacob[7][12];
      ARRAY2D(dfdu,7,13) = jacob[7][13];
      ARRAY2D(dfdu,7,14) = jacob[7][14];
      ARRAY2D(dfdu,8,0) = jacob[8][0];
      ARRAY2D(dfdu,8,1) = jacob[8][1];
      ARRAY2D(dfdu,8,2) = jacob[8][2];
      ARRAY2D(dfdu,8,3) = jacob[8][3];
      ARRAY2D(dfdu,8,4) = jacob[8][4];
      ARRAY2D(dfdu,8,5) = jacob[8][5];
      ARRAY2D(dfdu,8,6) = jacob[8][6];
      ARRAY2D(dfdu,8,7) = jacob[8][7];
      ARRAY2D(dfdu,8,8) = jacob[8][8];
      ARRAY2D(dfdu,8,9) = jacob[8][9];
      ARRAY2D(dfdu,8,10) = jacob[8][10];
      ARRAY2D(dfdu,8,11) = jacob[8][11];
      ARRAY2D(dfdu,8,12) = jacob[8][12];
      ARRAY2D(dfdu,8,13) = jacob[8][13];
      ARRAY2D(dfdu,8,14) = jacob[8][14];
      ARRAY2D(dfdu,9,0) = jacob[9][0];
      ARRAY2D(dfdu,9,1) = jacob[9][1];
      ARRAY2D(dfdu,9,2) = jacob[9][2];
      ARRAY2D(dfdu,9,3) = jacob[9][3];
      ARRAY2D(dfdu,9,4) = jacob[9][4];
      ARRAY2D(dfdu,9,5) = jacob[9][5];
      ARRAY2D(dfdu,9,6) = jacob[9][6];
      ARRAY2D(dfdu,9,7) = jacob[9][7];
      ARRAY2D(dfdu,9,8) = jacob[9][8];
      ARRAY2D(dfdu,9,9) = jacob[9][9];
      ARRAY2D(dfdu,9,10) = jacob[9][10];
      ARRAY2D(dfdu,9,11) = jacob[9][11];
      ARRAY2D(dfdu,9,12) = jacob[9][12];
      ARRAY2D(dfdu,9,13) = jacob[9][13];
      ARRAY2D(dfdu,9,14) = jacob[9][14];
      ARRAY2D(dfdu,10,0) = jacob[10][0];
      ARRAY2D(dfdu,10,1) = jacob[10][1];
      ARRAY2D(dfdu,10,2) = jacob[10][2];
      ARRAY2D(dfdu,10,3) = jacob[10][3];
      ARRAY2D(dfdu,10,4) = jacob[10][4];
      ARRAY2D(dfdu,10,5) = jacob[10][5];
      ARRAY2D(dfdu,10,6) = jacob[10][6];
      ARRAY2D(dfdu,10,7) = jacob[10][7];
      ARRAY2D(dfdu,10,8) = jacob[10][8];
      ARRAY2D(dfdu,10,9) = jacob[10][9];
      ARRAY2D(dfdu,10,10) = jacob[10][10];
      ARRAY2D(dfdu,10,11) = jacob[10][11];
      ARRAY2D(dfdu,10,12) = jacob[10][12];
      ARRAY2D(dfdu,10,13) = jacob[10][13];
      ARRAY2D(dfdu,10,14) = jacob[10][14];
      ARRAY2D(dfdu,11,0) = jacob[11][0];
      ARRAY2D(dfdu,11,1) = jacob[11][1];
      ARRAY2D(dfdu,11,2) = jacob[11][2];
      ARRAY2D(dfdu,11,3) = jacob[11][3];
      ARRAY2D(dfdu,11,4) = jacob[11][4];
      ARRAY2D(dfdu,11,5) = jacob[11][5];
      ARRAY2D(dfdu,11,6) = jacob[11][6];
      ARRAY2D(dfdu,11,7) = jacob[11][7];
      ARRAY2D(dfdu,11,8) = jacob[11][8];
      ARRAY2D(dfdu,11,9) = jacob[11][9];
      ARRAY2D(dfdu,11,10) = jacob[11][10];
      ARRAY2D(dfdu,11,11) = jacob[11][11];
      ARRAY2D(dfdu,11,12) = jacob[11][12];
      ARRAY2D(dfdu,11,13) = jacob[11][13];
      ARRAY2D(dfdu,11,14) = jacob[11][14];
      ARRAY2D(dfdu,12,0) = jacob[12][0];
      ARRAY2D(dfdu,12,1) = jacob[12][1];
      ARRAY2D(dfdu,12,2) = jacob[12][2];
      ARRAY2D(dfdu,12,3) = jacob[12][3];
      ARRAY2D(dfdu,12,4) = jacob[12][4];
      ARRAY2D(dfdu,12,5) = jacob[12][5];
      ARRAY2D(dfdu,12,6) = jacob[12][6];
      ARRAY2D(dfdu,12,7) = jacob[12][7];
      ARRAY2D(dfdu,12,8) = jacob[12][8];
      ARRAY2D(dfdu,12,9) = jacob[12][9];
      ARRAY2D(dfdu,12,10) = jacob[12][10];
      ARRAY2D(dfdu,12,11) = jacob[12][11];
      ARRAY2D(dfdu,12,12) = jacob[12][12];
      ARRAY2D(dfdu,12,13) = jacob[12][13];
      ARRAY2D(dfdu,12,14) = jacob[12][14];
      ARRAY2D(dfdu,13,0) = jacob[13][0];
      ARRAY2D(dfdu,13,1) = jacob[13][1];
      ARRAY2D(dfdu,13,2) = jacob[13][2];
      ARRAY2D(dfdu,13,3) = jacob[13][3];
      ARRAY2D(dfdu,13,4) = jacob[13][4];
      ARRAY2D(dfdu,13,5) = jacob[13][5];
      ARRAY2D(dfdu,13,6) = jacob[13][6];
      ARRAY2D(dfdu,13,7) = jacob[13][7];
      ARRAY2D(dfdu,13,8) = jacob[13][8];
      ARRAY2D(dfdu,13,9) = jacob[13][9];
      ARRAY2D(dfdu,13,10) = jacob[13][10];
      ARRAY2D(dfdu,13,11) = jacob[13][11];
      ARRAY2D(dfdu,13,12) = jacob[13][12];
      ARRAY2D(dfdu,13,13) = jacob[13][13];
      ARRAY2D(dfdu,13,14) = jacob[13][14];
      ARRAY2D(dfdu,14,0) = jacob[14][0];
      ARRAY2D(dfdu,14,1) = jacob[14][1];
      ARRAY2D(dfdu,14,2) = jacob[14][2];
      ARRAY2D(dfdu,14,3) = jacob[14][3];
      ARRAY2D(dfdu,14,4) = jacob[14][4];
      ARRAY2D(dfdu,14,5) = jacob[14][5];
      ARRAY2D(dfdu,14,6) = jacob[14][6];
      ARRAY2D(dfdu,14,7) = jacob[14][7];
      ARRAY2D(dfdu,14,8) = jacob[14][8];
      ARRAY2D(dfdu,14,9) = jacob[14][9];
      ARRAY2D(dfdu,14,10) = jacob[14][10];
      ARRAY2D(dfdu,14,11) = jacob[14][11];
      ARRAY2D(dfdu,14,12) = jacob[14][12];
      ARRAY2D(dfdu,14,13) = jacob[14][13];
      ARRAY2D(dfdu,14,14) = jacob[14][14];

      if (ijac == 1) {
          return 0;
      }

      /* Compute the derivative of the vector field  */
      /* with respect to the parameters.             */

      dfdpar[0][0] = 2.0*omega*(15.0*sin(phi[0])+14.0*sin(phi[1])+13.0*sin(phi
[2])+12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin
(phi[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*
sin(phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[1][0] = 2.0*omega*(14.0*sin(phi[0])+14.0*sin(phi[1])+13.0*sin(phi
[2])+12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin
(phi[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*
sin(phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[2][0] = 2.0*omega*(13.0*sin(phi[0])+13.0*sin(phi[1])+13.0*sin(phi
[2])+12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin
(phi[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*
sin(phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[3][0] = 2.0*omega*(12.0*sin(phi[0])+12.0*sin(phi[1])+12.0*sin(phi
[2])+12.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin
(phi[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*
sin(phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[4][0] = 2.0*omega*(11.0*sin(phi[0])+11.0*sin(phi[1])+11.0*sin(phi
[2])+11.0*sin(phi[3])+11.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin
(phi[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*
sin(phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[5][0] = 2.0*omega*(10.0*sin(phi[0])+10.0*sin(phi[1])+10.0*sin(phi
[2])+10.0*sin(phi[3])+10.0*sin(phi[4])+10.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin
(phi[7])+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*
sin(phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[6][0] = 2.0*omega*(9.0*sin(phi[0])+9.0*sin(phi[1])+9.0*sin(phi[2])
+9.0*sin(phi[3])+9.0*sin(phi[4])+9.0*sin(phi[5])+9.0*sin(phi[6])+8.0*sin(phi[7]
)+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi
[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[7][0] = 2.0*omega*(8.0*sin(phi[0])+8.0*sin(phi[1])+8.0*sin(phi[2])
+8.0*sin(phi[3])+8.0*sin(phi[4])+8.0*sin(phi[5])+8.0*sin(phi[6])+8.0*sin(phi[7]
)+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi
[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[8][0] = 2.0*omega*(7.0*sin(phi[0])+7.0*sin(phi[1])+7.0*sin(phi[2])
+7.0*sin(phi[3])+7.0*sin(phi[4])+7.0*sin(phi[5])+7.0*sin(phi[6])+7.0*sin(phi[7]
)+7.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi
[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[9][0] = 2.0*omega*(6.0*sin(phi[0])+6.0*sin(phi[1])+6.0*sin(phi[2])
+6.0*sin(phi[3])+6.0*sin(phi[4])+6.0*sin(phi[5])+6.0*sin(phi[6])+6.0*sin(phi[7]
)+6.0*sin(phi[8])+6.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(phi
[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[10][0] = 2.0*omega*(5.0*sin(phi[0])+5.0*sin(phi[1])+5.0*sin(phi[2]
)+5.0*sin(phi[3])+5.0*sin(phi[4])+5.0*sin(phi[5])+5.0*sin(phi[6])+5.0*sin(phi
[7])+5.0*sin(phi[8])+5.0*sin(phi[9])+5.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[11][0] = 2.0*omega*(4.0*sin(phi[0])+4.0*sin(phi[1])+4.0*sin(phi[2]
)+4.0*sin(phi[3])+4.0*sin(phi[4])+4.0*sin(phi[5])+4.0*sin(phi[6])+4.0*sin(phi
[7])+4.0*sin(phi[8])+4.0*sin(phi[9])+4.0*sin(phi[10])+4.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[12][0] = 2.0*omega*(3.0*sin(phi[0])+3.0*sin(phi[1])+3.0*sin(phi[2]
)+3.0*sin(phi[3])+3.0*sin(phi[4])+3.0*sin(phi[5])+3.0*sin(phi[6])+3.0*sin(phi
[7])+3.0*sin(phi[8])+3.0*sin(phi[9])+3.0*sin(phi[10])+3.0*sin(phi[11])+3.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[13][0] = 2.0*omega*(2.0*sin(phi[0])+2.0*sin(phi[1])+2.0*sin(phi[2]
)+2.0*sin(phi[3])+2.0*sin(phi[4])+2.0*sin(phi[5])+2.0*sin(phi[6])+2.0*sin(phi
[7])+2.0*sin(phi[8])+2.0*sin(phi[9])+2.0*sin(phi[10])+2.0*sin(phi[11])+2.0*sin(
phi[12])+2.0*sin(phi[13])+sin(phi[14]));
      dfdpar[14][0] = 2.0*omega*(sin(phi[0])+sin(phi[1])+sin(phi[2])+sin(phi[3]
)+sin(phi[4])+sin(phi[5])+sin(phi[6])+sin(phi[7])+sin(phi[8])+sin(phi[9])+sin(
phi[10])+sin(phi[11])+sin(phi[12])+sin(phi[13])+sin(phi[14]));
      ARRAY2D(dfdp,0,0) = dfdpar[0][0];
      ARRAY2D(dfdp,1,0) = dfdpar[1][0];
      ARRAY2D(dfdp,2,0) = dfdpar[2][0];
      ARRAY2D(dfdp,3,0) = dfdpar[3][0];
      ARRAY2D(dfdp,4,0) = dfdpar[4][0];
      ARRAY2D(dfdp,5,0) = dfdpar[5][0];
      ARRAY2D(dfdp,6,0) = dfdpar[6][0];
      ARRAY2D(dfdp,7,0) = dfdpar[7][0];
      ARRAY2D(dfdp,8,0) = dfdpar[8][0];
      ARRAY2D(dfdp,9,0) = dfdpar[9][0];
      ARRAY2D(dfdp,10,0) = dfdpar[10][0];
      ARRAY2D(dfdp,11,0) = dfdpar[11][0];
      ARRAY2D(dfdp,12,0) = dfdpar[12][0];
      ARRAY2D(dfdp,13,0) = dfdpar[13][0];
      ARRAY2D(dfdp,14,0) = dfdpar[14][0];

      return 0;
}

/*--------------------------------------------------------------------------*/
/* STPNT  Gives a starting pointing                                         */
/*--------------------------------------------------------------------------*/

int stpnt(integer ndim, doublereal t,
          doublereal *u, doublereal *par)
{
      doublereal phi[15];
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
      phi[7] = 0.0;
      phi[8] = 0.0;
      phi[9] = 0.0;
      phi[10] = 0.0;
      phi[11] = 0.0;
      phi[12] = 0.0;
      phi[13] = 0.0;
      phi[14] = 0.0;

      par[0] = omega;
      u[0] = phi[0];
      u[1] = phi[1];
      u[2] = phi[2];
      u[3] = phi[3];
      u[4] = phi[4];
      u[5] = phi[5];
      u[6] = phi[6];
      u[7] = phi[7];
      u[8] = phi[8];
      u[9] = phi[9];
      u[10] = phi[10];
      u[11] = phi[11];
      u[12] = phi[12];
      u[13] = phi[13];
      u[14] = phi[14];
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
