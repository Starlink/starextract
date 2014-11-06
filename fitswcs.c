/*
*				fitswcs.c
*
* Manage World Coordinate System data.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		13/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*
 *  History:
 *         05/09/2014 (PWD): Convert to AST as necessary.
 */

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include        "ast.h"
#include        "define.h"
#include        "globals.h"
#include	"fitswcs.h"

/******* copy_wcs ************************************************************
PROTO	wcsstruct *copy_wcs(wcsstruct *wcsin)
PURPOSE	Copy a WCS (World Coordinate System) structure.
INPUT	WCS structure to be copied.
OUTPUT	pointer to a copy of the input structure.
AUTHOR	E. Bertin (IAP)
VERSION	31/08/2002
 ***/
wcsstruct	*copy_wcs(wcsstruct *wcsin)
  {
  wcsstruct	*wcs = NULL;
  QMEMCPY(wcsin, wcs, wcsstruct, 1);
  return wcs;
  }


/******* end_wcs **************************************************************
PROTO	void end_wcs(wcsstruct *wcs)
PURPOSE	Free WCS (World Coordinate System) infos.
INPUT	WCS structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	24/05/2000
 ***/
void	end_wcs(wcsstruct *wcs)
  {
  if ( wcs )
      free( wcs );
  return;
  }


/******* raw_to_wcs ***********************************************************
PROTO	int raw_to_wcs(wcsstruct *, double *, double *)
PURPOSE	Convert raw (pixel) coordinates to WCS (World Coordinate System).
INPUT	WCS structure,
	Pointer to the array of input coordinates,
	Pointer to the array of output coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	raw_to_wcs(wcsstruct *wcs, double *pixpos, double *wcspos)
  {
  double xin[1], yin[1];
  double xout[1], yout[1];

  xin[0] = pixpos[0];
  yin[0] = pixpos[1];
  if ( wcs->astwcs != NULL ) {
      astTran2( wcs->astwcs, 1, xin, yin, 1, xout, yout );
      wcspos[0] = xout[0];
      wcspos[1] = yout[0];
      astNorm( wcs->astwcs, wcspos );
      wcspos[0] /= DEG;
      wcspos[1] /= DEG;
  }
  else {
      wcspos[0] = 0.0;
      wcspos[1] = 0.0;
  }
  if ( ! astOK )
    {
    astClearStatus;
    return RETURN_ERROR;
    }
  return RETURN_OK;
  }


/******* wcs_to_raw ***********************************************************
PROTO	int wcs_to_raw(wcsstruct *, double *, double *)
PURPOSE	Convert WCS (World Coordinate System) coords to raw (pixel) coords.
INPUT	WCS structure,
	Pointer to the array of input coordinates,
	Pointer to the array of output coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	wcs_to_raw(wcsstruct *wcs, double *wcspos, double *pixpos)
  {
  double xin[1], yin[1];
  double xout[1], yout[1];

  xin[0] = wcspos[0] * DEG;
  yin[0] = wcspos[1] * DEG;
  if ( wcs->astwcs != NULL ) {
      astTran2( wcs->astwcs, 1, xin, yin, 0, xout, yout );
      pixpos[0] = xout[0];
      pixpos[1] = yout[0];
  }
  else {
      pixpos[0] = 0.0;
      pixpos[1] = 0.0;
  }
  if ( ! astOK )
    {
    astClearStatus;
    return RETURN_ERROR;
    }
  return RETURN_OK;
  }

/******* red_to_raw **********************************************************
PROTO   int red_to_raw(wcsstruct *, double *, double *)
PURPOSE Convert reduced (World Coordinate System) coords to raw (pixel)
        coords.
INPUT   WCS structure,
        Pointer to the array of input (reduced) coordinates,
        Pointer to the array of output (pixel) coordinates.
OUTPUT  RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 23/10/2003
 ***/
int     red_to_raw(wcsstruct *wcs, double *redpos, double *pixpos)
  {
  /*  This concept doesn't make sense in NDF/AST (transforms using the
   *  CD/PC terms only), so just return dummy values and an error. */
  pixpos[0] = 0.0;
  pixpos[1] = 0.0;
  return RETURN_ERROR;
  }

/******* raw_to_red **********************************************************
PROTO   int raw_to_red(wcsstruct *, double *, double *)
PURPOSE Convert raw (pixel) coordinates to reduced WCS coordinates.
INPUT   WCS structure,
        Pointer to the array of input (pixel) coordinates,
        Pointer to the array of output (reduced) coordinates.
OUTPUT  RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 23/10/2003
 ***/
int     raw_to_red(wcsstruct *wcs, double *pixpos, double *redpos)
  {
  redpos[0] = 0.0;
  redpos[1] = 0.0;
  return RETURN_ERROR;
  }

/****** wcs jacobian *********************************************************
PROTO	double wcs_jacobian(wcsstruct *wcs, double *pixpos, double *jacob)
PURPOSE	Compute the local Jacobian matrix of the astrometric deprojection.
INPUT	WCS structure,
	Pointer to the array of local raw coordinates,
	Pointer to the jacobian array (output).
OUTPUT	Determinant over spatial coordinates (=pixel area), or -1.0 if mapping
	was unsuccesful.
NOTES   Memory must have been allocated (naxis*naxis*sizeof(double)) for the
        Jacobian array.
AUTHOR	E. Bertin (IAP)
VERSION	11/10/2007
 ***/
double	wcs_jacobian(wcsstruct *wcs, double *pixpos, double *jacob)
  {
   double	pixpos0[NAXIS], wcspos0[NAXIS], wcspos[NAXIS],
		dpos;
   int		i,j, lng,lat,naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    pixpos0[i] = pixpos[i];
  if (raw_to_wcs(wcs, pixpos0, wcspos0) == RETURN_ERROR)
    return -1.0;
  for (i=0; i<naxis; i++)
    {
    pixpos0[i] += 1.0;
    if (raw_to_wcs(wcs, pixpos0, wcspos) == RETURN_ERROR)
      return -1.0;
    pixpos0[i] -= 1.0;
    for (j=0; j<naxis; j++)
      {
      dpos = wcspos[j]-wcspos0[j];
      if (lng!=lat && j==lng)
        {
        if (dpos>180.0)
          dpos -= 360.0;
        else if (dpos<-180.0)
          dpos += 360.0;
        dpos *= cos(wcspos0[lat]*DEG);
        }
      jacob[j*naxis+i] = dpos;
      }
    }

  if (lng==lat)
    {
    lng = 0;
    lat = 1;
    }

  return fabs(jacob[lng+naxis*lng]*jacob[lat+naxis*lat]
		- jacob[lat+naxis*lng]*jacob[lng+naxis*lat]);
  }


/********************************* precess ***********************************/
/*
precess equatorial coordinates according to the equinox (from Ephemerides du
Bureau des Longitudes 1992). Epoch for coordinates should be J2000
(FK5 system).
*/
void	precess(double yearin, double alphain, double deltain,
		double yearout, double *alphaout, double *deltaout)

  {
   double	dzeta,theta,z, t1,t1t1, t2,t2t2,t2t2t2,
		cddsadz, cddcadz, cdd, sdd, adz, cdin,sdin,ct,st,caindz;

  alphain *= DEG;
  deltain *= DEG;

  t1 = (yearin - 2000.0)/1000.0;
  t2 = (yearout - yearin)/1000.0;
  t1t1 = t1*t1;
  t2t2t2 = (t2t2 = t2*t2)*t2;
  theta = (97171.735e-06 - 413.691e-06*t1 - 1.052e-06 * t1t1) * t2
	+ (-206.846e-06 - 1.052e-06*t1) * t2t2 - 202.812e-06 * t2t2t2;
  dzeta = (111808.609e-06 + 677.071e-06*t1 - 0.674e-06 * t1t1) * t2
	+ (146.356e-06 - 1.673e-06*t1) * t2t2 + 87.257e-06 * t2t2t2;
  z = (111808.609e-06 +677.071e-06*t1 - 0.674e-06 * t1t1) * t2
	+ (530.716e-06 + 0.320e-06*t1) * t2t2 + 88.251e-06 * t2t2t2;
  cddsadz = (cdin=cos(deltain)) * sin(alphain+dzeta);
  cddcadz = -(sdin=sin(deltain))*(st=sin(theta))
	+cdin*(ct=cos(theta))*(caindz=cos(alphain+dzeta));
  sdd = sdin*ct + cdin*st*caindz;
  cdd = cos(*deltaout = asin(sdd));
  adz = asin(cddsadz/cdd);
  if (cddcadz<0.0)
    adz = PI - adz;
  if (adz<0.0)
    adz += 2.0*PI;
  adz += z;
  *alphaout = adz/DEG;
  *deltaout /= DEG;

  return;
  }


/*********************************** j2b *************************************/
/*
conver equatorial coordinates from equinox and epoch J2000 to equinox and
epoch B1950 for extragalactic sources (from Aoki et al. 1983, after
inversion of their matrix and some custom arrangements).
*/
void    j2b(double yearobs, double alphain, double deltain,
	double *alphaout, double *deltaout)
  {
   int			i,j;
   double		a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6},
                        ap[3] = {1.245e-3, -1.580e-3, -0.659e-3},
                        m[6][6] = {
  { 0.9999256794678425,    0.01118148281196562,   0.004859003848996022,
   -2.423898417033081e-06,-2.710547600126671e-08,-1.177738063266745e-08},
  {-0.01118148272969232,   0.9999374849247641,   -2.717708936468247e-05,
    2.710547578707874e-08,-2.423927042585208e-06, 6.588254898401055e-11},
  {-0.00485900399622881,  -2.715579322970546e-05, 0.999988194643078,
    1.177738102358923e-08, 6.582788892816657e-11,-2.424049920613325e-06},
  {-0.0005508458576414713, 0.2384844384742432,   -0.4356144527773499,
    0.9999043171308133,    0.01118145410120206,   0.004858518651645554},
  {-0.2385354433560954,   -0.002664266996872802,  0.01225282765749546,
   -0.01118145417187502,   0.9999161290795875,   -2.717034576263522e-05},
  { 0.4357269351676567,   -0.008536768476441086,  0.002113420799663768,
   -0.004858518477064975, -2.715994547222661e-05, 0.9999668385070383}},
			a1[3], r[3], ro[3], r1[3], r2[3], v1[3], v[3];
   double		cai, sai, cdi, sdi, dotp, rmod, alpha, delta, t1;

/* Convert Julian years from J2000.0 to tropic centuries from B1950.0 */
  t1 = ((yearobs - 2000.0) + (MJD2000 - MJD1950)/365.25)*JU2TROP/100.0;
  alphain *= DEG;
  deltain *= DEG;
  cai = cos(alphain);
  sai = sin(alphain);
  cdi = cos(deltain);
  sdi = sin(deltain);
  r[0] = cdi*cai;
  r[1] = cdi*sai;
  r[2] = sdi;
  for (i=0; i<3; i++)
    v[i] = r2[i] = v1[i] = 0.0;
  for (j=0; j<6; j++)
    for (i=0; i<6; i++)
      if (j<3)
        r2[j] += m[j][i]*(i<3?r[i]:v[i-3]);
      else
        v1[j-3] += m[j][i]*(i<3?r[i]:v[i-3]);

  for (i=0; i<3; i++)
    r1[i] = r2[i]+v1[i]*ARCSEC*t1;

  dotp = 0.0;
  for (i=0; i<3; i++)
    {
    a1[i] = a[i]+ap[i]*ARCSEC*t1;
    dotp += a1[i]*(r1[i]+a1[i]);
    }
  dotp = 2.0/(sqrt(1+4.0*dotp)+1.0);
  rmod = 0.0;
  for (i=0; i<3; i++)
    {
    ro[i] = dotp*(r1[i]+a1[i]);
    rmod += ro[i]*ro[i];
    }
  rmod = sqrt(rmod);
  delta = asin(ro[2]/rmod);
  alpha = acos(ro[0]/cos(delta)/rmod);
  if (ro[1]<0)
    alpha = 2.0*PI - alpha;
  *alphaout = alpha/DEG;
  *deltaout = delta/DEG;

  return;
  }


/******************************** fmod_0_p360 *******************************/
/*
Fold input angle in the [0,+360[ domain.
*/
double  fmod_0_p360(double angle)
  {
  return angle>0.0? fmod(angle,360.0) : fmod(angle,360.0)+360.0;
  }


/******************************** fmod_m90_p90 *******************************/
/*
Fold input angle in the [-90,+90[ domain.
*/
double  fmod_m90_p90(double angle)
  {
  return angle>0.0? fmod(angle+90.0,180.0)-90.0 : fmod(angle-90.0,180.0)+90.0;
  }


/********************************* fcmp_0_p360 *******************************/
/*
Compare angles in the [0,+360[ domain: return 1 if anglep>anglem, 0 otherwise.
*/
int  fcmp_0_p360(double anglep, double anglem)
  {
   double dval = anglep - anglem;

  return (int)((dval>0.0 && dval<180.0) || dval<-180.0);
  }

