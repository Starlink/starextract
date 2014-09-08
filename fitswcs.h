/*
*				fitswcs.h
*
* Include file for fitswcs.c
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
*	Last modified:		14/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#define _FITSWCS_H_

/*-------------------------------- macros -----------------------------------*/

/*----------------------------- Internal constants --------------------------*/

#define		NAXIS	2		/* Max number of FITS axes */

#define		DEG	(PI/180.0)	/* 1 deg in radians */
#define		ARCMIN	(DEG/60.0)	/* 1 arcsec in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MAS	(ARCSEC/1000.0)	/* 1 mas in radians */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/

/*------------------------------- structures --------------------------------*/

typedef struct wcs
  {
  int		naxis;			/* Number of image axes */
  int		lat,lng;		/* longitude and latitude axes # */
  double	pixscale;		/* (Local) pixel scale */
  double	ap2000,dp2000;		/* J2000 coordinates of pole */
  double	ap1950,dp1950;		/* B1950 coordinates of pole */
  double	equinox;		/* Equinox of observations */
  double	epoch;			/* Epoch of observations (deprec.) */
  AstFrameSet   *astwcs;                /* AST FrameSet, the primary WCS */
  }	wcsstruct;

/*------------------------------- functions ---------------------------------*/

extern wcsstruct	*copy_wcs(wcsstruct *wcsin),
			*read_wcs(AstFrameSet *frmset);

extern double		fmod_0_p360(double angle),
			fmod_m90_p90(double angle),
			wcs_jacobian(wcsstruct *wcs, double *pixpos,
                                     double *jacob);

extern int		
			fcmp_0_p360(double anglep, double anglem),
			raw_to_wcs(wcsstruct *wcs,
                                   double *pixpos, double *wcspos),
			wcs_to_raw(wcsstruct *wcs,
                                   double *wcspos, double *pixpos),
                        raw_to_red(wcsstruct *wcs,
                                   double *pixpos, double *redpos),
                        red_to_raw(wcsstruct *wcs,
                                   double *redpos, double *pixpos);


extern void		
			end_wcs(wcsstruct *wcs),
			j2b(double yearobs, double alphain, double deltain,
                            double *alphaout, double *deltaout),
			precess(double yearin, double alphain, double deltain,
				double yearout,
				double *alphaout, double *deltaout);

#endif
