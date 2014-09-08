/*
*				field.c
*
* Handle field (image) structures.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SExtractor is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SExtractor is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		12/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*
*       History:
*                       27/10/98 (AJC)
*                          Use AFPRINTF not fprintf
*                       06/01/99 (PWD)
*                          Changed use of field->file member. This is
*                          used differently in NDF interface (was
*                          being closed!). 
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"assoc.h"
#include	"astrom.h"
#include	"back.h"
#include	"field.h"
#include	"filter.h"
#include	"fitswcs.h"
#include	"header.h"
#include	"interpolate.h"
#include        "adam_defs.h"

/********************************* newfield **********************************/
/*
Returns a pointer to a new field for extension ext, ready to go!
*/
picstruct	*newfield(char *filename, int flags, int ext)

  {
   picstruct	*field;
   int margin;

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(field, picstruct, 1);
  field->flags = flags;
  strcpy (field->filename, filename);
/* A short, "relative" version of the filename */
  if (!(field->rfilename = strrchr(field->filename, '/')))
    field->rfilename = field->filename;
  else
    field->rfilename++;

  sprintf(gstr, "Looking for %s", field->rfilename);
  NFPRINTF(OUTPUT, gstr);
/* Check the image exists and read important info (image size, etc...) */
  readimagehead(field);

  QPRINTF(OUTPUT, "----- %s %s%s\n",
        flags&FLAG_FIELD?   "Flagging  from:" :
       (flags&(RMS_FIELD|VAR_FIELD|WEIGHT_FIELD)?
                             "Weighting from:" :
       (flags&MEASURE_FIELD? "Measuring from:" :
                             "Detecting from:")),
        field->rfilename,
        gstr );
  QPRINTF(OUTPUT, "      \"%.20s\" / %s / %dx%d / %d bits\n",
        field->ident,
        field->headflag? "EXT. HEADER" : "no ext. header",
          field->width, field->height, field->bytepix*8 );

/* Check the astrometric system and do the setup of the astrometric stuff */
  if (prefs.world_flag && (flags & (MEASURE_FIELD|DETECT_FIELD)))
    initastrom(field);
  else
    field->pixscale=prefs.pixel_scale;

/* Gain and Saturation */
/* XXX Need NDF FITS header access.
  if (flags & (DETECT_FIELD|MEASURE_FIELD))
    {
    if (fitsread(field->tab->headbuf, prefs.gain_key, &field->gain,
	H_FLOAT, T_DOUBLE) != RETURN_OK)
      field->gain = prefs.gain;
    if (fitsread(field->tab->headbuf, prefs.satur_key, &field->satur_level,
	H_FLOAT, T_DOUBLE) !=RETURN_OK)
      field->satur_level = prefs.satur_level;
    }
*/

/* Background */
  if (flags & (DETECT_FIELD|MEASURE_FIELD|WEIGHT_FIELD|VAR_FIELD|RMS_FIELD))
    {
    field->ngamma = prefs.mag_gamma/log(10.0);

    field->backw = prefs.backsize[0]<field->width ? prefs.backsize[0]
						  : field->width;
    field->backh = prefs.backsize[1]<field->height ? prefs.backsize[1]
						   : field->height;
    field->nbackp = field->backw * field->backh;
    if ((field->nbackx = (field->width-1)/field->backw + 1) < 1)
      field->nbackx = 1;
    if ((field->nbacky = (field->height-1)/field->backh + 1) < 1)
      field->nbacky = 1;
    field->nback = field->nbackx * field->nbacky;
    field->nbackfx = field->nbackx>1 ? prefs.backfsize[0] : 1;
    field->nbackfy = field->nbacky>1 ? prefs.backfsize[1] : 1;
/*--  Set the back_type flag if absolute background is selected */
    if (((flags & DETECT_FIELD) && prefs.back_type[0]==BACK_ABSOLUTE)
	|| ((flags & MEASURE_FIELD) && prefs.back_type[1]==BACK_ABSOLUTE))
      field->back_type = BACK_ABSOLUTE;
    }

/* Add a comfortable margin for local background estimates */
  margin = (prefs.pback_type == LOCAL)? prefs.pback_size + prefs.mem_bufsize/4
					: 0;

  field->stripheight = prefs.mem_bufsize + margin;
  if (field->stripheight>field->height)
    field->stripheight = field->height;
/* Compute the image buffer size */
/* Basically, only one margin line is sufficient... */
  field->stripmargin = 1 + margin;
/* ...but : */
  if (prefs.filter_flag)
    {
/*-- If filtering is on, one should consider the height of the conv. mask */
/*-- + 1 line for detectinhg zero-weight neighbours */
    if (field->stripheight < thefilter->convh)
      field->stripheight = thefilter->convh;
    if (field->stripmargin < (margin = (thefilter->convh-1)/2+1))
      field->stripmargin = margin;
    }

  return field;
  }


/******************************* inheritfield *******************************/
/*
Make a copy of a field structure, e.g. for interpolation purposes.
*/
picstruct	*inheritfield(picstruct *infield, int flags)

  {
   picstruct	*field;

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(field, picstruct, 1);

/* Copy what is important and reset the remaining */
  *field = *infield;
  field->flags = flags;
  if (infield->wcs)
    field->wcs = copy_wcs(infield->wcs);
  field->interp_flag = 0;
  field->assoc = NULL;
  field->strip = NULL;
  field->fstrip = NULL;
  field->reffield = infield;
  field->file = 0;            /*  PWD: was NULL, file reused as position;*/

  return field;
  }


/********************************* endfield **********************************/
/*
Free and close everything related to a field structure.
*/
void	endfield(picstruct *field)

  {

  free(field->strip);
  free(field->fstrip);
  if (field->wcs)
    end_wcs(field->wcs);
  if (field->interp_flag)
    end_interpolate(field);
  endback(field);
  free(field);

  return;
  }

