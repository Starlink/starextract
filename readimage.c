/*
*				readimage.c.
*
* Read image data.
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
*	Last modified:		26/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*
*       History:
*	        	28/10/98 (AJC)
*                          In line with V2.0.15
*                       14/12/98 (PWD):
*                          Added USHORT and UBYTE support.
*                       17/12/98 (PWD):
*                          Changed to use NDF WCS component for astrometry.
*	                27/11/2003 (EB):
*                       23-NOV-2005 (TIMJ):
*                          Remove DAT__ROOT
*                       26/01/2006 (PWD):
*                          Changed to handle redundant axes in data and native
*                          data type of NDF.
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include        "ast.h"

#include	"define.h"
#include	"globals.h"
#include        "prefs.h"
#include	"check.h"
#include	"field.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include	"interpolate.h"
#include	"back.h"
#include	"astrom.h"
#include	"weight.h"

#include        "sae_par.h"
#include        "star/hds.h"
#include        "ndf.h"
#include        "merswrap.h"

/******************************* loadstrip ***********************************/
/*
Load a new strip of pixel data into the buffer.
*/
void	*loadstrip(picstruct *field, picstruct *wfield)
  {
   tabstruct	*tab;
   checkstruct	*check;
   int		y, w, flags, interpflag;
   PIXTYPE	*data, *wdata, *rmsdata;

  w = field->width;
  flags = field->flags;
  interpflag = (wfield && wfield->interp_flag);
  wdata = NULL;			/* To avoid gcc -Wall warnings */

  if (!field->y)
    {
/*- First strip */
     int	nbpix;

    nbpix = w*field->stripheight;

    if (flags ^ FLAG_FIELD)
      {
/*---- Allocate space for the frame-buffer */
      if (!(field->strip=(PIXTYPE *)malloc(field->stripheight*field->width
        *sizeof(PIXTYPE))))
        error(EXIT_FAILURE,"Not enough memory for the image buffer of ",
		field->rfilename);

      data = field->strip;
/*---- We assume weight data have been read just before */
      if (interpflag)
        wdata = wfield->strip;
      if (flags & BACKRMS_FIELD)
        for (y=0, rmsdata=data; y<field->stripheight; y++, rmsdata += w)
          backrmsline(field, y, rmsdata);
      else if (flags & INTERP_FIELD)
        copydata(field, 0, nbpix);
      else
        readdata(field, data, nbpix);
      if (flags & (WEIGHT_FIELD|RMS_FIELD|BACKRMS_FIELD|VAR_FIELD))
        weight_to_var(field, data, nbpix);
      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_IDENTICAL]))
        writecheck(check, data, nbpix);
      for (y=0; y<field->stripheight; y++, data += w)
        {
/*------ This is the only place where one can pick-up safely the current bkg */
        if (flags & (MEASURE_FIELD|DETECT_FIELD))
          subbackline(field, y, data);
/*------ Go to interpolation process */
        if (interpflag)
          {
          interpolate(field,wfield, data, wdata);
          wdata += w;
          }
/*------ Check-image stuff */
        if (prefs.check_flag)
          {
          if (flags & MEASURE_FIELD)
            {
            if ((check = prefs.check[CHECK_BACKGROUND]))
              writecheck(check, field->backline, w);
            if ((check = prefs.check[CHECK_SUBTRACTED]))
              writecheck(check, data, w);
            if ((check = prefs.check[CHECK_APERTURES]))
              writecheck(check, data, w);
            if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
              writecheck(check, data, w);
            if ((check = prefs.check[CHECK_SUBPCPROTOS]))
              writecheck(check, data, w);
            if ((check = prefs.check[CHECK_SUBPROFILES]))
              writecheck(check, data, w);
            if ((check = prefs.check[CHECK_SUBSPHEROIDS]))
              writecheck(check, data, w);
            if ((check = prefs.check[CHECK_SUBDISKS]))
              writecheck(check, data, w);
            }
          if ((flags&DETECT_FIELD) && (check=prefs.check[CHECK_BACKRMS]))
            {
            backrmsline(field, y, (PIXTYPE *)check->pix);
            writecheck(check, check->pix, w);
            }
          }
        }
      }
    else
      {
      if (!(field->fstrip=(FLAGTYPE *)malloc(field->stripheight*field->width
		*sizeof(FLAGTYPE))))
      error(EXIT_FAILURE,"Not enough memory for the flag buffer of ",
	field->rfilename);
      readidata(field, field->fstrip, nbpix);
      }

    field->ymax = field->stripheight;
    if (field->ymax < field->height)
      field->stripysclim = field->stripheight - field->stripmargin;
    }
  else
    {
/*- other strips */
    if (flags ^ FLAG_FIELD)
      {
      data = field->strip + field->stripylim*w;
/*---- We assume weight data have been read just before */
      if (interpflag)
        wdata = wfield->strip + field->stripylim*w;

/*---- copy to Check-image the "oldest" line before it is replaced */
      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_SUBOBJECTS]))
        writecheck(check, data, w);

      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_MASK]))
        writecheck(check, data, w);

      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_SUBMASK]))
        writecheck(check, data, w);

      if (flags & BACKRMS_FIELD)
        backrmsline(field, field->ymax, data);
      else if (flags & INTERP_FIELD)
        copydata(field, field->stripylim*w, w);
      else
        readdata(field, data, w);

      if (flags & (WEIGHT_FIELD|RMS_FIELD|BACKRMS_FIELD|VAR_FIELD))
        weight_to_var(field, data, w);

      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_IDENTICAL]))
        writecheck(check, data, w);
/*---- Interpolate and subtract the background at current line */
      if (flags & (MEASURE_FIELD|DETECT_FIELD))
        subbackline(field, field->ymax, data);
      if (interpflag)
        interpolate(field,wfield, data, wdata);
/*---- Check-image stuff */
      if (prefs.check_flag)
        {
        if (flags & MEASURE_FIELD)
          {
          if ((check = prefs.check[CHECK_BACKGROUND]))
            writecheck(check, field->backline, w);
          if ((check = prefs.check[CHECK_SUBTRACTED]))
            writecheck(check, data, w);
          if ((check = prefs.check[CHECK_APERTURES]))
            writecheck(check, data, w);
          if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
            writecheck(check, data, w);
          if ((check = prefs.check[CHECK_SUBPCPROTOS]))
            writecheck(check, data, w);
          if ((check = prefs.check[CHECK_SUBPROFILES]))
            writecheck(check, data, w);
          if ((check = prefs.check[CHECK_SUBSPHEROIDS]))
            writecheck(check, data, w);
          if ((check = prefs.check[CHECK_SUBDISKS]))
            writecheck(check, data, w);
          }
        if ((flags&DETECT_FIELD) && (check=prefs.check[CHECK_BACKRMS]))
          {
          backrmsline(field, field->ymax, (PIXTYPE *)check->pix);
          writecheck(check, check->pix, w);
          }
        }
      }
    else
      readidata(field, field->fstrip + field->stripylim*w, w);

    field->stripylim = (++field->ymin)%field->stripheight;
    if ((++field->ymax)<field->height)
      field->stripysclim = (++field->stripysclim)%field->stripheight;
    }

  return (flags ^ FLAG_FIELD)?
		  (void *)(field->strip + field->stripy*w)
		: (void *)(field->fstrip + field->stripy*w);
  }


/******************************** copydata **********************************/
/*
Copy image data from one field to the other.
*/
void	copydata(picstruct *field, int offset, int size)
  {
  memcpy(field->strip+offset, field->reffield->strip+offset,
		size*sizeof(PIXTYPE));
  return;
  }


/******************************** readdata **********************************/
/*
read and convert input data stream in PIXTYPE (float) format.
  field    is pointer to the image picstruct
    field->map     pointer to mapped data
    field->file    number of next pixel to be read (from 0)
    field->nel     number of pixels in image
  ptr      is pointer to the strip buffer
  size     is the number of pixels to read
*/
void	readdata(picstruct *field, PIXTYPE *ptr, int size)
  {
  int		i,left;
  PIXTYPE	bs,bz;

  bs = (PIXTYPE)field->bscale;
  bz = (PIXTYPE)field->bzero;

  left = field->nel - field->file;

  switch(field->bitpix)
    {
    case BP_BYTE:
        if ( field->bitsgn ) {
            for ( i=0; i< (size < left ? size : left ); i++ ) {
                *(ptr++) = (PIXTYPE) *(((char *)field->map) + field->file++)*bs+bz;
            }
        }
        else {
            for ( i=0; i< (size < left ? size : left ); i++ ) {
                *(ptr++) = (PIXTYPE) *(((BYTE *)field->map) + field->file++)*bs+bz;
            }
        }
        break;

    case BP_SHORT:
        if ( field->bitsgn ) {
            for ( i=0; i< (size < left ? size : left ); i++ ) {
                *(ptr++) = (PIXTYPE) *(((short *)field->map) + field->file++)*bs+bz;
            }
        }
        else {
            for ( i=0; i< (size < left ? size : left ); i++ ) {
                *(ptr++) = (PIXTYPE) *(((USHORT *)field->map) + field->file++)*bs+bz;
            }
        }
        break;

    case BP_LONG:
        if ( field->bitsgn ) {
            for ( i=0; i< (size < left ? size : left ); i++ ) {
                *(ptr++) = (PIXTYPE) *(((LONG *)field->map) + field->file++)*bs+bz;
            }
        }
        else {
            for ( i=0; i< (size < left ? size : left ); i++ ) {
                *(ptr++) = (PIXTYPE) *(((ULONG *)field->map) + field->file++)*bs+bz;
            }
        }
        break;

    case BP_FLOAT:
        for ( i=0; i< (size < left ? size : left ); i++ ) {
            *(ptr++) = (PIXTYPE) *(((float *)field->map) + field->file++)*bs+bz;
        }
        break;

    case BP_DOUBLE:
        for ( i=0; i< (size < left ? size : left ); i++ ) {
            *(ptr++) = (PIXTYPE) *(((double *)field->map) + field->file++)*bs+bz;
        }
        break;

    default:
        error(EXIT_FAILURE,"*FATAL ERROR*: unknown BITPIX type in ",
              "readdata()");
        break;
    }


/* Reset field->file if at end */
/* This prevents endfield from crashing if field->file is used for NDF */
  field->file = (field->file >= field->nel)?0:field->file;
  return;
  }


/******************************** readidata *********************************/
/*
read and convert input data stream in FLAGTYPE (unsigned int) format.
  field    is pointer to the image picstruct
    field->file   number of next pixel to be read (from 0)
    field->nel    number of pixels in image
  ptr      is pointer to the strip buffer
  size     is the number of pixels to read
*/
void	readidata(picstruct *field, FLAGTYPE *ptr, int size)
  {
  int		i, left;

  left = field->nel - field->file;
  for (i=0; i<(size<left?size:left); i++)
     *(ptr++) = *((FLAGTYPE *)field->map + field->file++);

/* Reset field->file if at end */
/* This prevents endfield from crashing if field->file is used for NDF */
  field->file = (field->file >= field->nel)?0:field->file;
  return;
  }

/******************************* readimagehead *******************************/
/*
extract some data from the FITS-file header
*/
void	readimagehead(picstruct *field)
  {
   int          status = SAI__OK;
   int          ndims, dims[NDF__MXDIM];
   char         type[20];
   void         *pntr[3];
   int          nel;
   int          placehldr;
   int          lbnd[NDF__MXDIM];
   int          ubnd[NDF__MXDIM];
   int          sigaxis[2];
   int          nsig;
   int          i;
   int          exists;


/* Open the file */
  field->file = 0;
  ndfOpen( NULL, field->filename, "READ", "OLD",
           &field->ndf, &placehldr, &status );

  if (status != SAI__OK)
    error(EXIT_FAILURE,"*Error*: Failed to open ", field->filename);

/*---------------------------- Basic keywords ------------------------------*/
  ndfType( field->ndf, "DATA", type, 20, &status );
  if (!(status == SAI__OK))
    error(EXIT_FAILURE,"*Error*: Failed to get data type. ", field->filename);

  field->bitsgn = 1;
  if (!strcmp(type,"_BYTE")) {
      field->bitpix = BP_BYTE;
  } else if (!strcmp(type,"_UBYTE")) {
    field->bitpix = BP_BYTE;
    field->bitsgn = 0;
  } else if (!strcmp(type,"_WORD")) {
    field->bitpix = BP_SHORT;
  } else if (!strcmp(type,"_UWORD")) {
    field->bitpix = BP_SHORT;
    field->bitsgn = 0;
  } else if (!strcmp(type,"_INTEGER")) {
    field->bitpix = BP_LONG;
  } else if (!strcmp(type,"_REAL")) {
    field->bitpix = BP_FLOAT;
  } else if (!strcmp(type,"_DOUBLE")) {
    field->bitpix = BP_DOUBLE;
  } else error(EXIT_FAILURE, "Sorry, I don't know that kind of data.", "");

  field->bytepix = (field->bitpix>0?field->bitpix:-field->bitpix)>>3;

  ndfDim( field->ndf, NDF__MXDIM, dims, &ndims, &status );
  sigaxis[0] = 0;
  sigaxis[1] = 1;
  if ( ndims != 2 ) {
      /* Look for 2D, but with insignificant dimensions */
      nsig = 0;
      for( i = 0; i < ndims; i++ ) {
          if( dims[ i ] > 1 ) {
              if ( nsig == 0 ) {
                  sigaxis[0] = i;
              }
              else {
                  sigaxis[1] = i;
              }
              nsig++;
          }
      }
      if ( nsig != 2 ) {
          error( EXIT_FAILURE, field->filename, " does NOT contain 2D data." );
      }
  }

  field->width = dims[sigaxis[0]];
  field->height = dims[sigaxis[1]];
  field->npix = (KINGSIZE_T)field->width*field->height;

  field->bscale = 1.0;
  field->bzero = 0.0;
  if (field->bitsgn && prefs.fitsunsigned_flag)
    field->bitsgn = 0;

  ndfCget( field->ndf, "TITLE", field->ident, MAXCHAR, &status );

  ndfBound( field->ndf, NDF__MXDIM, lbnd, ubnd, &ndims, &status );
  field->origin[0] = lbnd[sigaxis[0]];
  field->origin[1] = lbnd[sigaxis[1]];

  /*  Look for the FITS header and map, note never unmap, let
   *  annul clear this resource. */
  ndfXstat( field->ndf, "FITS", &exists, &status );
  field->fitsheadsize = 0;
  if ( exists ) {
      HDSLoc *fitsloc = NULL;
      size_t nhead = 0;
      char *fitshead;
      ndfXloc( field->ndf, "FITS", "READ", &fitsloc, &status );
      datMapV( fitsloc, "_CHAR*80", "READ", (void**) &fitshead,
               &nhead, &status );
      field->fitsheadsize = nhead;

      /*  Make sure we have an "END" card. NDFs may not add these. */
      if ( strncmp( fitshead+(80*(nhead-1)), "END     ", 8 ) != 0 ) {
          QMEMCPY2( fitshead, field->fitshead, char, (nhead + 1)*80,
                    nhead*80 );
          strncpy( field->fitshead+(80*nhead), "END     ", 8 );
          field->fitsheadsize = nhead + 1;
      }
      else {
          QMEMCPY2( fitshead, field->fitshead, char, nhead*80, nhead*80 );
      }
  }
  else {
      field->fitshead = "END                                         ";
      field->fitsheadsize = 1;
  }


/*----------------------------- Astrometry ---------------------------------*/
/* Presently, astrometry is done only on the measurement and detect images */
  if ( field->flags & ( MEASURE_FIELD | DETECT_FIELD ) ) {
    AstFrame     *baseframe;
    AstFrame     *frame;
    AstFrameSet  *fs;
    AstFrameSet  *wcsinfo;
    AstMapping   *map;
    AstSkyFrame  *template;
    struct wcs   *wcs;
    int          base;
    int          current;
    int          outperm[2];

    QCALLOC(wcs, struct wcs, 1);
    field->wcs = wcs;

    /* See if the NDF has a WCS component */
    status = SAI__OK;
    wcs->astwcs = NULL;
    ndfState( field->ndf, "WCS", &exists, &status );
    if ( exists ) {

        /* Get a pointer to the WCS FrameSet. */
        ndfGtwcs( field->ndf, &wcsinfo, &status );

        /* We can only deal with 2D celestial coordinate systems in which the
           longitude is the first axis and the latitude is the second
           axis. Create such a Frame (a SkyFrame) which we can use as a
           template for probing the WCS FrameSet. All attributes which are
           left unset (such as System) will act as wild cards and match any
           value in the WCS FrameSet.
        */
        template = astSkyFrame( " " );

        /* Allow the find to pick out a SkyFrame that's inside a nD frame,
         * most likely this matches the current frame dimensionality */
        astSetI( template, "MaxAxes", astGetI( wcsinfo, "Nout" ) );

        /* Search the WCS FrameSet for a Frame matching this template. This
           will match any SkyFrame, no matter what the axis order, system type,
           etc. It returns a FrameSet connecting the base Frame in the WCS
           FrameSet (GRID coords) to a Frame which has the axis order of the
           template but inherits attribute values from the matching (Sky)Frame
           in the WCS FrameSet.
        */
        fs = astFindFrame( wcsinfo, template, " " );
        template = astAnnul( template );

        /* A NULL pointer is returned if no SkyFrame is found in the WCS
           FrameSet. */
        if ( fs ) {

            /* Add the resulting SkyFrame into the WCS FrameSet, making it the
               new current Frame. */
            map = astGetMapping( fs, AST__BASE, AST__CURRENT );
            frame = astGetFrame( fs, AST__CURRENT );

            astAddFrame( wcsinfo, AST__BASE, map, frame );
            map = astAnnul( map );
            frame = astAnnul( frame );

            /*  Longtitude and latitude, zero based axes */
            wcs->lat = astGetI( wcsinfo, "LatAxis" ) - 1;
            wcs->lng = astGetI( wcsinfo, "LonAxis" ) - 1;

            /* Equinox and epoch. */
            wcs->equinox = astGetD( wcsinfo, "Equinox" );
            wcs->epoch = astGetD( wcsinfo, "Epoch" );
        }
        else {
            wcs->lat = 0;
            wcs->lng = 0;
        }

        /* The BASE coordinates should be 2D as well. Make sure that's the
         * case. Pick out the select sigaxes. */
        if ( ndims != 2 ) {
            outperm[0] = sigaxis[0] + 1;
            outperm[1] = sigaxis[1] + 1;
            baseframe = astGetFrame( fs, AST__BASE );
            frame = astPickAxes( baseframe, 2, outperm, &map );
            baseframe = astAnnul( baseframe );

            /* Now add this frame to the FrameSet and make it the base
             * one. Also reinstate the skyframe as the current frame.*/
            current = astGetI( wcsinfo, "Current" );
            astAddFrame( wcsinfo, AST__BASE, map, frame );
            base = astGetI( wcsinfo, "Current" );
            astSetI( wcsinfo, "Base", base );
            astSetI( wcsinfo, "Current", current );
            frame = astAnnul( frame );
            map = astAnnul( map );
        }

        /* Store the pointer to the WCS FrameSet */
        wcs->naxis = 2;
        wcs->astwcs = wcsinfo;

        /* Need a default pixel scale, only sensible for Sky domains. */
        field->pixscale = 1.0;
        if ( wcs->lat != wcs->lng ) {
            double point1[2], point2[2];
            double xin[2], yin[2], xout[2], yout[2];
            double dist;
            double xcen, ycen;

            /* Compute the scales the the sizes of a pixel near the centre of
             * the image. */
            xcen = 0.5 * ( (double) field->width );
            ycen = 0.5 * ( (double) field->height );
            xin[0] = xcen - 0.5;
            xin[1] = xcen + 0.5;
            yin[0] = yin[1] = ycen;

            /* Transform these image positions into sky coordinates. */
            astTran2( wcsinfo, 2, xin, yin, 1, xout, yout );

            /* And now get the distance between these positions in degrees. */
            point1[0] = xout[0];
            point1[1] = yout[0];
            point2[0] = xout[1];
            point2[1] = yout[1];
            dist = astDistance( wcsinfo, point1, point2 );
            if ( ! astOK ) astClearStatus;
            if ( dist != AST__BAD )  {
                field->pixscale = dist / DEG * 3600.0;
            }
        }

        if ( status != SAI__OK ) errAnnul( &status );
    }
    else {
      /*  No WCS */
      wcs->lat = 0;
      wcs->lng = 0;
      wcs->epoch = 2000.0;
      wcs->equinox = 2000.0;
    }

    /*-----------------------------------------------------------------------*/
  } /* end of MEASURE or DETECT field */

  /*-------------------------------------------------------------------------*/

  /* Map the NDF in the appropriate type */
  if (field->flags ^ FLAG_FIELD) {
      ndfMap( field->ndf, "DATA", type, "READ", pntr, &nel, &status );
  }
  else {
      /* All FLAG images are used with type "unsigned int" */
      ndfMap( field->ndf, "DATA", "_INTEGER", "READ", pntr, &nel, &status );
  }

  field->map = pntr[0];
  field->file = 0;
  field->nel = nel;

  return;
  }
