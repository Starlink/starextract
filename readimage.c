/*
 				readimage.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	ADAM SExtractor
*
*	Author:		A. Chipperfield (Starlink, CCLRC, RAL)
*                       P.W. Draper (Starlink, Durham University)
*                       E. Bertin (IAP, Leiden observatory & ESO)
*
*	Contents:	functions for input of image data.
*
*	Last modify:	06/05/99 (EB):
*	Last modify:	28/10/98 (AJC)
*                          In line with V2.0.15
*                       14/12/98 (PWD):
*                          Added USHORT and UBYTE support.
*                       17/12/98 (PWD):
*                          Changed to use NDF WCS component for astrometry.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#include	"interpolate.h"
#include	"back.h"
#include	"astrom.h"
#include	"weight.h"

#include        "sae_par.h"
#include        "dat_par.h"
#include        "ndf.h"
#include        "merswrap.h"

/******************************* readimagehead *******************************/
/*
extract some data from the FITS-file header
*/
void	readimagehead(picstruct *field)
  {
   int          status = SAI__OK;
   int		i,j,l, n;
   int          ndims, dims[NDF__MXDIM];
   int          there;
   char		*buf, st[80], str[80], *point;
   char         type[20];
   void         *pntr[3];
   int          nel;
   int          placehldr;
   int          lbnd[NDF__MXDIM];
   int          ubnd[NDF__MXDIM];

/* Open the file */
  field->file = 0;
  ndfOpen( DAT__ROOT, field->filename, "READ", "OLD",
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

/*  field->bytepix = (field->bitpix>0?field->bitpix:-field->bitpix)>>3;*/
  field->bytepix = 4;

  ndfDim( field->ndf, NDF__MXDIM, dims, &ndims, &status );
  if ( !(ndims==2))
    error( EXIT_FAILURE, field->filename, " does NOT contain 2D data." );

  field->width = dims[0];
  field->height = dims[1];
  field->npix = (KINGSIZE_T)field->width*field->height;

  field->bscale = 1.0;
  field->bzero = 0.0;
  if (field->bitsgn && prefs.fitsunsigned_flag)
    field->bitsgn = 0;

  ndfCget( field->ndf, "TITLE", field->ident, MAXCHAR, &status );

  ndfBound( field->ndf, NDF__MXDIM, lbnd, ubnd, &ndims, &status );
  field->origin[0] = lbnd[0];
  field->origin[1] = lbnd[1];

/*----------------------------- Astrometry ---------------------------------*/
/* Presently, astrometry is done only on the measurement and detect images */
  if ( field->flags & ( MEASURE_FIELD | DETECT_FIELD ) ) {
    astromstruct *as;
    AstFrameSet  *wcsinfo;
    AstFrame     *frame;
    int          exists;
    int          issky;

    QCALLOC(as, astromstruct, 1);
    field->astrom = as;

    /* See if the NDF has a WCS component */
    status = SAI__OK;
    field->astwcs = NULL;
    ndfState( field->ndf, "WCS", &exists, &status );
    if ( exists ) {

      /* Access it and record the identifier */
      ndfGtwcs( field->ndf, &wcsinfo, &status );
      as->naxis = 2;
      field->astwcs = wcsinfo;

      /* Make sure the WCS represents a celestial coordinate system
         (that's all we can deal with). */
      frame = astGetFrame( wcsinfo, AST__CURRENT );
      issky = astIsASkyFrame( frame );
      frame = astAnnul( frame );
      if ( ! issky ) {
        /*  No WCS */
        as->wcs_flag = 0;
      } else {
        as->wcs_flag = 1;
      }
      if ( status != SAI__OK ) errAnnul( &status );
    }
    else {
      /*  No WCS */
      as->wcs_flag = 0;
    }

    /*-----------------------------------------------------------------------*/
  } /* end if MEASURE or DETECT field */

  /*-------------------------------------------------------------------------*/

  /* Map the NDF in the appropriate type */
    if (field->flags ^ FLAG_FIELD)
      ndfMap( field->ndf, "DATA", "_REAL", "READ", pntr, &nel, &status );
    else
      ndfMap( field->ndf, "DATA", "_REAL", "READ", pntr, &nel, &status );
    field->map = pntr[0];
    field->file = 0;
    field->nel = nel;

  return;
  }

/******************************* loadstrip ***********************************/
/*
Load a new strip of pixel data into the buffer.
*/
void	*loadstrip(picstruct *field, picstruct *wfield)

  {
   checkstruct	*check;
   int		y, w, flags, interpflag;
   PIXTYPE	*data, *wdata, *rmsdata;

  w = field->width;
  flags = field->flags;
  interpflag = (wfield && wfield->interp_flag);

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
            if (check = prefs.check[CHECK_BACKGROUND])
              writecheck(check, field->backline, w);
            if (check = prefs.check[CHECK_SUBTRACTED])
              writecheck(check, data, w);
            if (check = prefs.check[CHECK_APERTURES])
              writecheck(check, data, w);
            if (check = prefs.check[CHECK_SUBPSFPROTOS])
              writecheck(check, data, w);
            if (check = prefs.check[CHECK_SUBPCPROTOS])
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
          if (check = prefs.check[CHECK_BACKGROUND])
            writecheck(check, field->backline, w);
          if (check = prefs.check[CHECK_SUBTRACTED])
            writecheck(check, data, w);
          if (check = prefs.check[CHECK_APERTURES])
            writecheck(check, data, w);
          if (check = prefs.check[CHECK_SUBPSFPROTOS])
            writecheck(check, data, w);
          if (check = prefs.check[CHECK_SUBPCPROTOS])
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
  for (i=0; i<(size<left?size:left); i++)
     *(ptr++) = *((PIXTYPE *)field->map + field->file++)*bs + bz;

/* Reset field->file if at end */
/* This prevents endfield from crashing if field->file is used for NDF */
  field->file = (field->file >= field->nel)?0:field->file;
  return;
  }


/******************************** readidata **********************************/
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
