/*
 				fitswrite.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	The LDAC Tools
*
*	Author:		E.BERTIN, DeNIS/LDAC
*
*	Contents:	low-level functions for writing LDAC FITS catalogs.
*
*	Last modify:	13/12/2002
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef	HAVE_CONFIG_H
#include "config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"fitscat_defs.h"
#include	"fitscat.h"

char	*lineout_buf, padbuf[FBSIZE];
int	lineout_size, nlineout;

/****** save_cat ***************************************************************
PROTO	void save_cat(catstruct *cat, char *filename)
PURPOSE	Save a FITS catalog with name filename.
INPUT	catalog structure,
	filename.
OUTPUT	-.
NOTES	Any preexisting file with name filename is overwritten.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	08/05/2002
 ***/
void	save_cat(catstruct *cat, char *filename)

  {
   tabstruct	*tab;
   int		i;

  strcpy(cat->filename, filename);
  if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);

  tab = cat->tab;
/*Go through each segment in the right order to save data*/
  for (i=0; i<cat->ntab; i++)
    {
/*-- Make sure that the tab header is primary or extension! */
    if (i)
      ext_head(tab);
    else
      prim_head(tab);
    save_tab(cat, tab);
    while (!((tab=tab->nexttab)->nseg));
    }

  if (close_cat(cat) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: Problem while closing", cat->filename);

  return;
  }


/****** save_tab **************************************************************
PROTO	void save_tab(catstruct *cat, tabstruct *tab)
PURPOSE	Save a FITS catalog with name filename.
INPUT	catalog structure,
	filename.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	13/12/2002
 ***/
void	save_tab(catstruct *cat, tabstruct *tab)

  {
   catstruct	*tabcat;
   keystruct	*key;
   tabstruct	*keytab;
   KINGSIZE_T	tabsize;
   KINGLONG	size;
   int		b,j,k,o, nbytes,nkey,nobj,spoonful,
		tabflag, larrayin,larrayout;
   char		*buf, *inbuf, *outbuf, *fptr,*ptr;
   int		esize;

/*  Make the table parameters reflect its content*/
  update_tab(tab);
/*  The header itself*/
  tabflag = update_head(tab)==RETURN_OK?1:0;
  QFTELL(cat->file, tab->headpos, cat->filename);
  QFWRITE(tab->headbuf, tab->headnblock*FBSIZE, cat->file, cat->filename);
/*  Allocate memory for the output buffer */
  tabsize = 0;
  if (tabflag)
    {
/*-- If segment is a binary table, save it row by row */
    QMALLOC(outbuf, char, (larrayout = tab->naxisn[0]));
    nkey = tab->nkey;
    tabsize = 0;
    for (j=tab->nseg; j--;)
      {
      update_tab(tab);
/*---- Scan keys to find the reference tab and other things*/
      keytab = NULL;
      key = tab->key;
      for (k=nkey; k--; key = key->nextkey)
        if (!key->ptr)
          {
          keytab = key->tab;
          tabcat = keytab->cat;
          }
/*---- If table contains some keys with no ptrs, we have to access a file */
      if (keytab)
        {
        QMALLOC(inbuf, char, (larrayin = keytab->naxisn[0]));
        if (open_cat(tabcat, READ_ONLY) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Cannot access ", tabcat->filename);
        QFSEEK(tabcat->file, keytab->bodypos, SEEK_SET, tabcat->filename);
        }
      nobj = tab->naxisn[1];
      for (o=0; o<nobj; o++)
        {
        if (keytab)
          QFREAD(inbuf, larrayin, tabcat->file, tabcat->filename);
        fptr = outbuf;
        for (k=nkey; k--; key = key->nextkey)
          {
          nbytes = key->nbytes;
          ptr = key->ptr? (char *)key->ptr+nbytes*o:inbuf+key->pos;
          for (b=nbytes; b--;)
            *(fptr++) = *(ptr++);
          if (bswapflag)
            if (key->ptr)
              {
              esize = t_size[key->ttype];
              swapbytes(fptr-nbytes, esize, nbytes/esize);
              }
          }
        QFWRITE(outbuf, larrayout, cat->file, cat->filename);
        }
      if (keytab)
        {
        free(inbuf);
        if (close_cat(tabcat) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Problem while closing",
		tabcat->filename);
        }
      tabsize += tab->tabsize;
      tab = tab->nexttab;
      }
    free(outbuf);
    }
  else
    {
/*-- If segment is not a binary table, save it ``as it is'' */
/*-- We use a limited-size buffer ``in case of'' */
    size = tabsize = tab->tabsize;
    if (tabsize)
      {
      if (tab->bodybuf)
        {
/*------ A body is present in memory and needs to be written */
        if (bswapflag)
          swapbytes(tab->bodybuf, tab->bytepix, tabsize/tab->bytepix);
        QFWRITE(tab->bodybuf, (size_t)tabsize, cat->file, cat->filename);
        if (bswapflag)
          swapbytes(tab->bodybuf, tab->bytepix, tabsize/tab->bytepix);
        }
      else
/*------ The body should be copied from the source tab */
        {
        tabcat = tab->cat;
        spoonful = size<DATA_BUFSIZE?size:DATA_BUFSIZE;
        QMALLOC(buf, char, spoonful);
        if (open_cat(tabcat, READ_ONLY) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Cannot access ", tabcat->filename);
        QFSEEK(tabcat->file, tab->bodypos, SEEK_SET, tabcat->filename);
        for (;size>0; size -= spoonful)
          {
          if (spoonful>size)
          spoonful = size;
          QFREAD(buf, spoonful, tabcat->file, tabcat->filename);
          QFWRITE(buf, spoonful, cat->file, cat->filename);
          }
        free(buf);
        if (close_cat(tabcat) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Problem while closing",
		tabcat->filename);
        }
      }
    }

/* FITS padding*/
  size = PADEXTRA(tabsize);
  if (size)
    QFWRITE(padbuf, (size_t)size, cat->file, cat->filename);

  return;
  }


/****** init_writeobj *********************************************************
PROTO	void init_writeobj(catstruct *cat, tabstruct *tab)
PURPOSE	Prepare the writing of individual sources in a FITS table
INPUT	catalog structure,
	table structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	13/12/2002
 ***/
void	init_writeobj(catstruct *cat, tabstruct *tab)

  {
   keystruct	*key;
   int		k;

/* Make the table parameters reflect its content*/
  update_tab(tab);
/* The header itself */
  if (update_head(tab) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: Not a binary table: ", tab->extname);

  QFTELL(cat->file, tab->headpos, cat->filename);
  QFWRITE(tab->headbuf, tab->headnblock*FBSIZE, cat->file, cat->filename);
/* Store the current position */
  QFTELL(cat->file, tab->bodypos, cat->filename);

/* Allocate memory for the output buffer  (or increase it if done already) */
  if (lineout_buf)
    {
    if (lineout_size < tab->naxisn[0])
      {
      QREALLOC(lineout_buf, char, tab->naxisn[0]);
      lineout_size = tab->naxisn[0];
      }
    }
  else
    {
    QMALLOC(lineout_buf, char, tab->naxisn[0]);
    lineout_size = tab->naxisn[0];
    }

  nlineout++;

/* Scan keys to check that memory has been allocated */
  key = tab->key;
  for (k=tab->nkey; k--; key = key->nextkey)
    if (!key->ptr)
      error(EXIT_FAILURE, "*Error*: no memory allocated for ", key->name);

/* No object written yet: initialize counter */
  tab->naxisn[1] = 0;

  return;
  }


/****** write_obj *************************************************************
PROTO	int write_obj(tabstruct *tab)
PURPOSE	Write one individual source in a FITS table
INPUT	Destination catalog
	Table structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	18/02/2000
 ***/
int	write_obj(tabstruct *tab)

  {
   keystruct	*key;
   char		*pin, *pout;
   int		b,k;
   int		esize;

  key = tab->key;
  pout = lineout_buf;
  for (k=tab->nkey; k--; key = key->nextkey)
    {
    pin = key->ptr;
    if (bswapflag)
      {
      esize = t_size[key->ttype];
      swapbytes(pin, esize, key->nbytes/esize);
      }
    for (b=key->nbytes; b--;)
      *(pout++) = *(pin++);
    }

  QFWRITE(lineout_buf, *tab->naxisn, tab->cat->file, tab->cat->filename);

  return ++tab->naxisn[1];
  }


/****** end_writeobj **********************************************************
PROTO	void end_writeobj(catstruct *cat, tabstruct *tab)
PURPOSE	End the writing of individual sources in a FITS table
INPUT	catalog structure,
	table structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	13/12/2002
 ***/
void	end_writeobj(catstruct *cat, tabstruct *tab)

  {
   keystruct	*key;
   OFF_T	pos;
   int		k, size;

/* Make the table parameters reflect its content*/
  key = tab->key;
  for (k=tab->nkey; k--; key = key->nextkey)
    key->nobj = tab->naxisn[1];
  update_tab(tab);
/* The header itself */
  if (update_head(tab) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: Not a binary table: ", tab->extname);

/*--FITS padding*/
  size = PADEXTRA(tab->tabsize);
  if (size)
    QFWRITE(padbuf, size, cat->file, cat->filename);
  QFTELL(cat->file, pos, cat->filename);
  QFSEEK(cat->file, tab->headpos, SEEK_SET, cat->filename);
  QFWRITE(tab->headbuf, FBSIZE*tab->headnblock, cat->file, cat->filename);
  QFSEEK(cat->file, pos, SEEK_SET, cat->filename);

  if (!(--nlineout))
    {
    free(lineout_buf);
    lineout_buf = NULL;
    lineout_size = 0;
    }

  return;
  }


/****** print_obj *************************************************************
PROTO	void print_obj(FILE *stream, tabstruct *tab)
PURPOSE	Print one individual source to the output stream.
INPUT	Output stream
	Table structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	13/06/97
 ***/
void	print_obj(FILE *stream, tabstruct *tab)

  {
   keystruct	*key;
   char		*ptr;
   int		i,k, esize;

  if (!(key = tab->key))
    error(EXIT_FAILURE, "*Error*: no key to print in table ", tab->extname);

  for (k=tab->nkey; k--; key = key->nextkey)
    {
    esize = t_size[key->ttype];
    ptr = key->ptr;
    for (i = key->nbytes/esize; i--; ptr += esize)
      switch(key->ttype)
        {
        case T_SHORT:
          fprintf(stream, *key->printf?key->printf:"%d", *(short *)ptr);
          if (i)
            putc(' ', stream);
          break;
        case T_LONG:
          fprintf(stream, *key->printf?key->printf:"%d", *(int *)ptr);
          if (i)
            putc(' ', stream);
          break;
        case T_FLOAT:
          fprintf(stream, *key->printf?key->printf:"%g", *(float *)ptr);
          if (i)
            putc(' ', stream);
          break;
        case T_DOUBLE:
          fprintf(stream, *key->printf?key->printf:"%f", *(double *)ptr);
          if (i)
            putc(' ', stream);
          break;
        case T_BYTE:
          if (key->htype==H_BOOL)
            {
            if (*ptr)
              fprintf(stream, "T");
            else
              fprintf(stream, "F");
            }
          else
            fprintf(stream, key->printf?key->printf:"%d", (int)*ptr);
          if (i)
            putc(' ', stream);
          break;
        case T_STRING:
          fprintf(stream, "%c", (int)*ptr);
          break;
        default:
          error(EXIT_FAILURE, "*FATAL ERROR*: Unknown FITS type in ",
		"show_keys()");
        }
    if (k)
      putc(' ', stream);
    }

  putc('\n', stream);

  return;
  }

