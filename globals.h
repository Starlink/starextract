 /*
 				globals.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global declarations.
*
*	Last modify:   25/05/99 (EB):
*                      02/11/98 (AJC)
*                        Add initialise flag
*	Last modify:   02/04/2003
*                        (EB): 2.3
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	"types.h"

/*----------------------- initialization flags ------------------------------
 * These are required for the ADAM version so that the need for initialization
 * (or not) of static and external variables can be flagged when an executable
 * image is run repeatedly.
*/
int                     initparcelout; /* initialization in parcelout */
int                     initprintf; /* Initialization in initcat */

/*----------------------- miscellaneous variables ---------------------------*/

sexcatstruct		thecat;
picstruct		thefield1,thefield2, thewfield1,thewfield2;
objstruct		flagobj;
obj2struct		flagobj2;
extern obj2struct	outobj2;
static obj2struct	*obj2 = &outobj2;
float			ctg[37], stg[37];
char			gstr[MAXCHAR];

/*------------------------------- functions ---------------------------------*/
extern void	allocparcelout(void),
		analyse(picstruct *, picstruct *, int, objliststruct *),
		blankit(char *, int),
                endcat(void),
                reendcat(void),
                closecheck(void),
		copydata(picstruct *, int, int),
		endfield(picstruct *),
		endobject(picstruct *, picstruct *, picstruct *, picstruct *,
			int, objliststruct *),
		examineiso(picstruct *, picstruct *, objstruct *,
			pliststruct *),
		flagcleancrowded(int, objliststruct *),
		freeparcelout(void),
		getnnw(void),
		initcat(void),
		reinitcat(picstruct *),
		initglob(void),
		makeit(void),
		mergeobject(objstruct *, objstruct *),
		neurinit(void),
		neurclose(void),
		neurresp(double *, double *),
		preanalyse(int, objliststruct *, int),
		readcatparams(char *),
		readdata(picstruct *, PIXTYPE *, int),
		readidata(picstruct *, FLAGTYPE *, int),
		readimagehead(picstruct *),
		readprefs(char *, char **, char **, int),
		scanimage(picstruct *, picstruct *, picstruct **, int,
			picstruct *, picstruct *),
		sexcircle(PIXTYPE *bmp, int, int, double, double, double,
			PIXTYPE),
		sexdraw(PIXTYPE *bmp, int, int, double, double, PIXTYPE),
		sexellips(PIXTYPE *bmp, int, int, double, double, double,
			double, double, PIXTYPE, int),
		sexmove(double, double),
		updateparamflags(void),
		useprefs(void),
		writecat(int, objliststruct *);

extern float	hmedian(float *, int);

extern int	addobj(int, objliststruct *, objliststruct *),
		belong(int, objliststruct *, int, objliststruct *),
		gatherup(objliststruct *, objliststruct *),
		parcelout(objliststruct *, objliststruct *);

extern void	*loadstrip(picstruct *, picstruct *);

extern char	*readfitshead(FILE *, char *, int *);

extern picstruct	*inheritfield(picstruct *infield, int flags),
			*newfield(char *, int , int);

