#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#include "sae_par.h"
#include "mers.h"

/*
 * adamprint - convert a C-style print arguments to ADAM MSG output
 *
 * Uses va_list arguments as expected by the format statement.
 *
 */
void adamprint( FILE *file, char *fmt, ... )
{
    va_list args;
    int status = SAI__OK;
    char *fmtcpy;
    int i;
    int j;
    int len;

    /*  Strip out any 4 character escape sequences from format string.
     *  These are added by SExtractor and shouldn't be passed through the 
     *  message system. */
    len = strlen( fmt );
    fmtcpy = (char *) malloc( (size_t) len + 1 );
    for ( i = 0, j = 0; i < len; i++, j++ ) {
        if ( fmt[i] == '\33' ) i += 4;
        fmtcpy[j] = fmt[i];
    }
    fmtcpy[j]= '\0';

    /* Write out the string via the message system */
    va_start( args, fmt );
    msgOutifv( MSG__QUIET, " ", fmtcpy, args, &status );
    va_end( args );
    free( fmtcpy );

}
