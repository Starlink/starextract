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

    /* Write out the string via the message system */
    va_start( args, fmt );
    msgOutifv( MSG__QUIET, " ", string, args, &status );
    va_end( args );

}
