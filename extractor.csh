#!/bin/csh
#+          
#  Name:
#     extractor.csh

#  Purpose:
#     Set up aliases for the EXTRACTOR package.

#  Type of Module:
#     C shell script.

#  Invocation:
#     source extractor.csh

#  Description:
#     This procedure defines an alias for each EXTRACTOR command. The 
#     string install_bin (upper-case) is replaced by the path of the 
#     directory containing the package executable files when the package
#     is installed.  The string help_dir (upper-case) is likewise replaced
#     by the path to the directory containing the help files.

#  Authors:
#     BLY: M.J. Bly (Starlink, RAL)
#     AJC: A.J. Chipperfield (Starlink, RAL)
#     {enter_new_authors_here}

#  History:
#     23-JUN-1995 (BLY):
#       Original Version.
#     12-DEC-1996 (BLY):
#       Cosmetic mods.
#     23-NOV-1998 (AJC):
#       Modify for extractor
#     {enter_changes_here}

#-

#  Prepare to run ADAM applications if this has not been done already.
#  ===================================================================
#
#  Here look to see if there is an ADAM_USER directory.  If there is not
#  check whether or not there is an adam file that is not a directory.
#  If there is, issue a warning and exit.  Otherwise create the required
#  directory.

if (-d ${HOME}/adam) then
   echo -n
else
   if (-f ${HOME}/adam) then
      echo "You have a file called adam in your home directory.  Please rename "
      echo "since adam must be a directory for ADAM files."
      exit
   else
      mkdir ${HOME}/adam
   endif
endif

#
#  Set up an environment variable pointing to the help library. 
#  This is refered to within the appliation interface files.

setenv EXTRACTOR_HELP INSTALL_HELP/extractor

#
#  Locate the installed binaries, scripts etc.

#setenv EXTRACTOR_DIR INSTALL_BIN
setenv EXTRACTOR_CONFIG INSTALL_BIN/config

#
#  Define symbols for the applications and scripts.
#  ===============================================

# eg:  alias command ${EXTRACTOR_BIN}/command
alias extract INSTALL_BIN/extractor

#
#  Now do the same with alternative names.
#  ======================================

# eg:  alias extractor_command ${EXTRACTOR_DIR}/command

#
#  Tell the user that EXTRACTOR commands are now available.
#  =======================================================

echo ""
echo "   EXTRACTOR commands are now available -- (Version PKG_VERS)"
echo " "
#echo "   Type extractorhelp for help on EXTRACTOR commands"
#echo " "

#
# end
