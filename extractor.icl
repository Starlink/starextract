{+          
{  Name:
{     extractor.icl

{  Purpose:
{     Set up commands for the EXTRACTOR package.

{  Type of Module:
{     ICL command file.

{  Invocation:
{     load extractor.icl

{  Description:
{     This procedure defines an command for each EXTRACTOR command, and 
{     also defines the help commands.

{  Authors:
{     BLY: M.J. Bly (Starlink, RAL)
{     AJC: A.J. Chipperfield (Starlink, RAL)
{     {enter_new_authors_here}

{  History:
{     12-DEC-1996 (BLY):
{       Original Version, based on `csh' model.
{     23-NOV-1998 (AJC):
{       Modified for EXTRACTOR
{     {enter_changes_here}

{-

{  Define main package help
 
defhelp extractor $EXTRACTOR_HELP 0
defstring extractorhelp help extractor
 
{  Basic command definitions.
{  eg:  define command $EXTRACTOR_DIR/extractor_mon command

{  Now do the same with alternative names.
{  eg:  define {pkg}_command $EXTRACTOR_DIR/extractor_mon command

define extractor $EXTRACTOR_DIR/extractor

{
{  Tell the user that EXTRACTOR commands are now available.

print " "
print "   EXTRACTOR commands are now available -- (Version PKG_VERS)"
print " "
{print "   Type `extractorhelp' or `help extractor' for help on EXTRACTOR commands"
{print " "

{
{  end
{.
