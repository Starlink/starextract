/*
 				param.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*                       P.W.DRAPER (STARLINK)
*
*	Contents:	parameter list for catalog data.
*
*	Last modify:	20/07/99
*                       12/06/01 PWD: Changed various _WORLD, degree
*                       units, coordinate outputs to use format "%15.8g"
*                       (from %15e). This gets to 0.01 precision in
*                       arcsec.
*                       19/07/01 PWD: added BKGSIG for GIM2D people.
*                       20/01/02 PWD: added X_PIXEL and Y_PIXEL for
*                                     NDF pixel coordinates.
*	Last modify:	16/12/2002
*                       (EB): 2.3
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

objstruct	outobj;
obj2struct	outobj2;

/*--------------------------------- initialization --------------------------*/
keystruct	objkey[] = {
  {"NUMBER", "Running object number",
	&outobj.number, H_INT, T_LONG, "%10d", ""},
  {"EXT_NUMBER", "FITS extension number",
	&outobj2.ext_number, H_INT, T_SHORT, "%3d", ""},
  {"FLUX_ISO", "Isophotal flux",
	&outobj2.flux_iso, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUXERR_ISO", "RMS error for isophotal flux",
	&outobj2.fluxerr_iso, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"MAG_ISO", "Isophotal magnitude",
	&outobj2.mag_iso, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR_ISO", "RMS error for isophotal magnitude",
	&outobj2.magerr_iso, H_FLOAT, T_FLOAT, "%8.4f", "mag"},

  {"FLUX_ISOCOR", "Corrected isophotal flux",
	&outobj2.flux_isocor, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUXERR_ISOCOR", "RMS error for corrected isophotal flux",
	&outobj2.fluxerr_isocor, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"MAG_ISOCOR", "Corrected isophotal magnitude",
	&outobj2.mag_isocor, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR_ISOCOR", "RMS error for corrected isophotal magnitude",
	&outobj2.magerr_isocor, H_FLOAT, T_FLOAT, "%8.4f", "mag"},

  {"FLUX_APER", "Flux vector within fixed circular aperture(s)",
	&outobj2.flux_aper, H_FLOAT, T_FLOAT, "%12g", "count", 1,
	&prefs.flux_apersize},
  {"FLUXERR_APER", "RMS error vector for aperture flux(es)",
	&outobj2.fluxerr_aper, H_FLOAT, T_FLOAT, "%12g", "count", 1,
	&prefs.fluxerr_apersize},
  {"MAG_APER", "Fixed aperture magnitude vector",
	&outobj2.mag_aper, H_FLOAT, T_FLOAT, "%8.4f", "mag", 1,
	&prefs.mag_apersize},
  {"MAGERR_APER", "RMS error vector for fixed aperture mag.",
	&outobj2.magerr_aper, H_FLOAT, T_FLOAT, "%8.4f", "mag", 1,
	&prefs.magerr_apersize},

  {"FLUX_AUTO", "Flux within a Kron-like elliptical aperture",
	&outobj2.flux_auto, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUXERR_AUTO", "RMS error for AUTO flux",
	&outobj2.fluxerr_auto, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"MAG_AUTO", "Kron-like elliptical aperture magnitude",
	&outobj2.mag_auto, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR_AUTO", "RMS error for AUTO magnitude",
	&outobj2.magerr_auto, H_FLOAT, T_FLOAT, "%8.4f", "mag"},

  {"FLUX_BEST", "Best of FLUX_AUTO and FLUX_ISOCOR",
	&outobj2.flux_best, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUXERR_BEST", "RMS error for BEST flux",
	&outobj2.fluxerr_best, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"MAG_BEST", "Best of MAG_AUTO and MAG_ISOCOR",
	&outobj2.mag_best, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR_BEST", "RMS error for MAG_BEST",
	&outobj2.magerr_best, H_FLOAT, T_FLOAT, "%8.4f", "mag"},

  {"FLUX_PROFILE", "Flux weighted by the FILTERed profile",
	&outobj2.flux_prof, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUXERR_PROFILE", "RMS error for PROFILE flux",
	&outobj2.fluxerr_prof, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"MAG_PROFILE", "Magnitude weighted by the FILTERed profile",
	&outobj2.mag_prof, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR_PROFILE", "RMS error for MAG_PROFILE",
	&outobj2.magerr_prof, H_FLOAT, T_FLOAT, "%8.4f", "mag"},

  {"FLUX_SOMFIT", "Flux derived from SOM fit",
	&outobj2.flux_somfit, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUXERR_SOMFIT", "RMS error for SOMFIT flux",
	&outobj2.fluxerr_somfit, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"MAG_SOMFIT", "Magnitude derived from SOM fit",
	&outobj2.mag_somfit, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR_SOMFIT", "Magnitude error derived from SOM fit",
	&outobj2.magerr_somfit, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"ERROR_SOMFIT", "Reduced Chi-square error of the SOM fit",
	&outobj2.stderr_somfit, H_FLOAT, T_FLOAT, "%10g", ""},
  {"VECTOR_SOMFIT", "Position vector of the winning SOM node",
	&outobj2.vector_somfit, H_FLOAT, T_FLOAT, "%5.2f", "", 1,
	&prefs.somfit_vectorsize},

  {"FLUX_GALFIT", "Flux derived from the galaxy fit",
	&outobj2.flux_galfit, H_FLOAT, T_FLOAT, "%12g", "count"},
/*
  {"FLUXERR_GALFIT", "RMS error for GALFIT flux",
	&outobj2.fluxerr_galfit, H_FLOAT, T_FLOAT, "%12g", "count"},
*/
  {"MAG_GALFIT", "Magnitude derived from galaxy fit",
	&outobj2.mag_galfit, H_FLOAT, T_FLOAT, "%8.4f", "mag"},

/*
  {"MAGERR_GALFIT", "Magnitude error derived from galaxy fit",
	&outobj2.magerr_galfit, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"ERROR_GALFIT", "Reduced Chi-square error of the galaxy fit",
	&outobj2.stderr_galfit, H_FLOAT, T_FLOAT, "%10g", ""},
*/
  {"GALDANG_IMAGE", "Galaxy disk position angle  from the galaxy fit",
	&outobj2.gdposang, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"GALDSCALE_IMAGE", "Galaxy disk-scale from the galaxy fit",
	&outobj2.gdscale, H_FLOAT, T_FLOAT, "%9.3f", "pixel"},
  {"GALDASPEC_IMAGE", "Galaxy disk aspect ratio from the galaxy fit",
	&outobj2.gdaspect, H_FLOAT, T_FLOAT, "%5.3f", ""},
  {"GALDE1_IMAGE", "Galaxy disk ellipticity #1 from the galaxy fit",
	&outobj2.gde1, H_FLOAT, T_FLOAT, "%6.4f", ""},
  {"GALDE2_IMAGE", "Galaxy disk ellipticity #2 from the galaxy fit",
	&outobj2.gde2, H_FLOAT, T_FLOAT, "%6.4f", ""},
  {"GALBRATIO_IMAGE", "Galaxy bulge ratio from the galaxy fit",
	&outobj2.gbratio, H_FLOAT, T_FLOAT, "%5.3f", ""},
  {"GALBANG_IMAGE", "Galaxy bulge position angle  from the galaxy fit",
	&outobj2.gbposang, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"GALBSCALE_IMAGE", "Galaxy bulge-scale from the galaxy fit",
	&outobj2.gbscale, H_FLOAT, T_FLOAT, "%9.3f", "pixel"},
  {"GALBASPEC_IMAGE", "Galaxy bulge aspect ratio from the galaxy fit",
	&outobj2.gbaspect, H_FLOAT, T_FLOAT, "%5.3f", ""},

  {"KRON_RADIUS", "Kron apertures in units of A or B",
	&outobj2.kronfactor, H_FLOAT, T_FLOAT, "%5.2f", ""},
  {"BACKGROUND", "Background at centroid position",
	&outobj.bkg, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"BKGSIG", "Local background standard deviation",
        &outobj.sigbkg, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"THRESHOLD", "Detection threshold above background",
	&outobj.dthresh, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"FLUX_MAX", "Peak flux above background",
	&outobj.peak, H_FLOAT, T_FLOAT, "%12g", "count"},
  {"ISOAREA_IMAGE", "Isophotal area above Analysis threshold",
	&outobj.npix, H_INT, T_LONG, "%9d", "pixel**2"},
  {"ISOAREAF_IMAGE", "Isophotal area (filtered) above Detection threshold",
	&outobj.fdnpix, H_INT, T_LONG, "%9d", "pixel**2"},

  {"XMIN_IMAGE", "Minimum x-coordinate among detected pixels",
	&outobj.xmin, H_INT, T_LONG, "%10d", "pixel"},
  {"YMIN_IMAGE", "Minimum y-coordinate among detected pixels",
	&outobj.ymin, H_INT, T_LONG, "%10d", "pixel"},
  {"XMAX_IMAGE", "Maximum x-coordinate among detected pixels",
	&outobj.xmax, H_INT, T_LONG, "%10d", "pixel"},
  {"YMAX_IMAGE", "Maximum y-coordinate among detected pixels",
	&outobj.ymax, H_INT, T_LONG, "%10d", "pixel"},

  {"XPEAK_IMAGE", "x-coordinate of the brightest pixel",
	&outobj.peakx, H_INT, T_LONG, "%10d", "pixel"},
  {"YPEAK_IMAGE", "y-coordinate of the brightest pixel",
	&outobj.peaky, H_INT, T_LONG, "%10d", "pixel"},
  {"XPEAK_WORLD", "World-x coordinate of the brightest pixel",
	&outobj2.peakxw, H_FLOAT, T_DOUBLE, "%15.8g", "deg"},
  {"YPEAK_WORLD", "World-y coordinate of the brightest pixel",
	&outobj2.peakyw, H_FLOAT, T_DOUBLE, "%15.8g", "deg"},

  {"ALPHAPEAK_SKY", "Right ascension of brightest pix (native)",
	&outobj2.peakalphas, H_FLOAT, T_DOUBLE, "%11.7f", "deg"},
  {"DELTAPEAK_SKY", "Declination of brightest pix (native)",
	&outobj2.peakdeltas, H_FLOAT, T_DOUBLE, "%+11.7f", "deg"},

  {"ALPHAPEAK_J2000", "Right ascension of brightest pix (J2000)",
	&outobj2.peakalpha2000, H_FLOAT, T_DOUBLE, "%11.7f", "deg"},
  {"DELTAPEAK_J2000", "Declination of brightest pix (J2000)",
	&outobj2.peakdelta2000, H_FLOAT, T_DOUBLE, "%+11.7f", "deg"},

  {"ALPHAPEAK_B1950", "Right ascension of brightest pix (B1950)",
	&outobj2.peakalpha1950, H_FLOAT, T_DOUBLE, "%11.7f", "deg"},
  {"DELTAPEAK_B1950", "Declination of brightest pix (B1950)",
	&outobj2.peakdelta1950, H_FLOAT, T_DOUBLE, "%+11.7f", "deg"},

  {"X_IMAGE", "Object position along x",
	&outobj2.sposx, H_FLOAT, T_FLOAT, "%10.3f", "pixel"},
  {"Y_IMAGE", "Object position along y",
	&outobj2.sposy, H_FLOAT, T_FLOAT, "%10.3f", "pixel"},
  {"X_IMAGE_DBL", "Object position along x (double precision)",
	&outobj2.posx, H_FLOAT, T_DOUBLE, "%10.3f", "pixel"},
  {"Y_IMAGE_DBL", "Object position along y (double precision)",
	&outobj2.posy, H_FLOAT, T_DOUBLE, "%10.3f", "pixel"},
  {"X_WORLD", "Barycenter position along world x axis",
	&outobj2.mxw, H_FLOAT, T_DOUBLE, "%15.8g", "deg"},
  {"Y_WORLD", "Barycenter position along world y axis",
	&outobj2.myw, H_FLOAT, T_DOUBLE, "%15.8g", "deg"},
  {"X_MAMA", "Barycenter position along MAMA x axis",
	&outobj2.mamaposx, H_FLOAT, T_DOUBLE, "%8.1f", "m**(-6)"},
  {"Y_MAMA", "Barycenter position along MAMA y axis",
	&outobj2.mamaposy, H_FLOAT, T_DOUBLE, "%8.1f", "m**(-6)"},

  {"X_PIXEL", "Object position along x in NDF pixel coordinates",
	&outobj2.ndfposx, H_FLOAT, T_FLOAT, "%10.3f", "pixel"},
  {"Y_PIXEL", "Object position along y in NDF pixel coordinates",
	&outobj2.ndfposy, H_FLOAT, T_FLOAT, "%10.3f", "pixel"},

  {"ALPHA_SKY", "Right ascension of barycenter (native)",
	&outobj2.alphas, H_FLOAT, T_DOUBLE, "%11.7f", "deg"},
  {"DELTA_SKY", "Declination of barycenter (native)",
	&outobj2.deltas, H_FLOAT, T_DOUBLE, "%+11.7f", "deg"},

  {"ALPHA_J2000", "Right ascension of barycenter (J2000)",
	&outobj2.alpha2000, H_FLOAT, T_DOUBLE, "%11.7f", "deg"},
  {"DELTA_J2000", "Declination of barycenter (J2000)",
	&outobj2.delta2000, H_FLOAT, T_DOUBLE, "%+11.7f", "deg"},

  {"ALPHA_B1950", "Right ascension of barycenter (B1950)",
	&outobj2.alpha1950, H_FLOAT, T_DOUBLE, "%11.7f", "deg"},
  {"DELTA_B1950", "Declination of barycenter (B1950)",
	&outobj2.delta1950, H_FLOAT, T_DOUBLE, "%+11.7f", "deg"},

  {"X2_IMAGE", "Variance along x",
	&outobj.mx2, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"Y2_IMAGE", "Variance along y",
	&outobj.my2, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"XY_IMAGE", "Covariance between x and y",
	&outobj.mxy, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"X2_WORLD", "Variance along X-WORLD (alpha)",
	&outobj2.mx2w, H_EXPO, T_DOUBLE, "%15g", "deg**2"},
  {"Y2_WORLD", "Variance along Y-WORLD (delta)",
	&outobj2.my2w, H_EXPO, T_DOUBLE, "%15g", "deg**2"},
  {"XY_WORLD", "Covariance between X-WORLD and Y-WORLD",
	&outobj2.mxyw, H_EXPO, T_DOUBLE, "%15g", "deg**2"},

  {"CXX_IMAGE", "Cxx object ellipse parameter",
	&outobj.cxx, H_EXPO, T_FLOAT, "%12e", "pixel**(-2)"},
  {"CYY_IMAGE", "Cyy object ellipse parameter",
	&outobj.cyy, H_EXPO, T_FLOAT, "%12e", "pixel**(-2)"},
  {"CXY_IMAGE", "Cxy object ellipse parameter",
	&outobj.cxy, H_EXPO, T_FLOAT, "%12e", "pixel**(-2)"},
  {"CXX_WORLD", "Cxx object ellipse parameter (WORLD units)",
	&outobj2.cxxw, H_EXPO, T_FLOAT, "%12e", "deg**(-2)"},
  {"CYY_WORLD", "Cyy object ellipse parameter (WORLD units)",
	&outobj2.cyyw, H_EXPO, T_FLOAT, "%12e", "deg**(-2)"},
  {"CXY_WORLD", "Cxy object ellipse parameter (WORLD units)",
	&outobj2.cxyw, H_EXPO, T_FLOAT, "%12e", "deg**(-2)"},

  {"A_IMAGE", "Profile RMS along major axis",
	&outobj.a, H_FLOAT, T_FLOAT, "%9.3f", "pixel"},
  {"B_IMAGE", "Profile RMS along minor axis",
	&outobj.b, H_FLOAT, T_FLOAT, "%9.3f", "pixel"},
  {"THETA_IMAGE", "Position angle (CCW/x)",
	&outobj.theta, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"A_WORLD", "Profile RMS along major axis (world units)",
	&outobj2.aw, H_FLOAT, T_FLOAT, "%12g", "deg"},
  {"B_WORLD", "Profile RMS along minor axis (world units)",
	&outobj2.bw, H_FLOAT, T_FLOAT, "%12g", "deg"},
  {"THETA_WORLD", "Position angle (CCW/world-x)",
	&outobj2.thetaw, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"THETA_SKY", "Position angle (east of north) (native)",
	&outobj2.thetas, H_FLOAT, T_FLOAT, "%+6.2f", "deg"},
  {"THETA_J2000", "Position angle (east of north) (J2000)",
	&outobj2.theta2000, H_FLOAT, T_FLOAT, "%+6.2f", "deg"},
  {"THETA_B1950", "Position angle (east of north) (B1950)",
	&outobj2.theta1950, H_FLOAT, T_FLOAT, "%+6.2f", "deg"},

  {"ERRX2_IMAGE", "Variance of position along x",
	&outobj.poserr_mx2, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"ERRY2_IMAGE", "Variance of position along y",
	&outobj.poserr_my2, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"ERRXY_IMAGE", "Covariance of position between x and y",
	&outobj.poserr_mxy, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"ERRX2_WORLD", "Variance of position along X-WORLD (alpha)",
	&outobj2.poserr_mx2w, H_EXPO, T_DOUBLE, "%15e", "deg**2"},
  {"ERRY2_WORLD", "Variance of position along Y-WORLD (delta)",
	&outobj2.poserr_my2w, H_EXPO, T_DOUBLE, "%15e", "deg**2"},
  {"ERRXY_WORLD", "Covariance of position X-WORLD/Y-WORLD",
	&outobj2.poserr_mxyw, H_EXPO, T_DOUBLE, "%15e", "deg**2"},

  {"ERRCXX_IMAGE", "Cxx error ellipse parameter",
	&outobj2.poserr_cxx, H_EXPO, T_FLOAT, "%12g", "pixel**(-2)"},
  {"ERRCYY_IMAGE", "Cyy error ellipse parameter",
	&outobj2.poserr_cyy, H_EXPO, T_FLOAT, "%12g", "pixel**(-2)"},
  {"ERRCXY_IMAGE", "Cxy error ellipse parameter",
	&outobj2.poserr_cxy, H_EXPO, T_FLOAT, "%12g", "pixel**(-2)"},
  {"ERRCXX_WORLD", "Cxx error ellipse parameter (WORLD units)",
	&outobj2.poserr_cxxw, H_EXPO, T_FLOAT, "%12g", "deg**(-2)"},
  {"ERRCYY_WORLD", "Cyy error ellipse parameter (WORLD units)",
	&outobj2.poserr_cyyw, H_EXPO, T_FLOAT, "%12g", "deg**(-2)"},
  {"ERRCXY_WORLD", "Cxy error ellipse parameter (WORLD units)",
	&outobj2.poserr_cxyw, H_EXPO, T_FLOAT, "%12g", "deg**(-2)"},

  {"ERRA_IMAGE", "RMS position error along major axis",
	&outobj2.poserr_a, H_FLOAT, T_FLOAT, "%8.4f", "pixel"},
  {"ERRB_IMAGE", "RMS position error along minor axis",
	&outobj2.poserr_b, H_FLOAT, T_FLOAT, "%8.4f", "pixel"},
  {"ERRTHETA_IMAGE", "Error ellipse position angle (CCW/x)",
	&outobj2.poserr_theta, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRA_WORLD", "World RMS position error along major axis",
	&outobj2.poserr_aw, H_FLOAT, T_FLOAT, "%12g", "pixel"},
  {"ERRB_WORLD", "World RMS position error along minor axis",
	&outobj2.poserr_bw, H_FLOAT, T_FLOAT, "%12g", "pixel"},
  {"ERRTHETA_WORLD", "Error ellipse pos. angle (CCW/world-x)",
	&outobj2.poserr_thetaw, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRTHETA_SKY", "Native error ellipse pos. angle (east of north)",
	&outobj2.poserr_thetas, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRTHETA_J2000", "J2000 error ellipse pos. angle (east of north)",
	&outobj2.poserr_theta2000, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRTHETA_B1950", "B1950 error ellipse pos. angle (east of north)",
	&outobj2.poserr_theta1950, H_FLOAT, T_FLOAT, "%5.1f", "deg"},

  {"MU_THRESHOLD", "Detection threshold above background",
	&outobj2.threshmu, H_FLOAT, T_FLOAT, "%8.4f", "mag * arcsec**(-2)"},
  {"MU_MAX", "Peak surface brightness above background",
	&outobj2.maxmu, H_FLOAT, T_FLOAT, "%8.4f", "mag * arcsec**(-2)"},
  {"ISOAREA_WORLD", "Isophotal area above Analysis threshold",
	&outobj2.npixw, H_FLOAT, T_FLOAT, "%12g", "deg**2"},
  {"ISOAREAF_WORLD", "Isophotal area (filtered) above Detection threshold",
	&outobj2.fdnpixw, H_FLOAT, T_FLOAT, "%12g", "deg**2"},
  {"ISO0", "Isophotal area at level 0",
	&outobj.iso[0], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO1", "Isophotal area at level 1",
	&outobj.iso[1], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO2", "Isophotal area at level 2",
	&outobj.iso[2], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO3", "Isophotal area at level 3",
	&outobj.iso[3], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO4", "Isophotal area at level 4",
	&outobj.iso[4], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO5", "Isophotal area at level 5",
	&outobj.iso[5], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO6", "Isophotal area at level 6",
	&outobj.iso[6], H_INT, T_LONG, "%8d", "pixel**2"},
  {"ISO7", "Isophotal area at level 7",
	&outobj.iso[7], H_INT, T_LONG, "%8d", "pixel**2"},

  {"RAD0", "Mean radius at brightness threshold 0",
	&outobj.rad[0], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD1", "Mean radius at brightness threshold 1",
	&outobj.rad[1], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD2", "Mean radius at brightness threshold 2",
	&outobj.rad[2], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD3", "Mean radius at brightness threshold 3",
	&outobj.rad[3], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD4", "Mean radius at brightness threshold 4",
	&outobj.rad[4], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD5", "Mean radius at brightness threshold 5",
	&outobj.rad[5], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD6", "Mean radius at brightness threshold 6",
	&outobj.rad[6], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD7", "Mean radius at brightness threshold 7",
	&outobj.rad[7], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD8", "Mean radius at brightness threshold 8",
	&outobj.rad[8], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD9", "Mean radius at brightness threshold 9",
	&outobj.rad[9], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD10", "Mean radius at brightness threshold 10",
	&outobj.rad[10], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD11", "Mean radius at brightness threshold 11",
	&outobj.rad[11], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD12", "Mean radius at brightness threshold 12",
	&outobj.rad[12], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD13", "Mean radius at brightness threshold 13",
	&outobj.rad[13], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD14", "Mean radius at brightness threshold 14",
	&outobj.rad[14], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},
  {"RAD15", "Mean radius at brightness threshold 15",
	&outobj.rad[15], H_FLOAT, T_FLOAT, "%8.4f", "arcsec"},

  {"FLAGS", "Extraction flags",
	&outobj.flag, H_INT, T_SHORT, "%3d", ""},
  {"IMAFLAGS_ISO", "FLAG-image flags OR'ed over the iso. profile",
	outobj.imaflag, H_INT, T_LONG, "%9u", "",
	1, &prefs.imaflag_size},
  {"NIMAFLAGS_ISO", "Number of flagged pixels entering IMAFLAGS_ISO",
	outobj.imanflag, H_INT, T_LONG, "%9d", "",
	1, &prefs.imanflag_size},

  {"FWHM_IMAGE", "FWHM assuming a gaussian core",
	&outobj.fwhm, H_FLOAT, T_FLOAT, "%8.2f", "pixel"},
  {"FWHM_WORLD", "FWHM assuming a gaussian core",
	&outobj2.fwhmw, H_FLOAT, T_FLOAT, "%12g", "deg"},
  {"ELONGATION", "A_IMAGE/B_IMAGE",
	&outobj2.elong, H_FLOAT, T_FLOAT, "%8.3f", ""},
  {"ELLIPTICITY", "1 - B_IMAGE/A_IMAGE",
	&outobj2.ellip, H_FLOAT, T_FLOAT, "%8.3f", ""},
  {"CLASS_STAR", "S/G classifier output",
	&outobj2.sprob, H_FLOAT, T_FLOAT, "%5.2f", ""},
  {"VIGNET", "Pixel data around detection",
	&outobj2.vignet, H_FLOAT, T_FLOAT, "%12g", "count", 2,
	prefs.vignetsize},
  {"VIGNET_SHIFT", "Pixel data around detection, corrected for shift",
	&outobj2.vigshift, H_FLOAT, T_FLOAT, "%12g", "count", 2,
	prefs.vigshiftsize},
  {"VECTOR_ASSOC", "ASSOCiated parameter vector",
	&outobj2.assoc, H_FLOAT, T_FLOAT, "%12g", "", 1,
	&prefs.assoc_size},
  {"NUMBER_ASSOC", "Number of ASSOCiated IDs",
	&outobj2.assoc_number, H_INT, T_LONG, "%10d", ""},

  {"THRESHOLDMAX", "Maximum threshold possible for detection",
	&outobj.dthresh, H_FLOAT, T_FLOAT, "%12g", "count"},

  {"FLUX_GROWTH", "Cumulated growth-curve",
	&outobj2.flux_growth, H_FLOAT, T_FLOAT, "%12g", "count", 1,
	&prefs.flux_growthsize},
  {"FLUX_GROWTHSTEP", "Step for growth-curves",
	&outobj2.flux_growthstep, H_FLOAT, T_FLOAT, "%10.3f", "pixel"},
  {"MAG_GROWTH", "Cumulated magnitude growth-curve",
	&outobj2.mag_growth, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	1, &prefs.mag_growthsize},
  {"MAG_GROWTHSTEP", "Step for growth-curves",
	&outobj2.mag_growthstep, H_FLOAT, T_FLOAT, "%10.3f", "pixel"},
  {"FLUX_RADIUS", "Fraction-of-light radii",
	&outobj2.flux_radius, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	1, &prefs.flux_radiussize},

  {"XPSF_IMAGE", "X coordinate from PSF-fitting",
	&outobj2.x_psf, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	1, &prefs.psf_xsize},
  {"YPSF_IMAGE", "Y coordinate from PSF-fitting",
	&outobj2.y_psf, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	1, &prefs.psf_ysize},
  {"XPSF_WORLD", "PSF position along world x axis",
	&outobj2.xw_psf, H_FLOAT, T_DOUBLE, "%15.8g", "deg",
	1, &prefs.psf_xwsize},
  {"YPSF_WORLD", "PSF position along world y axis",
	&outobj2.yw_psf, H_FLOAT, T_DOUBLE, "%15.8g", "deg",
	1, &prefs.psf_ywsize},

  {"ALPHAPSF_SKY", "Right ascension of the fitted PSF (native)",
	&outobj2.alphas_psf, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	1, &prefs.psf_alphassize},
  {"DELTAPSF_SKY", "Declination of the fitted PSF (native)",
	&outobj2.deltas_psf, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	1, &prefs.psf_deltassize},

  {"ALPHAPSF_J2000", "Right ascension of the fitted PSF (J2000)",
	&outobj2.alpha2000_psf, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	1, &prefs.psf_alpha2000size},
  {"DELTAPSF_J2000", "Declination of the fitted PSF (J2000)",
	&outobj2.delta2000_psf, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	1, &prefs.psf_delta2000size},

  {"ALPHAPSF_B1950", "Right ascension of the fitted PSF (B1950)",
	&outobj2.alpha1950_psf, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	1, &prefs.psf_alpha1950size},
  {"DELTAPSF_B1950", "Declination of the fitted PSF (B1950)",
	&outobj2.delta1950_psf, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	1, &prefs.psf_delta1950size},

  {"FLUX_PSF", "Flux from PSF-fitting",
	&outobj2.flux_psf, H_FLOAT, T_FLOAT, "%12g", "count",
	1, &prefs.psf_fluxsize},
  {"FLUXERR_PSF", "RMS flux error for PSF-fitting",
	&outobj2.fluxerr_psf, H_FLOAT, T_FLOAT, "%12g", "count",
	1, &prefs.psf_fluxerrsize},
  {"MAG_PSF", "Magnitude from PSF-fitting",
	&outobj2.mag_psf, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	1, &prefs.psf_magsize},
  {"MAGERR_PSF", "RMS magnitude error from PSF-fitting",
	&outobj2.magerr_psf, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	1, &prefs.psf_magsize},

  {"NITER_PSF", "Number of iterations for PSF-fitting",
	&outobj2.niter_psf, H_INT, T_SHORT, "%3d", ""},
  {"CHI2_PSF", "Reduced chi2 from PSF-fitting",
	&outobj2.chi2_psf, H_FLOAT, T_FLOAT, "%9g", ""},

  {"ERRX2PSF_IMAGE", "Variance of PSF position along x",
	&outobj2.poserrmx2_psf, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"ERRY2PSF_IMAGE", "Variance of PSF position along y",
	&outobj2.poserrmy2_psf, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"ERRXYPSF_IMAGE", "Covariance of PSF position between x and y",
	&outobj2.poserrmxy_psf, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"ERRX2PSF_WORLD", "Variance of PSF position along X-WORLD (alpha)",
	&outobj2.poserrmx2w_psf, H_EXPO, T_DOUBLE, "%15e", "deg**2"},
  {"ERRY2PSF_WORLD", "Variance of PSF position along Y-WORLD (delta)",
	&outobj2.poserrmy2w_psf, H_EXPO, T_DOUBLE, "%15e", "deg**2"},
  {"ERRXYPSF_WORLD", "Covariance of PSF position X-WORLD/Y-WORLD",
	&outobj2.poserrmxyw_psf, H_EXPO, T_DOUBLE, "%15e", "deg**2"},

  {"ERRCXXPSF_IMAGE", "Cxx PSF error ellipse parameter",
	&outobj2.poserrcxx_psf, H_EXPO, T_FLOAT, "%12g", "pixel**(-2)"},
  {"ERRCYYPSF_IMAGE", "Cyy PSF error ellipse parameter",
	&outobj2.poserrcyy_psf, H_EXPO, T_FLOAT, "%12g", "pixel**(-2)"},
  {"ERRCXYPSF_IMAGE", "Cxy PSF error ellipse parameter",
	&outobj2.poserrcxy_psf, H_EXPO, T_FLOAT, "%12g", "pixel**(-2)"},
  {"ERRCXXPSF_WORLD", "Cxx PSF error ellipse parameter (WORLD units)",
	&outobj2.poserrcxxw_psf, H_EXPO, T_FLOAT, "%12g", "deg**(-2)"},
  {"ERRCYYPSF_WORLD", "Cyy PSF error ellipse parameter (WORLD units)",
	&outobj2.poserrcyyw_psf, H_EXPO, T_FLOAT, "%12g", "deg**(-2)"},
  {"ERRCXYPSF_WORLD", "Cxy PSF error ellipse parameter (WORLD units)",
	&outobj2.poserrcxyw_psf, H_EXPO, T_FLOAT, "%12g", "deg**(-2)"},

  {"ERRAPSF_IMAGE", "PSF RMS position error along major axis",
	&outobj2.poserra_psf, H_FLOAT, T_FLOAT, "%8.4f", "pixel"},
  {"ERRBPSF_IMAGE", "PSF RMS position error along minor axis",
	&outobj2.poserrb_psf, H_FLOAT, T_FLOAT, "%8.4f", "pixel"},
  {"ERRTHTPSF_IMAGE", "PSF error ellipse position angle (CCW/x)",
	&outobj2.poserrtheta_psf, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRAPSF_WORLD", "World PSF RMS position error along major axis",
	&outobj2.poserraw_psf, H_FLOAT, T_FLOAT, "%12g", "pixel"},
  {"ERRBPSF_WORLD", "World PSF RMS position error along minor axis",
	&outobj2.poserrbw_psf, H_FLOAT, T_FLOAT, "%12g", "pixel"},
  {"ERRTHTPSF_WORLD", "PSF error ellipse pos. angle (CCW/world-x)",
	&outobj2.poserrthetaw_psf, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRTHTPSF_SKY", "Native PSF error ellipse pos. angle (east of north)",
	&outobj2.poserrthetas_psf, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRTHTPSF_J2000", "J2000 PSF error ellipse pos. angle (east of north)",
	&outobj2.poserrtheta2000_psf, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"ERRTHTPSF_B1950", "B1950 PSF error ellipse pos. angle (east of north)",
	&outobj2.poserrtheta1950_psf, H_FLOAT, T_FLOAT, "%5.1f", "deg"},

  {"X2PC_IMAGE", "PC variance along x",
	&outobj2.mx2_pc, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"Y2PC_IMAGE", "PC variance along y",
	&outobj2.my2_pc, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},
  {"XYPC_IMAGE", "PC covariance between x and y",
	&outobj2.mxy_pc, H_EXPO, T_DOUBLE, "%15e", "pixel**2"},

  {"APC_IMAGE", "PC profile RMS along major axis",
	&outobj2.a_pc, H_FLOAT, T_FLOAT, "%8.2f", "pixel"},
  {"BPC_IMAGE", "PC profile RMS along minor axis",
	&outobj2.b_pc, H_FLOAT, T_FLOAT, "%8.2f", "pixel"},
  {"THETAPC_IMAGE", "PC position angle (CCW/x)",
	&outobj2.theta_pc, H_FLOAT, T_FLOAT, "%5.1f", "deg"},
  {"PC", "Principal components",
	&outobj2.vector_pc, H_FLOAT, T_FLOAT, "%15e", "",
	1, &prefs.pc_vectorsize},
/*
	{"RETINOUT", T_FLOAT, &outobj.retinout, "%13g "},
*/
  {""}
  };

