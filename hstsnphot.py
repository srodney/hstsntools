#!/usr/bin/env python
# C.McCully, S.Rodney and B.Patel
#
"""
hstsnphot.py
Perform aperture photometry in Vega magnitudes an image or set of images.
Command-line operation for single-image, single-source photometry:

Syntax :  hstsnphot.py [options] image [xcoord ycoord]

Options: 
-h/--help      : print this help
-v             : verbose mode

For measuring a single point:
--forced       : forced photometry mode (position center is not allowed to float)
--wcs          : coordinates are WCS (decimal deg) instead of pixel coords
--AB           : report AB mags instead of Vega
--decliner     : source is a negative PSF (e.g. from a declining SN)
--upperlim     : report the mag corresponding to a 3-sigma flux upper limit
--noupperlim   : never report 3-sigma upper limits, even when flux is low
--snanadat     : record and report the magnitudes in SNANA-style OBS: lines

For measuring photometry on all Fake SNe in the image:
--fake         : Produces a .pdf showing plots of the recovered fake SN photometry
                 In this case, don't provide any x/y coordinates. 
"""

"""

WFC3/IR PSF Encircled Energy Fraction vs. Aperture Radius (arcsec)
wl(um):   0.7   0.8   0.9   1.0   1.1   1.2   1.3   1.4   1.5   1.6   1.7
 0.10   0.575 0.549 0.524 0.502 0.484 0.468 0.453 0.438 0.426 0.410 0.394
 0.15   0.736 0.714 0.685 0.653 0.623 0.596 0.575 0.558 0.550 0.539 0.531
 0.20   0.802 0.794 0.780 0.762 0.739 0.712 0.683 0.653 0.631 0.608 0.590
 0.25   0.831 0.827 0.821 0.813 0.804 0.792 0.776 0.756 0.735 0.708 0.679
 0.30   0.850 0.845 0.838 0.833 0.828 0.822 0.816 0.808 0.803 0.789 0.770
 0.40   0.878 0.876 0.869 0.859 0.850 0.845 0.841 0.838 0.840 0.836 0.832
 0.50   0.899 0.894 0.889 0.884 0.878 0.868 0.858 0.852 0.854 0.850 0.848
 0.60   0.916 0.913 0.904 0.897 0.893 0.889 0.883 0.875 0.870 0.863 0.859
 0.80   0.937 0.936 0.929 0.924 0.918 0.909 0.903 0.900 0.903 0.900 0.895
 1.00   0.951 0.951 0.946 0.941 0.935 0.930 0.925 0.920 0.917 0.912 0.909
 1.50   0.967 0.969 0.967 0.965 0.963 0.959 0.954 0.951 0.952 0.948 0.943
 2.00   0.974 0.977 0.976 0.975 0.973 0.972 0.969 0.967 0.970 0.967 0.963

WFC3/UVIS PSF Encircled Energy Fraction vs. Aperture Radius (arcsec)
wl(um):   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0   1.1
 0.08                                 0.633 0.587 0.561
 0.10   0.660 0.739 0.754 0.745 0.720 0.687 0.650 0.623 0.612 0.605
 0.12                                 0.741 0.713 0.685
 0.15   0.717 0.793 0.823 0.834 0.832 0.823 0.807 0.778 0.742 0.699
 0.16                                 0.823 0.816 0.791
 0.20   0.752 0.822 0.845 0.859 0.859 0.857 0.853 0.847 0.844 0.829
 0.25   0.781 0.844 0.864 0.875 0.877 0.874 0.870 0.867 0.868 0.864
 0.30   0.802 0.858 0.880 0.888 0.890 0.889 0.883 0.879 0.879 0.876
 0.40   0.831 0.880 0.899 0.911 0.910 0.907 0.906 0.904 0.900 0.894
 0.50   0.861 0.894 0.912 0.923 0.925 0.923 0.918 0.915 0.918 0.917
 0.60   0.884 0.906 0.922 0.932 0.934 0.933 0.931 0.927 0.927 0.923
 0.80   0.936 0.928 0.936 0.944 0.947 0.946 0.945 0.942 0.944 0.942
 1.00   0.967 0.946 0.948 0.954 0.955 0.955 0.955 0.952 0.955 0.952
 1.50   0.989 0.984 0.973 0.970 0.970 0.969 0.967 0.966 0.970 0.968
 2.00   0.994 0.992 0.989 0.985 0.980 0.977 0.976 0.975 0.978 0.976


To update or compute new aperture corrections:
  fobs = ftrue * EEfrac
  mobs = -2.5*log10( fobs ) + ZPT
  mtrue = -2.5*log10( ftrue ) + ZPT
 
  apcor = mobs - mtrue 
        = -2.5 * (  log10( fobs ) - log10( ftrue ) )
        = -2.5 * log10( fobs / ftrue )
  apcor = -2.5 * log10( EEfrac )

NOTE: use the encircled energy fraction relative to the 0.4" aperture EEfrac
"""

import os
import pyfits
from numpy import array, shape, nan_to_num

# A list of apertures we want phot to do:
APLIST = '2,3,4,5'
APSIZE = APLIST.split(',')

# dictionary of AB and Vega zero points for all filters
# AB zero points are included in case we want to use AB magnitudes in the future.
# These are the zero-points for a 0.4 arcsec aperture, from  http://www.stsci.edu/hst/wfc3/phot_zp_lbn/
# collected c.2012,June.  Note that the aperture corrections given below are computed relative to 
# these 0.4" aperture zeropoints, not to the infinite aperture zeropoints
ZPT_ACS_VEGA = {'F435W': 25.76695,
                'F475W': 26.16252,
                'F555W': 25.72747,
                'F606W': 26.40598,
                'F625W': 25.74339,
                'F775W': 25.27728,
                'F814W': 25.51994,
                'F850LP':24.32305 }

ZPT_ACS_AB = {'F435W': 25.65777,
              'F475W': 26.05923,
              'F555W': 25.7184,
              'F606W': 26.49113,
              'F625W': 25.90667,
              'F775W': 25.66504,
              'F814W': 25.94333,
              'F850LP':24.84245 }

ZPT_WFC3_IR_AB = {'F105W':26.0974,
                  'F110W':26.6424,
                  'F125W':26.0449,
                  'F140W':26.2608,
                  'F160W':25.7551,
                  'F098M':25.5041,
                  'F127M':24.4545,
                  'F139M':24.2880,
                  'F153M':24.2725 }

ZPT_WFC3_IR_VEGA = {'F105W':25.4523,
                    'F110W':25.8829,
                    'F125W':25.1439,
                    'F140W':25.1845,
                    'F160W':24.5037,
                    'F098M':24.9424,
                    'F127M':23.4932,
                    'F139M':23.2093,
                    'F153M':23.0188 }

ZPT_WFC3_UVIS_AB = {'F218W':22.7776,
                    'F225W':23.8629,
                    'F275W':23.9740,
                    'F336W':24.5377,
                    'F350LP':26.8413,
                    'F390W':25.2389,
                    'F438W':24.7097,
                    'F475W':25.5755,
                    'F555W':25.6890,
                    'F600LP':25.7681,
                    'F606W':25.9668,
                    'F625W':25.4230,
                    'F775W':24.7487,
                    'F814W':24.9912,
                    'F850LP':23.7234,
                    'F763M':24.2070,
                    'F845M':23.7811,
                    }

ZPT_WFC3_UVIS_VEGA = {'F218W':21.0878,
                      'F225W':22.2034,
                      'F275W':22.4757,
                      'F336W':23.3531,
                      'F350LP':26.6852,
                      'F390W':25.0240,
                      'F438W':24.8629,
                      'F475W':25.6739,
                      'F555W':25.7144,
                      'F600LP':25.4440,
                      'F606W':25.8843,
                      'F625W':25.2750,
                      'F775W':24.3679,
                      'F814W':24.5730,
                      'F850LP':23.20260,
                      'F763M':23.8283,
                      'F845M':23.2809,
                      }

# Here are the infinite aperture zeropoints, used by sextract.py
# for extended source photometry without aperture corrections
ZPT_WFC3_IR_AB_INF = {'F105W':26.2687,'F110W':26.8223,'F125W':26.2303,
                      'F140W':26.4524,'F160W':25.9463,'F098M':25.6674,
                      'F127M':24.6412,'F139M':24.4793,'F153M':24.4635 } 

ZPT_WFC3_IR_VEGA_INF = { 'F105W':25.6236,'F110W':26.0628,'F125W':25.3293,
                         'F140W':25.3761,'F160W':24.6949,'F098M':25.1057,
                         'F127M':23.6799,'F139M':23.4006,'F153M':23.2098 }

ZPT_ACS_WFC_AB_INF = { 'F850LP':24.750 }


# WFC3 IR and WFC3 UVIS F350LP Corrections were calculated by B. Patel.
# Patel measured several stars with different apertures. He calculated the correction to both
# 0.4" and to a nominally infinite aperture of 28" using the correct zeropoints respectively.
# The corrected magnitudes match well so we have adopted the aperture corrections to 0.4". The
# aperture uncertainty corresponds to the standard deviation in the recovered magnitudes for the stars.
# ACS aperture corrections are taken from Table 3 of Siranni et al. 2005.:

#Added by C.McCully
#Fractional Aperture Correction from Sirianni for ACS
#ACOR_FRAC = {
#    'F435W': [0.626, 0.775, 0.832, 0.859],
#    'F475W': [0.669, 0.801, 0.850, 0.876],
#    'F555W': [0.662, 0.795, 0.843, 0.871],
#    'F606W': [0.656, 0.796, 0.845, 0.871],
#    'F625W': [0.648, 0.795, 0.844, 0.868],
#    'F775W': [0.623, 0.780, 0.840, 0.863],
#    'F814W': [0.603, 0.764, 0.833, 0.858],
#    'F850LP':[0.538, 0.689, 0.779, 0.815],
#}
#Fractional Aperture Uncertainty
#If the uncertainty was reported as 0.000, I assume a 0.001 uncertainty.
#AERR_FRAC = {
#    'F435W': [0.011, 0.004, 0.002, 0.001],
#    'F475W': [0.006, 0.003, 0.002, 0.002],
#    'F555W': [0.003, 0.002, 0.001, 0.001],
#    'F606W': [0.005, 0.003, 0.001, 0.001],
#    'F625W': [0.004, 0.002, 0.001, 0.001],
#    'F775W': [0.005, 0.003, 0.001, 0.001],
#    'F814W': [0.006, 0.004, 0.002, 0.001],
#    'F850LP':[0.001, 0.003, 0.002, 0.001],
#}

APCOR = {
    # WFC3- IR 
    # aperture size = 0.18, 0.27, 0.36, 0.45 arcsec
    'F160W': array([0.486, 0.163, 0.026, -0.018]), #H
    'F140W': array([0.435, 0.135, 0.021, -0.017]), #Y
    'F125W': array([0.412, 0.112, 0.018, -0.018]), #J
    'F110W': array([0.387, 0.105, 0.018, -0.019]), #N
    'F105W': array([0.384, 0.105, 0.020, -0.022]), #M
    'F098M': array([0.347, 0.093, 0.019, -0.022]), #L 
    'F127M': array([0.432, 0.120, 0.018, -0.017]), #O
    'F139M': array([0.454, 0.134, 0.020, -0.016]), #P
    'F153M': array([0.485, 0.163, 0.025, -0.018]), #Q

    # ACS-WFC
    # aperture size = 0.10, 0.15, 0.20, 0.25 arcsec
    'F850LP':array([0.673, 0.404, 0.271,  0.222]), #Z
    'F775W': array([0.514, 0.270, 0.189,  0.160]), #X
    'F606W': array([0.458, 0.248, 0.183,  0.150]), #V
    'F625W': array([0.471, 0.249, 0.184,  0.154]), #R
    'F435W': array([0.509, 0.277, 0.200,  0.165]), #B
    'F475W': array([0.436, 0.241, 0.176,  0.144]), #G
    'F814W': array([0.549, 0.292, 0.198,  0.166]), #I

    # WFC3-UVIS
    # aperture size = 0.08, 0.12, 0.16, 0.20 arcsec
    'F350LP':array([0.931, 0.423, 0.193,  0.104]), #W
    'F555W': array([0.448, 0.249, 0.185,  0.150]), #
    'F225W': array([0.490, 0.219, 0.113,  0.066]), #S
    'F275W': array([0.342, 0.325, 0.308,  0.200]), #T
    'F336W': array([0.467, 0.208, 0.116,  0.071]), #U
    'F390W': array([0.384, 0.178, 0.110,  0.070]), #C
    'F763M': array([0.431, 0.240, 0.111,  0.064]), #7 
    'F845M': array([0.494, 0.280, 0.128,  0.068]), #9 
    }

APERR = {
    'F350LP':array([0.187, 0.128, 0.073, 0.034]),
    'F160W': array([0.031, 0.011, 0.005, 0.004]),
    'F140W': array([0.034, 0.011, 0.002, 0.002]),
    'F125W': array([0.032, 0.010, 0.005, 0.005]),
    'F105W': array([0.034, 0.012, 0.004, 0.004]),
    'F110W': array([0.018, 0.006, 0.004, 0.004]),
    'F098M': array([0.089, 0.027, 0.005, 0.005]), #L
    'F127M': array([0.071, 0.023, 0.005, 0.003]), #O
    'F139M': array([0.075, 0.024, 0.004, 0.004]), #P
    'F153M': array([0.035, 0.010, 0.002, 0.001]), #Q
    'F850LP':array([0.002, 0.005, 0.003, 0.001]),
    'F775W': array([0.009, 0.004, 0.001, 0.001]),
    'F606W': array([0.008, 0.004, 0.001, 0.001]),
    'F625W': array([0.007, 0.003, 0.001, 0.001]),
    'F435W': array([0.019, 0.006, 0.003, 0.001]),
    'F475W': array([0.010, 0.004, 0.003, 0.002]),
    'F814W': array([0.011, 0.006, 0.003, 0.001]),
    'F555W': array([0.005, 0.003, 0.001, 0.001]),
    'F225W': array([0.097, 0.077, 0.065, 0.059]), #S
    'F275W': array([0.1, 0.1, 0.05, 0.02]),       #T
    'F336W': array([0.047, 0.013, 0.007, 0.006]), #U
    'F390W': array([0.022, 0.008, 0.005, 0.005]), #C
    'F763M': array([0.05, 0.015, 0.01, 0.008]), #7 CRUDE ESTIMATES!
    'F845M': array([0.05, 0.015, 0.01, 0.008]), #9 CRUDE ESTIMATES!
    }

APDEFAULTSIZE = { 
    # Default aperture size, in pixels.
    # Basically we choose the smallest aperture that will have a small enough
    # uncertainty to not dominate the statistical uncertainty
    # This aperture will be reported in the final output file, but we will save the
    # results from the other aperture file in the full photometry file.
    'ACS':4.0,      # = 0.20" (native=0.05"/pix, we drizzle to 0.05"/pix)
    'WFC3_IR':3.0,  # = 0.27" (native=0.13"/pix, we drizzle to 0.09"/pix)
    'WFC3_UVIS':5.0 # = 0.20" (native=0.04"/pix, we drizzle to 0.04"/pix)
    }

# A dictionary with the size of the pixels in arcseconds for a given instrument
PIXSIZE = {
           'ACS':0.05, # native=0.05, we drizzle to 0.05"/pixel
           'WFC3_IR':0.09, # native=0.13, we drizzle to 0.09"/pixel
           'WFC3_UVIS':0.04 # native=0.04, we drizzle to 0.04"/pixel
           }
NATIVE_PIXSIZE = {
           'ACS':0.0495, # native=0.05, we drizzle to 0.05"/pixel
           'WFC3_IR':0.1282, # native=0.13, we drizzle to 0.09"/pixel
           'WFC3_UVIS':0.0396 # native=0.04, we drizzle to 0.04"/pixel
           }
# Conversion table from full filter names to single-letter abbreviations
FilterAlpha = { 'F225W':'S','F275W':'T','F336W':'U','F390W':'C',
                'F350LP':'W',
                'F435W':'B','F475W':'G','F606W':'V','F625W':'R',
                'F775W':'X','F814W':'I','F850LP':'Z',
                'F125W':'J','F160W':'H','F125W+F160W':'A',
                'F105W':'Y','F110W':'M','F140W':'N',
                'F098M':'L','F127M':'O','F139M':'P','F153M':'Q',
                'G141':'4','G102':'2','blank':'0',
                'F763M':'7','F845M':'9',
                }
SNbandlist =   ['H','J','W']
ACSbandlist =  ['B','G','V','R','X','I','Z']
IRbandlist =   ['H','J','Y','M','N']
UVISbandlist = ['S','T','U','W','C','7','9']
WFC3bandlist = IRbandlist + UVISbandlist
cambandlist = { 'a':ACSbandlist, 'u':UVISbandlist, 'i':IRbandlist }

#From Krist 2003 ISR ACS 2003-06 we find the fwhm to be about 2.5 +- 0.5.
#For WFC3 IR we use table 7.5 from the handbook to get a FWHM at 1.1 + 0.1
#For WFC UVIS we use table 6.7 from the handbook. Note that we have multiplied the WFC3 values
# by 1.5 and are therefore slightly larger than the mean because
# our undersampling with only 2 dither positions will broaden the fitted FWHM slightly.
psffwhm = {'ACS':2.5,'WFC3_IR':2.35,'WFC3_UVIS':2.7}

def filter2band( filter ):
    """ convert a full ACS or WFC3 filter string into its
    single-digit alphabetic code.  e.g:
      filter2band('F125W')  ==>  'J'
    """
    import exceptions
    if filter in FilterAlpha.keys():
        return( FilterAlpha[filter] )
    else :
        raise exceptions.RuntimeError(
            "Unknown filter %s"%filter )

def filter2alpha( filter ):
    return( filter2band( filter ) )

def band2filter( alpha ):
    """ convert a single-digit alphabetic bandpass code
    into a full ACS or WFC3 filter string. e.g:
      band2filter('J') ==> 'F125W'
    """
    import exceptions
    alpha = alpha.upper()
    for filter in FilterAlpha.keys():
        if FilterAlpha[filter] == alpha :
            return( filter )
    else :
        raise exceptions.RuntimeError(
            "Unknown filter %s"%alpha )

def alpha2filter( alpha ) :
    return( band2filter(alpha) )

def run(imfile, coordlist, WCS=False, outroot=None, poisson=True, verbose=True, debug=False, abmags=False, decliner=False, smallskyann=False,
        calgorithm='gauss', upperlim=None, snanadat=False, cbox = 5.0 ):
    """
    Use iraf.digiphot.apphot to collect aperture photometry in Vega magnitudes.

    Required Arguments:
        imfile (str):  fits image file name (cannot be fpacked)
        coord : the coordinate list for photometry
          as a str, this is the name of a text file in two columns giving x,y coordinates.
          as a python list or numpy array, this gives coordinate pairs for each target
             e.g.  coord=[ [1024,440], [502,680] ]

    Optional arguments:
        WCS : If input coordinates are in RA and DEC  (degrees), user must set WCS=True.
        outroot :  root of the output file names (outroot.mag, outroot.full, outroot.out)
            if unspecified, the root of  imfile  is used.
        poisson: Use a poisson noise model or not. If set to false, a constant noise model is
                used. Default is True
        calgorithm: Centering algorithm to pass to iraf.apphot. Choices are "gauss","centroid",
                    "none", or "ofilter". Default is "gauss"
        cbox: Size of centering box to pass to iraf.apphot. Default is 5.0
        upperlim:  None = when flux<3-sigma then report 3-sig upper limit
                   True = force computation of 3-sigma upper limit
                   False = disallow upper limits (i.e. always report measured flux)
        snanadat:  True = Report the mjd, flux and mag in a SNANA-style OBS: line
                   False = report the filename, mjd, source position, mag and errors.
        verbose: Default is True.
        debug: Start pdb. Default is False.

    Output products :  phot.out, phot.mag, phot.full
        phot.out : the raw output from apphot
        phot.full :  a detailed photometry file
        phot.mag :  header+one line summary output file :
              #image name, filter, xpos, ypos, magnitude, and errors

    Requires :  numpy, pyfits  (astLib is required if WCS=True)
    """
    if debug:
        import pdb; pdb.set_trace()
    from math import sqrt
    from numpy import nan, log10
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)

    if outroot == None : outroot = os.path.splitext(os.path.basename(imfile))[0]

    # Get info from header
    hdr = pyfits.getheader(imfile)
    exptime = hdr['EXPTIME']
    if 'EXPSTART' in hdr:
        mjdobs =  "%.1f"%hdr['EXPSTART']
    else : 
        mjdobs = '0.0'

    if 'DATE-OBS' in hdr : 
        dateobs = hdr['DATE-OBS']
    elif 'DATE' in hdr : 
        dateobs = hdr['DATE']
    else : 
        dateobs = 'unknown'

    if 'D001OUUN' in hdr : 
        if hdr['D001OUUN'].upper() == 'CPS' : 
            gain = exptime
            exptime=1.0
        elif hdr['D001OUUN'].upper() == 'COUNTS' : 
            gain = 1.0
    elif hdr['BUNIT'].upper() == 'ELECTRONS/S':
        gain = exptime
        exptime=1.0
    elif hdr['BUNIT'].upper() == 'ELECTRONS' :
        gain = 1.0
    
    else:
        print 'Problem determining units, check image'
        return
    

    if 'FILTER' in hdr:
        filt = hdr['FILTER']
    elif 'FILTER1' in hdr :
        filt = hdr['FILTER1']

    if filt.startswith('CLEAR'):
        filt = hdr['FILTER2']

    instrument = hdr['INSTRUME']
    if instrument == 'WFC3':
        instrument = instrument + '_' + hdr['DETECTOR']

    if abmags : 
        if instrument == 'WFC3_IR' :
            ZPT = ZPT_WFC3_IR_AB[filt]
        elif instrument == 'WFC3_UVIS' :
            ZPT = ZPT_WFC3_UVIS_AB[filt]
        elif instrument == 'ACS' :
            ZPT = ZPT_ACS_AB[filt]
        else :
            print "Can't handle instrument: %s; check the fits header." % instrument
            #Give up and die
            return
    else : 
        if instrument == 'WFC3_IR' :
            ZPT = ZPT_WFC3_IR_VEGA[filt]
        elif instrument == 'WFC3_UVIS' :
            ZPT = ZPT_WFC3_UVIS_VEGA[filt]
        elif instrument == 'ACS' :
            ZPT = ZPT_ACS_VEGA[filt]
        else :
            print "Can't handle instrument: %s; check the fits header." % instrument
            #Give up and die
            return


    # format for coord list:
    #    coordlist= [ [x0,y0], [x1,y1], ... ]
    #    or a two-column text file
    if type(coordlist) == str :
        # user provided name of a coordinate file
        fin = open(coordlist, 'r')
        coordlines = fin.readlines()
        fin.close()
        coordvals = array([ cline.split() for cline in coordlines ], dtype=float)
    elif len(shape(coordlist)) == 1 :
        # user provided something like coord=[x,y]
        coordvals = array([ coordlist ])
    else :
        # user provided something like coord= [[x1,y1],[x2,y2],[x3,y3]]
        coordvals = array(coordlist)

    # how many objects do we have ?
    numcoo = len(coordvals)

    if WCS:
        # If coords in wcs instead of x,y,
        # get wcs information from imfile for converting to xy.
        from astLib import astWCS
        wcsfits = astWCS.WCS(imfile)

    # (re)write a list of x,y positions
    coxyfile = '%s.xycoo' % outroot
    coxy = open(coxyfile, 'w')
    for coord in coordvals :
        if WCS :
            #Convert from RA and Dec to xy 
            # NOTE: wcsfits returns values based on a 0,0 origin, but the iraf phot 
            # packages expect a 1,1 origin. So we add 1 to each value
            xy = wcsfits.wcs2pix(coord[0], coord[1]) 
            xy[0] += 1 
            xy[1] += 1 
        else :
            xy = coord
        print >> coxy, "%10.2f  %10.2f" % (float(xy[0]), float(xy[1]))
    coxy.close()
    if verbose>1: print("XY coords written to %s" % coxyfile)

    """ iraf.digiphot.apphot.datapars :
    2013.09.06 SR: updated to use Poisson noise model (which includes sky noise) as 
    the default, instead of 'constant' """
    iraf.unlearn(iraf.apphot.phot)
    iraf.unlearn(iraf.datapars)
    iraf.datapars.scale = 1.0
    iraf.datapars.fwhmpsf = psffwhm[instrument]
    iraf.datapars.emission = not decliner
    iraf.datapars.sigma = 'INDEF'
    iraf.datapars.datamin = 'INDEF'
    iraf.datapars.datamax = 'INDEF'

    if poisson:
        iraf.datapars.noise = 'poisson'
    else:
        iraf.datapars.noise = 'constant'

    iraf.datapars.ccdread = ''
    #iraf.datapars.gain = ''
    iraf.datapars.readnoise = 0.0
    #iraf.datapars.exposure = ' '
    #iraf.datapars.airmass = ''
    #iraf.datapars.obstime = ''
    iraf.datapars.itime = exptime
    iraf.datapars.epadu = gain
    iraf.datapars.xairmass = 'INDEF'
    iraf.datapars.ifilter = 'INDEF'
    iraf.datapars.otime = 'INDEF'

    # iraf.digiphot.apphot.centerpars :
    iraf.unlearn(iraf.centerpars)

    iraf.centerpars.calgorithm = calgorithm
    iraf.centerpars.cbox = cbox

    iraf.centerpars.cthreshold = 0.0
    iraf.centerpars.minsnratio = 1.0
    iraf.centerpars.cmaxiter = 10.0
    iraf.centerpars.maxshift = 1.0
    iraf.centerpars.clean = False
    iraf.centerpars.rclean = 1.0
    iraf.centerpars.rclip = 2.0
    iraf.centerpars.kclean = 3.0
    iraf.centerpars.mkcenter = False

    # iraf.digiphot.apphot.fitskypars :
    iraf.unlearn(iraf.fitskypars)
    iraf.fitskypars.salgorithm = 'median'
    if smallskyann : 
        iraf.fitskypars.annulus = 8.0
        iraf.fitskypars.dannulus = 12.0
    else : 
        iraf.fitskypars.annulus = 25.0
        iraf.fitskypars.dannulus = 40.0
    iraf.fitskypars.skyvalue = 0.0
    iraf.fitskypars.smaxiter = 10.0
    iraf.fitskypars.sloclip = 0.0
    iraf.fitskypars.shiclip = 0.0
    iraf.fitskypars.snreject = 50.0
    iraf.fitskypars.sloreject = 3.0
    iraf.fitskypars.shireject = 3.0
    iraf.fitskypars.khist = 3.0
    iraf.fitskypars.binsize = 0.1
    iraf.fitskypars.smooth = False
    iraf.fitskypars.rgrow = 0.0
    iraf.fitskypars.mksky = False

    # iraf.digiphot.apphot.photpars :
    iraf.unlearn(iraf.photpars)
    iraf.photpars.weighting = 'constant'
    iraf.photpars.apertures = APLIST
    iraf.photpars.zmag = ZPT
    iraf.photpars.mkapert = False

    photparams = {
        'interac':False,
        'radplot':False,
        }

    magfile_out = outroot + '.out'
    magfile_full = outroot + '.full'
    magfile_phot = outroot + '.mag'

    if os.path.exists(magfile_out) :
        os.remove(magfile_out)

    # run photometry using the newly created coxyfile for providing input coordinates
    try:
        iraf.phot(image=imfile, skyfile='', coords=coxyfile, output=magfile_out,
                  verify=False, verbose=True, Stdout=1, **photparams)
    except iraf.IrafError, e:
        print("phot failed on %s with IRAF error :\n%s"%(imfile,e))

    if verbose > 1: print("Output from apphot written to %s" % magfile_out)
    f = open(magfile_out, 'r')
    maglines = f.readlines()
    f.close()

    full = open(magfile_full, 'w')
    magfile = open(magfile_phot, 'w')

    if snanadat : 
        magfile.write("VARLIST:  MJD  FLT FIELD  FLUXCAL  FLUXCALERR   MAG    MAGERR\n")
    else : 
        magfile.write("# Image MJD Filter XPos(Pix) YPos(Pix) Mag StatErr SysErr TotalErr \n")

    #get the photometric data from the apphot output
    for i in range(numcoo) :
        #x and y center
        xpos = maglines[75 + i * (4 + len(APSIZE)) + 1].split()[0]
        ypos = maglines[75 + i * (4 + len(APSIZE)) + 1].split()[1]

        #sky info
        msky = maglines[75 + i * (4 + len(APSIZE)) + 2].split()[0]
        skystdev = maglines[75 + i * (4 + len(APSIZE)) + 2].split()[1]
        nskypix = maglines[75 + i * (4 + len(APSIZE)) + 2].split()[3]

        # epadu : this is set to equal the exposure time when working with JHU products in cps.
        #   so that the daophot error calculation appropriately converts the measured cps flux 
        #   into counts when computing the flux and mag uncertainties.
        epadu = float(maglines[21].split()[3])

        if WCS :
            # Calculate RA,DEC pos of object
            radec = wcsfits.pix2wcs(float(xpos), float(ypos))
        else :
            radec = [0, 0]

        full.write("# Image Date MJD Instrument Filter Zeropoint X(Pix) Y(Pix) RA Dec Sky(Cts) SkyStd(Cts) \n")
        full.write(imfile + " " + dateobs + " " + mjdobs + " " + instrument + " " + filt + " " + "%.3f"%ZPT + " " + xpos + " " + ypos + " " + str(radec[0]) + " "
                   + str(radec[1]) + " " + msky + " " + skystdev + "\n")

        full.write("\n")

        #full.write("# Radius(Pix) Radius(Arcseconds) ApertureCorrection ApertureError Flux Mag MagError TotalError\n")
        full.write('# Rap(Pix) Rap(") ApCorr ApEr      Flux   FluxErrTot   Mag   MagErr MagErrTot\n' )

        #Get the flux and magnitude info
        for k in range(len(APSIZE)) :
            aper_rad = float(APSIZE[k])
            # calculate the aperture radius in arcsec
            aper_arc = aper_rad * PIXSIZE[instrument]

            ### need something here to handle INDEFs, should they arise, LGS 12/03/12
            # INDEF magnitudes now get recorded as NaNs as suggested by LGS (Implemented by CVM,SR)
            mag = maglines[75 + i * (4 + len(APSIZE)) + 4 + k].split()[4] 
            if mag=='INDEF': mag=nan
            else: mag = float(mag) - float(APCOR[filt][k])

            mag_err = maglines[75 + i * (4 + len(APSIZE)) + 4 + k].split()[5]
            if mag_err =='INDEF':
                mag_err=nan
                tot_err=nan
            else:
                mag_err = float(mag_err)
                tot_err = sqrt(mag_err**2 + APERR[filt][k]**2)

            # collect flux and compute the flux error
            area = float(maglines[75 + i * (4 + len(APSIZE)) + 4 + k].split()[2])
            try : flux = float(maglines[75 + i * (4 + len(APSIZE)) + 4 + k].split()[3])
            except ValueError : flux = 0.0
            try : flux_err = sqrt (flux / float(epadu) + float(area) * float(skystdev)**2 + float(area)**2 * float(skystdev)**2 / float(nskypix))
            except ValueError : flux_err = 0.0

            # apply the aperture correction to the flux and flux_err
            if flux>0: 
                flux *= 10**(0.4*APCOR[filt][k])
                flux_aperr = 1.0857 * flux * APERR[filt][k]
                flux_err *= 10**(0.4*APCOR[filt][k])
                flux_err = sqrt( flux_err**2 + flux_aperr**2 )
            else : 
                flux_err *= 10**(0.4*APCOR[filt][k])
                flux_aperr = 1.0857 * flux_err * APERR[filt][k]
                flux_err = sqrt( flux_err**2 + flux_aperr**2 )
                
            # Check if the source is below detectable flux levels,
            # and if so, report a 3-sigma upper limit on the flux
            if upperlim or (flux < 3*flux_err and upperlim!=0): 
                if verbose>2: print("Reporting 3-sigma flux limit for %.1f pix aperture."%aper_rad)
                flux = 0.0
                flux_err = 3 * flux_err
                mag = -2.5*log10( flux_err ) + ZPT + APCOR[filt][k]
                mag_err = -9
                tot_err = -9

            elif flux < 3*flux_err and upperlim==False and verbose>3: 
                print("WARNING: %.1f pix flux is below 3-sigma flux limit, but not reporting upper limits."%aper_rad)
                
            # write it out to the .full file
            full.write('  %3.1f      %5.3f  %6.3f %5.3f  %10.5f %9.5f   %6.3f %5.3f %5.3f\n'%(
                    aper_rad, aper_arc, APCOR[filt][k], APERR[filt][k], flux, flux_err, mag, mag_err, tot_err ) )

            if aper_rad == APDEFAULTSIZE[instrument]:
                # write lines for phot.mag output
                if snanadat : 
                    # Print a SNANA-style OBS line to the screen, e.g. 
                    # OBS: 56456.500  H  wol    0.000    8.630   25.160   -9.000
                    fluxcal = flux * 10**(0.4*(27.5-ZPT)) 
                    fluxcalerr = flux_err * 10**(0.4*(27.5-ZPT))
                    magline = "OBS: %8.2f   %s   %s %8.3f %8.3f    %8.3f %8.3f\nx,y = %s, %s     ra,dec= %.6f,%6f\n"%(
                        float(mjdobs), FilterAlpha[filt], imfile[:3], fluxcal, fluxcalerr, mag, tot_err, xpos, ypos, radec[0], radec[1] )
                else : 
                    magline = imfile +" "+mjdobs+" "+filt + " " + xpos + " " + ypos + \
                        " %.3f %.3f"%(mag,mag_err) + " %.3f %.3f\n"%(APERR[filt][k],tot_err)
                magfile.write( magline )

    full.close()
    magfile.close()
    if verbose > 2 :
        print("Full photometry output written to " + magfile_full)
        print("Summary photometry at default aperture size written to " + magfile_phot)

    if verbose : 
        # print mag summary to screen
        fin = open(magfile_phot, 'r')
        maglines = fin.readlines()
        fin.close()
        for mline in maglines :
            print(mline.strip())

    if snanadat : magdat=magline
    else : magdat = rdmagfile( magfile_phot )
    return( magdat )


def rdmagfile( magfile ):
    """ read data from a .mag summary photometry file into ndarrays."""
    from numpy import loadtxt

    magdat = loadtxt( magfile, dtype=str )

    if len(magdat.shape) == 1 :
        magdat = array( [ magdat ] )
    imfile = magdat[:,0]
    mjd = magdat[:,1].astype(float)
    filt = magdat[:,2]
    x = magdat[:,3].astype(float)
    y = magdat[:,4].astype(float)
    mag = magdat[:,5].astype(float)
    magerr = magdat[:,-1].astype(float)
    
    return( {'imfile':imfile,'mjd':mjd, 'filter':filt,'x':x,'y':y, 'mag':mag, 'magerr':magerr})


def rdfullfile( fullfile ):
    """ read data from a .full summary photometry file into ndarrays."""
    from numpy import loadtxt

    # TODO : make this more general, and robust against .full format changes
    # Read the header data line
    fin = open( fullfile, 'r' )
    hdr = fin.readline()
    hdrdat = fin.readline().strip().split()
    fin.close()
    imfile = hdrdat[0]
    obsdate = hdrdat[1]
    mjdobs = float(hdrdat[2])
    instrument = hdrdat[3]
    filt = hdrdat[4]
    zpt = float(hdrdat[5])
    x = float(hdrdat[6])
    y = float(hdrdat[7])
    ra = float(hdrdat[8])
    dec = float(hdrdat[9])
    skycts = float(hdrdat[10])
    skystdcts = float(hdrdat[11])

    # read the aperture photometry data grid
    fulldat = loadtxt( fullfile, skiprows=3 )
    rpix = fulldat[:,0]
    rarcsec = fulldat[:,1]
    apcor = fulldat[:,2]
    aperr = fulldat[:,3]
    flux = fulldat[:,4]
    fluxerr = fulldat[:,5]
    mag = fulldat[:,6]
    magerr= fulldat[:,7]
    toterr = fulldat[:,8]
    
    return( {'imfile':imfile,'mjd':mjdobs, 'instrument':instrument, 
             'filter':filt,'zpt':zpt,'x':x,'y':y, 'ra':ra,'dec':dec,
             'skycts':skycts,'skystdcts':skystdcts, 'rpix':rpix, 'rarcsec':rarcsec, 'apcor':apcor,
             'aperr':aperr,'flux':flux,'fluxerr':fluxerr,'mag':mag,'magerr':magerr,'toterr':toterr } )


def mksnanadat( snanadatfile=None,
                survey='HST', nickname='Supernova', snid='HST00hst',
                photversion='HST_CANDELS1',
                ra=0.0, dec=0.0, mwebv=0.008,
                specz=-9.0, speczerr=-9.0,
                photz=-9.0, photzerr=-9.0,
                debug=False ):
    """ read in data from all the .mag files in this
    directory, convert to fluxes and print out
    a SNANA .dat file with the full light curve.
    Header is populated with metadata info specified in the arguments,
    which default to some generic CANDELS values.
    """
    import glob
    from numpy import unique, argmin, argmax, isnan, abs, log10
    import util

    if debug : import pdb;  pdb.set_trace()

    if specz>0 : z, zerr = specz, speczerr
    elif photz>0: z, zerr = photz, photzerr
    else : z, zerr = specz, speczerr

    # Read in each .full file, convert it to the SNANA "OBS:" line format.
    fullfilelist = glob.glob( "*.full")
    Nobs = len(fullfilelist)
    bandlist, maglist, mjdlist, obslinelist, signoiselist = [], [], [], [], []
    for fullfile in fullfilelist :
        fulldat = rdfullfile( fullfile )
        mjd = fulldat['mjd']
        filt = fulldat['filter']
        band = FilterAlpha[filt] 
        field = os.path.basename(fulldat['imfile'])[:3]
        instrument = fulldat['instrument']
        ibestap = fulldat['rpix'].tolist().index( APDEFAULTSIZE[instrument] )

        # rescale flux to SNANA FLUXCAL units (FLUXCAL zpt=27.5)
        fluxcal = fulldat['flux'][ibestap] * 10**(0.4*(27.5-fulldat['zpt']))
        fluxcalerr = fulldat['fluxerr'][ibestap] * 10**(0.4*(27.5-fulldat['zpt']))
        mag = fulldat['mag'][ibestap]
        magerr = fulldat['toterr'][ibestap]

        if isnan(mag) or isnan(magerr): 
            # Magnitude is undefined (negative flux) 
            # so we report a 3-sigma upper limit, indicated by setting the magerr to -9
            fluxcal = 0.0 
            fluxcalerr = 3*fluxcalerr
            mag = -2.5*log10( fluxcalerr )  + 27.5
            magerr = -9.0

        # write out the SNANA 7-variable OBS line
        obsline = 'OBS: %9.3f %s %6s  %8.3f %8.3f  %8.3f %8.3f\n'%(
                    mjd, band, field, fluxcal, fluxcalerr, mag, magerr )
        obslinelist.append( obsline )
        bandlist.append( band ) 
        maglist.append( mag )
        mjdlist.append( mjd )
        signoiselist.append( fluxcal / abs(fluxcalerr) )

    # find the date of the highest S/N observation (gives SNANA a guess at the time of peak)
    mjdpk = mjdlist[ argmax( signoiselist ) ]
    # make a string listing all observation bands
    filtstr = ''.join( unique( bandlist )  )

    header = """
SURVEY: %s
NICKNAME: %s
SNID: %s
IAUC:    NULL
PHOTOMETRY_VERSION: %s
SNTYPE: -9
FILTERS: %s
RA:      %.6f  deg
DECL:    %.7f  deg
MAGTYPE: LOG10
MAGREF:  VEGA
MWEBV:   %.4f    MW E(B-V)
REDSHIFT_FINAL:   %.7f +- %.7f (CMB)
HOST_GALAXY_SPEC-Z: %.7f +- %.7f (CMB)
HOST_GALAXY_PHOTO-Z: %.7f +- %.7f (CMB)
REDSHIFT_STATUS: OK
SEARCH_PEAKMJD: %.3f

NOBS: %i
NVAR: 7
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR    MAG     MAGERR
    """% ( survey, nickname, snid, photversion, filtstr, ra, dec, mwebv,
           z, zerr, specz, speczerr, photz, photzerr, mjdpk, Nobs )

    footer = """
END:
"""

    # make a default output file name if needed
    if snanadatfile == None :
        snanadatfile = magfilelist[0].split('_')[0] + '_snana.dat'

    # write out the data
    fout = open(snanadatfile, 'w')
    print >> fout, header
    fout.writelines( obslinelist )
    print >> fout, footer
    fout.close()
    return( snanadatfile )


def fakephot( fakefilename, outfilename=None ):
    """
    collect photometry for all fakes in the image,
    make plots of input vs output positions and mags
    """
    import os
    fakeroot = os.path.splitext(os.path.basename( fakefilename ))[0]
    from pylab import rcParams, rcdefaults
    from math import sqrt

    rcParams['lines.linewidth']=1.5
    rcParams['font.size']=15    # text big enough for a half-page fig
    rcParams['axes.labelsize'] = 15      # fontsize of the x and y labels

    
    if outfilename==None : 
        outfilename = fakefilename.replace('.fits', '_fakephot.pdf')

    # construct coord and mag lists from the fits header keywords
    coordlist,xin,yin,magin = [],[],[],[]
    hdr = pyfits.getheader( fakefilename )
    for key in hdr.keys() : 
        if (key.startswith('FAKE') and key.endswith('X') ): xin.append(hdr[key])
        if (key.startswith('FAKE') and key.endswith('Y') ): yin.append(hdr[key])
        if (key.startswith('FAKE') and key.endswith('M') ): magin.append(hdr[key])

    for i in range(len(xin)):
        coordlist.append((xin[i],yin[i]))

    # collect forced photometry at fixed positions (for checking mags)
    forced = run( fakefilename, coordlist, WCS=False, outroot=fakeroot+'_forced', calgorithm='none',poisson=True, smallskyann=False )

    # collect photometry with recentering allowed (for checking positions)
    free = run( fakefilename, coordlist, WCS=False, outroot=fakeroot+'_free',poisson=True, smallskyann=False )

    # compute residuals
    dx = free['x']-xin
    dy = free['y']-yin
    dmag = forced['mag'] - magin
    dmag = nan_to_num( dmag )

    # make some plots
    from matplotlib.pyplot import figure, clf, savefig, autoscale
    # autoscale( enable=True, axis='both', tight=False )
    fig = figure(1, figsize=[8.5,11])
    ax1 = fig.add_subplot(3,2,1)
    ax1.plot( xin, dx, 'ro', ls=' ' )
    ax1.set_ylabel('$x_{out}-x_{in}$')
    ax1.set_xlabel('$x_{in}$')
    ax1.axhline(0, color='0.3', lw=0.7 )
    
    ax2 = fig.add_subplot(3,2,2)
    ax2.plot( yin, dy, 'ro', ls=' ' )
    ax2.set_ylabel('$y_{out}-y_{in}$')
    ax2.set_xlabel('$y_{in}$')
    ax2.axhline(0, color='0.3', lw=0.7 )

    ax3 = fig.add_subplot(3,2,3)
    ax3.plot( magin, forced['mag'], 'bs', ls=' ' )
    ax3.plot( [ max(magin), min(magin)], [ max(magin), min(magin)], 'k-' )
    ax3.errorbar( magin, forced['mag'], forced['magerr'], marker='s', color='b', ls=' ' )
    ax3.set_ylabel('$m_{out}$')
    ax3.set_xlabel('$m_{in}$')

    #bins,edges = histogram( dmag, bins=max(10, len(dmag)/10) )
    ax4 = fig.add_subplot(3,2,4)
    histout = ax4.hist( dmag, bins=max(5,sqrt(len(magin))), align='mid',histtype='stepfilled' )
    #ax4.plot( edges[:-1], bins, 'b-', drawstyle='steps-mid', ls='-' )
    ax4.set_ylabel('count')
    ax4.set_xlabel('$m_{out}-m_{in}$')
    ax4.set_ylim( tuple(array(ax4.get_ylim())*1.2) )

    ax5 = fig.add_subplot(3,2,5)
    #ax5.plot( xin, dmag, 'bs',ls=' ')
    ax5.errorbar( xin, dmag, forced['magerr'], color='b',marker='s',ls=' ')
    ax5.set_ylabel('$m_{out}-m_{in}$')
    ax5.set_xlabel('$x_{in}$')
    ax5.axhline(0, color='0.3', lw=0.7 )

    ax6 = fig.add_subplot(3,2,6)
    # ax6.plot( yin, dmag, 'bs',ls=' ')
    ax6.errorbar( yin, dmag, forced['magerr'], color='b',marker='s',ls=' ')
    ax6.set_ylabel('$m_{out}-m_{in}$')
    ax6.set_xlabel('$y_{in}$')
    ax6.axhline(0, color='0.3', lw=0.7 )

    fig.suptitle("Fake SN Photometry for %s"%fakefilename)
    fig.subplots_adjust( left=0.12, right=0.95, bottom=0.08, top=0.95, wspace=0.32, hspace=0.2 )

    savefig(outfilename)

    rcdefaults()
    return( dmag )


def get_fake_centroid_and_fluxcorr(filename,x,y,instrument,filt):
    """
    Locate the center of a fake psf and determine the flux correction
    needed to bring it into alignment with the hstsnphot zero points.

    INPUTS: The fake-SN psf image in filename, the expected x,y position
    of the center of the psf, the instrument and filter being modeled.
    RETURNS: xcentroid, ycentroid, fluxcorr
    """
    from pyraf import iraf

    if instrument == 'WFC3_IR' :
        ZPT = ZPT_WFC3_IR_VEGA[filt]
    elif instrument == 'WFC3_UVIS' :
        ZPT = ZPT_WFC3_UVIS_VEGA[filt]
    elif instrument.startswith('ACS') :
        ZPT = ZPT_ACS_VEGA[filt]

    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.unlearn(iraf.apphot.phot)
    iraf.unlearn(iraf.datapars)
    iraf.unlearn(iraf.centerpars)
    #Use the centroid algorithm right now as it seems more robust to geometric distortion.
    iraf.centerpars.calgorithm = 'centroid'
    iraf.centerpars.cbox = 5.0

    iraf.unlearn(iraf.fitskypars)
    iraf.unlearn(iraf.photpars)
    photparams = {
        'interac':False,
        'radplot':False,
        }
    iraf.datapars.readnoise = 0.0
    iraf.datapars.itime = 1.0
    iraf.datapars.epadu = 1.0

    # iraf.digiphot.apphot.fitskypars :
    iraf.unlearn(iraf.fitskypars)
    iraf.fitskypars.salgorithm = 'constant'
    iraf.fitskypars.skyvalue = 0.0

    # iraf.digiphot.apphot.photpars :
    iraf.unlearn(iraf.photpars)
    iraf.photpars.weighting = 'constant'
    iraf.photpars.apertures = APDEFAULTSIZE[instrument] * PIXSIZE[instrument] / NATIVE_PIXSIZE[instrument]
    iraf.photpars.zmag = ZPT
    iraf.photpars.mkapert = False

    #Write the coordinate file starting as position x and y
    coxyfile = 'centroid.xycoo'
    coxy = open(coxyfile, 'w')
    print >> coxy, "%10.2f  %10.2f" % (x,y)
    coxy.close()
    if os.path.exists('centroid.mag'): os.remove('centroid.mag')
    iraf.phot(image=filename, skyfile='', coords=coxyfile, output='centroid.mag',
              verify=False, verbose=True, Stdout=1, **photparams)
    f = open('centroid.mag', 'r')
    maglines = f.readlines()
    f.close()
    xcentroid = float(maglines[76].split()[0])
    ycentroid = float(maglines[76].split()[1])

    #Calculate the aperture correction necessary to match the observed aperture correction.
    meas_flux_frac = float(maglines[79].split()[1])

    expected_flux_frac = 10.0**(-0.4*APCOR[filt][APSIZE.index( '%1i' % APDEFAULTSIZE[instrument])])

    flux_corr = expected_flux_frac / meas_flux_frac
    return xcentroid,ycentroid,flux_corr



if __name__ == '__main__':
    import sys
    import getopt

    verbose=1
    forced = False
    wcs = False
    abmags=False
    dofake=False
    decliner=False
    upperlim=None
    snanadat=False
    smallskyann=False
    debug=False

    # read in arguments and options
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","quiet","help","forced","fake","wcs",'AB','ab','noupperlim','upperlim','snanadat','decliner','debug','smallskyann' ] )
    except getopt.GetoptError: 
        print "Error : incorrect option or missing argument."
        print __doc__
        sys.exit(1)
    for o, a in opt:
        if o in ["-h", "--help"]:
            print __doc__
            sys.exit(0)
        elif o == "-v" :
            verbose = 2
        elif o == "--quiet" :
            verbose = False
        elif o == "--verbose" :
            try: verbose = int(a)
            except ValueError : 
                print("Error : must specify verbosity level as an integer 1-10.")
                sys.exit(0)

        elif o == "--forced" :
            forced=True
        elif o == "--upperlim" :
            upperlim=True
        elif o == "--noupperlim" :
            upperlim=False
        elif o == "--fake" :
            dofake=True
        elif o == "--wcs" :
            wcs=True
        elif o == "--decliner" :
            decliner=True
        elif o == "--snanadat" :
            snanadat=True
        elif o in ["--AB",'--ab'] :
            abmags=True
        elif o == "--smallskyann" :
            smallskyann=True
        elif o == "--debug" :
            debug=True

    if dofake : 
        if len(arg) != 1 :
            print __doc__
            print("\n You provided too many or not enough arguments. Try again")
            sys.exit(0)
        filename=arg[0]
        if verbose : print("Running hstsnphot.fakephot on %s"%filename)
        fakephot( filename )
        if verbose : print("hstsnphot.fakephot finished")
        sys.exit(0)

    # parse the arguments : filename xcoord ycoord
    if len(arg) != 3 :
        print __doc__
        print("\n You provided too many or not enough arguments. Try again")
        sys.exit(0)
    
    filename = arg[0]
    xcoord = float(arg[1])
    ycoord = float(arg[2])
    coordlist = [xcoord, ycoord]
    if verbose > 1 : print( "    %s %.2f %.2f"%( filename, xcoord, ycoord ) )

    if forced: 
        # collect forced photometry at a fixed position
        run( filename, coordlist, WCS=wcs, calgorithm='none', cbox=0, poisson=True, abmags=abmags, upperlim=upperlim, decliner=decliner, snanadat=snanadat, verbose=verbose, debug=debug, smallskyann=smallskyann )
    else : 
        # collect photometry with recentering allowed (for checking positions)
        run( filename, coordlist, WCS=wcs, abmags=abmags, upperlim=upperlim, decliner=decliner, snanadat=snanadat, verbose=verbose, debug=debug, smallskyann=smallskyann )

def mag2fluxcal( mag, magerr=0 ):
    """ convert magnitudes into SNANA-style FLUXCAL units
    (fixed zero point of 27.5 for all bands) """
    from numpy import iterable, abs, array, zeros, any
    
    if not iterable( mag ) : 
        mag = array( [ mag ] )
        magerr = array( [ magerr ] )
    if not iterable( magerr ) : 
        magerr = zeros( len(mag) ) 

    fluxcal, fluxcalerr = [],[]
    for m,me in zip( mag, magerr) : 
        if me < 0 : 
            fluxcal.append( 0 ) 
            fluxcalerr.append( 10**(-0.4*(m-27.5)) )
        else : 
            fluxcal.append( 10**(-0.4*(m-27.5)) )
            fluxcalerr.append( 0.92103 * me * fluxcal[-1] )
    fluxcal = array( fluxcal )
    fluxcalerr = array( fluxcalerr )

    if len(mag)==1 : 
        fluxcal = fluxcal[0]
        fluxcalerr = fluxcalerr[0] 

    if any( magerr ) : return( fluxcal, fluxcalerr ) 
    else : return( fluxcal )
