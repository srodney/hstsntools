import exceptions
import sys
import os
import numpy as np
import glob
import pyfits
from matplotlib import pyplot as pl, cm
from matplotlib import rcParams
from matplotlib.patches import Circle
from hstsntools import imageops
from hstphot import hstphot
from PythonPhot import cntrd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
from astropy.table import Column, vstack

rcParams['text.usetex'] = True

_thisfile = sys.argv[0]
if 'ipython' in _thisfile:
    _thisfile = __file__
_thisdir = os.path.abspath(os.path.dirname(_thisfile))

def convert_deg_to_hmsdms(ra, dec):
    """ Convert from decimal degrees to sexagesimal coordinates.
    :param ra: R.A. in decimal degrees
    :param dec: Decl. in decimal degrees
    :return: two strings giving H:M:S and D:M:S
    """
    coord = SkyCoord( ra, dec, frame="icrs", unit=[u.deg,u.deg] )
    raout  = coord.ra.to_string( unit=u.hour, decimal=False, pad=True,
                                 sep=':', precision=3 )
    decout  = coord.dec.to_string( unit=u.degree, decimal=False, pad=True,
                                   alwayssign=True, sep=':', precision=2 )
    return raout, decout

def convert_hmsdms_to_deg(ra, dec):
    """ Convert from sexagesimal coordinates to decimal degrees.
    :param ra: R.A. in H:M:S
    :param dec: Decl. in D:M:S
    :return: two strings giving decimal degrees
    """
    coord = SkyCoord( ra, dec, frame="icrs", unit=[u.hour,u.deg] )
    raout  = coord.ra.to_string( unit=u.deg, decimal=True, pad=True,
                                 precision=6 )
    decout  = coord.dec.to_string( unit=u.deg, decimal=True, pad=True,
                                   alwayssign=False, precision=6 )
    return float(raout[0]), float(decout[0])

#def reformat_datfile(datfile='candels_sn_hostphot.txt'):


def get_all_positions_and_hostphot(indatfile, outdatfile,
                                   verbose=True, clobber=False):

    indat = ascii.read(indatfile, format='commented_header',
                       header_start=-1, data_start=0,
                       )
    nicknamelist = np.array([nick.lower() for nick in indat['nickname']])
    filterlist = ['f606w','f775w','f850l','f105w','f125w','f140w','f160w']
    for filtername in filterlist:
        mcol = Column( data=np.zeros(len(indat)),
                       name='sb_%s'%filtername.lower(),
                       dtype=float, format='%6.3f')
        merrcol = Column( data=np.zeros(len(indat)),
                          name='sberr_%s'%filtername.lower(),
                          dtype=float, format='%6.3f')
        indat.add_columns([mcol,merrcol])
    indat['z'].format = '%6.3f'
    indat['dz'].format = '%5.3f'
    indat['name'].format = '%13s'
    indat['grade'].format = '%3s'
    indat['decliner'].format = '%3i'

    for nickname in nicknamelist:
        isn = np.where(nicknamelist==nickname)[0]
        raHMS, decDMS, raDD, decDD = update_sn_coordinates(
            nickname, verbose=verbose, clobber=clobber)
        indat['RA_SN'][isn] = raHMS
        indat['DEC_SN'][isn] = decDMS
        indat.write(outdatfile, format='ascii.commented_header')
        hostphotout, sbdict, sberrdict = measure_host_photometry(
            nickname, datfile=outdatfile, verbose=True)
        for filtername in filterlist:
            if filtername in sbdict.keys():
                indat['sb_%s'%filtername][isn] = sbdict[filtername]
                indat['sberr_%s'%filtername][isn] = sberrdict[filtername]
            else :
                indat['sb_%s'%filtername][isn] = -9
                indat['sberr_%s'%filtername][isn] = -9
        indat.write(outdatfile, format='ascii.commented_header')
    print "all done"

def measure_host_photometry(nickname, datfile='candels_sn_hostphot.txt',
                            imsizearcsec=2.5, aparcsec=0.2, verbose=True):
    """ measure the photometry of the host galaxy at the site of the SN
    :param nickname:
    :param verbose:
    :return:  host galaxy surface brightness in a 0.2" radius aperture
    """
    indat = ascii.read(datfile, format='commented_header',
                       header_start=-1, data_start=0)
    nicknamelist = np.array([nick.lower() for nick in indat['nickname']])
    isn = np.where(nicknamelist==nickname)[0]
    raHMS = indat['RA_SN'][isn]
    decDMS = indat['DEC_SN'][isn]

    cwd = os.path.abspath('.')
    sndir = os.path.join(cwd, nickname.lower())
    if not os.path.isdir(sndir):
        raise exceptions.RuntimeError("No directory %s" % sndir)
    hostimlist = glob.glob(sndir+"/%s*e00_reg_dr?_sci.fits" % nickname.lower())

    if verbose :
        fig = pl.figure(2, figsize=[12, 4])
        fig.subplots_adjust(left=0.03, bottom=0.03, right=0.97, top=0.97)
        fig.clf()

    raDD, decDD = convert_hmsdms_to_deg(raHMS, decDMS)
    iax = 0
    sbdict = {}
    sberrdict = {}
    hstphotoutstr = ''
    for hostim in hostimlist:
        imname = os.path.basename(hostim)
        imnameparts = imname.split('_')
        filtername = imnameparts[1]
        if filtername.startswith('~'):
            continue

        # Measure the surface brightness at the SN position
        x,y = hstphot.radec2xy(hostim, raDD, decDD)
        hstphotout = hstphot.dophot(hostim, x, y, aparcsec=aparcsec,
                                    system='AB', printstyle='short',
                                    recenter=False)
        mag = float(hstphotout[0].split()[5])
        magerr = float(hstphotout[0].split()[6])
        surfacebrightness = mag + 2.5*np.log10(np.pi*aparcsec**2)
        hstphotoutstr += hstphotout[0] +'\n'
        sbdict[filtername] = surfacebrightness
        sberrdict[filtername] = magerr
        iax += 1

        # show the host galaxy image with the SN position highlighted
        if verbose:
            ax = fig.add_subplot(1, len(hostimlist)+1, iax)
            show_sn_on_image(nickname, imname, datfile=datfile,
                             imsizearcsec=imsizearcsec, circleradarcsec=0.2,
                             color='red')
            epoch = imnameparts[2]
            ax.text(0.0, 1.2, '%i: %s \n%s' % (iax-1, filtername, epoch),
                    transform=ax.transAxes, ha='left', va='bottom')
        continue
    if verbose:
        print hstphotoutstr
    ax = fig.add_subplot(1, len(hostimlist)+1, len(hostimlist)+1)
    combofile = '%s_composite_sub_masked.fits' % nickname.lower()
    show_sn_on_image(nickname, combofile, datfile=datfile,
                     imsizearcsec=imsizearcsec, circleradarcsec=0.2)
    fig.suptitle('Host Images for %s' % nickname)
    pl.draw()
    raw_input("Showing %s.   Hit return to continue." % nickname)
    return hstphotoutstr, sbdict, sberrdict

def show_sn_on_image(nickname, imfilename,
                     datfile='candels_sn_hostphot.txt',
                     imsizearcsec=2.5, circleradarcsec=0.2,
                     color='red'):
    indat = ascii.read(datfile, format='commented_header',
                       header_start=-1, data_start=0)
    nicknamelist = np.array([nick.lower() for nick in indat['nickname']])
    isn = np.where(nicknamelist==nickname)[0]
    raHMS = indat['RA_SN'][isn]
    decDMS = indat['DEC_SN'][isn]
    raDD, decDD = convert_hmsdms_to_deg(raHMS, decDMS)

    cwd = os.path.abspath('.')
    sndir = os.path.join(cwd, nickname.lower())
    ax = pl.gca()
    combofile = os.path.join(sndir, imfilename )
    combodat = pyfits.getdata(combofile)
    xsn, ysn = hstphot.radec2xy(combofile, raDD, decDD)
    pixscale = imageops.getpixscale(combofile)
    imsizepix = imsizearcsec / pixscale
    halfimsizepix = round(imsizepix/2)
    imdatcenter = combodat[
                  int(ysn-halfimsizepix):int(ysn+halfimsizepix),
                  int(xsn-halfimsizepix):int(xsn+halfimsizepix)]

    ax.imshow(imdatcenter, aspect='equal', interpolation='nearest',
              vmin=-0.05, vmax=0.15, cmap=cm.Greys_r)
    circleradpix = circleradarcsec / pixscale
    c = Circle((halfimsizepix - 1 + (xsn % 1), halfimsizepix - 1 + (ysn % 1)),
                radius=circleradpix, edgecolor=color, facecolor='None')
    ax.add_patch(c)
    ax.set_xticks([])
    ax.set_yticks([])
    return

def update_sn_coordinates(nickname, datfile='candels_sn_hostphot.txt',
                          verbose=True, clobber=False):
    """ display IR diff image, let the user select the epochs where the SN is
    visible, make a stack, get an updated measurement of the SN position
    :return:
    """
    cwd = os.path.abspath('.')
    sndir = os.path.join(cwd, nickname.lower())
    if not os.path.isdir(sndir):
        raise exceptions.RuntimeError("No directory %s" % sndir)
    diffimlist = glob.glob(sndir+"/*f1*sub_masked.fits")
    combofile = os.path.join(sndir, '%s_composite_sub_masked.fits' %
                             nickname.lower())
    indat = ascii.read(datfile, format='commented_header',
                       header_start=-1, data_start=0)
    nicknamelist = np.array([nick.lower() for nick in indat['nickname']])
    isn = np.where(nicknamelist==nickname)[0]
    raHMS = indat['RA_SN'][isn]
    decDMS = indat['DEC_SN'][isn]
    raDD, decDD = convert_hmsdms_to_deg(raHMS, decDMS)

    if not os.path.isfile(combofile) or clobber:
        fig = pl.figure(1, figsize=[12, 4])
        fig.subplots_adjust(left=0.03, bottom=0.03, right=0.97, top=0.97)
        if verbose:
            print("Constructing composite image %s from: %s" % (
                combofile,
                str([os.path.basename(diffim) for diffim in combinelist])))
        fig.clf()
        fig.suptitle('%s Recentering' % nickname)
        iax = 0
        if verbose:
            print "plotting diff images to select for composite image"
        for diffim in diffimlist:
            imname = os.path.basename(diffim)
            imnameparts = imname.split('_')
            filter = imnameparts[1]
            epoch = imnameparts[2]
            iax += 1

            ax = fig.add_subplot(1, len(diffimlist), iax)
            show_sn_on_image(nickname, diffim, datfile=datfile,
                             imsizearcsec=1.0, circleradarcsec=0.2)
            ax.text(0.0, 1.2, '%i: %s \n%s' % (iax, filter, epoch),
                    transform=ax.transAxes, ha='left', va='bottom')
        pl.draw()
        iwithsn = input("Recentering for %s\n" % nickname +
                        "Enter a comma-sep'd list of image numbers in which"
                        " the SN is bright enough for centroiding:"
                        "\n   ")
        combinelist = np.array(diffimlist)[list(iwithsn)]
        combofile = imageops.imaverage(combinelist, combofile, clobber=True)
    elif os.path.isfile(combofile) and not clobber:
        print "%s exists.  Not clobbering." % combofile

    fig = pl.figure(10, figsize=[4, 4])
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)
    show_sn_on_image(nickname, combofile, datfile=datfile,
                     imsizearcsec=1.0, circleradarcsec=0.2, color='cyan')

    # locate the sn centroid position
    xsn, ysn = hstphot.radec2xy(combofile, raDD, decDD)
    combodat = pyfits.getdata(combofile)
    xnew, ynew = cntrd.cntrd(combodat, xsn, ysn, 5)
    xnew += 1 # convert from iraf to python-style coordinates
    ynew += 1 # convert from iraf to python-style coordinates
    ranew, decnew = hstphot.xy2radec(combofile, xnew, ynew)
    ranew = ranew[0]
    decnew = decnew[0]
    raHMSnew, decDMSnew =  convert_deg_to_hmsdms(ranew, decnew)
    oldstr = '%-15s old(cyan): %11s %-11s  %10.6f %10.6f %9.3f %9.3f'%(
        nickname.lower(), raHMS[0], decDMS[0], raDD, decDD, xsn, ysn)
    newstr = '%-15s new(red) : %11s %-11s  %10.6f %10.6f %9.3f %9.3f'%(
        nickname.lower(), raHMSnew, decDMSnew, ranew, decnew, xnew, ynew)
    print oldstr
    print newstr

    # show the new position as a red circle
    pixscale = imageops.getpixscale(combofile)
    circleradpix = 0.2 / pixscale
    imsizepix = 1.0 / pixscale
    halfimsizepix = round(imsizepix/2)
    c = Circle((halfimsizepix - 1 + int(xnew)-int(xsn) + (xnew % 1),
                halfimsizepix - 1 + int(ynew)-int(ysn) + (ynew % 1)),
                radius=circleradpix, edgecolor='red', facecolor='None')
    ax.add_patch(c)
    fig.suptitle('%s : red=new position' % nickname)
    return raHMSnew, decDMSnew, ranew, decnew

from pytools import plotsetup, colorpalette as cp
def mk_histogram_figure(sndatfile='candels_sn_hostphot2.txt',
                        galdata=None,
                        galdatfilelist=['CANDELS.GOODSS.F160W.v1.PHOTOZ.CAT',
                                        'CANDELS.GOODSN.F160W.v1.PHOTOZ.CAT',
                                        'CANDELS.UDS.F160W.v1.photoz.cat',
                                        ],
                        insetstamps=False):
    fig = plotsetup.fullpaperfig()
    cwd = os.path.abspath('.')
    sndata = ascii.read(sndatfile, format='commented_header',
                          header_start=-1, data_start=0)
    igot160 = np.where(sndata['sb_f160w']>0)[0]
    sb160 = sndata['sb_f160w'][igot160]
    pIa = sndata['PIa'][igot160]

    pl.clf()
    ax = pl.axes([0.06,0.12,0.9,0.6])

    ax.hist( sb160, bins=np.arange(20,28,0.5), weights=pIa,
             color=cp.darkred, alpha=0.5, normed=True )
    ax.hist( sb160, bins=np.arange(20,28,0.5), weights=(1-pIa),
             color=cp.darkblue, alpha=0.5, normed=True )

    if galdata is None:
        galdatalist = [ascii.read(galdatfile) for galdatfile in galdatfilelist]
        galdata = vstack(galdatalist, join_type='outer')

    m160 = galdata['mag_f160w']
    zphot = galdata['Photo_z']
    iz12 = np.where((zphot>1) & (zphot<2.5) )[0]
    m160 = m160[iz12]

    ax.hist( m160, bins=np.arange(20,28,0.2), color=cp.black,
             alpha=0.3, normed=True, zorder=-1000 )
    fig.subplots_adjust(left=0.06, bottom=0.12, right=0.97, top=0.97)
    ax.set_xlabel('SN Host Galaxy Surface Brightness (mag arcsec$^{-2}$)')
    ax.set_ylabel('Normalized Count')
    ax.set_yticks([])

    ax.text(0.7, 0.92, '\\noindent Magnitude from total flux,\n'
                      'all galaxies with $1<z<2.5$',
            transform=ax.transAxes, ha='left', va='top', color=cp.darkgrey)
    ax.text(0.3, 0.92, '\\noindent CANDELS SN Host galaxies,\n'
                      'with $1<z<2.5$',
            transform=ax.transAxes, ha='left', va='top', color=cp.darkgrey)
    ax.text(0.15, 0.52, 'Type Ia SN',
            transform=ax.transAxes, ha='right', va='top', color=cp.darkred)
    ax.text(0.35, 0.62, 'CC SN',
            transform=ax.transAxes, ha='left', va='top', color=cp.darkblue)
    ax.plot([22.25, 22.75], [0.2, 0.3], marker=' ', ls='-', color=cp.darkblue,
            lw=1)
    ax.plot([21.75, 20.75], [0.13, 0.24], marker=' ', ls='-', color=cp.darkred,
            lw=1)

    pl.draw()

    if insetstamps:
        ax1 = pl.axes([0.1,0.76,0.18,0.2])
        show_sn_on_image('camille','camille_f160w_e00_reg_drz_sci.fits')

        ax2 = pl.axes([0.42,0.76,0.18,0.2])
        show_sn_on_image('gore','gore_f160w_e00_reg_drz_sci.fits')

        ax3 = pl.axes([0.72,0.76,0.18,0.2])
        show_sn_on_image('fairbanks','fairbanks_f160w_e00_reg_drz_sci.fits')

        #ax1.set_xticks([])
        #ax2.set_xticks([])
        #ax3.set_xticks([])
        #ax1.set_yticks([])
        #ax2.set_yticks([])
        #ax3.set_yticks([])

        imdat = pyfits.getdata('camille/camille_f160w_e00_reg_drz_sci.fits')
        h,w = imdat.shape
        #ax1.imshow( imdat[h/2-15:h/2+15,w/2-15:w/2+15], vmin=-0.05, vmax=0.15,
        #            cmap=cm.Greys_r, aspect='equal', interpolation='nearest')

    return galdata

