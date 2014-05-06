# 2014.04.29 
# S.Rodney
# HST Filter transmission curves: plotting and such

import numpy as np
from matplotlib import pylab as pl
import os

topdir = os.path.abspath( '.' )
try :
    sndataroot = os.environ['SNDATA_ROOT']
    os.chdir( sndataroot+'/filters/HST')

    w435, f435 = np.loadtxt( 'ACS_WFC_F435W.dat', unpack=True )
    w606, f606 = np.loadtxt( 'ACS_WFC_F606W.dat', unpack=True )
    w625, f625 = np.loadtxt( 'ACS_WFC_F625W.dat', unpack=True )
    w814, f814 = np.loadtxt( 'ACS_WFC_F814W.dat', unpack=True )

    w350, f350 = np.loadtxt( 'WFC3_UVIS_F350LP.dat', unpack=True )
    w606u, f606u = np.loadtxt( 'WFC3_UVIS_F606W.dat', unpack=True )
    w127, f127 = np.loadtxt( 'WFC3_IR_F127M.dat', unpack=True )
    w125, f125 = np.loadtxt( 'WFC3_IR_F125W.dat', unpack=True )
    w160, f160 = np.loadtxt( 'WFC3_IR_F160W.dat', unpack=True )
    w153, f153 = np.loadtxt( 'WFC3_IR_F153M.dat', unpack=True )
    w139, f139 = np.loadtxt( 'WFC3_IR_F139M.dat', unpack=True )
    w140, f140 = np.loadtxt( 'WFC3_IR_F140W.dat', unpack=True )

    os.chdir( sndataroot+'/filters/Bessell90')
    wB, fB = np.loadtxt( 'Bessell90_B.dat', unpack=True )
    wV, fV = np.loadtxt( 'Bessell90_V.dat', unpack=True )
    wR, fR = np.loadtxt( 'Bessell90_R.dat', unpack=True )
    wI, fI = np.loadtxt( 'Bessell90_I.dat', unpack=True ) 

except KeyError : 
    pass
finally : 
    os.chdir(topdir)

def filtername2datfile( filtername, camera=None):
    """ Given an abbreviated filter name, returns the name of the .dat file
    containing the transmission curve.
    """
    fname = filtername.upper()
    if fname.startswith('F1') : return( 'WFC3_IR_%s.dat'%fname )
    elif 'UV' in camera.upper():
        return( 'WFC3_UVIS_%s.dat'%fname )
    elif 'ACS' in camera.upper():
        return( 'ACS_WFC_%s.dat'%fname )
    elif fname=='F350LP' :
        return( 'WFC3_UVIS_%s.dat'%fname )
    else :
        print("Must specify a camera for filter %s."%fname)
        return(None)


def computeScaling( filt1, filt2, camera1=None, camera2=None ) :
    """determine the flux scaling factor that should be multiplied to
    filt1 to match the throughput of filt2.  This returns just a
    single number, effectively assuming the source SED is flat across
    the bandpass, so that we just need to correct for total
    throughput, not for the shape of the filter.
    """
    from scipy import integrate as scint

    if filt1.lower().startswith('f') :
        filt1 = filtername2datfile( filt1, camera=camera1 )
    if filt2.lower().startswith('f') :
        filt2 = filtername2datfile( filt2, camera=camera2 )
    if not filt1.endswith('.dat') or not filt2.endswith('.dat') :
        print("Must specify a filter name (e.g. F160W) or a .dat file.")
        return( None )

    # read in the transmission curves for filters 1 and 2 
    topdir = os.path.abspath( '.' )
    sndataroot = os.environ['SNDATA_ROOT']
    os.chdir( sndataroot+'/filters/HST')
    w1, f1 = np.loadtxt( filt1, unpack=True )
    w2, f2 = np.loadtxt( filt2, unpack=True )
    os.chdir( topdir )

    # integrate 
    int1 = scint.simps( f1, w1 )
    int2 = scint.simps( f2, w2 )
   
    # divide
    return( int2 / int1 )


def computeScaling2to1( filt1, filt2, filt3,
                        camera1=None, camera2=None, camera3=None) :
    """Determine the flux scaling factor for matching the sum of filt1+filt2
    to filt3.  This returns the value that should be multiplied to
    (filt1+filt2) to match the throughput of filt3.  This returns just a
    single number, effectively assuming the source SED is flat across
    the bandpass, so that we just need to correct for total
    throughput, not for the shape of the filter.
    """
    from scipy import integrate as scint

    if filt1.lower().startswith('f') :
        filt1 = filtername2datfile( filt1, camera=camera1 )
    if filt2.lower().startswith('f') :
        filt2 = filtername2datfile( filt2, camera=camera2 )
    if filt3.lower().startswith('f') :
        filt3 = filtername2datfile( filt3, camera=camera3 )
    if not (filt1.endswith('.dat') and filt2.endswith('.dat')
            and filt3.endswith('.dat') ):
        print("Must specify a filter name (e.g. F160W) or a .dat file.")
        return( None )

    # read in the transmission curves for filters
    topdir = os.path.abspath( '.' )
    sndataroot = os.environ['SNDATA_ROOT']
    os.chdir( sndataroot+'/filters/HST')
    w1, f1 = np.loadtxt( filt1, unpack=True )
    w2, f2 = np.loadtxt( filt2, unpack=True )
    w3, f3 = np.loadtxt( filt3, unpack=True )
    os.chdir( topdir )

    # integrate
    int1 = scint.simps( f1, w1 )
    int2 = scint.simps( f2, w2 )
    int3 = scint.simps( f3, w3 )

    # sum and divide
    return( int3 / (int1+int2) )



def plotmedbands( z = 2, day=5 ):
    from hstsntools import snana
    w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.dat', day=day )

    w1az = w1a * (1+z)
    f1az = f1a / f1a.max() / 2.
    clf()

    ax1 = subplot(3,1,1)
    plot(w125, f125, 'b--', label='F125W')
    plot(w127, f127, 'b-', label='F127M')
    plot(w1az, f1az, 'r-', label='_nolegend_')
    ax1.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
    ax1.set_xlim( 9000, 20000 )
    ax1.text(9500,0.2, 'SNIa\nz=%.1f\nt=%i'%(z,day), color='r',ha='left',va='bottom')
    setp(ax1.get_xticklabels(), visible=False)
    setp(ax1.get_yticklabels(), visible=False)


    ax2 = subplot(3,1,2, sharex=ax1, sharey=ax1)
    plot(w140, f140, 'g--',label='F140W')
    plot(w139, f139, 'g-',label='F139M')
    plot(w1az, f1az, 'r-', label='_nolegend_')
    ax2.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
    ax2.set_xlim( 9000, 20000 )
    setp(ax2.get_xticklabels(), visible=False)
    setp(ax2.get_yticklabels(), visible=False)
    ax2.set_ylabel('Flux / Transmission (arbitrary units)')

    ax3= subplot(3,1,3, sharex=ax1, sharey=ax1)
    plot(w160, f160, 'm--',label='F160W')
    plot(w153, f153, 'm-',label='F153M')
    plot(w1az, f1az, 'r-',label='_nolegend_')
    ax3.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
    setp(ax3.get_yticklabels(), visible=False)

    ax1.set_xlim( 9000, 20000 )
    ax1.set_xlabe

    l('observed wavelength (Angstroms)')

    fig = gcf()
    fig.subplots_adjust( wspace=0, hspace=0, left=0.05, bottom=0.12, right=0.95, top=0.95)



def plotbroadbandz(  zvals=[1,1.5,2.0], day=0 ):
    """ show how broad bands cover the SED at high z"""
    from hstsnpipe import tools
    from tools import snana
    w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.extrap.dat', day=day )
    print("SALT2")
    # w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/models/SALT2/SALT2.Guy10_UV2IR/salt2_template_0.dat', day=day )
    #w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/models/SALT2/SALT2.Guy10_UV2IR/salt2_template_1.dat', day=day )
    #wII, fII = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/non1a/SDSS-000018.DAT', day=0 )
    #wIb, fIb = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/non1a/SDSS-000020.DAT', day=0 )

    clf()

    i = 0
    for z in zvals: 
        i+=1
        w1az = w1a * (1+z)
        f1az = f1a / f1a.max() / 2.

        #wII = wII * (1+z)
        #fII = fII / fII.max() / 2.

        #wIb = wIb * (1+z)
        #fIb = fIb / fIb.max() / 2.

        ax = subplot(3,1,i)
        plot(w350, f350, 'b--', label='F350LP(W)')
        plot(w125, f125, 'g--', label='F125W(J)')
        plot(w160, f160, 'r--', label='F160W(H)')

        plot(w1az, f1az, 'k-', label='_nolegend_')
        #ax.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
        ax.set_xlim( 3000, 20000 )
        ax.text(0.98,0.95, 'z=%.1f'%(z), color='k',ha='right',va='top',transform=ax.transAxes)
        setp(ax.get_yticklabels(), visible=False)

        if i==1 : 
            top = ax.get_ylim()[1]
            ax.text(16000,top, 'F160W(H)', color='r',ha='center',va='bottom')
            ax.text(12500,top, 'F125W(J)', color='g',ha='center',va='bottom')
            ax.text(3500,top, 'F350LP(W)', color='b',ha='left',va='bottom')
        if i<3 : 
            setp(ax.get_xticklabels(), visible=False)
        if i==2 : 
            ax.set_ylabel('Flux or Transmission (arbitrary units)')
        if i==3 : 
            ax.set_xlabel('observed wavelength (Angstroms)')

    fig = gcf()
    fig.subplots_adjust( wspace=0, hspace=0, left=0.05, bottom=0.12, right=0.95, top=0.95)


def plotBVRI( ):
    """ show how broad ACS bands cover the SN SED """
    from hstsnpipe import tools
    from tools import snana
    w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.extrap.dat', day=0 )
    clf()

    f1a = f1a / f1a.max()
    plot(wB, fB, 'b--', label='B')
    plot(wV, fV, 'g--', label='V')
    plot(wR, fR, 'r--', label='R')
    plot(wI, fI, 'k--', label='I')

    plot(w435, f435, 'b-', label='F435W')
    plot(w606, f606, 'g-', label='F606W')
    plot(w625, f625, 'r-', label='F625W')
    plot(w814, f814, 'k-', label='F814W')

    plot(w1a, f1a, 'k-', label='_nolegend_')
    ax = gca()
    ax.set_xlim( 3000, 10000 )
    #setp(ax.get_yticklabels(), visible=False)
