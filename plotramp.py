#! /usr/bin/env python
# 2010.09.15
# S.Rodney
#
# Tools for examining HST WFC3 drz images
"""
Plot the up-the-ramp sampling for a list of pixel locations in a given
_drz image file.  Pixel locations (in the drz frame) are given as
series of comma-separated pairs or in a text file.  

Syntax :  plotramp.py drzfile [xylist.dat or x1,y1 x2,y2 ...]

"""

def main():
    import getopt
    import sys
    import os, subprocess
    import time

    verbose = True
    debug = False 
    clobber = False

    linfit=False

    # read in arguments and options
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","help","debug","clobber","linfit" ] )
    except getopt.GetoptError: 
        print "Error : incorrect option or missing argument."
        print __doc__
        return(-1)
    for o, a in opt:
        if o in ["-h", "--help"]:
            print __doc__
            return(0)
        elif o == "-v" :
            verbose = True
        elif o == "--verbose" :
            verbose = int(a)
        elif o == "--debug" :
            debug = True
        elif o == "--clobber" :
            clobber = True
        elif o == "--linfit" :
            linfit = True
    if debug: import pdb; pdb.set_trace()

    if len(arg) < 2 : 
        print __doc__
        sys.exit()
    imfile = arg[0]
    xylist = []
    for xy in arg[1:] :
        if xy.find(',')>0:
            xylist.append((float(xy.split(',')[0]),float(xy.split(',')[1])))
        else : 
            xylist = xy
            break

    # user provides drz.fits or dif.fits file name and coordinates
    # plots the ramp from the ima file for each contributing flt file
    multiramp( imfile, xylist, linfit=linfit )
    from pylab import show, savefig, ioff
    savefile = imfile[:-5]+"_plotramp.png"
    userin = raw_input("File name for saved figure [%s] : "%savefile)
    if len(userin): savefile = userin
    savefig( savefile )
    if os.uname()[0]=='Darwin':
        os.system('open %s'%savefile)


def multiramp( drzfile, xylist, linfit=False, debug=False) :
    """ plot a set of ramps from a list of x,y positions """
    from numpy import shape
    import pyfits
    xylistfile = ''
    if not len(shape(xylist))>1 :
        if not xylist.find(',')>0 : 
            xylistfile = xylist 
            xylist = [(0,0)]
        else : 
            xylist = [xy]

    # plot the ramp from every file at each position
    oldfile = ''
    markers = ['s','D','^','v','>','<','x','+','o','h']
    colors =  ['orchid','blue','green','darkorange','red','maroon','k']
    im,ic = -1,-1

    # import pdb; pdb.set_trace()
    if pyfits.getval(drzfile,'NEXTEND') > 1 : ext='SCI'
    else : ext=0

    for xy in xylist :
        xx,yy = xy
        # get the list of files and x,y positions
        filexylist = drztranback( drzfile, xx, yy, xylistfile=xylistfile, ext=ext)
        ic+=1 
        if ic>len(colors): ic=ic%len(colors)

        for imfile,xtran,ytran in filexylist:
            if imfile == drzfile : continue
            im+=1
            if im>len(markers): im=im%len(markers)
            imafile = imfile[:-9] + '_ima.fits'
            # add 5 to x,y coords to account for ref pixels
            sci,err,tsamp,nsamp = getrampdat( imafile, round(xtran+5), round(ytran+5) )
            plotramp(sci,err,tsamp,nsamp, plotfit=linfit, 
                     marker=markers[im], color=colors[ic],
                     label="%s [%i,%i]"%(imafile.rstrip('.fits'),round(xtran),round(ytran)) )
            

def fitramp(  sci,err,tsamp,nsamp ):
    from scipy import optimize as so
    param = [100., 0.1]
    minparam = so.fmin( chi2, param, args=(sci,err,tsamp,nsamp) )
    return( minparam )

def chi2( linparam, sci,err,tsamp,nsamp):
    #import pdb; pdb.set_trace()
    m,b = linparam
    fit = m * nsamp + b
    chi2sum = sum( (sci[1:]*tsamp[1:] - fit[1:])**2/(err[1:]*tsamp[1:])**2 )
    return( chi2sum )


def plotramp( sci,err,tsamp,nsamp, plotfit=False, **kwargs) :
    """ plot the up-the-ramp sampling sequence 
    for the given pixel, given in _ima coordinates.
    kwargs are passed to pylab.errorbar().
    """
    from pylab import plot, errorbar, legend, xlabel, ylabel, axes, axhline
    #sci, err, tsamp, nsamp  = getrampdat(imafile, x, y )
    if plotfit : 
        ax1 = axes( [0.1,0.35,0.85,0.6] )
        m,b = fitramp(  sci,err,tsamp,nsamp )
        fit = m * nsamp + b
        plot( nsamp, fit, ls='--', color='k', marker='')
        errorbar( nsamp, sci*tsamp, err*tsamp, **kwargs )
        ylabel('cumulative counts (sci*tsamp)')
        legend( loc='upper left')

        ax2 = axes( [0.1,0.1,0.85,0.25], sharex=ax1 )
        kwargs['ls']=' '
        errorbar( nsamp[1:], sci[1:]*tsamp[1:]-fit[1:], err[1:]*tsamp[1:],
                   **kwargs )
        axhline(ls='--',color='k')

        xlabel('SAMP NUMBER')
    else : 
        errorbar( nsamp, sci*tsamp, err*tsamp, **kwargs )
        xlabel('SAMP NUMBER')
        ylabel('cumulative counts (sci*tsamp)')
        legend( loc='upper left')

def getrampdat( imafile, x, y ):
    """ read the values of the science array from the 
    ima file at pixel x,y.  Return an array of pixel values
    as a function of sample number.
    """
    import pyfits
    from numpy import array

    ima = pyfits.open( imafile )
    nextend = ima[0].header['NEXTEND']
    NSAMP   = ima[0].header['NSAMP']
    
    tsamp = []
    nsamp = []
    sci = []    
    err = []
    dq = []
    for isamp in range(NSAMP-1,-1,-1) :
        sci.append( ima[ 1+isamp*5 ].data[y-1,x-1] )
        tsamp.append( ima[ 1+isamp*5 ].header['SAMPTIME'] )
        err.append( ima[ 2+isamp*5 ].data[y-1,x-1] )
        nsamp.append( ima[ 1+isamp*5 ].header['SAMPNUM'] )
    return( array(sci), array(err), array(tsamp), array(nsamp) )


def drztranback( drzfile, x=0, y=0, xylistfile='', ext='SCI',verbose=False ):
    """ 
    convert a set of coordinates from a drz image into
    the coordinate space of each of the contributing _flt files
    returns a list of tuples with (fltfile,x,y)
    """
    import os
    import pyfits
    from pyraf import iraf
    from iraf import stsdas
    from numpy import loadtxt
    stsdas.analysis()
    stsdas.dither()

    # get output (i.e. drizzled) image size
    nxout = pyfits.getval( drzfile, 'NAXIS1', ext=ext)
    nyout = pyfits.getval( drzfile, 'NAXIS2', ext=ext)
    scaleout = 3600*( abs(pyfits.getval(drzfile,'CD1_1',ext=ext)) +
                      abs(pyfits.getval(drzfile,'CD2_2',ext=ext)) )
    
    xscaleout = 7200*abs(pyfits.getval(drzfile,'CD1_1',ext=ext))
    yscaleout = 7200*abs(pyfits.getval(drzfile,'CD2_2',ext=ext))
    
    # get a list of contributing flt files
    fltfilelist = getfltlist( drzfile )

    scifile = drzfile[:-8]+"sci.fits"
    
    # translate the x,y coords back to _flt coordinates
    returnlist = []
    if xylistfile : 
        xlist,ylist = loadtxt( xylistfile, unpack=True )
        for x,y in zip(xlist,ylist) : 
            returnlist.append( (drzfile, x, y))
    else : 
        returnlist.append( (drzfile, x, y) )
    for fltfile in fltfilelist : 
        if verbose: print("translating %s to %s coords"%(drzfile,fltfile))
        # import pdb; pdb.set_trace()
        nxin = pyfits.getval( fltfile, 'NAXIS1', ext='SCI')
        nyin = pyfits.getval( fltfile, 'NAXIS2', ext='SCI')
        scalein = 3600*( abs(pyfits.getval(fltfile,'CD1_1',ext='SCI')) +
                    abs(pyfits.getval(fltfile,'CD2_2',ext='SCI')) )
        xscalein = 7200*abs(pyfits.getval(fltfile,'CD1_1',ext='SCI'))
        yscalein = 7200*abs(pyfits.getval(fltfile,'CD2_2',ext='SCI'))
        coeffile = fltfile[:-5] + '_coeffs1.dat'

        # wtranback is better than tranback
        iraf.flpr()
        iraf.flpr()
        iraf.unlearn( iraf.wtranback )
        if ext=='SCI' : extstr = '[1]'
        else : extstr = ''
        output = iraf.wtranback( x, y, nxin=nxin, nxout=nxout, 
                                 nyin=nyin, nyout=nyout, xylist=xylistfile,
                                 coeffs=coeffile,                                 
                                 geomode='wcs',
                                 refim=drzfile+extstr,
                                 inimage=fltfile+'[sci,1]',
                                 raref =pyfits.getval(drzfile,'CRVAL1',ext=ext),
                                 decref=pyfits.getval(drzfile,'CRVAL2',ext=ext),
                                 xrefpix=pyfits.getval(drzfile,'CRPIX1',ext=ext),
                                 yrefpix=pyfits.getval(drzfile,'CRPIX2',ext=ext),
                                 orient=pyfits.getval(drzfile,'ORIENTAT',ext=ext),
                                 Stdout=1)
        for line in output:
            if line.startswith(' Xin,Yin:') : 
                xin,yin = line.split()[1:3]
                returnlist.append( (fltfile,float(xin),float(yin)) )
    return( returnlist )


def getfltlist( drzfile ) :
    """get names of contributing _flt files from the _drz.fits header"""
    import os
    import pyfits
    import exceptions
    Nflt = pyfits.getval( drzfile, 'NDRIZIM', ext=0 )
    fltfilelist = []
    for i in range( Nflt ):
        fltfile = pyfits.getval( drzfile, 'D%03iDATA'%(i+1), ext=0 )
        fltfile = fltfile[:fltfile.find('[')]
        if not os.path.isfile( fltfile ):
            raise exceptions.RuntimeError( "%s (contributing to %s) is missing"%(fltfile,drzfile) )
        fltfilelist.append( fltfile )
    return( fltfilelist )


if __name__=='__main__':
    main()
