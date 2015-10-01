#! /usr/bin/env python
# S.Rodney 
# use the stsci iraf function drztranback to 
# convert drizzled coordinates


def drztranback( drzfile, x=0, y=0, fltlist=None, 
                 xylistfile='', verbose=True ):
    """ 
    convert a set of coordinates from a drz image into
    the coordinate space of each of the contributing _flt files
    returns a list of tuples with (fltfile,x,y)
    """
    import os
    import pyfits
    from pyraf import iraf
    from iraf import stsdas
    from numpy import loadtxt, iterable
    stsdas.analysis()
    stsdas.dither()

    # get output (i.e. drizzled) image size
    nxout = pyfits.getval( drzfile, 'NAXIS1')
    nyout = pyfits.getval( drzfile, 'NAXIS2')
    scaleout = 3600*( abs(pyfits.getval(drzfile,'CD1_1')) +
                     abs(pyfits.getval(drzfile,'CD2_2')) )

    xscaleout = 7200*abs(pyfits.getval(drzfile,'CD1_1'))
    yscaleout = 7200*abs(pyfits.getval(drzfile,'CD2_2'))

    # if needed, get a list of contributing flt files and/or coeff files
    if not fltlist : 
        fltlist = getfltlist( drzfile )

    # build the list of coordinates, starting with the original drzfile x,y coords 
    returnlist = []
    if xylistfile : 
        xlist,ylist = loadtxt( xylistfile, unpack=True )
    elif iterable(x) : 
        xlist, ylist = x, y
    else : 
        xlist,ylist = [float(x)],[float(y)]
    for xx,yy in zip(x,y) : 
        returnlist.append( (drzfile,1,xx,yy) )

    # write out a list of drz-frame x,y positions into a text file.
    # this is used as input to wtranback for all flt files
    xylistfile = util.naming.chsuffix( os.path.basename(drzfile), '_fake.xylist')
    fout = open(xylistfile, 'w') 
    for xx,yy in zip(xlist,ylist) : 
        print >> fout, "%15.5f %15.5f"%(xx,yy)
    fout.close()

    # translate the x,y coords back to _flt coordinates
    # TODO : probably need to allow for up to 2 coeff files 
    # for every flt file, for UVIS and ACS
    # for fltfile,coefffile in zip(fltlist,coefflist) : 
    for fltfile in fltlist: 
        if verbose: print("translating %s to %s coords"%(drzfile,fltfile))

        # find all the sci extensions
        flthdulist = pyfits.open( fltfile ) 
        extlist = [ i for i in range(len(flthdulist)) if flthdulist[i].name.lower().startswith('sci') ]

        for ext in extlist : 
            nxin = pyfits.getval( fltfile, 'NAXIS1', ext=ext)
            nyin = pyfits.getval( fltfile, 'NAXIS2', ext=ext)
            scalein = 3600*( abs(pyfits.getval(fltfile,'CD1_1',ext=ext)) +
                             abs(pyfits.getval(fltfile,'CD2_2',ext=ext)) )
            xscalein = 7200*abs(pyfits.getval(fltfile,'CD1_1',ext=ext))
            yscalein = 7200*abs(pyfits.getval(fltfile,'CD2_2',ext=ext))
        
            # slimmed down wtranback call from LS: 2011.04.28
            # 
            # NOTE: we use the xylistfile as input even when 
            # there is only one pair of coordinates
            iraf.wtranback.unlearn()
            iraf.flpr(); iraf.flpr()
            iraf.gflush(); iraf.gflush()

            #coeffile = os.path.join(snworkdir,fltfile[:-5] + '_coeffs1.dat')
            coeffile = fltfile[:-5] + '_coeffs1.dat'
            output = iraf.wtranback( 0, 0, nxin=nxin, nxout=nxout, 
                                     nyin=nyin, nyout=nyout, xylist=xylistfile,
                                     coeffs=coeffile,  
                                     geomode='wcs', refim=drzfile,
                                     inimage=fltfile+'[%i]'%ext, Stdout=1 )

            for line in output:
                if line.startswith(' Xin,Yin:') : 
                    xin,yin = map(float,line.split()[1:3])
                    if xin>0 and xin<nxin and yin>0 and yin<nyin :
                        returnlist.append( (fltfile,ext,float(xin),float(yin)) )
    return( returnlist )


def getfltlist( drzfile, mustexist=False ) :
    """ 
    extract a list of contributing flt files from the header of 
    the given drizzled image file. 
    """
    import os
    import pyfits
    import exceptions

    Nflt = pyfits.getval( drzfile, 'NDRIZIM', ext=0 )
    fltfilelist = []
    for i in range( Nflt ):
        fltfile = pyfits.getval( drzfile, 'D%03iDATA'%(i+1), ext=0 )
        fltfile = fltfile[:fltfile.find('[')]

        if mustexist and not os.path.isfile( fltfile ) : 
            raise exceptions.RuntimeError( "%s (contributing to %s) is missing"%(fltfile,drzfile) )

        if fltfile not in fltfilelist : 
            fltfilelist.append( fltfile )
    return( fltfilelist )
