"""
2010
S.Rodney

Simple image operations. 

"""

def imstamp( image, xc, yc, nx, ny, ext=0, fillpix=0, saveas=None, clobber=False ):
    """ 
    cut a postage stamp out of the given image
    centered on (xc,yc) with pixel size nx x ny.
    the 'image' may be the filename of a fits file
    or a numpy array.

    TODO: If the stamp extends beyond the image boundary 
    then it is filled with pixels of value 'fillpix'. 

    TODO : update the wcs keywords to correct values

    TODO : mef extension handling
    """
    from numpy import ndarray
    import pyfits
    import os

    # read in the image data 
    if isinstance( image, str ) :
        imdat = pyfits.getdata( image, ext=ext )
    elif isinstance( image, ndarray ):
        imdat = image
    naxis1,naxis2 = imdat.shape

    #  TODO:
    # check if stamp will overrun the boundary
    # make a blank stamp with value fillpix
    # paste in good pixels from image

    # extract the stamp and return it
    stamp = imdat[ (yc-ny/2):(yc+ny/2), (xc-nx/2):(xc+nx/2)]

    if not saveas : 
        return( stamp )

    if os.path.exists( saveas ) and not clobber : 
        print("%s exists. Not clobbering.")
        return( stamp )
    elif os.path.exists( saveas ) and clobber :
        os.remove( saveas )

    hdr = pyfits.getheader( image, ext=ext )
    hdu = pyfits.PrimaryHDU( data=stamp, header=hdr )
    hdulist = pyfits.HDUList( hdus=[hdu])
    hdulist.writeto( saveas )
    return( saveas )

def imsec( image, sec, ext=0 ):
    """ 
    cut a rectangular section out of the 
    given image defined by the
    pixel coordinates sec=[x1,x2,y1,y2]

    TODO : mef extension handling
    TODO : error handling
    """
    from numpy import ndarray
    from pyfits import getdata

    x1,x2,y1,y2 = sec

    # read in the image data 
    if isinstance( image, str ) :
        imdat = getdata( image, ext=ext )
    elif isinstance( image, ndarray ):
        imdat = image

    # extract the stamp and return it
    sec = imdat[ y1:y2, x1:x2]
    return( sec )


def drz2sciwht( drzfile, scifile=None, whtfile=None, 
                verbose=False ):
    """ split up a MEF file produced by multidrizzle
    into the _sci and _wht components """
    import pyfits
    import os
    from util import naming

    drzfile = os.path.abspath(drzfile)
    if not scifile : scifile = naming.chsuffix(drzfile, 'sci')
    if not whtfile : whtfile = naming.chsuffix(drzfile, 'wht')

    if verbose: 
        print("Splitting up mef drizzled image %s "%drzfile)
        print("into single-extension drizzle products:")
        print("    %s\n    %s"%(scifile, whtfile) )
    drz = pyfits.open( drzfile )
    sci = pyfits.PrimaryHDU( data = drz['SCI'].data, 
                             header= drz['SCI'].header ) 
    if 'D001PIXF' in drz[0].header.keys():
        sci.header.update('PIXFRAC', drz[0].header['D001PIXF'],
                          comment="multidrizzle 'drop' parameter")
    if 'FILTER' in drz[0].header.keys():
        sci.header.update(
            'FILTER', drz[0].header['FILTER'],
            comment="element selected in filter wheel")
    if 'EXPTIME' in drz[0].header.keys():
        sci.header.update(
            'EXPTIME', drz[0].header['EXPTIME'],
            comment="exposure time")
    if 'INSTRUME' in drz[0].header.keys():
        sci.header.update(
            'INSTRUME', drz[0].header['INSTRUME'],
            comment="instrument")
    if 'DETECTOR' in drz[0].header.keys():
        sci.header.update(
            'DETECTOR', drz[0].header['DETECTOR'],
            comment="detector")

    sci.writeto( scifile, clobber=True )

    wht = pyfits.PrimaryHDU( data = drz['WHT'].data, 
                             header= drz['WHT'].header ) 
    if 'D001PIXF' in drz[0].header.keys():
        wht.header.update('PIXFRAC', drz[0].header['D001PIXF'],
                          comment="multidrizzle 'drop' parameter")
    if 'FILTER' in drz[0].header.keys():
        wht.header.update(
            'FILTER', drz[0].header['FILTER'],
            comment="element selected in filter wheel")
    if 'EXPTIME' in drz[0].header.keys():
        wht.header.update(
            'EXPTIME', drz[0].header['EXPTIME'],
            comment="exposure time")
    if 'INSTRUME' in drz[0].header.keys():
        wht.header.update(
            'INSTRUME', drz[0].header['INSTRUME'],
            comment="instrument")
    if 'DETECTOR' in drz[0].header.keys():
        wht.header.update(
            'DETECTOR', drz[0].header['DETECTOR'],
            comment="detector")

    wht.writeto( whtfile, clobber=True )
    
    return( scifile, whtfile )


def imsubtract( image1, image2, outfile=None, 
                clobber=False, verbose=False, debug=False):
    """
     construct a simple subtraction:  image2 - image1
     guards against different sized data arrays by assuming 
     that the lower left pixel (0,0) is the anchor point.
    """
    import os
    import pyfits
    from numpy import ndarray
    import exceptions

    if debug  : import pdb; pdb.set_trace()

    if outfile : 
        if os.path.isfile( outfile ) and not clobber : 
            print("%s exists. Not clobbering."%outfile)
            return( outfile )

    # read in the images
    if not os.path.isfile( image1 ) :
        raise exceptions.RuntimeError(
            "The image file %s is not valid."%image1 )
    im1head = pyfits.getheader( image1 )
    im1data = pyfits.getdata( image1 )

    if not os.path.isfile( image2 ) :
        raise exceptions.RuntimeError(
            "The image file %s is not valid."%image2 )
    im2head = pyfits.getheader( image2 )
    im2data = pyfits.getdata( image2 )

    # sometimes multidrizzle drops a pixel. Unpredictable.
    nx2,ny2 = im2data.shape
    nx1,ny1 = im1data.shape
    if nx2>nx1 or ny2>ny1 : 
        im2data = im2data[:min(nx1,nx2),:min(ny1,ny2)]
        im1data = im1data[:min(nx1,nx2),:min(ny1,ny2)]
    elif nx2<nx1 or ny2<ny1 : 
        im1data = im1data[:min(nx1,nx2),:min(ny1,ny2)]
        im2data = im2data[:min(nx1,nx2),:min(ny1,ny2)]

    diffim =  im2data - im1data
    
    if not outfile :
        return( diffim )
    else : 
        im2head.update("SRCIM1",image1,"First source image = template for subtraction")
        im2head.update("SRCIM2",image2,"Second source image = search epoch image")
        outdir = os.path.split( outfile )[0]
        if outdir and not os.path.isdir(outdir): 
            os.makedirs( outdir )
        pyfits.writeto( outfile, diffim, 
                        header=im2head,
                        clobber=clobber )
        return( outfile )


def imWeightedAve( image1, image2, weight1, weight2, outfile, clobber=False, verbose=False):
    """
     construct a weighted average of image1 and image2:

     (weight1*image1 + weight2*image2) / (weight1+weight2)

     Mean image is written to outfile.
    """
    import os
    import pyfits
    from numpy import ndarray, nan_to_num
    import exceptions

    if os.path.isfile(outfile)  : 
        if clobber : 
            os.unlink( outfile )
        else : 
            print( "%s exists. Not clobbering."%outfile )
            return( outfile )
        
    # read in the sci and wht images 
    im1hdr = pyfits.getheader( image1 )
    im1 = pyfits.getdata( image1 )
    im2 = pyfits.getdata( image2 )
    wht1 = pyfits.getdata( weight1 )
    wht2 = pyfits.getdata( weight2 )

    meanim = nan_to_num( (wht1*im1 + wht2*im2)/(wht1+wht2) )
    
    # TODO : make a useful header
    outdir = os.path.dirname( outfile )
    if not os.path.isdir(outdir): 
        os.makedirs( outdir )
    pyfits.writeto( outfile, meanim, header=im1hdr )
    return( outfile )


def imaverage( imagelist, outfile, 
               clobber=False, verbose=False):
    """
     construct a simple average of all images in the list.
     Assumes all input images have identical dimensions
     Returns name of outfile
    """
    import os
    import pyfits
    from numpy import where, ones, zeros, array, ndarray, nan_to_num,float32
    import exceptions

    if os.path.exists(outfile)  : 
        if clobber : 
            os.unlink( outfile )
        else : 
            print( "%s exists. Not clobbering."%outfile )
            return( outfile )

    # make empty arrays for components of the weighted average
    naxis1 = pyfits.getval( imagelist[0], 'NAXIS1')
    naxis2 = pyfits.getval( imagelist[0], 'NAXIS2')
    sumarray = zeros( [naxis2,naxis1], dtype=float32 )
    ncombinearray = zeros( [naxis2,naxis1], dtype=float32 )
        
    # construct the weighted average and update header keywords
    outhdr = pyfits.getheader( imagelist[0] )
    i = 1
    # import pdb; pdb.set_trace()
    for imfile in imagelist : 
        imdat = pyfits.getdata( imfile )
        sumarray +=  imdat 
        ncombinearray += where( imdat != 0 , ones(imdat.shape), zeros(imdat.shape) ) 
        outhdr.update("SRCIM%02i"%i,imfile,"source image %i, used in average "%i )
        i+= 1
    outscidat = where( ncombinearray > 0 , sumarray/ncombinearray, zeros(sumarray.shape) ) 
   
    outdir = os.path.dirname( outfile )
    if outdir : 
        if not os.path.isdir(outdir): 
            os.makedirs( outdir )
    pyfits.writeto( outfile, outscidat, header=outhdr )

    return( outfile )


def weightedAverage( imagelist, whtlist, outfile, outwht,
                     clobber=False, verbose=False):
    """
     construct a weighted average :

     (weight1*image1 + weight2*image2 + ...) / (weight1+weight2+...)

     And a composite weight map : 
       (weight1 + weight2 + ...) / Nweight

     Mean image is written to outfile.
     Assumes all input images have identical dimensions
     Returns name of outfile and whtfile
    """
    import os
    import pyfits
    from numpy import where, ones, zeros, array, ndarray, nan_to_num,float32
    import exceptions

    if os.path.exists(outfile)  : 
        if clobber : 
            os.unlink( outfile )
        else : 
            print( "%s exists. Not clobbering."%outfile )
            return( outfile, outwht )

    if os.path.exists(outwht)  : 
        if clobber : 
            os.unlink( outwht )
        else : 
            print( "%s exists. Not clobbering."%outwht )
            return( outfile, outwht )

    # make empty arrays for components of the weighted average
    naxis1 = pyfits.getval( imagelist[0], 'NAXIS1')
    naxis2 = pyfits.getval( imagelist[0], 'NAXIS2')
    sumarray = zeros( [naxis2,naxis1], dtype=float32 )
    ncombinearray = zeros( [naxis2,naxis1], dtype=float32 )
    whtarray = ones( [naxis2,naxis1], dtype=float32 )
        
    # construct the weighted average and update header keywords
    outhdr = pyfits.getheader( imagelist[0] )
    i = 1
    # import pdb; pdb.set_trace()
    for imfile, whtfile in zip( imagelist, whtlist ) : 
        imdat = pyfits.getdata( imfile )
        whtdat = pyfits.getdata( whtfile )
        sumarray +=  whtdat * imdat 
        whtarray += whtdat
        ncombinearray += where( whtdat > 0 , ones(whtdat.shape), zeros(whtdat.shape) ) 
        outhdr.update("SRCIM%02i"%i,imfile,"source image %i, used in weighted average "%i )
        i+= 1
    # outscidat = nan_to_num( sumarray / whtarray ) 
    # outscidat = sumarray / whtarray
    outscidat = where( ncombinearray > 0 , sumarray / whtarray, zeros(sumarray.shape) ) 
    outwhtdat = where( ncombinearray > 0 , whtarray/ncombinearray, zeros(whtarray.shape) ) 
    
    outdir = os.path.dirname( outfile )
    if outdir : 
        if not os.path.isdir(outdir): 
            os.makedirs( outdir )
    pyfits.writeto( outfile, outscidat, header=outhdr )
    pyfits.writeto( outwht,  outwhtdat, header=outhdr )

    return( outfile, outwht )


def minCombine( imagelist, whtlist, outfile, outwht,
                clobber=False, verbose=False):
    """
     combined images in image list by taking the minimum value for all
     pixels with weight > 0.
     Min-Combined image is written to outfile.
     Assumes all input images have identical dimensions
     Returns name of outfile and whtfile

    import os
    import pyfits
    from numpy import where, ones, zeros, array, ndarray, nan_to_num,float32, min
    import exceptions

    if os.path.exists(outfile)  : 
        if clobber : 
            os.unlink( outfile )
        else : 
            print( "%s exists. Not clobbering."%outfile )
            return( outfile, outwht )

    if os.path.exists(outwht)  : 
        if clobber : 
            os.unlink( outwht )
        else : 
            print( "%s exists. Not clobbering."%outwht )
            return( outfile, outwht )

    # TODO : how to take the minimum only when wht>0 ??

    # make empty arrays for components of the weighted average
    naxis1 = pyfits.getval( imagelist[0], 'NAXIS1')
    naxis2 = pyfits.getval( imagelist[0], 'NAXIS2')
    ncombinearray = zeros( [naxis2,naxis1], dtype=float32 )
    # minarray = zeros( [naxis2,naxis1], dtype=float32 )
    # whtarray = zeros( [naxis2,naxis1], dtype=float32 )
        
    # construct the weighted average and update header keywords
    outhdr = pyfits.getheader( imagelist[0] )
    i = 1
    # import pdb; pdb.set_trace()
    for imfile, whtfile in zip( imagelist, whtlist ) : 
        imdat = pyfits.getdata( imfile )
        whtdat = pyfits.getdata( whtfile )

        if not minarray : 
            minarray = imdat
            whtarray = whtdat
            continue

        minarray =  where( whtdat > 0 , min(minarray,imdat), zeros(whtdat.shape) ) 

whtdat * imdat 
        whtarray += whtdat
        ncombinearray += where( whtdat > 0 , ones(whtdat.shape), zeros(whtdat.shape) ) 
        outhdr.update("SRCIM%02i"%i,imfile,"source image %i, used in weighted average "%i )
        i+= 1
    # outscidat = nan_to_num( sumarray / whtarray ) 
    # outscidat = sumarray / whtarray
    outscidat = where( ncombinearray > 0 , sumarray / whtarray, zeros(sumarray.shape) ) 
    outwhtdat = where( ncombinearray > 0 , whtarray/ncombinearray, zeros(whtarray.shape) ) 
    
    outdir = os.path.dirname( outfile )
    if outdir : 
        if not os.path.isdir(outdir): 
            os.makedirs( outdir )
    pyfits.writeto( outfile, outscidat, header=outhdr )
    pyfits.writeto( outwht,  outwhtdat, header=outhdr )

    return( outfile, outwht )
    """




def imaverage2( image1, image2, 
               outfile=None, clobber=False, verbose=False):
    """
     construct a simple average of image1 and image2:

     (image1 + image2) / 2.0 

     Mean image is written to outfile.
    """
    import os
    import pyfits
    from numpy import ndarray
    import exceptions

    if os.path.isfile(outfile)  : 
        if clobber : 
            os.unlink( outfile )
        else : 
            print( "%s exists. Not clobbering."%outfile )
            return( outfile )

    im1data = pyfits.getdata( image1 )
    im1head = pyfits.getheader( image1 )
    im2data = pyfits.getdata( image2 )

    # sometimes multidrizzle drops a pixel. Unpredictable.
    nx2,ny2 = im2data.shape
    nx1,ny1 = im1data.shape
    if nx2>nx1 or ny2>ny1 : 
        im2data = im2data[:min(nx1,nx2),:min(ny1,ny2)]
        im1data = im1data[:min(nx1,nx2),:min(ny1,ny2)]
    elif nx2<nx1 or ny2<ny1 : 
        im1data = im1data[:min(nx1,nx2),:min(ny1,ny2)]
        im2data = im2data[:min(nx1,nx2),:min(ny1,ny2)]

    meanim =  (im2data + im1data)/2.
    
    # TODO : make a useful header
    im1head.update("SRCIM1",image1,"First source image  for imaverage")
    im1head.update("SRCIM2",image2,"Second source image for imaverage")
    outdir = os.path.dirname( outfile )
    if outdir : 
        if not os.path.isdir(outdir): 
            os.makedirs( outdir )
    pyfits.writeto( outfile, meanim, header=im1head )
    return( outfile )


def imsum( image1, image2, outfile=None, clobber=False, verbose=False):
    """
     add together image1 and image2
    """
    import os
    import pyfits
    from numpy import ndarray
    import exceptions

    # read in the images
    if isinstance( image1, str ):
        if not os.path.isfile( image1 ) :
            raise exceptions.RuntimeError(
                "The image file %s is not valid."%image1 )
        im1head = pyfits.getheader( image1 )
        im1data = pyfits.getdata( image1 )
        im1head.update("SRCIM1",image1,"First source image  for imsum")
    elif isinstance( image1, ndarray ):
        im1data = image1
        im1head = None
    else : 
        raise  exceptions.RuntimeError(
            "Provide a fits file name or numpy array for each image.")
    
    if isinstance( image2, str ):
        if not os.path.isfile( image2 ) :
            raise exceptions.RuntimeError(
                "The image file %s is not valid."%image2 )
        im2data = pyfits.getdata( image2 )
        im1head.update("SRCIM2",image2,"Second source image for imsum")
    elif isinstance( image2, ndarray ):
        im2data = image2
    else : 
        raise  exceptions.RuntimeError(
            "Provide a fits file name or numpy array for each image.")


    # sometimes multidrizzle drops a pixel. Unpredictable.
    nx2,ny2 = im2data.shape
    nx1,ny1 = im1data.shape
    if nx2>nx1 or ny2>ny1 : 
        im2data = im2data[:min(nx1,nx2),:min(ny1,ny2)]
        im1data = im1data[:min(nx1,nx2),:min(ny1,ny2)]
    elif nx2<nx1 or ny2<ny1 : 
        im1data = im1data[:min(nx1,nx2),:min(ny1,ny2)]
        im2data = im2data[:min(nx1,nx2),:min(ny1,ny2)]

    sumim =  im2data + im1data
    
    if outfile :
        # TODO : make a useful header
        outdir = os.path.split( outfile )[0]
        if outdir : 
            if not os.path.isdir(outdir): 
                os.makedirs( outdir )
        pyfits.writeto( outfile, sumim, 
                        header=im1head,
                        clobber=clobber )
        return( outfile )
    else : 
        return( sumim )


def imscaleflux( image1, scalefactor, outfile=None, clobber=False, verbose=False):
    """
    multiply the flux in image by scalefactor and save to outfile
    """
    import os
    import pyfits
    from numpy import ndarray
    import exceptions

    # read in the image
    if isinstance( image1, str ):
        if not os.path.isfile( image1 ) :
            raise exceptions.RuntimeError(
                "The image file %s is not valid."%image1 )
        im1head = pyfits.getheader( image1 )
        im1data = pyfits.getdata( image1 )
        im1head.update("FLXSCALE",scalefactor,"Flux scaling factor")
    elif isinstance( image1, ndarray ):
        im1data = image1
        im1head = None
    else : 
        raise  exceptions.RuntimeError("Provide a fits file name or numpy array for each image.")
    
    scaledim =  scalefactor * im1data
    
    if outfile :
        outdir = os.path.split( outfile )[0]
        if outdir : 
            if not os.path.isdir(outdir): 
                os.makedirs( outdir )
        pyfits.writeto( outfile, scaledim, header=im1head,clobber=clobber )
        return( outfile )
    else : 
        return( scaledim )
   


def int2bin( n, returntype='bitnumbers', bits=16 ) :
    """
    convert the integer n into a binary representation.
    The return value is either a list or a string according
    to the 'returntype' setting:
      'bitnumbers' : list of bit numbers 
      'bitvalues'  : list of bit values 
      'string'     : string of 1's and 0's 
                    (formatted in quartets,
                     most significant bit on the left)
    examples:
      >> int2bitflags( 2304, returntype='bitvalues' )
        [256,2048]
      >> int2bitflags( 2304, returntype='bitnumbers' )
        [9,12]
      >> int2bitflags( 2304, returntype='string' )
        '0000 
    """
    import exceptions
    if returntype not in ['bitnumbers','bitvalues','string']:
        raise exceptions.RuntimeError(
            "returntype must be one of ['bitnumbers','bitvalues','string']")
    bitlist = []
    vallist = []
    bitstr  = ''
    for y in range(bits-1, -1, -1):
        if (n >> y & 1) :
            bitlist.append( y+1 )
            vallist.append( 2**y ) 
            bitstr +=  '1' 
        else : 
            bitstr +=  '0'
        if not y%4 :
            bitstr +=  ' '
    if returntype=='bitnumbers': return(sorted(bitlist))
    elif returntype=='bitvalues': return(sorted(vallist))
    elif returntype=='string': return( bitstr )








def wcscopy( donor, recipient, dext=1, rext=1, verbose=False ):
    """ Copy the WCS header keywords from the donor 
    fits file into the header of the recipient fits file. 
    dext and rext specify the fits extensions to use for
    the donor and recipient, respectively.
    """
    import pyfits
    donfits = pyfits.open( donor, mode='readonly' )
    recfits = pyfits.open( recipient, mode='update' )
    donhdr = donfits[dext].header
    donfits.close()
    rechdr = recfits[rext].header

    wcskeylist = [ 'WCSAXES','CRPIX1 ','CRPIX2 ','CRVAL1 ','CRVAL2 ','CTYPE1 ','CTYPE2 ','CD1_1  ','CD1_2  ','CD2_1  ','CD2_2  ','LTV1   ','LTV2   ','LTM1_1 ','LTM2_2 ','ORIENTAT','RA_APER','DEC_APER','PA_APER','VAFACTOR' ]    

    Nupdated = 0
    after = 'NAXIS2'
    for key in wcskeylist : 
        if key in donhdr : 
            value = donhdr.get( key )
            if key in rechdr : 
                rechdr[key] = value 
            else : 
                rechdr.update( key, value, after=after )
            after = key 
            Nupdated += 1
    recfits.flush()
    recfits.close()
    if verbose : print( 'Updated %i WCS keywords in %s'%(Nupdated, recipient) )
    return( None ) 
