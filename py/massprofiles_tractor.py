#!/usr/bin/env python

"""Fit the compact starbursts with The Tractor.

"""

from __future__ import division, print_function

import os
import sys
import logging
import argparse

import fitsio
import numpy as np
import matplotlib.pyplot as plt

from tractor.engine import Tractor, Image
from tractor.psfex import PsfEx, PixelizedPsfEx
from tractor.basics import NanoMaggies, PointSource, RaDecPos, LinearPhotoCal

from astrometry.util.util import Tan

class HizEA():
    """Class to handle the HizEA HST image cutouts."""
    def __init__(self, galaxy, scale, band):
        self.topdir = os.path.join(os.getenv('HIZEA_DATA'),'hst','aplus','cutouts')
        self.imfile = os.path.join(self.topdir,'{}-{}-{}.fits'.format(galaxy,scale,band))
        self.ivarfile = os.path.join(self.topdir,'{}-{}-{}-ivar.fits'.format(galaxy,scale,band))
        self.psffile = os.path.join(self.topdir,'{}-{}-{}.psf'.format(galaxy,scale,band))
        
        print('Reading image {}'.format(self.imfile))
        self.image = fitsio.read(self.imfile, ext=0)
        self.header = fitsio.read_header(self.imfile)

        print('Reading image {}'.format(self.ivarfile))
        self.ivar = fitsio.read(self.ivarfile, ext=0)

    def get_wcs(self):
        """Extract the WCS from the header."""
        wcs = Tan(self.imfile)
        return wcs

def fit_psf(galaxy,scale,band,verbose=0):

    # Set the debugging level.
    if verbose==0:
        lvl = logging.INFO
    else:
        lvl = logging.DEBUG
    logging.basicConfig(level=lvl,format='%(message)s',stream=sys.stdout)

    stampsize = 50

    pngprefix = 'qapsf-{}-{}-{}'.format(galaxy,scale,band)

    # Read the image and get a postage stamp of the galaxy (which should be
    # centered in the image)
    hizea = HizEA(galaxy,scale,band)
    wcs = hizea.get_wcs()
    W,H = wcs.imagew, wcs.imageh

    xpos,ypos = W/2, H/2
    ra,dec = wcs.pixelxy2radec(xpos, ypos)
    ix,iy = int(xpos), int(ypos)

    slc = (slice(max(iy-stampsize, 0), min(iy+stampsize+1, H)),
           slice(max(ix-stampsize, 0), min(ix+stampsize+1, W)))
    img = hizea.image[slc]
    
    # Instantiate the PSFEx model of the PSF at the center
    psf = PixelizedPsfEx(hizea.psffile).constantPsfAt(xpos,ypos)

    star = PointSource(RaDecPos(ra,dec), NanoMaggies(**{band: 1.0}))

    # Fit just the source RA,Dec,flux.
    tractor = Tractor([img], [star])
    tractor.freezeParam('images')
    
    for step in range(50):
        dlnp,X,alpha = tractor.optimize()
        if dlnp < 0.1:
            break

        print('Fit:', star)
        model_psfex = tractor.getModelImage(0)
        chi2_psfex = -2.0*tractor.getLogLikelihood()
        mag_psfex = NanoMaggies.nanomaggiesToMag(star.brightness)[0]

        #mn, mx = np.percentile((stamp-model_psfex)[ivarstamp>0],[1,95])
        sig = np.std((stamp-model_psfex)[ivarstamp>0])
        mn, mx = [-2.0*sig,5*sig]



    sys.exit(0)

    psf = hizea.get_psf()

    for istar, ps1star in enumerate(cat):
        ra, dec = (ps1star.ra, ps1star.dec)
        mag = ps1star.median[ps1band[band]] # r-band

        ok, xpos, ypos = wcs.radec2pixelxy(ra, dec)
        ix,iy = int(xpos), int(ypos)

        # create a little tractor Image object around the star
        slc = (slice(max(iy-stampsize, 0), min(iy+stampsize+1, H)),
               slice(max(ix-stampsize, 0), min(ix+stampsize+1, W)))

        # The PSF model 'const2Psf' is the one used in DR1: a 2-component
        # Gaussian fit to PsfEx instantiated in the image center.
        tim = im.get_tractor_image(slc=slc, const2psf=True)
        stamp = tim.getImage()
        ivarstamp = tim.getInvvar()

        # Initialize a tractor PointSource from PS1 measurements
        flux = NanoMaggies.magToNanomaggies(mag)
        star = PointSource(RaDecPos(ra,dec), NanoMaggies(**{band: flux}))

        # Fit just the source RA,Dec,flux.
        tractor = Tractor([tim], [star])
        tractor.freezeParam('images')

        print('2-component MOG:', tim.psf)
        tractor.printThawedParams()

        for step in range(50):
            dlnp,X,alpha = tractor.optimize()
            if dlnp < 0.1:
                break
        print('Fit:', star)
        model_mog = tractor.getModelImage(0)
        chi2_mog = -2.0*tractor.getLogLikelihood()
        mag_mog = NanoMaggies.nanomaggiesToMag(star.brightness)[0]


        # Generate a QAplot.
        if (istar>0) and (istar%(ncols)==0):
            irow = irow+1
        icol = 3*istar - 3*ncols*irow
        #print(istar, irow, icol, icol+1, icol+2)

        ax1 = plt.subplot2grid((nrows,3*ncols), (irow,icol), aspect='equal')
        ax1.axis('off')
        #ax1.imshow(stamp, **tim.ima)
        ax1.imshow(stamp,cmap='gray',interpolation='nearest',
                   origin='lower',vmin=mn,vmax=mx)
        ax1.text(0.1,0.9,'{:2d}'.format(istar+1),color='white',
                horizontalalignment='left',verticalalignment='top',
                transform=ax1.transAxes)

        ax2 = plt.subplot2grid((nrows,3*ncols), (irow,icol+1), aspect='equal')
        ax2.axis('off')
        #ax2.imshow(stamp-model_mog, **tim.ima)
        ax2.imshow(stamp-model_mog,cmap='gray',interpolation='nearest',
                   origin='lower',vmin=mn,vmax=mx)
        ax2.text(0.1,0.9,'MoG',color='white',
                horizontalalignment='left',verticalalignment='top',
                transform=ax2.transAxes)
        ax2.text(0.08,0.08,'{:.3f}'.format(mag_mog),color='white',
                 horizontalalignment='left',verticalalignment='bottom',
                 transform=ax2.transAxes)

        #ax2.set_title('{:.3f}, {:.2f}'.format(mag_psfex,chi2_psfex),fontsize=14)
        #ax2.set_title('{:.3f}, $\chi^{2}$={:.2f}'.format(mag_psfex,chi2_psfex))

        ax3 = plt.subplot2grid((nrows,3*ncols), (irow,icol+2), aspect='equal')
        ax3.axis('off')
        #ax3.imshow(stamp-model_psfex, **tim.ima)
        ax3.imshow(stamp-model_psfex,cmap='gray',interpolation='nearest',
                   origin='lower',vmin=mn,vmax=mx)
        ax3.text(0.1,0.9,'PSFEx',color='white',
                horizontalalignment='left',verticalalignment='top',
                transform=ax3.transAxes)
        ax3.text(0.08,0.08,'{:.3f}'.format(mag_psfex),color='white',
                 horizontalalignment='left',verticalalignment='bottom',
                 transform=ax3.transAxes)

        if istar==(nstar-1):
            break
    fig.savefig(pngprefix+'-stargrid.png',bbox_inches='tight')

def main():
    """
    Main routine.
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--galaxy', type=long, default='J1341-0321', metavar='', 
                        help='galaxy name')
    parser.add_argument('-s', '--scale', type=str, default='65mas', metavar='', 
                        help='pixel scale')
    parser.add_argument('-b', '--band', type=long, default='F160W', metavar='', 
                        help='bandpass to model')
    parser.add_argument('-v', '--verbose', dest='verbose', action='count', default=0,
                        help='Toggle on verbose output')
    args = parser.parse_args()

    fit_psf(galaxy=args.galaxy,scale=args.scale,band=args.band,
            verbose=args.verbose)
    
if __name__ == "__main__":
    main()
