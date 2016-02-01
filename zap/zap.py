# ZAP - Zurich Atmosphere Purge
#
#    Copyright (C) 2014  Kurt Soto, Simon Lilly
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import, division, print_function

import astropy.units as u
import logging
import numpy as np
import os
import sys

from astropy.io import fits
from astropy.wcs import WCS
from functools import wraps
from multiprocessing import cpu_count, Manager, Process
from scipy import ndimage as ndi
from scipy.stats import sigmaclip
from time import time

from .version import __version__

# Limits if the segments in Angstroms
SKYSEG = [0, 5400, 5850, 6440, 6750, 7200, 7700, 8265, 8602, 8731, 9275, 10000]

NCPU = cpu_count()
PY2 = sys.version_info[0] == 2

if not PY2:
    text_type = str
    string_types = (str,)
else:
    text_type = unicode
    string_types = (str, unicode)

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO,
                    stream=sys.stdout)
logger = logging.getLogger(__name__)


###############################################################################
################################### Top Level Functions #######################
###############################################################################


def process(musecubefits, outcubefits='DATACUBE_FINAL_ZAP.fits', clean=True,
            zlevel='median', cftype='weight', cfwidthSVD=100, cfwidthSP=50,
            pevals=[], nevals=[], optimizeType='normal', extSVD=None,
            skycubefits=None, svdoutputfits='ZAP_SVD.fits', mask=None,
            interactive=False):
    """ Performs the entire ZAP sky subtraction algorithm.

    Work on an input FITS file and optionally writes the product to an output
    FITS file.

    Parameters
    ----------

    musecubefits : str
        Input FITS file, containing a cube with data in the first extension.
    outcubefits : str
        Output FITS file, based on the input one to propagate all header
        information and other extensions. Default to `DATACUBE_FINAL_ZAP.fits`.
    clean : bool
        If True (default value), the NaN values are cleaned. Spaxels with more
        then 25% of NaN values are removed, the others are replaced with an
        interpolation from the neighbors.
    zlevel : str
        Method for the zeroth order sky removal: `none`, `sigclip` or `median`
        (default).
    cftype : str
        Method for the continuum filter: `median` or `weight` (default). For
        the `weight` method, a zeroth order sky is required (see `zlevel`).
    cfwidthSVD : int or float
        Window size for the continuum filter, for the SVD computation.
        Default to 100.
    cfwidthSP : int or float
        Window size for the continuum filter. Default to 50.
    optimizeType : str
        Optimization method to compute the number of eigenspectra used for each
        segment: `none`, `normal` (default), `enhanced`. If `none`, the number
        of eigenspectra must be specified with `nevals` or `pevals`, otherwise
        `normal` is used.
    pevals : list
        Allow to specify the percentage of eigenspectra used for each segment.
    nevals : list
        Allow to specify the number of eigenspectra used for each segment.
    extSVD : str
        Path of an input FITS file containing a SVD computed by the
        :func:`~zap.SVDoutput` function. Otherwise the SVD is computed.
    skycubefits : str
        Path for the optional output of the sky that is subtracted from the
        cube. This is simply the input cube minus the output cube.
    svdoutputfits : str
        Output FITS file. Default to `ZAP_SVD.fits`.
    interactive : bool
        If True, a :class:`~zap.zclass` object containing all information on
        the ZAP process is returned, and can be used to explore the
        eigenspectra and recompute the output (with the
        :meth:`~zap.zclass.reprocess` method). In this case, the output files
        are not saved (`outcubefits` and `skycubefits` are ignored). Default to
        False.

    """
    if not isinstance(musecubefits, string_types):
        raise TypeError('The process method only accepts a single datacube '
                        'filename.')

    # make sure it has the right extension
    outcubefits = outcubefits.split('.fits')[0] + '.fits'

    # check if outcubefits/skycubefits exists before beginning
    check_file_exists(outcubefits)
    check_file_exists(skycubefits)
    check_file_exists(svdoutputfits)

    # Check for consistency between weighted median and zlevel keywords
    if cftype == 'weight' and zlevel == 'none':
        raise ValueError('Weighted median requires a zlevel calculation')

    if optimizeType not in ('none', 'normal', 'enhanced'):
        raise ValueError('Invalid value for optimizeType')

    if extSVD is not None and mask is not None:
        raise ValueError('extSVD and mask parameters are incompatible: if mask'
                         ' must be used, then the SVD has to be recomputed')

    if mask is not None or cfwidthSVD != cfwidthSP:
        # In this case we have to run SVDoutput first to compute the SVD
        SVDoutput(musecubefits, svdoutputfits=svdoutputfits,
                  clean=clean, zlevel=zlevel, cftype=cftype,
                  cfwidth=cfwidthSVD, mask=mask)
        extSVD = svdoutputfits

    zobj = zclass(musecubefits)
    zobj._run(clean=clean, zlevel=zlevel, cfwidth=cfwidthSP, cftype=cftype,
              pevals=pevals, nevals=nevals, optimizeType=optimizeType,
              extSVD=extSVD)

    if interactive:
        # Return the zobj object without saving files
        return zobj

    if zobj.run_zlevel != 'extSVD' and svdoutputfits is not None:
        # Save SVD only if it was computed in _run, i.e. if an external SVD
        # was not given
        zobj.writeSVD(svdoutputfits=svdoutputfits)
    if skycubefits is not None:
        zobj.writeskycube(skycubefits=skycubefits)

    zobj.mergefits(outcubefits)


def SVDoutput(musecubefits, svdoutputfits='ZAP_SVD.fits', clean=True,
              zlevel='median', cftype='weight', cfwidth=100, mask=None):
    """ Performs the SVD decomposition of a datacube.

    This allows to use the SVD for a different datacube.

    Parameters
    ----------

    musecubefits : str
        Input FITS file, containing a cube with data in the first extension.
    svdoutputfits : str
        Output FITS file. Default to ZAP_SVD.fits
    clean : bool
        If True (default value), the NaN values are cleaned. Spaxels with more
        then 25% of NaN values are removed, the others are replaced with an
        interpolation from the neighbors.
    zlevel : str
        Method for the zeroth order sky removal: `none`, `sigclip` or `median`
        (default).
    cftype : str
        Method for the continuum filter: `median` or `weight` (default). For
        the `weight` method, a zeroth order sky is required (see `zlevel`).
    cfwidth : int or float
        Window size for the continuum filter, default to 300.
    mask : str
        Path of a FITS file containing a mask (1 for objects, 0 for sky).

    """
    logger.info('Running ZAP %s !', __version__)
    logger.info('Processing %s to compute the SVD', musecubefits)
    check_file_exists(svdoutputfits)

    # Check for consistency between weighted median and zlevel keywords
    if cftype == 'weight' and zlevel == 'none':
        raise ValueError('Weighted median requires a zlevel calculation')

    zobj = zclass(musecubefits)

    # clean up the nan values
    if clean:
        zobj._nanclean()

    # if mask is supplied, apply it
    if mask is not None:
        zobj._applymask(mask)

    # Extract the spectra that we will be working with
    zobj._extract()

    # remove the median along the spectral axis
    if zlevel.lower() != 'none':
        zobj._zlevel(calctype=zlevel)

    # remove the continuum level - this is multiprocessed to speed it up
    zobj._continuumfilter(cftype=cftype, cfwidth=cfwidth)

    # do the multiprocessed SVD calculation
    zobj._msvd()

    # write to file
    zobj.writeSVD(svdoutputfits=svdoutputfits)


def contsubfits(musecubefits, contsubfn='CONTSUB_CUBE.fits', cfwidth=100):
    """ A multiprocessed implementation of the continuum removal.

    This process distributes the data to many processes that then reassemble
    the data. Uses two filters, a small scale (less than the line spread
    function) uniform filter, and a large scale median filter to capture the
    structure of a variety of continuum shapes.

    """
    check_file_exists(contsubfn)
    hdu = fits.open(musecubefits)
    data = hdu[1].data
    stack = data.reshape(data.shape[0], (data.shape[1] * data.shape[2]))
    contarray = _continuumfilter(stack, 'median', cfwidth=cfwidth)

    # remove continuum features
    stack -= contarray
    hdu[1].data = stack.reshape(data.shape[0], data.shape[1], data.shape[2])
    hdu.writeto(contsubfn)
    hdu.close()


def nancleanfits(musecubefits, outfn='NANCLEAN_CUBE.fits', rejectratio=0.25,
                 boxsz=1):
    """
    Detects NaN values in cube and removes them by replacing them with an
    interpolation of the nearest neighbors in the data cube. The positions in
    the cube are retained in nancube for later remasking.

    Parameters
    ----------

    musecubefits : str
        Input FITS file, containing a cube with data in the first extension.
    outfn : str
        Output FITS file. Default to NANCLEAN_CUBE.fits
    rejectratio : float
        Defines a cutoff for the ratio of NAN to total pixels in a spaxel
         before the spaxel is avoided completely. Default to 0.25
    boxsz : int
        Defines the number of pixels around the offending NaN pixel.
        Default to 1, which looks for the 26 nearest neighbors which
        is a 3x3x3 cube.

    """
    check_file_exists(outfn)
    hdu = fits.open(musecubefits)
    cleancube = _nanclean(hdu[1].data, rejectratio=rejectratio, boxsz=boxsz)
    hdu[1].data = cleancube[0]
    hdu.writeto(outfn)
    hdu.close()


def check_file_exists(filename):
    if filename is not None and os.path.exists(filename):
        raise IOError('Output file "{0}" exists'.format(filename))


def timeit(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        t0 = time()
        res = func(*args, **kwargs)
        logger.info('%s - Time: %.2f sec.', func.__name__, time() - t0)
        return res
    return wrapped


###############################################################################
##################################### Process Steps ###########################
###############################################################################

class zclass(object):

    """ Main class to run each of the steps of ZAP.

    Attributes
    ----------

    cleancube : numpy.ndarray
        The final datacube after removing all of the residual features.
    contarray : numpy.ndarray
        A 2D array containing the subtracted continuum per spaxel.
    cube : numpy.ndarray
        The original cube with the zlevel subtraction performed per spaxel.
    especeval : list of (eigenspectra, eval)
        A list containing the full set of eigenspectra and eigenvalues
        generated by the SVD calculation that is used toy reconstruct the
        entire datacube.
    laxis : numpy.ndarray
        A 1d array containing the wavelength solution generated from the header
        parameters.
    wcs : astropy.wcs.WCS
        WCS object with the wavelength solution.
    lranges : list
        A list of the wavelength bin limits used in segmenting the sepctrum
        for SVD.
    nancube : numpy.ndarray
        A 3d boolean datacube containing True in voxels where a NaN value was
        replaced with an interpolation.
    nevals : numpy.ndarray
        A 1d array containing the number of eigenvalues used per segment to
        reconstruct the residuals.
    normstack : numpy.ndarray
        A normalized version of the datacube decunstructed into a 2d array.
    varlist : numpy.ndarray
        An array for each segment with the variance curve, calculated for the
        optimize method.
    pranges : numpy.ndarray
        The pixel indices of the bounding regions for each spectral segment.
    recon : numpy.ndarray
        A 2d array containing the reconstructed emission line residuals.
    run_clean : bool
        Boolean that indicates that the NaN cleaning method was used.
    run_zlevel : bool
        Boolean indicating that the zero level correction was used.
    stack : numpy.ndarray
        The datacube deconstructed into a 2d array for use in the the SVD.
    subespeceval : list of (eigenspectra, eval)
        The subset of eigenvalues and eigenspectra used to reconstruct the sky
        residuals.
    variancearray : numpy.ndarray
        A list of length nsegments containing variances calculated per spaxel
        used for normalization
    y,x : numpy.ndarray
        The position in the cube of the spaxels that are in the 2d
        deconstructed stack
    zlsky : numpy.ndarray
        A 1d array containing the result of the zero level subtraction

    """

    def __init__(self, musecubefits):
        """ Initialization of the zclass.

        Pulls the datacube into the class and trims it based on the known
        optimal spectral range of MUSE.

        """
        hdu = fits.open(musecubefits)
        self.cube = hdu[1].data
        self.header = hdu[1].header
        self.musecubefits = musecubefits
        hdu.close()

        # Workaround for floating points errors in wcs computation: if cunit is
        # specified, wcslib will convert in meters instead of angstroms, so we
        # remove cunit before creating the wcs object
        header = self.header.copy()
        unit = u.Unit(header.pop('CUNIT3'))
        self.wcs = WCS(header).sub([3])

        # Create Lambda axis
        wlaxis = np.arange(self.cube.shape[0])
        self.laxis = self.wcs.all_pix2world(wlaxis, 0)[0]
        if unit != u.angstrom:
            # Make sure lambda is in angstroms
            self.laxis = (self.laxis * unit).to(u.angstrom).value

        # NaN Cleaning
        self.run_clean = False
        self.nancube = None
        self._boxsz = 1
        self._rejectratio = 0.25

        # zlevel parameters
        self.run_zlevel = False
        self.zlsky = np.zeros_like(self.laxis)

        # Extraction results
        self.stack = None
        self.y = None
        self.x = None

        # Normalization Maps
        self.contarray = None
        self.variancearray = None
        self.normstack = None

        # identify the spectral range of the dataset
        laxmin = min(self.laxis)
        laxmax = max(self.laxis)

        # List of segmentation limits in the optical
        skyseg = np.array(SKYSEG)
        skyseg = skyseg[(skyseg > laxmin) & (skyseg < laxmax)]

        # segment limit in angstroms
        self.lranges = (np.vstack([np.append(laxmin - 10, skyseg),
                                   np.append(skyseg, laxmax + 10)])).T

        # segment limit in pixels
        laxis = self.laxis
        lranges = self.lranges
        pranges = []
        for i in range(len(lranges)):
            paxis = wlaxis[(laxis > lranges[i, 0]) & (laxis <= lranges[i, 1])]
            pranges.append((np.min(paxis), np.max(paxis) + 1))
        self.pranges = np.array(pranges)

        # eigenspace Subset
        self.especeval = []
        self.subespeceval = []

        # Reconstruction of sky features
        self.recon = None
        self.cleancube = None
        self.varlist = None  # container for variance curves

    @timeit
    def _run(self, clean=True, zlevel='median', cftype='weight',
             cfwidth=100, pevals=[], nevals=[], optimizeType='normal',
             extSVD=None):
        """ Perform all zclass to ZAP a datacube:

        - NaN re/masking,
        - deconstruction into "stacks",
        - zerolevel subraction,
        - continuum removal,
        - normalization,
        - singular value decomposition,
        - eigenvector selection,
        - residual reconstruction and subtraction,
        - data cube reconstruction.

        """
        logger.info('Running ZAP %s !', __version__)

        self.optimizeType = optimizeType

        # clean up the nan values
        if clean:
            self._nanclean()

        # Extract the spectra that we will be working with
        self._extract()

        # remove the median along the spectral axis
        if extSVD is None:
            if zlevel.lower() != 'none':
                self._zlevel(calctype=zlevel)
        else:
            self._externalzlevel(extSVD)

        # remove the continuum level - this is multiprocessed to speed it up
        self._continuumfilter(cfwidth=cfwidth, cftype=cftype)

        # do the multiprocessed SVD calculation
        if extSVD is None:
            self._msvd()
        else:
            self._externalSVD(extSVD)

        # choose some fraction of eigenspectra or some finite number of
        # eigenspectra
        if optimizeType != 'none' or (nevals == [] and pevals == []):
            self.optimize()
            self.chooseevals(nevals=self.nevals)
        else:
            self.chooseevals(pevals=pevals, nevals=nevals)

        # reconstruct the sky residuals using the subset of eigenspace
        self.reconstruct()

        # stuff the new spectra back into the cube
        self.remold()

    # Clean up the nan value spaxels
    def _nanclean(self):
        """
        Detects NaN values in cube and removes them by replacing them with an
        interpolation of the nearest neighbors in the data cube. The positions
        in the cube are retained in nancube for later remasking.
        """
        self.cube, self.nancube = _nanclean(
            self.cube, rejectratio=self._rejectratio, boxsz=self._boxsz)
        self.run_clean = True

    @timeit
    def _extract(self):
        """
        Deconstruct the datacube into a 2d array, since spatial information is
        not required, and the linear algebra routines require 2d arrays.

        The operation rejects any spaxel with even a single NaN value, since
        this would cause the linear algebra routines to crash.

        Adds the x and y data of these positions into the zclass

        """
        logger.debug('Extracting to 2D')
        # make a map of spaxels with NaNs
        badmap = (np.logical_not(np.isfinite(self.cube))).sum(axis=0)
        # get positions of those with no NaNs
        self.y, self.x = np.where(badmap == 0)
        # extract those positions into a 2d array
        self.stack = self.cube[:, self.y, self.x]
        logger.debug('%d valid spaxels', len(self.x))

    def _externalzlevel(self, extSVD):
        """Remove the zero level from the extSVD file."""
        logger.info('Using external zlevel')
        self.zlsky = fits.getdata(extSVD, 0)
        self.stack -= self.zlsky[:, np.newaxis]
        self.run_zlevel = 'extSVD'

    @timeit
    def _zlevel(self, calctype='median'):
        """
        Removes a 'zero' level from each spectral plane. Spatial information is
        not required, so it operates on the extracted stack.

        Operates on stack, leaving it with this level removed and adds the data
        'zlsky' to the class. zlsky is a spectrum of the zero levels.

        This zero level is currently calculated with a median.

        Experimental operations -

        - exclude top quartile
        - run in an iterative sigma clipped mode

        """

        self.run_zlevel = calctype

        if calctype != 'none':
            logger.info('Subtracting Zero Level')

            zlstack = self.stack

            if calctype == 'median':
                logger.info('Median zlevel calculation')
                func = _imedian
            elif calctype == 'sigclip':
                logger.info('Iterative Sigma Clipping zlevel calculation')
                func = _isigclip

            self.zlsky = np.hstack(parallel_map(func, zlstack, NCPU, axis=0))
            self.stack -= self.zlsky[:, np.newaxis]
        else:
            logger.info('Skipping zlevel subtraction')

    def _continuumfilter(self, cfwidth=100, cftype='weight'):
        """ A multiprocessed implementation of the continuum removal.

        This process distributes the data to many processes that then
        reassemble the data.  Uses two filters, a small scale (less than the
        line spread function) uniform filter, and a large scale median filter
        to capture the structure of a variety of continuum shapes.

        added to class
        contarray - the removed continuua
        normstack - "normalized" version of the stack with the continuua
            removed

        """
        logger.info('Applying Continuum Filter, cfwidth=%d', cfwidth)
        if cftype not in ('weight', 'median'):
            raise ValueError("cftype must be 'weight' or 'median', got {}"
                             .format(cftype))
        self._cftype = cftype
        self._cfwidth = cfwidth

        if cftype == 'median':
            weight = None
        elif cftype == 'weight':
            weight = np.abs(self.zlsky - (np.max(self.zlsky) + 1))

        # remove continuum features
        self.contarray = _continuumfilter(self.stack, cftype, weight=weight,
                                          cfwidth=cfwidth)
        self.normstack = self.stack - self.contarray

    @timeit
    def _msvd(self):
        """ Multiprocessed singular value decomposition.

        First the normstack is normalized per segment per spaxel by the
        variance.  Takes the normalized, spectral segments and distributes them
        to the individual svd methods.

        """
        logger.info('Calculating SVD')

        # normalize the variance in the segments
        nseg = len(self.pranges)
        self.variancearray = var = np.zeros((nseg, self.stack.shape[1]))

        for i in range(nseg):
            pmin, pmax = self.pranges[i]
            var[i, :] = np.var(self.normstack[pmin:pmax, :], axis=0)
            self.normstack[pmin:pmax, :] /= var[i, :]

        logger.debug('Beginning SVD on %d segments', nseg)
        indices = [x[0] for x in self.pranges[1:]]
        self.especeval = parallel_map(_isvd, self.normstack, indices, axis=0)

    def chooseevals(self, nevals=[], pevals=[]):
        """ Choose the number of eigenspectra/evals to use for reconstruction.

        User supplies the number of eigen spectra to be used (neval) or the
        percentage of the eigenspectra that were calculated (peval) from each
        spectral segment to be used.

        The user can either provide a single value to be used for all segments,
        or provide an array that defines neval or peval per segment.

        """
        nranges = len(self.especeval)
        nevals = np.atleast_1d(nevals)
        pevals = np.atleast_1d(pevals)
        nespec = np.array([self.especeval[i][0].shape[1]
                           for i in range(nranges)])

        # deal with no selection
        if len(nevals) == 0 and len(pevals) == 0:
            logger.info('Number of modes not selected')
            nevals = np.array([1])

        # deal with an input list
        if len(nevals) > 1:
            if len(nevals) != nranges:
                nevals = np.array([nevals[0]])
                logger.info('Chosen eigenspectra array does not correspond to '
                            'number of segments')
            else:
                logger.info('Choosing %s eigenspectra for segments', nevals)

        if len(pevals) > 1:
            if len(pevals) != nranges:
                pevals = np.array([pevals[0]])
                logger.info('Chosen eigenspectra array does not correspond to '
                            'number of segments')
            else:
                logger.info('Choosing %s%% of eigenspectra for segments',
                            pevals)
                nevals = (pevals * nespec / 100.).round().astype(int)

        # deal with single value entries
        if len(pevals) == 1:
            logger.info('Choosing %s%% of eigenspectra for all segments',
                        pevals)
            nevals = (pevals * nespec / 100.).round().astype(int)
        elif len(nevals) == 1:
            logger.info('Choosing %s eigenspectra for all segments', nevals)
            nevals = np.zeros(nranges, dtype=int) + nevals

        # take subset of the eigenspectra and put them in a list
        subespeceval = []
        for i in range(nranges):
            eigenspectra, evals = self.especeval[i]
            tevals = (evals[0:nevals[i], :]).copy()
            teigenspectra = (eigenspectra[:, 0:nevals[i]]).copy()
            subespeceval.append((teigenspectra, tevals))

        self.subespeceval = subespeceval
        self.nevals = nevals

    @timeit
    def reconstruct(self):
        """Reconstruct the residuals from a given set of eigenspectra and
        eigenvalues
        """

        logger.info('Reconstructing Sky Residuals')
        nseg = len(self.especeval)
        rec = [(eig[:, :, np.newaxis] * ev[np.newaxis, :, :]).sum(axis=1)
               for eig, ev in self.subespeceval]

        # rescale to correct variance
        for i in range(nseg):
            rec[i] *= self.variancearray[i, :]
        self.recon = np.concatenate(rec)

    # stuff the stack back into a cube
    def remold(self):
        """ Subtracts the reconstructed residuals and places the cleaned
        spectra into the duplicated datacube.
        """
        logger.info('Applying correction and reshaping data product')
        self.cleancube = self.cube.copy()
        self.cleancube[:, self.y, self.x] = self.stack - self.recon
        if self.run_clean:
            self.cleancube[self.nancube] = np.nan

    # redo the residual reconstruction with a different set of parameters
    def reprocess(self, pevals=[], nevals=[]):
        """
        A method that redoes the eigenvalue selection, reconstruction, and
        remolding of the data.
        """

        self.chooseevals(pevals=pevals, nevals=nevals)
        self.reconstruct()
        self.remold()

    @timeit
    def optimize(self):
        """ Function to optimize the number of components used to characterize
        the residuals.

        This function calculates the variance per segment with an increasing
        number of eigenspectra/eigenvalues. It then deterimines the point at
        which the second derivative of this variance curve reaches zero. When
        this occurs, the linear reduction in variance is attributable to the
        removal of astronomical features rather than emission line residuals.

        """
        logger.info('Optimizing')

        normstack = self.stack - self.contarray
        nseg = len(self.especeval)
        self.nevals = np.zeros(nseg, dtype=int)
        indices = [x[0] for x in self.pranges[1:]]
        self.varlist = parallel_map(_ivarcurve, normstack, indices, axis=0,
                                    especeval=self.especeval,
                                    variancearray=self.variancearray)

        if self.optimizeType == 'enhanced':
            logger.info('Enhanced Optimization')
        else:
            logger.info('Normal Optimization')

        for i in range(nseg):
            # optimize
            varlist = self.varlist[i]
            deriv = varlist[1:] - varlist[:-1]
            deriv2 = deriv[1:] - deriv[:-1]
            noptpix = varlist.size

            if self.optimizeType != 'enhanced':
                # statistics on the derivatives
                ind = int(.5 * (noptpix - 2))
                mn1 = deriv[ind:].mean()
                std1 = deriv[ind:].std() * 2
                mn2 = deriv2[ind:].mean()
                std2 = deriv2[ind:].std() * 2
                # look for crossing points. When they get within 1 sigma of
                # mean in settled region.
                cross1 = np.append([False], deriv >= (mn1 - std1))  # pad by 1 for 1st deriv
                cross2 = np.append([False, False], np.abs(deriv2) <= (mn2 + std2))  # pad by 2 for 2nd
                cross = np.logical_or(cross1, cross2)
            else:
                # statistics on the derivatives
                ind = int(.75 * (noptpix - 2))
                mn1 = deriv[ind:].mean()
                std1 = deriv[ind:].std()
                mn2 = deriv2[ind:].mean()
                std2 = deriv2[ind:].std()
                cross = np.append([False], deriv >= (mn1 - std1))  # pad by 1 for 1st deriv

            self.nevals[i] = np.where(cross)[0][0]

    # #########################################################################
    # #################################### Extra Functions ####################
    # #########################################################################

    def make_contcube(self):
        """ Remold the continuum array so it can be investigated.

        Takes the continuum stack and returns it into a familiar cube form.
        """
        contcube = self.cube.copy() * np.nan
        contcube[:, self.y, self.x] = self.contarray
        return contcube

    def _externalSVD(self, extSVD):
        logger.info('Calculating eigenvalues for input eigenspectra')
        hdu = fits.open(extSVD)
        nseg = len(self.pranges)

        # normalize the variance in the segments
        self.variancearray = np.zeros((nseg, self.stack.shape[1]))

        for i in range(nseg):
            pmin, pmax = self.pranges[i]
            self.variancearray[i, :] = np.var(self.normstack[pmin:pmax, :],
                                              axis=0)
            self.normstack[pmin:pmax, :] /= self.variancearray[i, :]

        especeval = []
        for i in range(nseg):
            eigenspectra = hdu[i + 1].data
            ns = self.normstack[self.pranges[i][0]:self.pranges[i][1]]
            evals = np.transpose(np.transpose(ns).dot(eigenspectra))
            especeval.append([eigenspectra, evals])

        self.especeval = especeval
        hdu.close()

    def _applymask(self, mask):
        """Apply a mask to the input data to provide a cleaner basis set.

        mask is >1 for objects, 0 for sky so that people can use sextractor.
        The file is read with ``astropy.io.fits.getdata`` which first tries to
        read the primary extension, then the first extension is no data was
        found before.

        """
        logger.info('Applying Mask for SVD Calculation from %s', mask)
        mask = fits.getdata(mask).astype(bool)
        nmasked = np.count_nonzero(mask)
        logger.info('Masking %d pixels (%d%%)', nmasked,
                    nmasked / np.prod(mask.shape) * 100)
        self.cube[:, mask] = np.nan

    ###########################################################################
    ##################################### Output Functions ####################
    ###########################################################################

    def writecube(self, outcubefits='DATACUBE_ZAP.fits'):
        """Write the processed datacube to an individual fits file."""

        check_file_exists(outcubefits)
        # fix up for writing
        outhead = _newheader(self)

        # create hdu and write
        outhdu = fits.PrimaryHDU(data=self.cleancube, header=outhead)
        outhdu.writeto(outcubefits)
        logger.info('Cube file saved to %s', outcubefits)

    def writeskycube(self, skycubefits='SKYCUBE_ZAP.fits'):
        """Write the processed datacube to an individual fits file."""

        check_file_exists(skycubefits)
        # fix up for writing
        outcube = self.cube - self.cleancube
        outhead = _newheader(self)

        # create hdu and write
        outhdu = fits.PrimaryHDU(data=outcube, header=outhead)
        outhdu.writeto(skycubefits)
        logger.info('Sky cube file saved to %s', skycubefits)

    def mergefits(self, outcubefits):
        """Merge the ZAP cube into the full muse datacube and write."""

        # make sure it has the right extension
        outcubefits = outcubefits.split('.fits')[0] + '.fits'
        check_file_exists(outcubefits)
        hdu = fits.open(self.musecubefits)
        hdu[1].header = _newheader(self)
        hdu[1].data = self.cleancube
        hdu.writeto(outcubefits)
        hdu.close()
        logger.info('Cube file saved to %s', outcubefits)

    def writeSVD(self, svdoutputfits='ZAP_SVD.fits'):
        """Write the SVD to an individual fits file."""

        check_file_exists(svdoutputfits)
        hdu = fits.HDUList([fits.PrimaryHDU(self.zlsky)])
        for i in range(len(self.pranges)):
            hdu.append(fits.ImageHDU(self.especeval[i][0]))
        # write for later use
        hdu.writeto(svdoutputfits)
        logger.info('SVD file saved to %s', svdoutputfits)

    def plotvarcurve(self, i=0, ax=None):
        if len(self.varlist) == 0:
            logger.warning('No varlist found. The optimize method must be '
                           'run first.')
            return

        # optimize
        deriv = (np.roll(self.varlist[i], -1) - self.varlist[i])[:-1]
        deriv2 = (np.roll(deriv, -1) - deriv)[:-1]
        noptpix = self.varlist[i].size

        if self.optimizeType == 'normal':
            # statistics on the derivatives
            mn1 = deriv[.5 * (noptpix - 2):].mean()
            std1 = deriv[.5 * (noptpix - 2):].std() * 2
            mn2 = deriv2[.5 * (noptpix - 2):].mean()
            std2 = deriv2[.5 * (noptpix - 2):].std() * 2
        else:
            # statistics on the derivatives
            mn1 = deriv[.75 * (noptpix - 2):].mean()
            std1 = deriv[.75 * (noptpix - 2):].std()
            mn2 = deriv2[.75 * (noptpix - 2):].mean()
            std2 = deriv2[.75 * (noptpix - 2):].std()

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(3, 1, figsize=[10, 15])

        ax1, ax2, ax3 = ax
        ax1.plot(self.varlist[i], linewidth=3)
        ax1.plot([self.nevals[i], self.nevals[i]],
                 [min(self.varlist[i]), max(self.varlist[i])])
        ax1.set_ylabel('Variance')

        ax2.plot(np.arange(deriv.size), deriv)
        ax2.plot([0, len(deriv)], [mn1, mn1], 'k')
        ax2.plot([0, len(deriv)], [mn1 - std1, mn1 - std1], '0.5')
        ax2.plot([self.nevals[i] - 1, self.nevals[i] - 1],
                 [min(deriv), max(deriv)])
        ax2.set_ylabel('d/dn Var')

        ax3.plot(np.arange(deriv2.size), np.abs(deriv2))
        ax3.plot([0, len(deriv2)], [mn2, mn2], 'k')
        ax3.plot([0, len(deriv2)], [mn2 + std2, mn2 + std2], '0.5')
        ax3.plot([self.nevals[i] - 2, self.nevals[i] - 2],
                 [min(deriv2), max(deriv2)])
        ax3.set_ylabel('(d^2/dn^2) Var')
        # ax3.set_xlabel('Number of Components')

        ax1.set_title('Segment {0}, {1} - {2} Angstroms'.format(
            i, self.lranges[i][0], self.lranges[i][1]))


###############################################################################
##################################### Helper Functions ########################
###############################################################################


def worker(f, i, chunk, out_q, err_q, kwargs):
    try:
        result = f(i, chunk, **kwargs)
    except Exception as e:
        err_q.put(e)
        return

    # output the result and task ID to output queue
    out_q.put((i, result))


def parallel_map(func, arr, indices, **kwargs):
    manager = Manager()
    out_q = manager.Queue()
    err_q = manager.Queue()
    jobs = []
    axis = kwargs.pop('axis', None)
    chunks = np.array_split(arr, indices, axis=axis)

    for i, chunk in enumerate(chunks):
        p = Process(target=worker, args=(func, i, chunk, out_q, err_q, kwargs))
        jobs.append(p)
        p.start()

    # gather the results
    for proc in jobs:
        proc.join()

    if not err_q.empty():
        # kill all on any exception from any one slave
        raise err_q.get()

    # Processes finish in arbitrary order. Process IDs double
    # as index in the resultant array.
    results = [None] * len(jobs)
    while not out_q.empty():
        idx, result = out_q.get()
        results[idx] = result

    return results


##### Continuum Filtering #####

@timeit
def _continuumfilter(stack, cftype, weight=None, cfwidth=300):
    if cftype == 'median':
        func = _icfmedian
        weight = None
    elif cftype == 'weight':
        func = _icfweight
    c = parallel_map(func, stack, NCPU, axis=1, weight=weight, cfwidth=cfwidth)
    return np.concatenate(c, axis=1)


def _icfweight(i, stack, weight=None, cfwidth=None):
    return np.array([wmedian(row, weight, cfwidth=cfwidth)
                     for row in stack.T]).T


def _icfmedian(i, stack, weight=None, cfwidth=None):
    ufilt = 3  # set this to help with extreme over/under corrections
    return ndi.median_filter(
        ndi.uniform_filter(stack, (ufilt, 1)), (cfwidth, 1))


def rolling_window(a, window):  # function for striding to help speed up
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def wmedian(spec, wt, cfwidth=100):
    """ Performs a weighted median filtering of a 1d spectrum

    Operates using a cumulative sum curve

    Parameters
    ----------

    spec : numpy.ndarray
        Input 1d spectrum to be filtered
    wt : numpy.ndarray
        A spectrum of equal length as the input array to provide the weights.
    cfwidth : int or float
        Window size for the continuum filter, for the SVD computation.
        Default to 100.

    """

    # ignore the warning (feature not a bug)
    old_settings = np.seterr(divide='ignore')
    spec = np.pad(spec, (cfwidth, cfwidth), 'constant', constant_values=0)
    wt = np.pad(wt, (cfwidth, cfwidth), 'constant',
                constant_values=(np.min(wt) / 1000., np.min(wt) / 1000.))

    # do some striding for speed
    swin = rolling_window(spec, cfwidth)  # create window container array
    wwin = rolling_window(wt, cfwidth)  # create window container array

    # sort based on data
    srt = np.argsort(swin, axis=-1)
    ind = np.ogrid[0:swin.shape[0], 0:swin.shape[1]]
    sdata = swin[ind[0], srt]
    swt = wwin[ind[0], srt]

    # calculate accumulated weights
    awt = np.cumsum(swt, axis=-1)

    # new weightsort for normalization and consideration of data
    nw = (awt - 0.5 * swt) / awt[:, -1][:, np.newaxis]

    # find the midpoint in the new weight sort
    s = np.argmin(np.abs(nw - 0.5), axis=-1)
    sl = np.arange(len(s))
    nws = nw[sl, s]
    nws1 = nw[sl, s - 1]

    f1 = (nws - 0.5) / (nws - nws1)
    f2 = (0.5 - nws1) / (nws - nws1)
    wmed = sdata[sl, s - 1] * f1 + sdata[sl, s] * f2
    width = cfwidth // 2
    wmed = wmed[width:-width - 1]
    np.seterr(old_settings['divide'])

    return wmed


# ### SVD #####

def _isvd(i, normstack):
    """
    Perform single value decomposition and Calculate PC amplitudes (projection)
    outputs are eigenspectra operates on a 2D array.

    eigenspectra = [nbins, naxes]
    evals = [naxes, nobj]
    data = [nbins, nobj]
    """

    inormstack = normstack.T
    U, s, V = np.linalg.svd(inormstack, full_matrices=0)
    eigenspectra = np.transpose(V)
    evals = inormstack.dot(eigenspectra)
    logger.info('Finished SVD Segment %d', i)
    return [eigenspectra, evals.T]


# ### OPTIMIZE #####


def _ivarcurve(i, istack, especeval=None, variancearray=None):
    """
    Reconstruct the residuals from a given set of eigenspectra and eigenvalues.

    this is a special version for caculating the variance curve. It adds the
    contribution of a single mode to an existing reconstruction.

    """

    iprecon = np.zeros_like(istack)
    eigenspectra, evals = especeval[i]
    ivarlist = []
    totalnevals = int(np.round(evals.shape[0] * 0.25))

    for nevals in range(totalnevals):
        if nevals % (totalnevals * .2) <= 1:
            logger.info('Seg %d: %d%% complete ',
                        i, int(nevals / (totalnevals - 1.) * 100.))

        eig = eigenspectra[:, nevals]
        ev = evals[nevals, :]
        # broadcast evals on evects and sum
        iprecon += (eig[:, np.newaxis] * ev[np.newaxis, :])

        icleanstack = istack - (iprecon * variancearray[i])
        # calculate the variance on the cleaned segment
        ivarlist.append(np.var(icleanstack))

    return np.array(ivarlist)


def _newheader(zobj):
    """Put the pertinent zap parameters into the header"""
    header = zobj.header.copy()
    header['COMMENT'] = 'These data have been ZAPped!'
    header.append(('ZAPvers', __version__, 'ZAP version'), end=True)
    # zlevel removal performed
    header.append(('ZAPzlvl', zobj.run_zlevel, 'ZAP zero level correction'))
    # Nanclean performed
    header['ZAPclean'] = (zobj.run_clean,
                          'ZAP NaN cleaning performed for calculation')
    # Continuum Filtering
    header['ZAPcftyp'] = (zobj._cftype, 'ZAP continuum filter type')
    header['ZAPcfwid'] = (zobj._cfwidth, 'ZAP continuum filter size')

    # number of segments
    nseg = len(zobj.pranges)
    header['ZAPnseg'] = (nseg, 'Number of segments used for ZAP SVD')

    # per segment variables
    for i in range(nseg):
        header['ZAPseg{0}'.format(i)] = (
            '{0}:{1}'.format(zobj.pranges[i][0], zobj.pranges[i][1] - 1),
            'spectrum segment (pixels)')
        header['ZAPnev{0}'.format(i)] = (zobj.nevals[i],
                                         'number of eigenvals/spectra used')

    return header


def _isigclip(i, istack):
    mn = []
    for col in istack:
        clipped, bot, top = sigmaclip(col, low=3, high=3)
        mn.append(clipped.mean())
    return np.array(mn)


def _imedian(i, istack):
    return np.median(istack, axis=1)


@timeit
def _nanclean(cube, rejectratio=0.25, boxsz=1):
    """
    Detects NaN values in cube and removes them by replacing them with an
    interpolation of the nearest neighbors in the data cube. The positions in
    the cube are retained in nancube for later remasking.

    """
    logger.info('Cleaning NaN values in the cube')
    cleancube = cube.copy()
    badcube = np.logical_not(np.isfinite(cleancube))        # find NaNs
    badmap = badcube.sum(axis=0)  # map of total nans in a spaxel

    # choose some maximum number of bad pixels in the spaxel and extract
    # positions
    badmask = badmap > (rejectratio * cleancube.shape[0])
    logger.info('Rejected %d spaxels with more than %.1f%% NaN pixels',
                np.count_nonzero(badmask), rejectratio * 100)

    # make cube mask of bad spaxels
    badcube &= (~badmask[np.newaxis, :, :])
    z, y, x = np.where(badcube)

    neighbor = np.zeros((z.size, (2 * boxsz + 1)**3))
    icounter = 0
    logger.info("Fixing %d remaining NaN pixels", len(z))

    # loop over samplecubes
    nz, ny, nx = cleancube.shape
    for j in range(-boxsz, boxsz + 1, 1):
        for k in range(-boxsz, boxsz + 1, 1):
            for l in range(-boxsz, boxsz + 1, 1):
                iz, iy, ix = z + l, y + k, x + j
                outsider = ((ix <= 0) | (ix >= nx - 1) |
                            (iy <= 0) | (iy >= ny - 1) |
                            (iz <= 0) | (iz >= nz - 1))
                ins = ~outsider
                neighbor[ins, icounter] = cleancube[iz[ins], iy[ins], ix[ins]]
                neighbor[outsider, icounter] = np.nan
                icounter = icounter + 1

    mn = np.ma.masked_invalid(neighbor)
    cleancube[z, y, x] = mn.mean(axis=1).filled(np.nan)
    return cleancube, badcube
