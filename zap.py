#ZAP - Zurich Atmosphere Purge
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

import numpy as np
from time import time
from scipy import ndimage as ndi
from multiprocessing import Pool
import multiprocessing
import matplotlib.pyplot as plt
import os
from scipy.stats import sigmaclip
from scipy import sparse
try:
    from astropy.io import fits as pyfits
except:
    import pyfits
import pdb
import numpy.ma as ma

##################################################################################################
################################### Top Level Functions ##########################################
##################################################################################################

def process(musecubefits, outcubefits='DATACUBE_FINAL_ZAP.fits', clean=False, zlevel='median',
            q=0, cftype='weight', cfwidth=100, pevals=[], nevals=[], optimize=True, silent=False,
            extSVD='', skycubefits='', rec_settings='', refilter = True,  enhanced_optimize=False):
    """
    Performs the entire ZAP sky subtraction algorithm on an input fits file and writes the
    product to an output fits file.

    """
    if type(musecubefits) == list:
        print 'The process method only accepts single datacubes.'
        print 'Multiple inputs for sky calculation can be performed with SVDoutput.'
        return

    outcubefits=outcubefits.split('.fits')[0]+'.fits' #make sure it has the right extension

    #check if outcubefits exists before beginning
    if os.path.exists(outcubefits):
        print 'output filename "{0}" exists'.format(outcubefits)
        return

    if os.path.exists(skycubefits):
        print 'output filename "{0}" exists'.format(skycubefits)
        return

    # Check for consistency between weighted median and zlevel keywords
    if cftype == 'weight' and zlevel == 'none':
        print 'weighted median requires a zlevel calculation'
        return
    
    zobj = zclass(musecubefits)

    if rec_settings != '':
        if rec_settings == 'filled' or rec_settings == 'extSVD':
            if extSVD == '':
                print 'Filled Field case requires external SVD'
                return
            enhanced_optimize = True
            refilter = True
            cfwidth = 20
            cftype = 'weight'
            print 'Using recommended settings for filled field case:'
            print "cfwidth = 20, cftype = 'weight', enhanced_optimize=True, refilter = True" 
        if rec_settings == 'sparse':
            print 'Using recommended settings for sparse field case:'
            print "cfwidth = 100, cftype = 'weight', enhanced_optimize=True, refilter = True" 
            refilter = True
            enhanced_optimize = True
            cfwidth = 100
            cftype = 'weight'

    if enhanced_optimize == True:
        optimize = True

    zobj._run(clean=clean, zlevel=zlevel, q=q, cfwidth=cfwidth, cftype=cftype,
              pevals=pevals, nevals=nevals, optimize=optimize, silent=silent, 
              enhanced_optimize=enhanced_optimize, extSVD=extSVD)
    
    if skycubefits != '':
        zobj.writeskycube(skycubefits=skycubefits)
    
    zobj.mergefits(outcubefits)
    

def interactive(musecubefits, clean=True, zlevel='median', q=0, cfwidth=100, cftype='weight', 
                pevals=[], nevals=[], optimize=False, silent=False, extSVD='',rec_settings='', 
                enhanced_optimize=False):
    """
    Performs the entire ZAP sky subtraction algorithm on an input datacube and header. A class
    containing all of the necessary data to examine the result and modify as desired.

    """

    # Check for consistency between weighted median and zlevel keywords
    if cftype == 'weight' and zlevel == 'none':
        print 'weighted median requires a zlevel calculation'
        return

    #create an instance
    zobj= zclass(musecubefits)

    if rec_settings != '':
        if rec_settings == 'filled' or rec_settings == 'extSVD':
            if extSVD == '':
                print 'Filled Field case requires external SVD'
                return
            enhanced_optimize = True
            cfwidth = 20
            cftype = 'weight'
            refilter = True
            print 'Using recommended settings for filled field case:'
            print "cfwidth = 20, cftype = 'weight', enhanced_optimization=True, refilter = True" 
        if rec_settings == 'sparse':
            print 'Using recommended settings for sparse field case:'
            print "cfwidth = 100, cftype = 'weight', enhanced_optimization=True, refilter = True" 
            enhanced_optimize = True
            cfwidth = 100
            cftype = 'weight'
            refilter = True

    if enhanced_optimize == True:
        optimize = True

    zobj._run(clean=clean, zlevel=zlevel, q=q, cfwidth=cfwidth, cftype=cftype,
              pevals=pevals, nevals=nevals, optimize=optimize,
              enhanced_optimize = enhanced_optimize,
              silent=silent, extSVD=extSVD)

    return zobj

def SVDoutput(musecubefitslst, svdfn='ZAP_SVD.fits', clean=True, zlevel='median', q=0,
              cftype='weight', cfwidth=100, mask=''):
    """
    Performs the SVD decomposition of the datacube for use in a different datacube.

    """

    if os.path.exists(svdfn):
        print 'output filename exists'
        return

    # Check for consistency between weighted median and zlevel keywords
    if cftype == 'weight' and zlevel == 'none':
        print 'weighted median requires a zlevel calculation'
        return

    if type(musecubefitslst) == list:
        print 'Combining multiple inputs'
        
        for i in range(len(musecubefitslst)):
    
            zobj= zclass(musecubefitslst[i])
            # clean up the nan values
            if clean != False:
                zobj._nanclean()
        
            #if mask is supplied, apply it
            if mask != '':
                zobj._applymask(mask[i])
        
                # Extract the spectra that we will be working with
            zobj._extract()
        
            #remove the median along the spectral axis
            if zlevel.lower() != 'none':
                zobj._zlevel(calctype=zlevel, q=q)
            
            if i == 0:
                combstack=zobj.stack
                zlstack=[zobj.zlsky]
            else:
                combstack=np.append(combstack, zobj.stack, axis=1)
                zlstack=np.append(zlstack,[zobj.zlsky], axis=0)
        zl=zlstack.mean(axis=0)

        zobj.zlsky = zl
        zobj.stack = combstack
    else:

        zobj= zclass(musecubefitslst)
        # clean up the nan values
        if clean != False:
            zobj._nanclean()
        
        #if mask is supplied, apply it
        if mask != '':
            zobj._applymask(mask)

        # Extract the spectra that we will be working with
        zobj._extract()

        #remove the median along the spectral axis
        if zlevel.lower() != 'none':
            zobj._zlevel(calctype=zlevel, q=q)

    #remove the continuum level - this is multiprocessed to speed it up
    print 'Continuum Filtering'
    zobj._continuumfilter(cftype=cftype, cfwidth=cfwidth)

    # do the multiprocessed SVD calculation
    print 'Calculating SVD'
    zobj._msvd()

    #write to file
    zobj.writeSVD(svdfn=svdfn)

def contsubfits(musecubefits, contsubfn='CONTSUB_CUBE.fits', cfwidth=100):
    """
    A multiprocessed implementation of the continuum removal. This process distributes the
    data to many processes that then reassemble the data. Uses two filters, a small scale
    (less than the line spread function) uniform filter, and a large scale median filter
    to capture the structure of a variety of continuum shapes.

    added to class
    contarray - the removed continuua
    normstack - "normalized" version of the stack with the continuua removed
    """
    
    if os.path.exists(contsubfn):
        print 'output filename "{0}" exists'.format(contsubfn)
        return

    hdu = pyfits.open(musecubefits)
    stack = hdu[1].data.reshape(hdu[1].data.shape[0],(hdu[1].data.shape[1]*hdu[1].data.shape[2]))

    contarray=_cfmedian(stack, cfwidth=cfwidth)

    #remove continuum features
    hdu[1].data = (stack - contarray).reshape(hdu[1].data.shape[0],
                                              hdu[1].data.shape[1],
                                              hdu[1].data.shape[2])
    
    hdu.writeto(contsubfn)
    

def nancleanfits(musecubefits, outfn='NANCLEAN_CUBE.fits', rejectratio=0.25, boxsz=1, silent=False):
    """
    Detects NaN values in cube and removes them by replacing them with an
    interpolation of the nearest neighbors in the data cube. The positions in the cube
    are retained in nancube for later remasking.
    """
    if os.path.exists(outfn):
        print 'output filename "{0}" exists'.format(outfn)
        return

    hdu = pyfits.open(musecubefits)
    cleancube= _nanclean(hdu[1].data, rejectratio=rejectratio, boxsz=boxsz, silent=silent)
    hdu[1].data = cleancube[0]
    hdu.writeto(outfn)


##################################################################################################
##################################### Process Steps ##############################################
##################################################################################################

class zclass:

    """
    The zclass retains all methods and attributes to run each of the steps of ZAP.

    Attributes:

      cleancube - The final datacube after removing all of the residual features.

      contarray - A 2d array containing the subtracted continuum per spaxel.

      cube - The original data cube with the zlevel subtraction performed per spaxel.

      especeval - A list containing the full set of eigenspectra and eigenvalues generated by the
                     SVD calculation that is used toy reconstruct the entire datacube.

      laxis - A 1d array containing the wavelength solution generated from the header
                     parameters.

      lparams - An array of parameters taken from the header to generate the wavelength solution.

      lranges - A list of the wavelength bin limits used in segmenting the sepctrum for SVD.

      lmin,lmax - the wavelength limits placed on the datacube

      nancube - A 3d boolean datacube containing True in voxels where a NaN value was replaced
                     with an interpolation.

      nevals - A 1d array containing the number of eigenvalues used per segment to reconstruct
                     the residuals.

      normstack - A normalized version of the datacube decunstructed into a 2d array.

      nsegments - The number of divisions in wavelength space that the cube is cut into in order
                     to perform the SVD.

      varlist - An array for each segment with the variance curve, calculated for the
                    optimize method.

      pranges - The pixel indices of the bounding regions for each spectral segment.

      recon - A 2d array containing the reconstructed emission line residuals.

      run_clean - Boolean that indicates that the NaN cleaning method was used.

      run_zlevel - Boolean indicating that the zero level correction was used.

      stack - The datacube deconstructed into a 2d array for use in the the SVD.

      subespeceval - The subset of eigenvalues and eigenspectra used to reconstruct the sky
                        residuals.

      variancearray - A list of length nsegments containing variances calculated per spaxel used
                        for normalization

      y,x - The position in the cube of the spaxels that are in the 2d deconstructed stack

      zlsky - A 1d array containing the result of the zero level subtraction



    """

    #setup the data structure
    def __init__(self, musecubefits):
        """
        Initialization of the zclass. Pulls the datacube into the class and trims it based
        on the known optimal spectral range of MUSE.
        """
        
        hdu = pyfits.open(musecubefits)
        cube = hdu[1].data
        header = hdu[1].header


        self.musecubefits = musecubefits

        self.header = header

        lparams = [header['NAXIS3'], header['CRVAL3'], header['CD3_3'], header['CRPIX3']]
        laxis=lparams[1] + lparams[2] * (np.arange(lparams[0]) + lparams[3] - 1)
        
        lmin = min(laxis)
        self.lmin = lmin
        lmax = max(laxis)
        self.lmax = lmax
        
        wlaxis = np.where(np.logical_and(laxis >= lmin, laxis <= lmax))[0]
        wlmin = min(wlaxis)
        wlmax = max(wlaxis)
        self._wlmin = wlmin
        self._wlmax = wlmax
        
        laxis = laxis[wlmin:wlmax+1]
        self.cubetrimb = cube[:wlmin,:,:] #save the trimmings
        self.cubetrimr = cube[wlmax+1:,:,:]
        
        cube = cube[wlmin:wlmax+1, :, :] #cut off the unusable bits
        self.cube = cube
        
        self.laxis = laxis
        
        #NaN Cleaning
        self.run_clean = False
        self.nancube = None
        self._boxsz = 1
        self._rejectratio = 0.25
        
        #zlevel parameters
        self.run_zlevel = False
        self.zlsky = np.zeros(laxis.shape)
        self.zlq = 0
        
        #Extraction results
        self.stack = np.array([])
        self.y = np.array([])
        self.x = np.array([])
        
        #Normalization Maps
        self.contarray = np.array([])
        self.variancearray = np.array([])
        self.normstack = np.array([])
        
        #List of segmentation limits in the optical
        skyseg=np.array([0., 4850., 5400., 5850., 6440.,
                         6750., 7200., 7575., 7700., 8265., 8602., 8731., 9275., 10000.])
        
        #identify the spectral range of the dataset
        laxmin = min(laxis)
        laxmax = max(laxis)
        
        #select the appropriate limits add a few Angstroms for padding
        lranges = np.transpose(np.vstack([
            np.append(laxmin-10,skyseg[np.logical_and(skyseg >= laxmin, 
                                                   skyseg <= laxmax)]), 
            np.append(skyseg[np.logical_and(skyseg >= laxmin, 
                                            skyseg <= laxmax)],laxmax+10)]))
        
        self.lranges=lranges
        self.nsegments = len(lranges)
        self.lparams = [header['NAXIS3'], header['CRVAL3'], header['CD3_3'], header['CRPIX3']]

        paxis=np.arange(len(laxis))
        pranges=[]
        for i in range(len(lranges)):
            lrangelogical= np.logical_and(laxis > lranges[i,0],laxis <= lranges[i,1])
            pranges.append((np.min(paxis[lrangelogical]),np.max(paxis[lrangelogical])+1))

        self.pranges=np.array(pranges)

        self.especeval = []

        #eigenspace Subset
        self.subespeceval = []

        #Reconstruction of sky features
        self.recon = np.array([])
        self.cleancube = np.array([])
        self.varlist = np.array([]) #container for variance curves
        hdu.close()


    def _run(self, clean=True, zlevel='median', q=0, cftype='weight', cfwidth=100, pevals=[], 
             nevals=[], optimize=False, silent=False, extSVD='', enhanced_optimize=False, refilter=True):

        """

        Perform all zclass to ZAP a datacube, including NaN re/masking, deconstruction into
        "stacks", zerolevel subraction, continuum removal, normalization, singular value
        decomposition, eigenvector selection, residual reconstruction and subtraction, and
        data cube reconstruction.

        Returns a "zclass" class that retains all of the data needed as the routine progresses.

        """

        t0=time()
        
        self.enhanced_optimize = enhanced_optimize
        
        print 'Preparing Data for eigenspectra calculation'
        # clean up the nan values
        if clean != False:
            self._nanclean()
        else:
            self.badcube = np.logical_not(np.isfinite(self.cube))
            badmap = self.badcube.sum(axis=0)
            self.cleancube = self.cube
            self.nancube = self.badcube
            self.spaxreject = badmap > (.99*(self.cleancube.shape)[0])
            

        # Extract the spectra that we will be working with
        self._extract()

        #remove the median along the spectral axis
        if extSVD == '':
            if zlevel.lower() != 'none':
                self._zlevel(calctype=zlevel, q=q)
        else:
            self._externalzlevel(extSVD)
            self.zlq=q

        #remove the continuum level - this is multiprocessed to speed it up
        print 'Applying Continuum Filter'
        self._continuumfilter(cfwidth=cfwidth, cftype=cftype)

        # do the multiprocessed SVD calculation
        if extSVD == '':
            self._msvd()
            #remove strong sources and use the previous SVD to help isolate sky components
            if refilter == True:
                print 'Re-applying Continuum Filter for object avoidance in eigenvalues'
                self._continuumfilter(cfwidth=20, cftype='weight')
        else:
            self._externalSVD(extSVD)
            
        
        # choose some fraction of eigenspectra or some finite number of eigenspectra
        if optimize == True or (nevals == [] and pevals == []):
            self.optimize()
        else:
            self.chooseevals(pevals=pevals, nevals=nevals)
            self.reconstruct() # reconstruct the sky residuals using the subset of eigenspace

        # stuff the new spectra back into the cube
        self.remold()

        t=time()-t0
        print 'Time to Fully processs: {0}'.format(int(t))


    ## Clean up the nan value spaxels
    def _nanclean(self, silent=False):
       """
       Detects NaN values in cube and removes them by replacing them with an
       interpolation of the nearest neighbors in the data cube. The positions in the cube
       are retained in nancube for later remasking.
       """
       boxsz = self._boxsz
       rejectratio = self._rejectratio
    
       cleancube = _nanclean(self.cube, rejectratio=rejectratio, boxsz=boxsz, silent=silent)
    
       self.run_clean = True
       self.cube = cleancube[0]
       self.nancube = cleancube[1]
       self.spaxreject = cleancube[2]

    def _extract(self, silent=False):
        """
        Deconstruct the datacube into a 2d array, since spatial information is not required, and
        the linear algebra routines require 2d arrays.
        
        The operation rejects any spaxel with even a single NaN value, since this would cause the
        linear algebra routines to crash.
        
        Adds the x and y data of these positions into the zclass
        
        """
        if silent == False:
            print 'Extracting to 2D'
        
        # make a map of spaxels with NaNs
        self.y, self.x = np.where(self.spaxreject == 0)    # get positions of acceptable spaxels
        self.stack = self.cube[:, self.y, self.x]  # extract those positions into a 2d array
        self.nanstack = self.nancube[:, self.y, self.x]
        
    def _externalzlevel(self, extSVD):
        svdhdu = pyfits.open(extSVD)
        #remove the same zero level
        print 'Using external zlevel'
        self.zlsky = svdhdu[0].data
        self.stack = self.stack - self.zlsky[:,np.newaxis]
        self.run_zlevel = 'extSVD'
        svdhdu.close()

    def _zlevel(self, calctype='median', q=0):
        """

        Removes a 'zero' level from each monochromatic layer. Spatial information 
        is not required,
        so it operates on the extracted stack.

        Operates on stack, leaving it with this level removed and adds the data 'zlsky' to the
        class. zlsky is a spectrum of the zero levels.

        This zero level is currently calculated with a median.

        Experimental operations -

        - exclude top quartile
        - run in an iterative sigma clipped mode

        """
        
        if calctype != 'none':

            print 'Subtracting Zero Level'
            t0=time()

            zlstack = np.ma.array(self.stack.copy(), mask=np.logical_not(np.isfinite(self.stack)))
            
            #choose the included quartiles
            q = int(q)
            if q > 3:
                q = 3
            self.zlq = q
            if q >= 1 and q <= 3:
                print 'Removing top {0} quartiles from zlevel calculation'.format(self.zlq)
                zlstack.sort(axis = 1)
                zlstack=zlstack[:, 0:zlstack.shape[1] * (4-self.zlq) * 0.25]

            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            jobs = []
            nseg = len(self.pranges) # choose some arbitrary number of processes


            #Do a simple median calculation
            if calctype == 'median':
                print 'Median zlevel calculation'

                #multiprocess the zlevel calculation, operating per segment
                for i in range(nseg):
                    p = multiprocessing.Process(target=_imedian, args=(i,
                                                                       zlstack,
                                                                       self.pranges,
                                                                       return_dict))
                    jobs.append(p)
                    p.start()
                
                #gather the results
                for proc in jobs:
                    proc.join()

            #Do an iterative sigma clip
            elif calctype == 'sigclip':
                print 'Iterative Sigma Clipping zlevel calculation'

                #multiprocess the zlevel calculation, operating per segment
                for i in range(nseg):
                    p = multiprocessing.Process(target=_isigclip, args=(i,
                                                                        zlstack,
                                                                        self.pranges,
                                                                        return_dict))
                    jobs.append(p)
                    p.start()

                for proc in jobs:
                    proc.join()

            zlsky=[]
            segs=return_dict.values()
            for seg in segs:
                zlsky=zlsky+seg
            self.zlsky=np.array(zlsky)

            self.stack = self.stack - self.zlsky[:,np.newaxis]
            t = time() - t0
            print 'Time to remove zero level: {0} s'.format(int(t))

        else:
            print 'Skipping zlevel subtraction'

        self.run_zlevel = calctype


    def _continuumfilter(self, cfwidth=100, cftype='weight', silent=False):

        """
        A multiprocessed implementation of the continuum removal. This process distributes the
        data to many processes that then reassemble the data. Uses two filters, a small scale
        (less than the line spread function) uniform filter, and a large scale median filter
        to capture the structure of a variety of continuum shapes.

        added to class
        contarray - the removed continuua
        normstack - "normalized" version of the stack with the continuua removed
        """

        self.stack[np.isnan(self.stack)] = 0
        
        self._cfwidth = cfwidth
        if cftype != 'weight':
            cftype = 'median'
        self._cftype = cftype
        
        
        if cftype == 'median':
            self.contarray=_cfmedian(self.stack, cfwidth=self._cfwidth)
        elif cftype == 'weight':
            self.contarray=_cfweight(self.stack, 1./self.zlsky**2, #np.abs(self.zlsky - (np.max(self.zlsky)+1)), 
                                     cfwidth=self._cfwidth)
        
        #remove continuum features
        self.normstack = self.stack - self.contarray
        
    def _msvd(self, silent=False):

        """
        Multiprocessed singular value decomposition.

        First the normstack is normalized per segment per spaxel by the variance.
        Takes the normalized, spectral segments and distributes them to the individual svd
        methods.

        """
        t0=time()
        
        #split the range
        
        nseg = len(self.pranges)
        
        #normalize the variance in the segments
        variancearray = np.zeros((nseg, self.stack.shape[1]))
        mnormstack = np.ma.array(self.normstack.copy(), mask = np.logical_not(np.isfinite(self.normstack)))
        
        for i in range(nseg):
            variancearray[i,:] = np.ma.var(mnormstack[
                self.pranges[i,0]:self.pranges[i,1], :], axis=0)
            self.normstack[self.pranges[i,0]:self.pranges[i,1], :] = self.normstack[
                self.pranges[i,0]:self.pranges[i,1], :] / variancearray[i,:]
        
        self.variancearray = variancearray
        
        print 'Beginning SVD on {0} segments'.format(nseg)
        
        #for receiving results of processes
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        
        jobs = []

        #multiprocess the svd calculation, operating per segment
        for i in range(nseg):
            p = multiprocessing.Process(target=_isvd,args=(i,
                                                           self.pranges[i],
                                                           self.normstack,
                                                           self.nanstack,
                                                           self.zlsky,
                                                           return_dict))
            jobs.append(p)
            p.start()

        #gather the results
        for proc in jobs:
            proc.join()
        
        especeval=[]
        for i in range(nseg):
            especeval.append(return_dict[i])
        
        t=time()-t0
        if silent == False:
            print 'Time to run svd : {0} s'.format(int(t))

        self.especeval = especeval

    def chooseevals(self, nevals=[], pevals=[]):
        """
        Choose the number of eigenspectra/evals to use for reconstruction

        user supplies the number of eigen spectra to be used (neval) or the percentage
        of the eigenspectra that were calculated (peval) from each spectral segment to be used.

        The user can either provide a single value to be used for all segments, or
        provide an array that defines neval or peval per segment.

        """
        if type(nevals) != list and type(nevals) != type(np.array([0])): # this can be cleaner
            nevals = [nevals]
        if type(pevals) != list and type(pevals) != type(np.array([0])):
            pevals = [pevals]
        nranges = len(self.especeval)
        nevals = np.array(nevals)
        pevals = np.array(pevals)
        nespec=[]
        for i in range(nranges):
            nespec.append((self.especeval[i][0]).shape[1])
        nespec=np.array(nespec)

        #deal with no selection
        if len(nevals) == 0 and len(pevals) == 0:
            print 'number of modes not selected'
            nevals=np.array([1])

        #deal with an input list
        if len(nevals) > 1:
            if len(nevals) != nranges:
                nevals = np.array([nevals[0]])
                print('Chosen eigenspectra array does not correspond to number of segments')
            else:
                print('Choosing {0} eigenspectra for segments'.format(nevals))

        if len(pevals) > 1:
            if len(pevals) != nranges:
                pevals = np.array([pevals[0]])
                print('Chosen eigenspectra array does not correspond to number of segments')
            else:
                print('Choosing {0}% of eigenspectra for segments'.format(pevals))
                nevals = (pevals*nespec/100.).round().astype(int)

        # deal with single value entries
        if len(pevals) == 1:
            print('Choosing {0}% of eigenspectra for all segments'.format(pevals))
            nevals = (pevals*nespec/100.).round().astype(int)
        elif len(nevals) == 1:
            print('Choosing {0} eigenspectra for all segments'.format(nevals))
            nevals = np.zeros(nranges, dtype=int) + nevals

        #take subset of the eigenspectra and put them in a list
        subespeceval=[]
        for i in range(nranges):
            eigenspectra, evals = self.especeval[i]
            tevals = (evals[0:nevals[i],:]).copy()
            teigenspectra = (eigenspectra[:,0:nevals[i]]).copy()
            subespeceval.append((teigenspectra,tevals))

        self.subespeceval = subespeceval
        self.nevals = nevals

    def reconstruct(self):

        """
        Multiprocessed residual reconstruction.

        Distributes the trimmed eigenspectra/eigenvalues to the reconstruction method.
        """

        print 'Reconstructing Sky Residuals'

        nseg=len(self.especeval)

        # take each ministack and run them independently
        pool = Pool(processes=nseg) # start pool

        #do the multiprocessing on each ministack and return eigenvalues/eigenvectors
        reconpieces = pool.map(_ireconstruct, self.subespeceval)

        pool.close()
        pool.join()

        #rescale to correct variance
        for i in range(nseg):
            reconpieces[i] = (reconpieces[i] * self.variancearray[i,:])
        self.recon=np.concatenate(reconpieces)
    
    #stuff the stack back into a cube
    def remold(self):
        """
        Subtracts the reconstructed residuals and places the cleaned spectra into the duplicated
        datacube.
        """
        print 'Applying correction and reshaping data product'
        self.cleancube = self.cube.copy()
        self.cleancube[:,self.y,self.x] = self.stack-self.recon
        #if self.run_clean == True:
        self.cleancube[self.nancube] = np.nan

    #redo the residual reconstruction with a different set of parameters
    def reprocess(self, pevals=[], nevals=[]):
        """
        A method that redoes the eigenvalue selection, reconstruction, and remolding of the
        data.
        """

        self.chooseevals(pevals=pevals, nevals=nevals)
        self.reconstruct()
        self.remold()


    def optimize(self):
        """
        Function to optimize the number of components used to characterize the residuals.

        This function calculates the variance per segment with an increasing number of
        eigenspectra/eigenvalues. It then deterimines the point at which the second derivative
        of this variance curve reaches zero. When this occurs, the linear reduction in variance
        is attributable to the removal of astronomical features rather than emission line
        residuals.

        """
        print 'Optimizing'

        t0 = time()

        nseg=len(self.especeval)

        #for receiving results of processes
        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        jobs = []

        #multiprocess the variance calculation, operating per segment
        for i in range(nseg):
            p = multiprocessing.Process(target=_ivarcurve,args=(i,
                                                                self.stack,
                                                                self.pranges,
                                                                self.especeval,
                                                                self.variancearray,
                                                                self.contarray,
                                                                self.zlsky,
                                                                return_dict))
            jobs.append(p)
            p.start()

        #gather the results
        for proc in jobs:
            proc.join()

        self.varlist = np.array(return_dict.values())
        self.nevals=np.zeros(nseg, dtype=int)

        for i in range(nseg):

            if self.enhanced_optimize == False:
                #optimize
                deriv = (np.roll(self.varlist[i],-1)-self.varlist[i])[:-1]
                deriv2 = (np.roll(deriv,-1)-deriv)[:-1]
                smderiv = ndi.uniform_filter(deriv,3)
                smderiv2 = ndi.uniform_filter(ndi.uniform_filter(np.abs(deriv2), 3), 3)
                
                noptpix=self.varlist[i].size
                
                #statistics on the derivatives
                mn1=deriv[.5 * (noptpix-2):].mean()
                std1=deriv[.5 * (noptpix-2):].std()*2
                mn2=deriv2[.5 * (noptpix-2):].mean()
                std2=deriv2[.5 * (noptpix-2):].std()*2
                #look for crossing points. When they get within 1 sigma of mean in settled region.
                cross1 = np.append([False], deriv >= (mn1 - std1)) #pad by 1 for 1st deriv
                cross2 = np.append([False, False], np.abs(deriv2) <= (mn2 + std2)) #pad by 2 for 2nd
                cross = np.logical_or(cross1,cross2)
            if self.enhanced_optimize == True:
                print 'Enhanced Optimization'

            #optimize
                deriv = (np.roll(self.varlist[i],-1)-self.varlist[i])[:-1]
                deriv2 = (np.roll(deriv,-1)-deriv)[:-1]
                smderiv = ndi.uniform_filter(deriv,3)
                smderiv2 = ndi.uniform_filter(ndi.uniform_filter(np.abs(deriv2), 3), 3)
                
                noptpix=self.varlist[i].size
                
                #statistics on the derivatives
                mn1=deriv[.75 * (noptpix-2):].mean()
                std1=deriv[.75 * (noptpix-2):].std()
                mn2=deriv2[.75 * (noptpix-2):].mean()
                std2=deriv2[.75 * (noptpix-2):].std()

                #cross = np.append([False, False], np.abs(smderiv2) <= (mn2 + std2)) #pad by 2 for 2nd
                cross = np.append([False], deriv >= (mn1 - std1)) #pad by 1 for 1st deriv

            self.nevals[i] = np.where(cross)[0][0]

        self.chooseevals(nevals=self.nevals)
        self.reconstruct()

        print 'Time to optimize: {0} s'.format(int(time()-t0))



    ##############################################################################################
    ###################################### Extra Functions #######################################
    ##############################################################################################

    def make_contcube(self):
        #remold the continuum array so it can be investigated

        """
        Takes the continuum stack and returns it into a familiar cube form.
        """
        contcube = self.cube.copy() * np.nan
        contcube[:, self.y, self.x] = self.contarray
        return contcube

    def _externalSVD(self, extSVD):

        print 'Calculating eigenvalues for input eigenspectra'
        svdhdu = pyfits.open(extSVD)
        nseg = len(self.pranges)

        #normalize the variance in the segments
        self.variancearray = np.zeros((nseg, self.stack.shape[1]))

        for i in range(nseg):
            self.variancearray[i,:] = np.var(self.normstack[
                self.pranges[i,0]:self.pranges[i,1], :], axis=0)
            self.normstack[self.pranges[i,0]:self.pranges[i,1], :] = self.normstack[
                self.pranges[i,0]:self.pranges[i,1], :] / self.variancearray[i,:]

        especeval=[]
        for i in range(nseg):
            eigenspectra = svdhdu[i+1].data
            ns = self.normstack[self.pranges[i][0]:self.pranges[i][1]]
            evals = np.transpose(np.transpose(ns).dot(eigenspectra))
            especeval.append([eigenspectra,evals])

        self.especeval = especeval

    def _recalc_evals(self, extSVD):

        print 'Recalculating eigenvalues for strong line object removal'

        #normalize the variance in the segments
        self.variancearray = np.zeros((nseg, self.stack.shape[1]))

        for i in range(nseg):
            self.variancearray[i,:] = np.var(self.normstack[
                self.pranges[i,0]:self.pranges[i,1], :], axis=0)
            self.normstack[self.pranges[i,0]:self.pranges[i,1], :] = self.normstack[
                self.pranges[i,0]:self.pranges[i,1], :] / self.variancearray[i,:]

        especeval=[]
        for i in range(nseg):
            #eigenspectra = svdhdu[i+1].data
            ns = self.normstack[self.pranges[i][0]:self.pranges[i][1]]
            evals = np.transpose(np.transpose(ns).dot(especeval[i][0]))
            especeval.append([eigenspectra,evals])

        self.especeval = especeval



    #apply a mask to the input data in order to provide a cleaner basis set
    def _applymask(self, mask):
        print 'Applying Mask for SVD Calculation file {0}'.format(mask)
        #mask is >1 for objects, 0 for sky so that people can use sextractor
        hmsk=pyfits.open(mask)
        bmask = hmsk[0].data.astype(bool)
        self.cube[:, bmask] = np.nan #self.cube / bmask[np.newaxis,:,:]

    ##############################################################################################
    ##################################### Output Functions #######################################
    ##############################################################################################
    def _cubetowrite(self):
        return np.concatenate((self.cubetrimb,self.cleancube,self.cubetrimr), axis=0)

    def _skycubetowrite(self):
        return np.concatenate((self.cubetrimb,self.cube-self.cleancube,self.cubetrimr), axis=0)

    def writecube(self, outcubefits='DATACUBE_ZAP.fits'):
        """
        write the processed datacube to an individual fits file.
        """

        if os.path.exists(outcubefits):
            print 'output filename exists'
            return

        #fix up for writing
        outcube = self._cubetowrite()
        outhead = _newheader(self)

        #create hdu and write
        outhdu = pyfits.PrimaryHDU(data=outcube,header=outhead)
        outhdu.writeto(outcubefits)

    def writeskycube(self, skycubefits='SKYCUBE_ZAP.fits'):
        """
        write the processed datacube to an individual fits file.
        """

        if os.path.exists(skycubefits):
            print 'output filename "{0}" exists'.format(outcubefits)
            return

        #fix up for writing
        outcube = self._skycubetowrite()
        outhead = _newheader(self)

        #create hdu and write
        outhdu = pyfits.PrimaryHDU(data=outcube,header=outhead)
        outhdu.writeto(skycubefits)
    
    def mergefits(self, outcubefits):
        """
        Merge the ZAP cube into the full muse datacube and write
        """

        outcubefits=outcubefits.split('.fits')[0]+'.fits' #make sure it has the right extension

        if os.path.exists(outcubefits):
            print 'output filename "{0}" exists'.format(outcubefits)
            return

        hdu = pyfits.open(self.musecubefits)
        hdu[1].header = _newheader(self)
        hdu[1].data = self._cubetowrite()

        hdu.writeto(outcubefits)

    def writeSVD(self, svdfn='ZAP_SVD.fits'):

        if os.path.exists(svdfn):
            print 'output filename "{0}" exists'.format(svdfn)
            return

        hdu = pyfits.HDUList()

        hdu.append(pyfits.PrimaryHDU(self.zlsky))

        for i in range(len(self.pranges)):
            hdu.append(pyfits.ImageHDU(self.especeval[i][0]))

        #write for later use
        hdu.writeto(svdfn)

    def plotvarcurve(self, i=0):
    
        if len(self.varlist) == 0:
            print 'No varlist found. The optimize method must be run first. \n'
            return

        if self.enhanced_optimize == False:        
            #optimize
            deriv = (np.roll(self.varlist[i],-1)-self.varlist[i])[:-1]
            deriv2 = (np.roll(deriv,-1)-deriv)[:-1]
            
            noptpix=self.varlist[i].size
            
            #statistics on the derivatives
            mn1=deriv[.5 * (noptpix-2):].mean()
            std1=deriv[.5 * (noptpix-2):].std()*2
            mn2=deriv2[.5 * (noptpix-2):].mean()
            std2=deriv2[.5 * (noptpix-2):].std()*2
        if self.enhanced_optimize == True:
            
            #optimize
            deriv = (np.roll(self.varlist[i],-1)-self.varlist[i])[:-1]
            deriv2 = (np.roll(deriv,-1)-deriv)[:-1]
            smderiv = ndi.uniform_filter(deriv,3)
            smderiv2 = ndi.uniform_filter(ndi.uniform_filter(np.abs(deriv2), 3), 3)
            
            noptpix=self.varlist[i].size
            
            #statistics on the derivatives
            mn1=deriv[.75 * (noptpix-2):].mean()
            std1=deriv[.75 * (noptpix-2):].std()
            mn2=deriv2[.75 * (noptpix-2):].mean()
            std2=deriv2[.75 * (noptpix-2):].std()
            
        fig=plt.figure(figsize=[10,15])
        ax= fig.add_subplot(3,1,1)
        plt.plot(self.varlist[i],linewidth=3)
        plt.plot([self.nevals[i],self.nevals[i]],[min(self.varlist[i]),max(self.varlist[i])])
        plt.ylabel('Variance')
        
        ax= fig.add_subplot(3,1,2)
        plt.plot(np.arange(deriv.size),deriv)
        plt.plot([0,len(deriv)],[mn1,mn1],'k')
        plt.plot([0,len(deriv)],[mn1-std1,mn1-std1],'0.5')
        plt.plot([self.nevals[i]-1,self.nevals[i]-1],[min(deriv),max(deriv)])
        plt.ylabel('d/dn Variance')

        ax= fig.add_subplot(3,1,3)
        plt.plot(np.arange(deriv2.size),np.abs(deriv2))
        plt.plot([0,len(deriv2)],[mn2,mn2],'k')
        plt.plot([0,len(deriv2)],[mn2+std2,mn2+std2],'0.5')
        plt.plot([self.nevals[i]-2,self.nevals[i]-2],[min(deriv2),max(deriv2)])
        plt.ylabel('(d^2/dn^2) Variance')
        plt.xlabel('Number of Components')

        plt.suptitle('Segment {0}, {1} - {2} Angstroms'.format(i, 
                                                               self.lranges[i][0],
                                                               self.lranges[i][1]))
        
        

##################################################################################################
##################################### Helper Functions ###########################################
##################################################################################################


##### Continuum Filtering #####
def rolling_window(a, window): #function for striding to help speed up
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def _icfweight(i, stack, wt, cfwidth, sprange, return_dict):
    istack = np.rollaxis(stack[:,sprange[0]:sprange[1]],1)
    result = np.rollaxis(np.array([wmedian(row, wt, cfwidth=cfwidth) for row in istack]),1)
    result = ndi.uniform_filter1d(result, 3, axis=0)
    return_dict[i] = result

def wmedian(spec, wt, cfwidth=100):

    # ignore the warning (feature not a bug)
    old_settings = np.seterr(divide='ignore') 
    spec = np.lib.pad(spec,(cfwidth,cfwidth),'constant', constant_values=0)
    wt = np.lib.pad(wt,(cfwidth,cfwidth),'constant', 
                  constant_values=(np.min(wt)/1000.,np.min(wt)/1000.))

    #do some striding for speed
    swin = rolling_window(spec, cfwidth) # create window container array
    wwin = rolling_window(wt, cfwidth) # create window container array
    
    #sort based on data
    srt = np.argsort(swin, axis=-1)
    ind = np.ogrid[0:swin.shape[0], 0:swin.shape[1]]
    sdata = swin[ind[0], srt]
    swt = wwin[ind[0], srt]
    
    #calculate accumulated weights
    awt = np.cumsum(swt,axis=-1)
    
    # new weightsort for normalization and consideration of data
    nw = (awt-0.5*swt)/awt[:,-1][:,np.newaxis]
    
    #find the midpoint in the new weight sort
    s = np.argmin(np.abs(nw-0.5),axis=-1)
    sl=np.arange(len(s))
    nws=nw[sl,s]
    nws1=nw[sl,s-1]

    f1 = (nws - 0.5)/(nws - nws1) 
    f2 = (0.5 - nws1)/(nws - nws1) 
    wmed = sdata[sl,s-1]*f1 + sdata[sl,s]*f2
    wmed = wmed[cfwidth/2:-cfwidth/2-1]
    np.seterr(old_settings['divide'])

    return wmed

def _cfweight(stack, weight, cfwidth=100, silent=False):

    """
    A multiprocessed implementation of the continuum removal. This process distributes the
    data to many processes that then reassemble the data. Uses two filters, a small scale
    (less than the line spread function) uniform filter, and a large scale median filter
    to capture the structure of a variety of continuum shapes.

    added to class
    contarray - the removed continuua
    normstack - "normalized" version of the stack with the continuua removed
    """

    t0=time()
    if silent != False:
        print 'Continuum Subtracting - weighted median filter method'
    nmedpieces=multiprocessing.cpu_count()

    #define bins

    edges = np.append(np.floor(
        np.arange(0, stack.shape[1], stack.shape[1]/np.float(nmedpieces))),
                      stack.shape[1])

    medianranges=np.array(zip(edges[0:-1],edges[1::])).astype(int)

    #for receiving results of processes
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    jobs = []

    #multiprocess the variance calculation, operating per segment
    for i in range(len(medianranges)):
        p = multiprocessing.Process(target=_icfweight, args=(i,
                                                             stack,
                                                             weight,
                                                             cfwidth,
                                                             medianranges[i],
                                                             return_dict))
        jobs.append(p)
        p.start()

    #gather the results
    for proc in jobs:
        proc.join()

    contarray=np.concatenate(return_dict, axis=1)

    if silent == False:
        print 'Time to Subtract : {0} s'.format(int(time()-t0))

    return contarray


def _icfmedian(i, stack, cfwidth, sprange, return_dict):  #for distributing data to the pool
    """
    Helper function to distribute data to Pool for SVD
    """
    ufilt=3 #set this to help with extreme over/under corrections
    result = ndi.median_filter(
        ndi.uniform_filter1d(stack[:,sprange[0]:sprange[1]], ufilt, axis=0), (cfwidth,1))
    return_dict[i]=result

def _cfmedian(stack, cfwidth=100, silent=False):

    """
    A multiprocessed implementation of the continuum removal. This process distributes the
    data to many processes that then reassemble the data. Uses two filters, a small scale
    (less than the line spread function) uniform filter, and a large scale median filter
    to capture the structure of a variety of continuum shapes.

    added to class
    contarray - the removed continuua
    normstack - "normalized" version of the stack with the continuua removed
    """

    t0=time()
    if silent != False:
        print 'Continuum Subtracting - median filter method'
    nmedpieces=multiprocessing.cpu_count()

    #define bins

    edges = np.append(np.floor(
        np.arange(0, stack.shape[1], stack.shape[1]/np.float(nmedpieces))),
                      stack.shape[1])

    medianranges=np.array(zip(edges[0:-1],edges[1::])).astype(int)

    #for receiving results of processes
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    jobs = []

    #multiprocess the variance calculation, operating per segment
    for i in range(len(medianranges)):
        p = multiprocessing.Process(target=_icfmedian, args=(i,
                                                             stack,
                                                             cfwidth,
                                                             medianranges[i],
                                                             return_dict))
        
        jobs.append(p)
        p.start()
    
    #gather the results
    for proc in jobs:
        proc.join()

    contarray = np.concatenate(return_dict, axis=1)

    if silent == False:
        print 'Time to Subtract : {0} s'.format(int(time()-t0))

    return contarray


##### SVD #####
def _isvd(i, prange, normstack, nanstack, zlsky, return_dict): #for distributing data to the pool
    """
    Perform single value decomposition and Calculate PC amplitudes (projection)
    outputs are eigenspectra operates on a 2D array.
    
    eigenspectra = [nbins, naxes]
    evals = [naxes, nobj]
    data = [nbins, nobj]
    """
    
    inormstack = normstack[prange[0]:prange[1],:]
    inanstack = nanstack[prange[0]:prange[1],:] #new
    inormstack[inanstack] = 0 #new
    inormstack[np.logical_not(np.isfinite(inormstack))] = 0 #new    
    inormstack = np.transpose(inormstack)
    izlsky=zlsky[prange[0]:prange[1]]
    
    #new
    #inormstack = sparse.dok_matrix(inormstack) #new    
    U,s,V = sparse.linalg.svds(inormstack, k=prange[1]-prange[0]-1) #10) #new
    eigenspectra = np.fliplr(np.transpose(V)) #put strongest first
        
    evals = (inormstack).dot(eigenspectra)
    evals = np.transpose(evals)
    
    print 'Finished SVD Segment'
    
    return_dict[i]= [eigenspectra, evals]


##### RECONSTRUCTION  #####
def _ireconstruct(iespeceval):
    """
    Reconstruct the residuals from a given set of eigenspectra and eigenvalues
    """

    eigenspectra, evals = iespeceval
    nrows=(evals.shape)[1]
    reconpiece=np.zeros([(eigenspectra.shape)[0],nrows]) # make container
    for i in np.arange(nrows): # this loop is FASTER than a fully vectorized one-liner command
        evalvect=(evals[:,i]) # choose single eval set
        reconpiece[:,i]=np.sum(eigenspectra * evalvect,axis=1) # broadcast evals on evects and sum

    return reconpiece

def _ipreconstruct(i, iespeceval, precon):
    """
    Reconstruct the residuals from a given set of eigenspectra and eigenvalues.

    this is a special version for caculating the variance curve. It adds the contribution of a
    single mode to an existing reconstruction.

    """

    eigenspectra, evals = iespeceval
    eigenspectra = eigenspectra[:,i]
    evals = evals[i,:]
    reconpiece = precon + (eigenspectra[:,np.newaxis] * evals[np.newaxis,:]) # broadcast evals on evects and sum

    return reconpiece

##### OPTIMIZE #####
def _ivarcurve(i, stack, pranges, especeval, variancearray,contarray, zlsky, return_dict):

    #segment the data
    istack = (stack-contarray)[pranges[i,0]:pranges[i,1],:]
    izlsky = zlsky[pranges[i,0]:pranges[i,1]][:,np.newaxis]*np.ones_like(istack)
    iprecon = np.zeros_like(istack)
    iespeceval = especeval[i]
    ivariancearray = variancearray[i]
    ivarlist = []
    totalnevals = int(np.round((iespeceval[1].shape[0])*0.25))

    for nevals in range(totalnevals):
        if nevals % (totalnevals * .1) <= 1:
            print 'Seg {0}: {1}% complete '.format(i, int(nevals/(totalnevals-1.)*100.))
        iprecon = _ipreconstruct(nevals, iespeceval, iprecon)
        icleanstack = istack - (iprecon * ivariancearray)
        av = np.average(icleanstack, weights=(izlsky-np.min(izlsky))**2)
        var = np.average((icleanstack - av)**2, weights=(izlsky-np.min(izlsky))**2)
        #ivarlist.append(np.var(icleanstack)) #calculate the variance on the cleaned segment
        ivarlist.append(var)
    
    ivarlist = np.array(ivarlist)
    return_dict[i] = ivarlist

def _newheader(zobj):

    header=zobj.header.copy()

    #put the pertinent zap parameters into the header
    header['COMMENT']='These data have been ZAPped!'

    # zlevel removal performed
    header.append(('ZAPzlvl',zobj.run_zlevel, 'ZAP zero level correction'), end=True)

    # zlevel removal performed
    header.append(('ZAPzlq',zobj.zlq, 'ZAP quartiles used for zero level correction'))

    # Nanclean performed
    header['ZAPclean'] = (zobj.run_clean, 'ZAP NaN cleaning performed for calculation')

    # Continuum Filtering
    header['ZAPcftyp'] = (zobj._cftype, 'ZAP continuum filter type')
    header['ZAPcfwid'] = (zobj._cfwidth, 'ZAP continuum filter size')

    # number of segments
    nseg=len(zobj.pranges)
    header['ZAPnseg'] = (nseg, 'Number of segments used for ZAP SVD')

    # per segment variables
    for i in range(nseg):
        header['ZAPseg{0}'.format(i)] = ('{0}:{1}'.format(zobj._wlmin+zobj.pranges[i][0],
                                                          zobj._wlmin+zobj.pranges[i][1]-1),
                                         'spectrum segment (pixels)')

        header['ZAPnev{0}'.format(i)] = (zobj.nevals[i], 'number of eigenvals/spectra used')

    return header

def _isigclip(i, stack, pranges, return_dict):

    istack = stack[pranges[i,0]:pranges[i,1],:]
    mn  = []

    for col in istack:
        col = col[np.isfinite(col)]
        clipped,bot,top=sigmaclip(col, low=3, high=3)
        mn.append(clipped.mean())

    return_dict[i]=mn

def _imedian(i, stack, pranges, return_dict):

    istack = stack[pranges[i,0]:pranges[i,1],:]

    medn = list(np.ma.median(istack, axis=1))

    return_dict[i]=medn


## Clean up the nan value spaxels
def _nanclean(cube, rejectratio=0.25, boxsz=1, silent=False):
    """
    Detects NaN values in cube and removes them by replacing them with an
    interpolation of the nearest neighbors in the data cube. The positions in the cube
    are retained in nancube for later remasking.
    """
    t0=time()
    
    cleancube = cube.copy()                #
    badcube = np.logical_not(np.isfinite(cleancube))        # find NaNs
    badmap = badcube.sum(axis=0)  # map of total nans in a spaxel
    spaxreject = badmap > (rejectratio*(cleancube.shape)[0])
    badcube = np.logical_or(badcube, spaxreject[np.newaxis,...])

# choose some maximum number of bad pixels in the spaxel and extract positions
    #y, x = np.where(badmap > (rejectratio*(cleancube.shape)[0]))
    nbadspax = spaxreject.sum()
    print "{0} spaxels rejected: > {1}% NaN pixels".format(nbadspax, rejectratio * 100)
    #y, x = np.where(badmap)
    
# make cube mask of bad spaxels
    #bcube = np.ones(cleancube.shape, dtype=bool)
    #bcube[:,y,x] = False

    #badcube=np.logical_and(badcube == True, bcube == True) # combine masking
    z, y, x = np.where(np.logical_xor(badcube, spaxreject[np.newaxis,...]))

    neighbor=np.zeros((z.size,(2*boxsz+1)**3))
    icounter=0
    print "fixing {0} pixels".format(len(z))

#loop over samplecubes
    for j in range(-boxsz,boxsz+1,1):
        for k in range(-boxsz,boxsz+1,1):
            for l in range(-boxsz,boxsz+1,1):
                iz, iy, ix = z+l, y+k, x+j
                outx=np.logical_or(ix <= 0, ix >= (cleancube.shape)[2]-1)
                outy=np.logical_or(iy <= 0, iy >= (cleancube.shape)[1]-1)
                outz=np.logical_or(iz <= 0, iz >= (cleancube.shape)[0]-1)
                outsider=np.where(np.logical_or(np.logical_or(outx,outy), outz))
                ix[ix < 0] = 0 ; ix[ix > (cleancube.shape)[2]-1] = (cleancube.shape)[2]-1
                iy[iy < 0] = 0 ; iy[iy > (cleancube.shape)[1]-1] = (cleancube.shape)[1]-1
                iz[iz < 0] = 0 ; iz[iz > (cleancube.shape)[0]-1] = (cleancube.shape)[0]-1
                neighbor[:,icounter]=cleancube[iz,iy,ix]
                neighbor[outsider,icounter] = np.nan
                icounter=icounter+1
    goodneighbor = np.isfinite(neighbor)
    neighbor[np.logical_not(goodneighbor)]=0
    nfix=goodneighbor.sum(axis=1)
    tfix= neighbor.sum(axis=1)
    neighborless = nfix == 0
    fix=np.zeros(z.size)
    fix[neighborless]=np.nan
    fix[np.logical_not(neighborless)] = \
        tfix[np.logical_not(neighborless)]/nfix[np.logical_not(neighborless)]
    cleancube[z,y,x]=fix

    t=time()-t0
    if silent == False:
        print 'Time to clean NaNs: {0} s'.format(int(t))

    return cleancube, badcube, spaxreject

