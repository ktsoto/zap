.. zap documentation master file, created by
   sphinx-quickstart on Mon Nov 25 09:46:49 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ZAP's documentation!
===============================

.. toctree::
   :maxdepth: 2

ZAP (the Zurich Atmosphere Purge) is a high precision sky subtraction tool which can be used as complete sky subtraction solution, or as an enhancement to previously sky-subtracted data.  **Currently, the best results come from applying ZAP to a datacube with no initial sky subtraction.** The method uses PCA to isolate the residual sky subtraction features and remove them from the observed datacube. Though the operation of ZAP is not dependent on perfect flatfielding of the data in a MUSE exposure, better results are obtained when these corrections are made ahead of time.


Examples
========

In its most hands-off form, ZAP can take an input fits datacube, operate on it, and output a final fits datacube.::

  from mpdaf_user import zap
  zap.process('INPUT.fits', 'OUTPUT.fits')

Care should be taken, however, since this case assumes a sparse field.

There are a number of options that can be passed to the code which we tabulate here, and describe in several different use cases.

===========  ===================  ==============================================================================
Keyword      Default              Function
===========  ===================  ==============================================================================
outcubefits  'DATACUBE_ZAP.fits'  the fits filename for the output datacube.

clean        True                 Interpolates over NaN values in the datacube to allow processing on all
                                  spectra. The NaN values are replaced in the final datacube.
                                  Any spaxel that includes a NaN pixel will hinder the calculation, so this
				  step is used to maximize the number of contributors. The NaN values are
				  reinserted into the final datacube.
	     	                  
zlevel       'median'             This option is used to define the method for determining the zeroth order
                                  subtraction of the sky. This is used to remove any systematic sky feature
                                  that exists over the whole field. Options are 'median', 'sigclip', and 'none'.
                                  The 'none' option should only be applied when ZAP is applied to enhance a
                                  previous sky subtraction.
	     	                  
cftype       'weight'             The type of filtering that is applied to remove the continuum.  'weight'
                                  refers to the weighted median, which uses the calculated zlevel sky as the
                                  weight. 'median' refers to a rolling median filter with a nested small
                                  uniform filter. The 'weight' option provides a better result, but is much
                                  slower (an additional 10 minutes on a single exposure) and should only be
                                  run in the complete sky subtraction case.
	     	                  
cfwidth      100                  Size of the filterbox used to remove the continuum features in order to
                                  sterilize the basis set used to calculate the eigenbasis.
	     	                  
optimize     'normal'             A flag used to call the optimization method. The possible options are
                                  'normal', 'enhanced', and 'none'.
	     	                  
nevals       []                   This option is used for the manual selection of the the number of eigenvectors
                                  to be used. If this is used, the pevals is ignored. Provide either a single
				  value that will be used for all of the segments, or a list of 11 values that
				  will be used for each of the segments.
	     	                  
pevals       []                   This option is used for the manual selection of the the number of eigenvectors
                                  to be used. This value is the percentage of the calculated
				  eigenspectra/eigenvalues used per spectral segment. Provide either a single
				  value that will be used for all of the segments, or a list of 11 values that
				  will be used for each of the segments.
	     	                  
extSVD       ''                   An optional parameter that allows the input of a externally calculated
                                  eigenbasis as well as a zlevel. This can be constructed from either a masked
				  version of a sparse field case, or an external sky frame.
	     	                  
mask         ''                   A 2D fits image to exclude regions that may contaminate the zlevel or
                                  eigenspectra. This image should be constructed from the datacube itself to
				  match the dimensionality. Sky regions should be marked as 0, and astronomical
				  sources should be identified with an integer greater than or equal to 1.

interactive  False                Setting this option to True will allow the user to pass out the zclass which
                                  contains all of the data and methods of ZAP. We describe this use below.
	     
===========  ===================  ==============================================================================

The code can handle datacubes trimmed in wavelength space. Since the code uses the correlation of segments of the emission line spectrum, it is best to trim the cube at specific wavelengths. The cube can include any connected subset of these segments. (for example 6400 - 8200 Angstroms) ::

  [4500, 5400]
  [5400, 5850]
  [5850, 6440]
  [6440, 6750]
  [6750, 7200]
  [7200, 7700]
  [7700, 8265]
  [8265, 8602]
  [8602, 8731]
  [8731, 9275]
  [9275, 9500]


Sparse Field Case
-----------------

As noted above, this case can be handled simply with the observed datacube itself, using: ::

  zap.process('INPUT.fits', 'OUTPUT.fits')

It is possible to add a mask to this operation: ::

  zap.process('INPUT.fits', 'OUTPUT.fits', mask='mask.fits')


Masked Processing
-----------------

Another option is to use a mask to isolate a sky within an exposure to pre-determine the zlevel and eigenspectra, which is then passed back into zap. This approach will allow the inclusion of a mask file, which is a 2d fits image matching the spatial dimensions of the input datacube. The values in the mask image will be 0 in the masked regions (such a where an extended object is) and 1 in the unmasked regions. Set this parameter with ``mask='maskfile.fits'`` ::

  zap.SVDoutput('INPUT.fits', svdfn='ZAP_SVD.fits', mask='maskfile.fits') # create SVD file
  zap.process('INPUT.fits', 'OUTPUT.fits', extSVD='ZAP_SVD.fits', cfwidth=50)

Filled Field Case
-----------------

This approach also can address the saturated field case and is robust in the case of strong emission lines, in this case the input is an offset sky observation. ::

  zap.SVDoutput('Offset_Field_CUBE.fits', svdfn='Offset_ZAP_SVD.fits')
  zap.process('INPUT.fits', 'OUTPUT.fits', extSVD='ZAP_SVD.fits', cfwidth=50)

The integration time of this frame does not need to be the same as the object exposure, but rather just a 2-3 minute exposure.

===================
Top Level Functions
===================

Aside from the main "full process", two functions are included that can be run outside of the entire zap process to facilitate some investigations.

**nan cleaning**

This function replaces the nan valued pixels with an average of the adjacent valid pixels. It can be called as below::

  zap.nancleanfits('INPUT.fits', outfn='NANCLEAN_CUBE.fits', rejectratio=0.25, boxsz=1)

"rejectratio" defines a cutoff for the ratio of pixels in a spaxel before the spaxel is avoided completely.

"boxsz" defines the number of pixels that defines the box around the offending nan pixel. With boxsz set to 1 the function looks for the nearest 26 neighbors which is a 3x3x3 cube.

This step is an intermediary step in the full ZAP process, but this function allows you to run it as a standalone step.

**continuum removal**

This function applies a filter on the datacube that removes most continuum features. This function allows for the enhancement of emission line characteristics. The filtering method has been enhanced by multiprocessing to produce a rapid result. It can be called as below::

  zap.contsubfits(musecubefits, contsubfn='CONTSUB_CUBE.fits', cfilter=100):

It applies a nested set of filters, one uniform of width 3 pixels and one median with a width defined by cfilter. Since it does not calculate a zlevel, it can only use the 'median' method.

================
Interactive mode
================

ZAP can also  be used interactively from within ipython using pyfits. ::

  from mpdaf_user import zap
  zclass = zap.process('INPUT.fits', interactive=True)

The run method operates on the datacube, and retains all of the data and methods necessary to
process a final data cube in a python class named zclass. You can elect to investigate the data product via the zclass, and even reprocess the cube with a different number of eigenspectra per region.  A workflow may go as follows: ::

  from mpdaf_user import zap
  from matplotlib import pyplot as plt

  zobj = zap.process('INPUT.fits', optimization='normal') #allow ZAP to run the optimize routine
  zobj.plotvarcurve(5) #plot the variance curves and the selection of the number of eigenspectra used

  #plot a spectrum extracted from the original cube
  plt.figure()
  plt.plot(zobj.cube[:,50:100,50:100].sum(axis=(1,2)), 'b', alpha=0.3)

  #plot a spectrum of the cleaned ZAP dataproduct
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2)), 'g')

  #choose just the first 3 spectra for all segmments 
  zobj.reprocess(nevals=3)

  #plot a spectrum extracted from the original cube
  plt.plot(zobj.cube[:,50:100,50:100].sum(axis=(1,2)), 'b', alpha=0.3)

  #plot a spectrum of the cleaned ZAP dataproduct
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2))), 'g')

  #choose some number of modes by hand
  zobj.reprocess(nevals=[2,5,2,4,6,7,9,8,5])

  #plot a spectrum
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2))), 'k')

  #Use the optimization algorithm to identify the best number of modes per segment
  zobj.optimize()

  #compare to the previous versions
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2))), 'r')

  #identify a pixel in the dispersion axis that shows a residual feature in the original
  plt.figure()
  plt.matshow(zobj.cube[2903,:,:])

  #compare this to the zap dataproduct
  plt.figure()
  plt.matshow(zobj.cleancube[2903,:,:])

  #write the processed cube as a single extension fits
  zobj.writecube('DATACUBE_ZAP.fits')

  #or merge the zap datacube into the original input datacube, replacing the data extension
  zobj.writefits(outcubefits='DATACUBE_FINAL_ZAP.fits')

======
ZCLASS
======
.. autoclass:: zap.zclass
   :members:

