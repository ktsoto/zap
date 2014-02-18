.. zap documentation master file, created by
   sphinx-quickstart on Mon Nov 25 09:46:49 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ZAP's documentation!
===============================

.. toctree::
   :maxdepth: 2

ZAP (the Zurich Atmosphere Purge) is a sky subtraction tool which removes residual sky emission line features.
The tool employs multiprocessing to operate quickly and robustly on a large MUSE datacube.

Examples
========

In its most hands-off form, ZAP can take an input fits datacube, operate on it, and output a
final fits datacube.::

  zap.process('INPUT.fits', 'OUTPUT.fits')


This function will modify the data and header in the first extension of the fits file

However, ZAP can also  be used interactively from within ipython using pyfits. ::

  import zap
  zclass = zap.interactive('INPUT.fits')

The run method operates on the datacube, and retains all of the data and methods necessary to
process a final data cube in a python class named zclass.

In both cases, the options are passed via the ipython command line via keywords:

========  =======  ==============================================================================
Keyword   Default  Function
========  =======  ==============================================================================
clean     True     Removes NaN values from the datacube to allow processing on allspectra.
                   The NaN values are replaced in the final datacube.

zlevel    True     Subtracts a systematic offset level from the data by calculating a
                   median per spectral channel.

cfilter   100      Size of the filterbox used to remove the continuum features in order to
                   sterilize the basis set used to calculate the eigenbasis.

nevals    False    Number of eigenspectra/eigenvaules used per spectral segment. If this
                   is used, the pevals is ignored. Provide either a single value that will be 
		   used for all of the segments, or a list of 9 values that will be used for 
		   each of the segments.

pevals    False    Percentage of the caclulated eigenspectra/eigenvaules used per spectral
                   segment. Provide either a single value that will be used for all of the 
		   segments, or a list of 9 values that will be used for each of the segments.

optimize  True     A flag used to call the optimization method. This automatically determines the
                   number of eigenspectra/eigenvalues to use per segment. 

========  =======  ==============================================================================

After running ZAP interactively, the user can elect to investigate the data product via the zclass, and even reprocess the cube with a different number of eigenspectra per region.  A workflow may go as follows: ::

  import zap
  from matplotlib import pyplot as plt
  
  zclass = zap.interactive('INPUT.fits', pevals=1) #choose 1% of modes per segment
  
  #investigate the dataproduct with pyplot
  plt.figure()
  
  #plot a spectrum extracted from the original cube
  plt.plot(zclass.cube[:,50:100,50:100].sum(axis=-1).sum(axis=-1), 'b', alpha=0.3)
  
  #plot a spectrum of the cleaned ZAP dataproduct
  plt.plot(zclass.cleancube[:,50:100,50:100].sum(axis=-1).sum(axis=-1), 'g')
  
  #choose some number of modes by hand
  zclass.reprocess(nevals=[2,5,2,4,6,7,9,8,5])

  #plot a spectrum
  plt.plot(zclass.cleancube[:,50:100,50:100].sum(axis=-1).sum(axis=-1), 'k')

  #Use the optimization algorithm to identify the number of modes per segment
  zclass.optimize()

  #compare to the previous versions
  plt.plot(zclass.cleancube[:,50:100,50:100].sum(axis=-1).sum(axis=-1), 'r')  

  #identify a pixel in the dispersion axis that shows a strong residual feature in the original
  plt.figure()
  plt.matshow(zclass.cube[,:,:])

  plt.figure()
  plt.matshow(zclass.cleancube[,:,:])

  #write the processed cube
  zclass.writecube('DATACUBE_ZAP.fits')

  #or merge the zap datacube into to whole inout datacube, replacing the data extension
  zclass.writecube('INPUT.fits','DATACUBE_FINAL_ZAP.fits')


ZCLASS
======
.. autoclass:: zap.zclass
   :members:


Algorithm description
=============================
ZAP is designed to remove the residual sky emission features after an initial Sky subtraction. To make this correction, several steps must be performed. Many of these are constructed to create a set of input spectra for singular value decomposition. By performing these prior steps we can better isolate the eigenspectra and eigenvalues that reconstruct the sky residuals without influencing the Astronomical objects. Below is a description of each of the processing steps.

The key point in this approach is that the preliminary steps are use to remove as many easily determined non-sky features from the spectra going into the SVD calculation. The sky features are then reconstructed and removed from the original datacube.


NaN cleaning
++++++++++++
Set by keyword "clean = True"

NaN (Not a Number) pixels interfere with the SVD (Singular Value Decomposition) by being unbounded. Since the eigenspectra can not reconstruct these values, the pixel must be either replaced with a float or the spaxel for which it is a member removed from the calculation. By Performing the NaN cleaning step, we choose the first of these options.

This algorithm identifies NaN pixels in the datacube and replaces them by interpolating with the 26 (3x3x3 - 1) nearest neighbor pixels in the 3D datacube. Any NaN values in this set of pixels is ingored and the mean of the finite values replaces the central pixel. ZAP retains the position of these NaN values so that the pixels can be converted back into NaN values in the final dataproduct.

The user has the option to exlude this step (run_clean = False), but should be advised that the presence of a single NaN pixel will exclude an entire spaxel from the entire calculation leaving it uncorrected.


Zero Level Subtraction
++++++++++++++++++++++
Remove the continuum from each spaxel via filtering.

Set by keyword "zlevel = True"

This processing step is performed to remove residual sky features that are consistent over the entire field. The "zero level" is determined from the median per spectral channel, producing a spectrum of the residual. This spectrum can be accessed in the interactive mode from the produced instance of the zclass (Described below)

.. warning::
   The user should be cautious in cases where a spectral feature, such as an emission line, covers the entire field. In these cases, the median calculation will erroneously remove this feature.  Several approaches are in development for handling this case.

Continuum Filtering
+++++++++++++++++++
Remove the continuum from each spaxel via filtering.

Filter size set by keyword "cfilt = 100"

The continuum filter is a process that removes astrophysical continuum features through the use of a combination of two sprectral filters. A uniform filter that is rougly on the scale of the spectral resolution (3 pixels) smooths any extreme variations on this small scale. The next filter is a large scale (100 pixels) median filter, which has the property of tracing a multitude of continuum shapes, including sharp edges in the spectrum. The result of these filters is then subtracted from the spectra leaving only object emission lines and sky residuals in the spectra.

Spectral Segmentation
+++++++++++++++++++++
This algorithm operates by segmenting the data into regions of the spectrum with strongly correlated features. This segmentation benefits the calculation in two ways. First, the correlated features are isolated from each other, making the dominant modes more pronounced, and therefore easier to choose when the algorithm reconstructs the residuals.  Second, this data segmentation allows us to implement multiprocessing on the calculations, which greatly reduces the calculation time.

Normalization
+++++++++++++
Normalize the spectra to be inserted into the SVD calculation.

The method of normalization has a strong effect on any SVD calculation. In this approach, we normalize to the variance for each spectrum within the operating spectral segments.


Singular Value Decomposition
++++++++++++++++++++++++++++
Calulate the eigenspectra and eigenvalues that characterize the residual emission line features.

Using the previous steps, we have produced a set of spectra that are prepared for the Singular Value Decomposition, which cacluated the eigenspectra and the eigenvalues that create the entire set of input spectra. To use these eigenspectra effectively, the eigenmodes must be chosen to identify only the contributions by the sky residuals.

Optimization
++++++++++++

