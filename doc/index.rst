.. zap documentation master file, created by
   sphinx-quickstart on Mon Nov 25 09:46:49 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ZAP's documentation!
===============================

.. toctree::
   :maxdepth: 2

ZAP (the Zurich Atmosphere Purge) is a high precision sky subtraction tool which can be used as complete sky subtraction solution, or as an enhancement to previously sky-subtracted MUSE integral field spectroscopic data.  The method uses PCA to isolate the residual sky subtraction features and remove them from the observed datacube. Though the operation of ZAP is not dependent on perfect flatfielding of the data in a MUSE exposure, better results are obtained when these corrections are made ahead of time. Future development will include expansion to more instruments.

Installation
============

Requirements
------------
ZAP requires the following packages:

* Numpy 1.6.0 or later
* Astropy v1.0 or later
* SciPy v0.13.3 or later

Many linear algebra operations are performed in ZAP, so it can be beneficial to use an alternative BLAS package. In the Anaconda distribution, the default BLAS comes with Numpy linked to OpenBlas, which can amount to a 20% speedup of ZAP.

Steps
-----

ZAP can be installed using pip ::
  pip install zap


Examples
========

In its most hands-off form, ZAP can take an input fits datacube, operate on it, and output a final fits datacube::

  import zap
  zap.process('INPUT.fits', outcubefits='OUTPUT.fits')

Care should be taken, however, since this case assumes a sparse field, and better results can be obtained by applying masks.

The main function is ``zap.process``:

There are a number of options that can be passed to the code which we describe here:

.. autofunction:: zap.process

The code can handle datacubes trimmed in wavelength space. Since the code uses the correlation of segments of the emission line spectrum, it is best to trim the cube at specific wavelengths. The cube can include any connected subset of these segments. (for example 6400 - 8200 Angstroms) ::

  [0,    5400]
  [5400, 5850]
  [5850, 6440]
  [6440, 6750]
  [6750, 7200]
  [7200, 7700]
  [7700, 8265]
  [8265, 8602]
  [8602, 8731]
  [8731, 9275]
  [9275, 10000]


Sparse Field Case
-----------------

This case specifically refers to the case where the sky can be measured in the sky frame itself, using::

  zap.process('INPUT.fits', outcubefits='OUTPUT.fits')

In both cases, the code will create a resulting processed datacube named ``DATACUBE_ZAP.fits`` and an SVD file named ``ZAP_SVD.fits`` in the current directory. While this can work well in the case of very faint sources, masks can improve the results.

For the sparse field case, a mask file can be included, which is a 2d fits image matching the spatial dimensions of the input datacube. Masks are defined to be >= 1 on astronomical sources and 0 at the position of the sky. Set this parameter with the ``mask`` keyword ::

  zap.process('INPUT.fits', outcubefits='OUTPUT.fits', mask='mask.fits')

Filled Field Case
-----------------

This approach also can address the saturated field case and is robust in the case of strong emission lines, in this case the input is an offset sky observation. To achieve this, we calculate the SVD on an external sky frame using the function ``zap.SVDoutput``

.. autofunction:: zap.SVDoutput

An example of running the code in this way is as follows::

  zap.SVDoutput('Offset_Field_CUBE.fits', svdfn='ZAP_SVD.fits', mask='mask.fits')
  zap.process('Source_cube.fits', outcubefits='OUTPUT.fits', extSVD='ZAP_SVD.fits', cfwidthSP=50)

The integration time of this frame does not need to be the same as the object exposure, but rather just a 2-3 minute exposure. Often residuals can be further reduced by changing `cfwidthSP` to a smaller value. However, this parameter should not be reduced to smaller than 15 pixels.

Extra Functions
===============

Aside from the main process, two functions are included that can be run outside of the entire zap process to facilitate some investigations.

.. autofunction:: zap.nancleanfits

.. autofunction:: zap.wmedian

		  
Command Line Interface
======================

ZAP can also be used from the command line::
  
  python -m zap INPUT_CUBE.fits

More information use of the command line interface can be found with the command ::
  
  python -m zap -h


Interactive mode
================

ZAP can also  be used interactively from within IPython ::

  import zap
  zobj = zap.process('INPUT.fits', interactive=True)

The run method operates on the datacube, and retains all of the data and
methods necessary to process a final data cube in a python class named
:class:`~zap.zclass`. You can elect to investigate the data product via the
:class:`~zap.zclass`, and even reprocess the cube with a different number of
eigenspectra per region.  A workflow may go as follows:

.. code-block:: python

  import zap
  from matplotlib import pyplot as plt

  # allow ZAP to run the optimize routine
  zobj = zap.process('INPUT.fits', optimization='normal', interactive=True)

  # plot the variance curves and the selection of the number of eigenspectra used
  zobj.plotvarcurve(5)

  # plot a spectrum extracted from the original cube
  plt.figure()
  plt.plot(zobj.cube[:,50:100,50:100].sum(axis=(1,2)), 'b', alpha=0.3)

  # plot a spectrum of the cleaned ZAP dataproduct
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2)), 'g')

  # choose just the first 3 spectra for all segmments
  zobj.reprocess(nevals=3)

  # plot a spectrum extracted from the original cube
  plt.plot(zobj.cube[:,50:100,50:100].sum(axis=(1,2)), 'b', alpha=0.3)

  # plot a spectrum of the cleaned ZAP dataproduct
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2))), 'g')

  # choose some number of modes by hand
  zobj.reprocess(nevals=[2,5,2,4,6,7,9,8,5,3,5])

  # plot a spectrum
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2))), 'k')

  # Use the optimization algorithm to identify the best number of modes per segment
  zobj.optimize()

  # compare to the previous versions
  plt.plot(zobj.cleancube[:,50:100,50:100].sum(axis=(1,2))), 'r')

  # identify a pixel in the dispersion axis that shows a residual feature in the original
  plt.figure()
  plt.matshow(zobj.cube[2903,:,:])

  # compare this to the zap dataproduct
  plt.figure()
  plt.matshow(zobj.cleancube[2903,:,:])

  # write the processed cube as a single extension fits
  zobj.writecube('DATACUBE_ZAP.fits')

  # or merge the zap datacube into the original input datacube, replacing the data extension
  zobj.writefits(outcubefits='DATACUBE_FINAL_ZAP.fits')

======
ZCLASS
======
.. autoclass:: zap.zclass
   :members:


