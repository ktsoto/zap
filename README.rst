ZAP (the Zurich Atmosphere Purge)
---------------------------------

Tired of sky subtraction residuals? ZAP them!

ZAP is a high precision sky subtraction tool which can be used as complete sky
subtraction solution, or as an enhancement to previously sky-subtracted MUSE data.
The method uses PCA to isolate the residual sky subtraction features and remove
them from the observed datacube. Future developments will include modification for
use on a variety of instruments.

The last stable release of ZAP can be installed simply with::

    pip install zap

Or into the user path with::

    pip install --user zap

Documentation : `zap.readthedocs.org <http://zap.readthedocs.org/en/latest/>`_

Read about the method here: `arXiv <http://arxiv.org/abs/1602.08037>`_

ZAP can be cited via the NASA ADS http://adsabs.harvard.edu/abs/2016arXiv160208037S

For questions on the use of zap please contact Kurt soto at kurt.soto@phys.ethz.ch
