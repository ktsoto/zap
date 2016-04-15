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

Please cite ZAP as::

\bibitem[Soto et al.(2016)]{2016MNRAS.458.3210S} Soto, K.~T., Lilly, S.~J., Bacon, R., Richard, J., \& Conseil, S.\ 2016, \mnras, 458, 3210 

The paper describing the method can be found here: http://adsabs.harvard.edu/abs/2016MNRAS.458.3210S

For questions on the use of zap please contact Kurt soto at kurt.soto@phys.ethz.ch
