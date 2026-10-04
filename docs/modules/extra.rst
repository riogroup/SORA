sora.extra
----------

ChiSquare Class
===============
.. autoclass:: sora.extra.ChiSquare
   :members:

``ChiSquare.get_nsigma()`` excludes non-finite chi-square trials when computing
the minimum and confidence intervals. It raises ``ValueError`` if no finite
trial is available, and leaves the stored trial data unchanged.

Plot ellipse
============
.. automodule:: sora.extra.plots
   :members:


Complementary functions
=======================
.. automodule:: sora.extra.utils
   :members:

