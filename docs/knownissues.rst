.. _Sec:KnownIssues:


Known Issues
============


Progress Bar Output
^^^^^^^^^^^^^^^^^^^

- Jupyter Notebooks running on Windows OS will experience missing progress bar outputs when setting `"threads"` with a value >1 in the `LightCurve.occ_lcfit()` and  the `Occultation.fit_ellipse()` functions. The issue arises from how the multiprocessing python module was implemented in Windows OS and how the SORA multithreading progress bars piping output was designed. Nonetheless, the progress bar output will still be printed in the execution terminal. [:issue:`79`]

Download Kernels from JPL
^^^^^^^^^^^^^^^^^^^^^^^^^

- The `getBSPfromJPL` function is no longer working due to changes in the JPL query service. Alternatives are being searched.
