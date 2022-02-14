.. seispy documentation master file, created by
   sphinx-quickstart on Sat Jul 18 19:55:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Seispy Documentation
==========================

.. image:: https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml/badge.svg
   :target: https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml
.. image:: https://github.com/xumi1993/seispy-doc.post/actions/workflows/deploy.yml/badge.svg
   :target: https://github.com/xumi1993/seispy-doc.post/actions/workflows/deploy.yml
.. image:: https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/seispy-feedstock?branchName=master
   :target: https://dev.azure.com/conda-forge/feedstock-builds/_build/latest?definitionId=13623&branchName=master
.. image:: https://img.shields.io/conda/vn/conda-forge/seispy.svg
   :target: https://anaconda.org/conda-forge/seispy

.. image:: https://img.shields.io/pypi/v/python-seispy
   :target: https://pypi.org/project/python-seispy
.. image:: https://img.shields.io/pypi/pyversions/python-seispy
.. image:: https://img.shields.io/github/license/xumi1993/seispy
.. image:: https://img.shields.io/github/v/tag/xumi1993/seispy
.. image:: https://img.shields.io/github/last-commit/xumi1993/seispy
.. image:: https://img.shields.io/github/commit-activity/m/xumi1993/seispy
.. image:: https://img.shields.io/github/languages/code-size/xumi1993/seispy
.. image:: https://img.shields.io/github/repo-size/xumi1993/seispy
.. image:: https://zenodo.org/badge/41006349.svg
   :target: https://zenodo.org/badge/latestdoi/41006349


.. image:: https://img.shields.io/github/stars/xumi1993/seispy?style=social
   :target: https://github.com/xumi1993/seispy
.. image:: https://img.shields.io/github/forks/xumi1993/seispy?style=social
   :target: https://github.com/xumi1993/seispy


Seispy is a Python module for processing seismological data and calculating Receiver Functions. The advanced functions are available to improve the Obspy.

I have been writing Seispy when I was a master student. At first, I wanted to calculate Receiver Functions in Python, but there is no suitable toolkit. Fortunately, The Obspy provided mounts of APIs for processing seismic data, so I ported codes for calculating Receiver Functions from Matlab to Python. Today increased functions have been added to Seispy to further process seismic data over than Obspy.

.. image:: _static/seispy_shortcut.png
   :height: 719
   :width: 1182
   :scale: 50
   :align: center


Libraries
------------

- ``seispy.distaz``: Calculate distance and azimuth (by the `lithospheric seismology program at USC <http://www.seis.sc.edu/software/distaz/>`_).
- ``seispy.geo``: Tiny codes of geophysics.
- ``seispy.decon``: Functions of deconvolution transferred from `iwbailey/processRFmatlab <https://github.com/iwbailey/processRFmatlab>`_ including

  - Iterative time domain deconvolution method (Ligorría and Ammon 1999 BSSA). 
  
  - Water level frequency domain deconvolution method (CJ. Ammon 1991 BSSA)
- ``seispy.rf``: Procedure for RF calculation. The functions of ``match_eq``, ``search_eq`` invoked ``obspy.core.UTCDateTime`` and ``obspy.clients`` from the `Obspy <https://docs.obspy.org/>`_.
- ``seispy.eq``: RF processing for each event, which invoked ``obspy.io.sac``, ``obspy.signal``, ``obspy.taup`` and ``obspy.core.Stream`` from the `Obspy <https://docs.obspy.org/>`_.
- ``seispy.hk``: H-k stacking for single station (Zhu and Kanamori 2000 JGR).
- ``seispy.rfani``: A joint method for crustal anisotropic calculation (Liu and Niu 2011 GJI).
- ``seispy.slantstack``: Slant stacking for single station (Tauzin et al., 2008)
- ``seispy.rfcorrect``: Subsequent process of RFs including moveout correction and time to depth conversion (1D and 3D) (see `Xu et al., 2018 EPSL <https://www.sciencedirect.com/science/article/pii/S0012821X17306921?via%3Dihub>`_)
- ``seispy.ccpprofile``: CCP stacking along a profile.
- ``seispy.ccp3d``: 3-D CCP stacking with extracting depth D410 and D660.


.. toctree::
   :maxdepth: 1

   usage/index
   examples/index
   notes/index
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


References
============

Ammon C J. The isolation of receiver effects from teleseismic P waveforms[J]. Bulletin-Seismological Society of America, 1991, 81(6): 2504-2510.

Krischer L, Megies T, Barsch R, et al. ObsPy: A bridge for seismology into the scientific Python ecosystem[J]. Computational Science & Discovery, 2015, 8(1): 014003.

Ligorría J P, Ammon C J. Iterative deconvolution and receiver-function estimation[J]. Bulletin of the seismological Society of America, 1999, 89(5): 1395-1400.

Xu M, Huang H, Huang Z, et al. Insight into the subducted Indian slab and origin of the Tengchong volcano in SE Tibet from receiver function analysis[J]. Earth and Planetary Science Letters, 2018, 482: 567-579.