<!-- .. seispy documentation master file, created by
   sphinx-quickstart on Sat Jul 18 19:55:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive. -->

# Seispy Documentation

-----------------------

[![CRV test](https://github.com/xumi1993/seispy/actions/workflows/test.yml/badge.svg?branch=dev)](https://github.com/xumi1993/seispy/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/xumi1993/seispy/branch/dev/graph/badge.svg?token=XN3E3N6S3V)](https://codecov.io/gh/xumi1993/seispy)

[![Upload Python Package](https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml/badge.svg)](https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml)
[![Deploy Seispy Docs](https://github.com/xumi1993/seispy-doc.post/actions/workflows/deploy.yml/badge.svg)](https://github.com/xumi1993/seispy-doc.post/actions/workflows/deploy.yml)
<a href="https://dev.azure.com/conda-forge/feedstock-builds/_build/latest?definitionId=13623&branchName=master">
  <img src="https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/seispy-feedstock?branchName=master">
</a> 

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/seispy.svg)](https://anaconda.org/conda-forge/seispy)
[![PyPI](https://img.shields.io/pypi/v/python-seispy)](https://pypi.org/project/python-seispy/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/python-seispy)](https://pypi.org/project/python-seispy/)
[![GitHub](https://img.shields.io/github/license/xumi1993/seispy)](https://github.com/xumi1993/seispy/blob/master/LICENSE.txt)
[![](https://img.shields.io/github/last-commit/xumi1993/seispy)](https://github.com/xumi1993/seispy)
[![](https://img.shields.io/github/commit-activity/m/xumi1993/seispy)](https://github.com/xumi1993/seispy)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/xumi1993/seispy)](https://github.com/xumi1993/seispy)
[![GitHub repo size](https://img.shields.io/github/repo-size/xumi1993/seispy)](https://github.com/xumi1993/seispy)
[![DOI](https://zenodo.org/badge/41006349.svg)](https://zenodo.org/badge/latestdoi/41006349)

[![GitHub stars](https://img.shields.io/github/stars/xumi1993/seispy?style=social)](https://github.com/xumi1993/seispy)
[![](https://img.shields.io/github/forks/xumi1993/seispy?style=social)](https://github.com/xumi1993/seispy)

Seispy is a Python module for processing seismological data and calculating Receiver Functions. The advanced functions are available to improve the Obspy.

I have been writing Seispy when I was a master student. At first, I wanted to calculate Receiver Functions in Python, but there is no suitable toolkit. Fortunately, The Obspy provided mounts of APIs for processing seismic data, so I ported codes for calculating Receiver Functions from Matlab to Python. Today increased functions have been added to Seispy to further process seismic data over than Obspy.

<img src="_static/seispy_shortcut.png"/> 

## Installation

::::{panels}
:container: full-width
----
:header: bg-myst-one
**Install via PyPI**
^^^^
```
pip install python-seispy
```
---
:header: bg-myst-two
**Install via conda**
^^^^
```
conda install seispy -c conda-forge
```
::::
Further details are available in the [Installation Guide](installation).


## Learning resources

:::{seealso}
Access [here](https://www.koushare.com/video/videodetail/14734) for video tutorial in Chinese
:::

::::{panels} 
:container: full-width
:column: col-lg-4 p-2
---
:header: bg-myst-two
**Getting Started**
^^^
- [Introduction to receiver function](usage/PRF_Process)
- [Procedure of RF calculation](usage/rf-sta)
- [Visual check RFs](usage/pickrf)
- [Time-to-depth & CCP stacking](usage/ccp)
---
:header: bg-myst-three
**Examples**
^^^
- [RF calculation & H-k stacking](examples/ex-prf)
- [CCP staking along profile](examples/ex-ccp)
- [3D CCP stacking](examples/ex-ccp3d)
---
:header: bg-myst-three
**Configuration**
^^^
- [RF Configure](notes/rf_cfg)
- [Hk Configure](notes/hk_cfg)
- [CCP Configure](notes/ccp_cfg)

::::

## Libraries

- `seispy.distaz`: Calculate distance and azimuth credited by the [lithospheric seismology program at USC](http://www.seis.sc.edu/software/distaz/), but `numpy.ndarray` operations are supported.
- `seispy.geo`: Tiny codes of geophysics.
- `seispy.decon`: Functions of deconvolution transferred from [iwbailey/processRFmatlab](https://github.com/iwbailey/processRFmatlab) including
  - Iterative time domain deconvolution method (Ligorría and Ammon 1999 BSSA). 
  - Water level frequency domain deconvolution method (CJ. Ammon 1991 BSSA)
- `seispy.rf`: Procedure for RF calculation. The functions of `match_eq`, `search_eq` invoked `obspy.core.UTCDateTime` and `obspy.clients` from the [Obspy](https://docs.obspy.org/).
- `seispy.eq`: RF processing for each event, which invoked `obspy.io.sac`, `obspy.signal`, `obspy.taup` and `obspy.core.Stream` from the [Obspy](https://docs.obspy.org/).
- `seispy.hk`: H-k stacking for single station (Zhu and Kanamori 2000 JGR).
- `seispy.rfani`: A joint method for crustal anisotropic calculation (Liu and Niu 2011 GJI).
- `seispy.slantstack`: Slant stacking for single station (Tauzin et al., 2008)
- `seispy.rfcorrect`: Subsequent process of RFs including moveout correction and time to depth conversion (1D and 3D) (see [Xu et al., 2018 EPSL](https://www.sciencedirect.com/science/article/pii/S0012821X17306921?via%3Dihub))
- `seispy.ccpprofile`: CCP stacking along a profile.
- `seispy.ccp3d`: 3-D CCP stacking with extracting depth D410 and D660.


:::{toctree}
:maxdepth: 1
:hidden:

installation
usage/index
examples/index
notes/index
modules
:::


## References


Ammon C J. The isolation of receiver effects from teleseismic P waveforms[J]. Bulletin-Seismological Society of America, 1991, 81(6): 2504-2510.

Liu H, Niu F. Estimating crustal seismic anisotropy with a joint analysis of radial and transverse receiver function data[J]. Geophysical Journal International, 2012, 188(1): 144-164.

Krischer L, Megies T, Barsch R, et al. ObsPy: A bridge for seismology into the scientific Python ecosystem[J]. Computational Science & Discovery, 2015, 8(1): 014003.

Ligorría J P, Ammon C J. Iterative deconvolution and receiver-function estimation[J]. Bulletin of the seismological Society of America, 1999, 89(5): 1395-1400.

Tauzin B, Debayle E, Wittlinger G. The mantle transition zone as seen by global Pds phases: No clear evidence for a thin transition zone beneath hotspots[J]. Journal of Geophysical Research: Solid Earth, 2008, 113(B8).

Xu M, Huang H, Huang Z, et al. Insight into the subducted Indian slab and origin of the Tengchong volcano in SE Tibet from receiver function analysis[J]. Earth and Planetary Science Letters, 2018, 482: 567-579.

Zhu, L., and Kanamori, H. Moho depth variation in southern California from teleseismic receiver functions[J]. J. Geophys. Res., 2000, 105( B2), 2969– 2980, doi:10.1029/1999JB900322.
