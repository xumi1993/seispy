# <img src="https://user-images.githubusercontent.com/7437523/128596331-dc5c5e40-93e1-4d9e-b92d-9c53fe51145a.png" width="500"/> 

[![Upload Python Package](https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml/badge.svg)](https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml)
[![Deploy Seispy Docs](https://github.com/xumi1993/seispy-doc.post/actions/workflows/deploy.yml/badge.svg)](https://github.com/xumi1993/seispy-doc.post/actions/workflows/deploy.yml)
<a href="https://dev.azure.com/conda-forge/feedstock-builds/_build/latest?definitionId=13623&branchName=master">
  <img src="https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/seispy-feedstock?branchName=master">
</a> 

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/seispy.svg)](https://anaconda.org/conda-forge/seispy)
[![PyPI](https://img.shields.io/pypi/v/python-seispy)](https://pypi.org/project/python-seispy/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/python-seispy)](https://pypi.org/project/python-seispy/)
[![GitHub](https://img.shields.io/github/license/xumi1993/seispy)]()
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/xumi1993/seispy)]()
[![](https://img.shields.io/github/last-commit/xumi1993/seispy)]()
[![](https://img.shields.io/github/commit-activity/m/xumi1993/seispy)]()
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/xumi1993/seispy)]()
[![GitHub repo size](https://img.shields.io/github/repo-size/xumi1993/seispy)]()
[![DOI](https://zenodo.org/badge/41006349.svg)](https://zenodo.org/badge/latestdoi/41006349)

[![GitHub stars](https://img.shields.io/github/stars/xumi1993/seispy?style=social)]()
[![](https://img.shields.io/github/forks/xumi1993/seispy?style=social)]()


Seispy is a Python module for processing seismological data and calculating Receiver Functions. The advanced functions are available to improve the Obspy.


# Installation
## Dependencies
  * [Python]() >= 3.6
  * [ObsPy](http://docs.obspy.org) >= 1.1.0
  * [NumPy](http://www.numpy.org/) >= 1.16
  * [SciPy](http://www.scipy.org/) >= 1.2.0
  * [matplotlib](https://matplotlib.org/) >= 3.0.0
  * [PyQt5](https://www.riverbankcomputing.com/software/pyqt/) >= 5.12.0
  * [scikits.bootstrap](https://github.com/cgevans/scikits-bootstrap) >= 1.0.0
  
## Installation
```
git clone https://github.com/xumi1993/seispy.git
python setup.py install
```

# Inclusion
## Libraries

-   `seispy.distaz`: Calculate distance and azimuth (by the
    [lithospheric seismology program at USC][]).
-   `seispy.geo`: Tiny codes of geophysics.
-   `seispy.decon`: Functions of deconvolution transferred from
    [iwbailey/processRFmatlab][] including
    -   Iterative time domain deconvolution method (Ligorría and Ammon
        1999 BSSA).
    -   Water level frequency domain deconvolution method (CJ. Ammon
        1991 BSSA)
-   `seispy.rf`: Procedure for RF calculation. The functions of
    `match_eq`, `search_eq` invoked `obspy.core.UTCDateTime` and
    `obspy.clients` from the [Obspy][].
-   `seispy.eq`: RF processing for each event, which invoked
    `obspy.io.sac`, `obspy.signal`, `obspy.taup` and `obspy.core.Stream`
    from the [Obspy][].
-   `seispy.rfcorrect`: Subsequent process of RFs including moveout
    correction and time to depth conversion (1D and 3D) (see [Xu et al., 2018 EPSL](https://www.sciencedirect.com/science/article/pii/S0012821X17306921?via%3Dihub))
-   `seispy.ccpprofile`: CCP stacking along a profile.
-   `seispy.ccp3d`: 3-D CCP stacking with extracting depth D410 and
    D660.
-   `seispy.get_cpt`: Convert color map from the `cpt` format to the `matplotlib.cmap` modified from [bouziot/get-cpt](https://github.com/bouziot/get-cpt) based on the GPLv3 license.

  [lithospheric seismology program at USC]: http://www.seis.sc.edu/software/distaz/
  [scikits-bootstrap]: https://github.com/cgevans/scikits-bootstrap
  [iwbailey/processRFmatlab]: https://github.com/iwbailey/processRFmatlab
  [Obspy]: https://docs.obspy.org/
  [Xu et al., 2018 EPSL]: https://www.sciencedirect.com/science/article/pii/S0012821X17306921?via%3Dihub


## Commands
### Receiver Functions
 * `prf`: Calculate PRFs for a station.
 * `pickrf`: Pick PRFs with virtual quality control after the calculation.
 * `plotrt`: Plot PRFs with R and T components order by back-azimuth.
 * `plotr`: Plot PRFs with R component order by back-azimuth.
 * `hk`: H-Kappa stacking for estimating Moho depth and crustal Vp/Vs.
 * `rf2depth`: Convert PRFs to depth axis.
 * `ccp_profile`: Stack PRFs along a profile with a CCP stacking method.
 * `ccp3d`: Stack PRFs with spaced bins.
 * `rfani`: Estimating crustal anisotropy with a joint method.

### Others
 * `ndk2dat`: Convert the GCMT catalog file ("ndk" format) to the list file for the `prf` command.
 * `updatecatalog`: Automatically update the GCMT catalog.
 * `setpar`: Set up the values in configure files.

