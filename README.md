# <img src="https://user-images.githubusercontent.com/7437523/128596331-dc5c5e40-93e1-4d9e-b92d-9c53fe51145a.png" width="500"/> 


[![License](https://img.shields.io/github/license/xumi1993/seispy)]()
[![](https://img.shields.io/github/last-commit/xumi1993/seispy)]()
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/xumi1993/seispy)]()
[![GitHub repo size](https://img.shields.io/github/repo-size/xumi1993/seispy)]()
[![Static Badge](https://img.shields.io/badge/DOI-10.1785%2F0220220288-pink)](https://doi.org/10.1785/0220220288)


[![CRV test](https://github.com/xumi1993/seispy/actions/workflows/test.yml/badge.svg?branch=dev)](https://github.com/xumi1993/seispy/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/xumi1993/seispy/branch/dev/graph/badge.svg?token=XN3E3N6S3V)](https://codecov.io/gh/xumi1993/seispy)
[![Upload Python Package](https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml/badge.svg)](https://github.com/xumi1993/seispy/actions/workflows/python-publish.yml)
[![Seispy docs](https://github.com/xumi1993/seispy/actions/workflows/doc_build.yml/badge.svg)](https://github.com/xumi1993/seispy/actions/workflows/doc_build.yml)
<a href="https://dev.azure.com/conda-forge/feedstock-builds/_build/latest?definitionId=13623&branchName=master">
  <img src="https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/seispy-feedstock?branchName=master">
</a> 

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/seispy/badges/version.svg)](https://anaconda.org/conda-forge/seispy)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/seispy.svg)](https://anaconda.org/conda-forge/seispy)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/seispy/badges/downloads.svg)](https://anaconda.org/conda-forge/seispy)

[![PyPI](https://img.shields.io/pypi/v/python-seispy)](https://pypi.org/project/python-seispy/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/python-seispy)](https://pypi.org/project/python-seispy/)

[![GitHub stars](https://img.shields.io/github/stars/xumi1993/seispy?style=social)]()
[![](https://img.shields.io/github/forks/xumi1993/seispy?style=social)]()


Seispy is a Python module for processing seismological data and calculating Receiver Functions. The advanced functions are available to improve the Obspy.

## Acknowledgements

For the use of the Seispy package, please cite as:

- Xu, M. & He, J. (2023). Seispy: Python Module for Batch Calculation and Postprocessing of Receiver Functions. Seismological Research Letters, 94 (2A): 935–943. [![Static Badge](https://img.shields.io/badge/DOI-10.1785%2F0220220288-pink)](https://doi.org/10.1785/0220220288)

For 3D time-difference correction, please also consider citing:

- Xu, M., Huang, H., Huang, Z., Wang, P., Wang, L., Xu, M., ... & Yuan, X. (2018). Insight into the subducted Indian slab and origin of the Tengchong volcano in SE Tibet from receiver function analysis. Earth and Planetary Science Letters, 482, 567-579.  [![Static Badge](https://img.shields.io/badge/DOI-10.1016%2Fj.epsl.2017.11.048-pink)
](https://doi.org/10.1785/0220220288)

- Xu, M., Huang, Z., Wang, L., Xu, M., Mi, N., & Yu, D. (2020). Lateral variation of the mantle transition zone beneath the Tibetan plateau: Insight into thermal processes during Indian–Asian collision. Physics of the Earth and Planetary Interiors, 301, 106452. [![Static Badge](https://img.shields.io/badge/DOI-10.1016%2Fj.pepi.2020.106452-pink)](https://doi.org/10.1016/j.pepi.2020.106452)

For 2D and 3D CCP stacking, please also consider citing:

- Xu, M., Huang, Z., Wang, L., Xu, M., Zhang, Y., Mi, N., ... & Yuan, X. (2020). Sharp lateral Moho variations across the SE Tibetan margin and their implications for plateau growth. Journal of Geophysical Research: Solid Earth, 125(5), e2019JB018117. [![Static Badge](https://img.shields.io/badge/DOI-10.1029%2F2019JB018117-pink)](https://doi.org/10.1029/2019JB018117)

## Installation

See [Seispy documentation](https://seispy.xumijian.me/installation.html) in detail.
 
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
 * `rfharmo`: Harmonic decomposition to extract constant component of RF and plot dip/anisotropic components.
 * `pickdepth`: Pick depth of stacked PRFs following `ccp3d`

### Others
 * `veltxt2mod`: Create 3D velocity model with `numpy.lib.npyio.NpzFile` format from a ASCII table file.
 * `downlod_catalog`: Download catalogs from FDSN web-service.
 * `gen_rayp_lib`: Generate a rap-parameter library with depth of source and epicentral distance.
 * `setpar`: Set up the values in configure files.

