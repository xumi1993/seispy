#!/usr/bin/env python
from setuptools import find_packages, setup
packages = find_packages()

VERSION = "1.1.5"
setup(name='seispy',
      version=VERSION,
      author='Mijian Xu',
      author_email='gomijianxu@gmail.com',
      license='GPLv3',
      packages=find_packages(),
      package_dir={'seispy': 'seispy'},
      package_data={'': ['data/*.vel']},
      install_requires=['obspy', 'pandas', 'numpy', 'scipy', 'matplotlib'],
      entry_points={'console_scripts': ['gen_rayp_lib=seispy.psrayp:gen_rayp_lib',
                                        'rf2depth=seispy.rf2depth_makedata:rf2depth',
                                        'plotrt=seispy.plotRT:main']},
      include_package_data=True,
      zip_safe=False
      )
