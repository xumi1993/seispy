#!/usr/bin/env python
from setuptools import find_packages, setup
packages = find_packages()

VERSION = "1.1.6"
setup(name='seispy',
      version=VERSION,
      author='Mijian Xu',
      author_email='gomijianxu@gmail.com',
      license='GPLv3',
      packages=find_packages(),
      package_dir={'seispy': 'seispy'},
      package_data={'': ['data/*']},
      install_requires=['pyerf', 'obspy', 'pandas', 'numpy', 'scipy', 'matplotlib'],
      entry_points={'console_scripts': ['gen_rayp_lib=seispy.psrayp:gen_rayp_lib',
                                        'rf2depth=seispy.rf2depth_makedata:rf2depth',
                                        'plotrt=seispy.plotRT:main',
                                        'plotr=seispy.plotR:main',
                                        'updatecatalog=seispy.updatecatalog:main',
                                        'ndk2dat=seispy.updatecatalog:ndk2dat',
                                        'lsnc=seispy.io:lsnc',
                                        'ccp_profile=seispy.ccp:ccp_profile']},
      include_package_data=True,
      zip_safe=False
      )
