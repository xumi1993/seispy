#!/usr/bin/env python
from setuptools import find_packages, setup
packages = find_packages()

VERSION = "1.1.20"
setup(name='seispy',
      version=VERSION,
      author='Mijian Xu',
      author_email='gomijianxu@gmail.com',
      license='GPLv3',
      packages=find_packages(),
      package_dir={'seispy': 'seispy'},
      package_data={'': ['data/*']},
      install_requires=[
                'netcdf4>=1.5.2',
                'pyerf>=1.0.1',
                'obspy>=1.2.0',
                'pandas>=1.0.0',
                'numpy>=1.19.0',
                'scipy>=1.1.0',
                'matplotlib>=3.2.0'],
      entry_points={'console_scripts': ['gen_rayp_lib=seispy.psrayp:gen_rayp_lib',
                                        'prf=seispy.rf:prf',
                                        'setpar=seispy.rf:setpar',
                                        'rf2depth=seispy.rf2depth_makedata:rf2depth',
                                        'plotrt=seispy.plotRT:main',
                                        'plotr=seispy.plotR:main',
                                        'updatecatalog=seispy.updatecatalog:main',
                                        'ndk2dat=seispy.updatecatalog:ndk2dat',
                                        'lsnc=seispy.io:lsnc',
                                        'ccp_profile=seispy.ccp:ccp_profile',
                                        'hk=seispy.hk:hk',
                                        'pickrf=seispy.pickui:main']},
      include_package_data=True,
      zip_safe=False
      )
