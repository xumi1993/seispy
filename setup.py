#!/usr/bin/env python
from setuptools import find_packages, setup
packages = find_packages()

with open("README.md", "r") as fh:
    long_description = fh.read()


VERSION = "1.3.0"
setup(name='python-seispy',
      version=VERSION,
      author='Mijian Xu',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author_email='gomijianxu@gmail.com',
      license='GPLv3',
      packages=find_packages(),
      package_dir={'seispy': 'seispy'},
      package_data={'': ['data/*']},
      install_requires=[
                'numpy>=1.19.0',
                'scipy>=1.9.1',
                'matplotlib>=3.5.0',
                'pandas>=1.0.0',
                'obspy>=1.2.1',
                'pyside6>=6.3.0',
                'scikits.bootstrap>=1.0.0'],
      entry_points={'console_scripts': ['gen_rayp_lib=seispy.psrayp:gen_rayp_lib',
                                        'prf=seispy.scripts:prf',
                                        'srf=seispy.scripts:srf',
                                        'setpar=seispy.rf:setpar',
                                        'rf2depth=seispy.rf2depth_makedata:rf2depth',
                                        'plotrt=seispy.scripts:plot_rt',
                                        'plotr=seispy.scripts:plot_r',
                                        'download_catalog=seispy.catalog:main',
                                        'ccp_profile=seispy.scripts:ccp_profile',
                                        'hk=seispy.hk:hk',
                                        'pickrf=seispy.pickrf.pickui:main',
                                        'pickdepth=seispy.pickdepth.pickdepthui:main',
                                        'rfani=seispy.scripts:rfani',
                                        'ccp3d=seispy.scripts:ccp3d',
                                        'rfharmo=seispy.scripts:rfharmo',
                                        'get_pierce_points=seispy.scripts:get_pierce_points',
                                        'veltxt2mod=seispy.modcreator:veltxt2mod']},
      zip_safe=False,
      classifiers=['Programming Language :: Python',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10',
                   'Programming Language :: Python :: 3.11']
      )
