#!/usr/bin/env python
from setuptools import find_packages, setup
packages=find_packages()

VERSION = "1.0"
# url='https://github.com/trichter/rf', 
#      scripts = ['PlotR.py', 'PyCCP.py', 'PyCCP_Boot.py', 'makeRF4station.py'], 
setup(name='seispy',
      version=VERSION,
      author='Mijian Xu',
      author_email='gomijianxu@gmail.com',
      license='GNU',
      packages=find_packages(),
      package_dir={'seispy': 'seispy'},
      install_requires=['obspy'],
      include_package_data=True,
      zip_safe=False
      )
