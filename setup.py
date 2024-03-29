######!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# Import of the lib pyswath
import pyswath

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='pyswath',
	#version='0.1.1',
	version=pyswath.__version__,
	description='package that provides tools to extract swath profiles from a raster map',
	long_descritpion=open('README.rst').read(),
	url='https://github.com/robertxa/pyswath',
	download_url='https://github.com/robertxa/pyswath/archive/master.zip',
	author='Xavier Robert',
	author_email='xavier.robert@univ-grenoble-alpes.fr',
	license='GPL-V3.0',
	packages=find_packages(),
	#include_package_data=True,	# What is the use of it ?
	install_requires=[
	      'rasterstats>=0.11.0',
	      'shapely',
	      'numpy',
	      'utm',
	      'gdal',
	      'progress'
	],
	classifiers=[
		"Operating System :: OS Independent",
		"Topic :: Scientific/Engineering :: GIS"
	],
	include_package_data=True,
	zip_safe=False)
      