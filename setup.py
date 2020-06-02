#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import setuptools

def readme():
    with open('README.md') as f:
        return f.read()

setuptools.setup(name='quickpbsa',
                 version='2020.0.1',
                 description='Fast and Complete Photobleaching Step Analysis',
                 long_description=readme(),
                 classifiers=[
                              'Development Status :: 3 - Alpha',
                              'GPLv3 License',
                              'Programming Language :: Python :: 3',
                              'Topic :: Scientific/Engineering',
                              ],
                url='https://github.com/JohnDieSchere/quickpbsa',
                author='Johan Hummert',
                author_email='johndieschere@gmail.com',
                license='GPLv3',
                packages=setuptools.find_packages(),
                install_requires=[
                                  'numpy',
                                  'scipy',
                                  'pandas',
                                  'sympy'
                                  ],
                extras_require={
                                'Trace Extraction': ['tifffile', 'matplotlib']
                                },
                zip_safe=False)