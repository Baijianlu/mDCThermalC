# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 10:09:30 2019

@author: Tao.Fan
This script used to setup the scripts
"""
import os
import sys
import numpy

try:
    from setuptools import setup, Extension
    use_setuptools = True
    print("setuptools is used.")
except ImportError:
    from distutils.core import setup, Extension
    use_setuptools = False
    print("distutils is used.")
    
extension_mDCThermalC = Extension('mdcthermalc._extern', sources = [])
ext_modules_mDCThermalC = [extension_mDCThermalC]
packages_mDCThermalC = ['mdcthermalc']
scripts_mDCThermalC = ['Scripts/mDCThermalC']

if __name__ == '__main__':

    version_nums = [None, None, None]
    with open("mdcthermalc/version.py") as f:
        for line in f:
            if "__version__" in line:
                for i, num in enumerate(line.split()[2].strip('\"').split('.')):
                    version_nums[i] = int(num)
                break


    if None in version_nums:
        print("Failed to get version number in setup.py.")
        raise

    version_number = ".".join(["%d" % n for n in version_nums])

    if use_setuptools:
        setup(name='mdcthermalc',
              version=version_number,
              description='This is the mDCThermalC module.',
              author='Tao Fan',
              author_email='Tao.Fan@skoltech.ru',
              url='http://atztogo.github.io/phonopy/',
              packages=packages_mDCThermalC,
              install_requires=['numpy', 'scipy', 'pymatgen'],
              provides=['mdcthermalc'],
              scripts=scripts_mDCThermalC,
              ext_modules=ext_modules_mDCThermalC)
    else:
        setup(name='mdcthermalc',
              version=version_number,
              description='This is the mDCThermalC module.',
              author='Tao Fan',
              author_email='Tao.Fan@skoltech.ru',
              url='http://atztogo.github.io/phonopy/',
              packages=packages_mDCThermalC,
              requires=['numpy', 'scipy', 'pymatgen'],
              provides=['mdcthermalc'],
              scripts=scripts_mDCThermalC,
              ext_modules=ext_modules_mDCThermalC)
