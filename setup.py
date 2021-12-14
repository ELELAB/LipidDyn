#!/usr/bin/env python

#    setup.py: installation script for the LypidDyn pipeline
#    Copyright (C) 2020, SimoneScrima <simonescrima@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

setup(
    name="LipidDyn",
    version="0.3",
    packages=find_packages(include=['LipidDyn',
                                    'LipidDyn.*']),
    scripts=['bin/LipidDyn',
             'bin/profiler',
             'bin/ordpar',
             'bin/diffusion',
             'bin/dmaps'],
    install_requires=[ 'matplotlib',
                       'numpy',
                       'MDAnalysis',
                       'progressbar',
                       'seaborn',
                       'pyyaml',
                       'lipyphilic'],
    include_package_data=True,
    zip_safe = False,
    author="Simone Scrima, Matteo Tiberti, Matteo Lambrughi,Elena Papaleo",
    author_email="simonescrima@gmail.com",
    description="Scripts to automate pipeline of analysis for Lipid Simulations",
    keywords="Lipid pipeline analysis",
    url="https://github.com/ELELAB/LipidDyn",
)

