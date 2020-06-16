# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:02:21 2020

@author: rjovelin
"""


from setuptools import setup


__version__ = "1.0.0"

# get the long description from the readme
with open("README.md") as infile:
    content = infile.read().rstrip()

setup(
	name = "barcodex",
	version = __version__,
	author = "Richard Jovelin",
	author_email = "richard.jovelin@oicr.on.ca",
	description = ("A tool for extracting Unique Molecular Identifiers (UMIs) \
                from single or paired-end read sequences"),
	license = "MIT License",
	keywords = "computational genomics",
	url = "https://github.com/oicr-gsi/barcodex",
	py_modules = ['barcodex'],
	long_description = content,
	classifiers = [
	"Development Status :: 3 - Alpha",
	"Intended Audience :: Science/Research",
	"Intended Audience :: Developers",
	"License :: OSI Approved :: MIT License",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 3.6",
	"Topic :: Software Development",
	"Topic :: Scientific/Engineering",
	"Operating System :: POSIX",
	"Operating System :: Unix",
	"Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows"
    ],
	install_requires = ["regex>=2.4.104"],
    python_requires=">=3.6",
)
    
    