#!/usr/bin/env python
import os

from setuptools import setup

# https://packaging.python.org/guides/single-sourcing-package-version/
version = {}
with open(os.path.join("src", "einsteinpy", "__init__.py")) as fp:
    exec(fp.read(), version)


# http://blog.ionelmc.ro/2014/05/25/python-packaging/
setup(    
    version=version["__version__"]
)
