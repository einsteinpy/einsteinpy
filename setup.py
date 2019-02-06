#!/usr/bin/env python
import os

from setuptools import find_packages, setup

# https://packaging.python.org/guides/single-sourcing-package-version/
version = {}
with open(os.path.join("src", "einsteinpy", "__init__.py")) as fp:
    exec(fp.read(), version)


# http://blog.ionelmc.ro/2014/05/25/python-packaging/
setup(
    name="einsteinpy",
    version=version["__version__"],
    description="Python package for Orbital Mechanics",
    author="Shreyas Bapat",
    author_email="bapat.shreyas@gmail.com",
    url="https://einsteinpy.readthedocs.io",
    download_url="https://github.com/einsteinpy/einsteinpy",
    license="MIT",
    keywords=["general relativity", "geodesics", "relativity", "tensors"],
    python_requires=">=3.6",
    install_requires=["numpy", "astropy", "matplotlib", "scipy>=1.0"],
    extras_require={
        "dev": [
            "black ; python_version>='3.6'",
            "coverage",
            "isort",
            "pytest",
            "pytest-cov<2.6.0",
            "pycodestyle",
            "sphinx",
            # "sphinx_rtd_theme",  # Use https://github.com/Juanlu001/sphinx_rtd_theme/archive/js-head.zip
            "nbsphinx",
            "ipython>=5.0",
            "jupyter-client",
            "ipykernel",
            "ipywidgets",
        ]
    },
    packages=find_packages("src"),
    package_dir={"": "src"},
    entry_points={"console_scripts": ["einsteinpy = einsteinpy.cli:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    long_description=open("README.rst", encoding="utf-8").read(),
    include_package_data=True,
    zip_safe=False,
)
