#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#
# fluctmatch --- https://github.com/tclick/python-fluctmatch
# Copyright (c) 2013-2017 The fluctmatch Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the New BSD license.
#
# Please cite your use of fluctmatch in published work:
#
# Timothy H. Click, Nixon Raj, and Jhih-Wei Chu.
# Calculation of Enzyme Fluctuograms from All-Atom Molecular Dynamics
# Simulation. Meth Enzymology. 578 (2016), 327-342,
# doi:10.1016/bs.mie.2016.05.024.
#

import io
from glob import glob
from os.path import (
    basename,
    dirname,
    join,
    splitext
)

from setuptools import (
    find_packages,
    setup
)

from pathlib import Path

version_ns = {}
init_py = Path(__file__).parent / "src" / "fluctmatch" / "__init__.py"
exec(init_py.read_text(), version_ns)

def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ).read()

setup(
    name="fluctmatch",
    version=version_ns["__version__"],
    license="BSD license",
    description="Structure Mechanics Statistical Learning of ENM parameters by fluctuation matching",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="Nixon Raj",
    author_email="nixon.raj@nationwidechildrens.org",
    url="https://github.com/nixnmtm/SMSL",
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Utilities",
    ],
    python_requires=">=3.9",
    keywords=[
        # eg: "keyword1", "keyword2", "keyword3",
    ],
    install_requires=[
        "click",
        "MDAnalysis",
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
    ],
    extras_require={
        # eg:
        #   "rst": ["docutils>=0.11"],
        #   ":python_version=="2.6"": ["argparse"],
    },
    entry_points={
        "console_scripts": [
            "fluctmatch = fluctmatch.cli:main",
        ]
    },
)
