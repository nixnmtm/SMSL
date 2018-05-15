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
from __future__ import (
    absolute_import,
    division,
    print_function,
)

import io
import re
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


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ).read()


setup(
    name="fluctmatch",
    version="3.3.0",
    license="BSD license",
    description="Elastic network model using fluctuation matching.",
    long_description="%s\n%s" % (
        re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub("", read("README.rst")),
        re.sub(":[a-z]+:`~?(.*?)`", r"``\1``", read("CHANGELOG.rst"))
    ),
    author="Timothy Click",
    author_email="tclick@nctu.edu.tw",
    url="https://github.com/tclick/python-fluctmatch",
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        # uncomment if you test on these interpreters:
        # "Programming Language :: Python :: Implementation :: IronPython",
        # "Programming Language :: Python :: Implementation :: Jython",
        # "Programming Language :: Python :: Implementation :: Stackless",
        "Topic :: Utilities",
    ],
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
