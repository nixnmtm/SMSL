# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
    unicode_literals,
)

from MDAnalysis.coordinates.base import (_Readermeta, _Writermeta, IOBase)
from future.utils import (
    with_metaclass,
)


class TopologyReaderBase(with_metaclass(_Readermeta, IOBase)):
    def read(self):  # pragma: no cover
        """Read the file"""
        raise NotImplementedError("Override this in each subclass")


class TopologyWriterBase(with_metaclass(_Writermeta, IOBase)):
    def write(self, selection):  # pragma: no cover
        # type: (object) -> object
        """Write selection at current trajectory frame to file.

        Parameters
        ----------
        selection : AtomGroup
             group of atoms to be written

        """
        raise NotImplementedError("Override this in each subclass")
