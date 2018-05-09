# -*- coding: utf-8 -*-
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
from future.builtins import (
    super, )

from collections import defaultdict

from MDAnalysis.core import (
    groups,
    topologyattrs,
)


class _Beads(topologyattrs.AtomAttr):
    """Underlying group of atoms for each bead"""
    attrname = "_beads"
    singular = "_bead"
    target_classes = [groups.Atom, groups.Residue, groups.Segment]
    per_object = "atom"
    transplants = defaultdict(list)

    def __init__(self, values, guessed=False):
        super().__init__(values, guessed)


class XplorTypes(topologyattrs.Atomtypes):
    """String types for atoms used for XPLOR-PSF."""
    attrname = "xplortypes"
    singular = "xplortype"
    target_classes = [groups.Atom, groups.Residue, groups.Segment]
    per_object = "atom"
    transplants = defaultdict(list)

    def __init__(self, values, guessed=False):
        super().__init__(values, guessed)


class NumTypes(topologyattrs.Atomtypes):
    """Number types for atoms used for PSF."""
    attrname = "numtypes"
    singular = "numtype"
    target_classes = [groups.Atom, groups.Residue, groups.Segment]
    per_object = "atom"
    transplants = defaultdict(list)

    def __init__(self, values, guessed=False):
        super().__init__(values, guessed)
