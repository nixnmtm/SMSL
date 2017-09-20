# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from collections import defaultdict

from future.builtins import (
    super,
)

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
