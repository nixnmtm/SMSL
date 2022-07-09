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
    super,
    zip,
)

from collections import OrderedDict

from MDAnalysis.core import topologyattrs
from fluctmatch.models.base import ModelBase


class Water(ModelBase):
    """Create a universe consisting of the water oxygen.
    """
    model = "WATER"
    describe = "c.o.m./c.o.g. of whole water molecule"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["OW"] = "name OH2"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self.atoms.select_atoms("name OW").types = 1

    def _add_bonds(self):
        self._topology.add_TopologyAttr(topologyattrs.Bonds([]))
        self._generate_from_topology()


class Tip3p(ModelBase):
    """Create a universe containing all three water atoms."""
    model = "TIP3P"
    describe = "All-atom watter"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["OW"] = "name OW MW"
        self._mapping["HW1"] = "name HW1"
        self._mapping["HW2"] = "name HW2"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self.atoms.select_atoms("name OW").types = 1
        self.atoms.select_atoms("name HW1").types = 2
        self.atoms.select_atoms("name HW2").types = 3

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name OW").ix,
                s.atoms.select_atoms("name HW1").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name OW").ix,
                s.atoms.select_atoms("name HW2").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name HW1").ix,
                s.atoms.select_atoms("name HW2").ix)
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Dma(ModelBase):
    """Create a universe for N-dimethylacetamide.
    """
    model = "DMA"
    describe = "c.o.m./c.o.g. of C1, N, C2, and C3 of DMA"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["C1"] = "resname DMA and name C1 H1*"
        self._mapping["N"] = "resname DMA and name C N O"
        self._mapping["C2"] = "resname DMA and name C2 H2*"
        self._mapping["C3"] = "resname DMA and name C3 H3*"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self.atoms.select_atoms("name C1").types = 4
        self.atoms.select_atoms("name N").types = 5
        self.atoms.select_atoms("name C2").types = 6
        self.atoms.select_atoms("name C3").types = 7

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C1").ix,
                s.atoms.select_atoms("name N").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C2").ix,
                s.atoms.select_atoms("name N").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C3").ix,
                s.atoms.select_atoms("name N").ix)
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()
