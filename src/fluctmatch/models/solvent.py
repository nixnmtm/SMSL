# -*- coding: utf-8 -*-
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
from . import universe
from .selection import *


class Water(universe._Universe):
    """Create a universe consisting of the water oxygen.
    """
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["OW"] = "segid SOL WAT and name OW HW* MW"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)
        self.atoms.select_atoms("name OW").types = 1

    def _add_bonds(self):
        pass


class Tip3p(universe._Universe):
    """Create a universe containing all three water atoms."""
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["OW"] = "segid SOL WAT and name OW MW"
        self._mapping["HW1"] = "segid SOL WAT and name HW1"
        self._mapping["HW2"] = "segid SOL WAT and name HW2"

        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)
        self.atoms.select_atoms("name OW").types = 1
        self.atoms.select_atoms("name HW1").types = 2
        self.atoms.select_atoms("name HW2").types = 3

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name OW").ix,
                         s.atoms.select_atoms("name HW1").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name OW").ix,
                         s.atoms.select_atoms("name HW2").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name HW1").ix,
                         s.atoms.select_atoms("name HW2").ix)
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Dma(universe._Universe):
    """Create a universe for N-dimethylacetamide.
    """
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["C1"] = "resname DMA and name C1 H1*"
        self._mapping["N"] = "resname DMA and name C N O"
        self._mapping["C2"] = "resname DMA and name C2 H2*"
        self._mapping["C3"] = "resname DMA and name C3 H3*"

        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)
        self.atoms.select_atoms("name C1").types = 4
        self.atoms.select_atoms("name N").types = 5
        self.atoms.select_atoms("name C2").types = 6
        self.atoms.select_atoms("name C3").types = 7

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C1").ix,
                         s.atoms.select_atoms("name N").ix)
        ])
        bonds.extend([
            _ for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C2").ix,
                         s.atoms.select_atoms("name N").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C3").ix,
                         s.atoms.select_atoms("name N").ix)
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()
