# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from collections import OrderedDict

from MDAnalysis.core import topologyattrs
from future.builtins import (
    super,
    zip,
)

from fluctmatch.models.base import ModelBase


class Water(ModelBase):
    """Create a universe consisting of the water oxygen.
    """
    model = "WATER"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["OW"] = "water"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self.atoms.select_atoms("name OW").types = 1

    def _add_bonds(self):
        pass


class Tip3p(ModelBase):
    """Create a universe containing all three water atoms."""
    model = "TIP3P"
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


class Dma(ModelBase):
    """Create a universe for N-dimethylacetamide.
    """
    model = "DMA"
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
