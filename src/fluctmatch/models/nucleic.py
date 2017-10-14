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

from collections import OrderedDict

from MDAnalysis.core import topologyattrs
from future.builtins import (
    super,
    zip,
)

from fluctmatch.models.base import ModelBase
from fluctmatch.models.selection import *


class Nucleic3(ModelBase):
    """A universe consisting of the phosphate, sugar, and base of the nucleic acid.
    """
    model = "NUCLEIC3"
    describe = "Phosohate, sugar, and nucleotide of nucleic acid"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = "nucleicphosphate"
        self._mapping["C4'"] = "hnucleicsugar"
        self._mapping["C5"] = "hnucleicbase"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name P").ix,
                         s.atoms.select_atoms("name C4'").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C4'").ix,
                         s.atoms.select_atoms("name C5").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C4'").ix[:-1],
                         s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Nucleic4(ModelBase):
    """A universe consisting of the phosphate, C4', C3', and base of the nucleic acid.
    """
    model = "NUCLEIC4"
    describe = "Phosphate, C3', C4', and c.o.m./c.o.g. of C4/C5 of nucleic acid"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = "nucleicphosphate"
        self._mapping["C4'"] = "sugarC4"
        self._mapping["C3'"] = "sugarC3"
        self._mapping["C5"] = "nucleiccenter"
        _nucl = "hnucleicbase"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)

        # Update the masses and charges
        nucl_base = self.atu.select_atoms(_nucl).split("residue")
        self.atoms.select_atoms("name C5").masses = np.array([
            _.total_mass()
            for _ in nucl_base
        ])

        try:
            self.atoms.select_atoms("name C5").charges = np.array([
                _.total_charge()
                for _ in nucl_base
            ])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name P").ix,
                         s.atoms.select_atoms("name C4'").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C4'").ix,
                         s.atoms.select_atoms("name C3'").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C4'").ix,
                         s.atoms.select_atoms("name C5").ix)
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(s.atoms.select_atoms("name C3'").ix[:-1],
                         s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()
