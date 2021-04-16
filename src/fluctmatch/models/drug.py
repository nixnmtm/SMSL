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
from future.builtins import (
    super,
    zip,
)

from collections import OrderedDict

from MDAnalysis.core import topologyattrs
from fluctmatch.models.base import ModelBase


class PolarN3(ModelBase):
    """Create a universe consisting of N, C and polar heavy atoms in N3 biomimmetic drug.
    """
    model = "POLARN3"
    describe = "c.o.m./c.o.g. of N, C, and polar sidechains of protein"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["N1"] = "resname N3P and name N2"
        self._mapping["O1"] = "resname N3P and name O42"
        self._mapping["CB1"] = "resname N3P and name C6 O1"

        self._mapping["N2"] = "resname N3P and name N11"
        self._mapping["O2"] = "resname N3P and name O11"
        self._mapping["CB2"] = "resname N3P and name CB11"

        self._mapping["N3"] = "resname N3P and name N12"
        self._mapping["O3"] = "resname N3P and name O12"
        self._mapping["CB3"] = "resname N3P and name CG1 CG2"

        self._mapping["N4"] = "resname N3P and name N13"
        self._mapping["O4"] = "resname N3P and name O13"
        self._mapping["CB4"] = "resname N3P and name CD1 CD2"

        self._mapping["N5"] = "resname N3P and name N5"
        self._mapping["O5"] = "resname N3P and name O7"
        self._mapping["CB5"] = "resname N3P and name C26 C27 C28 C29 N6 O8"

        self._mapping["O6"] = "resname N3P and name O"
        self._mapping["CB6"] = "resname N3P and name C1 C2 C7 C8 C9 C10"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_masses()
        self._set_charges()

    def _add_bonds(self):
        bonds = []

        for i in range(1, 6):

            if i == 6:
                bonds.extend([
                    _ for s in self.segments for _ in zip(
                        s.atoms.select_atoms(f"name O6").ix,
                        s.atoms.select_atoms(f"name CB6").ix)
                ])
            else:
                bonds.extend([
                    _ for s in self.segments for _ in zip(
                        s.atoms.select_atoms(f"name N{i}").ix,
                        s.atoms.select_atoms(f"name O{i}").ix)
                ])
                bonds.extend([
                    _ for s in self.segments for _ in zip(
                        s.atoms.select_atoms(f"name N{i}").ix,
                        s.atoms.select_atoms(f"name CB{i}").ix)
                ])
                bonds.extend([
                    _ for s in self.segments for _ in zip(
                        s.atoms.select_atoms(f"name O{i}").ix,
                        s.atoms.select_atoms(f"name CB{i}").ix)
                ])
            if i < 5:
                bonds.extend([
                    _ for s in self.segments for _ in zip(
                        s.atoms.select_atoms(f"name O{i}").ix,
                        s.atoms.select_atoms(f"name N{i+1}").ix)
                ])
            if i == 5:
                bonds.extend([
                    _ for s in self.segments for _ in zip(
                        s.atoms.select_atoms(f"name O{i}").ix,
                        s.atoms.select_atoms(f"name O{i + 1}").ix)
                ])

        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()

    def _set_masses(self):
        # sharing the mass of atom CA to atoms O and N as in protein Polar CG model for amino acids in N3 drug
        CA_map = OrderedDict()
        CA_map["CA1"] = 0.
        CA_map["CA2"] = 0.5 * self.atu.select_atoms("resname N3P and name CA11").total_mass()
        CA_map["CA3"] = 0.5 * self.atu.select_atoms("resname N3P and name CA12").total_mass()
        CA_map["CA4"] = 0.5 * self.atu.select_atoms("resname N3P and name CA13").total_mass()
        CA_map["CA5"] = 0.

        for i in range(1, 7):
            if i == 6:
                self.atoms.select_atoms(f"name O{i}").masses = self.atu.select_atoms(
                    self._mapping[f"O{i}"]).total_mass()/self.atu.segments.n_segments
                self.atoms.select_atoms(f"name CB{i}").masses = self.atu.select_atoms(
                    self._mapping[f"CB{i}"]).total_mass()/self.atu.segments.n_segments
            else:
                self.atoms.select_atoms(f"name N{i}").masses = (self.atu.select_atoms(
                    self._mapping[f"N{i}"]).total_mass() + CA_map[f"CA{i}"])/self.atu.segments.n_segments
                self.atoms.select_atoms(f"name O{i}").masses = (self.atu.select_atoms(
                    self._mapping[f"O{i}"]).total_mass() + CA_map[f"CA{i}"])/self.atu.segments.n_segments
                self.atoms.select_atoms(f"name CB{i}").masses = self.atu.select_atoms(
                    self._mapping[f"CB{i}"]).total_mass()/self.atu.segments.n_segments

    def _set_charges(self):
        # sharing the mass of atom CA to atoms O and N as in protein Polar CG model for amino acids in N3 drug
        CA_map = OrderedDict()
        CA_map["CA1"] = 0.
        CA_map["CA2"] = 0.5 * self.atu.select_atoms("resname N3P and name CA11").total_charge()
        CA_map["CA3"] = 0.5 * self.atu.select_atoms("resname N3P and name CA12").total_charge()
        CA_map["CA4"] = 0.5 * self.atu.select_atoms("resname N3P and name CA13").total_charge()
        CA_map["CA5"] = 0.

        for i in range(1, 7):
            if i == 6:
                self.atoms.select_atoms(f"name O{i}").charges = self.atu.select_atoms(
                    self._mapping[f"O{i}"]).total_charge()/self.atu.segments.n_segments
                self.atoms.select_atoms(f"name CB{i}").charges = self.atu.select_atoms(
                    self._mapping[f"CB{i}"]).total_charge()/self.atu.segments.n_segments
            else:
                self.atoms.select_atoms(f"name N{i}").charges = (self.atu.select_atoms(
                    self._mapping[f"N{i}"]).total_charge() + CA_map[f"CA{i}"])/self.atu.segments.n_segments
                self.atoms.select_atoms(f"name O{i}").charges = (self.atu.select_atoms(
                    self._mapping[f"O{i}"]).total_charge() + CA_map[f"CA{i}"])/self.atu.segments.n_segments
                self.atoms.select_atoms(f"name CB{i}").charges = self.atu.select_atoms(
                    self._mapping[f"CB{i}"]).total_charge()/self.atu.segments.n_segments
