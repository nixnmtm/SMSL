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
from fluctmatch.models.selection import *


class Nucleic3(ModelBase):
    """A universe consisting of the phosphate, sugar, and base of the nucleic acid.
    """
    model = "NUCLEIC3"
    describe = "Phosohate, sugar, and nucleotide of nucleic acid"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = "nucleicphosphate and not name H*"
        self._mapping["C4'"] = "hnucleicsugar and not name H*"
        self._mapping["C5"] = "hnucleicbase and not name H*"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_masses()
        self._set_chargess()

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name P").ix,
                s.atoms.select_atoms("name C4'").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix,
                s.atoms.select_atoms("name C5").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix[:-1],
                s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()

    def _set_masses(self):
        p_atu = self.atu.select_atoms("nucleicphosphate").split("residue")
        sugar_atu = self.atu.select_atoms("hnucleicsugar").split("residue")
        nucl_atu = self.atu.select_atoms("hnucleicbase").split("residue")

        self.atoms.select_atoms("name P").masses = np.array(
            [_.total_mass() for _ in p_atu])
        self.atoms.select_atoms("name C4'").masses = np.array(
            [_.total_mass() for _ in sugar_atu])
        self.atoms.select_atoms("name C5").masses = np.array(
            [_.total_mass() for _ in nucl_atu])

    def _set_chargess(self):
        p_atu = self.atu.select_atoms("nucleicphosphate").split("residue")
        sugar_atu = self.atu.select_atoms("hnucleicsugar").split("residue")
        nucl_atu = self.atu.select_atoms("hnucleicbase").split("residue")

        try:
            self.atoms.select_atoms("name P").charges = np.array(
                [_.total_charge() for _ in p_atu])
            self.atoms.select_atoms("name C4'").charges = np.array(
                [_.total_charge() for _ in sugar_atu])
            self.atoms.select_atoms("name C5").charges = np.array(
                [_.total_charge() for _ in nucl_atu])
        except AttributeError:
            pass


class Nucleic4(ModelBase):
    """A universe consisting of the phosphate, C4', C3', and base of the nucleic acid.
    """
    model = "NUCLEIC4"
    describe = "Phosphate, C2', C4', and c.o.m./c.o.g. of C4/C5 of nucleic acid"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = "nucleicphosphate and not name H*"
        self._mapping["C4'"] = "name C4'"
        self._mapping["C2'"] = "name C2'"
        self._mapping["C5"] = "nucleiccenter and not name H*"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_masses()
        self._set_chargess()

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name P").ix,
                s.atoms.select_atoms("name C4'").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix,
                s.atoms.select_atoms("name C2'").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix,
                s.atoms.select_atoms("name C5").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix[:-1],
                s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()

    def _set_masses(self):
        p_atu = self.atu.select_atoms("nucleicphosphate").split("residue")
        sugar_atu = self.atu.select_atoms("sugarC4").split("residue")
        sugar2_atu = self.atu.select_atoms("sugarC2").split("residue")
        nucl_atu = self.atu.select_atoms("hnucleicbase").split("residue")

        self.atoms.select_atoms("name P").masses = np.array(
            [_.total_mass() for _ in p_atu])
        self.atoms.select_atoms("name C4'").masses = np.array(
            [_.total_mass() for _ in sugar_atu])
        self.atoms.select_atoms("name C2'").masses = np.array(
            [_.total_mass() for _ in sugar2_atu])
        self.atoms.select_atoms("name C5").masses = np.array(
            [_.total_mass() for _ in nucl_atu])

    def _set_chargess(self):
        p_atu = self.atu.select_atoms("nucleicphosphate").split("residue")
        sugar_atu = self.atu.select_atoms("sugarC4").split("residue")
        sugar2_atu = self.atu.select_atoms("sugarC2").split("residue")
        nucl_atu = self.atu.select_atoms("hnucleicbase").split("residue")

        try:
            self.atoms.select_atoms("name P").charges = np.array(
                [_.total_charge() for _ in p_atu])
            self.atoms.select_atoms("name C4'").charges = np.array(
                [_.total_charge() for _ in sugar_atu])
            self.atoms.select_atoms("name C2'").charges = np.array(
                [_.total_charge() for _ in sugar2_atu])
            self.atoms.select_atoms("name C5").charges = np.array(
                [_.total_charge() for _ in nucl_atu])
        except AttributeError:
            pass


class Nucleic6(ModelBase):
    """A universe accounting for six sites involved with hydrogen bonding.
    """
    model = "NUCLEIC6"
    describe = "Phosphate, C2', C4', and 3 sites on the nucleotide"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = "name P H5T"
        self._mapping["C4'"] = "name C4'"
        self._mapping["C2'"] = "name C2'"
        self._mapping["H1"] = ("(resname ADE DA* RA* and name N6) or "
                               "(resname OXG GUA DG* RG* and name O6) or "
                               "(resname CYT DC* RC* and name N4) or "
                               "(resname THY URA DT* RU* and name O4)")
        self._mapping["H2"] = (
            "(resname ADE DA* RA* OXG GUA DG* RG* and name N1) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name N3)")
        self._mapping["H3"] = (
            "(resname ADE DA* RA* and name H2) or "
            "(resname OXG GUA DG* RG* and name N2) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name O2)")

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_chargess()
        self._set_masses()

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name P").ix,
                s.atoms.select_atoms("name C4'").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix,
                s.atoms.select_atoms("name C2'").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C2'").ix,
                s.atoms.select_atoms("name H1").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name H1").ix,
                s.atoms.select_atoms("name H2").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name H2").ix,
                s.atoms.select_atoms("name H3").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name C4'").ix[:-1],
                s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()

    def _set_chargess(self):
        self.atoms.charges = 0.

    def _set_masses(self):
        self.atoms.masses = self.atu.select_atoms("name C4'")[0].mass
