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
        self._set_charges()

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
            if s.residues.n_residues > 1
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

    def _set_charges(self):
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
        self._set_charges()

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

    def _set_charges(self):
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
        self._mapping["P"] = "name P O5'"
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
            "(resname ADE DA* RA* and name N3) or "
            "(resname OXG GUA DG* RG* and name N2) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name O2)")

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_charges()
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

    def _set_charges(self):
        self.atoms.charges = 0.

    def _set_masses(self):
        self.atoms.masses = 1.0  # Nix


# Nix Test
class DNA6S(ModelBase):
    """DNA 6 site with mass - A universe accounting for six sites involved with hydrogen bonding and mass determined from
    atoms selected for each site (not constant 1).
    """
    model = "DNA6S"
    describe = "Phosphate, RB, R, and 3 sites on the Base"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = ("(resname DC5 DA5 DT5 DG5 and name O5') or "
                              "(name P)")
        self._mapping["R"] = "name C1' C2' C3' C4' O4'"
        self._mapping["RB"] = (
            "(resname ADE DA* RA* OXG GUA DG* RG* and name N9) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name N1)")
        self._mapping["B1"] = ("(resname ADE DA* RA* and name N6) or "
                               "(resname OXG GUA DG* RG* and name O6) or "
                               "(resname CYT DC* RC* and name N4) or "
                               "(resname THY URA DT* RU* and name O4)")
        self._mapping["B2"] = (
            "(resname ADE DA* RA* OXG GUA DG* RG* and name N1) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name N3)")
        self._mapping["B3"] = (
            "(resname ADE DA* RA* and name C2) or "
            "(resname OXG GUA DG* RG* and name N2) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name O2)")

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_charges()
        self._set_masses()

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name P").ix,
                s.atoms.select_atoms("name R").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name R").ix,
                s.atoms.select_atoms("name RB").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name RB").ix,
                s.atoms.select_atoms("name B1").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name RB").ix,
                s.atoms.select_atoms("name B2").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name RB").ix,
                s.atoms.select_atoms("name B3").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name R").ix[:-1],
                s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()

    def _set_charges(self):
        self.atoms.charges = 0.

    def _set_masses(self):
        P_atu = self.atu.select_atoms(self._mapping["P"]).split("residue")
        R_atu = self.atu.select_atoms(self._mapping["R"]).split("residue")
        RB_atu = self.atu.select_atoms(self._mapping["RB"]).split("residue")
        B1_atu = self.atu.select_atoms(("(resname ADE DA* RA* and name N6) or "
                                        "(resname OXG GUA DG* RG* and name O6) or "
                                        "(resname CYT DC* RC* and name N4) or "
                                        "(resname THY URA DT* RU* and name O4)")).split("residue")
        B2_atu = self.atu.select_atoms(("(resname ADE DA* RA* OXG GUA DG* RG* and name N1) or "
                                        "(resname CYT DC* RC* THY URA DT* RU* and name N3)")).split("residue")
        B3_atu = self.atu.select_atoms(("(resname ADE DA* RA* and name C2) or "
                                        "(resname OXG GUA DG* RG* and name N2) or "
                                        "(resname CYT DC* RC* THY URA DT* RU* and name O2)")).split("residue")

        self.atoms.select_atoms("name P").masses = np.array(
            [_.total_mass() for _ in P_atu])
        self.atoms.select_atoms("name R").masses = np.array(
            [_.total_mass() for _ in R_atu])
        self.atoms.select_atoms("name RB").masses = np.array(
            [_.total_mass() for _ in RB_atu])
        self.atoms.select_atoms("name B1").masses = np.array(
            [_.total_mass() for _ in B1_atu])
        self.atoms.select_atoms("name B2").masses = np.array(
            [_.total_mass() for _ in B2_atu])
        self.atoms.select_atoms("name B3").masses = np.array(
            [_.total_mass() for _ in B3_atu])


class RNA7S(ModelBase):
    """RNA 7 site with mass - A universe accounting for seven sites involved with hydrogen bonding and mass determined
    from atoms selected for each site (not constant 1).
    """
    model = "RNA7S"
    describe = "Phosphate, R, RO, RB, and 3 sites on the Base"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["P"] = ("(resname DC5 DA5 DT5 DG5 and name O5') or "
                              "(name P)")
        self._mapping["R"] = "name C1' C2' C3' C4' O4'"
        self._mapping["RO"] = "name O2'"
        self._mapping["RB"] = (
            "(resname ADE DA* RA* OXG GUA DG* RG* and name N9) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name N1)")
        self._mapping["B1"] = ("(resname ADE DA* RA* and name N6) or "
                               "(resname OXG GUA DG* RG* and name O6) or "
                               "(resname CYT DC* RC* and name N4) or "
                               "(resname THY URA DT* RU* and name O4)")
        self._mapping["B2"] = (
            "(resname ADE DA* RA* OXG GUA DG* RG* and name N1) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name N3)")
        self._mapping["B3"] = (
            "(resname ADE DA* RA* and name C2) or "
            "(resname OXG GUA DG* RG* and name N2) or "
            "(resname CYT DC* RC* THY URA DT* RU* and name O2)")

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_charges()
        self._set_masses()

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name P").ix,
                s.atoms.select_atoms("name R").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name R").ix,
                s.atoms.select_atoms("name RB").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name R").ix,
                s.atoms.select_atoms("name RO").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name RB").ix,
                s.atoms.select_atoms("name B1").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name RB").ix,
                s.atoms.select_atoms("name B2").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name RB").ix,
                s.atoms.select_atoms("name B3").ix)
        ])
        bonds.extend([
            _ for s in self.segments for _ in zip(
                s.atoms.select_atoms("name R").ix[:-1],
                s.atoms.select_atoms("name P").ix[1:])
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()

    def _set_charges(self):
        self.atoms.charges = 0.

    def _set_masses(self):
        P_atu = self.atu.select_atoms(self._mapping["P"]).split("residue")
        R_atu = self.atu.select_atoms(self._mapping["R"]).split("residue")
        RO_atu = self.atu.select_atoms(self._mapping["RO"]).split("residue")
        RB_atu = self.atu.select_atoms(self._mapping["RB"]).split("residue")
        B1_atu = self.atu.select_atoms(("(resname ADE DA* RA* and name N6) or "
                                        "(resname OXG GUA DG* RG* and name O6) or "
                                        "(resname CYT DC* RC* and name N4) or "
                                        "(resname THY URA DT* RU* and name O4)")).split("residue")
        B2_atu = self.atu.select_atoms(("(resname ADE DA* RA* OXG GUA DG* RG* and name N1) or "
                                        "(resname CYT DC* RC* THY URA DT* RU* and name N3)")).split("residue")
        B3_atu = self.atu.select_atoms(("(resname ADE DA* RA* and name C2) or "
                                        "(resname OXG GUA DG* RG* and name N2) or "
                                        "(resname CYT DC* RC* THY URA DT* RU* and name O2)")).split("residue")

        self.atoms.select_atoms("name P").masses = np.array(
            [_.total_mass() for _ in P_atu])
        self.atoms.select_atoms("name R").masses = np.array(
            [_.total_mass() for _ in R_atu])
        self.atoms.select_atoms("name RO").masses = np.array(
            [_.total_mass() for _ in RO_atu])
        self.atoms.select_atoms("name RB").masses = np.array(
            [_.total_mass() for _ in RB_atu])
        self.atoms.select_atoms("name B1").masses = np.array(
            [_.total_mass() for _ in B1_atu])
        self.atoms.select_atoms("name B2").masses = np.array(
            [_.total_mass() for _ in B2_atu])
        self.atoms.select_atoms("name B3").masses = np.array(
            [_.total_mass() for _ in B3_atu])
