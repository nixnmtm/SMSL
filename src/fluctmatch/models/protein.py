# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

from . import universe
from .selection import *


class Calpha(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["CA"] = "calpha"
        self._mapping["ions"] = "bioion"

        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)

        # Update the masses and charges
        ca_atu = self.atu.select_atoms("protein").split("residue")
        self.atoms.select_atoms("calpha").masses = np.array([_.total_mass() for _ in ca_atu])

        try:
            self.atoms.select_atoms("calpha").charges = np.array([_.total_charge() for _ in ca_atu])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("calpha").ix,
                                                            s.atoms.select_atoms("calpha").ix[1:])])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Caside(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["CA"] = "calpha"
        self._mapping["CB"] = "hsidechain"
        self._mapping["ions"] = "bioion"
        _back = "hbackbone"

        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)

        # Update the masses and charges
        ca_atu = self.atu.select_atoms(_back).split("residue")
        cb_atu = self.atu.select_atoms(self._mapping["CB"]).split("residue")
        self.atoms.select_atoms("calpha").masses = np.array([_.total_mass() for _ in ca_atu])
        self.atoms.select_atoms("cbeta").masses = np.array([_.total_mass() for _ in cb_atu])

        try:
            self.atoms.select_atoms("calpha").charges = np.array([_.total_charge() for _ in ca_atu])
        except AttributeError:
            pass
        try:
            self.atoms.select_atoms("cbeta").charges = np.array([_.total_charge() for _ in cb_atu])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("calpha").ix,
                                                            s.atoms.select_atoms("calpha").ix[1:])])
        bonds.extend([(r.atoms.select_atoms("calpha").ix[0], r.atoms.select_atoms("cbeta").ix[0]) for r in self.residues
                      if r.resname != "GLY" and r.resname in selection.ProteinSelection.prot_res])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Ncsc(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["N"] = "amine"
        self._mapping["CB"] = "hsidechain"
        self._mapping["C"] = "carboxyl"
        self._mapping["ions"] = "bioion"
        _back = "hcalpha"

        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)

        # Update the masses and charges
        ca_atu = self.atu.select_atoms(_back).split("residue")
        cb_atu = self.atu.select_atoms(self._mapping["CB"]).split("residue")
        self.atoms.select_atoms("name N").masses += 0.5 * np.array([_.total_mass() for _ in ca_atu])
        self.atoms.select_atoms("name C").masses += 0.5 * np.array([_.total_mass() for _ in ca_atu])
        self.atoms.select_atoms("cbeta").masses = np.array([_.total_mass() for _ in cb_atu])

        try:
            self.atoms.select_atoms("name N").charges += 0.5 * np.array([_.total_charge() for _ in ca_atu])
        except AttributeError:
            pass
        try:
            self.atoms.select_atoms("name C").charges += 0.5 * np.array([_.total_charge() for _ in ca_atu])
        except AttributeError:
            pass
        try:
            self.atoms.select_atoms("cbeta").charges = np.array([_.total_charge() for _ in cb_atu])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("name N").ix,
                                                            s.atoms.select_atoms("name C").ix)])
        bonds.extend([(r.atoms.select_atoms("name N").ix[0], r.atoms.select_atoms("cbeta").ix[0]) for r in self.residues
                      if r.resname != "GLY" and r.resname in selection.ProteinSelection.prot_res])
        bonds.extend([(r.atoms.select_atoms("cbeta").ix[0], r.atoms.select_atoms("name C").ix[0]) for r in self.residues
                      if r.resname != "GLY" and r.resname in selection.ProteinSelection.prot_res])
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("name C").ix,
                                                            s.atoms.select_atoms("name N").ix[1:])])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()
