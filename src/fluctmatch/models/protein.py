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

import itertools
from collections import OrderedDict

from MDAnalysis.core import (
    topology,
    topologyattrs,
)
from MDAnalysis.topology import base as topbase
from future.builtins import (
    dict,
    zip,
)
from future.utils import viewitems

from fluctmatch.models import topattrs
from fluctmatch.models.base import ModelBase
from fluctmatch.models.selection import *


class Calpha(ModelBase):
    """Create a universe defined by the protein C-alpha.
    """
    model = "CALPHA"
    describe = "C-alpha of a protein"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["CA"] = "calpha"
        self._mapping["ions"] = "bioion"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)

        # Update the masses and charges
        ca_atu = self.atu.select_atoms("protein").split("residue")
        self.atoms.select_atoms("calpha").masses = np.array([
            _.total_mass()
            for _ in ca_atu
        ])

        try:
            self.atoms.select_atoms("calpha").charges = np.array([
                _.total_charge()
                for _ in ca_atu
            ])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(
                s.atoms.select_atoms("calpha").ix,
                s.atoms.select_atoms("calpha").ix[1:]
            )
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Caside(ModelBase):
    """Create a universe consisting of the C-alpha and sidechains of a protein.
    """
    model = "CASIDE"
    describe = "C-alpha and sidechain (c.o.m./c.o.g.) of protein"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["CA"] = "calpha"
        self._mapping["CB"] = "hsidechain and not name H*"
        self._mapping["ions"] = "bioion"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)

        # Update the masses and charges
        ca_atu = self.atu.select_atoms("hbackbone").split("residue")
        cb_atu = self.atu.select_atoms("hsidechain").split("residue")
        self.atoms.select_atoms("calpha").masses = np.array([
            _.total_mass()
            for _ in ca_atu
        ])
        self.atoms.select_atoms("cbeta").masses = np.array([
            _.total_mass()
            for _ in cb_atu
        ])

        try:
            self.atoms.select_atoms("calpha").charges = np.array([
                _.total_charge()
                for _ in ca_atu
            ])
        except AttributeError:
            pass
        try:
            self.atoms.select_atoms("cbeta").charges = np.array([
                _.total_charge()
                for _ in cb_atu
            ])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(
                s.atoms.select_atoms("calpha").ix,
                s.atoms.select_atoms("calpha").ix[1:]
            )
        ])
        bonds.extend([
            (
                r.atoms.select_atoms("calpha").ix[0],
                r.atoms.select_atoms("cbeta").ix[0]
            )
            for r in self.residues
            if (
                r.resname != "GLY" and
                r.resname in selection.ProteinSelection.prot_res
            )
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Ncsc(ModelBase):
    """Create a universe consisting of the amine, carboxyl, and sidechain regions.
    """
    model = "NCSC"
    describe = "c.o.m./c.o.g. of N, O, and sidechain of protein"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["N"] = "protein and name N"
        self._mapping["CB"] = "hsidechain and not name H*"
        self._mapping["O"] = "protein and name O OT1 OT2 OXT"
        self._mapping["ions"] = "bioion"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)

        # Update the masses and charges
        n_atu = self.atu.select_atoms("amine").split("residue")
        o_atu = self.atu.select_atoms("carboxyl").split("residue")
        ca_atu = self.atu.select_atoms("hcalpha").split("residue")
        cb_atu = self.atu.select_atoms("hsidechain").split("residue")

        ca_masses = 0.5 * np.array([_.total_mass() for _ in ca_atu])
        self.atoms.select_atoms("name N").masses = np.array([
            _.total_mass()
            for _ in n_atu
        ]) + ca_masses
        self.atoms.select_atoms("cbeta").masses = np.array([
            _.total_mass()
            for _ in cb_atu
        ])
        self.atoms.select_atoms("name O").masses = np.array([
            _.total_mass()
            for _ in o_atu
        ]) + ca_masses

        try:
            ca_charges = 0.5 * np.array([_.total_charge() for _ in ca_atu])
            self.atoms.select_atoms("name N").charges = np.array([
                _.total_mass()
                for _ in n_atu
            ]) + ca_charges
            self.atoms.select_atoms("name O").charges = np.array([
                _.total_charge()
                for _ in o_atu
            ]) + ca_charges
            self.atoms.select_atoms("cbeta").charges = np.array([
                _.total_charge()
                for _ in cb_atu
            ])
        except AttributeError:
            pass

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(
                s.atoms.select_atoms("name N").ix,
                s.atoms.select_atoms("name O").ix
            )
        ])
        bonds.extend([
            (
                r.atoms.select_atoms("name N").ix[0],
                r.atoms.select_atoms("cbeta").ix[0]
            )
            for r in self.residues
            if (
                r.resname != "GLY" and
                r.resname in selection.ProteinSelection.prot_res
            )
        ])
        bonds.extend([
            (
                r.atoms.select_atoms("cbeta").ix[0],
                r.atoms.select_atoms("name O").ix[0]
            )
            for r in self.residues
            if (
                r.resname != "GLY" and
                r.resname in selection.ProteinSelection.prot_res
            )
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(
                s.atoms.select_atoms("name C").ix,
                s.atoms.select_atoms("name N").ix[1:]
            )
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()


class Polar(ModelBase):
    """Create a universe consisting of the amine, carboxyl, and polar regions.
    """
    model = "POLAR"
    describe = "c.o.m./c.o.g. of N, C, and polar sidechains of protein"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["N"] = "protein and name N"
        self._mapping["CB"] = dict(
            ALA="name CB",
            ARG="name NH*",
            ASN="name OD1 ND2",
            ASP="name OD*",
            CYS="name SG",
            GLN="name OE1 NE2",
            GLU="name OE*",
            HIS="name CG ND1 CD2 CE1 NE2",
            HSD="name CG ND1 CD2 CE1 NE2",
            HSE="name CG ND1 CD2 CE1 NE2",
            HSP="name CG ND1 CD2 CE1 NE2",
            ILE="name CG1 CG2 CD",
            LEU="name CD1 CD2",
            LYS="name NZ",
            MET="name SD",
            PHE="name CG CD* CE* CZ",
            PRO="name CG",
            SER="name OG",
            THR="name OG1",
            TRP="name CG CD* NE CE* CZ* CH",
            TYR="name CG CD* CE* CZ OH",
            VAL="name CG1 CG2",
        )
        self._mapping["O"] = "protein and name O OT1 OT2 OXT"
        self._mapping["ions"] = "bioion"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)

        # Update the masses and charges
        n_atu = self.atu.select_atoms("amine").split("residue")
        o_atu = self.atu.select_atoms("carboxyl").split("residue")
        ca_atu = self.atu.select_atoms("hcalpha").split("residue")
        cb_atu = self.atu.select_atoms("hsidechain").split("residue")

        ca_masses = 0.5 * np.array([_.total_mass() for _ in ca_atu])
        self.atoms.select_atoms("name N").masses = np.array([
            _.total_mass()
            for _ in n_atu
        ]) + ca_masses
        self.atoms.select_atoms("cbeta").masses = np.array([
            _.total_mass()
            for _ in cb_atu
        ])
        self.atoms.select_atoms("name O").masses = np.array([
            _.total_mass()
            for _ in o_atu
        ]) + ca_masses

        try:
            ca_charges = 0.5 * np.array([_.total_charge() for _ in ca_atu])
            self.atoms.select_atoms("name N").charges = np.array([
                _.total_mass()
                for _ in n_atu
            ]) + ca_charges
            self.atoms.select_atoms("name O").charges = np.array([
                _.total_charge()
                for _ in o_atu
            ]) + ca_charges
            self.atoms.select_atoms("cbeta").charges = np.array([
                _.total_charge()
                for _ in cb_atu
            ])
        except AttributeError:
            pass

    def _apply_map(self, mapping):
        """Apply the mapping scheme to the beads.

        Parameters
        ----------
        mapping : dict
            Mapping definitions per bead/

        Returns
        -------
        :class:`~MDAnalysis.core.topology.Topology` defining the new universe.
        """
        # Allocate arrays
        _beads = []
        atomnames = []
        atomids = []
        resids = []
        resnames = []
        segids = []
        charges = []
        masses = []

        residues = self.atu.atoms.split("residue")
        select_residues = enumerate(
            itertools.product(residues, viewitems(mapping))
        )
        for i, (res, (name, selection)) in select_residues:
            if name != "CB":
                bead = res.select_atoms(selection)
            else:
                bead = res.select_atoms(
                    selection.get(res, "hsidechain and not name H*")
                )
            if bead:
                _beads.append(bead)
                atomnames.append(name)
                atomids.append(i)
                resids.append(bead.resids[0])
                resnames.append(bead.resnames[0])
                segids.append(bead.segids[0].split("_")[-1])
                try:
                    charges.append(bead.total_charge())
                except AttributeError:
                    charges.append(0.)
                masses.append(bead.total_mass())

        _beads = np.array(_beads)
        n_atoms = len(_beads)

        # Atom
        _beads = topattrs._Beads(_beads)
        vdwradii = np.zeros_like(atomids)
        vdwradii = topologyattrs.Radii(vdwradii)
        atomids = topologyattrs.Atomids(np.asarray(atomids))
        atomnames = topologyattrs.Atomnames(
            np.asarray(atomnames, dtype=np.object)
        )
        atomtypes = topologyattrs.Atomtypes(np.asarray(np.arange(n_atoms)+100))
        charges = topologyattrs.Charges(np.asarray(charges))
        masses = topologyattrs.Masses(np.asarray(masses))

        # Residue
        # resids, resnames
        segids = np.asarray(segids, dtype=np.object)
        resids = np.asarray(resids)
        resnames = np.asarray(resnames, dtype=np.object)
        residx, (new_resids, new_resnames, new_segids) = topbase.change_squash(
            (resids,), (resids, resnames, segids)
        )

        # transform from atom:Rid to atom:Rix
        residueids = topologyattrs.Resids(new_resids)
        residuenums = topologyattrs.Resnums(new_resids.copy())
        residuenames = topologyattrs.Resnames(new_resnames)

        # Segment
        segidx, (perseg_segids,) = topbase.change_squash(
            (new_segids,), (new_segids,)
        )
        segids = topologyattrs.Segids(perseg_segids)

        # Setup topology
        top = topology.Topology(len(atomids), len(new_resids), len(segids),
                                attrs=[_beads, atomids, atomnames, atomtypes,
                                       charges, masses, vdwradii, residueids,
                                       residuenums, residuenames, segids],
                                atom_resindex=residx,
                                residue_segindex=segidx)
        return top

    def _add_bonds(self):
        bonds = []
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(
                s.atoms.select_atoms("name N").ix,
                s.atoms.select_atoms("name O").ix
            )
        ])
        bonds.extend([
            (
                r.atoms.select_atoms("name N").ix[0],
                r.atoms.select_atoms("cbeta").ix[0]
            )
            for r in self.residues
            if (
                r.resname != "GLY" and
                r.resname in selection.ProteinSelection.prot_res
            )
        ])
        bonds.extend([
            (
                r.atoms.select_atoms("cbeta").ix[0],
                r.atoms.select_atoms("name O").ix[0]
            )
            for r in self.residues
            if (
                r.resname != "GLY" and
                r.resname in selection.ProteinSelection.prot_res
            )
        ])
        bonds.extend([
            _
            for s in self.segments
            for _ in zip(
                s.atoms.select_atoms("name C").ix,
                s.atoms.select_atoms("name N").ix[1:]
            )
        ])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()
