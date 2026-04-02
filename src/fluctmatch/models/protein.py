# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# SMSL - https://github.com/nixnmtm/SMSL
#
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.builtins import dict, zip
from collections import OrderedDict

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topology, topologyattrs
from MDAnalysis.topology import base as topbase

from fluctmatch.models.base import ModelBase
from fluctmatch.models.selection import *


def _safe_total_charge(ag):
    try:
        return ag.total_charge()
    except AttributeError:
        return 0.0


class Calpha(ModelBase):
    """Create a universe defined by the protein C-alpha."""
    model = "CALPHA"
    describe = "C-alpha of a protein"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["CA"] = "calpha"
        self._mapping["ions"] = "bioion"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_masses()
        self._set_charges()

    def _add_bonds(self):
        bonds = []
        bonds.extend(
            [
                _
                for s in self.segments
                for _ in zip(
                    s.atoms.select_atoms("name CA").ix,
                    s.atoms.select_atoms("name CA").ix[1:],
                )
            ]
        )
        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)

    def _set_masses(self):
        ca_masses = [
            res.atoms.total_mass()
            for res in self.atu.select_atoms("protein").residues
            if res.atoms.select_atoms("name CA").n_atoms > 0
        ]
        self.atoms.select_atoms("name CA").masses = np.asarray(ca_masses)

        ion_masses = [
            ion.total_mass()
            for ion in self.atu.select_atoms("bioion").split("residue")
            if ion.n_atoms > 0
        ]
        if self.atoms.select_atoms("name ions").n_atoms > 0:
            self.atoms.select_atoms("name ions").masses = np.asarray(ion_masses)

    def _set_charges(self):
        try:
            ca_charges = [
                res.atoms.total_charge()
                for res in self.atu.select_atoms("protein").residues
                if res.atoms.select_atoms("name CA").n_atoms > 0
            ]
            self.atoms.select_atoms("name CA").charges = np.asarray(ca_charges)

            ion_charges = [
                ion.total_charge()
                for ion in self.atu.select_atoms("bioion").split("residue")
                if ion.n_atoms > 0
            ]
            if self.atoms.select_atoms("name ions").n_atoms > 0:
                self.atoms.select_atoms("name ions").charges = np.asarray(ion_charges)
        except AttributeError:
            pass


class Caside(ModelBase):
    """Create a universe consisting of the C-alpha and sidechains of a protein."""
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
        self._set_masses()
        self._set_charges()

    def _add_bonds(self):
        bonds = []
        bonds.extend(
            [
                _
                for s in self.segments
                for _ in zip(
                    s.atoms.select_atoms("name CA").ix,
                    s.atoms.select_atoms("name CA").ix[1:],
                )
            ]
        )
        bonds.extend(
            [
                (
                    r.atoms.select_atoms("name CA").ix[0],
                    r.atoms.select_atoms("name CB").ix[0],
                )
                for r in self.residues
                if r.atoms.select_atoms("name CA").n_atoms > 0
                and r.atoms.select_atoms("name CB").n_atoms > 0
            ]
        )
        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)

    def _set_masses(self):
        ca_masses = []
        cb_masses = []

        for r in self.atu.select_atoms("protein").residues:
            ca_bead = r.atoms.select_atoms("name CA")
            cb_bead = r.atoms.select_atoms("name CB")

            bb = r.atoms.select_atoms("hbackbone")
            sc = r.atoms.select_atoms("hsidechain")

            if ca_bead.n_atoms > 0 and bb.n_atoms > 0:
                ca_masses.append(bb.total_mass())

            if cb_bead.n_atoms > 0 and sc.n_atoms > 0:
                cb_masses.append(sc.total_mass())

        self.atoms.select_atoms("name CA").masses = np.asarray(ca_masses)
        self.atoms.select_atoms("name CB").masses = np.asarray(cb_masses)

        ion_masses = [
            ion.total_mass()
            for ion in self.atu.select_atoms("bioion").split("residue")
            if ion.n_atoms > 0
        ]
        if self.atoms.select_atoms("name ions").n_atoms > 0:
            self.atoms.select_atoms("name ions").masses = np.asarray(ion_masses)

    def _set_charges(self):
        try:
            ca_charges = []
            cb_charges = []

            for r in self.atu.select_atoms("protein").residues:
                ca_bead = r.atoms.select_atoms("name CA")
                cb_bead = r.atoms.select_atoms("name CB")

                bb = r.atoms.select_atoms("hbackbone")
                sc = r.atoms.select_atoms("hsidechain")

                if ca_bead.n_atoms > 0 and bb.n_atoms > 0:
                    ca_charges.append(bb.total_charge())

                if cb_bead.n_atoms > 0 and sc.n_atoms > 0:
                    cb_charges.append(sc.total_charge())

            self.atoms.select_atoms("name CA").charges = np.asarray(ca_charges)
            self.atoms.select_atoms("name CB").charges = np.asarray(cb_charges)

            ion_charges = [
                ion.total_charge()
                for ion in self.atu.select_atoms("bioion").split("residue")
                if ion.n_atoms > 0
            ]
            if self.atoms.select_atoms("name ions").n_atoms > 0:
                self.atoms.select_atoms("name ions").charges = np.asarray(ion_charges)
        except AttributeError:
            pass


class Polar(ModelBase):
    """Create a universe consisting of the amine, carboxyl, and polar regions."""
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
        self._set_masses()
        self._set_charges()

    def _cb_selection(self, resname):
        return self._mapping["CB"].get(resname, None)

    def _apply_map(self, mapping):
        beads = []
        atomnames = []
        atomids = []
        resids = []
        resnames = []
        segids = []
        charges = []
        masses = []

        i = 0

        for res in self.atu.select_atoms("protein").residues:
            # N bead
            n_bead = res.atoms.select_atoms(mapping["N"])
            if n_bead.n_atoms > 0:
                beads.append(n_bead)
                atomnames.append("N")
                atomids.append(i)
                resids.append(n_bead.resids[0])
                resnames.append(n_bead.resnames[0])
                segids.append(n_bead.segids[0].split("_")[-1])
                charges.append(0.0)
                masses.append(0.0)
                i += 1

            # CB bead
            cb_sel = self._cb_selection(res.resname)
            if cb_sel is not None:
                cb_bead = res.atoms.select_atoms(cb_sel)
                if cb_bead.n_atoms > 0:
                    beads.append(cb_bead)
                    atomnames.append("CB")
                    atomids.append(i)
                    resids.append(cb_bead.resids[0])
                    resnames.append(cb_bead.resnames[0])
                    segids.append(cb_bead.segids[0].split("_")[-1])
                    charges.append(0.0)
                    masses.append(0.0)
                    i += 1

            # O bead
            o_bead = res.atoms.select_atoms(mapping["O"])
            if o_bead.n_atoms > 0:
                beads.append(o_bead)
                atomnames.append("O")
                atomids.append(i)
                resids.append(o_bead.resids[0])
                resnames.append(o_bead.resnames[0])
                segids.append(o_bead.segids[0].split("_")[-1])
                charges.append(0.0)
                masses.append(0.0)
                i += 1

        # ions
        for ion in self.atu.select_atoms(mapping["ions"]).split("residue"):
            if ion.n_atoms > 0:
                beads.append(ion)
                atomnames.append("ion")
                atomids.append(i)
                resids.append(ion.resids[0])
                resnames.append(ion.resnames[0])
                segids.append(ion.segids[0].split("_")[-1])
                charges.append(0.0)
                masses.append(0.0)
                i += 1

        beads = np.asarray(beads, dtype=object)
        n_atoms = len(beads)

        vdwradii = topologyattrs.Radii(np.zeros(n_atoms, dtype=float))
        atomids = topologyattrs.Atomids(np.asarray(atomids))
        atomnames = topologyattrs.Atomnames(np.asarray(atomnames, dtype=object))
        atomtypes = topologyattrs.Atomtypes(np.asarray(np.arange(n_atoms) + 100))
        charges = topologyattrs.Charges(np.asarray(charges, dtype=float))
        masses = topologyattrs.Masses(np.asarray(masses, dtype=float))

        segids = np.asarray(segids, dtype=object)
        resids = np.asarray(resids)
        resnames = np.asarray(resnames, dtype=object)

        residx, (new_resids, new_resnames, new_segids) = topbase.change_squash(
            (resids,), (resids, resnames, segids)
        )

        residueids = topologyattrs.Resids(new_resids)
        residuenums = topologyattrs.Resnums(new_resids.copy())
        residuenames = topologyattrs.Resnames(new_resnames)

        segidx, (perseg_segids,) = topbase.change_squash((new_segids,), (new_segids,))
        segids = topologyattrs.Segids(perseg_segids)

        top = topology.Topology(
            len(atomids),
            len(new_resids),
            len(segids),
            attrs=[
                atomids,
                atomnames,
                atomtypes,
                charges,
                masses,
                vdwradii,
                residueids,
                residuenums,
                residuenames,
                segids,
            ],
            atom_resindex=residx,
            residue_segindex=segidx,
        )
        return top

    def _add_bonds(self):
        bonds = []

        for s in self.segments:
            n_ix = s.atoms.select_atoms("name N").ix
            o_ix = s.atoms.select_atoms("name O").ix
            bonds.extend(zip(n_ix, o_ix))
            if len(n_ix) > 1 and len(o_ix) > 1:
                bonds.extend(zip(o_ix[:-1], n_ix[1:]))

        for r in self.residues:
            n = r.atoms.select_atoms("name N")
            cb = r.atoms.select_atoms("name CB")
            o = r.atoms.select_atoms("name O")

            if n.n_atoms > 0 and cb.n_atoms > 0:
                bonds.append((n.ix[0], cb.ix[0]))
            if cb.n_atoms > 0 and o.n_atoms > 0:
                bonds.append((cb.ix[0], o.ix[0]))

        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)

    def _set_masses(self):
        n_masses = []
        cb_masses = []
        o_masses = []

        for r in self.atu.select_atoms("protein").residues:
            n_bead = r.atoms.select_atoms("protein and name N")
            o_bead = r.atoms.select_atoms("protein and name O OT1 OT2 OXT")

            n = r.atoms.select_atoms("amine")
            ca = r.atoms.select_atoms("hcalpha")
            o = r.atoms.select_atoms("carboxyl")

            cb_sel = self._cb_selection(r.resname)
            cb = r.atoms.select_atoms(cb_sel) if cb_sel is not None else r.atoms[:0]

            if n_bead.n_atoms > 0 and n.n_atoms > 0 and ca.n_atoms > 0:
                n_masses.append(n.total_mass() + 0.5 * ca.total_mass())

            if cb.n_atoms > 0:
                cb_masses.append(cb.total_mass())

            if o_bead.n_atoms > 0 and o.n_atoms > 0 and ca.n_atoms > 0:
                o_masses.append(o.total_mass() + 0.5 * ca.total_mass())

        self.atoms.select_atoms("name N").masses = np.asarray(n_masses)
        self.atoms.select_atoms("name CB").masses = np.asarray(cb_masses)
        self.atoms.select_atoms("name O").masses = np.asarray(o_masses)

        ion_masses = [
            ion.total_mass()
            for ion in self.atu.select_atoms("bioion").split("residue")
            if ion.n_atoms > 0
        ]
        if self.atoms.select_atoms("name ion").n_atoms > 0:
            self.atoms.select_atoms("name ion").masses = np.asarray(ion_masses)

    def _set_charges(self):
        try:
            n_charges = []
            cb_charges = []
            o_charges = []

            for r in self.atu.select_atoms("protein").residues:
                n_bead = r.atoms.select_atoms("protein and name N")
                o_bead = r.atoms.select_atoms("protein and name O OT1 OT2 OXT")

                n = r.atoms.select_atoms("amine")
                ca = r.atoms.select_atoms("hcalpha")
                o = r.atoms.select_atoms("carboxyl")

                cb_sel = self._cb_selection(r.resname)
                cb = r.atoms.select_atoms(cb_sel) if cb_sel is not None else r.atoms[:0]

                if n_bead.n_atoms > 0 and n.n_atoms > 0 and ca.n_atoms > 0:
                    n_charges.append(n.total_charge() + 0.5 * ca.total_charge())

                if cb.n_atoms > 0:
                    cb_charges.append(cb.total_charge())

                if o_bead.n_atoms > 0 and o.n_atoms > 0 and ca.n_atoms > 0:
                    o_charges.append(o.total_charge() + 0.5 * ca.total_charge())

            self.atoms.select_atoms("name N").charges = np.asarray(n_charges)
            self.atoms.select_atoms("name CB").charges = np.asarray(cb_charges)
            self.atoms.select_atoms("name O").charges = np.asarray(o_charges)

            ion_charges = [
                ion.total_charge()
                for ion in self.atu.select_atoms("bioion").split("residue")
                if ion.n_atoms > 0
            ]
            if self.atoms.select_atoms("name ion").n_atoms > 0:
                self.atoms.select_atoms("name ion").charges = np.asarray(ion_charges)
        except AttributeError:
            pass