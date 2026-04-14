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

def _trueprotein_residues(u):
    return list(u.select_atoms("trueprotein").residues)

def _cap_before(u, res):
    idx = res.resindex
    if idx > 0:
        prev_res = u.residues[idx - 1]
        if prev_res.resname == "ACE":
            return prev_res
    return None

def _cap_after(u, res):
    idx = res.resindex
    if idx < len(u.residues) - 1:
        next_res = u.residues[idx + 1]
        if next_res.resname == "NME":
            return next_res
    return None

def _alpha_code(i):
    """0 -> A, 1 -> B, ..., 25 -> Z, 26 -> AA, ..."""
    letters = []
    i = int(i)
    while True:
        i, rem = divmod(i, 26)
        letters.append(chr(ord("A") + rem))
        if i == 0:
            break
        i -= 1
    return "".join(reversed(letters))

def _first_molnum(bead):
    try:
        molnums = bead.atoms.molnums
        if molnums is not None and len(molnums) > 0:
            return int(molnums[0])
    except Exception:
        pass
    return None

def _first_chainid(bead):
    try:
        chainids = bead.atoms.chainIDs
        if chainids is not None and len(chainids) > 0:
            chainid = str(chainids[0]).strip()
            if chainid:
                return chainid
    except Exception:
        pass
    return None

def _protein_group_key(bead):
    """
    Prefer molnum for GROMACS/TPR systems.
    Fall back to short chainID only if it looks real.
    """
    molnum = _first_molnum(bead)
    if molnum is not None:
        return ("molnum", molnum)

    chainid = _first_chainid(bead)
    if chainid and len(chainid) <= 2 and " " not in chainid and chainid.upper() not in {"P", "PROTEIN"}:
        return ("chain", chainid)

    try:
        segid = str(bead.segids[0]).strip()
        if segid and " " not in segid and len(segid) <= 4 and segid.upper() not in {"PROT", "SYS"}:
            return ("segid", segid)
    except Exception:
        pass

    return ("fallback", "protein")

def _ion_group_key(bead):
    molnum = _first_molnum(bead)
    if molnum is not None:
        return ("ion_molnum", molnum)

    try:
        return ("ion_resid", int(bead.resids[0]))
    except Exception:
        return ("ion_fallback", "ion")

def _safe_segid(name, bead, segid_map=None, default="SYS"):
    """
    Return 4-character PSF-safe SEGIDs.

    Protein beads: PROA, PROB, ...
    Ion beads:     IONA, IONB, ...
    Ligands:       residue-name based, e.g. OM, ADP
    """
    if segid_map is None:
        segid_map = {}

    # proteins
    if name in {"CA", "CB", "N", "O"}:
        key = _protein_group_key(bead)
        if key not in segid_map:
            segid_map[key] = f"PRO{_alpha_code(sum(v.startswith('PRO') for v in segid_map.values()))}"[:4]
        return segid_map[key]

    # ions
    if name in {"ion", "ions"}:
        key = _ion_group_key(bead)
        if key not in segid_map:
            segid_map[key] = f"ION{_alpha_code(sum(v.startswith('ION') for v in segid_map.values()))}"[:4]
        return segid_map[key]

    # ligands and others
    try:
        resname = str(bead.resnames[0]).strip().upper()
        if resname:
            return resname[:4]
    except Exception:
        pass

    return default[:4]

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

    def _iter_beads(self, universe=None):
        u = self.atu if universe is None else universe

        for res in _trueprotein_residues(u):
            ca = res.atoms.select_atoms("calpha")

            cap_b = _cap_before(u, res)
            cap_a = _cap_after(u, res)

            if cap_b is not None:
                ca = ca + cap_b.atoms
            if cap_a is not None:
                ca = ca + cap_a.atoms

            if len(ca) > 0:
                yield "CA", ca.unique

        for ion in u.select_atoms("bioion").split("residue"):
            if ion.n_atoms > 0:
                yield "ions", ion

    def _add_bonds(self):
        bonds = []

        for s in self.segments:
            residues = list(s.residues)
            for r1, r2 in zip(residues[:-1], residues[1:]):
                if (r2.resid - r1.resid) != 1:
                    continue

                ca1 = r1.atoms.select_atoms("name CA")
                ca2 = r2.atoms.select_atoms("name CA")

                if ca1.n_atoms > 0 and ca2.n_atoms > 0:
                    bonds.append((ca1.ix[0], ca2.ix[0]))

        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)

    def _set_masses(self):
        ca_masses = []

        for res in _trueprotein_residues(self.atu):
            atoms = res.atoms

            cap_b = _cap_before(self.atu, res)
            cap_a = _cap_after(self.atu, res)

            if cap_b is not None:
                atoms = atoms + cap_b.atoms
            if cap_a is not None:
                atoms = atoms + cap_a.atoms

            ca_masses.append(atoms.total_mass())

        self.atoms.select_atoms("name CA").masses = np.asarray(ca_masses)

    def _set_charges(self):
        try:
            ca_charges = []

            for res in _trueprotein_residues(self.atu):
                atoms = res.atoms

                cap_b = _cap_before(self.atu, res)
                cap_a = _cap_after(self.atu, res)

                if cap_b is not None:
                    atoms = atoms + cap_b.atoms
                if cap_a is not None:
                    atoms = atoms + cap_a.atoms

                ca_charges.append(_safe_total_charge(atoms))

            self.atoms.select_atoms("name CA").charges = np.asarray(ca_charges)

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

    def _iter_beads(self, universe=None):
        u = self.atu if universe is None else universe

        for res in _trueprotein_residues(u):
            # ---- CA (with caps) ----
            ca = res.atoms.select_atoms("calpha")

            cap_b = _cap_before(u, res)
            cap_a = _cap_after(u, res)

            if cap_b is not None:
                ca = ca + cap_b.atoms
            if cap_a is not None:
                ca = ca + cap_a.atoms

            if len(ca) > 0:
                yield "CA", ca.unique

            # ---- CB (NO caps!) ----
            cb = res.atoms.select_atoms(self._mapping["CB"])
            if cb.n_atoms > 0:
                yield "CB", cb

        for ion in u.select_atoms("bioion").split("residue"):
            if ion.n_atoms > 0:
                yield "ions", ion

    def _apply_map(self, mapping):
        bead_items = list(self._iter_beads(self.atu))

        atomnames = []
        atomids = []
        resids = []
        resnames = []
        segids = []
        charges = []
        masses = []

        segid_map = {}

        for i, (name, bead) in enumerate(bead_items):
            atomnames.append(name)
            atomids.append(i)
            resids.append(bead.resids[0])
            resnames.append(bead.resnames[0])
            segids.append(_safe_segid(name, bead, segid_map=segid_map))
            charges.append(0.0)
            masses.append(0.0)

        n_atoms = len(bead_items)

        vdwradii = topologyattrs.Radii(np.zeros(n_atoms, dtype=float))
        atomids = topologyattrs.Atomids(np.asarray(atomids))
        atomnames = topologyattrs.Atomnames(np.asarray(atomnames, dtype=object))
        atomtypes = topologyattrs.Atomtypes(np.asarray(np.arange(n_atoms) + 100))
        charges = topologyattrs.Charges(np.asarray(charges, dtype=float))
        masses = topologyattrs.Masses(np.asarray(masses, dtype=float))

        segids = np.asarray(segids, dtype=object)
        resids = np.asarray(resids)
        resnames = np.asarray(resnames, dtype=object)

        # IMPORTANT: distinguish residues by both resid and segid
        residx, (new_resids, new_resnames, new_segids) = topbase.change_squash(
            (resids, segids), (resids, resnames, segids)
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

        # intra-residue CA-CB
        for r in self.residues:
            ca = r.atoms.select_atoms("name CA")
            cb = r.atoms.select_atoms("name CB")

            if ca.n_atoms > 0 and cb.n_atoms > 0:
                bonds.append((ca.ix[0], cb.ix[0]))

        # inter-residue CA(i)-CA(i+1), only for consecutive residue IDs
        for s in self.segments:
            residues = list(s.residues)
            for r1, r2 in zip(residues[:-1], residues[1:]):
                if (r2.resid - r1.resid) != 1:
                    continue

                ca1 = r1.atoms.select_atoms("name CA")
                ca2 = r2.atoms.select_atoms("name CA")

                if ca1.n_atoms > 0 and ca2.n_atoms > 0:
                    bonds.append((ca1.ix[0], ca2.ix[0]))

        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)
    
    def _set_masses(self):
        ca_masses = []
        cb_masses = []

        for r in _trueprotein_residues(self.atu):
            # ---- CA (with caps) ----
            atoms = r.atoms

            cap_b = _cap_before(self.atu, r)
            cap_a = _cap_after(self.atu, r)

            if cap_b is not None:
                atoms = atoms + cap_b.atoms
            if cap_a is not None:
                atoms = atoms + cap_a.atoms

            ca_masses.append(atoms.total_mass())

            # ---- CB (pure sidechain) ----
            cb = r.atoms.select_atoms(self._mapping["CB"])
            if cb.n_atoms > 0:
                cb_masses.append(cb.total_mass())

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

            for r in _trueprotein_residues(self.atu):
                # CA bead gets full residue + terminal cap(s)
                atoms = r.atoms

                cap_b = _cap_before(self.atu, r)
                cap_a = _cap_after(self.atu, r)

                if cap_b is not None:
                    atoms = atoms + cap_b.atoms
                if cap_a is not None:
                    atoms = atoms + cap_a.atoms

                ca_charges.append(_safe_total_charge(atoms))

                # CB bead stays pure sidechain
                cb = r.atoms.select_atoms(self._mapping["CB"])
                if cb.n_atoms > 0:
                    cb_charges.append(_safe_total_charge(cb))

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

    def _n_bead_atoms(self, res):
        atoms = res.atoms.select_atoms("amine") + res.atoms.select_atoms("hcalpha")

        cap = _cap_before(self.atu, res)
        if cap is not None:
            atoms = atoms + cap.atoms  # includes CH3

        return atoms.unique

    def _o_bead_atoms(self, res):
        atoms = res.atoms.select_atoms("carboxyl") + res.atoms.select_atoms("hcalpha")

        cap = _cap_after(self.atu, res)
        if cap is not None:
            atoms = atoms + cap.atoms  # includes CH3

        return atoms.unique

    def _iter_beads(self, universe=None):
        u = self.atu if universe is None else universe

        for res in _trueprotein_residues(u):
            n_atoms = self._n_bead_atoms(res)
            if n_atoms.n_atoms > 0:
                yield "N", n_atoms

            cb_sel = self._cb_selection(res.resname)
            if cb_sel is not None:
                cb = res.atoms.select_atoms(cb_sel)
                if cb.n_atoms > 0:
                    yield "CB", cb

            o_atoms = self._o_bead_atoms(res)
            if o_atoms.n_atoms > 0:
                yield "O", o_atoms

        for ion in u.select_atoms(self._mapping["ions"]).split("residue"):
            if ion.n_atoms > 0:
                yield "ion", ion

    def _apply_map(self, mapping):
        bead_items = list(self._iter_beads(self.atu))
        from collections import Counter
        atomnames = []
        atomids = []
        resids = []
        resnames = []
        segids = []
        charges = []
        masses = []

        segid_map = {}

        for i, (name, bead) in enumerate(bead_items):
            atomnames.append(name)
            atomids.append(i)
            resids.append(bead.resids[0])
            resnames.append(bead.resnames[0])
            segids.append(_safe_segid(name, bead, segid_map=segid_map))
            charges.append(0.0)
            masses.append(0.0)

        n_atoms = len(bead_items)

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

        # intra-residue bonds
        for r in self.residues:
            n = r.atoms.select_atoms("name N")
            cb = r.atoms.select_atoms("name CB")
            o = r.atoms.select_atoms("name O")

            if n.n_atoms > 0 and o.n_atoms > 0:
                bonds.append((n.ix[0], o.ix[0]))
            if n.n_atoms > 0 and cb.n_atoms > 0:
                bonds.append((n.ix[0], cb.ix[0]))
            if cb.n_atoms > 0 and o.n_atoms > 0:
                bonds.append((cb.ix[0], o.ix[0]))

        # inter-residue O(i)-N(i+1), only for consecutive residue IDs
        for s in self.segments:
            residues = list(s.residues)
            for r1, r2 in zip(residues[:-1], residues[1:]):
                if (r2.resid - r1.resid) != 1:
                    continue

                o = r1.atoms.select_atoms("name O")
                n = r2.atoms.select_atoms("name N")

                if o.n_atoms > 0 and n.n_atoms > 0:
                    bonds.append((o.ix[0], n.ix[0]))

        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)

    def _set_masses(self):
        n_masses, cb_masses, o_masses = [], [], []

        for r in _trueprotein_residues(self.atu):
            n_atoms = self._n_bead_atoms(r)
            o_atoms = self._o_bead_atoms(r)

            cb_sel = self._cb_selection(r.resname)
            cb = r.atoms.select_atoms(cb_sel) if cb_sel is not None else r.atoms[:0]

            if n_atoms.n_atoms > 0:
                n_masses.append(n_atoms.total_mass())

            if cb.n_atoms > 0:
                cb_masses.append(cb.total_mass())

            if o_atoms.n_atoms > 0:
                o_masses.append(o_atoms.total_mass())

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
            n_charges, cb_charges, o_charges = [], [], []

            for r in _trueprotein_residues(self.atu):
                n_atoms = self._n_bead_atoms(r)
                o_atoms = self._o_bead_atoms(r)

                cb_sel = self._cb_selection(r.resname)
                cb = r.atoms.select_atoms(cb_sel) if cb_sel is not None else r.atoms[:0]

                if n_atoms.n_atoms > 0:
                    n_charges.append(_safe_total_charge(n_atoms))

                if cb.n_atoms > 0:
                    cb_charges.append(_safe_total_charge(cb))

                if o_atoms.n_atoms > 0:
                    o_charges.append(_safe_total_charge(o_atoms))

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