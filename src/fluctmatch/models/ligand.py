# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

from collections import OrderedDict

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topology, topologyattrs
from MDAnalysis.topology import base as topbase

from fluctmatch.models.base import ModelBase


def _safe_total_charge(ag):
    try:
        return ag.total_charge()
    except AttributeError:
        return 0.0

def _is_psf_safe_segid(segid):
        if segid is None:
            return False
        segid = str(segid).strip()
        if not segid:
            return False
        if " " in segid:
            return False
        if len(segid) > 8:
            return False
        return True

def _safe_segid(name, bead, default="SYS"):
    """
    Preserve an existing SEGID if it is already PSF-safe.
    Otherwise, synthesize one based on bead/model identity.
    """
    raw = None
    try:
        if bead.segids is not None and len(bead.segids) > 0:
            raw = str(bead.segids[0]).strip()
    except Exception:
        raw = None

    # keep an existing segid only if it looks PSF-safe
    if _is_psf_safe_segid(raw) and raw.upper() not in {
        "SYSTEM", 
        "PROTEININWATER", 
        "PROTEIN IN WATER"
        }:
        return raw

    # fallback: assign a clean segid by bead/model type
    if name in {"CA", "CB", "N", "O"}:
        return "PROT"
    if name in {"ion", "ions"}:
        return "ION"

    # for ligands, use residue name if available
    try:
        resname = str(bead.resnames[0]).strip()
        if resname:
            return resname[:8]
    except Exception:
        pass

    return default

class LigandBase(ModelBase):
    """Base class for coarse-grained ligand models."""
    resname = None
    _mapping = OrderedDict()
    _bonded = []

    def _iter_beads(self, universe=None):
        u = self.atu if universe is None else universe

        for res in u.residues:
            if self.resname is not None and res.resname != self.resname:
                continue

            for bead_name, sel in self._mapping.items():
                bead = res.atoms.select_atoms(sel)
                if bead.n_atoms > 0:
                    yield bead_name, bead

    def _apply_map(self, mapping):
        bead_items = list(self._iter_beads(self.atu))

        atomnames = []
        atomids = []
        resids = []
        resnames = []
        segids = []
        charges = []
        masses = []

        for i, (name, bead) in enumerate(bead_items):
            atomnames.append(name)
            atomids.append(i)
            resids.append(bead.resids[0])
            resnames.append(bead.resnames[0])
            segids.append(_safe_segid(name, bead))
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
        bead_index = {a.name: a.ix for a in self.atoms}

        for a_name, b_name in self._bonded:
            a_ix = bead_index.get(a_name, None)
            b_ix = bead_index.get(b_name, None)
            if a_ix is not None and b_ix is not None:
                bonds.append((a_ix, b_ix))

        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)

    def _set_masses(self):
        masses = []
        for _, bead in self._iter_beads(self.atu):
            masses.append(bead.total_mass())
        self.atoms.masses = np.asarray(masses)

    def _set_charges(self):
        try:
            charges = []
            for _, bead in self._iter_beads(self.atu):
                charges.append(bead.total_charge())
            self.atoms.charges = np.asarray(charges)
        except AttributeError:
            pass


class ADP(LigandBase):
    """Coarse-grained ADP model using nucleic6-like base decomposition."""
    model = "ADP"
    describe = "ADP with explicit base/sugar/phosphate pharmacophore-like beads"
    resname = "ADP"
    _mapping = OrderedDict([
        # adenine split
        ("A1", "name C5 N7 C8 H8 N9"),
        ("A2", "name N1 C2 H2 N3 C4 C6 N6 H61 H62"),
        # sugar split
        ("S1", "name C4' H4' O4' C1' H1'"),
        ("S2", "name C2' H2'' O2' H2' C3' H3' O3' H3T"),
        # phosphate side
        ("P1", "name C5' H5' H5'' O5' PA O1A O2A O3A"),
        ("P2", "name PB O1B O2B O3B"),
    ])
    _bonded = [
        ("A1", "A2"),
        ("A1", "S1"),
        ("S1", "S2"),
        ("S1", "P1"),
        ("P1", "P2"),
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_masses()
        self._set_charges()


class OM(LigandBase):
    """Coarse-grained Omecamtiv Mecarbil model."""
    model = "OM"
    describe = "16-bead pharmacophore-oriented OM model"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._mapping["OM1"] = "resname 2OW and name C18 H19 H20 H21"
        self._mapping["OM2"] = "resname 2OW and name C3 C4 C5 C19 C17 H2 H3 H22"
        self._mapping["OM3"] = "resname 2OW and name N5 H2"
        self._mapping["OM4"] = "resname 2OW and name N4 H18"
        self._mapping["OM5"] = "resname 2OW and name O3"
        self._mapping["OM6"] = "resname 2OW and name N3 H17"
        self._mapping["OM7"] = "resname 2OW and name C15 C2 C14"
        self._mapping["OM8"] = "resname 2OW and name C13 C1 C6"
        self._mapping["OM9"] = "resname 2OW and name F"
        self._mapping["OM10"] = "resname 2OW and name C12"
        self._mapping["OM11"] = "resname 2OW and name N2 C7 C11"
        self._mapping["OM12"] = "resname 2OW and name N1 C20 C10"
        self._mapping["OM13"] = "resname 2OW and name O2"
        self._mapping["OM14"] = "resname 2OW and name O1"
        self._mapping["OM15"] = "resname 2OW and name C8"

        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        self._set_masses()
        self._set_charges()

    def _add_bonds(self):
        bonds = []

        bead_names = [f"OM{i}" for i in range(1, 17)]
        bead_ix = {}
        for name in bead_names:
            ag = self.atoms.select_atoms(f"name {name}")
            if ag.n_atoms > 0:
                bead_ix[name] = ag.ix[0]

        bonded_pairs = [
            ("OM1", "OM2"),
            ("OM2", "OM3"),
            ("OM2", "OM4"),
            ("OM4", "OM5"),
            ("OM5", "OM6"),
            ("OM6", "OM7"),
            ("OM7", "OM8"),
            ("OM8", "OM9"),
            ("OM8", "OM10"),
            ("OM10", "OM11"),
            ("OM11", "OM12"),
            ("OM12", "OM13"),
            ("OM13", "OM14"),
            ("OM14", "OM15"),
            ("OM15", "OM16"),
        ]

        for a, b in bonded_pairs:
            if a in bead_ix and b in bead_ix:
                bonds.append((bead_ix[a], bead_ix[b]))

        if bonds:
            self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        mda.Universe.__init__(self, self._topology)