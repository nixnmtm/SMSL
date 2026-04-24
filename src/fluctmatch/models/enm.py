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

from collections import defaultdict

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib.distances import distance_array

from fluctmatch.fluctmatch import utils as fmutils
from fluctmatch.models.base import (
    ModelBase,
    rename_universe,
)


class _EnmBase(ModelBase):
    """Shared ENM construction logic."""

    model = "ENM_BASE"
    describe = "Elastic network model base"

    def __init__(self, *args, **kwargs):
        self._rmin = kwargs.pop("rmin", 0.0)
        self._rmax = kwargs.pop("rmax", 10.0)
        super().__init__(*args, **kwargs)
        self._initialize(*args, **kwargs)

    def __repr__(self):
        message = "<CG Universe with {} beads".format(self.atoms.n_atoms)
        try:
            message += " and {:d} bonds".format(len(self._topology.bonds.values))
        except AttributeError:
            pass
        finally:
            message += ">"
        return message

    def _attach_coords_from_atu(self):
        coords = self.atu.atoms.positions.copy()
        self.load_new(coords[np.newaxis, :, :], format=MemoryReader)

    def _initialize(self, *args, **kwargs):
        self.__dict__.update(self.atu.__dict__)

        rename_universe(self)

        charges = kwargs.get("charges", False)
        if not charges:
            self._topology.add_TopologyAttr(
                topologyattrs.Charges(np.zeros(self.atoms.n_atoms))
            )

        self._topology.add_TopologyAttr(
            topologyattrs.Atomtypes(np.arange(self.atoms.n_atoms) + 1)
        )
        self._topology.add_TopologyAttr(topologyattrs.Angles([]))
        self._topology.add_TopologyAttr(topologyattrs.Dihedrals([]))
        self._topology.add_TopologyAttr(topologyattrs.Impropers([]))

        # Rebuild topology view once after topology attrs are updated
        mda.Universe.__init__(self, self._topology)

        self._add_bonds()

        # FINAL step: attach coordinates after all topology rebuilds
        self._attach_coords_from_atu()

        if kwargs.get("guess_angles", False):
            self._add_angles()
            self._add_dihedrals()
            self._add_impropers()

    def _average_positions(self):
        return fmutils.AverageStructure(self.atu.atoms).run().result

    def _candidate_pairs_from_distance(self, positions):
        distmat = distance_array(positions, positions, backend="OpenMP")

        if self._rmin > 0.0:
            a0, a1 = np.where((distmat >= self._rmin) & (distmat <= self._rmax))
        else:
            a0, a1 = np.where((distmat > self._rmin) & (distmat <= self._rmax))

        pairs = []
        for i, j in zip(a0, a1):
            if j > i:
                pairs.append((int(i), int(j), float(distmat[i, j])))
        return pairs

    def _skeleton_bonds(self):
        try:
            return set(
                tuple(sorted((int(bond[0].ix), int(bond[1].ix))))
                for bond in self.atu.bonds
            )
        except Exception:
            return set()

    def _finalize_bonds(self, bond_set):
        self._topology.add_TopologyAttr(topologyattrs.Bonds(sorted(bond_set)))
        mda.Universe.__init__(self, self._topology)

    def _select_enm_bonds(self, positions):
        """Override in subclasses."""
        raise NotImplementedError

    def _add_bonds(self):
        positions = self._average_positions()
        bond_set = self._select_enm_bonds(positions)
        self._finalize_bonds(bond_set)

    @property
    def rmin(self):
        return self._rmin

    @rmin.setter
    def rmin(self, distance):
        self._rmin = distance

    @property
    def rmax(self):
        return self._rmax

    @rmax.setter
    def rmax(self, distance):
        self._rmax = distance


class Enm(_EnmBase):
    model = "ENM"
    describe = "Elastic network model"

    def _select_enm_bonds(self, positions):
        skeleton_bonds = self._skeleton_bonds()
        candidate_pairs = self._candidate_pairs_from_distance(positions)

        rmax_bonds = set((i, j) for i, j, _ in candidate_pairs)
        return skeleton_bonds.union(rmax_bonds)


class EnmSparse(_EnmBase):
    model = "ENM_SPARSE"
    describe = "Elastic network model with reduced redundant long-range contacts"

    def __init__(self, *args, **kwargs):
        self._local_seq_sep = kwargs.pop("local_seq_sep", 2)
        self._short_contact_cutoff = kwargs.pop("short_contact_cutoff", 4.5)
        self._max_long_range_per_node = kwargs.pop("max_long_range_per_node", 10)
        super().__init__(*args, **kwargs)

    def _node_metadata(self):
        atoms = self.atu.atoms

        segids = []
        resids = []

        for atom in atoms:
            try:
                segid = str(atom.segid).strip()
            except Exception:
                segid = ""
            segids.append(segid)

            try:
                resid = int(atom.resid)
            except Exception:
                resid = int(atom.ix) + 1
            resids.append(resid)

        return np.array(segids, dtype=object), np.array(resids, dtype=int)

    def _select_enm_bonds(self, positions):
        skeleton_bonds = self._skeleton_bonds()
        candidate_pairs = self._candidate_pairs_from_distance(positions)
        segids, resids = self._node_metadata()

        local_bonds = set()
        short_bonds = set()
        long_range_candidates = []

        for i, j, d in candidate_pairs:
            pair = (i, j)

            same_seg = (segids[i] == segids[j])
            seq_sep = abs(resids[i] - resids[j])

            # Keep all local sequence-neighbor contacts on same segment
            if same_seg and seq_sep <= self._local_seq_sep:
                local_bonds.add(pair)
                continue

            # Keep all very short contacts regardless of sequence separation
            if d <= self._short_contact_cutoff:
                short_bonds.add(pair)
                continue

            # Everything else is considered pruneable long-range contact
            long_range_candidates.append((i, j, d))

        # Cap long-range contacts per node by shortest distance
        per_node = defaultdict(list)
        for i, j, d in long_range_candidates:
            per_node[i].append((d, j))
            per_node[j].append((d, i))

        kept_long_range = set()
        for i, neigh in per_node.items():
            neigh = sorted(neigh, key=lambda x: x[0])
            for d, j in neigh[: self._max_long_range_per_node]:
                kept_long_range.add(tuple(sorted((i, j))))

        return skeleton_bonds.union(local_bonds).union(short_bonds).union(kept_long_range)

    @property
    def local_seq_sep(self):
        return self._local_seq_sep

    @property
    def short_contact_cutoff(self):
        return self._short_contact_cutoff

    @property
    def max_long_range_per_node(self):
        return self._max_long_range_per_node