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
# The original code is from Richard J. Gowers.
# https://github.com/richardjgowers/MDAnalysis-coarsegraining
#
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)
from future.builtins import super, zip
from future.utils import (
    raise_with_traceback,
    viewitems,
    with_metaclass,
)

import abc
import itertools
import logging
import string

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import (
    topology,
    topologyattrs,
)
from MDAnalysis.lib.util import asiterable
from MDAnalysis.topology import base as topbase
from MDAnalysis.topology import guessers
from fluctmatch import (_DESCRIBE, _MODELS)
from fluctmatch.models import (
    topattrs,
    trajectory,
)

logger = logging.getLogger(__name__)


class _ModelMeta(abc.ABCMeta):
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = asiterable(classdict['model'])
        except KeyError:
            pass
        else:
            for f in fmt:
                f = f.upper()
                _MODELS[f] = cls
        try:
            description = asiterable(classdict['describe'])
        except KeyError:
            pass
        else:
            for f, d in zip(fmt, description):
                f = f.upper()
                _DESCRIBE[f] = d


class ModelBase(with_metaclass(_ModelMeta, mda.Universe)):
    """Base class for creating coarse-grain models.

    Parameters
    ----------
    topology : filename or Topology object
        A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
        define the list of atoms. If the file includes bond information,
        partial charges, atom masses, ... then these data will be available to
        MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
        topology) is always required. Alternatively, an existing
        :class:`MDAnalysis.core.topology.Topology` instance may also be given.
    extended
        Renames the residues and atoms according to the extended CHARMM PSF
        format. Standard CHARMM PSF limits the residue and atom names to four
        characters, but the extended CHARMM PSF permits eight characters. The
        residues and atoms are renamed according to the number of segments
        (1: A, 2: B, etc.) and then the residue number or atom index number.
     xplor
        Assigns the atom type as either a numerical or an alphanumerical
        designation. CHARMM normally assigns a numerical designation, but the
        XPLOR version permits an alphanumerical designation with a maximum
        size of 4. The numerical form corresponds to the atom index number plus
        a factor of 100, and the alphanumerical form will be similar the
        standard CHARMM atom name.
    topology_format
        Provide the file format of the topology file; ``None`` guesses it from
        the file extension [``None``] Can also pass a subclass of
        :class:`MDAnalysis.topology.base.TopologyReaderBase` to define a custom
        reader to be used on the topology file.
    format
        Provide the file format of the coordinate or trajectory file; ``None``
        guesses it from the file extension. Note that this keyword has no
        effect if a list of file names is supplied because the "chained" reader
        has to guess the file format for each individual list member.
        [``None``] Can also pass a subclass of
        :class:`MDAnalysis.coordinates.base.ProtoReader` to define a custom
        reader to be used on the trajectory file.
    guess_bonds : bool, optional
        Once Universe has been loaded, attempt to guess the connectivity
        between atoms.  This will populate the .bonds .angles and .dihedrals
        attributes of the Universe.
    vdwradii : dict, optional
        For use with *guess_bonds*. Supply a dict giving a vdwradii for each
        atom type which are used in guessing bonds.
    is_anchor : bool, optional
        When unpickling instances of
        :class:`MDAnalysis.core.groups.AtomGroup` existing Universes are
        searched for one where to anchor those atoms. Set to ``False`` to
        prevent this Universe from being considered. [``True``]
    anchor_name : str, optional
        Setting to other than ``None`` will cause
        :class:`MDAnalysis.core.groups.AtomGroup` instances pickled from the
        Universe to only unpickle if a compatible Universe with matching
        *anchor_name* is found. Even if *anchor_name* is set *is_anchor* will
        still be honored when unpickling.
    in_memory
        After reading in the trajectory, transfer it to an in-memory
        representations, which allow for manipulation of coordinates.
    in_memory_step
        Only read every nth frame into in-memory representation.

    Attributes
    ----------
    trajectory
        currently loaded trajectory reader;
    dimensions
        current system dimensions (simulation unit cell, if set in the
        trajectory)
    atoms, residues, segments
        master Groups for each topology level
    bonds, angles, dihedrals
        master ConnectivityGroups for each connectivity type
    """

    def __init__(self, *args, **kwargs):
        """Initialise like a normal MDAnalysis Universe but give the mapping and
        com keywords.

        Mapping must be a dictionary with atom names as keys.
        Each name must then correspond to a selection string,
        signifying how to split up a single residue into many beads.
        eg:
        mapping = {"CA":"protein and name CA",
                   "CB":"protein and not name N HN H HT* H1 H2 H3 CA HA* C O OXT
                   OT*"}
        would split residues into 2 beads containing the C-alpha atom and the
        sidechain.
        """
        # Coarse grained Universe
        # Make a blank Universe for myself.
        super().__init__()

        self._com = kwargs.pop("com", True)

        # Atomistic Universe
        try:
            self.atu = mda.Universe(*args, **kwargs)
        except (IOError, OSError, ValueError):
            raise_with_traceback(RuntimeError("Failed to create a universe."))

    def __repr__(self):
        message = "<CG Universe with {} beads".format(len(self.atoms))
        try:
            message += " and {:d} bonds".format(
                len(self._topology.bonds.values))
        except AttributeError as exc:
            pass
        finally:
            message += ">"
        return message

    def _initialize(self, *args, **kwargs):
        try:
            mapping = kwargs.pop("mapping")
        except KeyError:
            raise ValueError("CG mapping has not been defined.")

        # Fake up some beads
        self._topology = self._apply_map(mapping)
        self._generate_from_topology()
        self._add_bonds()
        if kwargs.get("guess_angles", True):
            self._add_angles()
            self._add_dihedrals()
            self._add_impropers()

        # This replaces load_new in a traditional Universe
        try:
            self.trajectory = trajectory._Trajectory(
                self.atu, mapping, n_atoms=self.atoms.n_atoms, com=self._com)
        except (IOError, TypeError) as exc:
            raise_with_traceback(
                RuntimeError("Unable to open {}".format(
                    self.atu.trajectory.filename)))

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
            itertools.product(residues, viewitems(mapping)))
        for i, (res, (name, selection)) in select_residues:
            bead = res.select_atoms(selection)
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
        # _beads = topattrs._Beads(_beads)
        vdwradii = np.zeros_like(atomids)
        vdwradii = topologyattrs.Radii(vdwradii)
        atomids = topologyattrs.Atomids(np.asarray(atomids))
        atomnames = topologyattrs.Atomnames(
            np.asarray(atomnames, dtype=np.object))
        atomtypes = topologyattrs.Atomtypes(
            np.asarray(np.arange(n_atoms) + 100))
        charges = topologyattrs.Charges(np.asarray(charges))
        masses = topologyattrs.Masses(np.asarray(masses))

        # Residue
        # resids, resnames
        segids = np.asarray(segids, dtype=np.object)
        resids = np.asarray(resids, dtype=np.int32)
        resnames = np.asarray(resnames, dtype=np.object)
        residx, (new_resids, new_resnames,
                 perres_segids) = topbase.change_squash(
                     (resids, resnames, segids), (resids, resnames, segids))

        # transform from atom:Rid to atom:Rix
        residueids = topologyattrs.Resids(new_resids)
        residuenums = topologyattrs.Resnums(new_resids.copy())
        residuenames = topologyattrs.Resnames(new_resnames)

        # Segment
        segidx, perseg_segids = topbase.squash_by(perres_segids)[:2]
        segids = topologyattrs.Segids(perseg_segids)

        # Setup topology
        top = topology.Topology(
            len(atomids),
            len(new_resids),
            len(segids),
            attrs=[
                atomids, atomnames, atomtypes, charges, masses,
                vdwradii, residueids, residuenums, residuenames, segids
            ],
            atom_resindex=residx,
            residue_segindex=segidx)
        return top

    @abc.abstractmethod
    def _add_bonds(self):
        pass

    def _add_angles(self):
        try:
            angles = guessers.guess_angles(self.bonds)
            self._topology.add_TopologyAttr(topologyattrs.Angles(angles))
            self._generate_from_topology()
        except AttributeError:
            pass

    def _add_dihedrals(self):
        try:
            dihedrals = guessers.guess_dihedrals(self.angles)
            self._topology.add_TopologyAttr(topologyattrs.Dihedrals(dihedrals))
            self._generate_from_topology()
        except AttributeError:
            pass

    def _add_impropers(self):
        try:
            impropers = guessers.guess_improper_dihedrals(self.angles)
            self._topology.add_TopologyAttr(
                (topologyattrs.Impropers(impropers)))
            self._generate_from_topology()
        except AttributeError:
            pass

    def _set_masses(self):
        pass

    def _set_charges(self):
        pass

    @property
    def cguniverse(self):
        """Convert a :class:`~MDAnalysis.AtomGroup` to a :class:`~MDAnalysis.Universe`.

        Returns
        -------
        :class:`~MDAnalysis.Universe`
        """
        return Merge(self)


def Merge(*args):
    """Combine multiple coarse-grain systems into one.

    Parameters
    ----------
    args : iterable of either :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`

    Returns
    -------
    :class:`~MDAnalysis.Universe`
        A merged universe.
    """
    import multiprocessing as mp
    from MDAnalysis.coordinates.memory import MemoryReader

    logger.warning("This might take a while depending upon the number of "
                   "trajectory frames.")
    if not all([
            u.universe.trajectory.n_frames == args[0]
            .universe.trajectory.n_frames for u in args
    ]):
        logger.error("The trajectories are not the same length.")
        raise ValueError("The trajectories are not the same length.")
    ag = [_.atoms for _ in args]
    unit_cell = np.mean([_._unitcell for _ in args[0].trajectory], axis=0)
    universe = mda.Merge(*ag)

    if args[0].universe.trajectory.n_frames > 1:
        traj = (_.trajectory for _ in args)
        coordinates = [
            np.concatenate([_.positions for _ in ts], axis=0)
            for ts in zip(*traj)
        ]
        coordinates = np.array(coordinates)
        if universe.atoms.n_atoms != coordinates.shape[1]:
            logger.error(
                "The number of sites does not match the number of coordinates."
            )
            raise RuntimeError(
                "The number of sites does not match the number of coordinates."
            )
        logger.info("The new universe has {1} beads in {0} frames.".format(
            *coordinates.shape))

        universe.load_new(coordinates, format=MemoryReader)
        logger.warning(
            "The new trajectory will is assigned an average unit cell "
            "for the entire trajectory. This is currently a limitation "
            "implemented by MDAnalysis.")
        universe.trajectory.ts._unitcell = unit_cell
    return universe


def rename_universe(universe):
    """Rename the atoms and residues within a universe.

    Standardizes naming of the universe by renaming atoms and residues based
    upon the number of segments. Atoms are labeled as 'A001', 'A002', 'A003',
    ..., 'A999' for the first segment, and 'B001', 'B002', 'B003', ..., 'B999'
    for the second segment. Residues are named in a similar fashion according to
    their segment.

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe`
        A collection of atoms in a universe.

    Returns
    -------
    :class:`~MDAnalysis.Universe`
        The universe with renamed residues and atoms.
    """
    logger.info("Renaming atom names and atom types within the universe.")
    atomnames = np.array([
        "{}{:0>3d}".format(lett, i)
        for lett, segment in zip(string.ascii_uppercase, universe.segments)
        for i, _ in enumerate(segment.atoms, 1)
    ])
    resnames = np.array([
        "{}{:0>3d}".format(lett, i)
        for lett, segment in zip(string.ascii_uppercase, universe.segments)
        for i, _ in enumerate(segment.residues, 1)
    ])

    universe._topology.add_TopologyAttr(topologyattrs.Atomnames(atomnames))
    universe._topology.add_TopologyAttr(topologyattrs.Resnames(resnames))
    if not np.issubdtype(universe.atoms.types.dtype, np.int):
        universe._topology.add_TopologyAttr(topologyattrs.Atomtypes(atomnames))
    universe._generate_from_topology()
    return universe
