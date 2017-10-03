# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import abc
import itertools
import string

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core import (
    topology,
    topologyattrs,
)
from MDAnalysis.lib.util import asiterable
from MDAnalysis.topology import base as topbase
from MDAnalysis.topology import guessers
from future.builtins import (
    super,
    zip,
)
from future.utils import (
    raise_with_traceback,
    viewitems,
    with_metaclass,
)

from . import (
    topattrs,
    trajectory,
)
from .. import _MODELS


class _ModelMeta(type):
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
        Renames the residues and atoms according to the extended CHARMM PSF format.
        Standard CHARMM PSF limits the residue and atom names to four characters,
        but the extended CHARMM PSF permits eight characters. The residues and
        atoms are renamed according to the number of segments (1: A, 2: B, etc.)
        and then the residue number or atom index number.
     xplor
        Assigns the atom type as either a numerical or an alphanumerical
        designation. CHARMM normally assigns a numerical designation, but the
        XPLOR version permits an alphanumerical designation with a maximum
        size of 4. The numerical form corresponds to the atom index number plus a
        factor of 100, and the alphanumerical form will be similar the standard
        CHARMM atom name.
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
        """Initialise like a normal MDAnalysis Universe but give the mapping and com keywords.

        Mapping must be a dictionary with atom names as keys.
        Each name must then correspond to a selection string,
        signifying how to split up a single residue into many beads.
        eg:
        mapping = {"CA":"protein and name CA",
                   "CB":"protein and not name N HN H HT* H1 H2 H3 CA HA* C O OXT OT*"}
        Would split residues into 2 beads containing the C-alpha atom and the sidechain.
        """
        # Coarse grained Universe
        # Make a blank Universe for myself.
        super().__init__()

        self._com = kwargs.pop("com", True)
        self._extended = kwargs.pop("extended", True)
        self._xplor = kwargs.pop("xplor", True)

        # Atomistic Universe
        try:
            self.atu = mda.Universe(*args, **kwargs)
        except (IOError, OSError, ValueError):
            raise_with_traceback(RuntimeError("Failed to create a universe."))

    def __repr__(self):
        message = "<CG Universe with {} beads".format(len(self.atoms._beads))
        try:
            message += " and {:d} bonds".format(len(self._topology.bonds.values))
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
            self.trajectory = trajectory._Trajectory(self.atu, mapping, n_atoms=self.atoms.n_atoms, com=self._com)
        except (IOError, TypeError) as exc:
            raise_with_traceback(RuntimeError("Unable to open {}".format(self.atu.trajectory.filename)))

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
        for i, (res, (name, selection)) in enumerate(itertools.product(residues, viewitems(mapping))):
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
        _beads = topattrs._Beads(_beads)
        vdwradii = np.zeros_like(atomids)
        vdwradii = topologyattrs.Radii(vdwradii)
        atomids = topologyattrs.Atomids(np.asarray(atomids))
        atomnames = topologyattrs.Atomnames(np.asarray(atomnames))
        atomtypes = topologyattrs.Atomtypes(np.asarray(np.arange(n_atoms)+100))
        charges = topologyattrs.Charges(np.asarray(charges))
        masses = topologyattrs.Masses(np.asarray(masses))

        # Residue
        # resids, resnames
        segids = np.asarray(segids)
        resids = np.asarray(resids)
        resnames = np.asarray(resnames)
        residx, (new_resids, new_resnames, new_segids) = topbase.change_squash((resids,), (resids, resnames, segids))

        # transform from atom:Rid to atom:Rix
        residueids = topologyattrs.Resids(new_resids)
        residuenums = topologyattrs.Resnums(new_resids.copy())
        residuenames = topologyattrs.Resnames(new_resnames)

        # Segment
        segidx, (perseg_segids,) = topbase.change_squash((new_segids,), (new_segids,))
        segids = topologyattrs.Segids(perseg_segids)

        # Setup topology
        top = topology.Topology(len(atomids), len(new_resids), len(segids),
                                attrs=[_beads, atomids, atomnames, atomtypes,
                                       charges, masses, vdwradii, residueids,
                                       residuenums, residuenames, segids],
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
            self._topology.add_TopologyAttr((topologyattrs.Impropers(impropers)))
            self._generate_from_topology()
        except AttributeError:
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
    args : list or tuple of either :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`

    Returns
    -------
    :class:`~MDAnalysis.Universe`
        A merged universe.
    """
    from MDAnalysis.coordinates.memory import MemoryReader
    from MDAnalysis.analysis.base import AnalysisFromFunction

    print("This might take a while depending upon the number of trajectory frames.")
    if not np.all([issubclass(u.__class__, mda.Universe) for u in args]):
        raise TypeError("The universes must all be derived from MDAnalysis.Universe.")
    if not np.all([u.trajectory.n_frames == args[0].trajectory.n_frames for u in args]):
        raise ValueError("The trajectories are not the same length.")
    ag = [_.atoms for _ in args]
    universe = mda.Merge(*ag)

    atoms = universe.atoms
    attrs = universe._topology
    residx, (new_resids, new_resnames, new_segids) = topbase.change_squash((atoms.resids,),
                                                                           (atoms.resids, atoms.resnames, atoms.segids))
    # Update segment list
    # transform from atom:Rid to atom:Rix
    residueids = topologyattrs.Resids(new_resids)
    residuenums = topologyattrs.Resnums(new_resids.copy())
    residuenames = topologyattrs.Resnames(new_resnames)

    # Segment
    segidx, (perseg_segids, ) = topbase.change_squash((new_segids,), (new_segids,))
    segids = topologyattrs.Segids(perseg_segids)
    universe._topology = topology.Topology(attrs.n_atoms, attrs.n_residues, len(segids),
                                           attrs=[attrs.indices, attrs.ids, attrs.names, attrs.types,
                                                  attrs.charges, attrs.masses, attrs.radii, residueids,
                                                  residuenums, residuenames, segids],
                                           atom_resindex=residx,
                                           residue_segindex=segidx)
    universe._generate_from_topology()

    if args[0].trajectory.n_frames > 1:
        coordinates = [
            AnalysisFromFunction(lambda u: u.positions.copy(), u).run().results
            for u in ag
        ]
        coordinates = np.concatenate(coordinates, axis=1)
        if universe.atoms.n_atoms != coordinates.shape[1]:
            raise RuntimeError("The number of sites does not match the number of coordinates.")
        print("The new universe has {1} beads in {0} frames.".format(*coordinates.shape))
        universe.load_new(coordinates, format=MemoryReader)
    return universe


def rename_universe(universe):
    """Rename the atoms and residues within a universe.

    Standardizes naming of the universe by renaming atoms and residues based upon the
    number of segments. Atoms are labeled as 'A001', 'A002', 'A003', ..., 'A999' for
    the first segment, and 'B001', 'B002', 'B003', ..., 'B999' for the second segment.
    Residues are named in a similar fashion according to their segment.

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe`
        A collection of atoms in a universe.

    Returns
    -------
    :class:`~MDAnalysis.Universe`
        The universe with renamed residues and atoms.
    """
    atomnames = np.array(["{}{:0>3d}".format(lett, i)
                          for lett, segment in zip(string.ascii_uppercase, universe.segments)
                          for i, _ in enumerate(segment.atoms, 1)])
    resnames = np.array(["{}{:0>3d}".format(lett, i)
                         for lett, segment in zip(string.ascii_uppercase, universe.segments)
                         for i, _ in enumerate(segment.residues, 1)])

    universe._topology.add_TopologyAttr(topologyattrs.Atomnames(atomnames))
    universe._topology.add_TopologyAttr(topologyattrs.Resnames(resnames))
    if not np.issubdtype(universe.atoms.types.dtype, np.int):
        universe._topology.add_TopologyAttr(topologyattrs.Atomtypes(atomnames))
    universe._generate_from_topology()
    return universe
