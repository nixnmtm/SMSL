# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import MDAnalysis as mda
import MDAnalysis.analysis.base as analysis
import numpy as np


def average_structure(*args, **kwargs):
    """Calculate the average structure from a trajectory.

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
    :return: coordinates
    """
    universe = mda.Universe(*args, **kwargs)
    system = analysis.AnalysisFromFunction(lambda atoms: atoms.positions, universe.trajectory, universe.atoms).run()
    positions = np.mean(system.results, axis=0)
    return positions


def average_bonds(*args, **kwargs):
    """Calculate the average bond lengths from a trajectory.

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
    :return: bond lengths
    """
    universe = mda.Universe(*args, **kwargs)
    system = analysis.AnalysisFromFunction(lambda bonds: bonds.bonds(), universe.trajectory, universe.bonds).run()
    bonds = np.mean(system.results, axis=0)
    return bonds

def bond_fluctuation(*args, **kwargs):
    """Calculate the average bond lengths from a trajectory.

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
    :return: bond lengths
    """
    universe = mda.Universe(*args, **kwargs)
    system = analysis.AnalysisFromFunction(lambda bonds: bonds.bonds(), universe.trajectory, universe.bonds).run()
    bonds = np.std(system.results, axis=0)
    return bonds
