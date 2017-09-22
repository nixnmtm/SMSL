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
import os

from future.builtins import (
    dict,
)
from future.utils import (
    with_metaclass,
)


class FluctMatch(with_metaclass(abc.ABCMeta, object)):
    """Base class for fluctuation matching."""
    parameters = dict()
    target = dict()

    def __init__(self, *args, **kwargs):
        """Initialization of fluctuation matching.

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
        outdir
            Output directory
        temperature
            Temperature (in K)
        rmin
            Minimum distance to consider for bond lengths.
        rmax
            Maximum distance to consider for bond lengths.
        """
        self.outdir = kwargs.get("outdir", os.curdir)
        self.temperature = kwargs.get("temperature", 300.0)
        self.args = args
        self.kwargs = kwargs

    @abc.abstractmethod
    def initialize(self, restart=False):
        """Create an elastic network model from a basic coarse-grain model.

        :param restart: reinitialize tye object
        """
        pass

    @abc.abstractmethod
    def run(self, nma_exec=None, tol=1.e-4, n_cycles=250):
        """Perform a self-consistent fluctuation matching

        :param nma_exec: executable file for normal mode analysis
        :param tol: error tolerance
        :param n_cycles: number of fluctuation matching cycles
        """
        pass

    @abc.abstractmethod
    def restart(self, nma_exec=None, tol=1.e-4, n_cycles=250):
        """Restart a simulation.

        :param nma_exec: executable file for normal mode analysis
        :param tol: error tolerance
        :param n_cycles: number of fluctuation matching cycles
        """
        pass
