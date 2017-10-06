# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)
from future.utils import PY2

import os
from os import path

import MDAnalysis as mda
import MDAnalysis.analysis.base as analysis
import numpy as np
import pandas as pd
from MDAnalysis.coordinates import memory
from future.builtins import super
from future.utils import native_str

if PY2:
    FileNotFoundError = OSError


class AverageStructure(analysis.AnalysisBase):
    """Calculate the average structure of a trajectory.
    """
    def __init__(self, atomgroup, **kwargs):
        """
        Parameters
        ----------
        atomgroup : :class:`~MDAnalysis.Universe.AtomGroup`
            An AtomGroup
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        verbose : bool, optional
            Turn on verbosity
        """
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self._ag = atomgroup

    def _prepare(self):
        self.result = []

    def _single_frame(self):
        self.result.append(self._ag.positions)

    def _conclude(self):
        self.result = np.mean(self.result, axis=0)


class BondStats(analysis.AnalysisBase):
    """Calculate either the average bond length or the fluctuation in bond lengths.

    """
    def __init__(self, atomgroup, func="mean", **kwargs):
        """
        Parameters
        ----------
        atomgroup : :class:`~MDAnalysis.Universe.AtomGroup`
            An AtomGroup
        func : {"mean", "std"}, optional
            Calculate either the mean or the standard deviation of the bonds
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        verbose : bool, optional
            Turn on verbosity
        """
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self._ag = atomgroup
        if func == "mean":
            self._func = np.mean
        elif func == "std":
            self._func = np.std
        else:
            raise AttributeError("func must either be 'mean' or 'std'")

    def _prepare(self):
        self.result = []

    def _single_frame(self):
        self.result.append(self._ag.bonds.bonds())

    def _conclude(self):
        self.result = self._func(self.result, axis=0)
        bonds = pd.concat([
            pd.Series(self._ag.bonds.atom1.names),
            pd.Series(self._ag.bonds.atom2.names),
            pd.Series(self.result),
        ], axis=1)
        bonds.columns = ["I", "J", "r_IJ"]
        self.result = bonds.copy(deep=True)


def write_charmm_files(universe, outdir=os.curdir, prefix="cg", write_traj=True, **kwargs):
    """Write CHARMM coordinate, topology PSF, stream, and topology RTF files.

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
        A collection of atoms in a universe or AtomGroup with bond definitions.
    outdir : str
        Location to write the files.
    prefix : str
        Prefix of filenames
    write_traj : bool
        Write the trajectory to disk.
    charmm_version
        Version of CHARMM for formatting (default: 41)
    extended
        Use the extended format.
    cmap
        Include CMAP section.
    cheq
        Include charge equilibration.
    title
        Title lines at the beginning of the file.
    """
    from MDAnalysis.core import (topologyattrs,)

    filename = path.join(outdir, prefix)
    filenames = dict(
        psf_file=".".join((filename, "psf")),
        xplor_psf_file=".".join((filename, "xplor", "psf")),
        crd_file=".".join((filename, "cor")),
        stream_file=".".join((filename, "stream")),
        topology_file=".".join((filename, "rtf")),
        traj_file=".".join((filename, "dcd")),
    )

    # Write required CHARMM input files.
    print("Writing {}...".format(filenames["topology_file"]))
    with mda.Writer(native_str(filenames["topology_file"]), **kwargs) as rtf:
        rtf.write(universe)
    print("Writing {}...".format(filenames["stream_file"]))
    with mda.Writer(native_str(filenames["stream_file"]), **kwargs) as stream:
        stream.write(universe)
    print("Writing {}...".format(filenames["psf_file"]))
    with mda.Writer(
        native_str(filenames["psf_file"]), **kwargs) as psf:
        psf.write(universe)

    # Write the new trajectory in Gromacs XTC format.
    if write_traj:
        kwargs["start"] = 1
        kwargs["step"] = 1
        print("Writing the trajectory {}...".format(filenames["traj_file"]))
        print("This may take a while depending upon the size and length of the trajectory.")
        with mda.Writer(native_str(filenames["traj_file"]), universe.atoms.n_atoms, dt=1.0, **kwargs) as trj:
            for ts in universe.trajectory:
                trj.write(ts)

    # Write an XPLOR version of the PSF
    atomtypes = topologyattrs.Atomtypes(universe.atoms.names)
    universe._topology.add_TopologyAttr(topologyattr=atomtypes)
    universe._generate_from_topology()
    print("Writing {}...".format(filenames["xplor_psf_file"]))
    with mda.Writer(native_str(filenames["xplor_psf_file"]), **kwargs) as psf:
        psf.write(universe)

    # Calculate the average coordinates, average bond lengths, and
    # fluctuations of bond lengths from the trajectory.
    if universe.trajectory.n_frames > 1:
        print("Determining the average structure of the trajectory. ")
        print("Note: This could take a while depending upon the side of your trajectory.")
        try:
            u = mda.Universe(filenames["psf_file"], filenames["traj_file"])
        except FileNotFoundError:
            u = mda.Universe(universe.filename, universe.trajectory.filename)
        positions = AverageStructure(u.atoms).run().result
        avg_universe = mda.Universe(
            filenames["psf_file"],
            [positions, ],
            format=memory.MemoryReader,
            order="fac"
        )
        print("Writing {}...".format(filenames["crd_file"]))
        with mda.Writer(native_str(filenames["crd_file"]), dt=1.0, **kwargs) as crd:
            crd.write(avg_universe.atoms)
    else:
        print("Writing {}...".format(filenames["crd_file"]))
        with mda.Writer(native_str(filenames["crd_file"]), dt=1.0, **kwargs) as crd:
            crd.write(universe.atoms)
