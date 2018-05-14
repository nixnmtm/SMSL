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
from future.builtins import (dict, super)
from future.utils import (PY2, native_str)

import copy
import logging
import os
import subprocess
import tempfile
import textwrap
from os import path

import click
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.analysis.base as analysis
from MDAnalysis.coordinates import memory
from MDAnalysis.lib import util as mdutil
from fluctmatch.fluctmatch.data import charmm_split

if PY2:
    FileNotFoundError = OSError

logger = logging.getLogger(__name__)


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
        func : {"mean", "std", "both"}, optional
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
            self._func = (np.mean, )
        elif func == "std":
            self._func = (np.std, )
        elif func == "both":
            self._func = (np.mean, np.std)
        else:
            raise AttributeError("func must either be 'mean' or 'std'")

    def _prepare(self):
        self.result = []

    def _single_frame(self):
        self.result.append(self._ag.bonds.bonds())

    def _conclude(self):
        self.result = [func(self.result, axis=0) for func in self._func]
        bonds = [
            pd.concat(
                [
                    pd.Series(self._ag.bonds.atom1.names),
                    pd.Series(self._ag.bonds.atom2.names),
                    pd.Series(_),
                ],
                axis=1) for _ in self.result
        ]
        for _ in bonds:
            _.columns = ["I", "J", "r_IJ"]
        self.result = copy.deepcopy(bonds)


def write_charmm_files(universe,
                       outdir=os.getcwd(),
                       prefix="cg",
                       write_traj=True,
                       **kwargs):
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
    from MDAnalysis.core import (
        topologyattrs, )

    # Attempt to create the necessary subdirectory
    try:
        os.makedirs(outdir)
    except OSError:
        pass

    filename = path.join(outdir, prefix)
    filenames = dict(
        psf_file=".".join((filename, "psf")),
        xplor_psf_file=".".join((filename, "xplor", "psf")),
        crd_file=".".join((filename, "cor")),
        stream_file=".".join((filename, "stream")),
        topology_file=".".join((filename, "rtf")),
        traj_file=".".join((filename, "dcd")),
    )

    n_atoms = universe.atoms.n_atoms
    n_bonds = len(universe.bonds)
    n_angles = len(universe.angles)
    n_dihedrals = len(universe.dihedrals)
    n_impropers = len(universe.impropers)
    logger.warning("The system has {:d} atoms, {:d} bonds, {:d} angles, {:d} "
                   "dihedrals, and {:d} impropers. Depending upon "
                   "the size of the system, file writing may take a while and "
                   "have a large file size.".format(n_atoms, n_bonds, n_angles,
                                                    n_dihedrals, n_impropers))

    # Write required CHARMM input files.
    with mda.Writer(native_str(filenames["topology_file"]), **kwargs) as rtf:
        logger.info("Writing {}...".format(filenames["topology_file"]))
        rtf.write(universe)
    with mda.Writer(native_str(filenames["stream_file"]), **kwargs) as stream:
        logger.info("Writing {}...".format(filenames["stream_file"]))
        stream.write(universe)
    with mda.Writer(native_str(filenames["psf_file"]), **kwargs) as psf:
        logger.info("Writing {}...".format(filenames["psf_file"]))
        psf.write(universe)

    # Write the new trajectory in Gromacs XTC format.
    if write_traj:
        with mda.Writer(
                native_str(filenames["traj_file"]),
                universe.atoms.n_atoms,
                remarks="Written by fluctmatch.") as trj:
            logger.info("Writing the trajectory {}...".format(
                filenames["traj_file"]))
            logger.warning("This may take a while depending upon the size and "
                           "length of the trajectory.")
            universe.trajectory.rewind()
            with click.progressbar(universe.trajectory) as bar:
                for ts in bar:
                    trj.write(ts)

    # Write an XPLOR version of the PSF
    atomtypes = topologyattrs.Atomtypes(universe.atoms.names)
    universe._topology.add_TopologyAttr(topologyattr=atomtypes)
    universe._generate_from_topology()
    with mda.Writer(native_str(filenames["xplor_psf_file"]), **kwargs) as psf:
        logger.info("Writing {}...".format(filenames["xplor_psf_file"]))
        psf.write(universe)

    # Calculate the average coordinates from the trajectory.
    if universe.trajectory.n_frames > 1:
        logger.info("Determining the average structure of the trajectory. ")
        logger.warning("Note: This could take a while depending upon the "
                       "size of your trajectory.")
        positions = AverageStructure(universe.atoms).run().result
        positions = positions.reshape((*positions.shape, 1))
        avg_universe = universe.copy()
        avg_universe.load_new(
            positions, format=memory.MemoryReader, order="acf")
        with mda.Writer(
                native_str(filenames["crd_file"]), dt=1.0, **kwargs) as crd:
            logger.info("Writing {}...".format(filenames["crd_file"]))
            crd.write(avg_universe.atoms)
    else:
        with mda.Writer(
                native_str(filenames["crd_file"]), dt=1.0, **kwargs) as crd:
            logger.info("Writing {}...".format(filenames["crd_file"]))
            crd.write(universe.atoms)


def split_gmx(info, data_dir=path.join(os.getcwd(), "data"), **kwargs):
    """Create a subtrajectory from a Gromacs trajectory.

    Parameters
    ----------
    info : :class:`collections.namedTuple`
        Contains information about the data subdirectory and start and
        stop frames
    data_dir : str, optional
        Location of the main data directory
    topology : str, optional
        Topology filename (e.g., tpr gro g96 pdb brk ent)
    trajectory : str, optional
        A Gromacs trajectory file (e.g., xtc trr)
    index : str, optional
        A Gromacs index file (e.g., ndx)
    outfile : str, optional
        A Gromacs trajectory file (e.g., xtc trr)
    logfile : str, optional
        Log file for output of command
    system : int
        Atom selection from Gromacs index file (0 = System, 1 = Protein)
    """
    # Trajectory splitting information
    subdir, start, stop = info
    subdir = path.join(data_dir, "{}".format(subdir))
    gromacs_exec = mdutil.which("gmx")

    # Attempt to create the necessary subdirectory
    try:
        os.makedirs(subdir)
    except OSError:
        pass

    # Various filenames
    topology = kwargs.get("topology", "md.tpr")
    trajectory = kwargs.get("trajectory", path.join(os.curdir, "md.xtc"))
    index = kwargs.get("index")
    outfile = kwargs.get("outfile", "aa.xtc")
    logfile = kwargs.get("logfile", "split.log")

    if index is not None:
        command = [
            "gmx",
            "-s",
            topology,
            "-f",
            trajectory,
            "-n",
            index,
            "-o",
            path.join(subdir, outfile),
            "-b",
            "{:d}".format(start),
            "-e",
            "{:d}".format(stop),
        ]
    else:
        command = [
            gromacs_exec,
            "trjconv",
            "-s",
            topology,
            "-f",
            trajectory,
            "-o",
            path.join(subdir, outfile),
            "-b",
            "{:d}".format(start),
            "-e",
            "{:d}".format(stop),
        ]
    fd, fpath = tempfile.mkstemp(text=True)
    with mdutil.openany(fpath, "w") as temp:
        print(kwargs.get("system", 0), file=temp)
    with mdutil.openany(fpath, "r") as temp, \
        mdutil.openany(path.join(subdir, logfile), mode="w") as log:
        subprocess.check_call(
            command, stdin=temp, stdout=log, stderr=subprocess.STDOUT)
    os.remove(fpath)


def split_charmm(info, data_dir=path.join(os.getcwd(), "data"), **kwargs):
    """Create a subtrajectory from a CHARMM trajectory.

    Parameters
    ----------
    info : :class:`collections.namedTuple`
        Contains information about the data subdirectory and start and
        stop frames
    data_dir : str, optional
        Location of the main data directory
    toppar : str, optional
        Directory containing CHARMM topology/parameter files
    trajectory : str, optional
        A CHARMM trajectory file (e.g., dcd)
    outfile : str, optional
        A CHARMM trajectory file (e.g., dcd)
    logfile : str, optional
        Log file for output of command
    charmm_version : int
        Version of CHARMM
    """
    # Trajectory splitting information
    subdir, start, stop = info
    subdir = path.join(data_dir, "{}".format(subdir))
    charmm_exec = mdutil.which("charmm")

    # Attempt to create the necessary subdirectory
    try:
        os.makedirs(subdir)
    except OSError:
        pass

    # Various filenames
    version = kwargs.get("charmm_version", 41)
    toppar = kwargs.get("toppar",
                        "/opt/local/charmm/c{:d}b1/toppar".format(version))
    trajectory = kwargs.get("trajectory", path.join(os.curdir, "md.dcd"))
    outfile = path.join(subdir, kwargs.get("outfile", "aa.dcd"))
    logfile = kwargs.get("logfile", "split.log")
    inpfile = path.join(subdir, "split.inp")

    with mdutil.openany(inpfile, "w") as charmm_input:
        charmm_inp = charmm_split.split_inp.format(
            toppar=toppar,
            trajectory=trajectory,
            outfile=outfile,
            version=version,
            start=start,
            stop=stop,
        )
        charmm_inp = textwrap.dedent(charmm_inp[1:])
        print(charmm_inp, file=charmm_input)
    command = [
        charmm_exec,
        "-i",
        inpfile,
        "-o",
        path.join(subdir, logfile),
    ]
    subprocess.check_call(command)
