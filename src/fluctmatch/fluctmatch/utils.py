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

import logging
import os
import subprocess
import tempfile
import textwrap
from os import path

import click
import math
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
        self._nframes = atomgroup.universe.trajectory.n_frames

    def _prepare(self):
        self.result = np.zeros_like(self._ag.positions)

    def _single_frame(self):
        self.result += self._ag.positions

    def _conclude(self):
        self.result /= self._nframes


class BondAverage(analysis.AnalysisBase):
    """Calculate the average bond length.

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
        self._nframes = atomgroup.universe.trajectory.n_frames

    def _prepare(self):
        self.result = np.zeros_like(self._ag.bonds.bonds())

    def _single_frame(self):
        self.result += self._ag.bonds.bonds()

    def _conclude(self):
        self.result = np.rec.fromarrays(
            [
                self._ag.bonds.atom1.names,
                self._ag.bonds.atom2.names,
                self.result / self._nframes
            ],
            names=["I", "J", "r_IJ"]
        )
        self.result = pd.DataFrame.from_records(self.result)


class BondStd(analysis.AnalysisBase):
    """Calculate the fluctuation in bond lengths.

    """

    def __init__(self, atomgroup, average, **kwargs):
        """
        Parameters
        ----------
        atomgroup : :class:`~MDAnalysis.Universe.AtomGroup`
            An AtomGroup
        average : float or ""lass:`~numpy.array`
            Average bond length
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
        self._nframes = atomgroup.universe.trajectory.n_frames
        self._average = average

    def _prepare(self):
        self.result = np.zeros_like(self._ag.bonds.bonds())

    def _single_frame(self):
        self.result += np.square(self._ag.bonds.bonds() - self._average)

    def _conclude(self):
        self.result = np.rec.fromarrays(
            [
                self._ag.bonds.atom1.names,
                self._ag.bonds.atom2.names,
                np.sqrt(self.result / self._nframes)
            ],
            names=["I", "J", "r_IJ"]
        )
        self.result = pd.DataFrame.from_records(self.result)


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
        universe.trajectory.rewind()
        with mda.Writer(
            native_str(filenames["traj_file"]),
            universe.atoms.n_atoms,
            istart=1.0,
                remarks="Written by fluctmatch.") as trj:
            logger.info("Writing the trajectory {}...".format(
                filenames["traj_file"]))
            logger.warning("This may take a while depending upon the size and "
                           "length of the trajectory.")
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
    logger.info("Determining the average structure of the trajectory. ")
    logger.warning("Note: This could take a while depending upon the "
                   "size of your trajectory.")
    positions = AverageStructure(universe.atoms).run().result
    positions = positions.reshape((*positions.shape, 1))

    # Create a new universe.
    topologies = ("names", "resids", "resnums", "resnames", "segids")
    avg_universe = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=universe.residues.n_residues,
        n_segments=universe.segments.n_segments,
        atom_resindex=universe.atoms.resindices,
        residue_segindex=universe.residues.segindices,
        trajectory=True)
    for _ in topologies:
        avg_universe.add_TopologyAttr(_)
    avg_universe.atoms.names = universe.atoms.names
    avg_universe.residues.resids = universe.residues.resids
    avg_universe.residues.resnums = universe.residues.resnums
    avg_universe.residues.resnames = universe.residues.resnames
    avg_universe.segments.segids = universe.segments.segids
    avg_universe.load_new(positions, order="acf")

    # avg_universe.load_new(
    #     positions, )
    with mda.Writer(
            native_str(filenames["crd_file"]), dt=1.0, **kwargs) as crd:
        logger.info("Writing {}...".format(filenames["crd_file"]))
        crd.write(avg_universe.atoms)


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
    outfile = path.join(subdir, kwargs.get("outfile", "aa.xtc"))
    logfile = path.join(subdir, kwargs.get("logfile", "split.log"))

    if index is not None:
        command = [
            "gmx",
            "trjconv",
            "-s",
            topology,
            "-f",
            trajectory,
            "-n",
            index,
            "-o",
            outfile,
            "-b",
            "{:d}".format(start),
            "-e",
            "{:d}".format(stop),
            "-ndec", 5
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
            outfile,
            "-b",
            "{:d}".format(start),
            "-e",
            "{:d}".format(stop),
            "-ndec", 5
        ]
    fd, fpath = tempfile.mkstemp(text=True)
    with mdutil.openany(fpath, "w") as temp:
        print(kwargs.get("system", 0), file=temp)
    with mdutil.openany(fpath, "r") as temp, \
        mdutil.openany(logfile, mode="w") as log:
        logger.info("Writing trajectory to {}".format(outfile))
        logger.info("Writing Gromacs output to {}".format(logfile))
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
    logfile = path.join(subdir, kwargs.get("logfile", "split.log"))
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
