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
    with mda.Writer(filenames["topology_file"], **kwargs) as rtf:
        logger.info("Writing {}...".format(filenames["topology_file"]))
        rtf.write(universe)
    with mda.Writer(filenames["stream_file"], **kwargs) as stream:
        logger.info("Writing {}...".format(filenames["stream_file"]))
        stream.write(universe)
    with mda.Writer(filenames["psf_file"], **kwargs) as psf:
        logger.info("Writing {}...".format(filenames["psf_file"]))
        psf.write(universe, xplor=False, atom_type_source="types")

    # Write the new trajectory.
    if write_traj:
        universe.trajectory.rewind()
        with mda.Writer(
            filenames["traj_file"],
            universe.atoms.n_atoms,
            istart=1,
            remarks="Written by fluctmatch.",
        ) as trj:
            logger.info("Writing the trajectory {}...".format(
                filenames["traj_file"]))
            logger.warning("This may take a while depending upon the size and "
                        "length of the trajectory.")
            with click.progressbar(universe.trajectory) as bar:
                for _ in bar:
                    trj.write(universe.atoms)

    #universe._generate_from_topology()
    with mda.Writer(filenames["xplor_psf_file"], **kwargs) as psf:
        logger.info("Writing {}...".format(filenames["xplor_psf_file"]))
        psf.write(universe, xplor=True, atom_type_source="names")

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
            filenames["crd_file"], dt=1.0, **kwargs) as crd:
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
    precision = kwargs.get("precision", 5)

    if index is not None:
        command = [
            gromacs_exec,
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
            "-ndec",
            "{:d}".format(precision),
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
            "-ndec",
            "{:d}".format(precision),
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

def split_mda(info, data_dir=path.join(os.getcwd(), "data"), frame_based=True, **kwargs):
    """Create a subtrajectory from a trajectory using MDAnalysis.

    Parameters
    ----------
    info : tuple
        Tuple of (subdir, start, stop), where subdir is the output folder
        label and start/stop define the selected interval.
    data_dir : str, optional
        Location of the main data directory.
    frame_based : bool, optional
        If True, select using 1-based frame indices:
        start <= ts.frame + 1 <= stop

        If False, select using trajectory time:
        start <= ts.time <= stop

    Notes
    -----
    Frame-based selection is safer for concatenated trajectories because
    trajectory time may reset between runs.
    """

    subdir, start, stop = info
    subdir = path.join(data_dir, "{}".format(subdir))

    try:
        os.makedirs(subdir)
    except OSError:
        pass

    topology = kwargs.get("topology", "md.tpr")
    trajectory = kwargs.get("trajectory", path.join(os.curdir, "md.xtc"))
    outfile = path.join(subdir, kwargs.get("outfile", "aa.xtc"))
    logfile = path.join(subdir, kwargs.get("logfile", "split.log"))
    precision = kwargs.get("precision", 5)

    u = mda.Universe(topology, trajectory)
    ag = u.atoms
    n_written = 0

    logger.info(
        "split_mda starting | outfile=%s | start=%s | stop=%s | frame_based=%s",
        outfile, start, stop, frame_based
    )

    with mda.Writer(outfile, ag.n_atoms, precision=precision) as W:
        if frame_based:
            for ts in u.trajectory[start - 1:stop]:
                W.write(ag)
                n_written += 1
        else:
            for ts in u.trajectory:
                current_value = ts.time
                if current_value < start:
                    continue
                if current_value > stop:
                    break
                W.write(ag)
                n_written += 1

    with mdutil.openany(logfile, mode="w") as log:
        print("Frames written: {}".format(n_written), file=log)
        print("Start: {}".format(start), file=log)
        print("Stop: {}".format(stop), file=log)
        print("Frame based: {}".format(frame_based), file=log)

    logger.info(
        "split_mda completed | outfile=%s | frames_written=%d",
        outfile, n_written
    )
