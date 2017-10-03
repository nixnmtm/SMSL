# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import os
from os import path

import MDAnalysis as mda
import MDAnalysis.analysis.base as analysis
import numpy as np
import pandas as pd
from MDAnalysis.coordinates import memory
from future.utils import native_str


def average_structure(universe):
    """Calculate the average structure of a trajectory.

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
        A collection of atoms in a universe or AtomGroup.

    Returns
    -------
    :class:`~numpy.ndarray`
        Average coordinates of the trajectory.
    """
    system = analysis.AnalysisFromFunction(lambda atoms: atoms.positions, universe.trajectory, universe.atoms).run()
    positions = np.mean(system.results, axis=0)
    return positions


def bond_stats(universe, func="mean"):
    """Calculate the average bond lengths from a trajectory.

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
        A collection of atoms in a universe or AtomGroup with bond information.
    func : {"mean", "std"}
        Determine either the mean or the standard deviation of the bond lengths.

    Returns
    -------
    :class:`~pandas.DataFrame`
        Average bond lengths of the trajectory.
    """
    system = analysis.AnalysisFromFunction(lambda bonds: bonds.bonds(), universe.trajectory, universe.bonds).run()
    if func == "mean":
        callback = np.mean
    elif func == "std":
        callback = np.std
    else:
        raise KeyError("The 'func' keyword must either be 'mean' or 'std'."
                       "")
    bonds = pd.concat([
        pd.Series(universe.bonds.atom1.names),
        pd.Series(universe.bonds.atom2.names),
        pd.Series(callback(system.results, axis=0)),
    ], axis=1)
    bonds.columns = ["I", "J", "r_IJ"]
    return bonds.set_index(["I", "J"])


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
        crd_file=".".join((filename, "crd")),
        stream_file=".".join((filename, "stream")),
        topology_file=".".join((filename, "rtf")),
        traj_file=".".join((filename, "xtc")),
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

    # Calculate the average coordinates, average bond lengths, and
    # fluctuations of bond lengths from the trajectory.
    print("Determining the average structure of the trajectory. ")
    print("Note: This could take a while depending upon the side of your trajectory.")
    positions = average_structure(universe)
    avg_universe = mda.Universe(
        filenames["psf_file"],
        [positions, ],
        format=memory.MemoryReader,
        order="fac"
    )
    print("Writing {}...".format(filenames["crd_file"]))
    with mda.Writer(native_str(filenames["crd_file"]), dt=1.0, **kwargs) as crd:
        crd.write(avg_universe.atoms)

    # Write the new trajectory in Gromacs XTC format.
    if write_traj:
        print("Writing the trajectory {}.".format(filenames["traj_file"]))
        print("This may take a while depending upon the size and length of the trajectory.")
        with mda.Writer(native_str(filenames["traj_file"]), n_atoms=universe.atoms.n_atoms, dt=1.0, **kwargs) as trj:
            for ts in universe.trajectory:
                trj.write(ts)

    # Write an XPLOR version of the PSF
    atomtypes = topologyattrs.Atomtypes(universe.atoms.names)
    universe._topology.add_TopologyAttr(topologyattr=atomtypes)
    universe._generate_from_topology()
    print("Writing {}...".format(filenames["xplor_psf_file"]))
    with mda.Writer(native_str(filenames["xplor_psf_file"]), **kwargs) as psf:
        psf.write(universe)
