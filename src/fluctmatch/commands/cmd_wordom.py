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

import click
import MDAnalysis as mda
from fluctmatch.cli import pass_context
from fluctmatch.fluctmatch import utils as fmutils


@click.command("write-charmm", short_help="Write various CHARMM files.")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Gromacs topology file (e.g., tpr gro g96 pdb brk ent)",
)
@click.option(
    "-f",
    "trajectory",
    metavar="FILE",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-o",
    "outfile",
    metavar="FILE",
    type=click.Path(
        exists=False,
        file_okay=False,
        resolve_path=True
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-b",
    "start",
    metavar="FRAME",
    default=0,
    type=click.INT,
    help="Start time of trajectory",
)
@click.option(
    "-e",
    "stop",
    metavar="FRAME",
    default=-1,
    type=click.INT,
    help="Stop time of total trajectory",
)
@click.option(
    "-dt",
    "step",
    metavar="STEP",
    default=1,
    type=click.INT,
    help="Interval between frames"
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Show progress of output"
)
def cli(topology, trajectory, outfile, start, stop, step, verbose):
    universe = mda.Universe(topology, trajectory)
    trajectory = universe.trajectory
    with mda.Writer(outfile, n_atoms=universe.atoms.n_atoms, dt=trajectory.dt) as trj:
        if verbose:
            with click.progressbar(
                length=trajectory.n_frames,
                label="Trajectory frames written"
            ) as bar:
                    for ts in trajectory[start:stop:step]:
                        trj.write(ts)
                        bar.update(ts.frame)
        else:
            for ts in trajectory[start:stop:step]:
                trj.write(ts)
