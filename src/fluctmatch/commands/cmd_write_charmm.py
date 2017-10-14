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
from fluctmatch.fluctmatch import utils as fmutils


@click.command("write_charmm", short_help="Write various CHARMM files.")
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
    "outdir",
    metavar="DIR",
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default="cg",
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-c",
    "--charmm",
    "charmm_version",
    metavar="VERSION",
    default=41,
    type=click.INT,
    help="CHARMM version",
)
@click.option(
    "--extended / --standard",
    "extended",
    default=True,
    help="Output using the extended or standard columns",
)
@click.option(
    "--nb / --no-nb",
    "nonbonded",
    default=True,
    help="Include nonbonded section in CHARMM parameter file",
)
@click.option(
    "--cmap / --no-cmap",
    "cmap",
    default=True,
    help="Include CMAP section in CHARMM PSF file",
)
@click.option(
    "--cheq / --no-cheq",
    default=True,
    help="Include charge equilibrium section in CHARMM PSF file",
)
@click.option(
    "--write / --no-write",
    "write_traj",
    default=True,
    help="Convert the trajectory file",
)
def cli(
    topology, trajectory, outdir, prefix, charmm_version,
    extended, cmap, cheq, nonbonded, write_traj,
):
    kwargs = dict(
        prefix=prefix,
        outdir=outdir,
        charmm_version=charmm_version,
        extended=extended,
        cmap=cmap,
        cheq=cheq,
        nonbonded=nonbonded,
        write_traj=write_traj,
    )
    universe = mda.Universe(topology, trajectory)
    fmutils.write_charmm_files(universe, **kwargs)
