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

import os
from os import path

import click
from MDAnalysis.lib.util import which
from fluctmatch.fluctmatch import charmmfluctmatch


@click.command("run_fm", short_help="Run fluctuation matching.")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.tpr"),
    show_default=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Gromacs topology file (e.g., tpr gro g96 pdb brk ent)",
)
@click.option(
    "-f",
    "trajectory",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.xtc"),
    show_default=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-o",
    "outdir",
    metavar="DIR",
    default=os.getcwd(),
    show_default=True,
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Directory",
)
@click.option(
    "-e",
    "--exec",
    "nma_exec",
    metavar="FILE",
    envvar="CHARMMEXEC",
    default=which("charmm"),
    show_default=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="CHARMM executable file",
)
@click.option(
    "-t",
    "--temperature",
    metavar="TEMP",
    type=click.FLOAT,
    default=300.0,
    show_default=True,
    help="Temperature of simulation",
)
@click.option(
    "-n",
    "--ncycles",
    "n_cycles",
    metavar="NCYCLES",
    type=click.IntRange(1, None, clamp=True),
    default=250,
    show_default=True,
    help="Number of simulation cycles",
)
@click.option(
    "--tol",
    metavar="TOL",
    type=click.FLOAT,
    default=1.e-4,
    show_default=True,
    help="Tolerance level between simulations",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default="fluctmatch",
    show_default=True,
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-c",
    "--charmm",
    "charmm_version",
    metavar="VERSION",
    default=41,
    show_default=True,
    type=click.IntRange(27, None, clamp=True),
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
    "--resid / --no-resid",
    "resid",
    default=True,
    help="Include segment IDs in internal coordinate files",
)
@click.option(
    "--restart",
    is_flag=True,
    help="Restart simulation",
)
def cli(
        topology,
        trajectory,
        outdir,
        nma_exec,
        temperature,
        n_cycles,
        tol,
        prefix,
        charmm_version,
        extended,
        resid,
        nonbonded,
        restart,
):
    kwargs = dict(
        prefix=prefix,
        outdir=outdir,
        temperature=temperature,
        charmm_version=charmm_version,
        extended=extended,
        resid=resid,
        nonbonded=nonbonded,
    )
    cfm = charmmfluctmatch.CharmmFluctMatch(topology, trajectory, **kwargs)
    cfm.initialize(restart=restart)
    cfm.run(nma_exec=nma_exec, tol=tol, n_cycles=n_cycles)
