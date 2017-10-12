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

import click
from MDAnalysis.lib.util import which
from fluctmatch.cli import pass_context
from fluctmatch.fluctmatch import charmmfluctmatch


@click.command("run-fm", short_help="Run fluctuation matching.")
@click.option(
    "-s",
    "topology",
    default=path.join(os.getcwd(), "md.tpr"),
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
    default=path.join(os.getcwd(), "md.xtc"),
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
    default=os.getcwd(),
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-e",
    "--exec",
    "nma_exec",
    envvar="CHARMMEXEC",
    default=which("charmm"),
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="CHARMM executable file",
)
@click.option(
    "-t",
    "--temperature",
    type=click.FLOAT,
    default=300.0,
    help="Temperature of simulation",
)
@click.option(
    "--rmin",
    type=click.FLOAT,
    default=0.0,
    help="Minimum distance between bonds",
)
@click.option(
    "--rmax",
    type=click.FLOAT,
    default=10.0,
    help="Maximum distance between bonds",
)
@click.option(
    "-n",
    "--ncycles",
    "n_cycles",
    type=click.INT,
    default=250,
    help="Number of simulation cycles",
)
@click.option(
    "--tol",
    type=click.FLOAT,
    default=1.e-4,
    help="Tolerance level between simulations",
)
@click.option(
    "-p",
    "--prefix",
    default="fluctmatch",
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-c",
    "--charmm",
    "charmm_version",
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
    "--resid / --no-resid",
    "resid",
    default=True,
    help="Include segment IDs in internal coordinate files",
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
    "--restart",
    is_flag=True,
    help="Restart simulation",
)
@pass_context
def cli(
    topology, trajectory, outdir, nma_exec, temperature,
    n_cylces, tol, rmin, rmax, prefix, charmm_version,
    extended, resid, cmap, cheq, nonbonded, restart,
):
    kwargs = dict(
        prefix=prefix,
        outdir=outdir,
        temperature=temperature,
        rmin=rmin,
        rmax=rmax,
        charmm_version=charmm_version,
        extended=extended,
        resid=resid,
        cmap=cmap,
        cheq=cheq,
        nonbonded=nonbonded,
    )
    cfm = charmmfluctmatch.CharmmFluctMatch(topology, trajectory, **kwargs)
    cfm.initialize(restart=restart)
    cfm.run(nma_exec=nma_exec, tol=tol, n_cycles=n_cylces)
