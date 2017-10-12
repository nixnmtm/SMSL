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
from future.utils import viewkeys
from fluctmatch import _MODELS
from fluctmatch.models.core import modeller
from fluctmatch.cli import pass_context
from fluctmatch.fluctmatch.utils import write_charmm_files


@click.command("convert", short_help="Convert from all-atom to coarse-grain model.")
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
    "--outdir",
    default=os.getcwd(),
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Directory",
)
@click.option(
    "-p",
    "--prefix",
    default="cg",
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-m",
    "--models",
    type=click.Choice(viewkeys(_MODELS)),
    multiple=True,
    help="Model(s) to convert to",
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
    "--com / --cog",
    "com",
    default=True,
    help="Use either center of mass or center of geometry",
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
    "--write / --no-write",
    "write_traj",
    default=True,
    help="Convert the trajectory file",
)
@pass_context
def cli(
    topology, trajectory, outdir, prefix, models, charmm_version,
    com, extended, resid, cmap, cheq, nonbonded, write_traj,
):
    kwargs = dict(
        prefix=prefix,
        outdir=outdir,
        charmm_version=charmm_version,
        extended=extended,
        resid=resid,
        cmap=cmap,
        cheq=cheq,
        nonbonded=nonbonded,
        write_traj=write_traj,
    )
    universe = modeller(topology, trajectory, com=com, dt=1.0, model=models)
    write_charmm_files(universe, **kwargs)
