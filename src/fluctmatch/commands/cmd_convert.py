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
from future.utils import (viewkeys, iteritems)
from fluctmatch import (_DESCRIBE, _MODELS)
from fluctmatch.models.core import modeller
from fluctmatch.fluctmatch.utils import write_charmm_files


@click.command("convert", short_help="Convert from all-atom to coarse-grain model.")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.tpr"),
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=True
    ),
    help="Gromacs topology file (e.g., tpr gro g96 pdb brk ent)",
)
@click.option(
    "-f",
    "trajectory",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.xtc"),
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=True
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-o",
    "--outdir",
    metavar="DIR",
    default=os.getcwd(),
    type=click.Path(
        exists=False,
        file_okay=False,
        resolve_path=True
    ),
    help="Directory",
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
    "--rmin",
    metavar="DIST",
    type=click.FLOAT,
    default=0.0,
    help="Minimum distance between bonds",
)
@click.option(
    "--rmax",
    metavar="DIST",
    type=click.FLOAT,
    default=10.0,
    help="Maximum distance between bonds",
)
@click.option(
    "-m",
    "--model",
    metavar="MODEL",
    type=click.Choice(viewkeys(_MODELS)),
    multiple=True,
    help="Model(s) to convert to",
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
@click.option(
    "-l",
    "--list",
    "model_list",
    is_flag=True,
    help="List available models with their descriptions"
)
def cli(
    topology, trajectory, outdir, prefix, rmin, rmax, model, charmm_version,
    com, extended, resid, cmap, cheq, nonbonded, write_traj, model_list,
):
    if model_list:
        for k, v in iteritems(_DESCRIBE):
            print("{:20}{}".format(k, v))
        return

    kwargs = dict(
        outdir=outdir,
        prefix=prefix,
        rmin=rmin,
        rmax=rmax,
        charmm_version=charmm_version,
        extended=extended,
        resid=resid,
        cmap=cmap,
        cheq=cheq,
        nonbonded=nonbonded,
        write_traj=write_traj,
    )
    universe = modeller(topology, trajectory, com=com, model=model, **kwargs)
    write_charmm_files(universe, **kwargs)
