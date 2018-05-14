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
from future.utils import (viewkeys, iteritems)

import logging
import logging.config
import os
from os import path

import click
from fluctmatch import (_DESCRIBE, _MODELS)
from fluctmatch.models.core import modeller
from fluctmatch.fluctmatch.utils import write_charmm_files


@click.command(
    "convert", short_help="Convert from all-atom to coarse-grain model.")
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
    "--outdir",
    metavar="DIR",
    show_default=True,
    default=os.getcwd(),
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Directory",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default="cg",
    show_default=True,
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "--rmin",
    metavar="DIST",
    type=click.FLOAT,
    default=0.0,
    show_default=True,
    help="Minimum distance between bonds",
)
@click.option(
    "--rmax",
    metavar="DIST",
    type=click.FLOAT,
    default=10.0,
    show_default=True,
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
    show_default=True,
    type=click.IntRange(27, None, clamp=True),
    help="CHARMM version",
)
@click.option(
    "--com / --cog",
    "com",
    default=True,
    show_default=True,
    help="Use either center of mass or center of geometry",
)
@click.option(
    "--extended / --standard",
    "extended",
    default=True,
    help="Output using the extended or standard columns",
)
@click.option(
    "--no-nb",
    "nonbonded",
    is_flag=True,
    help="Include nonbonded section in CHARMM parameter file",
)
@click.option(
    "--no-resid",
    "resid",
    is_flag=True,
    help="Include segment IDs in internal coordinate files",
)
@click.option(
    "--no-cmap",
    "cmap",
    is_flag=True,
    help="Include CMAP section in CHARMM PSF file",
)
@click.option(
    "--no-cheq",
    "cheq",
    is_flag=True,
    help="Include charge equilibrium section in CHARMM PSF file",
)
@click.option(
    "--uniform",
    "mass",
    is_flag=True,
    help="Set uniform mass of beads to 1.0",
)
@click.option(
    "--write",
    "write_traj",
    is_flag=True,
    help="Convert the trajectory file",
)
@click.option(
    "-l",
    "--list",
    "model_list",
    is_flag=True,
    help="List available models with their descriptions")
def cli(
        topology,
        trajectory,
        outdir,
        prefix,
        rmin,
        rmax,
        model,
        charmm_version,
        com,
        extended,
        resid,
        cmap,
        cheq,
        nonbonded,
        mass,
        write_traj,
        model_list,
):
    logging.config.dictConfig({
        "version"                 : 1,
        "disable_existing_loggers": False,  # this fixes the problem
        "formatters"              : {
            "standard": {
                "class" : "logging.Formatter",
                "format": "%(name)-12s %(levelname)-8s %(message)s",
            },
            "detailed": {
                "class"  : "logging.Formatter",
                "format" : "%(asctime)s %(name)-15s %(levelname)-8s %(processName)-10s %(message)s",
                "datefmt": "%m-%d-%y %H:%M",
            },
        },
        "handlers"                : {
            "console": {
                "class"    : "logging.StreamHandler",
                "level"    : "INFO",
                "formatter": "standard",
            },
            "file"   : {
                "class"    : "logging.FileHandler",
                "filename" : path.join(outdir, "convert.log"),
                "level"    : "DEBUG",
                "mode"     : "w",
                "formatter": "detailed",
            }
        },
        "root"                    : {
            "level"   : "DEBUG",
            "handlers": ["console", "file"]
        },
    })
    logger = logging.getLogger(__name__)

    if model_list:
        for k, v in iteritems(_DESCRIBE):
            print("{:20}{}".format(k, v))
        return

    kwargs = dict()
    universe = modeller(topology, trajectory, com=com, model=model, **kwargs)

    kwargs.update(dict(
        outdir=outdir,
        prefix=prefix,
        rmin=rmin,
        rmax=rmax,
        charmm_version=charmm_version,
        extended=extended,
        resid=not resid,
        cmap=not cmap,
        cheq=not cheq,
        nonbonded=not nonbonded,
        write_traj=write_traj,
    ))
    if mass:
        logger.info("Setting all bead masses to 1.0.")
        universe.atoms.mass = 1.0
    write_charmm_files(universe, **kwargs)
