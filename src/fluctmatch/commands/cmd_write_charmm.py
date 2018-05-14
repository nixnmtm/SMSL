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

import logging
import logging.config

import click
import MDAnalysis as mda
from fluctmatch.fluctmatch import utils as fmutils


@click.command("write_charmm", short_help="Write various CHARMM files.")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    help="Gromacs topology file (e.g., tpr gro g96 pdb brk ent)",
)
@click.option(
    "-f",
    "trajectory",
    metavar="FILE",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-o",
    "outdir",
    metavar="DIR",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Trajectory file (e.g. xtc trr dcd)",
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
    "--no-nb",
    "nonbonded",
    is_flag=True,
    help="Include nonbonded section in CHARMM parameter file",
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
        charmm_version,
        extended,
        cmap,
        cheq,
        nonbonded,
        write_traj,
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
                "format" : "%(asctime)s %(name)-15s %(levelname)-8s %(message)s",
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
                "filename" : path.join(outdir, "write_charmm.log"),
                "level"    : "INFO",
                "mode"     : "w",
                "formatter": "detailed",
            }
        },
        "root"                    : {
            "level"   : "INFO",
            "handlers": ["console", "file"]
        },
    })
    logger = logging.getLogger(__name__)

    kwargs = dict(
        outdir=outdir,
        prefix=prefix,
        charmm_version=charmm_version,
        extended=extended,
        cmap=not cmap,
        cheq=not cheq,
        nonbonded=not nonbonded,
        write_traj=write_traj,
    )
    universe = mda.Universe(topology, trajectory)
    fmutils.write_charmm_files(universe, **kwargs)
