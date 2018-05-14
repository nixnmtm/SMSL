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
import os
from os import path

import click
from MDAnalysis.lib.util import which
from fluctmatch.analysis import thermodynamics


@click.command("thermo", short_help="Calculate thermodynamic properties.")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    default="fluctmatch.xplor.psf",
    show_default=True,
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=False,
    ),
    help="Topology file (e.g., tpr gro g96 pdb brk ent psf)",
)
@click.option(
    "-f",
    "trajectory",
    metavar="FILE",
    default="cg.dcd",
    show_default=True,
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=False,
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-d",
    "datadir",
    metavar="DIR",
    default=path.join(os.getcwd(), "data"),
    show_default=True,
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Directory",
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
    "-c",
    "--charmm",
    "charmm_version",
    metavar="VERSION",
    default=41,
    show_default=True,
    type=click.IntRange(27, None, clamp=True),
    help="CHARMM version",
)
def cli(datadir, outdir, topology, trajectory, nma_exec, temperature,
        charmm_version):
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
                "filename" : path.join(outdir, "thermo.log"),
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

    # Attempt to create the necessary subdirectory
    try:
        os.makedirs(outdir)
    except OSError:
        pass

    thermodynamics.create_thermo_tables(
        datadir,
        outdir,
        topology=topology,
        trajectory=trajectory,
        temperature=temperature,
        nma_exec=nma_exec,
        charmm_version=charmm_version)
