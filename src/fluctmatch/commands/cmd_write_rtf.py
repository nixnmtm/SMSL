# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# pysca --- https://github.com/tclick/python-pysca
# Copyright (c) 2015-2017 The pySCA Development Team and contributors
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

from future.utils import (native_str)

import logging
import logging.config
import os
from os import path

import click
import MDAnalysis as mda
from fluctmatch.topology import RTF


@click.command(
    "write_rtf", short_help="Create an RTF file from a structure file.")
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
    default=path.join(os.getcwd(), "md.pdb"),
    show_default=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Trajectory file (e.g. xtc trr dcd, crd, cor, pdb, tpr)",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    show_default=True,
    default=path.join(os.getcwd(), "write_rtf.log"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Log file",
)
@click.option(
    "-o",
    "--outfile",
    metavar="FILE",
    show_default=True,
    default=path.join(os.getcwd(), "md.rtf"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="CHARMM topology file",
)
@click.option(
    "--no-decl",
    "decl",
    is_flag=True,
    help="Include declaration section in CHARMM topology file",
)
@click.option(
    "--uniform",
    "mass",
    is_flag=True,
    help="Set uniform mass of beads to 1.0",
)
def cli(
        topology,
        trajectory,
        logfile,
        outfile,
        decl,
        mass,
):
    logging.config.dictConfig({
        "version": 1,
        "disable_existing_loggers": False,  # this fixes the problem
        "formatters": {
            "standard": {
                "class": "logging.Formatter",
                "format": "%(name)-12s %(levelname)-8s %(message)s",
            },
            "detailed": {
                "class": "logging.Formatter",
                "format":
                "%(asctime)s %(name)-15s %(levelname)-8s %(message)s",
                "datefmt": "%m-%d-%y %H:%M",
            },
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "standard",
            },
            "file": {
                "class": "logging.FileHandler",
                "filename": logfile,
                "level": "INFO",
                "mode": "w",
                "formatter": "detailed",
            }
        },
        "root": {
            "level": "INFO",
            "handlers": ["console", "file"]
        },
    })
    logger = logging.getLogger(__name__)

    kwargs = dict()
    universe = mda.Universe(topology, trajectory)

    if mass:
        logger.info("Setting all bead masses to 1.0.")
        universe.atoms.mass = 1.0

    with mda.Writer(native_str(outfile), **kwargs) as rtf:
        logger.info("Writing {}...".format(outfile))
        rtf.write(universe, decl=not decl)
