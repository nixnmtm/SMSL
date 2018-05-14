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
from os import path

import click
import MDAnalysis as mda


@click.command("wordom", short_help="Convert between coordinate file types.")
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
    "outfile",
    metavar="FILE",
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-b",
    "start",
    metavar="FRAME",
    default=0,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Start time of trajectory",
)
@click.option(
    "-e",
    "stop",
    metavar="FRAME",
    type=click.IntRange(0, None, clamp=True),
    help="Stop time of total trajectory",
)
@click.option(
    "-dt",
    "step",
    metavar="STEP",
    default=1,
    show_default=True,
    type=click.IntRange(1, None, clamp=True),
    help="Interval between frames")
@click.option("-v", "--verbose", is_flag=True, help="Show progress of output")
def cli(topology, trajectory, outfile, start, stop, step, verbose):
    # Setup logger
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
                "filename" : path.join(path.dirname(outfile), "wordom.log"),
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

    logger.info("Loading the universe.")
    universe = mda.Universe(topology, trajectory)
    trajectory = universe.trajectory
    with mda.Writer(
            outfile, n_atoms=universe.atoms.n_atoms, dt=trajectory.dt) as trj:
        logger.info("Converting from {} to {}".format(trajectory.format,
                                                      trj.format))
        if verbose:
            with click.progressbar(
                    length=trajectory.n_frames,
                    label="Trajectory frames written") as bar:
                for ts in trajectory[start:stop:step]:
                    trj.write(ts)
                    bar.update(ts.frame)
        else:
            for ts in trajectory[start:stop:step]:
                trj.write(ts)
        logger.debug("Trajectory conversion successful.")
