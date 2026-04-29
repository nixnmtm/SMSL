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

import functools
import logging
import logging.config
import multiprocessing as mp
import os
from os import path

import click
from MDAnalysis.lib import util as mdutil
from fluctmatch.fluctmatch import utils

_CONVERT = {
    "MDA": utils.split_mda,
    "GMX": utils.split_gmx,
}

# IF MDA IS USED remember
# TODO: split_mda currently selects windows using ts.time, while the CLI
# arguments -b/-e/-w are documented as frame-based. This works for the
# current trajectory because 1 frame = 1 ps, but should be refactored to
# use ts.frame (or the CLI should be renamed to explicit time units) so
# behavior is consistent for trajectories saved at other time intervals.


@click.command(
    "splittraj", short_help="Split a trajectory using Gromacs or CHARMM, "
                            "the trajectory input expects to have atleast 5 floating point precision (i.e.) -ndec=5")
@click.option(
    "--type",
    "program",
    type=click.Choice(_CONVERT.keys()),
    default="MDA",
    help="Split using MDAnalysis or Gromacs")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.tpr"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Gromacs topology file (e.g., tpr gro g96 pdb brk ent)",
)
@click.option(
    "-f",
    "trajectory",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.xtc"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "--data",
    metavar="DIR",
    default=path.join(os.getcwd(), "data"),
    type=click.Path(
        exists=False,
        file_okay=False,
        writable=True,
        readable=True,
        resolve_path=True,
    ),
    help="Directory to write data.",
)
@click.option(
    "-n",
    "index",
    metavar="FILE",
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Gromacs index file (e.g. ndx)",
)
@click.option(
    "-o",
    "outfile",
    metavar="FILE",
    default="aa.xtc",
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=False,
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    show_default=True,
    default="splittraj.log",
    type=click.Path(exists=False, file_okay=True, resolve_path=False),
    help="Log file",
)
@click.option(
    "-t",
    "--system",
    metavar="NDXNUM",
    default=0,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="System selection based upon Gromacs index file",
)
@click.option(
    "-b",
    "start",
    metavar="FRAME",
    default=1,
    show_default=True,
    type=click.IntRange(1, None, clamp=True),
    help="Start time of trajectory",
)
@click.option(
    "-e",
    "stop",
    metavar="FRAME",
    default=10000,
    show_default=True,
    type=click.IntRange(1, None, clamp=True),
    help="Stop time of total trajectory",
)
@click.option(
    "-w",
    "window_size",
    metavar="WINSIZE",
    default=10000,
    show_default=True,
    type=click.IntRange(2, None, clamp=True),
    help="Size of each subtrajectory",
)
@click.option(
    "-p",
    "precision",
    metavar="PRECISION",
    default=5,
    show_default=True,
    type=click.IntRange(1, None, clamp=True),
    help="XTC number of decimal precision",
)
@click.option(
    "--n-jobs",
    "n_jobs",
    metavar="N",
    default=32,
    show_default=True,
    type=click.IntRange(1, None, clamp=True),
    help="Number of worker processes to use for trajectory splitting.",
)




def cli(program, topology, trajectory, data, index, outfile, logfile,
        system, start, stop, window_size, precision, n_jobs):
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
                "filename": path.join(os.getcwd(), logfile),
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

    if program == "GMX" and mdutil.which("gmx") is None:
        logger.error("Gromacs 5.0+ is required. "
                     "If installed, please ensure that it is in your path.")
        raise OSError("Gromacs 5.0+ is required. "
                      "If installed, please ensure that it is in your path.")

    half_size = window_size // 2
    beg = start - half_size if start >= window_size else start
    values = zip(
        range(beg, stop + 1, half_size),
        range(beg + window_size - 1, stop + 1, half_size))
    values = [(i + 1, x, y) for i, (x, y) in enumerate(values)]

    func = functools.partial(
        _CONVERT[program],
        data_dir=data,
        topology=topology,
        trajectory=trajectory,
        index=index,
        outfile=outfile,
        logfile=logfile,
        system=system,
        precision=precision,
    )

    logger.info("first 10 split values: %s", values[:10])
    logger.info("total split windows: %d", len(values))

    # Run multiple instances simultaneously
    with mp.Pool(processes=n_jobs) as pool:
        pool.map(func, values)
    logger.info("splittraj finished successfully | total_windows=%d", len(values))
