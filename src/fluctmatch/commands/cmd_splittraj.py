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

import functools
import multiprocessing as mp
import os
from os import path

from future.builtins import (
    dict,
    range,
    zip,
)
from future.utils import viewkeys

import click
from fluctmatch.fluctmatch import utils

_CONVERT = dict(
    GMX=utils.split_gmx,
    CHARMM=utils.split_charmm,
)


@click.command(
    "splittraj", short_help="Split a trajectory using Gromacs or CHARMM.")
@click.option(
    "--type",
    "program",
    type=click.Choice(viewkeys(_CONVERT)),
    default="GMX",
    help="Split using an external MD program")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    default=path.join(os.getcwd(), "md.tpr"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Gromacs topology file (e.g., tpr gro g96 pdb brk ent)",
)
@click.option(
    "--toppar",
    metavar="DIR",
    default=path.join(os.getcwd(), "toppar"),
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Location of CHARMM topology/parameter files",
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
    "logfile",
    metavar="FILE",
    default="split.log",
    show_default=True,
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=False,
    ),
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
def cli(program, toppar, topology, trajectory, data, index, outfile, logfile,
        system, start, stop, window_size):
    half_size = window_size // 2
    beg = start - half_size if start >= window_size else start
    values = zip(
        range(beg, stop + 1, half_size),
        range(beg + window_size - 1, stop + 1, half_size))
    values = [((y // half_size) - 1, x, y) for x, y in values]

    func = functools.partial(
        _CONVERT[program],
        data_dir=data,
        topology=topology,
        toppar=toppar,
        trajectory=trajectory,
        index=index,
        outfile=outfile,
        logfile=logfile,
        system=system,
    )

    # Run multiple instances simultaneously
    pool = mp.Pool()
    pool.map_async(func, values)
    pool.close()
    pool.join()
