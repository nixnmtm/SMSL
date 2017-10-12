# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import collections
import functools
from concurrent import futures
import os
from os import path

from future.builtins import (
    range,
    zip,
)

import click
from fluctmatch.cli import pass_context
from fluctmatch.fluctmatch import utils


@click.command("splittraj", short_help="Split a trajectory using Gromacs or CHARMM.")
@click.option(
    "--gmx / --charmm",
    "gromacs",
    default=True,
    help="Split using either Gromacs or CHARMM"
)
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
    "-d",
    "toppar",
    default=path.join(os.getcwd(), "toppar"),
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Location of CHARMM topology/parameter files",
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
    "-n",
    "index",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Gromacs index file (e.g. ndx)",
)
@click.option(
    "-o",
    "outfile",
    default="aa.xtc",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Trajectory file (e.g. xtc trr dcd)",
)
@click.option(
    "-l",
    "logfile",
    default="split.log",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Log file",
)
@click.option(
    "-t",
    "--system",
    default=0,
    type=click.INT,
    help="System selection based upon Gromacs index file",
)
@click.option(
    "-b",
    "start",
    default=1,
    type=click.INT,
    help="Start time of trajectory",
)
@click.option(
    "-e",
    "stop",
    default=10000,
    type=click.INT,
    help="Stop time of total trajectory",
)
@click.option(
    "-w",
    "window_size",
    default=10000,
    type=click.INT,
    help="Size of each subtrajectory",
)
@pass_context
def cli(ctx, gromacs, toppar, topology, trajectory, index, outfile, logfile, system, start, stop, window_size):
    half_ws = window_size // 2
    Info = collections.namedtuple("Info", "subdir start stop")
    values = zip(range(start, stop+1, half_ws), range(start+window_size, stop+1, half_ws))
    values = (Info(subdir=(y // half_ws) - 1, start=x, stop=y) for x, y in values)

    root0, ext0 = path.split(trajectory)
    root1, ext1 = path.splitext(outfile)
    if gromacs:
        if ext0 == ".dcd" or ext1 == ".dcd":
            raise IOError("The trajectory input and output files must "
                          "not end with '.dcd'.")
        func = functools.partial(
            utils.split_gmx,
            data_dir=ctx.data,
            topology=topology,
            trajectory=trajectory,
            index=index,
            outfile=outfile,
            logfile=logfile,
            system=system,
        )
    else:
        if ext0 != ".dcd" and ext1 != ".dcd":
            raise IOError("The trajectory input and output files must "
                          "end with '.dcd'.")
        func = functools.partial(
            utils.split_charmm,
            data_dir=ctx.data,
            topology=toppar,
            trajectory=trajectory,
            outfile=outfile,
            logfile=logfile,
        )

    # Run multiple instances simultaneously
    with futures.ProcessPoolExecutor() as pool:
        data = [_ for _ in pool.submit(func, values)]
        for _ in futures.as_completed(data):
            _.result()
