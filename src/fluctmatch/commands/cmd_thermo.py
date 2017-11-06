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
from MDAnalysis.lib.util import which
from fluctmatch.analysis import thermodynamics


@click.command("thermo", short_help="Calculate thermodynamic properties.")
@click.option(
    "-s",
    "topology",
    metavar="FILE",
    default="fluctmatch.xplor.psf",
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
    default="cg.xtc",
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
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Directory",
)
@click.option(
    "-o",
    "outdir",
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
    "-e",
    "--exec",
    "nma_exec",
    metavar="FILE",
    envvar="CHARMMEXEC",
    default=which("charmm"),
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=True
    ),
    help="CHARMM executable file",
)
@click.option(
    "-t",
    "--temperature",
    metavar="TEMP",
    type=click.FLOAT,
    default=300.0,
    help="Temperature of simulation",
)
@click.option(
    "-c",
    "--charmm",
    "charmm_version",
    metavar="VERSION",
    default=41,
    type=click.IntRange(27, None, clamp=True),
    help="CHARMM version",
)
def cli(
    datadir, outdir, topology, trajectory, nma_exec, temperature, charmm_version
):
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
        charmm_version=charmm_version
    )
