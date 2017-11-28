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

import os
from os import path

import click
import MDAnalysis as mda
import pandas as pd
from future.builtins import (dict, open, zip)


@click.command(
    "table_convert",
    short_help="Transform an ENM IC table name to corresponding atoms."
)
@click.option(
    "-s1",
    "top1",
    metavar="FILE",
    default=path.join(os.getcwd(), "cg.xplor.psf"),
    show_default=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Topology file",
)
@click.option(
    "-s2",
    "top1",
    metavar="FILE",
    default=path.join(os.getcwd(), "fluctmatch.xplor.psf"),
    show_default=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Topology file",
)
@click.option(
    "-f1",
    "coord1",
    metavar="FILE",
    default=path.join(os.getcwd(), "cg.cor"),
    show_default=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Coordinate file",
)
@click.option(
    "-f2",
    "coord2",
    metavar="FILE",
    default=path.join(os.getcwd(), "fluctmatch.cor"),
    show_default=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Coordinate file",
)
@click.option(
    "-t",
    "--table",
    metavar="FILE",
    default=path.join(os.getcwd(), "kb.txt"),
    show_default=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True
    ),
    help="Coordinate file",
)
@click.option(
    "-o",
    "--outfile",
    metavar="OUTFILE",
    default=path.join(os.getcwd(), "kb_aa.txt"),
    show_default=True,
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=True
    ),
    help="Table file",
)
def cli(top1, top2, coord1, coord2, table, outfile):
    cg = mda.Universe(top1, coord1)
    fluctmatch = mda.Universe(top2, coord2)
    convert = dict(zip(fluctmatch.atoms.names, cg.atoms.names))
    resnames = dict(zip(cg.residue.resnums, cg.residues.resnames))

    with open(table, "rb") as tbl:
        constants = pd.read_csv(
            tbl,
            header=0,
            skipinitialspace=True,
            delim_whitespace=True,
        )

    # Transform assigned bead names to an all-atom designation.
    constants["I"] = constants["I"].apply(
        lambda x: convert[x]
    )
    constants["J"] = constants["J"].apply(
        lambda x: convert[x]
    )

    # Create lists of corresponding residues
    resnI = constants["resI"].apply(lambda x: resnames[x]).to_frame()
    resnJ = constants["resJ"].apply(lambda x: resnames[x]).to_frame()

    # Concatenate the columns
    cols = ["segidI", "resI", "resnI", "I", "segidJ", "resJ", "resnJ", "J"]
    constants = pd.concat([constants, resnI, resnJ], axis=1)
    constants.set_index(cols, inplace=True)

    with open(outfile, "wb") as output:
        constants = constants.to_csv(
            header=True,
            index=True,
            sep=" ",
            float_format="%.4f",
        )
        output.write(constants.encode())
