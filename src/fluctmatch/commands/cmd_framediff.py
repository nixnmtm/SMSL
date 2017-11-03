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
from MDAnalysis.lib.util import openany

from fluctmatch.analysis import paramtable


@click.command(
    "framediff",
    short_help="Calculate differences between consecutive frames."
)
@click.option(
    "-o",
    "--outdir",
    metavar="OUTDIR",
    default=os.getcwd(),
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True,
    ),
    help="Directory",
)
@click.option(
    "-r",
    "--ressep",
    metavar="RESSEP",
    default=3,
    type=click.INT,
    help="Separation between residues (I,I+n)"
)
@click.argument(
    "table",
    metavar="TABLE",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True,
    ),
)
def cli(outdir, ressep, table):
    """Calculate the differences between consecutive frames.

    Parameters
    ----------
    outdir : str
        Output directory for files
    ressep : int
        Separation between residues
    table : str
        Filename of table
    """
    table_1 = paramtable.ParamTable(ressep=ressep)
    table_1.from_file(table)

    d_table = table_1.diff(axis=1).dropna(axis=1)
    d_perres = table_1.per_residue.diff(axis=1).dropna(axis=1)
    d_interactions = table_1.interactions.diff(axis=1).dropna(axis=1)

    filename = path.join(outdir, "dframe_coupling.txt")
    with openany(filename, mode="w") as output:
        print("Writing frame differences to {}".format(filename))
        d_table.to_csv(
            output,
            header=True,
            index=True,
            sep=" ",
            float_format="%.4f",
        )

    filename = path.join(outdir, "dframe_perres.txt")
    with openany(filename, mode="w") as output:
        print("Writing per residue frame differences to {}".format(filename))
        d_perres.to_csv(
            output,
            header=True,
            index=True,
            sep=" ",
            float_format="%.4f",
        )

    filename = path.join(outdir, "dframe_interactions.txt")
    with openany(filename, mode="w") as output:
        print(
            "Writing residue-residue frame differences to {}".format(filename)
        )
        d_interactions.to_csv(
            output,
            header=True,
            index=True,
            sep=" ",
            float_format="%.4f",
        )
