# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
from future.utils import native_str
from MDAnalysis.lib.util import (filename, openany)

from fluctmatch.analysis import paramtable


@click.command("table", short_help="Create a table from individual parameter files.")
@click.option(
    "-d",
    "--datadir"
    "data_dir",
    default=path.join(os.getcwd(), "data"),
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Data directory",
)
@click.option(
    "-o",
    "--outdir",
    default=os.getcwd(),
    type=click.Path(
        exists=True,
        file_okay=False,
        resolve_path=True
    ),
    help="Directory",
)
@click.option(
    "-p",
    "--prefix",
    default="cg",
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-t",
    "--type",
    "tbltype",
    default="Kb",
    type=click.Choice(["Kb", "b0"]),
    help="Force constant or equilibrium distance",
)
@click.option(
    "-r",
    "--ressep",
    default=3,
    type=click.INT,
    help="Separation between residues (I,I+n)"
)
@click.option(
    "--verbose",
    is_flag=True,
)

def table(data_dir, outdir, prefix, tbltype, ressep, verbose):
    pt = paramtable.ParamTable(data_dir, prefix=prefix, tbltype=tbltype, ressep=ressep)
    pt.run(verbose=verbose)

    # Write the various tables to different files.
    with openany(path.join(outdir, filename(tbltype.lower(), ext="txt"))) as table:
        pt.table.to_csv(table, header=True, index=True, sep=native_str(" "), float_format="%.4f")
    if tbltype == "Kb":
        with openany(path.join(outdir, filename(".".join((tbltype.lower()), "perres"), ext="txt"))) as table:
            pt.per_residue.to_csv(table, header=True, index=True, sep=native_str(" "), float_format="%.4f")
        with openany(path.join(outdir, filename(".".join((tbltype.lower()), "interactions"), ext="txt"))) as table:
            pt.interactions.to_csv(table, header=True, index=True, sep=native_str(" "), float_format="%.4f")
