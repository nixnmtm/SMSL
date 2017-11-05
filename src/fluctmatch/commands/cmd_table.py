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
from future.builtins import open
from future.utils import native_str
from MDAnalysis.lib.util import filename

from fluctmatch.analysis import paramtable


@click.command("table", short_help="Create a table from individual parameter files.")
@click.option(
    "-d",
    "--datadir",
    "data_dir",
    metavar="DIR",
    default=path.join(os.getcwd(), "data"),
    type=click.Path(
        exists=False,
        file_okay=False,
        resolve_path=True
    ),
    help="Data directory",
)
@click.option(
    "-o",
    "--outdir",
    metavar="OUTDIR",
    default=os.getcwd(),
    type=click.Path(
        exists=False,
        file_okay=False,
        resolve_path=True
    ),
    help="Directory",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default="cg",
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-t",
    "--type",
    "tbltype",
    metavar="TYPE",
    default="Kb",
    type=click.Choice(["Kb", "b0"]),
    help="Force constant or equilibrium distance",
)
@click.option(
    "-r",
    "--ressep",
    metavar="RESSEP",
    default=3,
    type=click.INT,
    help="Separation between residues (I,I+n)"
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
)
def cli(data_dir, outdir, prefix, tbltype, ressep, verbose):
    pt = paramtable.ParamTable(prefix=prefix, tbltype=tbltype, ressep=ressep)
    pt.run(verbose=verbose)

    # Write the various tables to different files.
    fn = path.join(outdir, filename(tbltype.lower(), ext="txt", keep=True))
    pt.write(fn)

    if tbltype == "Kb":
        fn = path.join(outdir, filename("perres", ext="txt"))
        with open(fn, mode="wb") as output:
            table = pt.per_residue.to_csv(
                header=True,
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
                encoding="utf-8",
            )
            output.write(table.encode())

        fn = path.join(outdir, filename("interactions", ext="txt"))
        with open(fn, mode="wb") as output:
            table = pt.interactions.to_csv(
                header=True,
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
                encoding="utf-8",
            )
            output.write(table.encode())
