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
from io import TextIOWrapper
from os import path

import click
from future.builtins import open
from future.utils import native_str
from MDAnalysis.lib.util import (filename, openany)

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
    pt = paramtable.ParamTable(data_dir, prefix=prefix, tbltype=tbltype, ressep=ressep)
    pt.run(verbose=verbose)

    # Write the various tables to different files.
    fn = path.join(outdir, filename(tbltype.lower(), ext="txt", keep=True))
    with open(
        fn, "wb"
    ) as table, TextIOWrapper(table, encoding="utf-8") as buf:
        pt.table.to_csv(
            buf,
            header=True,
            index=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
        )
    if tbltype == "Kb":
        fn = path.join(outdir, filename("perres", ext="txt"))
        with open(
            fn, "wb"
        ) as table, TextIOWrapper(table, encoding="utf-8") as buf:
            pt.table.to_csv(
                buf,
                header=True,
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
            )
        fn = path.join(outdir, filename("interactions", ext="txt"))
        with open(
            fn, "wb"
        ) as table, TextIOWrapper(table, encoding="utf-8") as buf:
            pt.table.to_csv(
                buf,
                header=True,
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
            )
