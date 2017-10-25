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
from future.utils import native_str
from MDAnalysis.lib.util import openany

from fluctmatch.analysis.paramstats import ParamStats
from fluctmatch.analysis.paramtable import ParamTable


@click.command(
    "diff",
    short_help="Calculate differences between two tables."
)
@click.option(
    "-s",
    "--stats",
    is_flag=True,
    help="Calculate statistics of tables",
)
@click.option(
    "-g",
    "--hist",
    is_flag=True,
    help="Calculate histograms of tables",
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
@click.option(
    "-t",
    "--type",
    "tbltype",
    metavar="TYPE",
    default="Kb",
    type=click.Choice(["Kb", "b0"]),
    help="Force constant or equilibrium distance",
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
def cli(stats, hist, outdir, ressep, tbltype, table):
    pt = ParamTable(ressep=ressep)
    pt.from_file(table)
    ps = ParamStats(pt)

    if stats:
        filename = path.join(
            outdir,
            "_".join(tbltype.lower(), "table", "stats.txt")
        )
        with openany(filename, "w") as stat_file:
            print("Writing table statistics to {}".format(filename))
            info = ps.table_stats()
            info.to_csv(
                stat_file,
                header=True,
                index=True,
                float_format="%.4f",
                sep=native_str(" ")
            )

        if tbltype == "Kb":
            filename = path.join(outdir, "interaction_stats.txt")
            with openany(filename, "w") as stat_file:
                print(
                    "Writing residue-residue statistics to {}".format(filename)
                )
                ps._table._ressep = 0
                info = ps.interaction_stats()
                info.to_csv(
                    stat_file,
                    header=True,
                    index=True,
                    float_format="%.4f",
                    sep=native_str(" ")
                )

            filename = path.join(outdir, "residue_stats.txt")
            with openany(filename, "w") as stat_file:
                print(
                    "Writing individual residue statistics "
                    "to {}".format(filename)
                )
                ps._table._ressep = ressep
                info = ps.residue_stats()
                info.to_csv(
                    stat_file,
                    header=True,
                    index=True,
                    float_format="%.4f",
                    sep=native_str(" ")
                )

    if hist:
        filename = path.join(
            outdir,
            "_".join(tbltype.lower(), "table", "hist.txt")
        )
        with openany(filename, "w") as stat_file:
            print("Writing table histogram to {}".format(filename))
            info = ps.table_hist()
            info.to_csv(
                stat_file,
                index=True,
                float_format="%.4f",
                sep=native_str(" ")
            )

        if tbltype == "Kb":
            filename = path.join(outdir, "interaction_hist.txt")
            with openany(filename, "w") as stat_file:
                print(
                    "Writing residue-residue histogram to {}".format(filename)
                )
                ps._table._ressep = 0
                info = ps.interaction_hist()
                info.to_csv(
                    stat_file,
                    index=True,
                    float_format="%.4f",
                    sep=native_str(" ")
                )

            filename = path.join(outdir, "residue_hist.txt")
            with openany(filename, "w") as stat_file:
                print(
                    "Writing individual residue histogram "
                    "to {}".format(filename)
                )
                ps._table._ressep = ressep
                info = ps.residue_hist()
                info.to_csv(
                    stat_file,
                    index=True,
                    float_format="%.4f",
                    sep=native_str(" ")
                )
