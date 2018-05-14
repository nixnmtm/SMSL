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
from future.builtins import open
from future.utils import native_str

import logging
import logging.config
import os
from os import path

import click
from fluctmatch.analysis.paramstats import ParamStats
from fluctmatch.analysis.paramtable import ParamTable


@click.command("stats", short_help="Calculate statistics of a table.")
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
    show_default=True,
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
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Number of residues to exclude in I,I+r")
@click.option(
    "-t",
    "--type",
    "tbltype",
    metavar="TYPE",
    default="Kb",
    show_default=True,
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
    # Setup logger
    logging.config.dictConfig({
        "version"                 : 1,
        "disable_existing_loggers": False,  # this fixes the problem
        "formatters"              : {
            "standard": {
                "class" : "logging.Formatter",
                "format": "%(name)-12s %(levelname)-8s %(message)s",
            },
            "detailed": {
                "class"  : "logging.Formatter",
                "format" : "%(asctime)s %(name)-15s %(levelname)-8s %(processName)-10s %(message)s",
                "datefmt": "%m-%d-%y %H:%M",
            },
        },
        "handlers"                : {
            "console": {
                "class"    : "logging.StreamHandler",
                "level"    : "INFO",
                "formatter": "standard",
            },
            "file"   : {
                "class"    : "logging.FileHandler",
                "filename" : path.join(outdir, "stats.log"),
                "level"    : "DEBUG",
                "mode"     : "w",
                "formatter": "detailed",
            }
        },
        "root"                    : {
            "level"   : "DEBUG",
            "handlers": ["console", "file"]
        },
    })
    logger = logging.getLogger(__name__)

    logger.info("Reading {}".format(table))
    pt = ParamTable(ressep=ressep)
    pt.from_file(table)
    ps = ParamStats(pt)

    if stats:
        filename = path.join(outdir, "_".join((tbltype.lower(), "table",
                                               "stats.txt")))
        with open(filename, mode="wb") as stat_file:
            logger.info("Writing table statistics to {}".format(filename))
            info = ps.table_stats().to_csv(
                header=True,
                index=True,
                float_format="%.4f",
                sep=" ",
                encoding="utf-8",
            )
            stat_file.write(info.encode())
            logger.debug("Table successfully written.")

        if tbltype == "Kb":
            filename = path.join(outdir, "interaction_stats.txt")
            with open(filename, mode="wb") as stat_file:
                logger.info("Writing residue-residue statistics to {}".format(
                    filename))
                ps._table._ressep = 0
                info = ps.interaction_stats().to_csv(
                    header=True,
                    index=True,
                    sep=native_str(" "),
                    float_format=native_str("%.4f"),
                    encoding="utf-8",
                )
                stat_file.write(info.encode())
                logger.debug("Table successfully written.")

            filename = path.join(outdir, "residue_stats.txt")
            with open(filename, mode="wb") as stat_file:
                logger.info(
                    "Writing individual residue statistics to {}".format(
                        filename))
                ps._table._ressep = ressep
                info = ps.residue_stats().to_csv(
                    header=True,
                    index=True,
                    sep=native_str(" "),
                    float_format=native_str("%.4f"),
                    encoding="utf-8",
                )
                stat_file.write(info.encode())
                logger.debug("Table successfully written.")

    if hist:
        filename = path.join(outdir, "_".join((tbltype.lower(), "table",
                                               "hist.txt")))
        with open(filename, mode="wb") as stat_file:
            logger.info("Writing table histogram to {}".format(filename))
            info = ps.table_hist().to_csv(
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
                encoding="utf-8",
            )
            stat_file.write(info.encode())
            logger.debug("Table successfully written.")

        if tbltype == "Kb":
            filename = path.join(outdir, "interaction_hist.txt")
            with open(filename, mode="wb") as stat_file:
                logger.info(
                    "Writing residue-residue histogram to {}".format(filename))
                ps._table._ressep = 0
                info = ps.interaction_hist().to_csv(
                    index=True,
                    sep=native_str(" "),
                    float_format=native_str("%.4f"),
                    encoding="utf-8",
                )
                stat_file.write(info.encode())
                logger.debug("Table successfully written.")

            filename = path.join(outdir, "residue_hist.txt")
            with open(filename, mode="wb") as stat_file:
                logger.info(
                    "Writing individual residue histogram to {}".format(
                        filename))
                ps._table._ressep = ressep
                info = ps.residue_hist().to_csv(
                    index=True,
                    sep=native_str(" "),
                    float_format=native_str("%.4f"),
                    encoding="utf-8",
                )
                stat_file.write(info.encode())
                logger.debug("Table successfully written.")
