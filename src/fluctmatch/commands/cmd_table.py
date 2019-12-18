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
from MDAnalysis.lib.util import filename
from fluctmatch.analysis import paramtable


@click.command(
    "table", short_help="Create a table from individual parameter files.")
@click.option(
    "-d",
    "--datadir",
    "data_dir",
    metavar="DIR",
    default=path.join(os.getcwd(), "data"),
    show_default=True,
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Data directory",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    show_default=True,
    default=path.join(os.getcwd(), "table.log"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Log file",
)
@click.option(
    "-o",
    "--outdir",
    metavar="OUTDIR",
    default=os.getcwd(),
    show_default=True,
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Directory",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default="cg",
    show_default=True,
    type=click.STRING,
    help="Prefix for filenames",
)
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
@click.option(
    "-r",
    "--ressep",
    metavar="RESSEP",
    default=3,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Number of residues to exclude in I,I+r")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
)
@click.option(
    "-s",
    "--start",
    metavar="START",
    default=None,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Start from given window directory")
@click.option(
    "-e",
    "--end",
    metavar="END",
    default=None,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="End at given window directory")
def cli(data_dir, logfile, outdir, prefix, tbltype, ressep, verbose, start, end):
    pt = paramtable.ParamTable(
        prefix=prefix,
        tbltype=tbltype,
        ressep=ressep,
        datadir=data_dir,
        start=start,
        end=end
    )
    # Setup logger
    logging.config.dictConfig({
        "version": 1,
        "disable_existing_loggers": False,  # this fixes the problem
        "formatters": {
            "standard": {
                "class": "logging.Formatter",
                "format": "%(name)-12s %(levelname)-8s %(message)s",
            },
            "detailed": {
                "class": "logging.Formatter",
                "format":
                "%(asctime)s %(name)-15s %(levelname)-8s %(message)s",
                "datefmt": "%m-%d-%y %H:%M",
            },
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "standard",
            },
            "file": {
                "class": "logging.FileHandler",
                "filename": logfile,
                "level": "INFO",
                "mode": "w",
                "formatter": "detailed",
            }
        },
        "root": {
            "level": "INFO",
            "handlers": ["console", "file"]
        },
    })
    logger = logging.getLogger(__name__)

    pt.run(verbose=verbose)

    # Write the various tables to different files.
    fn = path.join(outdir, filename(tbltype.lower(), ext="txt", keep=True))
    pt.write(fn)

    if tbltype == "Kb":
        fn = path.join(outdir, filename("perres", ext="txt"))
        with open(fn, mode="wb") as output:
            logger.info("Writing per-residue data to {}.".format(fn))
            table = pt.per_residue.to_csv(
                header=True,
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
                encoding="utf-8",
            )
            output.write(table.encode())
            logger.info("Table successfully written.")

        fn = path.join(outdir, filename("interactions", ext="txt"))
        with open(fn, mode="wb") as output:
            logger.info("Writing interactions to {}.".format(fn))
            table = pt.interactions.to_csv(
                header=True,
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
                encoding="utf-8",
            )
            output.write(table.encode())
            logger.info("Table successfully written.")
