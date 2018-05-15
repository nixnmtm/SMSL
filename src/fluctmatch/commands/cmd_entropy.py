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
from fluctmatch.analysis import entropy


@click.command(
    "entropy", short_help="Calculate the Shannon entropy of residues.")
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
    help="Separation between residues (I,I+n)")
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
                "filename": path.join(outdir, "entropy.log"),
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

    logger.info("Loading {}".format(table))
    ent_table = entropy.Entropy(table, ressep=ressep)

    filename = path.join(outdir, "coupling.entropy.txt")
    with open(filename, mode="wb") as output:
        logger.info("Writing coupling entropy to {}".format(filename))
        ent = ent_table.coupling_entropy().to_csv(
            index=True,
            header=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(ent.encode())
        logger.info("Table written successfully.")

    filename = path.join(outdir, "relative.entropy.txt")
    with open(filename, mode="wb") as output:
        logger.info("Writing relative entropy to {}".format(filename))
        ent = ent_table.relative_entropy().to_csv(
            index=True,
            header=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(ent.encode())
        logger.info("Table written successfully.")

    filename = path.join(outdir, "windiff.entropy.txt")
    with open(filename, mode="wb") as output:
        logger.info(
            "Writing entropy for window difference to {}".format(filename))
        ent = ent_table.windiff_entropy().to_csv(
            index=True,
            header=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(ent.encode())
        logger.info("Table written successfully.")
