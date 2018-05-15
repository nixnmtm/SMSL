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
import numpy as np
from fluctmatch.analysis import paramtable

logger = logging.getLogger(__name__)


@click.command(
    "normdiff", short_help="Normalize the differences between t and t - dt/2.")
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    show_default=True,
    default=path.join(os.getcwd(), "normdiff.log"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Log file",
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
@click.argument(
    "kb",
    metavar="kb_table",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True,
    ),
)
@click.argument(
    "b0",
    metavar="b0_table",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True,
    ),
)
def cli(logfile, outdir, ressep, kb, b0):
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

    resgrp = ['segidI', 'resI']

    logger.info("Loading coupling strength table.")
    kb_table = paramtable.ParamTable(ressep=ressep)
    kb_table.from_file(kb)
    kb_table = kb_table.table.copy(deep=True)
    logger.info("Table loaded successfully.")

    logger.info("Loading distance table.")
    b0_table = paramtable.ParamTable(ressep=ressep)
    b0_table.from_file(b0)
    b0_table = b0_table.table.copy(deep=True)
    logger.info("Table loaded successfully.")

    idx = (kb_table == 0.0)
    maxkb = np.maximum(kb_table[kb_table.columns[1:].values],
                       kb_table[kb_table.columns[:-1]].values)

    maxkb[maxkb == 0.0] = np.NaN
    kb_table = kb_table.diff(axis=1).dropna(axis=1) / maxkb
    kb_table.columns = kb_table.columns.astype(np.int32) - 1

    count = kb_table[kb_table.abs() > 0.].groupby(level=resgrp).count()
    kb_table = kb_table.groupby(level=resgrp).sum() / count
    kb_table.fillna(0., inplace=True)

    # Calculate the r.m.s.d. of equilibrium distances between sites.
    b0_table = b0_table.diff(axis=1).dropna(axis=1) / maxkb
    b0_table[idx] = 0.0
    b0_table = 0.5 * b0_table.pow(2).groupby(level=resgrp).sum()
    b0_table = b0_table.apply(np.sqrt)

    filename = path.join(outdir, "normed_kb.txt")
    with open(filename, mode="wb") as output:
        logger.info("Writing normed coupling strengths to {}".format(filename))
        kb_table = kb_table.to_csv(
            header=True,
            index=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(kb_table.encode())
        logger.info("Table written successfully.")

    filename = path.join(outdir, "normed_b0.txt")
    with open(filename, mode="wb") as output:
        logger.info("Writing normed distances to {}".format(filename))
        b0_table = b0_table.to_csv(
            header=True,
            index=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(b0_table.encode())
        logger.info("Table written successfully.")
