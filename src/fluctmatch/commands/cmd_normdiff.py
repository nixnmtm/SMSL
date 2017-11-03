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
import numpy as np
from MDAnalysis.lib.util import openany

from fluctmatch.analysis import paramtable


@click.command(
    "normdiff",
    short_help="Normalize the differences between t and t - dt/2."
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
def cli(outdir, ressep, kb, b0):
    resgrp = ['segidI', 'resI']

    kb_table = paramtable.ParamTable(ressep=ressep)
    kb_table.from_file(kb)
    kb_table = kb_table.table.copy(deep=True)

    b0_table = paramtable.ParamTable(ressep=ressep)
    b0_table.from_file(b0)
    b0_table = b0_table.table.copy(deep=True)

    idx = (kb_table == 0.0)
    maxkb = np.maximum(
        kb_table[kb_table.columns[1:].values],
        kb_table[kb_table.columns[:-1]].values
    )

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

    with openany("normed_kb.txt", "w") as output:
        kb_table.to_csv(
            output,
            header=True,
            index=True,
            sep=" ",
            float_format="%.4f",
        )

    with openany("normed_b0.txt", "w") as output:
        b0_table.to_csv(
            output,
            header=True,
            index=True,
            sep=" ",
            float_format="%.4f",
        )
