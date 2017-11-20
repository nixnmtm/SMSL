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

from fluctmatch.analysis import entropy


@click.command(
    "entropy",
    short_help="Calculate the Shannon entropy of residues."
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
    ent_table = entropy.Entropy(table, ressep=ressep)

    filename = path.join(outdir, "coupling.entropy.txt")
    with open(filename, mode="wb") as output:
        ent = ent_table.coupling_entropy().to_csv(
            index=True,
            header=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(ent.encode())

    filename = path.join(outdir, "relative.entropy.txt")
    with open(filename, mode="wb") as output:
        ent = ent_table.relative_entropy().to_csv(
            index=True,
            header=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(ent.encode())

    filename = path.join(outdir, "windiff.entropy.txt")
    with open(filename, mode="wb") as output:
        ent = ent_table.windiff_entropy().to_csv(
            index=True,
            header=True,
            sep=native_str(" "),
            float_format=native_str("%.4f"),
            encoding="utf-8",
        )
        output.write(ent.encode())
