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
from future.utils import native_str, PY3
from future.builtins import open
from six.moves import cPickle

import logging
import logging.config
from os import path

import click
import numpy as np
import pandas as pd
from fluctmatch.analysis.paramtable import ParamTable
from fluctmatch.analysis import (
    fluctsca, )


@click.command(
    "sca",
    short_help="Statistical coupling analysis (SCA) on coupling strength")
@click.option(
    "-n",
    "--ntrials",
    metavar="NTRIALS",
    default=100,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Number of random iterations")
@click.option(
    "--std",
    metavar="STDDEV",
    default=2,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Number of std. deviations for beyond second eigenmode")
@click.option(
    "-k",
    "--kpos",
    metavar="KPOS",
    default=0,
    type=click.IntRange(0, None, clamp=True),
    help="Number of eigenmodes [default: auto]")
@click.option(
    "-p",
    "--pcut",
    default=0.95,
    show_default=True,
    type=np.float,
    help="Cutoff value for sector selection")
@click.option(
    "-r",
    "--ressep",
    metavar="RESSEP",
    default=3,
    show_default=True,
    type=click.IntRange(0, None, clamp=True),
    help="Number of residues to exclude in I,I+r")
@click.option(
    "-o",
    "--output",
    default=path.join(path.curdir, "scafluct.db"),
    show_default=True,
    type=click.Path(
        exists=False,
        file_okay=True,
        resolve_path=True,
    ),
    help="Output filename")
@click.option(
    "-s",
    "--subset",
    metavar="SEGID RES RES",
    type=(native_str, click.IntRange(1, None, clamp=True),
          click.IntRange(1, None, clamp=True)),
    multiple=True,
    help="Subset of a system (SEGID FIRST LAST)")
@click.option(
    "--all",
    "transformation",
    flag_value="all",
    default=True,
)
@click.option(
    "--bb",
    "transformation",
    flag_value="backbone",
)
@click.option(
    "--sc",
    "transformation",
    flag_value="sidechain",
)
@click.option(
    "--bbsc",
    "transformation",
    flag_value="bbsc",
)
@click.argument(
    "filename",
    type=click.Path(
        exists=True,
        file_okay=True,
        resolve_path=True,
    ))
def cli(ntrials, std, kpos, pcut, ressep, output, subset, transformation,
        filename):
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
                "filename" : path.join(path.dirname(filename), "fluctsca.log"),
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

    # Load the table, separate by I,I+r, and if requested, create a subset.
    logger.info("Loading parameter table {}".format(
        click.format_filename(filename)))
    table = ParamTable(ressep=ressep)
    table.from_file(click.format_filename(filename))
    kb = table._separate(table.table)
    _index = kb.index.names
    if transformation == "backbone":
        logger.debug("Calculating backbone-backbone interactions.")
        kb.reset_index(inplace=True)
        kb = kb[(kb["I"] == "N") | (kb["I"] == "CA") | (kb["I"] == "O")]
        kb = kb[(kb["J"] == "N") | (kb["J"] == "CA") | (kb["J"] == "O")]
        table.table = kb.set_index(_index)
    elif transformation == "sidechain":
        logger.debug("Calculating sidechain-sidechain interactions.")
        kb.reset_index(inplace=True)
        kb = kb[(kb["I"] == "CB") & (kb["J"] == "CB")]
        table.table = kb.set_index(_index)
    elif transformation == "bbsc":
        logger.debug("Calculating backbone-sidechain interactions.")
        kb.reset_index(inplace=True)
        tmp1 = kb[(kb["I"] == "CB") & ((kb["J"] == "N") | (kb["J"] == "O"))]
        tmp2 = kb[(kb["J"] == "CB") & ((kb["I"] == "N") | (kb["I"] == "O"))]
        table.table = pd.concat([tmp1, tmp2], axis=0).set_index(_index)
    else:
        logger.debug("Accounting for all interactions.")
    kb = table.per_residue

    D_info = dict(
        kb=kb,
        ressep=ressep,
    )

    if subset:
        segid, start, stop = subset[0]
        logger.info("Using a subset of {} between {:d} and {:d}".format(
            segid, start, stop))
        kb = kb.loc[segid].loc[start:stop]
        D_info["kb"] = kb.copy(deep=True)
        D_info["subset"] = subset[0]

    # Calculate eigenvalues and eigenvectors for the time series with sign correction.
    U, Lsca, Vt = fluctsca.svd(kb)

    logger.info("Using {:d} random trials.".format(ntrials))
    Lrand = fluctsca.randomize(kb, ntrials=ntrials)
    Ucorrel = kb.values.dot(kb.T.values)
    Vcorrel = kb.values.T.dot(kb.values)
    D_sca = dict(
        U=U,
        Lsca=Lsca,
        Vt=Vt,
        Lrand=Lrand,
        Ucorrel=Ucorrel,
        Vcorrel=Vcorrel,
        ntrials=ntrials)

    # Determine the number of eigenmodes if kpos = 0
    _kpos = fluctsca.chooseKpos(Lsca, Lrand, stddev=std) if kpos == 0 else kpos
    logger.info("Selecting {:d} eigenmodes".format(_kpos))
    Ucorr = fluctsca.correlate(U, Lsca, kmax=_kpos)
    Vcorr = fluctsca.correlate(Vt.T, Lsca, kmax=_kpos)

    # Calculate IC sectors
    logger.debug("Calculating the ICA for the residues.")
    Uica, Wres = fluctsca.rotICA(U, kmax=_kpos)
    Uics, Uicsize, Usortedpos, Ucutoff, Uscaled_pd, Upd = fluctsca.icList(
        Uica, _kpos, Ucorrel, p_cut=pcut)

    logger.debug("Calculating the ICA for the windows.")
    Vica, Wtime = fluctsca.rotICA(Vt.T, kmax=_kpos)
    Utica = Wtime.dot(U[:, :_kpos].T).T
    Vrica = Wres.dot(Vt[:_kpos]).T
    Vics, Vicsize, Vsortedpos, Vcutoff, Vscaled_pd, Vpd = fluctsca.icList(
        Vica, _kpos, Vcorrel, p_cut=pcut)
    D_sector = dict(
        std=std,
        kpos=_kpos,
        Ucorr=Ucorr,
        Vcorr=Vcorr,
        Uica=Uica,
        Wres=Wres,
        Vica=Vica,
        Wtime=Wtime,
        Uics=Uics,
        Uicsize=Uicsize,
        Usortedpos=Usortedpos,
        Ucutoff=Ucutoff,
        Uscaled_pd=Uscaled_pd,
        Upd=Upd,
        Vics=Vics,
        Vicsize=Vicsize,
        Vsortedpos=Vsortedpos,
        Vcutoff=Vcutoff,
        Vscaled_pd=Vscaled_pd,
        Vpd=Vpd,
        Utica=Utica,
        Vrica=Vrica)

    D = dict(info=D_info, sca=D_sca, sector=D_sector)
    with open(click.format_filename(output), mode="wb") as dbf:
        logger.info("Saving data to {}".format(click.format_filename(output)))
        if PY3:
            logger.warning(
                "Note: The saved file will be incompatible with Python 2.")
        cPickle.dump(D, dbf, protocol=cPickle.HIGHEST_PROTOCOL)
