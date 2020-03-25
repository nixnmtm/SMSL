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

import logging
import logging.config
import os
from os import path

import click
from MDAnalysis.lib.util import which
from fluctmatch.fluctmatch import tune_topology
from fluctmatch.fluctmatch import atomfluc_rmsf
from fluctmatch.fluctmatch import charmmfluctmatch


@click.command("tune_fluc", short_help="Run Atomic Fluc")
@click.option(
    "--rcut",
    "--rcut",
    metavar="RMSF CUT",
    type=click.FLOAT,
    default=0.,
    show_default=True,
    help="RMSF Difference (abs(afQHA-afNMA)) cutoff",
)
@click.option(
    "--dcut",
    "--dcut",
    metavar="DISTANCE CUT",
    type=click.FLOAT,
    default=6.5,
    show_default=True,
    help="Trimming the topology to distance cut off, "
         "usually atleast 1 Ang less than the overall rmax",
)
@click.option(
    "--kbcut",
    "--kbcut",
    metavar="Force Constant CUT",
    type=click.FLOAT,
    default=0.,
    show_default=True,
    help="Trimming the bonds only if kb greater than kbcut",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    show_default=True,
    default=path.join(os.getcwd(), "charmmfm.log"),
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="Log file",
)
@click.option(
    "-o",
    "outdir",
    metavar="DIR",
    default=os.getcwd(),
    show_default=True,
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help="Directory",
)
@click.option(
    "-e",
    "--exec",
    "nma_exec",
    metavar="FILE",
    envvar="CHARMMEXEC",
    default=which("charmm"),
    show_default=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help="CHARMM executable file",
)
@click.option(
    "-t",
    "--temperature",
    metavar="TEMP",
    type=click.FLOAT,
    default=300.0,
    show_default=True,
    help="Temperature of simulation",
)

@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default="fluctmatch",
    show_default=True,
    type=click.STRING,
    help="Prefix for filenames",
)
@click.option(
    "-c",
    "--charmm",
    "charmm_version",
    metavar="VERSION",
    default=41,
    show_default=True,
    type=click.IntRange(27, None, clamp=True),
    help="CHARMM version",
)
@click.option(
    "--extended / --standard",
    "extended",
    default=True,
    help="Output using the extended or standard columns",
)
@click.option(
    "--nb / --no-nb",
    "nonbonded",
    default=True,
    help="Include nonbonded section in CHARMM parameter file",
)
@click.option(
    "--resid / --no-resid",
    "resid",
    default=True,
    help="Include segment IDs in internal coordinate files",
)
@click.option(
    "--restart",
    is_flag=True,
    help="Restart Fluctuation Matching using the trimmed topology",
)
@click.option(
    "-n",
    "--ncycles",
    "n_cycles",
    metavar="NCYCLES",
    type=click.IntRange(1, None, clamp=True),
    default=300,
    show_default=True,
    help="Number of simulation cycles",
)
@click.option(
    "-m",
    "--minval",
    "low_bound",
    metavar="LOW_BOUND",
    type=click.FLOAT,
    default=0.,
    show_default=True,
    help="Lower bound Kb value. Less than this values is set 0"
)
@click.option(
    "--tol",
    metavar="TOL",
    type=click.FLOAT,
    default=1.e-3,
    show_default=True,
    help="Fluct difference tolerance between simulations",
)
def cli(rcut,
        dcut,
        kbcut,
        logfile,
        outdir,
        nma_exec,
        temperature,
        prefix,
        charmm_version,
        extended,
        resid,
        nonbonded,
        restart,
        n_cycles,
        tol,
        low_bound,
):
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

    kwargs = dict(

        prefix=prefix,
        outdir=outdir,
        temperature=temperature,
        charmm_version=charmm_version,
        extended=extended,
        resid=resid,
        nonbonded=nonbonded,

    )

    logger.info("Tuning FM topology Started")
    # Run atomic fluctuations of NMA and QHA using CHARMM
    af = atomfluc_rmsf.AtomicFluctuations(**kwargs)
    af.run_atomic_fluct(charmm_exec=nma_exec)
    logger.info("Computing fluctuation difference")
    fluct_diff = af.get_rmsf_diff_qha_nma()

    kwargs.update(dict(
        rcut=rcut,
        dcut=dcut,
        kbcut=kbcut)
    )

    # Tuning Topology based on distance cutoff
    logger.info("Trimming bondlists based on dcut")
    tunefm = tune_topology.TopologyTuning(**kwargs)
    tunefm.tune_topology(fluct_diff)

    # Run FM in trimmed folder

    writedir = os.path.join(outdir, "trimmed")

    kwargs.update(
        dict(restart=restart,
             tol=tol,
             n_cycles=n_cycles,
             low_bound=low_bound,
             outdir=writedir)
    )

    #  Run Fluctuation Matching for the trimmed topolgy untill bonds converge
    cfm = charmmfluctmatch.CharmmFluctMatch(**kwargs)
    logger.info("Initializing the parameters.")
    cfm.initialize(nma_exec=nma_exec, restart=restart)
    logger.info("Running fluctuation matching with new trimmed topology")
    cfm.run(nma_exec=nma_exec, tol=tol, n_cycles=n_cycles, low_bound=low_bound)
    logger.info("Fluctuation matching successfully completed.")
