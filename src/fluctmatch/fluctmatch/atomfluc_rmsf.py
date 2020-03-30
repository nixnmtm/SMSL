# -*- coding: utf-8 -*-
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

import pandas as pd
import numpy as np
import os
from os import path
import textwrap
import logging
import subprocess
from MDAnalysis.lib import util
from fluctmatch.fluctmatch.data import charmm_afnma, charmm_afqha

from future.utils import (
    raise_with_traceback,
)

logger = logging.getLogger(__name__)


class AtomicFluctuations(object):
    """
    Compute atomic fluctuations from normal modes of Quasi harmonic and Normal mode Analyses


    """
    def __init__(self, topology, trajectory, **kwargs):
        self.kwargs = kwargs
        self.outdir = self.kwargs.get("outdir", os.getcwd())
        self.prefix = self.kwargs.get("prefix", "fluctmatch")
        self.temperature = self.kwargs.get("temperature", 300.0)
        self.trajectory = trajectory

        self.outdir = self.outdir
        self.filenames = dict(
            qha_input=path.join(self.outdir, "afqha.inp"),
            nma_input=path.join(self.outdir, "afnma.inp"),

            qha_log=path.join(self.outdir, "afqha.log"),
            nma_log=path.join(self.outdir, "afnma.log"),

            fixed_prm=path.join(self.outdir, ".".join((self.prefix, "prm"))),
            psf_file=path.join(self.outdir, ".".join((self.prefix, "psf"))),
            xplor_psf_file=path.join(self.outdir, ".".join((self.prefix,
                                                            "xplor", "psf"))),
            crd_file=path.join(self.outdir, ".".join((self.prefix, "cor"))),
            topology_file=path.join(self.outdir, ".".join((self.prefix,
                                                           "rtf"))),

            fit_traj_file=path.join(self.outdir, "qha.fit.dcd"),
            qha_avg_cor=path.join(self.outdir, "qha.avg.cor"),
            qhacor_file=path.join(self.outdir, "qha_mode1.cor"),
            nma_nmv_file=path.join(self.outdir, "nma_nmv.crd"),
            quasi_nmv_file=path.join(self.outdir, "quasi_nmv.crd"),

            nma_crd=path.join(self.outdir, ".".join((self.prefix, "mini",
                                                     "cor"))),
            nma_vib=path.join(self.outdir, ".".join((self.prefix, "vib"))),
            nmacor_file=path.join(self.outdir, "nma_mode1.cor"),
            traj_file=self.trajectory,
        )

    def run_atomic_fluct(self, charmm_exec=None):

        # Find CHARMM executable
        charmm_exec = (os.environ.get("CHARMMEXEC", util.which("charmm"))
                       if charmm_exec is None else charmm_exec)
        if charmm_exec is None:
            logger.exception(
                "Please set CHARMMEXEC with the location of your CHARMM "
                "executable file or add the charmm path to your PATH "
                "environment.")
            raise_with_traceback(
                OSError(
                    "Please set CHARMMEXEC with the location of your CHARMM "
                    "executable file or add the charmm path to your PATH "
                    "environment."))

        # Write CHARMM input files and run atomic fluctuations
        if not path.exists(self.filenames["qha_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = (
                "dimension chsize 1000000" if version >= 36 else "")
            with open(self.filenames["qha_input"], mode="wb") as charmm_file:
                logger.info("Writing CHARMM QHA input file.")
                charmm_qha_inp = charmm_afqha.afqha.format(
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames)
                charmm_qha_inp = textwrap.dedent(charmm_qha_inp[1:])
                charmm_file.write(charmm_qha_inp.encode())
            logger.info("Running QHA Atomic Fluctuations.")
            with open(self.filenames["qha_log"], "w") as log_file:
                subprocess.check_call(
                    [charmm_exec, "-i", self.filenames["qha_input"]],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                )

        if not path.exists(self.filenames["nma_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = (
                "dimension chsize 1000000" if version >= 36 else "")
            with open(self.filenames["nma_input"], mode="wb") as charmm_file:
                logger.info("Writing NMA CHARMM input file.")
                charmm_nma_inp = charmm_afnma.afnma.format(
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames)
                charmm_nma_inp = textwrap.dedent(charmm_nma_inp[1:])
                charmm_file.write(charmm_nma_inp.encode())

            logger.info("Running NMA Atomic Fluctuations.")
            with open(self.filenames["nma_log"], "w") as log_file:
                subprocess.check_call(
                    [charmm_exec, "-i", self.filenames["nma_input"]],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                )

    @staticmethod
    def get_atom_fluc_rmsf(file):
        """
        Read fluctuation coordinates

        Return:
            RMSF of atomic fluctuations as pd.Series
        """
        fluc = pd.read_csv(file, skiprows=5, delim_whitespace=True, header=None,
                           index_col=0, usecols=[0, 1, 2, 3, 4, 5, 6, 7])
        fluc.columns = ["resid", "CGid", "I", "deltaX", "deltaY", "deltaZ", "segid"]
        fluc["posfluc"] = np.sqrt(fluc[["deltaX", "deltaY", "deltaZ"]].sum(axis=1))
        return fluc

    def get_rmsf_diff_qha_nma(self):
        """

        Returns
        -------
        'pd.Dataframe' with rmsf difference between QHA and NMA rmsf fluctuations

        """
        nmafluc = self.get_atom_fluc_rmsf(self.filenames["nmacor_file"]).set_index(["I"])["posfluc"]
        qhafluc = self.get_atom_fluc_rmsf(self.filenames["qhacor_file"]).set_index(["I"])["posfluc"]
        diff = qhafluc-nmafluc
        return diff
