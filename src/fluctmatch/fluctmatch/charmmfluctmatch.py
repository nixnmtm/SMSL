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
# Nixon Raj, Timothy Click, Haw Yang, Jhih-Wei Chua
# Mechanical couplings of protein backbone and side chains exhibit
# scale-free network properties and specific hotspots for function
# Computational and Structural Biotechnology Journal, Volume 19, 2021, Pages 5309-5320
# https://doi.org/10.1016/j.csbj.2021.09.004.
#

"""Fluctuation matching using CHARMM.

Notes
------
For CHARMM to work with the fluctuation matching code, it must be
recompiled with some modifications to the source code. `ATBMX`, `MAXATC`,
`MAXCB` (located in dimens.fcm [c35] or dimens_ltm.src [c39]) must
be increased. `ATBMX` determines the number of bonds allowed per
atom, `MAXATC` describes the maximum number of atom types, and `MAXCB`
determines the maximum number of bond parameters in the CHARMM parameter
file. Additionally, `CHSIZE` may need to be increased if using an earlier
version (< c36).

Modified by Nixon Raj:
---------------------
Fluctuation matching sensitivity increased by adding bondwise convergence
instead of overall Root MEan Squared Error. Please look into the the "relative difference" part of the code.
In brief,
1. each bond is check for convergence, if certain bonds are taking long time for convergence then
2. the percentage of the number of bonds yet to converged is taken into account.
3. If this percentage is less that 0.3% of the total number of bonds in the system, then
4. Relative difference is implemented:
    4.1. Relative difference between (Fluct_diff: difference of previous step and target) and (tol:tolerance)
5. If relative difference is very low (< 0.00005) then these bonds are skipped as it is almost converged.

NOTE: The bonds that causes these long time convergence issues are mostly found in the flexible regions of the protein

"""

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)
from future.builtins import (
    dict,
    open,
    range,
    super,
)
from future.utils import (
    PY2,
    native_str,
    raise_with_traceback,
)

import copy
import logging
import os
import subprocess
import textwrap
import time
from os import path

import numpy as np
import pandas as pd
from scipy import constants
import MDAnalysis as mda
from MDAnalysis.lib import util
from MDAnalysis.coordinates.core import reader
from fluctmatch.fluctmatch import base as fmbase
from fluctmatch.fluctmatch import utils as fmutils
from fluctmatch.fluctmatch.data import (
    charmm_init,
    charmm_nma,
    charmm_thermo,
)
from fluctmatch.intcor import utils as icutils
from fluctmatch.parameter import utils as prmutils

if PY2:
    FileNotFoundError = IOError

logger = logging.getLogger(__name__)


class CharmmFluctMatch(fmbase.FluctMatch):
    """Fluctuation matching using CHARMM."""
    bond_def = ["I", "J"]
    error_hdr = ["step", "Kb_rms", "fluct_rms", "b0_rms"]

    def __init__(self, *args, **kwargs):
        """Initialization of fluctuation matching using the CHARMM program.

        Parameters
        ----------
        topology : filename or Topology object
            A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file;
            used to define the list of atoms. If the file includes bond
            information, partial charges, atom masses, ... then these data will
            be available to MDAnalysis. A "structure" file (PSF, PDB or GRO, in
            the sense of a topology) is always required. Alternatively, an
            existing :class:`MDAnalysis.core.topology.Topology` instance may
            also be given.
        topology_format
            Provide the file format of the topology file; ``None`` guesses it
            from the file extension [``None``] Can also pass a subclass of
            :class:`MDAnalysis.topology.base.TopologyReaderBase` to define a
            custom reader to be used on the topology file.
        format
            Provide the file format of the coordinate or trajectory file;
            ``None`` guesses it from the file extension. Note that this keyword
            has no effect if a list of file names is supplied because the
            "chained" reader has to guess the file format for each individual
            list member. [``None``] Can also pass a subclass of
            :class:`MDAnalysis.coordinates.base.ProtoReader` to define a custom
            reader to be used on the trajectory file.
        guess_bonds : bool, optional
            Once Universe has been loaded, attempt to guess the connectivity
            between atoms.  This will populate the .bonds .angles and .dihedrals
            attributes of the Universe.
        vdwradii : dict, optional
            For use with *guess_bonds*. Supply a dict giving a vdwradii for each
            atom type which are used in guessing bonds.
        is_anchor : bool, optional
            When unpickling instances of
            :class:`MDAnalysis.core.groups.AtomGroup` existing Universes are
            searched for one where to anchor those atoms. Set to ``False`` to
            prevent this Universe from being considered. [``True``]
        anchor_name : str, optional
            Setting to other than ``None`` will cause
            :class:`MDAnalysis.core.groups.AtomGroup` instances pickled from
            the Universe to only unpickle if a compatible Universe with matching
            *anchor_name* is found. Even if *anchor_name* is set *is_anchor*
            will still be honored when unpickling.
        in_memory
            After reading in the trajectory, transfer it to an in-memory
            representations, which allow for manipulation of coordinates.
        in_memory_step
            Only read every nth frame into in-memory representation.
        outdir
            Output directory
        temperature
            Temperature (in K)
        rmin
            Minimum distance to consider for bond lengths.
        rmax
            Maximum distance to consider for bond lengths.
        charmm_version
            Version of CHARMM for formatting (default: 41)
        extended
            Use the extended format.
        title
            Title lines at the beginning of the file.
        resid
            Include segment IDs in the internal coordinate files.
        nonbonded
            Include the nonbonded section in the parameter file.
        """
        super().__init__(*args, **kwargs)
        self.dynamic_params = dict()
        self.filenames = dict(
            init_input=path.join(self.outdir, "fluctinit.inp"),
            init_log=path.join(self.outdir, "fluctinit.log"),
            init_avg_ic=path.join(self.outdir, "init.average.ic"),
            init_fluct_ic=path.join(self.outdir, "init.fluct.ic"),
            avg_ic=path.join(self.outdir, "average.ic"),
            fluct_ic=path.join(self.outdir, "fluct.ic"),
            dynamic_prm=path.join(self.outdir, "{}.dist.prm".format(
                self.prefix)),
            fixed_prm=path.join(self.outdir, ".".join((self.prefix, "prm"))),
            psf_file=path.join(self.outdir, ".".join((self.prefix, "psf"))),
            xplor_psf_file=path.join(self.outdir, ".".join((self.prefix,
                                                            "xplor", "psf"))),
            crd_file=path.join(self.outdir, ".".join((self.prefix, "cor"))),
            stream_file=path.join(self.outdir, ".".join((self.prefix,
                                                         "stream"))),
            topology_file=path.join(self.outdir, ".".join((self.prefix,
                                                           "rtf"))),
            nma_crd=path.join(self.outdir, ".".join((self.prefix, "mini",
                                                     "cor"))),
            nma_vib=path.join(self.outdir, ".".join((self.prefix, "vib"))),
            charmm_input=path.join(self.outdir, ".".join((self.prefix,
                                                          "inp"))),
            charmm_log=path.join(self.outdir, ".".join((self.prefix, "log"))),
            error_data=path.join(self.outdir, "error.dat"),
            thermo_input=path.join(self.outdir, "thermo.inp"),
            thermo_log=path.join(self.outdir, "thermo.log"),
            thermo_data=path.join(self.outdir, "thermo.dat"),
            traj_file=self.args[1] if len(self.args) > 1 else path.join(self.outdir, "cg.dcd"),
            bond_convergence=path.join(self.outdir, "bond_convergence.txt")
        )

        # Boltzmann constant
        self.BOLTZ = self.temperature * (constants.k * constants.N_A /
                                         (constants.calorie * constants.kilo))

        # Bond factor mol^2-Ang./kcal^2
        self.KFACTOR = 0.02

        # Self consistent error information.
        self.error = pd.DataFrame(
            np.zeros((1, len(self.error_hdr)), dtype=np.int),
            columns=self.error_hdr,
        )

    def _create_ic_table(self, universe, data):
        data.set_index(self.bond_def, inplace=True)
        table = icutils.create_empty_table(universe.atoms)
        hdr = table.columns
        table.set_index(self.bond_def, inplace=True)
        table.drop(
            [
                "r_IJ",
            ], axis=1, inplace=True)
        table = pd.concat([table, data["r_IJ"]], axis=1)
        return table.reset_index()[hdr]

    def initialize(self, nma_exec=None, restart=False):
        """Create an elastic network model from a basic coarse-grain model.

        Parameters
        ----------
        nma_exec : str
            executable file for normal mode analysis
        restart : bool, optional
            Reinitialize the object by reading files instead of doing initial
            calculations.
        """
        self.restart = restart
        if not self.restart:
            # Write CHARMM input file.
            if not path.exists(self.filenames["init_input"]):
                version = self.kwargs.get("charmm_version", 41)
                dimension = (
                    "dimension chsize 1000000" if version >= 36 else "")
                with open(
                    self.filenames["init_input"], mode="wb") as charmm_file:
                    logger.info("Writing CHARMM input file.")
                    charmm_inp = charmm_init.init.format(
                        flex="flex" if version else "",
                        version=version,
                        dimension=dimension,
                        **self.filenames)
                    charmm_inp = textwrap.dedent(charmm_inp[1:])
                    charmm_file.write(charmm_inp.encode())

            charmm_exec = (os.environ.get("CHARMMEXEC", util.which("charmm"))
                           if nma_exec is None else nma_exec)
            with open(self.filenames["init_log"], "w") as log_file:
                subprocess.check_call(
                    [charmm_exec, "-i", self.filenames["init_input"]],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                )

            # Write the parameter files.
            with reader(self.filenames["init_fluct_ic"]) as icfile:
                std_bonds = icfile.read().set_index(self.bond_def)
            with reader(self.filenames["init_avg_ic"]) as icfile:
                avg_bonds = icfile.read().set_index(self.bond_def)
            target = pd.concat([std_bonds["r_IJ"], avg_bonds["r_IJ"]], axis=1)
            target.reset_index(inplace=True)
            logger.info("Calculating the initial CHARMM parameters...")
            universe = mda.Universe(
                self.filenames["xplor_psf_file"], self.filenames["crd_file"]
            )
            self.target = prmutils.create_empty_parameters(universe, **self.kwargs)
            target.columns = self.target["BONDS"].columns
            self.target["BONDS"] = target.copy(deep=True)
            self.parameters = copy.deepcopy(self.target)
            self.parameters["BONDS"]["Kb"] = (
                self.BOLTZ / self.parameters["BONDS"]["Kb"].apply(np.square))
            self.dynamic_params = copy.deepcopy(self.parameters)
            with mda.Writer(self.filenames["fixed_prm"], **self.kwargs) as prm:
                logger.info("Writing {}...".format(
                    self.filenames["fixed_prm"]))
                prm.write(self.parameters)
            with mda.Writer(self.filenames["dynamic_prm"],
                            **self.kwargs) as prm:
                logger.info("Writing {}...".format(
                    self.filenames["dynamic_prm"]))
                prm.write(self.dynamic_params)
        else:
            print("FM Restarted")
            if not path.exists(self.filenames["fixed_prm"]):
                self.initialize(nma_exec, restart=False)
            try:
                # Read the parameter files.
                logger.info("Loading parameter and internal coordinate files.")
                with reader(self.filenames["fixed_prm"]) as fixed:
                    self.parameters.update(fixed.read())
                with reader(self.filenames["dynamic_prm"]) as dynamic:
                    self.dynamic_params.update(dynamic.read())

                # Read the initial internal coordinate files.
                with reader(self.filenames["init_avg_ic"]) as init_avg:
                    avg_table = init_avg.read().set_index(
                        self.bond_def)["r_IJ"]

                with reader(self.filenames["init_fluct_ic"]) as init_fluct:
                    fluct_table = (init_fluct.read().set_index(
                        self.bond_def)["r_IJ"])
                table = pd.concat([fluct_table, avg_table], axis=1)

                # Set the target fluctuation values.
                logger.info("Files loaded successfully...")
                self.target = copy.deepcopy(self.parameters)
                self.target["BONDS"].set_index(self.bond_def, inplace=True)
                cols = self.target["BONDS"].columns
                table.columns = cols
                self.target["BONDS"] = table.copy(deep=True).reset_index()

            except (FileNotFoundError, IOError):
                raise_with_traceback(
                    (IOError("Some files are missing. Unable to restart.")))

    def run(self, nma_exec=None, tol=1.e-3, n_cycles=300, low_bound=0.):
        """Perform a self-consistent fluctuation matching.

        Parameters
        ----------
        nma_exec : str
            executable file for normal mode analysis
        tol : float, optional
            fluct difference tolerance
        n_cycles : int, optional
            number of fluctuation matching cycles
        low_bound  : float, optional
            lowest Kb values to reduce noise
        """
        # Find CHARMM executable
        charmm_exec = (os.environ.get("CHARMMEXEC", util.which("charmm"))
                       if nma_exec is None else nma_exec)
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

        # Read the parameters
        if not self.parameters:
            try:
                self.initialize(nma_exec, restart=True)
            except IOError:
                raise_with_traceback(
                    (IOError("Some files are missing. Unable to restart.")))

        # Write CHARMM input file.
        if not path.exists(self.filenames["charmm_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = ("dimension chsize 1000000" if version >= 36 else "")
            with open(
                self.filenames["charmm_input"], mode="wb") as charmm_file:
                logger.info("Writing CHARMM input file.")
                charmm_inp = charmm_nma.nma.format(
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames)
                charmm_inp = textwrap.dedent(charmm_inp[1:])
                charmm_file.write(charmm_inp.encode())

        # Set the indices for the parameter tables.
        self.target["BONDS"].set_index(self.bond_def, inplace=True)
        bond_values = self.target["BONDS"].columns

        # Check for restart.
        try:
            if os.stat(self.filenames["error_data"]).st_size > 0:
                with open(self.filenames["error_data"], "rb") as data:
                    error_info = pd.read_csv(
                        data,
                        header=0,
                        skipinitialspace=True,
                        delim_whitespace=True)
                    if not error_info.empty:
                        self.error["step"] = error_info["step"].values[-1]
            else:
                raise FileNotFoundError
        except (FileNotFoundError, OSError):
            with open(self.filenames["error_data"], "wb") as data:
                np.savetxt(
                    data, [
                        self.error_hdr,
                    ],
                    fmt=native_str("%15s"),  # Nix
                    delimiter=native_str(""))
        self.error["step"] += 1

        # Initiate an all true index data, for preserving bond convergence
        if not self.restart:
            temp = ~self.target["BONDS"]["Kb"].isna()
            temp = temp.reset_index()
            self.converge_bnd_list = temp.iloc[:, 2]

        # Start self-consistent iteration for Fluctuation Matching
        # Run simulation
        logger.info(f"Starting fluctuation matching--{n_cycles} iterations to run")
        if low_bound != 0.:
            logger.info(f"Lower bound after 75% iteration is set to {low_bound}")
        st = time.time()
        fdiff = []
        for i in range(n_cycles):
            ct = time.time()
            self.error["step"] = i + 1
            with open(self.filenames["charmm_log"], "w") as log_file:
                subprocess.check_call(
                    [charmm_exec, "-i", self.filenames["charmm_input"]],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                )
            self.dynamic_params["BONDS"].set_index(self.bond_def, inplace=True)
            self.parameters["BONDS"].set_index(self.bond_def, inplace=True)

            # Read the average bond distance.
            with reader(self.filenames["avg_ic"]) as icavg:
                avg_ic = icavg.read().set_index(self.bond_def)["r_IJ"]

            # Read the bond fluctuations.
            with reader(self.filenames["fluct_ic"]) as icfluct:
                fluct_ic = icfluct.read().set_index(self.bond_def)["r_IJ"]

            vib_ic = pd.concat([fluct_ic, avg_ic], axis=1)
            vib_ic.columns = bond_values
            logger.info(f"Checking for bondlist convergence")
            fluct_diff = np.abs(vib_ic[bond_values[0]] - self.target["BONDS"][bond_values[0]])
            fdiff.append(fluct_diff)
            fluct_diff = fluct_diff.reset_index()
            tmp = self.parameters["BONDS"][bond_values[0]].reset_index()

            if not self.restart:
                self.converge_bnd_list &= ((fluct_diff.iloc[:, 2] > tol) & (tmp.iloc[:, 2] > 0))
            else:
                if i == 0:
                    self.converge_bnd_list = ((fluct_diff.iloc[:, 2] > tol) & (tmp.iloc[:, 2] > 0))
                else:
                    self.converge_bnd_list &= ((fluct_diff.iloc[:, 2] > tol) & (tmp.iloc[:, 2] > 0))

            # Calculate the r.m.s.d. between fluctuation and distances
            # compared with the target values.
            vib_error = self.target["BONDS"] - vib_ic
            vib_error = vib_error.apply(np.square).mean(axis=0)
            vib_error = np.sqrt(vib_error)
            self.error[self.error.columns[-2:]] = vib_error.T.values

            # Calculate the new force constant.
            optimized = vib_ic.apply(np.reciprocal).apply(np.square)
            target = self.target["BONDS"].apply(np.reciprocal).apply(np.square)
            optimized -= target
            optimized *= self.BOLTZ * self.KFACTOR

            # update  bond list
            vib_ic[bond_values[0]] = (self.parameters["BONDS"][bond_values[0]]
                                      - optimized[bond_values[0]])
            vib_ic[bond_values[0]] = (
                vib_ic[bond_values[0]].where(vib_ic[bond_values[0]] >= 0., 0.))  # set negative to zero

            if low_bound > 0. and i > int(n_cycles * 0.75):
                logger.info(f"Fluctuation matching cycle {i}: low bound is {low_bound}")
                vib_ic[bond_values[0]] = (vib_ic[bond_values[0]].where(vib_ic[bond_values[0]] >= low_bound, 0.))

            # r.m.s.d. between previous and current force constant
            diff = self.dynamic_params["BONDS"] - vib_ic
            diff = diff.apply(np.square).mean(axis=0)
            diff = np.sqrt(diff)
            self.error[self.error.columns[1]] = diff.values[0]

            # Update the parameters and write to file.
            self.parameters["BONDS"][bond_values[0]] = vib_ic[bond_values[0]]
            self.dynamic_params["BONDS"][bond_values[0]] = vib_ic[bond_values[0]]
            self.dynamic_params["BONDS"][bond_values[1]] = vib_ic[bond_values[1]]

            self.parameters["BONDS"].reset_index(inplace=True)
            self.dynamic_params["BONDS"].reset_index(inplace=True)
            with mda.Writer(self.filenames["fixed_prm"], **self.kwargs) as prm:
                prm.write(self.parameters)
            with mda.Writer(self.filenames["dynamic_prm"],
                            **self.kwargs) as prm:
                prm.write(self.dynamic_params)

            # Update the error values.
            with open(self.filenames["error_data"], "ab") as error_file:
                np.savetxt(
                    error_file,
                    self.error,
                    fmt=native_str("%15d%15.6f%15.6f%15.6f", ),  # Nix
                    delimiter=native_str(""),
                )
            logger.info("Fluctuation matching cycle {} completed in {:.6f}".format(
                i, time.time() - ct))
            logger.info(f"{self.converge_bnd_list.sum()} not converged out of {len(self.converge_bnd_list)}")

            if self.converge_bnd_list.sum() <= len(self.converge_bnd_list.values.tolist()) * 0.003:
                # if bonds to converge is less than 0.3% of total bonds, use relative difference as criteria
                # as it takes more than 100 iterations for these 0.3%  bonds to converge.
                relative_diff = (fluct_diff.iloc[:, 2] - tol) / tol

                ### To know the late converged bonds uncomment the below 5 lines ###

                # late_converged = pd.DataFrame()
                # indx = self.converge_bnd_list[self.converge_bnd_list].index.values
                # late_converged = pd.concat([fluct_diff.loc[indx], relative_diff.loc[indx]], axis=1)
                # late_converged.columns = ["I", "J", "fluct_diff_Kb", "relative_diff_kb"]
                # print(late_converged)

                self.converge_bnd_list = self.converge_bnd_list & (relative_diff > 5)
                if self.converge_bnd_list.sum() == 0:
                    logger.info("Checking relative difference: All bonds converged, exiting")
                    break
        fluct_conv = pd.concat(fdiff, axis=1).round(6)
        fluct_conv.columns = [j for j in range(1, i + 2)]
        fluct_conv.to_csv(self.filenames["bond_convergence"])
        logger.info("Fluctuation matching completed in {:.6f}".format(
            time.time() - st))
        self.target["BONDS"].reset_index(inplace=True)

    def calculate_thermo(self, nma_exec=None):
        """Calculate the thermodynamic properties of the trajectory.

        Parameters
        ----------
        nma_exec : str
            executable file for normal mode analysis
        """
        # Find CHARMM executable
        charmm_exec = (os.environ.get("CHARMMEXEC", util.which("charmm"))
                       if nma_exec is None else nma_exec)
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

        if not path.exists(self.filenames["thermo_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = ("dimension chsize 500000 maxres 3000000"
                         if version >= 36 else "")
            with open(
                self.filenames["thermo_input"], mode="wb") as charmm_file:
                logger.info("Writing CHARMM input file.")
                charmm_inp = charmm_thermo.thermodynamics.format(
                    trajectory=path.join(self.outdir, self.args[-1]),
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames)
                charmm_inp = textwrap.dedent(charmm_inp[1:])
                charmm_file.write(charmm_inp.encode())

        # Calculate thermodynamic properties of the trajectory.
        with open(self.filenames["thermo_log"], "w") as log_file:
            logger.info("Running thermodynamic calculation.")
            subprocess.check_call(
                [charmm_exec, "-i", self.filenames["thermo_input"]],
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )
            logger.info("Calculations completed.")

        header = ("SEGI  RESN  RESI     Entropy    Enthalpy     "
                  "Heatcap     Atm/res   Ign.frq")
        columns = np.array(header.split())
        columns[:3] = np.array(["segidI", "RESN", "resI"])
        thermo = []

        # Read log file
        with open(self.filenames["thermo_log"], "rb") as log_file:
            logger.info("Reading CHARMM log file.")
            for line in log_file:
                if line.find(header) < 0:
                    continue
                break
            for line in log_file:
                if len(line.strip().split()) == 0:
                    break
                thermo.append(line.strip().split())

        # Create human-readable table
        thermo = pd.DataFrame(thermo, columns=columns)
        thermo.drop(["RESN", "Atm/res", "Ign.frq"], axis=1, inplace=True)
        thermo.set_index(["segidI", "resI"], inplace=True)
        thermo = thermo.astype(np.float)

        # Write data to file
        with open(self.filenames["thermo_data"], "wb") as data_file:
            logger.info("Writing thermodynamics data file.")
            thermo = thermo.to_csv(
                index=True,
                sep=native_str(" "),
                float_format=native_str("%.4f"),
                encoding="utf-8")
            data_file.write(thermo.encode())
