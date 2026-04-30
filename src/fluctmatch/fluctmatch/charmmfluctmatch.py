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
-----
This workflow relies on a correctly configured CHARMM build capable of
handling large coarse-grained fluctuation-matching calculations.
Use `charmmsize_increase_before_configuration.py` together with this
workflow to modify the CHARMM source code before configuration when
larger problem sizes are required. Please ensure the necessary source
changes are applied and that CHARMM is recompiled before use.

Workflow
--------
The initialization step computes target coarse-grained average bond
distances and fluctuations from the mapped trajectory and writes them to
`init.average.ic` and `init.fluct.ic`. During each fluctuation-matching
cycle, CHARMM normal mode analysis evaluates the current elastic network
model and produces `average.ic` and `fluct.ic`. The model parameters are
then iteratively updated so that the ENM-derived fluctuations converge
toward the target coarse-grained fluctuations derived from the trajectory.

Modified by Nixon Raj
---------------------
This implementation extends the original fluctuation-matching procedure
with bondwise convergence tracking, tail-convergence filtering, robust
fresh/restart/resume/auto execution modes, conservative checkpoint-based
recovery, and per-cycle timestamped progress logging.

Briefly,
1. Convergence is monitored bondwise rather than relying only on the
   overall root-mean-squared error.
2. At each cycle, every bond is evaluated individually against the
   fluctuation-difference tolerance.
3. The number and fraction of bonds that remain unconverged are tracked
   throughout the run using an active unconverged-bond mask.
4. When the remaining unconverged bonds fall below 0.3% of the total
   number of bonds in the system, a tail-convergence rule is applied.
5. In this tail region, a relative difference is evaluated as
   `(fluct_diff - tol) / tol`, where `fluct_diff` is the bondwise
   difference between the current fluctuation and the target fluctuation,
   and `tol` is the convergence tolerance.
6. Bonds with only a very small relative excess above tolerance are
   treated as effectively converged and removed from the active
   unconverged set.
7. A fresh run initializes the fluctuation-matching state from the
   trajectory-derived `init.average.ic` and `init.fluct.ic` files and
   begins the cycle count from step 1.
8. Restart explicitly rebuilds the fluctuation-matching state from the
   saved initialization files and begins the cycle count again from
   step 1.
9. Resume explicitly reloads the latest checkpointed run state,
   including the current parameter files, bond mask, and error history,
   and continues from the most conservative valid completed step.
10. Auto mode selects the execution path per window by prioritizing
    resume when a valid checkpointed state is available, otherwise
    falling back to restart if initialization files are present, and
    otherwise starting fresh.
11. Resume recovery is tolerant to minor checkpoint and `error.dat`
    step mismatches by continuing from the minimum consistent completed
    step rather than discarding recoverable progress.
12. Each completed fluctuation-matching cycle is logged with detailed
    timing information together with a wall-clock timestamp, making
    runtime progress easier to monitor on long HPC jobs.

Note
----
The bonds that remain longest in the active unconverged set are often
associated with highly flexible regions of the protein.

"""

import copy
import json
import logging
import os
import subprocess
import textwrap
import time
import uuid
from os import path

import numpy as np
import pandas as pd
from scipy import constants
import MDAnalysis as mda
from MDAnalysis.lib import util
from MDAnalysis.coordinates.core import reader
from fluctmatch.fluctmatch import base as fmbase
from fluctmatch.fluctmatch.data import (
    charmm_init,
    charmm_nma,
    charmm_thermo,
)
from fluctmatch.intcor import utils as icutils
from fluctmatch.parameter import utils as prmutils

from datetime import datetime, timedelta

logger = logging.getLogger(__name__)

def _fmt_timing(parts):
    return " | ".join(f"{k}={v:.3f}s" for k, v in parts.items())


class CharmmFluctMatch(fmbase.FluctMatch):
    """Fluctuation matching using CHARMM."""
    bond_def = ["I", "J"]
    error_hdr = ["step", "Kb_rms", "fluct_rms", "b0_rms"]
    TAIL_FRACTION = 0.003
    TAIL_RELATIVE_DIFF_THRESHOLD = 1e-3

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dynamic_params = dict()
        self.filenames = dict(
            init_input=path.join(self.outdir, "fluctinit.inp"),
            init_log=path.join(self.outdir, "fluctinit.log"),
            init_avg_ic=path.join(self.outdir, "init.average.ic"),
            init_fluct_ic=path.join(self.outdir, "init.fluct.ic"),
            avg_ic=path.join(self.outdir, "average.ic"),
            fluct_ic=path.join(self.outdir, "fluct.ic"),
            dynamic_prm=path.join(self.outdir, "{}.dist.prm".format(self.prefix)),
            fixed_prm=path.join(self.outdir, ".".join((self.prefix, "prm"))),
            psf_file=path.join(self.outdir, ".".join((self.prefix, "psf"))),
            xplor_psf_file=path.join(self.outdir, ".".join((self.prefix, "xplor", "psf"))),
            crd_file=path.join(self.outdir, ".".join((self.prefix, "cor"))),
            stream_file=path.join(self.outdir, ".".join((self.prefix, "stream"))),
            topology_file=path.join(self.outdir, ".".join((self.prefix, "rtf"))),
            nma_crd=path.join(self.outdir, ".".join((self.prefix, "mini", "cor"))),
            nma_vib=path.join(self.outdir, ".".join((self.prefix, "vib"))),
            charmm_input=path.join(self.outdir, ".".join((self.prefix, "inp"))),
            charmm_log=path.join(self.outdir, ".".join((self.prefix, "log"))),
            error_data=path.join(self.outdir, "error.dat"),
            thermo_input=path.join(self.outdir, "thermo.inp"),
            thermo_log=path.join(self.outdir, "thermo.log"),
            thermo_data=path.join(self.outdir, "thermo.dat"),
            traj_file=self.args[1] if len(self.args) > 1 else path.join(self.outdir, "cg.dcd"),
            bond_convergence=path.join(self.outdir, "bond_convergence.txt"),
            bond_mask=path.join(self.outdir, "bond_mask.txt"),
            checkpoint=path.join(self.outdir, "fm_checkpoint.json"),
        )

        self.BOLTZ = self.temperature * (
            constants.k * constants.N_A / (constants.calorie * constants.kilo)
        )
        self.KFACTOR = 0.02
        self.error = pd.DataFrame(
            np.zeros((1, len(self.error_hdr)), dtype=int),
            columns=self.error_hdr,
        )

        self.auto = False
        self.restart = False
        self.resume = False

    def _create_ic_table(self, universe, data):
        data.set_index(self.bond_def, inplace=True)
        table = icutils.create_empty_table(universe.atoms)
        hdr = table.columns
        table.set_index(self.bond_def, inplace=True)
        table.drop(["r_IJ"], axis=1, inplace=True)
        table = pd.concat([table, data["r_IJ"]], axis=1)
        return table.reset_index()[hdr]

    def _run_charmm(self, charmm_exec):
        with open(self.filenames["charmm_log"], "w") as log_file:
            subprocess.check_call(
                [charmm_exec, "-i", self.filenames["charmm_input"]],
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )

    def _read_ic_series(self, filename):
        with reader(filename) as icf:
            table = icf.read()
        return table.set_index(self.bond_def)["r_IJ"]

    def _read_cycle_ic(self):
        avg_ic = self._read_ic_series(self.filenames["avg_ic"])
        fluct_ic = self._read_ic_series(self.filenames["fluct_ic"])
        return avg_ic, fluct_ic

    def _build_vib_ic(self, fluct_ic, avg_ic, bond_values):
        vib_ic = pd.concat([fluct_ic, avg_ic], axis=1)
        vib_ic.columns = bond_values
        return vib_ic

    def _bond_index(self):
        return self.target["BONDS"].set_index(self.bond_def).index

    def _initialize_convergence_mask(self):
        kb = self.target["BONDS"].set_index(self.bond_def)["Kb"]
        self.unconverged_bond_list = (~kb.isna()).astype(bool).copy()

    def _update_convergence_mask(self, fluct_diff_series, current_kb_series, tol, i):
        mask = ((fluct_diff_series > tol) & (current_kb_series > 0)).astype(bool)
        if not self.unconverged_bond_list.index.equals(mask.index):
            self.unconverged_bond_list = self.unconverged_bond_list.reindex(mask.index)
            self.unconverged_bond_list = self.unconverged_bond_list.fillna(True).astype(bool)
        if not self.restart:
            self.unconverged_bond_list &= mask
        else:
            if i == 0:
                self.unconverged_bond_list = mask.copy()
            else:
                self.unconverged_bond_list &= mask

    def _compute_convergence_and_error(self, vib_ic, bond_values, tol, i):
        fluct_diff_series = np.abs(
            vib_ic[bond_values[0]] - self.target["BONDS"][bond_values[0]]
        )
        current_kb_series = self.parameters["BONDS"][bond_values[0]]
        self._update_convergence_mask(fluct_diff_series, current_kb_series, tol, i)

        vib_error = self.target["BONDS"] - vib_ic
        vib_error = vib_error.apply(np.square).mean(axis=0)
        vib_error = np.sqrt(vib_error)
        self.error[self.error.columns[-2:]] = vib_error.T.values
        return fluct_diff_series

    def _optimize_bonds(self, vib_ic, bond_values, low_bound, n_cycles, i):
        optimized = vib_ic.apply(np.reciprocal).apply(np.square)
        target = self.target["BONDS"].apply(np.reciprocal).apply(np.square)
        optimized -= target
        optimized *= self.BOLTZ * self.KFACTOR

        vib_ic[bond_values[0]] = (
            self.parameters["BONDS"][bond_values[0]] - optimized[bond_values[0]]
        )
        vib_ic[bond_values[0]] = vib_ic[bond_values[0]].where(
            vib_ic[bond_values[0]] >= 0.0, 0.0
        )

        if low_bound > 0.0 and i > int(n_cycles * 0.75):
            logger.info("Fluctuation matching cycle %d: low bound is %s", i, low_bound)
            vib_ic[bond_values[0]] = vib_ic[bond_values[0]].where(
                vib_ic[bond_values[0]] >= low_bound, 0.0
            )

        diff = self.dynamic_params["BONDS"] - vib_ic
        diff = diff.apply(np.square).mean(axis=0)
        diff = np.sqrt(diff)
        self.error[self.error.columns[1]] = diff.values[0]
        return vib_ic

    def _apply_updated_parameters(self, vib_ic, bond_values):
        self.parameters["BONDS"][bond_values[0]] = vib_ic[bond_values[0]]
        self.dynamic_params["BONDS"][bond_values[0]] = vib_ic[bond_values[0]]
        self.dynamic_params["BONDS"][bond_values[1]] = vib_ic[bond_values[1]]
        self.parameters["BONDS"].reset_index(inplace=True)
        self.dynamic_params["BONDS"].reset_index(inplace=True)

    def _write_parameter_files(self):
        with mda.Writer(self.filenames["fixed_prm"], **self.kwargs) as prm:
            prm.write(self.parameters)
        with mda.Writer(self.filenames["dynamic_prm"], **self.kwargs) as prm:
            prm.write(self.dynamic_params)

    def _append_error_row(self):
        with open(self.filenames["error_data"], "ab") as error_file:
            np.savetxt(
                error_file,
                self.error,
                fmt="%15d%15.6f%15.6f%15.6f",
                delimiter="",
            )

    def _write_error_header(self):
        with open(self.filenames["error_data"], "wb") as data:
            np.savetxt(
                data,
                [self.error_hdr],
                fmt="%15s",
                delimiter="",
            )

    def _get_error_last_step(self):
        try:
            if os.stat(self.filenames["error_data"]).st_size > 0:
                with open(self.filenames["error_data"], "rb") as data:
                    error_info = pd.read_csv(
                        data,
                        header=0,
                        skipinitialspace=True,
                        sep=r"\s+",
                        engine="python",
                    )
                if not error_info.empty:
                    return int(error_info["step"].values[-1])
        except (FileNotFoundError, OSError, pd.errors.EmptyDataError, KeyError, ValueError):
            pass
        return 0

    def _write_json(self, filename, payload):
        with open(filename, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)

    def _read_json(self, filename):
        with open(filename, "r") as fh:
            return json.load(fh)

    def _save_bond_mask(self, filename=None):
        bond_mask_file = filename or self.filenames["bond_mask"]
        pd.Series(self.unconverged_bond_list.astype(bool)).to_csv(
            bond_mask_file, index=False, header=["active"]
        )

    def _load_bond_mask(self, filename):
        mask_df = pd.read_csv(filename)
        if "active" in mask_df.columns:
            mask = mask_df["active"]
        else:
            mask = mask_df.iloc[:, 0]

        mask = mask.astype(bool).reset_index(drop=True)
        bond_index = self._bond_index()
        expected = len(bond_index)
        if mask.size != expected:
            raise IOError(
                "Bond mask length mismatch: expected {} got {}".format(
                    expected, mask.size
                )
            )
        mask.index = bond_index
        self.unconverged_bond_list = mask

    def _build_initial_state_from_init_ic(self):
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
            self.BOLTZ / self.parameters["BONDS"]["Kb"].apply(np.square)
        )
        self.dynamic_params = copy.deepcopy(self.parameters)
        self._initialize_convergence_mask()

    def _load_target_from_init_files(self):
        with reader(self.filenames["init_avg_ic"]) as init_avg:
            avg_table = init_avg.read().set_index(self.bond_def)["r_IJ"]
        with reader(self.filenames["init_fluct_ic"]) as init_fluct:
            fluct_table = init_fluct.read().set_index(self.bond_def)["r_IJ"]
        table = pd.concat([fluct_table, avg_table], axis=1)

        self.target = copy.deepcopy(self.parameters)
        self.target["BONDS"].set_index(self.bond_def, inplace=True)
        cols = self.target["BONDS"].columns
        table.columns = cols
        self.target["BONDS"] = table.copy(deep=True).reset_index()

    def _load_parameter_files(self, fixed_prm, dynamic_prm):
        logger.info("Loading parameter files: %s | %s", fixed_prm, dynamic_prm)
        with reader(fixed_prm) as fixed:
            self.parameters.update(fixed.read())
        with reader(dynamic_prm) as dynamic:
            self.dynamic_params.update(dynamic.read())

    def _load_restart_state(self):
        for required in (
            self.filenames["init_avg_ic"],
            self.filenames["init_fluct_ic"],
            self.filenames["xplor_psf_file"],
            self.filenames["crd_file"],
        ):
            if not path.exists(required):
                raise IOError("Required restart file is missing: {}".format(required))

        self.run_uuid = str(uuid.uuid4()) # since we are restrating, new runid created
        self._build_initial_state_from_init_ic()
        logger.info(
            "Restart state rebuilt successfully from init.average.ic and init.fluct.ic."
        )

    def _write_cycle_checkpoint(self, cycle_idx, tol, n_cycles, low_bound):
        self._save_bond_mask(self.filenames["bond_mask"])
        payload = {
            "run_uuid": getattr(self, "run_uuid", None),
            "last_completed_step": int(cycle_idx),
            "n_bonds": int(self.target["BONDS"].shape[0]),
            "tol": float(tol),
            "n_cycles_requested": int(n_cycles),
            "low_bound": float(low_bound),
            "files": {
                "fixed_prm": path.basename(self.filenames["fixed_prm"]),
                "dynamic_prm": path.basename(self.filenames["dynamic_prm"]),
                "bond_mask": path.basename(self.filenames["bond_mask"]),
                "error_data": path.basename(self.filenames["error_data"]),
                "init_avg_ic": path.basename(self.filenames["init_avg_ic"]),
                "init_fluct_ic": path.basename(self.filenames["init_fluct_ic"]),
            },
        }
        self._write_json(self.filenames["checkpoint"], payload)

    def _load_resume_state(self):
        if not path.exists(self.filenames["checkpoint"]):
            raise IOError("Checkpoint file is missing. Unable to resume.")

        state = self._read_json(self.filenames["checkpoint"])
        for required in (
            self.filenames["fixed_prm"],
            self.filenames["dynamic_prm"],
            self.filenames["bond_mask"],
            self.filenames["error_data"],
            self.filenames["init_avg_ic"],
            self.filenames["init_fluct_ic"],
        ):
            if not path.exists(required):
                raise IOError("Required resume state file is missing: {}".format(required))

        checkpoint_step = int(state["last_completed_step"])
        error_last_step = self._get_error_last_step()

        if checkpoint_step <= 0 or error_last_step <= 0:
            raise IOError(
                "Resume validation failed: invalid step values checkpoint={} error.dat={}".format(
                    checkpoint_step, error_last_step
                )
            )

        resume_step = min(checkpoint_step, error_last_step)

        if checkpoint_step != error_last_step:
            logger.warning(
                "Resume step mismatch detected: checkpoint step=%d, error.dat last step=%d. "
                "Using conservative resume_step=%d.",
                checkpoint_step,
                error_last_step,
                resume_step,
            )

        state["resume_step"] = int(resume_step)

        self._load_parameter_files(self.filenames["fixed_prm"], self.filenames["dynamic_prm"])
        self._load_target_from_init_files()
        self._load_bond_mask(self.filenames["bond_mask"])

        expected = int(state.get("n_bonds", self.target["BONDS"].shape[0]))
        current = int(self.target["BONDS"].shape[0])
        if expected != current:
            raise IOError(
                "Resume validation failed: checkpoint n_bonds {} does not match current target {}".format(
                    expected, current
                )
            )

        self.run_uuid = state.get("run_uuid", getattr(self, "run_uuid", None))
        logger.info(
            "Resume validation successful: checkpoint step=%d error.dat last step=%d",
            checkpoint_step,
            error_last_step,
        )
        return state

    def _can_restart(self):
        required = (
            self.filenames["init_avg_ic"],
            self.filenames["init_fluct_ic"],
            self.filenames["xplor_psf_file"],
            self.filenames["crd_file"],
        )
        for required_file in required:
            if not path.exists(required_file):
                return False
            try:
                if os.path.getsize(required_file) == 0:
                    return False
            except OSError:
                return False
        return True


    def _can_resume(self):
        required = (
            self.filenames["checkpoint"],
            self.filenames["fixed_prm"],
            self.filenames["dynamic_prm"],
            self.filenames["bond_mask"],
            self.filenames["error_data"],
            self.filenames["init_avg_ic"],
            self.filenames["init_fluct_ic"],
        )
        for required_file in required:
            if not path.exists(required_file):
                return False
            try:
                if os.path.getsize(required_file) == 0:
                    return False
            except OSError:
                return False

        try:
            state = self._read_json(self.filenames["checkpoint"])
            checkpoint_step = int(state["last_completed_step"])
            error_last_step = self._get_error_last_step()
            if checkpoint_step <= 0 or error_last_step <= 0:
                return False
        except Exception:
            return False

        return True


    def _resolve_requested_mode(self, restart=False, resume=False, auto=False):
        selected = sum(bool(x) for x in (restart, resume, auto))
        if selected > 1:
            raise ValueError("restart, resume, and auto are mutually exclusive.")
        if resume:
            return "resume"
        if restart:
            return "restart"
        if auto:
            return "auto"
        return "fresh"


    def _resolve_auto_mode(self):
        if self._can_resume():
            return "resume"
        if self._can_restart():
            return "restart"
        return "fresh"  

    def _compute_relative_diff(self, fluct_diff_series, tol):
        fluct_diff_series = pd.to_numeric(fluct_diff_series, errors="coerce")
        return (fluct_diff_series - tol) / tol

    def _log_tail_convergence_debug(self, fluct_diff, relative_diff, tol):
        n_remaining = int(self.unconverged_bond_list.sum())
        n_total = self.unconverged_bond_list.size

        logger.info("Remaining unconverged bonds: %d / %d", n_remaining, n_total)

        if n_remaining > n_total * self.TAIL_FRACTION:
            return False

        logger.info("Tail convergence region reached.")

        mask_tail = self.unconverged_bond_list & (
            relative_diff > self.TAIL_RELATIVE_DIFF_THRESHOLD
        )
        skipped_mask = self.unconverged_bond_list & (~mask_tail)

        if int(skipped_mask.sum()) > 0:
            skipped_df = pd.DataFrame({
                "fluct_diff": fluct_diff[skipped_mask],
                "relative_diff": relative_diff[skipped_mask],
            })
            skipped_df.index.names = self.bond_def

            logger.info(
                "Tail relative-difference rule skipped %d bonds:\n%s",
                int(skipped_mask.sum()),
                skipped_df.to_string()
            )

            skipped_file = path.join(self.outdir, "tail_skipped_bonds.txt")
            skipped_df.to_csv(skipped_file, sep=" ")
            logger.info("Skipped tail bonds written to %s", skipped_file)

        self.unconverged_bond_list = mask_tail

        if self.unconverged_bond_list.sum() == 0:
            logger.info("Checking relative difference: All bonds converged, exiting")
            return True

        return False

    def _reset_run_outputs(self):
        self.error.loc[:, :] = 0
        self._write_error_header()
        if path.exists(self.filenames["checkpoint"]):
            os.remove(self.filenames["checkpoint"])
        if path.exists(self.filenames["bond_mask"]):
            os.remove(self.filenames["bond_mask"])
        if path.exists(self.filenames["bond_convergence"]):
            os.remove(self.filenames["bond_convergence"])

    def _write_bond_convergence_cycle(self, fluct_diff, cycle_idx):
        cycle_column = str(cycle_idx)
        cycle_data = fluct_diff.round(6).rename(cycle_column)
        cycle_data.index.names = self.bond_def

        if path.exists(self.filenames["bond_convergence"]):
            convergence = pd.read_csv(
                self.filenames["bond_convergence"],
                index_col=list(range(len(self.bond_def))),
            )
            convergence.index.names = self.bond_def
            convergence = convergence.reindex(cycle_data.index)
            convergence[cycle_column] = cycle_data
        else:
            convergence = cycle_data.to_frame()

        convergence.to_csv(self.filenames["bond_convergence"])

    # Initial average and fluctuation values are calculated 
    # Using ic dyna aver and ic dyna fluct in CHARMM

    def initialize(self, nma_exec=None):
        """Prepare fresh fluctuation-matching initialization state."""

        if not path.exists(self.filenames["init_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = ("dimension chsize 8000000" if version >= 36 else "")
            with open(self.filenames["init_input"], mode="wb") as charmm_file:
                logger.info("Writing CHARMM input file.")
                charmm_inp = charmm_init.init.format(
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames
                )
                charmm_inp = textwrap.dedent(charmm_inp[1:])
                charmm_file.write(charmm_inp.encode())

        charmm_exec = (
            os.environ.get("CHARMMEXEC", util.which("charmm"))
            if nma_exec is None else nma_exec
        )
        with open(self.filenames["init_log"], "w") as log_file:
            subprocess.check_call(
                [charmm_exec, "-i", self.filenames["init_input"]],
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )

        for required in (self.filenames["init_avg_ic"], self.filenames["init_fluct_ic"]):
            if not path.exists(required):
                raise IOError(
                    "CHARMM initialization did not produce required file: {}".format(required)
                )
            if os.path.getsize(required) == 0:
                raise IOError(
                    "CHARMM initialization produced an empty IC file: {}".format(required)
                )
        self.run_uuid = str(uuid.uuid4())
        self._build_initial_state_from_init_ic()
        
        with mda.Writer(self.filenames["fixed_prm"], **self.kwargs) as prm:
            logger.info("Writing %s...", self.filenames["fixed_prm"])
            prm.write(self.parameters)
        with mda.Writer(self.filenames["dynamic_prm"], **self.kwargs) as prm:
            logger.info("Writing %s...", self.filenames["dynamic_prm"])
            prm.write(self.dynamic_params)


    # Fluctuation Matching starts here by calculating 
    # NMA (Using vibran nmode in CHARMM) from ENM and 
    # fluctuations are matched with CG initial values untill convergence 
    def run(
        self,
        nma_exec=None,
        tol=1.e-3,
        n_cycles=300,
        low_bound=0.0,
        restart=False,
        resume=False,
        auto=False,
    ):     
        """Perform self-consistent fluctuation matching.

        Parameters
        ----------
        nma_exec : str
            Executable file for normal mode analysis.
        tol : float, optional
            Fluctuation-difference tolerance.
        n_cycles : int, optional
            Total number of fluctuation-matching cycles desired for the run.
            When `resume=True`, the code continues from `last_completed_step + 1`
            up to `n_cycles`, rather than running `n_cycles` additional cycles.
        low_bound : float, optional
            Lowest Kb value retained to reduce noise.
        restart : bool, optional
            Rebuild the initial fluctuation-matching state from
            `init.average.ic` and `init.fluct.ic`, reset run outputs, and start
            the fluctuation-matching cycle again from step 1.
        resume : bool, optional
            Resume from the last completed checkpointed fluctuation-matching
            cycle by loading the current parameter files, `bond_mask.txt`,
            and the target reconstructed from `init.average.ic` and
            `init.fluct.ic`.
        """

        charmm_exec = (
            os.environ.get("CHARMMEXEC", util.which("charmm"))
            if nma_exec is None else nma_exec
        )
        if charmm_exec is None:
            logger.exception(
                "Please set CHARMMEXEC with the location of your CHARMM executable file or add the charmm path to your PATH environment."
            )
            raise OSError(
                    "Please set CHARMMEXEC with the location of your CHARMM executable file or add the charmm path to your PATH environment."
            )

        requested_mode = self._resolve_requested_mode(
            restart=restart, resume=resume, auto=auto
        )

        if requested_mode == "resume":
            state = self._load_resume_state()
            last_completed_step = int(state.get("resume_step", state["last_completed_step"]))
            self.restart = False
            self.resume = True
            self.auto = False
            logger.info("Run mode selected: resume")

        elif requested_mode == "restart":
            self._load_restart_state()
            self._reset_run_outputs()
            last_completed_step = 0
            self.restart = True
            self.resume = False
            self.auto = False
            logger.info("Run mode selected: restart")

        elif requested_mode == "auto":
            self.auto = True
            self.restart = False
            self.resume = False
            resolved_mode = self._resolve_auto_mode()
            logger.info("Auto run mode resolved to: %s", resolved_mode)

            if resolved_mode == "resume":
                try:
                    state = self._load_resume_state()
                    last_completed_step = int(state.get("resume_step", state["last_completed_step"]))
                    self.resume = True
                    logger.info("Auto mode using resume state.")
                except Exception as exc:
                    logger.warning(
                        "Auto mode could not load resume state (%s). Falling back.",
                        exc
                    )
                    if self._can_restart():
                        self._load_restart_state()
                        self._reset_run_outputs()
                        last_completed_step = 0
                        self.restart = True
                        logger.info("Auto mode falling back to restart state.")
                    else:
                        if not self.parameters:
                            try:
                                self.initialize(nma_exec)
                            except IOError:
                                raise IOError("Some files are missing. Unable to initialize.")
                        self._reset_run_outputs()
                        self._initialize_convergence_mask()
                        last_completed_step = 0
                        logger.info("Auto mode falling back to fresh state.")

            elif resolved_mode == "restart":
                self._load_restart_state()
                self._reset_run_outputs()
                last_completed_step = 0
                self.restart = True
                logger.info("Auto mode using restart state.")

            else:
                if not self.parameters:
                    try:
                        self.initialize(nma_exec)
                    except IOError:
                        raise IOError("Some files are missing. Unable to initialize.")
                self._reset_run_outputs()
                self._initialize_convergence_mask()
                last_completed_step = 0
                logger.info("Auto mode using fresh state.")

        else:  # start fresh
            if not self.parameters:
                try:
                    self.initialize(nma_exec)
                except IOError:
                    raise IOError("Some files are missing. Unable to initialize.")
            self._reset_run_outputs()
            self._initialize_convergence_mask()
            last_completed_step = 0
            self.restart = False
            self.resume = False
            self.auto = False
            logger.info("Run mode selected: fresh")

        # Check whether the requested cylce is already completed and exit early if so
        if last_completed_step >= n_cycles:
            logger.info(
                f"Requested {n_cycles} cycles already completed for this window "
                f"(last_completed_step={last_completed_step}). Skipping."
            )
            return

        if not path.exists(self.filenames["charmm_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = ("dimension chsize 8000000" if version >= 36 else "")
            with open(self.filenames["charmm_input"], mode="wb") as charmm_file:
                logger.info("Writing CHARMM input file.")
                charmm_inp = charmm_nma.nma.format(
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames
                )
                charmm_inp = textwrap.dedent(charmm_inp[1:])
                charmm_file.write(charmm_inp.encode())

        self.target["BONDS"].set_index(self.bond_def, inplace=True)
        bond_values = self.target["BONDS"].columns

        logger.info(f"Starting fluctuation matching --> {n_cycles - last_completed_step} iterations to run")
        if low_bound != 0.0:
            logger.info("Lower bound after 75%% iteration is set to %s", low_bound)

        st = time.perf_counter()

        for cycle_idx in range(last_completed_step + 1, n_cycles + 1):
            cycle_t0 = time.perf_counter()
            timings = {}
            self.error["step"] = cycle_idx

            t0 = time.perf_counter()
            self._run_charmm(charmm_exec)
            timings["charmm"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            self.dynamic_params["BONDS"].set_index(self.bond_def, inplace=True)
            self.parameters["BONDS"].set_index(self.bond_def, inplace=True)
            timings["set_index"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            avg_ic, fluct_ic = self._read_cycle_ic()
            timings["read_ic_total"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            vib_ic = self._build_vib_ic(fluct_ic, avg_ic, bond_values)
            timings["concat"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            fluct_diff = self._compute_convergence_and_error(
                vib_ic=vib_ic,
                bond_values=bond_values,
                tol=tol,
                i=cycle_idx - 1,
            )
            self._write_bond_convergence_cycle(fluct_diff, cycle_idx)
            timings["convergence"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            vib_ic = self._optimize_bonds(
                vib_ic=vib_ic,
                bond_values=bond_values,
                low_bound=low_bound,
                n_cycles=n_cycles,
                i=cycle_idx - 1,
            )
            timings["optimization"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            self._apply_updated_parameters(vib_ic, bond_values)
            timings["update_tables"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            self._write_parameter_files()
            timings["write_prm"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            self._append_error_row()
            timings["write_error"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            self._write_cycle_checkpoint(cycle_idx, tol, n_cycles, low_bound)
            timings["checkpoint"] = time.perf_counter() - t0

            timings["cycle_total"] = time.perf_counter() - cycle_t0

            logged_time_dt = datetime.now()
            estimated_next_end = logged_time_dt + timedelta(seconds=timings["cycle_total"])

            logger.info(
                "Fluctuation matching cycle %d timing | %s | logged_time=%s | est_next_cycle_end=%s",
                cycle_idx,
                _fmt_timing(timings),
                logged_time_dt.strftime("%Y-%m-%d %H:%M:%S"),
                estimated_next_end.strftime("%Y-%m-%d %H:%M:%S"),
            )
            logger.info(
                "%d not converged out of %d",
                self.unconverged_bond_list.sum(), self.unconverged_bond_list.size
            )

            relative_diff = self._compute_relative_diff(fluct_diff, tol)
            if self._log_tail_convergence_debug(fluct_diff, relative_diff, tol):
                self._write_cycle_checkpoint(cycle_idx, tol, n_cycles, low_bound)
                break

        logger.info(
            "Fluctuation matching completed in %.3fs",
            time.perf_counter() - st
        )
        self.target["BONDS"].reset_index(inplace=True)

    def calculate_thermo(self, nma_exec=None):
        """Calculate the thermodynamic properties of the trajectory."""
        charmm_exec = (
            os.environ.get("CHARMMEXEC", util.which("charmm"))
            if nma_exec is None else nma_exec
        )
        if charmm_exec is None:
            logger.exception(
                "Please set CHARMMEXEC with the location of your CHARMM executable file or add the charmm path to your PATH environment."
            )
            raise OSError(
                    "Please set CHARMMEXEC with the location of your CHARMM executable file or add the charmm path to your PATH environment."
                    )

        if not path.exists(self.filenames["thermo_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = ("dimension chsize 500000 maxres 3000000" if version >= 36 else "")
            with open(self.filenames["thermo_input"], mode="wb") as charmm_file:
                logger.info("Writing CHARMM input file.")
                charmm_inp = charmm_thermo.thermodynamics.format(
                    trajectory=path.join(self.outdir, self.args[-1]),
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames
                )
                charmm_inp = textwrap.dedent(charmm_inp[1:])
                charmm_file.write(charmm_inp.encode())

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

        thermo = pd.DataFrame(thermo, columns=columns)
        thermo.drop(["RESN", "Atm/res", "Ign.frq"], axis=1, inplace=True)
        thermo.set_index(["segidI", "resI"], inplace=True)
        thermo = thermo.astype(float)

        with open(self.filenames["thermo_data"], "wb") as data_file:
            logger.info("Writing thermodynamics data file.")
            thermo = thermo.to_csv(
                index=True,
                sep=" ",
                float_format="%.4f",
                encoding="utf-8",
            )
            data_file.write(thermo.encode())
