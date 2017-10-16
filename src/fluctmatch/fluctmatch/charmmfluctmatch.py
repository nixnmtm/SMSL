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
"""Fluctuation matching using CHARMM.

Notes
------
For CHARMM to work with the fluctuation matching code, it must be
recompiled with some modifications to the source code. `ATBMX` and
`MAXCB` (located in dimens.fcm [c35] or dimens_ltm.src [c41]) must
be increased. `ATBMX` determines the number of bonds allowed per
atom and `MAXCB` determines the maximum number of bond parameters
in the CHARMM parameter file. Additionally, `CHSIZE` may need to
be increased if using an earlier version (< c36).

"""

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import copy
import os
import subprocess
import textwrap
import time
from os import path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from MDAnalysis.coordinates.core import reader
from future.builtins import (
    dict,
    super,
)
from future.utils import (
    PY2,
    native_str,
    raise_with_traceback,
)
from scipy import constants

from fluctmatch.fluctmatch import base as fmbase
from fluctmatch.fluctmatch import utils as fmutils
from fluctmatch.fluctmatch.data import (
    charmm_nma,
)
from fluctmatch.intcor import utils as icutils
from fluctmatch.parameter import utils as prmutils

if PY2:
    FileNotFoundError = IOError


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
            init_avg_ic=path.join(self.outdir, "init.average.ic"),
            init_fluct_ic=path.join(self.outdir, "init.fluct.ic"),
            avg_ic=path.join(self.outdir, "average.ic"),
            fluct_ic=path.join(self.outdir, "fluct.ic"),
            dynamic_prm=path.join(
                self.outdir, "{}.dist.prm".format(self.prefix)
            ),
            fixed_prm=path.join(self.outdir, ".".join((self.prefix, "prm"))),
            psf_file=path.join(self.outdir, ".".join((self.prefix, "psf"))),
            xplor_psf_file=path.join(
                self.outdir, ".".join((self.prefix, "xplor", "psf"))
            ),
            crd_file=path.join(self.outdir, ".".join((self.prefix, "cor"))),
            stream_file=path.join(
                self.outdir, ".".join((self.prefix, "stream"))
            ),
            topology_file=path.join(
                self.outdir, ".".join((self.prefix, "rtf"))
            ),
            nma_crd=path.join(
                self.outdir, ".".join((self.prefix, "mini", "cor"))
            ),
            nma_vib=path.join(self.outdir, ".".join((self.prefix, "vib"))),
            charmm_input=path.join(self.outdir, ".".join((self.prefix, "inp"))),
            charmm_log=path.join(self.outdir, ".".join((self.prefix, "log"))),
            error_data=path.join(self.outdir, "error.dat"),
        )

        # Boltzmann constant
        self.BOLTZ = self.temperature * (
            constants.k * constants.N_A / (constants.calorie * constants.kilo)
        )

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
        table.drop(["r_IJ", ], axis=1, inplace=True)
        table = pd.concat([table, data["r_IJ"]], axis=1)
        return table.reset_index()[hdr]

    def initialize(self, restart=False):
        """Create an elastic network model from a basic coarse-grain model.

        Parameters
        ----------
        restart : bool, optional
            Reinitialize the object by reading files instead of doing initial
            calculations.
        """
        if not restart:
            universe = mda.Universe(*self.args, **self.kwargs)

            # Create and write initial internal coordinate files.
            print("Determining the average bond distances...")
            avg_bonds, std_bonds = fmutils.BondStats(
                universe, func="both"
            ).run().result
            print("Writing {}...".format(self.filenames["init_avg_ic"]))
            with mda.Writer(
                self.filenames["init_avg_ic"],
                **self.kwargs
            ) as table:
                avg_table = self._create_ic_table(universe, avg_bonds)
                table.write(avg_table)

            print("Determining the fluctuation of bond distances...")
            print("Writing {}...".format(self.filenames["init_fluct_ic"]))
            with mda.Writer(
                self.filenames["init_fluct_ic"],
                **self.kwargs
            ) as table:
                fluct_table = self._create_ic_table(universe, std_bonds)
                table.write(fluct_table)

            # Write the parameter files.
            print("Calculating the initial CHARMM parameters...")
            target = pd.concat([std_bonds, avg_bonds], axis=1).reset_index()
            self.target.update(
                prmutils.create_empty_parameters(universe, **self.kwargs)
            )
            target.columns = self.target["BONDS"].columns
            self.target["BONDS"] = target
            self.parameters = copy.deepcopy(self.target)
            self.parameters["BONDS"]["Kb"] = (
                self.BOLTZ / self.parameters["BONDS"]["Kb"].apply(np.square)
            )
            self.dynamic_params = copy.deepcopy(self.parameters)
            print("Writing {}...".format(self.filenames["fixed_prm"]))
            with mda.Writer(self.filenames["fixed_prm"], **self.kwargs) as prm:
                prm.write(self.parameters)
            print("Writing {}...".format(self.filenames["dynamic_prm"]))
            with mda.Writer(
                self.filenames["dynamic_prm"],
                **self.kwargs
            ) as prm:
                prm.write(self.dynamic_params)
        else:
            try:
                # Read the parameter files.
                print("Loading parameter and internal coordinate files.")
                with reader(self.filenames["fixed_prm"]) as fixed:
                    self.parameters.update(fixed.read())
                with reader(self.filenames["dynamic_prm"]) as dynamic:
                    self.dynamic_params.update(dynamic.read())

                # Read the initial internal coordinate files.
                with reader(self.filenames["init_avg_ic"]) as init_avg:
                    avg_table = init_avg.read().set_index(self.bond_def)["r_IJ"]
                with reader(self.filenames["init_fluct_ic"]) as init_fluct:
                    fluct_table = (
                        init_fluct.read().set_index(self.bond_def)["r_IJ"]
                    )
                table = pd.concat([fluct_table, avg_table], axis=1)

                # Set the target fluctuation values.
                print("Files loaded successfully...")
                self.target = copy.deepcopy(self.parameters)
                self.target["BONDS"].set_index(self.bond_def, inplace=True)
                table.columns = self.target["BONDS"].columns
                self.target["BONDS"] = table.copy(deep=True).reset_index()
            except (FileNotFoundError, IOError):
                raise_with_traceback(
                    (IOError("Some files are missing. Unable to restart."))
                )

    def run(self, nma_exec=None, tol=1.e-4, n_cycles=250):
        """Perform a self-consistent fluctuation matching.

        Parameters
        ----------
        nma_exec : str
            executable file for normal mode analysis
        tol : float, optional
            error tolerance
        n_cycles : int, optional
            number of fluctuation matching cycles
        """
        try:
            charmm_exec = (
                os.environ["CHARMMEXEC"] if nma_exec is None else nma_exec
            )
        except KeyError:
            raise_with_traceback(
                OSError("Please set CHARMMEXEC with the location of your "
                        "CHARMM executable file.")
            )

        # Read the parameters
        if not self.parameters:
            try:
                self.initialize(restart=True)
            except IOError:
                raise_with_traceback(
                    (IOError("Some files are missing. Unable to restart."))
                )

        # Write CHARMM input file.
        if not path.exists(self.filenames["charmm_input"]):
            version = self.kwargs.get("charmm_version", 41)
            dimension = (
                "dimension chsize 500000 maxres 3000000"
                if version >= 36
                else ""
            )
            with util.openany(
                self.filenames["charmm_input"],
                mode="w"
            ) as charmm_file:   # type: Optional[IO[str]]
                charmm_inp = charmm_nma.nma.format(
                    temperature=self.temperature,
                    flex="flex" if version else "",
                    version=version,
                    dimension=dimension,
                    **self.filenames)
                charmm_inp = textwrap.dedent(charmm_inp[1:])
                print(charmm_inp, file=charmm_file)

        # Set the indices for the parameter tables.
        self.target["BONDS"].set_index(self.bond_def, inplace=True)
        bond_values = self.target["BONDS"].columns

        # Check for restart.
        try:
            if os.stat(self.filenames["error_data"]).st_size > 0:
                with util.openany(self.filenames["error_data"], "r") as data:
                    error_info = pd.read_csv(
                        data,
                        header=0,
                        skipinitialspace=True,
                        delim_whitespace=True
                    )
                    if not error_info.empty:
                        self.error["step"] = error_info["step"].values[-1]
            else:
                raise FileNotFoundError
        except (FileNotFoundError, OSError):
            with util.openany(self.filenames["error_data"], "w") as data:
                np.savetxt(
                    data,
                    [self.error_hdr,],
                    fmt=native_str("%10s"),
                    delimiter=native_str("")
                )
        self.error["step"] += 1

        # Run simulation
        print("Starting fluctuation matching")
        st = time.time()
        while (self.error["step"] <= n_cycles).bool():
            subprocess.check_call([
                "charmm",
                "-i", self.filenames["charmm_input"],
                "-o", self.filenames["charmm_log"]
            ])
            self.dynamic_params["BONDS"].set_index(self.bond_def, inplace=True)
            self.parameters["BONDS"].set_index(self.bond_def, inplace=True)

            # Read the average bond distance.
            with reader(self.filenames["avg_ic"]) as intcor:
                avg_ic = intcor.read().set_index(self.bond_def)["r_IJ"]

            # Read the bond fluctuations.
            with reader(self.filenames["fluct_ic"]) as intcor:
                fluct_ic = intcor.read().set_index(self.bond_def)["r_IJ"]

            vib_ic = pd.concat([fluct_ic, avg_ic], axis=1)
            vib_ic.columns = bond_values

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
            vib_ic[bond_values[0]] = (
                self.parameters["BONDS"][bond_values[0]] -
                optimized[bond_values[0]]
            )
            vib_ic[bond_values[0]] = (
                vib_ic[bond_values[0]].where(vib_ic[bond_values[0]] >= 0., 0.)
            )

            # r.m.s.d. between previous and current force constant
            diff = self.dynamic_params["BONDS"] - vib_ic
            diff = diff.apply(np.square).mean(axis=0)
            diff = np.sqrt(diff)
            self.error[self.error.columns[1]] = diff.values[0]

            # Update the parameters and write to file.
            self.parameters["BONDS"][bond_values[0]] = (
                vib_ic[bond_values[0]].copy(deep=True)
            )
            self.dynamic_params["BONDS"] = vib_ic.copy(deep=True)
            self.parameters["BONDS"].reset_index(inplace=True)
            self.dynamic_params["BONDS"].reset_index(inplace=True)
            with mda.Writer(self.filenames["fixed_prm"], **self.kwargs) as prm:
                prm.write(self.parameters)
            with mda.Writer(
                self.filenames["dynamic_prm"],
                **self.kwargs
            ) as prm:
                prm.write(self.dynamic_params)

            # Update the error values.
            with util.openany(self.filenames["error_data"], "a") as error_file:
                np.savetxt(
                    error_file,
                    self.error,
                    fmt=native_str("%10d%10.6f%10.6f%10.6f",),
                    delimiter=native_str(""),
                )

            if (self.error[self.error.columns[1]] < tol).bool():
                break
            self.error["step"] += 1

        print(
            "Fluctuation matching completed in {:.6f}".format(time.time() - st)
        )
        self.target["BONDS"].reset_index(inplace=True)
