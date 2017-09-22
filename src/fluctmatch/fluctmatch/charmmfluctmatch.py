# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import os
import subprocess
from os import path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.coordinates import memory
from future.builtins import (
    dict,
    open,
    super,
)
from future.utils import (
    native_str,
    raise_with_traceback,
)
from scipy import constants

from fluctmatch.intcor import IC
from fluctmatch.intcor import utils as icutils
from fluctmatch.models import enm
from fluctmatch.parameter import PRM
from fluctmatch.parameter import utils as prmutils
from . import base as fmbase
from . import utils as fmutils
from .data import (
    charmm36_nma,
    charmm_nma,
)


class CharmmFluctMatch(fmbase.FluctMatch):
    """Fluctuation matching using CHARMM."""
    charmm39 = False
    dynamic_params = dict()
    bond_def = ("I", "J")
    error_hdr = ("step", "r.m.s.d.", "Kb_err", "b0_err")

    def __init__(self, *args, **kwargs):
        """Initialization of fluctuation matching using the CHARMM program.

        Parameters
        ----------
        topology : filename or Topology object
            A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
            define the list of atoms. If the file includes bond information,
            partial charges, atom masses, ... then these data will be available to
            MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
            topology) is always required. Alternatively, an existing
            :class:`MDAnalysis.core.topology.Topology` instance may also be given.
        extended
            Renames the residues and atoms according to the extended CHARMM PSF format.
            Standard CHARMM PSF limits the residue and atom names to four characters,
            but the extended CHARMM PSF permits eight characters. The residues and
            atoms are renamed according to the number of segments (1: A, 2: B, etc.)
            and then the residue number or atom index number.
         xplor
            Assigns the atom type as either a numerical or an alphanumerical
            designation. CHARMM normally assigns a numerical designation, but the
            XPLOR version permits an alphanumerical designation with a maximum
            size of 4. The numerical form corresponds to the atom index number plus a
            factor of 100, and the alphanumerical form will be similar the standard
            CHARMM atom name.
        topology_format
            Provide the file format of the topology file; ``None`` guesses it from
            the file extension [``None``] Can also pass a subclass of
            :class:`MDAnalysis.topology.base.TopologyReaderBase` to define a custom
            reader to be used on the topology file.
        format
            Provide the file format of the coordinate or trajectory file; ``None``
            guesses it from the file extension. Note that this keyword has no
            effect if a list of file names is supplied because the "chained" reader
            has to guess the file format for each individual list member.
            [``None``] Can also pass a subclass of
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
            :class:`MDAnalysis.core.groups.AtomGroup` instances pickled from the
            Universe to only unpickle if a compatible Universe with matching
            *anchor_name* is found. Even if *anchor_name* is set *is_anchor* will
            still be honored when unpickling.
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
        charmm39
            Use CHARMM39 RTF file format.
        charmm36
            Use CHARMM36 PRM file format
        extended
            Use the extended format.
        cmap
            Include CMAP section.
        cheq
            Include charge equilibration.
        title
            Title lines at the beginning of the file.
        resid
            Include segment IDs in the internal coordinate files.
        nonbonded
            Include the nonbonded section in the parameter file.
        """
        super().__init__(*args, **kwargs)
        self.filenames = dict(
            init_avg_ic=path.join(self.outdir, "init.avg.ic"),
            init_fluct_ic=path.join(self.outdir, "init.fluct.ic"),
            avg_ic=path.join(self.outdir, "avg.ic"),
            fluct_ic=path.join(self.outdir, "fluct.ic"),
            dynamic_prm=path.join(self.outdir, "fluctmatch.dist.prm"),
            fixed_prm=path.join(self.outdir, "fluctmatch.prm"),
            psf_file=path.join(self.outdir, "fluctmatch.psf"),
            xplor_psf_file=path.join(self.outdir, "fluctmatch.xplor.psf"),
            crd_file=path.join(self.outdir, "fluctmatch.crd"),
            stream_file=path.join(self.outdir, "fluctmatch.stream"),
            topology_file=path.join(self.outdir, "fluctmatch.rtf"),
            nma_crd=path.join(self.outdir, "fluctmatch.vib.crd"),
            nma_vib=path.join(self.outdir, "fluctmatch.vib"),
            charmm_input=path.join(self.outdir, "fluctmatch.inp"),
            charmm_log=path.join(self.outdir, "fluctmatch.log"),
            error_data=path.join(self.outdir, "error.dat"),
        )

        # Boltzmann constant
        self.BOLTZ = self._temperature * (constants.k * constants.N_A / (constants.calorie * constants.kilo))

        # Bond factor mol^2-Ang./kcal^2
        self.KFACTOR = 0.02

        # Self consistent error information.
        self.error = pd.DataFrame(columns=self.error_hdr)

    def _restart(self):
        """Check whether to restart from a previous run."""
        # Load target data
        if not self.target:
            try:


    def initialize(self, restart=False):
        """Create an elastic network model from a basic coarse-grain model.

        :param restart: reinitialize tye object
        """
        if not restart:
            # Initialize variables and load the universe.
            title = self.kwargs.get("title")
            extended = self.kwargs.get("extended", True)
            nonbonded = self.kwargs.get("nonbonded", False)
            cheq = self.kwargs.get("cheq", True)
            cmap = self.kwargs.get("cmap", True)
            charmm39 = self.kwargs.get("charmm39", True)
            charmm36 = self.kwargs.get("charmm36", True)
            universe = enm.Enm(*self.args, **self.kwargs)

            # Write required CHARMM input files.
            with mda.Writer(self.filenames["topology_file"], title=title, charmm39=charmm39) as rtf:
                rtf.write(universe)
            with mda.Writer(self.filenames["stream_file"]) as stream:
                stream.write(universe)
            with mda.Writer(self.filenames["psf_file"], title=title, extended=extended, cheq=cheq, cmap=cmap) as psf:
                psf.write(universe)

            # Calculate the average coordinates, average bond lengths, and fluctuations of bond lengths from the trajectory.
            positions = fmutils.average_structure(*self.args, **self.kwargs)
            avg_universe = mda.Universe(
                self.filenames["psf_file"],
                [positions, ],
                format=memory.MemoryReader, order="fac"
            )
            with mda.Writer(self.filenames["crd_file"], title=title, dt=1.0) as crd:
                crd.write(avg_universe.atoms)

            # Create and write initial internal coordinate files.
            avg_bonds = fmutils.average_bonds(*args, **self.kwargs).set_index(self.bond_def)
            std_bonds = fmutils.bond_fluctuation(*args, **self.kwargs).set_index(self.bond_def)
            avg_table = icutils.create_empty_table(universe.atoms)
            hdr = avg_table.columns
            avg_table.set_index(self.bond_def, inplace=True)
            avg_table["r_IJ"] = avg_bonds["r_IJ"]
            avg_table = avg_table.reset_index()[hdr]
            fluct_table = icutils.create_empty_table(universe.atoms)
            fluct_table["r_IJ"] = std_bonds["r_IJ"]
            fluct_table = fluct_table.reset_index()[hdr]
            with mda.Writer(self.filenames["init_avg_ic"], title=title, extended=extended, resid=resid) as table:
                table.write(avg_table)
            with mda.Writer(self.filenames["init_fluct_ic"], title=title, extended=extended, resid=resid) as table:
                table.write(fluct_table)

            # Write the parameter files.
            if self.charmm39 and not charmm36:
                charmm36 = True

            target = pd.concat([std_bonds, avg_bonds], axis=1)
            self.target.update(prmutils.create_empty_parameters(universe))
            target.columns = self.target["BONDS"].columns
            self.target["BONDS"] = target
            self.parameters.update(self.target)
            self.parameters["BONDS"]["Kb"] = self.BOLTZ / np.square(self.parameters["BONDS"]["Kb"])
            with mda.Writer(self.filenames["fixed_prm"], title=title, charmm36=charmm36, nonbonded=nonbonded) as prm:
                prm.write(self.parameters)

            # Write an XPLOR version of the PSF
            universe.atoms.types = universe.atoms.names
            with mda.Writer(self.filenames["xplor_psf_file"], title=title, extended=extended, cheq=cheq,
                            cmap=cmap) as psf:
                psf.write(universe)
        else:
            try:
                # Read the parameter files.
                with PRM.ParamReader(self.filenames["fixed_prm"]) as fixed:
                    self.parameters.update(fixed.read())
                with PRM.ParamReader(self.filenames["dynamic_prm"]) as dynamic:
                    self.dynamic_params.update(dynamic.read())

                # Read the initial internal coordinate files.
                with IC.IntcorReader(self.filenames["init_avg_ic"]) as init_avg:
                    avg_table = init_avg.read().set_index(self.bond_def)["r_IJ"]
                with IC.IntcorReader(self.filenames["init_fluct_ic"]) as init_fluct
                    fluct_table = init_fluct.read().set_index(self.bond_def)["r_IJ"]
                table = pd.concat([fluct_table, avg_table], axis=1)

                # Set the target fluctuation values.
                self.target.update(self.parameters)
                self.target["BONDS"].set_index(self.bond_def, inplace=True)
                table.columns = self.target["BONDS"].columns
                self.target["BONDS"] = table.copy(deep=True)
                self.target["BONDS"].reset_index(inplace=True)
            except IOError:
                raise_with_traceback((IOError("Some files are missing. Unable to restart.")))

    def run(self, nma_exec=None, tol=1.e-4, n_cycles=250):
        """Perform a self-consistent fluctuation matching

        :param nma_exec: executable file for normal mode analysis
        :param tol: error tolerance
        :param n_cycles: number of fluctuation matching cycles
        """
        try:
            charmm_exec = os.environ["CHARMMEXEC"] if nma_exec is None else nma_exec
        except KeyError:
            raise_with_traceback(OSError("Please set CHARMMEXEC with the location of your CHARMM executable file."))


        # Read the parameters
        if not self.parameters:
            try:
                self.initialize(restart=True)
            except IOError:
                raise_with_traceback((IOError("Some files are missing. Unable to restart.")))

        # Write CHARMM input file.
        if not path.exists(self.filenames["charmm_inp"]):
            with open(self.filenames["charmm_inp"], "wt") as charmm_file:
                charmm_inp = (
                    charmm36_nma.nma.format(temperature=self.temperature, **self.filenames)
                    if self.charmm39
                    else charmm_nma.nma.format(temperature=self.temperature, **self.filenames)
                )
                print(charmm_inp, file=charmm_file)

        # Set the indices for the parameter tables.
        self.target["BONDS"].set_index(self.bond_def, inplace=True)

        # Check for restart.
        if path.exists(self.filenames["error_data"]):
            with open(self.filenames["error_data"], "r") as data:
                error_info = pd.read_csv(data, header=0, skipinitialspace=True, delim_whitespace=True)
                self.error["step"] = error_info["step"].size + 1
        else:
            with open(self.filenames["error_data"], "a") as data:
                np.savetxt(data, self.error_hdr, fmt=native_str("%10s"))
            self.error["step"] = 1

        # Run simulation
        while self.error["step"] <= n_cycles or self.error["r.m.s.d."] > tol:
            subprocess.check_call([charmm_exec, "-i", charmm_file.name], stdout=self.filenames["charmm_log"])
            self.dynamic_params["BONDS"].set_index(self.bond_def, inplace=True)
            self.parameters["BONDS"].set_index(self.bond_def, inplace=True)

            # Read the average bond distance.
            with IC.IntcorReader(self.filenames["avg_ic"]) as intcor:
                avg_ic = intcor.read().set_index(self.bond_def)["r_IJ"]

            # Read the bond fluctuations.
            with IC.IntcorReader(self.filenames["fluct_ic"]) as intcor:
                fluct_ic = intcor.read().set_index(self.bond_def)["r_IJ"]

            vib_ic = pd.concat([fluct_ic, avg_ic], axis=1)
            vib_ic.columns = self.dynamic_params["BONDS"].columns

            vib_error = np.square(self.target["BONDS"] - vib_ic).mean(axis=0)
            vib_error = np.sqrt(vib_error)
            vib_error.columns = self.error_hdr[2:]
            error[header[2:]] - vib_error.copy()


            # Calculate the new force constant.
            optimized = 1. / np.square(self.target["BONDS"]["Kb"])
            optimized -= 1. / np.square(vib_ic["Kb"])
            vib_ic["Kb"] = self.parameters["BONDS"]["Kb"] - (optimized * self.BOLTZ * self.KFACTOR)

            # r.m.s.d. between previous and current force constant
            diff = np.square(self.parameters["BONDS"]["Kb"] - vib_ic["Kb"])
            error["r.m.s.d."] = np.sqrt(diff.mean(axis=0)).values

            # Update the parameters and write to file.
            self.parameters["BONDS"]["Kb"] = vib_ic["Kb"].copy(deep=True)
            self.dynamic_params["BONDS"] = vib_ic.copy(deep=True)
            self.parameters["BONDS"].reset_index(inplace=True)
            self.dynamic_params["BONDS"].reset_index(inplace=True)
            with mda.Writer(self.filenames["fixed_prm"]) as prm:
                prm.write(self.parameters)
            with mda.Writer(self.filenames["dynamic_prm"]) as prm:
                prm.write(self.dynamic_params)

            # Update the error values.
            with open(self.filenames["error_data"], "a") as error_file
                np.savetxt(error_file, self.error, fmt=native_str("%10d%10.6f%10.6f%10.6f"))
            self.error["step"] += 1
