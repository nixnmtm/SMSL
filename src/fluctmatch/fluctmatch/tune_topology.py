
from future.utils import (PY2, native_str)


import os.path as path
import os
import logging
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.core import reader
from MDAnalysis.core import topologyattrs
from fluctmatch.parameter import PRM
from fluctmatch.intcor import IC
from fluctmatch.topology import RTF, PSFParser, STR

logger = logging.getLogger(__name__)


class TopologyTuning(object):

    """
    Tune Topology based on the fluctuation difference, distance cutoff and force constant.

    Returns:
        Writes necessary topology and parameter files
        Returns the new universe after trimming

    """

    def __init__(self, **kwargs):

        self.outdir = kwargs.get("outdir", os.getcwd())
        self.prefix = kwargs.get("prefix", "fluctmatch")
        self.rmsf_cutoff = kwargs.pop("rcut", 0.)
        self.dist_cutoff = kwargs.pop("dcut", 6.5)
        self.kb_cutoff = kwargs.pop("kbcut", 0.)

        self.pair_idx = ["I", "J"]

        self.filenames = dict(
            init_avg_ic=path.join(self.outdir, "init.average.ic"),
            init_fluct_ic=path.join(self.outdir, "init.fluct.ic"),
            dynamic_prm=path.join(self.outdir, "{}.dist.prm".format(self.prefix)),
            fixed_prm=path.join(self.outdir, ".".join((self.prefix, "prm"))),
            psf_file=path.join(self.outdir, ".".join((self.prefix, "psf"))),
            xplor_psf_file=path.join(self.outdir, ".".join((self.prefix, "xplor", "psf"))),
            stream_file=path.join(self.outdir, ".".join((self.prefix, "stream"))),
            topology_file=path.join(self.outdir, ".".join((self.prefix, "rtf"))),
            traj_file=path.join(self.outdir, "cg.dcd"),
        )

        self.universe = mda.Universe(self.filenames["xplor_psf_file"], self.filenames["traj_file"])

        # Attempt to create the necessary subdirectory
        try:
            self.writedir = os.path.join(self.outdir, "trimmed")
            os.makedirs(self.writedir)
        except FileExistsError:
            # directory already exists
            pass

        self.write_filenames = dict(
            init_avg_ic=path.join(self.writedir, "init.average.ic"),
            init_fluct_ic=path.join(self.writedir, "init.fluct.ic"),
            dynamic_prm=path.join(self.writedir, "{}.dist.prm".format(self.prefix)),
            fixed_prm=path.join(self.writedir, ".".join((self.prefix, "prm"))),
            psf_file=path.join(self.writedir, ".".join((self.prefix, "psf"))),
            xplor_psf_file=path.join(self.writedir, ".".join((self.prefix, "xplor", "psf"))),
            stream_file=path.join(self.writedir, ".".join((self.prefix, "stream"))),
            topology_file=path.join(self.writedir, ".".join((self.prefix, "rtf"))),
        )

    def trim_get_new_bondlist_write_prm(self, atomfluc_diff, **kwargs):

        # READ fixed and dynamic parameter files and initial ic files
        params = dict()
        with reader(self.filenames["fixed_prm"]) as fixed:
            params.update(fixed.read())

        dyn_params = dict()
        with reader(self.filenames["dynamic_prm"]) as dynamic:
            dyn_params.update(dynamic.read())

        with reader(self.filenames["init_fluct_ic"]) as icfile:
            std_bonds = icfile.read().set_index(self.pair_idx)

        with reader(self.filenames["init_avg_ic"]) as avgfile:
            avg_bonds = avgfile.read().set_index(self.pair_idx)

        # Get the bonds to be trimmed from fixed parameter file
        params["BONDS"].set_index(self.pair_idx, inplace=True)
        logger.info(f"Before trimming NBONDS: {params['BONDS'].shape}")
        count = 0
        bonds2trim = []
        for i in params["BONDS"].index.values:
            if atomfluc_diff.loc[i[0]] > self.rmsf_cutoff and atomfluc_diff.loc[i[1]] > self.rmsf_cutoff:
                if params["BONDS"].loc[i, 'Kb'] > self.kb_cutoff and params["BONDS"].loc[i, 'b0'] > self.dist_cutoff:
                    resI = self.universe.atoms.select_atoms(f"type {i[0]}").resids[0]
                    resJ = self.universe.atoms.select_atoms(f"type {i[1]}").resids[0]
                    if resI > resJ + 2 or resJ > resI + 2:
                        bonds2trim.append(i)
                        count += 1
        logger.info(f"Totally {count} bonds found")
        params["BONDS"].drop(bonds2trim, inplace=True)
        logger.info(f"After trimming NBONDS: {params['BONDS'].shape}")
        params["BONDS"].reset_index(inplace=True)

        # Write both fixed and dynamics
        with mda.Writer(self.write_filenames["fixed_prm"], **kwargs) as fixed_param:
            fixed_param.write(params)
        new_bondlist = params["BONDS"].set_index(self.pair_idx).index.values

        # Make NBONDS in dynamic parameters to match fixed parameter NBONDS
        # as the # of bonds selected from parameters and dynamic_parameter
        # could be diffrent due to its update average bond details in every iteration of FM

        dyn_params["BONDS"].set_index(self.pair_idx, inplace=True)
        dyn_params["BONDS"] = dyn_params["BONDS"].loc[new_bondlist]
        dyn_params["BONDS"].reset_index(inplace=True)
        assert dyn_params['BONDS'].shape[0] == params['BONDS'].shape[0], "Shape Mismatch."

        with mda.Writer(self.write_filenames["dynamic_prm"], **kwargs) as dynamic_param:
            dynamic_param.write(dyn_params)

        """ for matching the shape of initial and current parameter files for force constant calculation, 
            the atom pairs trimmed in parameter files should also be trimmed in ic files """

        cols = np.asarray([
            "segidI", "resI", "I", "segidJ", "resJ", "J", "segidK", "resK", "K",
            "segidL", "resL", "L", "r_IJ", "T_IJK", "P_IJKL", "T_JKL", "r_KL"
        ])
        std_bonds = std_bonds.loc[new_bondlist]
        std_bonds.reset_index(inplace=True)
        std_bonds = std_bonds[cols]
        avg_bonds = avg_bonds.loc[new_bondlist]
        avg_bonds.reset_index(inplace=True)
        avg_bonds = avg_bonds[cols]
        with mda.Writer(self.write_filenames["init_fluct_ic"], **kwargs) as fluct_init:
            fluct_init.write(std_bonds)
        with mda.Writer(self.write_filenames["init_avg_ic"], **kwargs) as avg_init:
            avg_init.write(avg_bonds)

        return new_bondlist, bonds2trim

    def get_atm_idx_aft_trim(self, atompairs):
        _map = dict(zip(self.universe.atoms.names, self.universe.atoms.indices))
        return [(_map[b[0]], _map[b[1]]) for b in atompairs]

    def trim_write_topology_files(self, **kwargs):
        """Write CHARMM coordinate, topology PSF, stream, and topology RTF files.

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or AtomGroup with bond definitions.
        outdir : str
            Location to write the files.
        prefix : str
            Prefix of filenames
        charmm_version
            Version of CHARMM for formatting (default: 41)
        extended
            Use the extended format.
        cmap
            Include CMAP section.
        cheq
            Include charge equilibration.
        title
            Title lines at the beginning of the file.
        """

        n_atoms = self.universe.atoms.n_atoms
        n_bonds = len(self.universe.bonds)
        n_angles = len(self.universe.angles)
        n_dihedrals = len(self.universe.dihedrals)
        n_impropers = len(self.universe.impropers)
        logger.warning("The system has {:d} atoms, {:d} bonds, {:d} angles, {:d} "
                       "dihedrals, and {:d} impropers. Depending upon "
                       "the size of the system, file writing may take a while and "
                       "have a large file size.".format(n_atoms, n_bonds, n_angles,
                                                        n_dihedrals, n_impropers))
        #     for k, v in filenames.items():
        #         if path.isfile(v):
        #             head, tail = path.split(v)
        #             file = "/".join((head, "1"+tail))
        #             os.system(f'mv {v} {file}')

        # Write required CHARMM input files.
        with mda.Writer(native_str(self.write_filenames["topology_file"]), **kwargs) as rtf:
            logger.info("Writing {}...".format(self.write_filenames["topology_file"]))
            rtf.write(self.universe)
        with mda.Writer(native_str(self.write_filenames["psf_file"]), **kwargs) as psf:
            logger.info("Writing {}...".format(self.write_filenames["psf_file"]))
            psf.write(self.universe)
        with mda.Writer(native_str(self.write_filenames["stream_file"]), **kwargs) as stream:
            logger.info("Writing {}...".format(self.write_filenames["stream_file"]))
            stream.write(self.universe)

        # Write an XPLOR version of the PSF
        atomtypes = topologyattrs.Atomtypes(self.universe.atoms.names)
        self.universe._topology.add_TopologyAttr(topologyattr=atomtypes)
        self.universe._generate_from_topology()
        with mda.Writer(native_str(self.write_filenames["xplor_psf_file"]), **kwargs) as psf:
            logger.info("Writing {}...".format(self.write_filenames["xplor_psf_file"]))
            psf.write(self.universe)

        # Create a new universe.
        topologies = ("names", "resids", "resnums", "resnames", "segids")
        trimmed_universe = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=self.universe.residues.n_residues,
            n_segments=self.universe.segments.n_segments,
            atom_resindex=self.universe.atoms.resindices,
            residue_segindex=self.universe.residues.segindices,
            trajectory=True)
        for _ in topologies:
            trimmed_universe.add_TopologyAttr(_)
        trimmed_universe.atoms.names = self.universe.atoms.names
        trimmed_universe.residues.resids = self.universe.residues.resids
        trimmed_universe.residues.resnums = self.universe.residues.resnums
        trimmed_universe.residues.resnames = self.universe.residues.resnames
        trimmed_universe.segments.segids = self.universe.segments.segids
        return trimmed_universe

    def tune_topology(self, fluctdiff):
        """

        Parameters
        ----------
        fluctdiff: pd.Series
            RMS fluctuation diffrence between QHA & NMA

        Returns
        -------
            Writes new topology and parameter files and return the trimmed universe
        """
        new_bonds, bonds_trimmed = self.trim_get_new_bondlist_write_prm(fluctdiff)
        bonds = topologyattrs.Bonds(self.get_atm_idx_aft_trim(new_bonds))
        self.universe._topology.add_TopologyAttr(bonds)
        self.universe._generate_from_topology()
        return self.trim_write_topology_files()



