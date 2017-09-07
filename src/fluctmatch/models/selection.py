# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import numpy as np
from MDAnalysis.core import selection
from future.builtins import (
    super,
)


class BackboneSelection(selection.BackboneSelection):
    """Contains all heavy atoms within a protein backbone including the terminal carboxyl oxygens.
    """
    token = "backbone"
    oxy_atoms = ["OXT", "OT1", "OT2"]

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.bb_atoms, self.oxy_atoms]))
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class HBackboneSelection(BackboneSelection):
    """Includes all atoms found within a protein backbone including hydrogens.
    """
    token = "hbackbone"
    hbb_atoms = np.array(["H", "HN", "H1", "H2", "H3", "HT1", "HT2", "HT3", "HA", "HA1", "HA2", "1HA", "2HA"])

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.bb_atoms, self.oxy_atoms, self.hbb_atoms]))
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class CalphaSelection(selection.ProteinSelection):
    """Contains only the alpha-carbon of a protein.
    """
    token = "calpha"
    calpha = np.array(["CA"])

    def apply(self, group):
        mask = np.in1d(group.names, self.calpha)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class HCalphaSelection(CalphaSelection):
    """Contains the alpha-carbon and alpha-hydrogens of a protein.
    """
    token = "hcalpha"
    hcalpha = np.array(["HA", "HA1", "HA2", "1HA", "2HA"])

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.calpha, self.hcalpha]))
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class CbetaSelection(selection.ProteinSelection):
    """Contains only the beta-carbon of a protein.
    """
    token = "cbeta"
    cbeta = np.array(["CB"])

    def apply(self, group):
        mask = np.in1d(group.names, self.cbeta)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class AmineSelection(selection.ProteinSelection):
    """Contains atoms within the amine group of a protein.
    """
    token = "amine"
    amine = np.array(["N", "HN", "H", "H1", "H2", "H3", "HT1", "HT2", "HT3"])

    def apply(self, group):
        mask = np.in1d(group.names, self.amine)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class CarboxylSelection(selection.ProteinSelection):
    """Contains atoms within the carboxyl group of a protein.
    """
    token = "carboxyl"
    carboxyl = np.array(["C", "O", "OXT", "OT1", "OT2"])

    def apply(self, group):
        mask = np.in1d(group.names, self.carboxyl)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class HSidechainSelection(HBackboneSelection):
    """Includes hydrogens on the protein sidechain.
    """
    token = "hsidechain"

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.bb_atoms, self.oxy_atoms, self.hbb_atoms]), invert=True)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class BioIonSelection(selection.Selection):
    """Contains atoms commonly found in proteins.
    """
    token = "bioion"
    ion_atoms = np.array(["MG", "CAL", "MN", "FE", "CU", "ZN", "AG"])

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        mask = np.in1d(group.resnames, self.ion_atoms)
        return group[mask].unique


class AdditionalNucleicSelection(selection.NucleicSelection):
    """Contains additional nucleic acid residues."""
    nucl_res = np.concatenate(
        (nucl_res,
         ["OXG", "HPX"]),
        axis=0
    )

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        mask = np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class HNucleicSugarSelection(AdditionalNucleicSelection, selection.NucleicSugarSelection):
    """Contains the additional atoms definitions for the sugar.
    """
    token = "hnucleicsugar"
    sug_atoms = np.concatenate(
        (sug_atoms,
         np.array(["H1'", "O1'", "O2'", "H2'", "H2''", "O3'", "H3'", "H3T", "H4'"])),
        axis=0
    )

    def apply(self, group):
        mask = np.in1d(group.names, self.sug_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class NucleicPhosphateSelection(AdditionalNucleicSelection):
    """Contains the nucleic phosphate group including the C5'.
    """
    token = "nucleicphosphate"
    phos_atoms = np.array(["P", "O1P", "O2P", "O5'", "C5'", "H5'", "H5''", "H5T"])

    def apply(self, group):
        mask = np.in1d(group.names, self.phos_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class NucleicC3Selection(AdditionalNucleicSelection):
    """Contains the definition for the C3' region.
    """
    token = "C3'"
    c3_atoms = np.array(["C3'", "O3'", "H3'", "H3T", "C2'", "H2'", "O2'", "H2''"])

    def apply(self, group):
        mask = np.in1d(group.names, self.c3_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class NucleicC43Selection(AdditionalNucleicSelection):
    """Contains the definition for the C4' region.
    """
    token = "C4'"
    c3_atoms = np.array(["C4'", "O4'", "H4'", "C1'", "H1'", "C2'", "H2'", "O2'", "H2''"])

    def apply(self, group):
        mask = np.in1d(group.names, self.c3_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class HBaseSelection(AdditionalNucleicSelection, selection.BaseSelection):
    """Contains additional atoms on the base region of the nucleic acids.
    """
    token = "hnucleicbase"
    base_atoms = np.concatenate((
        base_atoms,
        ["O8", "H8", "H21", "H22", "H2", "O6", "H6", "H61", "H62", "H41", "H42", "H5",
         "H51", "H52", "H53", "H3"]
    ), axis=0)

    def apply(self, group):
        mask = np.in1d(group.names, self.base_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class BaseCenterSelection(AdditionalNucleicSelection):
    """Contains the central atoms (C4 and C5) on the base of the nuleic acid.
    """
    token = "nucleiccenter"
    center_atoms = np.array(["C4", "C5"])

    def apply(self, group):
        mask = np.in1d(group.names, self.center_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique
