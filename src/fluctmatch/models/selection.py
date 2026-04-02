# -*- coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# fluctmatch --- https://github.com/nixnmtm/SMSL
# Copyright (c) 2013-2017 The fluctmatch Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the New BSD license.
#
# Please cite your use of fluctmatch in published work:

import numpy as np
from MDAnalysis.core import selection


def _arr(x):
    return np.atleast_1d(np.asarray(x, dtype=object))


def _cat(*x):
    return np.concatenate([_arr(i) for i in x])


# ======================
# BASIC SELECTIONS
# ======================

class BioIonSelection(selection.Selection):
    token = "bioion"
    ion_atoms = _arr(["MG", "CAL", "MN", "FE", "CU", "ZN", "AG"])

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        names = _arr(group.names)
        return group[np.isin(names, self.ion_atoms)].unique


class WaterSelection(selection.Selection):
    token = "water"
    water_atoms = _arr(["OW", "HW1", "HW2", "MW"])

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        names = _arr(group.names)
        return group[np.isin(names, self.water_atoms)].unique


# ======================
# PROTEIN
# ======================

class CalphaSelection(selection.Selection):
    token = "calpha"
    calpha = _arr(["CA"])

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        names = _arr(group.names)
        return group[np.isin(names, self.calpha)].unique


class HCalphaSelection(CalphaSelection):
    token = "hcalpha"
    hcalpha = _arr(["HA", "HA1", "HA2", "1HA", "2HA"])

    def apply(self, group):
        names = _arr(group.names)
        mask = np.isin(names, _cat(self.calpha, self.hcalpha))
        return group[mask].unique

class BackboneSelection(selection.Selection):
    token = "backbone"
    bb_atoms = _arr(["N", "CA", "C", "O"])
    oxy_atoms = _arr(["OXT", "OT1", "OT2"])

    def apply(self, group):
        names = _arr(group.names)
        mask = np.isin(names, _cat(self.bb_atoms, self.oxy_atoms))
        return group[mask].unique


class HBackboneSelection(BackboneSelection):
    token = "hbackbone"
    hbb_atoms = _arr([
        "H", "HN", "H1", "H2", "H3", "HT1", "HT2", "HT3",
        "HA", "HA1", "HA2", "1HA", "2HA"
    ])

    def apply(self, group):
        names = _arr(group.names)
        mask = np.isin(names, _cat(self.bb_atoms, self.oxy_atoms, self.hbb_atoms))
        return group[mask].unique

class CbetaSelection(selection.ProteinSelection):
    token = "cbeta"
    cbeta = _arr(["CB"])

    def apply(self, group):
        names = _arr(group.names)
        mask = np.isin(names, self.cbeta)
        return group[mask].unique


class HSidechainSelection(HBackboneSelection):
    token = "hsidechain"

    def apply(self, group):
        names = _arr(group.names)

        mask = np.isin(
            names,
            _cat(
                self.bb_atoms, self.oxy_atoms, self.hbb_atoms,
                BioIonSelection.ion_atoms
            ),
            invert=True
        )

        return group[mask].unique

class AmineSelection(selection.Selection):
    token = "amine"
    amine = _arr(["N", "HN", "H", "H1", "H2", "H3", "HT1", "HT2", "HT3"])

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        names = _arr(group.names)
        return group[np.isin(names, self.amine)].unique


class CarboxylSelection(selection.Selection):
    token = "carboxyl"
    carboxyl = _arr(["C", "O", "OXT", "OT1", "OT2"])

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

    def apply(self, group):
        names = _arr(group.names)
        return group[np.isin(names, self.carboxyl)].unique


# ======================
# NUCLEIC
# ======================

class AdditionalNucleicSelection(selection.NucleicSelection):
    token = "nucleic"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.nucl_res = _cat(self.nucl_res, ["OXG", "HPX", "DC35"])

    def apply(self, group):
        return group[np.isin(_arr(group.resnames), self.nucl_res)].unique


class HNucleicSugarSelection(AdditionalNucleicSelection, selection.NucleicSugarSelection):
    token = "hnucleicsugar"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.sug_atoms = _cat(self.sug_atoms, [
            "H1'", "O1'", "O2'", "H2'", "H2''",
            "O3'", "H3'", "H3T", "H4'"
        ])

    def apply(self, group):
        names = _arr(group.names)
        resnames = _arr(group.resnames)

        mask = np.isin(names, self.sug_atoms)
        mask &= np.isin(resnames, self.nucl_res)
        return group[mask].unique


class HBaseSelection(AdditionalNucleicSelection, selection.BaseSelection):
    token = "hnucleicbase"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.base_atoms = _cat(self.base_atoms, [
            "O8", "H8", "H21", "H22", "H2",
            "O6", "H6", "H61", "H62",
            "H41", "H42", "H5", "H51", "H52", "H53",
            "H3", "H7"
        ])

    def apply(self, group):
        names = _arr(group.names)
        resnames = _arr(group.resnames)

        mask = np.isin(names, self.base_atoms)
        mask &= np.isin(resnames, self.nucl_res)
        return group[mask].unique


class NucleicPhosphateSelection(AdditionalNucleicSelection):
    token = "nucleicphosphate"
    phos_atoms = _arr(["P", "O1P", "O2P", "O5'", "C5'", "H5'", "H5''", "O5T", "H5T"])

    def apply(self, group):
        names = _arr(group.names)
        resnames = _arr(group.resnames)

        mask = np.isin(names, self.phos_atoms)
        mask &= np.isin(resnames, self.nucl_res)
        return group[mask].unique


class NucleicC2Selection(AdditionalNucleicSelection):
    token = "sugarC2"
    c3_atoms = _arr(["C1'", "H1'", "C2'", "O2'", "H2'", "H2''"])

    def apply(self, group):
        names = _arr(group.names)
        resnames = _arr(group.resnames)

        mask = np.isin(names, self.c3_atoms)
        mask &= np.isin(resnames, self.nucl_res)
        return group[mask].unique


class NucleicC4Selection(AdditionalNucleicSelection):
    token = "sugarC4"
    c3_atoms = _arr(["C3'", "O3'", "H3'", "H3T", "C4'", "O4'", "H4'"])

    def apply(self, group):
        names = _arr(group.names)
        resnames = _arr(group.resnames)

        mask = np.isin(names, self.c3_atoms)
        mask &= np.isin(resnames, self.nucl_res)
        return group[mask].unique


class BaseCenterSelection(AdditionalNucleicSelection):
    token = "nucleiccenter"
    center_atoms = _arr(["C4", "C5"])

    def apply(self, group):
        names = _arr(group.names)
        resnames = _arr(group.resnames)

        mask = np.isin(names, self.center_atoms)
        mask &= np.isin(resnames, self.nucl_res)
        return group[mask].unique