# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import numpy as np

from MDAnalysis.core import selection


class BackboneSelection(selection.BackboneSelection):
    token = "backbone"
    oxy_atoms = ["OXT", "OT1", "OT2"]

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.bb_atoms, self.oxy_atoms]))
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class HBackboneSelection(BackboneSelection):
    token = "hbackbone"
    hbb_atoms = np.array(["H", "HN", "H1", "H2", "H3", "HT1", "HT2", "HT3", "HA", "HA1", "HA2", "1HA", "2HA"])

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.bb_atoms, self.oxy_atoms, self.hbb_atoms]))
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class CalphaSelection(selection.ProteinSelection):
    token = "calpha"
    calpha = np.array(["CA"])

    def apply(self, group):
        mask = np.in1d(group.names, self.calpha)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class HCalphaSelection(CalphaSelection):
    token = "hcalpha"
    hcalpha = np.array(["HA", "HA1", "HA2", "1HA", "2HA"])

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.calpha, self.hcalpha]))
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class CbetaSelection(selection.ProteinSelection):
    token = "cbeta"
    cbeta = np.array(["CB"])

    def apply(self, group):
        mask = np.in1d(group.names, self.cbeta)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class AmineSelection(selection.ProteinSelection):
    token = "amine"
    amine = np.array(["N", "HN", "H", "H1", "H2", "H3", "HT1", "HT2", "HT3"])

    def apply(self, group):
        mask = np.in1d(group.names, self.amine)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class CarboxylSelection(selection.ProteinSelection):
    token = "carboxyl"
    carboxyl = np.array(["C", "O", "OXT", "OT1", "OT2"])

    def apply(self, group):
        mask = np.in1d(group.names, self.carboxyl)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class HSidechainSelection(HBackboneSelection):
    token = "hsidechain"

    def apply(self, group):
        mask = np.in1d(group.names, np.concatenate([self.bb_atoms, self.oxy_atoms, self.hbb_atoms]), invert=True)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class BioIonSelection(selection.Selection):
    token = "bioion"
    ion_atoms = np.array(["MG", "CAL", "MN", "FE", "CU", "ZN", "AG"])

    def __init__(self, parser, tokens):
        pass

    def apply(self, group):
        mask = np.in1d(group.resnames, self.ion_atoms)
        return group[mask].unique
