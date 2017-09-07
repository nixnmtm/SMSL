# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.builtins import (
    super,
    zip,
)
from collections import OrderedDict

import numpy as np

from . import universe


class SolventIons(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["ION"] = "name LI K NA F CL BR"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)
        resnames = np.unique(self.residues.names)
        restypes = {k: v for k, v in zip (resnames, np.arange(resnames.size)+10)}
        self.atoms.types = np.asarray([restypes[atom.resname] for atom in self.atoms])

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _add_bonds(self):
        pass


class OtherIons(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["ions"] = "name MG CAL MN FE CU ZN AG"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)
        resnames = np.unique(self.residues.names)
        restypes = {k: v for k, v in zip (resnames, np.arange(resnames.size)+20)}
        self.atoms.types = np.asarray([restypes[atom.resname] for atom in self.atoms])

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _add_bonds(self):
        pass


class NobleAtoms(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)
        resnames = np.unique(self.residues.names)
        restypes = {k: v for k, v in zip (resnames, np.arange(resnames.size)+40)}
        self.atoms.types = np.asarray([restypes[atom.resname] for atom in self.atoms])

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _add_bonds(self):
        pass
