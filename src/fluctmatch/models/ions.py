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


from . import universe
from .selection import *


class SolventIons(universe._Universe):
    """Include ions within the solvent.
    """
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["ION"] = "name LI LIT K NA F CL BR I"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        resnames = np.unique(self.residues.names)
        restypes = {
            k: v
            for k, v in zip(resnames, np.arange(resnames.size)+10)
        }
        self.atoms.types = np.asarray([
            restypes[atom.resname]
            for atom in self.atoms
        ])

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _add_bonds(self):
        pass


class BioIons(universe._Universe):
    """Select ions normally found within biological systems.
    """
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["ions"] = "bioion"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        resnames = np.unique(self.residues.names)
        restypes = {
            k: v
            for k, v in zip(resnames, np.arange(resnames.size)+20)
        }
        self.atoms.types = np.asarray([
            restypes[atom.resname]
            for atom in self.atoms
        ])

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _add_bonds(self):
        pass


class NobleAtoms(universe._Universe):
    """Select atoms column VIII of the periodic table.
    """
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        resnames = np.unique(self.residues.names)
        restypes = {
            k: v
            for k, v in zip(resnames, np.arange(resnames.size)+40)
        }
        self.atoms.types = np.asarray([
            restypes[atom.resname]
            for atom in self.atoms
        ])

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _add_bonds(self):
        pass
