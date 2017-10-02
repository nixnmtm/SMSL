# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from collections import OrderedDict

from future.builtins import (
    super,
    zip,
)

from .base import ModelBase
from .selection import *


class SolventIons(ModelBase):
    """Include ions within the solvent.
    """
    model = "SOLVENTIONS"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["ION"] = "name LI LIT K NA F CL BR I"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        resnames = np.unique(self.residues.resnames)
        restypes = {
            k: v
            for k, v in zip(resnames, np.arange(resnames.size)+10)
        }
        self.atoms.types = np.asarray([
            restypes[atom.resname]
            for atom in self.atoms
        ])

    def _add_bonds(self):
        pass


class BioIons(ModelBase):
    """Select ions normally found within biological systems.
    """
    model = "BIOIONS"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mapping["ions"] = "bioion"

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        resnames = np.unique(self.residues.resnames)
        restypes = {
            k: v
            for k, v in zip(resnames, np.arange(resnames.size)+20)
        }
        self.atoms.types = np.asarray([
            restypes[atom.resname]
            for atom in self.atoms
        ])

    def _add_bonds(self):
        pass


class NobleAtoms(ModelBase):
    """Select atoms column VIII of the periodic table.
    """
    model = "NOBLE"
    _mapping = OrderedDict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        kwargs["guess_bonds"] = False
        kwargs["mapping"] = self._mapping
        self._initialize(*args, **kwargs)
        resnames = np.unique(self.residues.resnames)
        restypes = {
            k: v
            for k, v in zip(resnames, np.arange(resnames.size)+40)
        }
        self.atoms.types = np.asarray([
            restypes[atom.resname]
            for atom in self.atoms
        ])

    def _add_bonds(self):
        pass
