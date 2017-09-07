# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.utils import viewvalues
from future.builtins import (
    super,
    zip,
)

from collections import OrderedDict

from MDAnalysis.core import topologyattrs
from . import universe
from .selection import *


class Nucleic4(universe._Universe):
    _mapping = OrderedDict()

    def __init__(self, topfn, crdfn, com=True, extended=True, xplor=True, **kwargs):
        super().__init__(topfn, crdfn, com, extended, xplor, **kwargs)
        self._mapping["P"] = "(nucleic or resname OXG) and name P O1P O2P O5' H5T"
        self._mapping["C4'"] = "(nucleic or resname OXG) and name C4' O4' H4' C5' H5'*"
        self._mapping["C1'"] = "(nucleic or resname OXG) and name C1' H1' C2' O2' H2'*"
        self._mapping["C5"] = "(nucleic or resname OXG) and name N4 C4 C5"
        _sugar = "(nucleic or resname OXG) and name C3' O3' H3'* H3T"
        _nucl = "(nucleic or resname OXG) and not ({} or {} or {} or {} or {})".format(_sugar, *viewvalues(self._mapping))

        kwargs["mapping"] = self._mapping
        self._initialize(topfn, crdfn, **kwargs)

        # Update the masses and charges
        sug_atu = self.atu.select_atoms(_sugar).split("residue")
        nucl_atu = self.atu.select_atoms(_nucl).split("residue")
        self.atoms.select_atoms("name C4"").masses += 0.5 * np.array([_.total_mass() for _ in sug_atu])
        self.atoms.select_atoms("name C1"").masses += 0.5 * np.array([_.total_mass() for _ in sug_atu])
        self.atoms.select_atoms("name C5").masses = np.array([_.total_mass() for _ in nucl_atu])

        if hasattr(sug_atu, "total_charge") and hasattr(nucl_atu, "total_charge"):
            self.atoms.select_atoms("name C1"").charges += 0.5 * np.array([_.total_charge() for _ in sug_atu])
            self.atoms.select_atoms("name C4"").charges += 0.5 * np.array([_.total_charge() for _ in sug_atu])
            self.atoms.select_atoms("name C5").charges = np.array([_.total_charge() for _ in nucl_atu])

    def _add_bonds(self):
        bonds = []
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("name P").ix,
                                                            s.atoms.select_atoms("name C4"").ix)])
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("name C4"").ix,
                                                            s.atoms.select_atoms("name C1"").ix)])
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("name C1"").ix,
                                                            s.atoms.select_atoms("name C5").ix)])
        bonds.extend([_ for s in self.segments for _ in zip(s.atoms.select_atoms("name C4"").ix[:-1],
                                                            s.atoms.select_atoms("name P").ix[1:])])
        self._topology.add_TopologyAttr(topologyattrs.Bonds(bonds))
        self._generate_from_topology()
