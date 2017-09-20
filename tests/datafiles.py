# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from os import path

__all__ = [
    "PDB",          # PDB
    "PDB_prot", "PDB_dna",
    "TIP3P", "TIP4P", "IONS", "DMA",
    "TPR", "XTC",   # Gromacs
]

from pkg_resources import resource_filename

PDB = resource_filename(__name__, path.join("data", "trex1.pdb"))
PDB_prot = resource_filename(__name__, path.join("data", "protein.pdb"))
PDB_dna = resource_filename(__name__, path.join("data", "dna.pdb"))
TPR = resource_filename(__name__, path.join("data", "trex1.tpr"))
XTC = resource_filename(__name__, path.join("data", "trex1.xtc"))
TIP3P = resource_filename(__name__, path.join("data", "spc216.gro"))
TIP4P = resource_filename(__name__, path.join("data", "tip4p.gro"))
IONS = resource_filename(__name__, path.join("data", "ions.pdb"))
DMA = resource_filename(__name__, path.join("data", "dma.gro"))
