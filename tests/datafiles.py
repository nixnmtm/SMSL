# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# fluctmatch --- https://github.com/tclick/python-fluctmatch
# Copyright (c) 2013-2017 The fluctmatch Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the New BSD license.
#
# Please cite your use of fluctmatch in published work:
#
# Timothy H. Click, Nixon Raj, and Jhih-Wei Chu.
# Calculation of Enzyme Fluctuograms from All-Atom Molecular Dynamics
# Simulation. Meth Enzymology. 578 (2016), 327-342,
# doi:10.1016/bs.mie.2016.05.024.
#
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from os import path

__all__ = [
    "PDB",  # PDB
    "PDB_prot",
    "PDB_dna",
    "TIP3P",
    "TIP4P",
    "IONS",
    "DMA",
    "TPR",
    "XTC",  # Gromacs
    "NCSC",
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
NCSC = resource_filename(__name__, path.join("data", "ncsc.pdb"))
