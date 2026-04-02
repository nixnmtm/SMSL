# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# SMSL - https://github.com/nixnmtm/SMSL
#
#

from pathlib import Path

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

_DATA = Path(__file__).resolve().parent / "data"

PDB = str(_DATA / "trex1.pdb")
PDB_prot = str(_DATA / "protein.pdb")
PDB_dna = str(_DATA / "dna.pdb")
TPR = str(_DATA / "trex1.tpr")
XTC = str(_DATA / "trex1.xtc")
TIP3P = str(_DATA / "spc216.gro")
TIP4P = str(_DATA / "tip4p.gro")
IONS = str(_DATA / "ions.pdb")
DMA = str(_DATA / "dma.gro")
NCSC = str(_DATA / "ncsc.pdb")