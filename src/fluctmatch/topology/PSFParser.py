# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import time
from os import environ

import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from future.builtins import (
    dict,
    open,
)
from future.utils import (
    native_str,
)

from . import base


class PSFWriter(base.TopologyWriterBase):
    """PSF writer that implements the CHARMM PSF topology format.

    Requires the following attributes to be present:
    - ids
    - names
    - types
    - masses
    - charges
    - resids
    - resnames
    - segids
    - bonds

    .. versionchanged:: 3.0.0
       Uses numpy arrays for bond, angle, dihedral, and improper outputs.
    """
    format = "PSF"
    units = dict(time=None, length=None)
    _fmt = dict(STD="%8d %-8s %-8d %-8s %-8s %4d %14.6f%14.6f%8d",
                STD_XPLOR="{%8d %-8s %-8d %-8s %-8s %-4s %14.6f%14.6f%8d",
                STD_XPLOR_C36="%8d %-8s %-8d %-8s %-8s %-6s %14.6f%14.6f%8d",
                EXT="%10d %-8s %-8d %-8s %-8s %4d %14.6f%14.6f%8d",
                EXT_XPLOR="%10d %-8s %-8d %-8s %-8s %-4s %14.6f%14.6f%8d",
                EXT_XPLOR_C36="%10d %-8s %-8d %-8s %-8s %-6s %14.6f%14.6f%8d")

    def __init__(self, filename, extended=True, cmap=True, cheq=True, title=None, **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        extended : bool
             extended format
        cmap : bool
             include CMAP section
        cheq : bool
             include charge equilibration
        title : str
             title lines at beginning of the file
        """

        self.filename = util.filename(filename, ext="psf")
        self.extended = extended
        self.cmap = cmap
        self.cheq = cheq
        self.title = ("* Created by fluctmatch on {}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())),
                      "* User: {}".format(environ["USER"])) if title is None else title

        self.col_width = 10 if self.extended else 8
        self.sect_hdr = "{:>10d} !{}\n" if self.extended else "{:>8d} !{}"
        self.sect_hdr2 = "{:>10d}{:>10d} !{}\n" if self.extended else "{:>8d}{:>8d} !{}"
        self.sections = (("bonds", "NBOND: bonds", 8),
                         ("angles", "NTHETA: angles", 9),
                         ("dihedrals", "NPHI: dihedrals", 8),
                         ("impropers", "NIMPHI: impropers", 8),
                         ("donors", "NDON: donors", 8),
                         ("acceptors", "NACC: acceptors", 8))

    def write(self, universe):
        # type: (object) -> object
        """Write universe to PSF format.

        Parameters
        ----------
        universe : Universe
             universe to be written

        """
        self.universe = universe
        self.xplor = not np.issubdtype(universe.atoms.types.dtype, np.int)
        header = "PSF"
        if self.extended:
            header += " EXT"
        if self.cmap:
            header += " CMAP"
        if self.cheq:
            header += " CHEQ"
        if self.xplor:
            header += " XPLOR"

        self._fmtkey = "EXT" if self.extended else "STD"
        self._fmtkey += "_XPLOR" if self.xplor else ""
        if not np.issubdtype(universe.atoms.types.dtype, np.int):
            self._fmtkey += "_C36" if np.any(np.where([len(_) for _ in universe.atoms.types.astype(np.unicode)])[0] > 4) else ""


        with open(self.filename, "wb") as psffile:
            psffile.write((header + "\n").encode())
            n_title = len(self.title)
            psffile.write(self.sect_hdr.format(n_title, "NTITLE").encode())
            for _ in self.title:
                psffile.write((_ + "\n").encode())
            self._write_atoms(psffile)
            for section in self.sections:
                self._write_sec(psffile, section)
            self._write_other(psffile)


    def _write_atoms(self, psffile):
        """Write atom section in a Charmm PSF file.

        Normal (standard) and extended (EXT) PSF format are
        supported.


        CHARMM Format from ``source/psffres.src``:

        no CHEQ::
         standard format:
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6) CHEQ
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR,CHEQ
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A6,1X,2G14.6,I8)  XPLOR,c36
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A6,1X,2G14.6,I8,2G14.6)  XPLOR,c36,CHEQ
          expanded format EXT:
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6) CHEQ
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR,CHEQ
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8) XPLOR,c36
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6) XPLOR,c36,CHEQ
        """
        fmt = self._fmt[self._fmtkey]
        psffile.write(self.sect_hdr.format(self.universe.atoms.n_atoms, "NATOM").encode())
        atoms = self.universe.atoms
        lines = [atoms.ids+1, atoms.segids, atoms.resids, atoms.resnames,
                 atoms.names, atoms.types, atoms.charges, atoms.masses,
                 np.zeros_like(atoms.ids)]
        lines = np.concatenate([pd.DataFrame(_) for _ in lines], axis=1)

        if self.cheq:
            fmt += "%10.6f%18s"
            cheq = [np.zeros((atoms.n_atoms, 1)), np.full((atoms.n_atoms, 1), "-0.301140E-02")]
            cheq = np.concatenate([pd.DataFrame(_) for _ in cheq], axis=1)
            lines = np.concatenate((lines, cheq), axis=1)
        np.savetxt(psffile, lines, fmt=native_str(fmt))
        psffile.write("\n".encode())

    def _write_sec(self, psffile, section_info):
        attr, header, n_perline = section_info
        if not hasattr(self.universe, attr):
            psffile.write((self.sect_hdr.format(0, header) + "\n").encode())
            return

        values = np.asarray(getattr(self.universe, attr).to_indices()) + 1
        values = values.astype(np.unicode)
        n_rows, n_cols = values.shape
        n_values = n_perline // n_cols
        if n_rows % n_values > 0:
            n_extra = n_values - (n_rows % n_values)
            values = np.concatenate((values, np.full((n_extra, n_cols), "", dtype=np.unicode)), axis=0)
        values = values.reshape((values.shape[0] // n_values, n_perline))
        psffile.write(self.sect_hdr.format(n_rows, header).encode())
        np.savetxt(psffile, values, fmt=native_str("%{:d}s".format(self.col_width)), delimiter=native_str(""))
        psffile.write("\n".encode())

    def _write_other(self, psffile):
        n_atoms = self.universe.atoms.n_atoms
        n_cols = 8
        dn_cols = n_atoms % n_cols
        missing = n_cols - dn_cols if dn_cols > 0 else dn_cols

        # NBB
        nbb = np.full(n_atoms, "0", dtype=np.unicode)
        if dn_cols > 0:
            nbb = np.concatenate([nbb, np.empty(missing, dtype=np.unicode)], axis=0)
        nbb = nbb.reshape((nbb.size // n_cols, n_cols))
        psffile.write((self.sect_hdr.format(0, "NBB") + "\n").encode())
        np.savetxt(psffile, nbb, fmt=native_str("%{:d}s".format(self.col_width)), delimiter=native_str(""))
        psffile.write("\n".encode())

        # NGRP NST2
        psffile.write(self.sect_hdr2.format(1, 0, "NGRP NST2").encode())
        line = np.zeros(3, dtype=np.int)
        line = line.reshape((1, line.size))
        np.savetxt(psffile, line, fmt=native_str("%{:d}d".format(self.col_width)), delimiter=native_str(""))
        psffile.write("\n".encode())

        # MOLNT
        if self.cheq:
            line = np.full(n_atoms, "1", dtype=np.unicode)
            if dn_cols > 0:
                line = np.concatenate([line, np.zeros(missing, dtype=np.unicode)], axis=0)
            line = line.reshape((line.size // n_cols, n_cols))
            psffile.write(self.sect_hdr.format(1, "MOLNT").encode())
            np.savetxt(psffile, line, fmt=native_str("%{:d}s".format(self.col_width)), delimiter=native_str(""))
            psffile.write("\n".encode())
        else:
            psffile.write(self.sect_hdr.format(0, "MOLNT").encode())
            psffile.write("\n".encode())

        # NUMLP NUMLPH
        psffile.write(self.sect_hdr2.format(0, 0, "NUMLP NUMLPH").encode())
        psffile.write("\n".encode())

        # NCRTERM: cross-terms
        psffile.write(self.sect_hdr.format(0, "NCRTERM: cross-terms").encode())
        psffile.write("\n".encode())
