# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import logging
import time
from os import environ

import numpy as np
import pandas as pd
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Charges,
    Resids,
    Resnums,
    Resnames,
    Segids,
)
from MDAnalysis.lib import util
from MDAnalysis.topology import PSFParser
from MDAnalysis.topology.base import change_squash
from future.builtins import (
    dict,
    range,
)
from future.utils import (
    native_str,
)

from MDAnalysis.core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Charges,
    Resids,
    Resnums,
    Resnames,
    Segids,
    Bonds,
    Angles,
    Dihedrals,
    Impropers
)
from MDAnalysis.core.topology import Topology
from fluctmatch.topology import base

logger = logging.getLogger("MDAnalysis.topology.PSF")


# Changed the segid squash_by to change_squash to prevent segment ID sorting.
class PSF36Parser(PSFParser.PSFParser):
    """Read topology information from a CHARMM/NAMD/XPLOR PSF_ file.

    Creates a Topology with the following Attributes:
    - ids
    - names
    - types
    - masses
    - charges
    - resids
    - resnames
    - segids
    - bonds
    - angles
    - dihedrals
    - impropers

    .. _PSF: http://www.charmm.org/documentation/c35b1/struct.html
    """
    format = 'PSF'

    def parse(self):
        """Parse PSF file into Topology

        Returns
        -------
        MDAnalysis *Topology* object
        """
        # Open and check psf validity
        with util.openany(self.filename, 'r') as psffile:
            header = next(psffile)
            if not header.startswith("PSF"):
                err = ("{0} is not valid PSF file (header = {1})"
                       "".format(self.filename, header))
                logger.error(err)
                raise ValueError(err)
            header_flags = header[3:].split()

            if "NAMD" in header_flags:
                self._format = "NAMD"        # NAMD/VMD
            elif "EXT" in header_flags:
                self._format = "EXTENDED"    # CHARMM
            else:
                self._format = "STANDARD"    # CHARMM
            if "XPLOR" in header_flags:
                self._format += "_XPLOR"

            next(psffile)
            title = next(psffile).split()
            if not (title[1] == "!NTITLE"):
                err = "{0} is not a valid PSF file".format(psffile.name)
                logger.error(err)
                raise ValueError(err)
            # psfremarks = [psffile.next() for i in range(int(title[0]))]
            for _ in range(int(title[0])):
                next(psffile)
            logger.debug("PSF file {0}: format {1}"
                         "".format(psffile.name, self._format))

            # Atoms first and mandatory
            top = self._parse_sec(
                psffile, ('NATOM', 1, 1, self._parseatoms))
            # Then possibly other sections
            sections = (
                #("atoms", ("NATOM", 1, 1, self._parseatoms)),
                (Bonds, ("NBOND", 2, 4, self._parsesection)),
                (Angles, ("NTHETA", 3, 3, self._parsesection)),
                (Dihedrals, ("NPHI", 4, 2, self._parsesection)),
                (Impropers, ("NIMPHI", 4, 2, self._parsesection)),
                #("donors", ("NDON", 2, 4, self._parsesection)),
                #("acceptors", ("NACC", 2, 4, self._parsesection))
            )

            try:
                for attr, info in sections:
                    next(psffile)
                    top.add_TopologyAttr(
                        attr(self._parse_sec(psffile, info)))
            except StopIteration:
                # Reached the end of the file before we expected
                pass

        return top

    def _parseatoms(self, lines, atoms_per, numlines):
        """Parses atom section in a Charmm PSF file.

        Normal (standard) and extended (EXT) PSF format are
        supported. CHEQ is supported in the sense that CHEQ data is simply
        ignored.


        CHARMM Format from ``source/psffres.src``:

        CHEQ::
          II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)

          standard format:
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR
          expanded format EXT:
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR

        no CHEQ::
          II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)

         standard format:
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR
          expanded format EXT:
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR

        NAMD PSF

        space separated, see release notes for VMD 1.9.1, psfplugin at
        http://www.ks.uiuc.edu/Research/vmd/current/devel.html :

        psfplugin: Added more logic to the PSF plugin to determine cases where the
        CHARMM "EXTended" PSF format cannot accomodate long atom types, and we add
        a "NAMD" keyword to the PSF file flags line at the top of the file. Upon
        reading, if we detect the "NAMD" flag there, we know that it is possible
        to parse the file correctly using a simple space-delimited scanf() format
        string, and we use that strategy rather than holding to the inflexible
        column-based fields that are a necessity for compatibility with CHARMM,
        CNS, X-PLOR, and other formats. NAMD and the psfgen plugin already assume
        this sort of space-delimited formatting, but that's because they aren't
        expected to parse the PSF variants associated with the other programs. For
        the VMD PSF plugin, having the "NAMD" tag in the flags line makes it
        absolutely clear that we're dealing with a NAMD-specific file so we can
        take the same approach.

        """
        # how to partition the line into the individual atom components
        atom_parsers = dict(
            STANDARD="I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2F14.6,I8",
            STANDARD_XPLOR="'(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2F14.6,I8",
            EXTENDED="I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2F14.6,I8",
            EXTENDED_XPLOR="I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2F14.6,I8",
            NAMD="I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2F14.6,I8",
        )
        atom_parser = util.FORTRANReader(atom_parsers[self._format])
        # once partitioned, assigned each component the correct type
        set_type = lambda x: (int(x[0]) - 1, x[1] or "SYSTEM", int(x[2]), x[3],
                              x[4], x[5], float(x[6]), float(x[7]))

        # Oli: I don't think that this is the correct OUTPUT format:
        #   psf_atom_format = "   %5d %4s %4d %4s %-4s %-4s %10.6f      %7.4f%s\n"
        # It should be rather something like:
        #   psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
        #                     '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'

        # source/psfres.src (CHEQ and now can be used for CHEQ EXTended), see comments above
        #   II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
        #  (I8,1X,A4, 1X,A4,  1X,A4,  1X,A4,  1X,I4,  1X,2G14.6,     I8,   2G14.6)
        #   0:8   9:13   14:18   19:23   24:28   29:33   34:48 48:62 62:70 70:84 84:98

        # Allocate arrays
        atomids = np.zeros(numlines, dtype=np.int32)
        segids = np.zeros(numlines, dtype=object)
        resids = np.zeros(numlines, dtype=np.int32)
        resnames = np.zeros(numlines, dtype=object)
        atomnames = np.zeros(numlines, dtype=object)
        atomtypes = np.zeros(numlines, dtype=object)
        charges = np.zeros(numlines, dtype=np.float32)
        masses = np.zeros(numlines, dtype=np.float64)

        for i in range(numlines):
            try:
                line = lines()
            except StopIteration:
                err = ("{0} is not valid PSF file"
                       "".format(self.filename))
                logger.error(err)
                raise ValueError(err)
            try:
                vals = atom_parser.read(line)
            except ValueError:
                # last ditch attempt: this *might* be a NAMD/VMD
                # space-separated "PSF" file from VMD version < 1.9.1
                try:
                    atom_parser = util.FORTRANReader(atom_parsers['NAMD'])
                    vals = atom_parser.read(line)
                    logger.warn("Guessing that this is actually a"
                                " NAMD-type PSF file..."
                                " continuing with fingers crossed!")
                    logger.debug("First NAMD-type line: {0}: {1}"
                                 "".format(i, line.rstrip()))
                except ValueError:
                    atom_parser = util.FORTRANReader(atom_parsers[self._format].replace("A6", "A4"))
                    vals = atom_parser.read(line)
                    logger.warn("Guessing that this is actually a"
                                " pre CHARMM36 PSF file..."
                                " continuing with fingers crossed!")
                    logger.debug("First NAMD-type line: {0}: {1}"
                                 "".format(i, line.rstrip()))

            atomids[i] = vals[0]
            segids[i] = vals[1] if vals[1] else "SYSTEM"
            resids[i] = vals[2]
            resnames[i] = vals[3]
            atomnames[i] = vals[4]
            atomtypes[i] = vals[5]
            charges[i] = vals[6]
            masses[i] = vals[7]

        # Atom
        atomids = Atomids(atomids - 1)
        atomnames = Atomnames(atomnames)
        atomtypes = Atomtypes(atomtypes)
        charges = Charges(charges)
        masses = Masses(masses)

        # Residue
        # resids, resnames
        residx, (new_resids, new_resnames, perres_segids) = change_squash(
            (resids, resnames, segids),
            (resids, resnames, segids))
        # transform from atom:Rid to atom:Rix
        residueids = Resids(new_resids)
        residuenums = Resnums(new_resids.copy())
        residuenames = Resnames(new_resnames)

        # Segment
        segidx, (perseg_segids,) = change_squash((perres_segids,), (perres_segids,))
        segids = Segids(perseg_segids)

        top = Topology(len(atomids), len(new_resids), len(segids),
                       attrs=[atomids, atomnames, atomtypes,
                              charges, masses,
                              residueids, residuenums, residuenames,
                              segids],
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top


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

    Parameters
    ----------
    filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
         name of the output file or a stream
    n_atoms : int, optional
        The number of atoms in the output trajectory.
    extended
         extended format
    cmap
         include CMAP section
    cheq
         include charge equilibration
    title
         title lines at beginning of the file
    charmm_version
        Version of CHARMM for formatting (default: 41)
    """
    format = "PSF"
    units = dict(time=None, length=None)
    _fmt = dict(
        STD="%8d %-4s %-4d %-4s %-4s %4d %14.6f%14.6f%8d",
        STD_XPLOR="{%8d %4s %-4d %-4s %-4s %-4s %14.6f%14.6f%8d",
        STD_XPLOR_C35="%4d %-4s %-4d %-4s %-4s %-4s %14.6f%14.6f%8d",
        EXT="%10d %-8s %8d %-8s %-8s %4d %14.6f%14.6f%8d",
        EXT_XPLOR="%10d %-8s %-8d %-8s %-8s %-6s %14.6f%14.6f%8d",
        EXT_XPLOR_C35="%10d %-8s %-8d %-8s %-8s %-4s %14.6f%14.6f%8d"
    )

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext="psf")
        self._extended = kwargs.get("extended", True)
        self._cmap = kwargs.get("cmap", True)
        self._cheq = kwargs.get("cheq", True)
        self._version = kwargs.get("charmm_version", 41)
        self._universe = None
        self._fmtkey = "EXT" if self._extended else "STD"

        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        user = environ["USER"]
        self._title = kwargs.get(
            "title",
            (
                "* Created by fluctmatch on {date}".format(date=date),
                "* User: {user}".format(user=user),
            )
        )
        if not util.iterable(self._title):
            self._title = util.asiterable(self._title)

        self.col_width = 10 if self._extended else 8
        self.sect_hdr = "{:>10d} !{}" if self._extended else "{:>8d} !{}"
        self.sect_hdr2 = "{:>10d}{:>10d} !{}" if self._extended else "{:>8d}{:>8d} !{}"
        self.sections = (("bonds", "NBOND: bonds", 8),
                         ("angles", "NTHETA: angles", 9),
                         ("dihedrals", "NPHI: dihedrals", 8),
                         ("impropers", "NIMPHI: impropers", 8),
                         ("donors", "NDON: donors", 8),
                         ("acceptors", "NACC: acceptors", 8))

    def write(self, universe):
        """Write universe to PSF format.

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or atomgroup with bond definitions.
        """
        self._universe = universe
        xplor = not np.issubdtype(universe.atoms.types.dtype, np.int)

        header = "PSF"
        if self._extended:
            header += " EXT"
        if self._cheq:
            header += " CHEQ"
        if xplor:
            header += " XPLOR"
        if self._cmap:
            header += " CMAP"
        header += "\n"

        if xplor:
            self._fmtkey += "_XPLOR"
            if self._version < 36:
                self._fmtkey += "_C35"

        with util.openany(self.filename, "w") as psffile:
            print(header, file=psffile)
            n_title = len(self._title)
            print(self.sect_hdr.format(n_title, "NTITLE"), file=psffile)
            for _ in self._title:
                print(_, file=psffile)
            print(file=psffile)
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
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A6,1X,2G14.6,I8)  XPLOR
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A6,1X,2G14.6,I8,2G14.6)  XPLOR,CHEQ
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR,c35
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR,c35,CHEQ
          expanded format EXT:
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6) CHEQ
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8) XPLOR
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6) XPLOR,CHEQ
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR,c35
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR,c35,CHEQ
        """
        fmt = self._fmt[self._fmtkey]
        print(self.sect_hdr.format(self._universe.atoms.n_atoms, "NATOM"), file=psffile)
        atoms = self._universe.atoms
        lines = (
            np.arange(atoms.n_atoms) + 1,
            atoms.segids,
            atoms.resids,
            atoms.resnames,
            atoms.names,
            atoms.types,
            atoms.charges,
            atoms.masses,
            np.zeros_like(atoms.ids)
        )
        lines = pd.concat([pd.DataFrame(_) for _ in lines], axis=1)

        if self._cheq:
            fmt += "%10.6f%18s"
            cheq = (
                np.zeros_like(atoms.masses),
                np.full_like(atoms.names.astype(np.object), "-0.301140E-02")
            )
            cheq = pd.concat([pd.DataFrame(_) for _ in cheq], axis=1)
            lines = pd.concat([lines, cheq], axis=1)
        np.savetxt(psffile, lines, fmt=native_str(fmt))
        print(file=psffile)

    def _write_sec(self, psffile, section_info):
        attr, header, n_perline = section_info
        if not hasattr(self._universe, attr):
            print(self.sect_hdr.format(0, header), file=psffile)
            print("\n", file=psffile)
            return
        if len(getattr(self._universe, attr).to_indices()) < 2:
            print(self.sect_hdr.format(0, header), file=psffile)
            print("\n", file=psffile)
            return

        values = np.asarray(getattr(self._universe, attr).to_indices()) + 1
        values = values.astype(np.object)
        n_rows, n_cols = values.shape
        n_values = n_perline // n_cols
        if n_rows % n_values > 0:
            n_extra = n_values - (n_rows % n_values)
            values = np.concatenate((values, np.full((n_extra, n_cols), "", dtype=np.object)), axis=0)
        values = values.reshape((values.shape[0] // n_values, n_perline))
        print(self.sect_hdr.format(n_rows, header), file=psffile)
        np.savetxt(psffile, values, fmt=native_str("%{:d}s".format(self.col_width)), delimiter=native_str(""))
        print(file=psffile)

    def _write_other(self, psffile):
        n_atoms = self._universe.atoms.n_atoms
        n_cols = 8
        dn_cols = n_atoms % n_cols
        missing = n_cols - dn_cols if dn_cols > 0 else dn_cols

        # NNB
        nnb = np.full(n_atoms, "0", dtype=np.object)
        if dn_cols > 0:
            nnb = np.concatenate([nnb, np.empty(missing, dtype=np.object)], axis=0)
        nnb = nnb.reshape((nnb.size // n_cols, n_cols))
        print(self.sect_hdr.format(0, "NNB") + "\n", file=psffile)
        np.savetxt(psffile, nnb, fmt=native_str("%{:d}s".format(self.col_width)), delimiter=native_str(""))
        print(file=psffile)

        # NGRP NST2
        print(self.sect_hdr2.format(1, 0, "NGRP NST2"), file=psffile)
        line = np.zeros(3, dtype=np.int)
        line = line.reshape((1, line.size))
        np.savetxt(psffile, line, fmt=native_str("%{:d}d".format(self.col_width)), delimiter=native_str(""))
        print(file=psffile)

        # MOLNT
        if self._cheq:
            line = np.full(n_atoms, "1", dtype=np.object)
            if dn_cols > 0:
                line = np.concatenate([line, np.zeros(missing, dtype=np.object)], axis=0)
            line = line.reshape((line.size // n_cols, n_cols))
            print(self.sect_hdr.format(1, "MOLNT"), file=psffile)
            np.savetxt(psffile, line, fmt=native_str("%{:d}s".format(self.col_width)), delimiter=native_str(""))
            print(file=psffile)
        else:
            print(self.sect_hdr.format(0, "MOLNT"), file=psffile)
            print(file=psffile)

        # NUMLP NUMLPH
        print(self.sect_hdr2.format(0, 0, "NUMLP NUMLPH"), file=psffile)
        print("\n", file=psffile)

        # NCRTERM: cross-terms
        print(self.sect_hdr.format(0, "NCRTERM: cross-terms"), file=psffile)
        print("\n", file=psffile)
