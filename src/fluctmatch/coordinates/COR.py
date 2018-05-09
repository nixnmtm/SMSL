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
"""COR structure files in fluctmatch --- :mod:`MDAnalysis.coordinates.CRD`
===========================================================================

Read and write coordinates in CHARMM CARD coordinate format (suffix
"cord"). The CHARMM "extended format" is handled automatically.
"""
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import itertools
import warnings
from io import TextIOWrapper

import numpy as np
from MDAnalysis.coordinates import CRD
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib import util
from future.builtins import (
    open,
    range,
    str,
    super,
    zip,
)


class CORReader(CRD.CRDReader):
    """CRD reader that implements the standard and extended CRD coordinate formats

    .. versionchanged:: 0.11.0
       Now returns a ValueError instead of FormatError.
       Frames now 0-based instead of 1-based.
    """
    format = "COR"
    units = {"time": None, "length": "Angstrom"}


class CORWriter(CRD.CRDWriter):
    """COR writer that implements the CHARMM CRD EXT coordinate format.

    This class supercedes the original class provided by MDAnalysis by
    writing only the EXT format regardless of the number of atoms in the
    universe.

    Requires the following attributes:
    - resids
    - resnames
    - names
    - chainIDs
    - tempfactors

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """
    format = "COR"
    units = {"time": None, "length": "Angstrom"}

    fmt = dict(
        # crdtype = "extended"
        # fortran_format = "(2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)"
        ATOM_EXT=("{serial:10d}{totRes:10d}  {resname:<8.8s}  {name:<8.8s}"
                  "{pos[0]:20.10f}{pos[1]:20.10f}{pos[2]:20.10f}  "
                  "{chainID:<8.8s}  {resSeq:<8d}{tempfactor:20.10f}\n"),
        NUMATOMS_EXT="{0:10d} EXT\n",
        # crdtype = "standard"
        # fortran_format = "(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)"
        ATOM=("{serial:5d}{totRes:5d} {resname:<4.4s} {name:<4.4s}"
              "{pos[0]:10.5f}{pos[1]:10.5f}{pos[2]:10.5f} "
              "{chainID:<4.4s} {resSeq:<4d}{tempfactor:10.5f}\n"),
        TITLE="* FRAME {frame} FROM {where}",
        NUMATOMS="{0:5d}\n",
    )

    def __init__(self, filename, **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        """
        super().__init__(filename, **kwargs)
        self.filename = util.filename(filename, ext="cor")
        self.crd = None

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        write(selection,frame=FRAME)

        selection         MDAnalysis AtomGroup
        frame             optionally move to frame FRAME
        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single PDB (?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.positions  # can write from selection == Universe (Issue 49)

        n_atoms = len(atoms)
        # Detect which format string we"re using to output (EXT or not)
        # *len refers to how to truncate various things,
        # depending on output format!
        at_fmt = self.fmt["ATOM_EXT"]
        serial_len = 10
        resid_len = 8
        totres_len = 10

        # Check for attributes, use defaults for missing ones
        attrs = {}
        missing_topology = []
        for attr, default in (
            ("resnames", itertools.cycle(("UNK", ))),
                # Resids *must* be an array because we index it later
            ("resids", np.ones(n_atoms, dtype=np.int)),
            ("names", itertools.cycle(("X", ))),
            ("tempfactors", itertools.cycle((0.0, ))),
        ):
            try:
                attrs[attr] = getattr(atoms, attr)
            except (NoDataError, AttributeError):
                attrs[attr] = default
                missing_topology.append(attr)
        # ChainIDs - Try ChainIDs first, fall back to Segids
        try:
            attrs["chainIDs"] = atoms.segids
        except (NoDataError, AttributeError):
            # try looking for segids instead
            try:
                attrs["chainIDs"] = atoms.chainIDs
            except (NoDataError, AttributeError):
                attrs["chainIDs"] = itertools.cycle(("", ))
                missing_topology.append(attr)
        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "".format(miss=", ".join(missing_topology)))

        with open(self.filename, "wb") as crd:
            # Write Title
            crd.write(self.fmt["TITLE"].format(
                frame=frame, where=u.trajectory.filename).encode())
            crd.write("\n".encode())
            crd.write("*\n".encode())

            # Write NUMATOMS
            crd.write(self.fmt["NUMATOMS_EXT"].format(n_atoms).encode())

            # Write all atoms

            current_resid = 1
            resids = attrs["resids"]
            for i, pos, resname, name, chainID, resid, tempfactor in zip(
                    range(n_atoms), coor, attrs["resnames"], attrs["names"],
                    attrs["chainIDs"], attrs["resids"], attrs["tempfactors"]):
                if not i == 0 and resids[i] != resids[i - 1]:
                    current_resid += 1

                # Truncate numbers
                serial = int(str(i + 1)[-serial_len:])
                resid = int(str(resid)[-resid_len:])
                current_resid = int(str(current_resid)[-totres_len:])

                crd.write(
                    at_fmt.format(
                        serial=serial,
                        totRes=current_resid,
                        resname=resname,
                        name=name,
                        pos=pos,
                        chainID=chainID,
                        resSeq=resid,
                        tempfactor=tempfactor).encode())
