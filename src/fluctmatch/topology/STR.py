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

import textwrap
import time
from os import (environ)

import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from future.builtins import (
    dict,
)
from future.utils import (
    native_str,
)
from future.utils import (
    raise_with_traceback,
)

from fluctmatch.topology import base as topbase


class STRWriter(topbase.TopologyWriterBase):
    """Write a stream file to define internal coordinates within CHARMM.

    Parameters
    ----------
    filename : str
        The filename to which the internal coordinate definitions are written.
    n_atoms : int, optional
        The number of atoms in the output trajectory.
    title
        A header section written at the beginning of the stream file.
        If no title is given, a default title will be written.
    """
    format = "STREAM"
    units = dict(time=None, length="Angstrom")

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext=native_str("stream"))
        self._version = kwargs.get("charmm_version", 41)

        width = 4 if self._version < 36 else 8
        if self._version >= 36:
            self.fmt = (
                """
                IC EDIT
                DIST %-{width}s %{width}d %-{width}s %-{width}s %{width}d %-{width}s%{width}.1f
                END
                """.format(width=width)
            )
        else:
            self.fmt = (
                """
                IC EDIT
                DIST BYNUM %{width}d BYNUM %{width}d %{width}.1f
                END
                """.format(width=width)
            )

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

    def write(self, universe):
        """Write the bond information to a CHARMM-formatted stream file.

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or atomgroup with bond
            definitions.
        """
        # Create the table
        try:
            dist = np.zeros_like(universe.atoms.bonds.bonds(), dtype=np.float)
            if self._version >= 36:
                a1, a2 = universe.atoms.bonds.atom1, universe.atoms.bonds.atom2
                data = (
                    a1.segids, a1.resids, a1.names,
                    a2.segids, a2.resids, a2.names,
                    dist,
                )
                data = pd.concat([pd.DataFrame(_) for _ in data], axis=1)
            else:
                data = universe._topology.bonds.values
                data = np.concatenate((data, dist[:, np.newaxis]), axis=1)
        except AttributeError:
            raise_with_traceback(AttributeError("No bonds were found."))

        # Write the data to the file.
        with util.openany(self.filename, "w") as stream_file:
            for _ in self._title:
                print(_, file=stream_file)
            np.savetxt(
                stream_file,
                data,
                fmt=native_str(textwrap.dedent(self.fmt[1:]))
            )
            print("\nRETURN", file=stream_file)
