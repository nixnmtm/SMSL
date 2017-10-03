# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
    open,
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
        A header section written at the beginning of the stream file. If no title
        is given, a default title will be written.
    """
    format = "STREAM"
    units = dict(time=None, length="Angstrom")

    def __init__(self, filename, **kwargs):
        width = 4 if kwargs.get("charmm_version", 41) < 36 else 8
        self.fmt = (
            """
            IC EDIT
            DIST %-{width}s %{width}d %-{width}s %-{width}s %{width}d %-{width}s%{width}.1f
            END
            """.format(width=width)
        )
        self.filename = util.filename(filename, ext=native_str("stream"))

        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        user = environ["USER"]
        self._title = kwargs.get(
            "title",
            (
                "* Created by fluctmatch on {date}".format(date=date),
                "* User: {user}".format(user=user),
            )
        )
        if issubclass(type(self._title), str) or issubclass(type(self._title), np.unicode):
            self._title = (self._title,)

    def write(self, universe):
        """Write the bond information to a CHARMM-formatted stream file.

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or atomgroup with bond definitions.
        """
        # Create the table
        try:
            a1, a2 = universe.atoms.bonds.atom1, universe.atoms.bonds.atom2
            data = (
                a1.segids, a1.resids, a1.names,
                a2.segids, a2.resids, a2.names,
                np.zeros(a1.n_atoms, dtype=np.float)
            )
            data = pd.concat([pd.DataFrame(_) for _ in data], axis=1)
        except AttributeError:
            raise_with_traceback(AttributeError("No bonds were found."))

        # Write the data to the file.
        with open(self.filename, "wb") as stream_file:
            for _ in self._title:
                stream_file.write((_ + "\n").encode())
            np.savetxt(stream_file, data, fmt=native_str(textwrap.dedent(self.fmt)), delimiter=native_str(""))
