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

from . import base


class STRWriter(base.TopologyWriterBase):
    """Write a stream file to define internal coordinates within CHARMM.

    Parameters
    ----------
    filename : str
        The filename to which the internal coordinate definitions are written.
    n_atoms : int, optional
        The number of atoms in the output trajectory.
    title : str or list of str, optional
        A header section written at the beginning of the stream file. If no title
        is given, a default title will be written.
    """
    format = "STREAM"
    units = dict(time=None, length="Angstrom")

    def __init__(self, filename, n_atoms=None, title=None):
        self.filename = util.filename(filename, ext=native_str("stream"))
        self.n_atoms = n_atoms
        if title is None:
            self.title = [
                "* Created by fluctmatch on {}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())),
                "* User: {}".format(environ["USER"]),
            ]
        elif issubclass(type(title), str) or issubclass(type(title), unicode):
            self.title = (title, )
        else:
            self.title = title
        self.fmt = "IC EDIT\nDIST %4s %10d %4s %4s %10d %4s%5.1f\nEND\n"

    def write(self, universe):
        """Write the bond information to a CHARMM-formatted stream file.

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or atomgroup with bond definitions.
        """
        if self.n_atoms is not None and self.n_atoms != universe.atoms.n_atoms:
            raise IOError(
                "The number of atoms from object initialization do not match the number of atoms "
                "from the universe [{:d} != {:d}]".format(self.n_atoms, universe.atoms.n_atoms)
            )

        # Create the table
        try:
            a1, a2 = universe.atoms.bonds.atom1, universe.atoms.bonds.atom2
            data = [
                a1.segids, a1.resids, a1.names,
                a2.segids, a2.resids, a2.names,
                np.zeros(a1.n_atoms, dtype=np.float)
            ]
            data = np.concatenate([pd.DataFrame(_) for _ in data], axis=1)

        except AttributeError:
            raise_with_traceback(AttributeError("No bonds were found."))

        # Write the data to the file.
        with open(self.filename, "wb") as stream_file:
            for _ in self.title:
                stream_file.write((_ + "\n").encode())
            np.savetxt(stream_file, data, fmt=native_str(self.fmt))
