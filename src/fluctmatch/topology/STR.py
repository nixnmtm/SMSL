# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.utils import (
    native_str,
)
from future.builtins import (
    dict,
    open,
    super,
)

from os import (path, environ)
import time

import numpy as np
import pandas as pd

from MDAnalysis.lib import util
from . import base


class STRWriter(base.TopologyWriterBase):
    format = "STREAM"
    units = dict(time=None, length="Angstrom")

    def __init__(self, filename, title=None, **kwargs):
        self.filename = util.filename(filename, ext="str")
        super().__init__(self.filename, **kwargs)

        self.title = ("* Created by fluctmatch on {}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())),
                      "* User: {}\n".format(environ["USER"])) if title is None else title
        self.fmt = "IC EDIT\nDIST %4s %10d %4s %4s %10d %4s%5.1f\nEND\n"

    def write(self, universe):
        # type: (object) -> object
        if not hasattr(universe.atoms, "bonds"):
            print("No bonds were found.")
            return

        # Create the table
        a1, a2 = universe.atoms.bonds.atom1, universe.atoms.bonds.atom2
        data = [a1.segids, a1.resids, a1.names, a2.segids, a2.resids, a2.names, np.zeros(a1.n_atoms, dtype=np.float)]
        data = np.concatenate([pd.DataFrame(_) for _ in data], axis=1)

        # Write the data to the file.
        with open(self.filename, "wb") as strfile:
            for _ in self.title:
                print(_.encode(), file=strfile)
            np.savetxt(strfile, data, fmt=native_str(self.fmt))
