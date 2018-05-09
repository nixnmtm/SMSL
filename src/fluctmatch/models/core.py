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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.utils import (
    viewkeys,
    raise_with_traceback,
)

import warnings

from fluctmatch import _MODELS
from fluctmatch.models import *
from fluctmatch.models.base import Merge


def modeller(*args, **kwargs):
    """Create coarse-grain model from universe selection.

    Parameters
    ----------
    args
    kwargs

    Returns
    -------
    A coarse-grain model
    """
    models = kwargs.pop("model", [
        "ncsc",
    ])
    models = [_.upper() for _ in models]
    if "ENM" in models:
        warnings.warn(
            "ENM model detected. All other models are being ignored.")
        universe = _MODELS["ENM"](*args, **kwargs)
        return universe

    try:
        universe = [_MODELS[model](*args, **kwargs) for model in models]
    except KeyError:
        msg = ("{0} is not an available model. "
               "Please try {1}".format(model, viewkeys(_MODELS)))
        raise_with_traceback(KeyError(msg))

    universe = Merge(*universe) if len(universe) > 1 else universe[0]
    return universe
