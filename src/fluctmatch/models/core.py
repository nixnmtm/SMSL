# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals
)

import warnings

from future.utils import (
    viewkeys,
    raise_with_traceback,
)

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
    models = kwargs.pop("model", ["ncsc",])
    models = [_.upper() for _ in models]
    if "ENM" in models:
        warnings.warn("ENM model detected. All other models are being ignored.")
        universe = _MODELS["ENM"](*args, **kwargs)
        return universe

    universe = []
    try:
        for model in models:
            universe.append(_MODELS[model](*args, **kwargs))
    except KeyError:
        msg = (
            "{0} is not an available model. "
            "Please try {1}".format(model, viewkeys(_MODELS))
        )
        raise_with_traceback(KeyError(msg))

    universe = Merge(*universe)
    return universe
