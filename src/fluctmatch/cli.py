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
"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mfluctmatch` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``fluctmatch.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``fluctmatch.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import os
import sys
from os import path
import click

CONTEXT_SETTINGS = dict(
    auto_envvar_prefix='COMPLEX',
    help_option_names=['-h', '--help'],
)


class Context(object):
    """Context manager for click command-line interface."""
    def __init__(self):
        self.verbose = False
        self.home = os.getcwd()

    def log(self, msg, *args):
        """Logs a message to stderr."""
        if args:
            msg %= args
        click.echo(msg, file=sys.stderr)

    def vlog(self, msg, *args):
        """Logs a message to stderr only if verbose is enabled."""
        if self.verbose:
            self.log(msg, *args)


pass_context = click.make_pass_decorator(Context, ensure=True)
cmd_folder = path.abspath(
    path.join(path.dirname(__file__), 'commands')
)


class ComplexCLI(click.MultiCommand):
    """Complex command-line options with subcommands for fluctmatch.
    """
    def list_commands(self, ctx):
        """List available commands.

        Parameters
        ----------
        ctx : :object:`Context`
            click context

        Returns
        -------
            List of available commands
        """
        rv = []
        for filename in os.listdir(cmd_folder):
            if (filename.endswith(".py") and filename.startswith("cmd_")):
                rv.append(filename[4:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        """Run the selected command

        Parameters
        ----------
        ctx : :class:`Context`
            click context
        name : str
            command name

        Returns
        -------
            The chosen command if present
        """
        try:
            if sys.version_info[0] == 2:
                name = name.encode("ascii", "replace")
            mod = __import__(
                "fluctmatch.commands.cmd_" + name,
                None,
                None,
                ["cli"]
            )
        except ImportError:
            return
        return mod.cli


@click.command(cls=ComplexCLI, context_settings=CONTEXT_SETTINGS)
@click.version_option()
def main():
    pass
