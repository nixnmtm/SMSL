# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# SMSL - https://github.com/nixnmtm/SMSL
#


import os
import sys

import click

try:
    import parmed as pmd
except ImportError:
    pmd = None


@click.command(
    "gro2psf",
    short_help="Convert a GROMACS GRO/TOP system into a PSF using ParmEd.",
)
@click.option(
    "-g",
    "--gro",
    "gro_file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Input full-system GRO file.",
)
@click.option(
    "-p",
    "--top",
    "top_file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Input master topology TOP file that includes all split ITPs.",
)
@click.option(
    "-o",
    "--out",
    "out_file",
    type=click.Path(dir_okay=False, resolve_path=True),
    default="system.psf",
    show_default=True,
    help="Output PSF filename.",
)
def cli(gro_file, top_file, out_file):
    """Convert GRO + TOP/ITP topology into a PSF."""
    if pmd is None:
        raise click.ClickException(
            "ParmEd is not installed. Install with: conda install -c conda-forge parmed"
        )

    gro_file = os.path.abspath(gro_file)
    top_file = os.path.abspath(top_file)
    out_file = os.path.abspath(out_file)

    if not os.path.isfile(gro_file):
        raise click.ClickException("GRO file not found: {}".format(gro_file))

    if not os.path.isfile(top_file):
        raise click.ClickException("TOP file not found: {}".format(top_file))

    top_dir = os.path.dirname(top_file)
    top_basename = os.path.basename(top_file)

    old_cwd = os.getcwd()
    try:
        # ParmEd resolves relative includes from cwd
        os.chdir(top_dir)

        click.echo("Loading topology: {}".format(top_file))
        click.echo("Loading coordinates: {}".format(gro_file))

        struct = pmd.load_file(top_basename, xyz=gro_file)

        click.echo("Loaded {} atoms".format(len(struct.atoms)))
        click.echo("Loaded {} bonds".format(len(struct.bonds)))

        click.echo("Writing PSF: {}".format(out_file))
        struct.save(out_file, overwrite=True)

        click.echo("Done.")
    except Exception as exc:
        raise click.ClickException(
            "ERROR while converting GRO/TOP to PSF:\n{}".format(exc)
        )
    finally:
        os.chdir(old_cwd)