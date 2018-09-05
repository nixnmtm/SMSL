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

thermodynamics = ("""
    * Calculate thermodynamic qunatities for an ENM structure.
    *

    {dimension}

    set version {version}
    bomlev -5 ! This is for CHARMM 39

    ! Additional data
    set ndcd    1
    set temp    {temperature}
    set fileu   10

    ! Open CHARMM topology and parameter file
    read rtf  card name "{topology_file}"
    read para card {flex} name "{fixed_prm}"

    ! Open PSF and coordinate files
    if @version .ge. 39 then
        read psf card name "{xplor_psf_file}"
    else
        read psf card name "{psf_file}"
    endif
    read coor card name "{crd_file}"
    coor copy comp

    stream "{stream_file}"

    skip all excl bond
    update inbfrq 0

    ! Load the trajectory
    open read  unit @fileu file name {trajectory}

    ! Gather information from the first trajectory assuming that all trajectories
    ! are similar.
    traj query unit @fileu

    calc nmod = 3*?NATOM
    vibran nmodes @nmod
        coor dyna sele all end nopr first @fileu nunit @ndcd begin ?START -
            nskip ?SKIP orient sele all end
        quasi first @fileu nunit @ndcd nskip ?SKIP begin ?START sele all end -
            temp @temp thermo resi
    end
    calc ts = ?stot * @temp
    close unit @fileu

    stop
    """)
