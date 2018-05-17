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

split_inp = ("""
    * Create a subtrajectory from a larger CHARMM trajectory.
    * This is for <= c35.
    *

    set version {version}
    if @version .ge. 36 then
        dimension chsize 500000 maxres 3000000
    endif
    bomlev -5 ! This is for CHARMM 39

    set toppar {toppar}
    set begin {start}
    set end {stop}

    if @version .lt. 36 then
        read rtf  card name @toppar/top_all27_prot_na.rtf
        read para card {flex} name @toppar/par_all27_prot_na.prm

        read rtf  card name @toppar/top_all27_lipid.rtf append
        read para card {flex} name @toppar/par_all27_lipid.prm append

        read rtf  card name @toppar/top_all35_sugar.rtf append
        read para card {flex} name @toppar/par_all35_sugar.prm append

        stream @toppar/stream/toppar_water_ions.str
    else
        read rtf  card name @toppar/top_all36_prot.rtf
        read para card {flex} name @toppar/par_all36_prot.prm

        read rtf  card name @toppar/top_all36_na.rtf append
        read para card {flex} name @toppar/par_all36_na.prm append

        read rtf  card name @toppar/top_all36_carb.rtf append
        read para card {flex} name @toppar/par_all36_carb.prm append

        read rtf  card name @toppar/top_all36_lipid.rtf append
        read para card {flex} name @toppar/par_all36_lipid.prm append

        read rtf  card name @toppar/top_all36_cgenff.rtf append
        read para card {flex} name @toppar/par_all36_cgenff.prm append

        stream @toppar/toppar_water_ions.str
    endif

    read psf card name {psf}

    set iu    20
    set ou    30

    ! Load the trajectory
    open read  unit @iu file name {trajectory}
    open write unit @ou file name {outfile}

    ! Gather information from the first trajectory assuming that all trajectories
    ! are similar.
    traj query unit @iu

    if @start .gt. 1 then
        calc nother @start - 1
    else
        set nother @start
    endif
    traj first @iu skip ?SKIP begin @start stop @stop -
        iwrite @ou noth @nother

    close @iu
    close @ou
    stop
    """)
