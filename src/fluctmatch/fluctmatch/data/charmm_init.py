# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# fluctmatch --- https://github.com/tclick/python-fluctmatch
# Copyright (c) 2015-2017 The pySCA Development Team and contributors
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
import numpy as np
import pandas as pd


init = ("""
    * Calculate the initial equilibrium distance and force constant for an ENM
    * structure. The original CHARMM code was written by Prof. Jhih-Wei Chu.
    * 
    
    {dimension}

    set version {version}
    bomlev -5 ! This is for CHARMM 39

    ! Additional information
    set fileu   10

    ! Open CHARMM topology and parameter file
    read rtf card name "{topology_file}"

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

    ! Determine "experimental" normal mode analysis information and output to a
    ! file using the "extended" format, which will be used for all file output.
    ioformat extended
    ic fill
    
    ! Load the trajectories
    open read unit @fileu file name {traj_file}
    
    ! Gather information from the first trajectory assuming that all 
    ! trajectories are similar.
    ! Calculate the beginning and final times for a trajectory sequence.
    traj query unit @fileu
    
    ! Calculate internal coordinate average movement
    ic dyna aver first @fileu nunit 1 skip ?SKIP
    write ic card resid name {init_avg_ic}
    * Internal coordinate averages
    
    ! Calculate fluctuations in internal coordinate movement
    ic dyna fluc first @fileu nunit 1 skip ?SKIP
    write ic card resid name {init_fluct_ic}
    * Internal coordinate fluctuations
    
    close unit @fileu
    
    stop
""")
