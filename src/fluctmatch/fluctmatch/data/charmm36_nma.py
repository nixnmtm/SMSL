# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

nma = (
"""* Normal mode analysis of structure for parameter fitting. The original CHARMM
* script was written by Prof. Jhih-Wei Chu
*

bomlev -5 ! This is for CHARMM 39

! Additional information
set temp    {temperature}
set fileu   10
set fluctu  20
set vibu    30

! Open CHARMM topology and parameter file
read rtf  card name "{topology_file}"
read para card name "{fixed_prm}" flex

! Open PSF and coordinate files
read psf  card name "{xplor_psf_file}" xplor
read coor card name "{crd_file}"
coor copy comp

skip all excl bond
update inbfrq 0

ener

! Minimize structure using steepest descent and ABNR
mini   sd nstep 100
mini abnr nstep 2000

coor orie rms mass
scalar wmain copy mass

ioformat extended
write coor card name "{nma_crd}"

stream "{stream_file}"

ic fill
write ic unit @fileu card resid name "{avg_ic}"

calc nmode   ?natom * 3

set nmodes   @nmode
set type     temp

open write unit @fluctu card name "{fluct_ic}"
open write unit @vibu   card name "{nma_vib}"

! Perform normal mode analysis at desired temperature for vibrational normal
! modes
vibran nmode @nmodes
    diag fini
    fluc ic @type @temp tfre 0.0 mode 7 thru @nmodes
    ic save
    ic write unit @fluctu resid
    write normal card mode 1 thru @nmodes unit @vibu
end

stop
"""
)


