#
# Released under the New BSD license.
#
# Please cite your use of fluctmatch in published work:
#
# Timothy H. Click, Nixon Raj, and Jhih-Wei Chu.
# Calculation of Enzyme Fluctuograms from All-Atom Molecular Dynamics
# Simulation. Meth Enzymology. 578 (2016), 327-342,
# doi:10.1016/bs.mie.2016.05.024.
# This code authored by: Nixon
#
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

afnma = ("""
    * Author: Nixon
    * Get the magnitude of atomic fluctuations from all normal modes
    *

    {dimension}

    set version {version}

    bomlev -2


    ! Additional information
    set temp    {temperature}

    ! Read topology file and parameter file
    read rtf  card name "{topology_file}"
    read para card {flex} name "{fixed_prm}"

    ! Read psf and coordinate file
    if @version .ge. 39 then
        read psf card name "{xplor_psf_file}"
    else
        read psf card name "{psf_file}"
    endif
    read coor card name "{crd_file}"

    skip all excl bond
    update inbfrq 0

    ener

    ! Minimize structure using steepest descent and ABNR
    mini   sd nstep 100
    mini abnr nstep 2000

    coor copy comp

    calc nmode   ?natom * 3

    set nmodes   @nmode
    set type     temp

    ! Perform normal mode analysis at desired temperature for vibrational normal
    ! mode and get the atomic fluctuations of first low frequency mode

    vibran nmode @nmodes
        diag fini
        print norm @type @temp tfre 0.0
        !open write unit 61 card name "{nma_nmv_file}"
        !write norm card mode 1 thru @nmodes unit 61
        fluc atom mode 7 thru 7 @type @temp tfre 0.0

    end

    print coor comp
    write coor comp card name "{nmacor_file}"

    stop

""")
