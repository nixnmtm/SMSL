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

afqha = ("""
    * Compute QHA and get atomic fluctuation from the eigenmodes
    * Author: Nixon Raj N
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
    
    ! Read trajectory file
    open read unit 51 file name "{traj_file}"
    
    ! open write file for rms fit 
    open unit 61 write unform name "{fit_traj_file}"
    
    
    !define core sele segid S1A .and. ( resid 1 : 2 .or. resid 15 : 18 .or. resid 25 : 31 .or. resid 33 : 40 -
    !                                  .or. resid 48 : 50 .or. resid 52 .or. resid 70 : 73 .or. resid 83 : 90 .or. resid 102 : 105 -
    !                                  .or. resid 109 .or. resid 116 : 124 .or. resid 135 : 139 .or. resid 142 .or. resid 160 : 162 .or. -
    !                                  resid 171 : 183 .or. resid 186 : 195 .or. resid 204 : 210 .or. resid 212 : 216 ) end
    
    ! Calculate the average structure of trajectory and save to comparsion set
    coor dyna comp firstu 51 nunit 1 sele all end orient mass sele all end
    
    ! Write the average structure
    open write unit 75 card name "{qha_avg_cor}"
    write coor card unit 75
    *
    close unit 75
    
    
    ! Mass-weighted fit trajectory to average structure
    
    merge coor firstu 51 nunit 1 output 61 begin 1 skip 1 nfile 10000 sele all end orie mass
    
    close unit 51
    
    
    read coor card name "{qha_avg_cor}"
    coor copy comp
    
    open unit 51 read unform name "{fit_traj_file}"
    
    set nmod ?natom
    calc nmod = @nmod * 3
    
    set nmodes @nmod
    
    nbonds inbfrq 0
    
    ! Atomic fluctuations of first low frequency mode calculated
    vibran nmod @nmod
        quasi temp @temp tfre 0.0 firstu 51 nunit 1 begin 1 skip 1
        open write unit 71 card name "{quasi_nmv_file}"
        write norm card mode 1 thru @nmodes unit 71
        fluc atom mode @nmodes thru @nmodes temp 300 tfre 0.0 
        print norm temp @temp tfre 0.0
    end
    
    print coor comp
    
    write coor comp card name "{qhacor_file}"
    
    close unit 51
        
    stop

    """)
