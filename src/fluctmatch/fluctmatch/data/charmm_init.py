init = ("""
    * Calculate the initial equilibrium distance and force constant for an ENM
    * structure. The original CHARMM code was written by Prof. Jhih-Wei Chu.
    *

    {dimension}

    set version {version}
    bomlev -5

    set fileu   10

    read rtf card name "{topology_file}"

    read psf card name "{xplor_psf_file}"
    read coor card name "{crd_file}"
    coor copy comp

    stream "{stream_file}"

    skip all excl bond
    update inbfrq 0

    ioformat extended
    ic fill

    open read unform unit @fileu file name "{traj_file}"

    traj query unit @fileu

    ic dyna aver first @fileu nunit 1 skip ?SKIP
    write ic card resid name "{init_avg_ic}"
    * Internal coordinate averages

    ic dyna fluc first @fileu nunit 1 skip ?SKIP
    write ic card resid name "{init_fluct_ic}"
    * Internal coordinate fluctuations

    close unit @fileu

    stop
""")