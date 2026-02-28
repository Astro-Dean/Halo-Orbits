import numpy as np
from cr3bp_halo_family import CR3BP_Halo
from plot_orbits import plot_family, plot_manifold

# Guess from NASA JPL Three-Body Periodic Orbit Catalog
# https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
# ID 1030 on Date: February 24th, 2026 @ 10:25 pm EST
ic_guess = np.array([
    -1.6967233092486720E+0,
    3.1756007279014768E-23,
    1.0136396075463166E-2,
    -5.8932494929994690E-13,
    1.2795419587406101E+0,
    -1.3980383346242373E-14
])

tf_guess = 6.2391471092964244E+0
system = 'earth_moon'

test_run = CR3BP_Halo(system=system, ic=ic_guess, tf=tf_guess)
family = test_run.npc()
family = test_run.palc(family, max_members=300)

plot_family(family=family, system=system, Lpt=3)
plot_manifold(system=system, family=family, member=200, Lpt=3, prop_time_scale=2, show_earth=True)
