import numpy as np
from cr3bp_halo_family import CR3BP_Halo
from plot_orbits import plot_family, plot_manifold

# Guess from NASA JPL Three-Body Periodic Orbit Catalog
# https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
# ID 819 on Date: February 28th, 2026 @ 5:39 pm EST
ic_guess = np.array([
    9.7748175708167062E-1,
    -1.2844955661062785E-21,
    1.5560134012085195E-3,
    2.3268901724697762E-16,
    1.7973750978536162E-2,
    -4.2927254818154733E-17
])

tf_guess = 3.0378505750205589E+0
system = 'jupiter_europa'

test_run = CR3BP_Halo(system=system, ic=ic_guess, tf=tf_guess, tol=1e-12)
family = test_run.npc(max_members=20, step=1e-4)
family = test_run.palc(family, max_members=95, step=1e-4)

plot_family(family=family, system=system, Lpt=1)
plot_manifold(system=system, family=family, member=35, prop_time_scale=1)
