import numpy as np
from cr3bp_halo_family import CR3BP_Halo
from plot_orbits import plot_family, plot_manifold

# Guess from NASA JPL Three-Body Periodic Orbit Catalog
# https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
# ID 1147 on Date: February 24th, 2026 @ 9:41 pm EST
ic_guess = np.array([
    8.2339081983651485E-1,
    -1.9017764504099543E-28,
    9.8941366235910004E-4,
    -2.3545391932685812E-15,
    1.2634272983881797E-1,
    2.2367029429442455E-16
])

tf_guess = 2.7430007981241529E+0
system = 'earth_moon'

test_run = CR3BP_Halo(system=system, ic=ic_guess, tf=tf_guess)
family = test_run.npc()
family = test_run.palc(family, max_members=100)

plot_family(family=family, system=system, Lpt=1)
plot_manifold(system=system, family=family, member=27, prop_time_scale=1.5)