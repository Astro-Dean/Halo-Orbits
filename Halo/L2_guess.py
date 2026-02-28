import numpy as np
from cr3bp_halo_family import CR3BP_Halo
from plot_orbits import plot_family, plot_manifold

# Guess from NASA JPL Three-Body Periodic Orbit Catalog
# https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
# ID 1524 on Date: February 24th, 2026 @ 9:53 pm EST
ic_guess = np.array([
    1.1808985497899205E+0,
    -2.5444988241150091E-26,
    1.0295054075242347E-4,
    3.3765359485568778E-15,
    -1.5585631393981156E-1,
    5.5263881873244218E-18
])

tf_guess = 3.4155308065628454E+0
system = 'earth_moon'

test_run = CR3BP_Halo(system=system, ic=ic_guess, tf=tf_guess)
family = test_run.npc()
family = test_run.palc(family, max_members=15)

plot_family(family=family, system=system, Lpt=2)
plot_manifold(system=system, family=family, member=17, Lpt=2)
