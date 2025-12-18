import numpy as np
from halo_orbit_tools import HaloOrbitFamily

mu = 1.215053e-2

halo = HaloOrbitFamily(mu=mu, Lpt=1, branch=1)
halo.plot_halo_family(mu=mu, Lpt=1, branch=1)
