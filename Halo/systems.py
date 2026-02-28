"""
This file holds both gravitational mass parameters for different CR3BP systems 
as well as the dimensionalization values for distance, velocity and time.

Values were found from these sources:

Pitjeva, E. V. (2015). Determination of the Value of the Heliocentric Gravitational Constant (GM⊙) from Modern Observations of Planets and Spacecraft. Journal of Physical and Chemical Reference Data, 44(3). https://doi.org/10.1063/1.4921980
Planetary Satellite Physical Parameters. (n.d.). Ssd.jpl.nasa.gov. https://ssd.jpl.nasa.gov/sats/phys_par/
IAU Working Group Numerical Standards for Fundamental Astronomy. (n.d.). Iau-A3.Gitlab.io. https://iau-a3.gitlab.io/NSFA/NSFA_cbe.html#GME2009
Wang, S., Koon, Lo, M., Marsden, J., & Ross, S. (n.d.). Dynamical Systems, the Three-Body Problem and Space Mission Design. Retrieved February 24, 2026, from https://ross.aoe.vt.edu/books/Ross_3BodyProblem_Book_2022.pdf

Created on January 5th, 2026
Author: Dean Hall
"""

from __future__ import annotations
import re

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

@dataclass(frozen=True, slots=True)
class Dimensional:
    """Dimensionalization constants for a specific CR3BP primary-secondary pair."""
    name: str
    primary: str
    secondary: str

    gm1: float
    gm2: float

    L: float | None = None
    V: float | None = None
    T: float | None = None

    @property
    def mu(self) -> float:
        """CR3BP mass parameter mu = gm2 / (gm1 + gm2)."""
        return self.gm2 / (self.gm1 + self.gm2)

_SUN_JUPITER = Dimensional(
    name =  "Sun-Jupiter",
    primary = "Sun",
    secondary = 'Jupiter',
    gm1 = 1.32712440042e11,
    gm2 = 1.26686534e8,
    L = 7.784e8,
    V = 13.102,
    T = 3.733e8,
)

_SUN_EARTHMOON = Dimensional(
    name =  "Sun-EarthMoon",
    primary = "Sun",
    secondary = 'Earth + Moon',
    gm1 = 132712e6,
    gm2 = 3.986004418e5 + 4.902800118e3,
    L = 1.496e8,
    V = 29.784,
    T = 3.147e7,
)

_EARTH_MOON = Dimensional(
    name =  "Earth-Moon",
    primary = "Earth",
    secondary = 'Moon',
    gm1 = 3.986004418e5,
    gm2 = 4.902800118e3,
    L = 3.850e5,
    V = 1.025,
    T = 2.361e6,
)

_MARS_PHOBOS = Dimensional(
    name =  "Mars-Phobos",
    primary = "Mars",
    secondary = 'Phobos',
    gm1 = 4.282837362e4,
    gm2 = 0.7087e-3,
    L = 9.830e3,
    V = 2.144,
    T = 2.749e4,
)

_JUPITER_IO = Dimensional(
    name =  "Jupiter-Io",
    primary = "Jupiter",
    secondary = 'Io',
    gm1 = 1.26686534e8,
    gm2 = 5.95991547e3,
    L = 4.218e5,
    V = 17.390,
    T = 1.524e5,
)

_JUPITER_EUROPA = Dimensional(
    name =  "Jupiter-Europa",
    primary = "Jupiter",
    secondary = 'Europa',
    gm1 = 1.26686534e8,
    gm2 = 3.20271210e3,
    L = 6.711e5,
    V = 13.780,
    T = 3.060e5,
)

_JUPITER_GANYMEDE = Dimensional(
    name =  "Jupiter-Ganymede",
    primary = "Jupiter",
    secondary = 'Ganymede',
    gm1 = 1.26686534e8,
    gm2 = 9.88783275e3,
    L = 1.070e6,
    V = 10.909,
    T = 6.165e5,
)

_JUPITER_CALLISTO = Dimensional(
    name =  "Jupiter-Callisto",
    primary = "Jupiter",
    secondary = 'Callisto',
    gm1 = 1.26686534e8,
    gm2 = 7.17928340e3,
    L = 1.883e6,
    V = 8.226,
    T = 1.438e6,
)

_SATURN_MIMAS = Dimensional(
    name =  "Saturn-Mimas",
    primary = "Saturn",
    secondary = 'Mimas',
    gm1 = 3.7931187e7,
    gm2 = 2.50349e0,
    L = 1.856e5,
    V = 14.367,
    T = 8.117e4,
)

_SATURN_TITAN = Dimensional(
    name =  "Saturn-Titan",
    primary = "Saturn",
    secondary = 'Titan',
    gm1 = 3.7931187e7,
    gm2 = 8.97813710e3,
    L = 1.222e6,
    V = 5.588,
    T = 1.374e6,
)

_NEPTUNE_TRITON = Dimensional(
    name =  "Neptune-Triton",
    primary = "Neptune",
    secondary = 'Triton',
    gm1 = 6.836529e6,
    gm2 = 1.42849546e3,
    L = 3.548e5,
    V = 4.402,
    T = 5.064e5,
)

_PLUTO_CHARON = Dimensional(
    name =  "Pluto-Charon",
    primary = "Pluto",
    secondary = 'Charon',
    gm1 = 8.71e2,
    gm2 = 1.061e2,
    L = 1.941e4,
    V = 0.222,
    T = 5.503e5,
)

_SYSTEMS: Mapping[str, Dimensional] = MappingProxyType({
    "sun_jupiter": _SUN_JUPITER,
    "sun_earthmoon": _SUN_EARTHMOON,
    "earth_moon": _EARTH_MOON,
    "mars_phobos": _MARS_PHOBOS,
    "jupiter_io": _JUPITER_IO,
    "jupiter_europa": _JUPITER_EUROPA,
    "jupiter_ganymede": _JUPITER_GANYMEDE,
    "jupiter_callisto": _JUPITER_CALLISTO,
    "saturn_mimas": _SATURN_MIMAS,
    "saturn_titan":_SATURN_TITAN,
    "neptune_triton":_NEPTUNE_TRITON,
    "pluto_charon":_PLUTO_CHARON
})

def _canon(key: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")

def get_system_dim(key: str):
    k = _canon(key)
    try:
        return _SYSTEMS[k]
    except KeyError as e:
        options = ", ".join(sorted(_SYSTEMS.keys()))
        raise KeyError(f"Unknown CR3BP system '{key}'. Options: {options}") from e

# Example usage:
# sys = get_system_dim("Mars-Phobos")
# print(sys.mu, sys.L)
