# 3D Periodic Halo Orbits
A computational framework for generating and analyzing periodic halo orbits in the CR3BP environment.

## Overview
This is my own independent research project exploring the generation of periodic halo orbit families around Lagrangian libration points. In its current state, Richardson's third-order analytical approximation for periodic halo orbits is used for initial conditions at two different z amplitudes which is utilized via a numerical continuation method to find the halo orbit family.

**Status**: This is currently a work in progress. As I learned more and more about Python, CR3BP, and the mathematics involved, expect significant changes.

## Current Capabilities
- Full Halo orbit family at $L_{1}$ including NRHOs and partial Halo orbit family @ $L_{2}$ excluding NRHOs
- Numerical continuation method
- Single shooting method

## Future Plans
- Restructure coding for better readability and maintainability
- Increase robustness of family generation via Pseudo Arc-length continuation
- Improve error handling by creating more user friendly error syntax
- Add visualization tools for orbit families
- Create examples and tutorials
- Write comprehensive mathematical documentation

## References
This work is been made possible by various people who spent their time working in this realm. Here are some papers that I've used so far and will be used to help advance this project:

- Howell, K.C., Pernicka, H.J. Numerical determination of Lissajous trajectories in the restricted three-body problem. Celestial Mechanics 41, 107–124 (1987). https://doi.org/10.1007/BF01238756
- Richardson, D.L. Analytic construction of periodic orbits about the collinear points. Celestial Mechanics 22, 241–253 (1980). https://doi.org/10.1007/BF01229511
- Breakwell, J.V., Brown, J.V. The ‘Halo’ family of 3-dimensional periodic orbits in the Earth-Moon restricted 3-body problem. Celestial Mechanics 20, 389–404 (1979). https://doi.org/10.1007/BF0123040
- Connor Howell, K. Three-dimensional, periodic, ‘halo’ orbits. Celestial Mechanics 32, 53–71 (1984). https://doi.org/10.1007/BF01358403
- Haapala, A. F. (n.d.). Trajectory design using periapse maps and invariant manifolds. Purdue e-Pubs. https://docs.lib.purdue.edu/dissertations/AAI1490654/
- W. Koon, M. Lo, J. Marsden, S. Ross, "Dynamical Systems, The Three-Body Problem, and Space Mission Design", 2006
- Howell, K.C., Breakwell, J.V. Almost rectilinear halo orbits. Celestial Mechanics 32, 29–52 (1984). https://doi.org/10.1007/BF01358402
