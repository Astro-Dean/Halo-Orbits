# 3D Periodic Halo Orbits
A computational framework for generating and analyzing periodic halo orbits in the CR3BP environment.

## Overview
This is my own independent research project exploring the generation of periodic halo orbit families around Lagrangian libration points. In its current state, initial guesses fron NASA JPL's Three-Body Periodic Orbit Catalog are used with single shooting as well as Natural Parameter and Pseudo Arc-Length Continuation to generate full Halo Orbits in different systems. Four tests files exist shows generate of halo orbits at $L_1$, $L_2$, and $L_3$ in the Earth-Moon system as well as a test at $L_1$ in the Jupiter-Europa system.

**Status**: This is currently a work in progress. As I learned more and more about Python, CR3BP, and the mathematics involved, expect significant changes. 

## Current Capabilities
- Full Halo orbit family at Collinear Lagrangian Points including NRHOs
- Natural Parameter and Pseudo Arc-length Continuation
- Single shooting method
- Orbit Family and Invariant manifold plotting with options
- Restructured coding for better readability and maintainability
- Create database of different CR3BP systems

## Future Plans
- Create examples and tutorials
- Write comprehensive mathematical documentation
- Increase robustness of shooting method (most likely with multiple shooting)

## References
This work is been made possible by various people who spent their time working in this realm. Here are some papers that I've used so far and will be used to help advance this project:

- Howell, K.C., Pernicka, H.J. Numerical determination of Lissajous trajectories in the restricted three-body problem. Celestial Mechanics 41, 107–124 (1987). https://doi.org/10.1007/BF01238756
- Richardson, D.L. Analytic construction of periodic orbits about the collinear points. Celestial Mechanics 22, 241–253 (1980). https://doi.org/10.1007/BF01229511
- Breakwell, J.V., Brown, J.V. The ‘Halo’ family of 3-dimensional periodic orbits in the Earth-Moon restricted 3-body problem. Celestial Mechanics 20, 389–404 (1979). https://doi.org/10.1007/BF0123040
- Connor Howell, K. Three-dimensional, periodic, ‘halo’ orbits. Celestial Mechanics 32, 53–71 (1984). https://doi.org/10.1007/BF01358403
- Haapala, A. F. (n.d.). Trajectory design using periapse maps and invariant manifolds. Purdue e-Pubs. https://docs.lib.purdue.edu/dissertations/AAI1490654/
- Howell, K.C., Breakwell, J.V. Almost rectilinear halo orbits. Celestial Mechanics 32, 29–52 (1984). https://doi.org/10.1007/BF01358402
- Parker, J., & Anderson, R. (2013). LOW-ENERGY LUNAR TRAJECTORY DESIGN. https://descanso.jpl.nasa.gov/monograph/series12/LunarTraj--Overall.pdf
- Wang, S., Koon, Lo, M., Marsden, J., & Ross, S. (n.d.). Dynamical Systems, the Three-Body Problem and Space Mission Design. Retrieved February 24, 2026, from https://ross.aoe.vt.edu/books/Ross_3BodyProblem_Book_2022.pdf
