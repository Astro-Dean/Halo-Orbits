# 3D Periodic Halo Orbits
This is my own independent research project that I do during my free time. 

This is currently a work in progress. As I learned more and more about Python, CR3BP, and the mathematics involved, expect significant changes. Most changes will come as mass updates.


The generation of halo orbits comes from the use of an analytical approximation via Richardson's third order solver which gives a good initial condition to develop halo orbits with a low $A_{z}$ amplitude. Using two consecutive orbits, both corrected via single shooting, we can create both the southern and northern halo orbit family about $L_{1}$ and $L_{2}$ Lagrangian Libration points using numerical continuation.

At this moment, the algorithm is capable of getting NRHOs from $L_{1}$ and close to NHROs from $L_{2}$. I will be adding psuedo arc-length continuation some time in the future to try to reach NRHOs for $L_{2}$

This work is been made possible by various people who spent thier time working in this realm. Here are some papers that I've used so far and will be used to help advance this project:

- Howell, K.C., Pernicka, H.J. Numerical determination of Lissajous trajectories in the restricted three-body problem. Celestial Mechanics 41, 107–124 (1987). https://doi.org/10.1007/BF01238756
- Richardson, D.L. Analytic construction of periodic orbits about the collinear points. Celestial Mechanics 22, 241–253 (1980). https://doi.org/10.1007/BF01229511
- Breakwell, J.V., Brown, J.V. The ‘Halo’ family of 3-dimensional periodic orbits in the Earth-Moon restricted 3-body problem. Celestial Mechanics 20, 389–404 (1979). https://doi.org/10.1007/BF0123040
- Connor Howell, K. Three-dimensional, periodic, ‘halo’ orbits. Celestial Mechanics 32, 53–71 (1984). https://doi.org/10.1007/BF01358403
- Haapala, A. F. (n.d.). Trajectory design using periapse maps and invariant manifolds. Purdue e-Pubs. https://docs.lib.purdue.edu/dissertations/AAI1490654/
- W. Koon, M. Lo, J. Marsden, S. Ross, "Dynamical Systems, The Three-Body Problem, and Space Mission Design", 2006
- Howell, K.C., Breakwell, J.V. Almost rectilinear halo orbits. Celestial Mechanics 32, 29–52 (1984). https://doi.org/10.1007/BF01358402
