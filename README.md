# 3D Periodic Halo Orbits
The generation of halo orbits comes from the use of an analytical approximation via Richardson's third order solver which gives a good initial condition to develop halo orbits with a low $A_{z}$ amplitude. Using two consecutive orbits, both corrected via single shooting, we can create both the southern and northern halo orbit family about $L_{1}$ and $L_{2}$ Lagrangian Libration points using numerical continuation.

This is currently a work in progress.

At this moment, the algorithm is capable of getting NHROs from $L_{1}$ and close to NHROs from $L_{2}$. I will be adding psuedo arc-length continuation some time in the future to try to reach NHROs for $L_{2}$
