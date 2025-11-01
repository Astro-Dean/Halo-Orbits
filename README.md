# 3D Periodic Halo Orbits
This is my own independent research project that I do during my free time. 


This is currently a work in progress. As I learned more and more about Python, CR3BP, and the mathematics involved, expect significant changes. Most changes will come as mass updates.


The generation of halo orbits comes from the use of an analytical approximation via Richardson's third order solver which gives a good initial condition to develop halo orbits with a low $A_{z}$ amplitude. Using two consecutive orbits, both corrected via single shooting, we can create both the southern and northern halo orbit family about $L_{1}$ and $L_{2}$ Lagrangian Libration points using numerical continuation.

At this moment, the algorithm is capable of getting NRHOs from $L_{1}$ and close to NHROs from $L_{2}$. I will be adding psuedo arc-length continuation some time in the future to try to reach NRHOs for $L_{2}$
