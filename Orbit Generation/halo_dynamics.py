import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class PotentialDerivatives:
    # First-order PDEs
    Ux: float
    Uy: float
    Uz: float
    # Second-order PDEs
    Uxx: float
    Uyy: float
    Uzz: float
    Uxy: float
    Uxz: float
    Uyz: float

class Halo_Orbit_Dynamics:
    def __init__(self, mu):
        self.mu = mu
    
    def distances(self, state):
        x, y, z = state[:3]
        mu = self.mu
        mu2 = 1 - mu
        r13 = np.sqrt((x+mu)**2  + y**2 + z**2)
        r23 = np.sqrt((x-mu2)**2 + y**2 + z**2)
        return r13, r23
    
    def potentials(self, state):
        
        x, y, z = state[:3]
        mu = self.mu
        mu2 = 1 - mu

        # Distances from primaries
        r13, r23 = self.distances(state)

        # First-Order PDE of Psuedo-potential function U(x,y,z)
        Ux = x - (mu2)*(x+mu)/r13**3 - mu*(x-mu2)/r23**3
        Uy = y - mu2*y/r13**3        - mu*y/r23**3
        Uz =   - mu2*z/r13**3        - mu*z/r23**3

        # Second-Order PDE of psuedo-potential function U(x,y,z)
        Uxx = 1 - mu2/r13**3 + 3*mu2*(x+mu)**2/r13**5 - mu/r23**3 + 3*mu*(x-mu2)/r23**5
        Uyy = 1 - mu2/r13**3 + 3*mu2*y**2/r13**5   - mu/r23**3 + 3*mu*y**2/r23**5
        Uzz =   - mu2/r13**3 + 3*mu2*z**2/r13**5   - mu/r23**3 + 3*mu*z**2/r23**5
        
        # Mixed Second-Order of psuedo-potential function U(x,y,z)
        Uxy = 3*(mu2*(mu+x)*y/r13**5 + mu*(x-mu2)*y/r23**5)
        Uxz = 3*(mu2*(x+mu)*z/r13**5 + mu*(x-mu2)*z/r23**5)
        Uyz = 3*y*z*(mu2/r13**5 + mu/r23**5)

        return PotentialDerivatives(Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz)

    def A_matrix(self, state):
        derivs = self.potentials(state)
        Atl = np.zeros((3,3))
        Atr = np.eye(3)
        Abr = np.array([
            [ 0, 1, 0],
            [-1, 0, 0],
            [ 0, 0, 0]
        ])
        Abl = np.array([
            [derivs.Uxx, derivs.Uxy, derivs.Uxz],
            [derivs.Uxy, derivs.Uyy, derivs.Uyz],
            [derivs.Uxz, derivs.Uyz, derivs.Uzz]
        ])
        
        A = np.block([
            [Atl,   Atr],
            [Abl, 2*Abr]
        ])
        return A


    def cr3bp(self, t, state):
        x, y, z, dx, dy, dz = state

        derivs = self.potentials(state)
        Ux, Uy, Uz = derivs.Ux, derivs.Uy, derivs.Uz

        ddx = 2*dy + Ux
        ddy = Uy - 2*dx
        ddz = Uz

        return [dx, dy, dz, ddx, ddy, ddz]
    
    def cr3bp_stm(self, t, state):
        X = state[:6]
        Phi = state[6:].reshape(6,6)
        dstate = self.cr3bp(t, X)

        A = self.A_matrix(X)

        dPhi = A @ Phi

        return np.concatenate((dstate, dPhi.flatten()))
    
    def jacobi_constant(self, state):
        x, y, z, dx, dy, dz = state[:6]
        r13, r23 = self.distances(state)
        mu = self.mu
        mu2 = 1 - mu
        U = mu2/r13 + mu/r23 + 1/2*(x**2 + y**2)
        C = 2*U - (dx**2 + dy**2 + dz**2)
        return C
    
    def energy(self, state):
        x, y, z, dx, dy, dz = state[:6]
        mu = self.mu
        mu2 = 1 - mu
        r13, r23 = self.distances(state)

        U = mu2/r13 + mu/r23 + 1/2*(x**2 + y**2)
        E = (1/2)*(dx**2 + dy**2 + dz**2) - U
        
        return E
    
    def libration_points(self):
        mu = self.mu
        def get_L_points(mu):
            """
            Gets the collinear Lagrangian points L1, L2, and L3
            Using Newton-Raphson Method to find points
            
            Parameters
            ---------------
            mu: float
                Normalized mass parameter
            
            Returns:
            ---------------
            sol_x: ndarray
                Non-dimensionalized x-axis locations of collinear lagrangian points L1, L2, and L3
            """
            def L1(x):
                return x - (1 - mu) / (x + mu)**2 + mu / (x - 1 + mu)**2
            def dL1(x):
                return 2*(1 - mu) / (x + mu)**3 - 2*mu / (x+mu-1)**3 + 1
            
            def L2(x):
                return -x + (1 - mu) / (x + mu)**2 + mu / (x - 1 + mu)**2
            def dL2(x):
                return -2*(1 + mu) / (x + mu)**3 - 2*mu / (x+mu-1)**3 - 1
            
            def L3(x):
                return -x - (1 - mu) / (x + mu)**2 - mu / (x - 1 + mu)**2
            def dL3(x):
                return 2*(1 + mu) / (x + mu)**3 + 2*mu / (x+mu-1)**3 - 1

            iter = 100
            error = 1e-13
            xL1 = NR_method( 0.8, L1, dL1, error, iter)
            xL2 = NR_method( 1.2, L2, dL2, error, iter)
            xL3 = NR_method(-1.0, L3, dL3, error, iter)

            xL45 = 1/2 - mu
            yL45 = np.sqrt(3)/2

            Lp = np.array([
                [ xL1,     0, 0],
                [ xL2,     0, 0],
                [ xL3,     0, 0],
                [xL45,  yL45, 0],
                [xL45, -yL45, 0]
            ])

            return Lp

        def NR_method(guess, f, f_prime, error, iter):
            x = np.zeros((int(iter),1))
            x[0] = guess
            for idx, _ in enumerate(x):
                x[idx+1,0] = x[idx, 0] - f(x[idx, 0])/f_prime(x[idx, 0])
                if np.abs(f(x[idx+1, 0])) < error:
                    return x[idx+1, 0]
                if idx == (iter-1):
                    print("Newton-Raphson couldn't converge")
    
        return get_L_points(mu)
