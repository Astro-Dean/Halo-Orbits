import numpy as np

def lib_gen(mu):
    """
    Finds and returns collinear lagrangian point locations with L4 and L5
    
    Parameters:
    -----------
        mu : float
            Gravitational mass parameter for CR3BP system
    
    Returns:
    -----------
        xlocs : ndarray
            x-axis locations of lagrangian points in the system
        ylocs : ndarray
            y-axis locations of lagrangian points in the system
    """
    def L1(x):
        return x - (1 - mu) / (x + mu)**2 + mu / (x - 1 + mu)**2
    def dL1(x):
        return 2*(1 - mu) / (x + mu)**3 - 2*mu / (x+mu-1)**3 + 1
    
    def L2(x):
        return -x + (1 - mu) / (x + mu)**2 + mu / (x - 1 + mu)**2
    def dL2(x):
        return -2*(1 - mu) / (x + mu)**3 - 2*mu / (x+mu-1)**3 - 1
    
    def L3(x):
        return -x - (1 - mu) / (x + mu)**2 - mu / (x - 1 + mu)**2
    def dL3(x):
        return 2*(1 - mu) / (x + mu)**3 + 2*mu / (x+mu-1)**3 - 1

    
    def newton(f, f_prime, guess, error = 1e-13, iter = 100):
        """Newton's method used to find collinear x locations"""
        x = np.zeros(int(iter) + 1)
        x[0] = guess
        for idx in range(iter):
            x[idx+1] = x[idx] - f(x[idx])/f_prime(x[idx])
            if np.abs(f(x[idx+1])) < error:
                return x[idx+1]
            if idx == (iter-1):
                print("Newton-Raphson couldn't converge")

    
    xL1 = newton(L1, dL1, 0)
    xL2 = newton(L2, dL2, 0)
    xL3 = newton(L3, dL3, 0)
    
    xlocs = np.array([xL1, xL2, xL3, 1/2 - mu, 1/2 - mu], dtype = np.float64)
    ylocs = np.array([0, 0, 0, np.sqrt(3)/2, -np.sqrt(3)/2], dtype = np.float64)
    return xlocs, ylocs

