import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy import cos, sin, sqrt

def distances(x, y, z, mu):
    r12 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r13 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    return r12, r13

def get_jacobi_constant(state, mu):
    x, y, z, dx, dy, dz = state[:6]
    r12, r13 = distances(x, y, z, mu)
    U = 0.5*(x**2 + y**2) + (1-mu)/r12 + mu/r13
    return 2*U - (dx**2 + dy**2 + dz**2)

def cr3bp_accel(t, X, mu):
    x, y, z, dx, dy, dz = X[:6]
    r12, r13 = distances(x, y, z, mu)

    # Potentials First Derivates
    Ux = x - (1 - mu)*(x + mu)/r12**3 - mu*(x + mu - 1)/r13**3
    Uy = y - (1 - mu)*y/r12**3 - mu*y/r13**3
    Uz = -(1 - mu)*z/r12**3 - mu*z/r13**3

    ddx =  2 * dy + Ux
    ddy = -2 * dx + Uy
    ddz =  Uz

    return dx, dy, dz, ddx, ddy, ddz

def cr3bp_stm(t, X, mu):
    x, y, z, dx, dy, dz = X[:6]
    Phi = X[6:].reshape(6,6)
    r12, r13 = distances(x, y, z, mu)
    
    # Potentials First Derivates
    Ux = x - (1 - mu)*(x + mu)/r12**3 - mu*(x + mu - 1)/r13**3
    Uy = y - (1 - mu)*y/r12**3 - mu*y/r13**3
    Uz = -(1 - mu)*z/r12**3 - mu*z/r13**3

    # Potentials Non-Zero Jacobian Elements 
    Uxx = 1 - (1-mu)/r12**3 + 3*(1-mu)*(x+mu)**2/r12**5 - mu/r13**3 + 3*mu*(x+mu-1)**2/r13**5
    Uyy = 1 - (1-mu)/r12**3 + 3*(1-mu)*y**2/r12**5      - mu/r13**3 + 3*mu*y**2/r13**5
    Uzz =   - (1-mu)/r12**3 + 3*(1-mu)*z**2/r12**5      - mu/r13**3 + 3*mu*z**2/r13**5
    Uxy =  3*((1-mu)*(x+mu)*y/r12**5 + mu*(x-1+mu)*y/r13**5)
    Uyx = Uxy
    Uxz =  3*z*((1 - mu)*(x + mu)/r12**5 + mu*(x - (1 - mu))/r13**5)
    Uyz =  3*y*z*((1 - mu)/r12**5 + mu/r13**5)
    Uzy = Uyz

    Omega = np.array([
        [ 0, 1, 0],
        [-1, 0, 0],
        [ 0, 0, 0]
    ])

    identity = np.eye(3)

    UXX = np.array([
        [Uxx, Uxy, Uxz],
        [Uyx, Uyy, Uyz],
        [Uxz, Uzy, Uzz]
    ])

    F = np.block([[np.zeros((3,3)), identity], [UXX, 2*Omega]])

    ddx =  2 * dy + Ux
    ddy = -2 * dx + Uy
    ddz =  Uz

    dPhi = F @ Phi

    dX = np.array([dx, dy, dz, ddx, ddy, ddz])

    dPhi_flatten = dPhi.flatten()

    return np.concatenate((dX, dPhi_flatten))

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
    sol_L1 = NR_method(0, L1, dL1, error, iter)
    sol_L2 = NR_method(0, L2, dL2, error, iter)
    sol_L3 = NR_method(0, L3, dL3, error, iter)

    return np.array([sol_L1, sol_L2, sol_L3])

def NR_method(guess, f, f_prime, error, iter):
    x = np.zeros(int(iter) + 1)
    x[0] = guess
    for idx in range(iter):
        x[idx+1] = x[idx] - f(x[idx])/f_prime(x[idx])
        if np.abs(f(x[idx+1])) < error:
            return x[idx+1]
        if idx == (iter-1):
            print("Newton-Raphson couldn't converge")

def get_jacobi(state, mu):
    x, y, z, dx, dy, dz = state[:6]
    r12, r13 = distances(x, y, z, mu)
    U = 0.5*(x**2 + y**2) + (1 - mu)/r12 + mu/r13
    C = 2*U - (dx**2 + dy**2 + dz**2)

    return C

def gammaL(Lpt, mu):
    Ls = get_L_points(mu)
    mu1 = -mu
    mu2 = 1 - mu
    if Lpt == 1:
        gamma = mu2 - Ls[0]
    elif Lpt == 2:
        gamma = mu2 - Ls[1]
    elif Lpt == 3:
        gamma = Ls[2] - mu1
    else:
        print("Please choose a collinear libration point.")
    return gamma

def halo_initial_state(mu, Lpt, Azlp, n):
    Az = Azlp

    gamma = gammaL(Lpt, mu)
    if Lpt == 1:
        won = 1
        primary = 1-mu
    elif Lpt == 2:
        won = -1
        primary = 1-mu
    elif Lpt == 3:
        won = 1
        primary = -mu
    else:
        print("Please choose a collinear libration point.")
    
    c = np.zeros(4)

    if Lpt == 3:
        for N in range(1,4):
            c[N] = (1/gamma**3) * (1 - mu + (-primary * gamma**(N+1)) / ((1 + gamma)**(N+2)))
    else:
        for N in range(1,4):
             c[N] = (1/gamma**3) * ((won**(N+1)) * mu + ((-1)**(N+1)) * (primary * gamma**(N+2)) / ((1 + (-won)*gamma)**(N+2)))
    polylambda = [1, 0, (c[1]-2), 0, -(c[1]-1)*(1+2*c[1])]
    lamb = np.roots(polylambda)

    if Lpt == 3:
        lamb = abs(lamb[2])
    else:
        lamb = abs(lamb[0])
    
    k = 2*lamb/(lamb**2 + 1 - c[1])
    delta = lamb**2 - c[1]

    d1 = ((3*lamb**2)/k)*(k*(6*lamb**2 -1)-2*lamb)
    d2 = ((8*lamb**2)/k)*(k*(11*lamb**2-1)-2*lamb)

    a21 = (3*c[2]*(k**2-2))/(4*(1+2*c[1]))
    a22 = 3*c[2]/(4*(1+2*c[1]))
    a23 = -(3*c[2]*lamb/(4*k*d1))*(3*(k**3)*lamb - 6*k*(k-lamb) + 4)
    a24 = -(3*c[2]*lamb/(4*k*d1))*(2 + 3*k*lamb)

    b21 = -(3*c[2]*lamb/(2*d1))*(3*k*lamb-4)
    b22 = 3*c[2]*lamb/d1
    d21 = -c[2]/(2*lamb**2)

    a31 = -(9*lamb/(4*d2))*(4*c[2]*(k*a23 - b21) + k*c[3]*(4 + k**2)) + ((9*lamb**2 + 1 -c[1])/(2*d2))*(3*c[2]*(2*a23 - k*b21) + c[3]*(2 + 3*k**2))
    a32 = -(1/d2)*( (9*lamb/4)*(4*c[2]*(k*a24 - b22) + k*c[3]) + 1.5*(9*lamb**2 + 1 - c[1])*( c[2]*(k*b22 + d21 - 2*a24) - c[3]) )

    b31 = (.375/d2)*( 8*lamb*(3*c[2]*(k*b21 - 2*a23) - c[3]*(2 + 3*k**2)) + (9*lamb**2 + 1 + 2*c[1])*(4*c[2]*(k*a23 - b21) + k*c[3]*(4 + k**2)) )
    b32 = (1/d2)*( 9*lamb*(c[2]*(k*b22 + d21 - 2*a24) - c[3]) + .375*(9*lamb**2 + 1 + 2*c[1])*(4*c[2]*(k*a24 - b22) + k*c[3]) )

    d31 = (3/(64*lamb**2))*(4*c[2]*a24 + c[3])
    d32 = (3/(64*lamb**2))*(4*c[2]*(a23 - d21) + c[3]*(4 + k**2))

    s1 = (1/(2*lamb*(lamb*(1+k**2) - 2*k)))*( 1.5*c[2]*(2*a21*(k**2 - 2)-a23*(k**2 + 2) - 2*k*b21) - .375*c[3]*(3*k**4 - 8*k**2 + 8) )
    s2 = (1/(2*lamb*(lamb*(1+k**2) - 2*k)))*( 1.5*c[2]*(2*a22*(k**2 - 2)+a24*(k**2 + 2) + 2*k*b22 + 5*d21) + .375*c[3]*(12 - k**2) )

    a1 = -1.5*c[2]*(2*a21+ a23 + 5*d21) - .375*c[3]*(12-k**2)
    a2 = 1.5*c[2]*(a24-2*a22) + 1.125*c[3]

    l1 = a1 + 2*(lamb**2)*s1
    l2 = a2 + 2*(lamb**2)*s2

    tau1 = 0
    deltan =-n
    Ax = sqrt((-delta - l2*Az**2)/l1)

    x = a21*Ax**2 + a22*Az**2 - Ax*cos(tau1) + (a23*Ax**2 - a24*Az**2)*cos(2*tau1) + (a31*Ax**3 - a32*Ax*Az**2)*cos(3*tau1)
    y = k*Ax*sin(tau1) + (b21*Ax**2 - b22*Az**2)*sin(2*tau1) + (b31*Ax**3 - b32*Ax*Az**2)*sin(3*tau1)
    z = deltan*Az*cos(tau1) + deltan*d21*Ax*Az*(cos(2*tau1) - 3) + deltan * (d32*Az*Ax**2 - d31*Az**3)*cos(3*tau1)
    vx = lamb*Ax*sin(tau1) - 2*lamb*(a23*Ax**2 - a24*Az**2)*sin(2*tau1) - 3*lamb*(a31*Ax**3 - a32*Ax*Az**2)*sin(3*tau1)
    vy = lamb*(k*Ax*cos(tau1) + 2*(b21*Ax**2 - b22*Az**2)*cos(2*tau1) + 3*(b31*Ax**3 - b32*Ax*Az**2)*cos(3*tau1))
    vz = - lamb*deltan*Az*sin(tau1) - 2*lamb*deltan*d21*Ax*Az*sin(2*tau1) - 3*lamb*deltan*(d32*Az*Ax**2 - d31*Az**3)*sin(3*tau1)

    r0 = np.array([(primary + gamma*(-won + x))/gamma, -y, z])*gamma
    v0 = np.array([vx, vy, vz])*gamma
    x0 = np.concatenate((r0,v0))

    return x0

def event(t, x0, *args):
    return x0[1]

def y_crossing(t_guess, x0, mu):
    sol = solve_ivp(cr3bp_accel, (0, t_guess), x0, args=(mu,), rtol=1e-10, atol=1e-10, events=[event])
    tc = sol.t_events[0][1]
    Xc = sol.y_events[0][1]
    return tc, Xc

def single_shooting(t_guess, x0, mu):
    max_attempt = 50
    dx = 1
    dz = 1
    for _ in range(max_attempt):
        if abs(dx) < 1e-8 and abs(dz) < 1e-8:
            return x0, tc
        x0 = np.concatenate((x0[:6], np.eye(6).flatten()))
        tc, Xc = y_crossing(t_guess, x0[:6], mu)
        sol = solve_ivp(cr3bp_stm, (0, tc), x0, args=(mu,), rtol=1e-10, atol=1e-10)
        x, y, z, dx, dy, dz = Xc[:6]
        phi = sol.y[6:, -1].reshape(6,6)
        _, _, _, ddx, ddy, ddz = cr3bp_accel(tc, Xc[:6], mu)
        D1 = np.array([
            [phi[3,0], phi[3,4]],
            [phi[5,0], phi[5,4]]
        ])

        D2 = D1 - (1/dy)*np.outer([ddx, ddz], [phi[1,0], phi[1,4]])
        D3 = np.linalg.inv(D2) @ np.array([-dx, -dz])
        x0[0] += D3[0]
        x0[4] += D3[1]
    raise ValueError("Single Shooting Method couldn't converge.")

def plot_halos(mu, Lpt=1, A_min=1e-4, A_max = 0.6, n_pts=30, t_guess=1.5):
    A_vals = np.linspace(A_min, A_max, n_pts)
    Lx = get_L_points(mu)[Lpt-1]

    fig, ax = plt.subplots(subplot_kw={"projection":"3d"})

    for idx, Az in enumerate(A_vals):
        try:
            state0 = halo_initial_state(mu, Lpt, Az, -1)
            X0 = np.concatenate((state0, np.eye(6).flatten()))
            X0, tc = single_shooting(t_guess, X0, mu)
            sol = solve_ivp(cr3bp_accel, (0, 2*tc), X0[:6], args=(mu,), rtol=1e-10, atol=1e-10, events=[event])

            x, y, z = sol.y[:3, :]
            if idx % 10 == 0:
                ax.plot(x, y, z, color="black", alpha=0.4)
            else:
                ax.plot(x, y, z, color="blue", alpha=0.2)
            
        except Exception as e:
            print(f"Failed for Az={Az:.3f}: {e}")
    
    ax.set_xlabel("X (non-dim)")
    ax.set_ylabel("Y (non-dim)")
    ax.set_zlabel("Z (non-dim)")
    ax.set_title(f"Halo Family @ $L_{Lpt}$ (μ= {mu})")
    ax.scatter(1-mu, 0, 0, color="grey", label="Moon")
    ax.scatter(Lx, 0, 0, color="red", marker="*", label=r"$L_{1}$")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    mu = 1.215057e-2
    plot_halos(mu, Lpt=1, t_guess=2, n_pts=100)

