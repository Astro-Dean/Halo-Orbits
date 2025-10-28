import numpy as np
from numpy import sqrt, cos, sin
from halo_dynamics import Halo_Orbit_Dynamics
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class HaloOrbitFamily:
    def __init__(self, mu, Lpt, branch):
        self.mu = mu
        self.Lpt = Lpt
        self.branch = branch
    
    def third_order_richardson(self, Az):
        n = self.branch
        mu = self.mu
        Lpt = self.Lpt

        orbit = Halo_Orbit_Dynamics(mu)

        def gammaL(Lpt, mu):
            Ls = orbit.libration_points()
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

        gamma = gammaL(self, Lpt, mu)
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

        omega = 1 + s1*Ax**2 + s2*Az**2
        T0 = 2*np.pi / omega

        state0 = np.array([x, y, z, vx, vy, vz])

        return state0, T0
    
    def y_crossing(self, t_guess, x0, mu):
        def event(t, x0, *args): return x0[1]
        cr3bp_accel = Halo_Orbit_Dynamics.cr3bp
        sol = solve_ivp(cr3bp_accel, (0, t_guess), x0, args=(mu,), rtol=1e-10, atol=1e-10, events=[event])
        tc = sol.t_events[0][1]
        Xc = sol.y_events[0][1]
        return tc, Xc
    
    def single_shooting(self, t_guess, x0, mu, case): 
        max_attempt = 50 
        dx = 1 
        dz = 1 
        y = 1
        cr3bp_stm = Halo_Orbit_Dynamics.cr3bp_stm
        cr3bp_accel = Halo_Orbit_Dynamics.cr3bp
        for _ in range(max_attempt): 
            
            if case == 2:
                F = np.array([y, dx, dz])
                if np.linalg.norm(F, np.inf) < 1e-10: return x0, tc
            else:
                if abs(dx) < 1e-8 and abs(dz) < 1e-8: return x0, tc
            
            x0 = np.concatenate((x0[:6], np.eye(6).flatten())) 
            tc, Xc = self.y_crossing(t_guess, x0[:6], mu)
            sol = solve_ivp(cr3bp_stm, (0, tc), x0, args=(mu,), rtol=1e-12, atol=1e-12, method='DOP853') 
            x, y, z, dx, dy, dz = Xc[:6] 
            phi = sol.y[6:, -1].reshape(6,6) 
            _, _, _, ddx, ddy, ddz = cr3bp_accel(tc, Xc[:6], mu)
            if case == 1:
                D1 = np.array([ [phi[3,0], phi[3,4]], [phi[5,0], phi[5,4]] ]) 
                D2 = D1 - (1/dy)*np.outer([ddx, ddz], [phi[1,0], phi[1,4]]) 
                D3 = np.linalg.inv(D2) @ np.array([-dx, -dz]) 
                x0[0] += D3[0] 
                x0[4] += D3[1]
            elif case == 2:
                D1 = np.array([
                    [phi[1,2], phi[1,4], dy],
                    [phi[3,2], phi[3,4], ddx],
                    [phi[5,2], phi[5,4], ddz]
                ])
                F = np.array([y, dx, dz])
                delta = -np.linalg.solve(D1,F)
                if abs(dy) < 1e-8:
                    scale = .5
                else:
                    scale = 1
                x0[2]   += scale*delta[0]     # z0
                x0[4]   += scale*delta[1]     # ydot0
                t_guess += scale*delta[2]     # T/2
            else:
                raise ValueError("case value must 1 or 2.")
            
        raise ValueError("Single Shooting Method couldn't converge.")

    def halo_family(self, mu, Lpt, branch, n_steps):
        t_max = 5
        case = Lpt
        ds = 1.0
        
        if Lpt == 1:
            amps = [1e-4*branch, 1e-3*branch]
        elif Lpt == 2:
            amps = [1e-4*-branch, 1e-1*-branch]
        else:
            raise ValueError("Please choose lagrangian Point 1 or 2")
        
        family_states = []
        family_periods = []

        for i in range(2):
            state0 = self.third_order_richardson(mu, Lpt, amps[i], -1)
            x0 = np.concatenate((state0, np.eye(6).flatten()))
            X0, tc = self.single_shooting(t_max, x0, mu, case)
            family_states.append(X0[:6])
            family_periods.append(tc)
        
        for k in range(2, n_steps+2):
            X_prev1, X_prev2 = family_states[-2], family_states[-1]
            X_predict = X_prev2 + ds * (X_prev2 - X_prev1)

            x_predict = np.concatenate((X_predict, np.eye(6).flatten()))

            try:
                X_corr, tc_corr = self.single_shooting(t_max, x_predict, mu, case)
                family_states.append(X_corr[:6])
                family_periods.append(tc_corr)
                print(f"Step {k}: success, ds={ds:.3f}")
                if k > 90:
                    print(X_corr[:6])
            except Exception as e:
                ds *= 0.5
                print(f"Step {k}: correction failed — {e}, lowering ds to {ds:.3f}")
        
        return family_states, family_periods
    
    def plot_halo_family(self, mu, Lpt, branch):
        cr3bp_stm = Halo_Orbit_Dynamics.cr3bp_stm
        get_L_points = Halo_Orbit_Dynamics.libration_points
        
        if Lpt == 1:
            Lp = r"L$_{1}$"
            n_steps = 1500 # Gets slower with increasing steps but will succeed without changes to ds
        elif Lpt == 2:
            Lp = r"L$_{2}$"
            n_steps = 86 # Gets very slow after step 88 which has a failure when ds = 1.0
        else:
            raise ValueError("Please choose lagrangian Point 1 or 2")
        
        states, periods = self.halo_family(mu, Lpt, branch, n_steps)

        colors = plt.cm.plasma(np.linspace(0, 1, n_steps+2))
        fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
        plt.rc('text.latex', preamble=r'\usepackage{textgreek}')
        for i, (X0, tc) in enumerate(zip(states, periods)):
            sol = solve_ivp(cr3bp_stm, (0, 2*tc), np.concatenate((X0, np.eye(6).flatten())),
                        args=(mu,), rtol=1e-12, atol=1e-12, method='DOP853')
            x, y, z = sol.y[:3, :]
            ax.plot(x, y, z, color=colors[i], lw=1.5)
        ax.scatter(1 - mu, 0, 0, color='gray', s=40, label='Secondary')
        ax.scatter(get_L_points(mu)[Lpt-1], 0, 0, color='red', marker='*', s=60, label=Lp)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(rf"Halo Orbit Family about {Lp}")
        ax.legend()
        plt.show()