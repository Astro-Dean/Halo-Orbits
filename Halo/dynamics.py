import numpy as np
from scipy.integrate import solve_ivp

class main_dynamics:
    """
    Holds both the propagation of the CR3BP system and the dynamical system itself
    """

    def __init__(self, mu, ic, tf, t_eval=None, track_cross=0, stm_bool=False, tol=1e-12, model="DOP853"):
        self.mu = mu
        self.ic = ic
        self.tf = tf
        self.t_eval = t_eval
        self.track_cross = track_cross
        self.stm_bool = stm_bool
        self.tol = tol
        self.model = model
        self.prop_results = None

    def propagate(self):
        """Propagate trajectory using solve_ivp"""
        
        # Prepare initial conditions
        if self.stm_bool:
            if self.ic.size == 6:
                ic_with_stm = np.concatenate([self.ic, np.eye(6).flatten()])
            elif self.ic.size == 42:
                ic_with_stm = self.ic
            else:
                raise ValueError("Initial conditions must be length 6 or 42")
        else:
            if self.ic.size != 6:
                raise ValueError("When stm_bool is False, ic must be length 6")
            ic_with_stm = self.ic
        
        # Setup event function
        def y_event_crossing(t, y):
            return y[1]
        
        if self.track_cross == 0:
            events = None
        elif self.track_cross == 1:
            y_event_crossing.terminal = True
            y_event_crossing.direction = 1
            events = y_event_crossing
        else:
            events = None
        
        # Integrate
        t_start = 1e-6 if self.track_cross == 1 else 0.0
        self.prop_results = solve_ivp(
            fun=self.cr3bp,
            t_span=(t_start, self.tf),
            y0=ic_with_stm,
            method=self.model,
            events=events,
            rtol=self.tol,
            atol=self.tol,
            dense_output=False
        )
        
        return self.prop_results

    def cr3bp(self, t, state):
        """
        CR3BP equations of motion with optional STM

        Parameters:
            t: float
                time (required by solve_ivp)
            state: ndarray
                state vector (6 or 42 elements)

        Returns:
            dstate/dt: ndarray
                derivative vector (same length as state)
        """
        x, y, z, dx, dy, dz = state[:6]

        # Compute potential derivatives
        Ux, Uy, Uz, *_ = self.potentials(state[:6])

        ddx = 2*dy + Ux
        ddy = -2*dx + Uy
        ddz = Uz

        dX = np.array([dx, dy, dz, ddx, ddy, ddz])

        if self.stm_bool and len(state) > 6:
            Phi = state[6:].reshape(6, 6)

            A = self.get_A(state)

            dPhi = A @ Phi
            return np.concatenate([dX, dPhi.flatten()])
        else:
            return dX

    def potentials(self, state):
        """
        Compute CR3BP potential derivatives

        Parameters:
            state: 6-element state vector

        Returns:
            Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz
        """
        mu1 = 1 - self.mu
        mu2 = self.mu

        x, y, z = state[:3]

        r1, r2 = self.distances(state)

        Ux = x - mu1*(x + mu2)/r1**3 - mu2*(x - mu1)/r2**3
        Uy = y - mu1*y/r1**3 - mu2*y/r2**3
        Uz = -mu1*z/r1**3 - mu2*z/r2**3

        Uxx = 1 - mu1/r1**3 - mu2/r2**3 + 3*mu1*(x + mu2)**2/r1**5 + 3*mu2*(x - mu1)**2/r2**5
        Uyy = 1 - mu1/r1**3 - mu2/r2**3 + 3*mu1*y**2/r1**5 + 3*mu2*y**2/r2**5
        Uzz = -mu1/r1**3 - mu2/r2**3 + 3*mu1*z**2/r1**5 + 3*mu2*z**2/r2**5

        Uxy = 3*(mu1*(x + mu2)/r1**5 + mu2*(x - mu1)/r2**5)*y
        Uxz = 3*(mu1*(x + mu2)*z/r1**5 + mu2*(x - mu1)*z/r2**5)
        Uyz = 3*(mu1*y*z/r1**5 + mu2*y*z/r2**5)

        return Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz

    def get_A(self, state):
            Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz = self.potentials(state[:6])
            A = np.array([
                [0,    0,    0,    1,  0,  0],
                [0,    0,    0,    0,  1,  0],
                [0,    0,    0,    0,  0,  1],
                [Uxx,  Uxy,  Uxz,  0,  2,  0],
                [Uxy,  Uyy,  Uyz, -2,  0,  0],
                [Uxz,  Uyz,  Uzz,  0,  0,  0]
            ])
            return A

    def distances(self, state):
        """Compute distances to primaries"""
        x, y, z = state[:3]
        r1 = np.sqrt((x + self.mu)**2 + y**2 + z**2)
        r2 = np.sqrt((x - 1 + self.mu)**2 + y**2 + z**2)
        return r1, r2

    def store_info(self):
        """Store propagation results"""
        prop_info = {}

        t = self.prop_results.t
        Y = self.prop_results.y

        prop_info["t"] = t
        prop_info["states"] = Y
        prop_info["final_state"] = Y[:, -1]

        if self.stm_bool:
            Phi_hist = Y[6:, :].T.reshape(len(t), 6, 6)
            prop_info["stm"] = Phi_hist
            prop_info["Phi_tf"] = Phi_hist[-1]

        if hasattr(self.prop_results, 't_events') and self.prop_results.t_events:
            if len(self.prop_results.t_events) > 0:
                prop_info["t_events"] = self.prop_results.t_events[0]
                prop_info["y_events"] = self.prop_results.y_events[0]
        else:
            prop_info["t_events"] = None
            prop_info["y_events"] = None

        prop_info["success"] = bool(self.prop_results.success)
        prop_info["message"] = self.prop_results.message
        prop_info["nfev"] = self.prop_results.nfev

        return prop_info
