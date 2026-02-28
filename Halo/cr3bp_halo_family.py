"""
Halo family generation code

Created by Dean Hall
Date: February 9th, 2026
"""


import numpy as np
from scipy.integrate import solve_ivp
from systems import get_system_dim
from dynamics import main_dynamics
import scipy.linalg

class CR3BP_Halo(main_dynamics):
    """
    # Halo Family Generation
    This class allows for the generation of a halo orbit family in the CR3BP environment from the Lagrangian points L1 and L2 
    using Natural Parameter Continuation and Psuedo Arclength Continuation to solve for orbits along a solution curve. The only requirements
    for the class is the selection of a system, the initial guess state, and the assumed orbital period of the initial guess.

    
    """
    def __init__(self, system, ic, tf, t_eval=None, track_cross=0, stm_bool=False, tol=1e-12, model="DOP853"):
        """
        Initialization of values necessary for functionality of Halo family generation

        Parameters:
        -------------
        system : str
            Name of CR3BP system
        ic : ndarray
            Initial guess state (6 dimensions)
        tf : float
            Initial guess period 
        
        """
        self.tf_full = tf
        tf_half = tf / 2  # Period guess halved for halo convergence
        super().__init__(get_system_dim(system).mu, ic, tf_half, t_eval, track_cross, stm_bool, tol, model)
        self.prop_results = None
    
    def single_shooting(self, max_iter=50, tol=1e-12, damping=1.0, palc_args=None):
        # Validate/fix initial conditions
        if abs(self.ic[1]) > 1e-6:
            print(f"Warning: y0={self.ic[1]:.3e}, setting to 0")
            self.ic[1] = 0.0
        if abs(self.ic[3]) > 1e-6:
            print(f"Warning: vx0={self.ic[3]:.3e}, setting to 0")
            self.ic[3] = 0.0
        if abs(self.ic[5]) > 1e-6:
            print(f"Warning: vz0={self.ic[5]:.3e}, setting to 0")
            self.ic[5] = 0.0
        
        self.stm_bool = True

        # Free variables: [x0, z0, dy0, T]
        X = np.array([self.ic[0], self.ic[2], self.ic[4], self.tf])
        X_prev_iter = X.copy()
        
        using_palc = palc_args is not None
        
        if using_palc:
            print("Starting standard single shooter with PALC")
        else:
            print("Starting standard single shooter (no PALC)")
        
        print(f"Free vars: [x0, z0, dy0, T]")
        print(f"Constraints: [y(T)=0, vx(T)=0, vz(T)=0]")
        print(f"Initial X: {X}")
        print(f"Initial IC: {self.ic}")
        print(f"Initial JC: {self._compute_jacobi_constant(self.ic):.6f}\n")

        best_error = np.inf
        best_X = X.copy()
        self.DF = None

        for iteration in range(max_iter):
            # Update state from X
            self.ic[0] = X[0]
            self.ic[1] = 0.0
            self.ic[2] = X[1]
            self.ic[3] = 0.0
            self.ic[4] = X[2]
            self.ic[5] = 0.0
            self.tf = X[3]

            # Propagate
            try:
                self.propagate()
                prop_info = self.store_info()
            except Exception as e:
                print(f"Integration error: {e}")
                X = X_prev_iter.copy()
                damping *= 0.5
                if damping < 0.01:
                    return None, False
                continue
            
            if not prop_info['success']:
                print(f"Integration failed")
                X = X_prev_iter.copy()
                damping *= 0.5
                if damping < 0.01:
                    return None, False
                continue
            
            # Extract final state and STM
            x_final = prop_info['final_state'][:6]
            Phi = prop_info['Phi_tf']
            
            # Constraint vector: F = [y(T), vx(T), vz(T)]
            F = np.array([
                x_final[1],  # y(T)
                x_final[3],  # vx(T)
                x_final[5]   # vz(T)
            ])
            
            # Add PALC constraint if using PALC
            if using_palc:
                palc_constraint = np.dot(X - palc_args['X_prev'], palc_args['tangent']) - palc_args['delta_s']
                F_augmented = np.append(F, palc_constraint)
            else:
                F_augmented = F
            
            # Check convergence
            error = np.linalg.norm(F)
            print(f"Iteration {iteration}: |F| = {error:.3e}, damping = {damping:.3f}")
            
            if error < best_error:
                best_error = error
                best_X = X.copy()
            
            if error < tol:
                print(f"\nConverged in {iteration} iterations!")
                print(f"Final IC: {self.ic}")
                print(f"Half period: {self.tf:.6f}")
                print(f"Full period: {2*self.tf:.6f}")
                print(f"JC: {self._compute_jacobi_constant(self.ic):.6f}")
                print(f"z₀: {self.ic[2]:.6f}")

                Phi_half = prop_info['Phi_tf']
                S = np.diag([1,-1,1,-1,1,-1])
                M = S @ np.linalg.inv(Phi_half) @ S @ Phi_half
                stability_index = self.stability(M)

                result = {
                    'ic': self.ic.copy(),
                    'tf': self.tf,
                    'states': prop_info['states'][:6, :].T,
                    't': prop_info['t'],
                    'stm': prop_info['stm'],
                    'stability': stability_index,
                    'free_vars': X.copy(),
                    'DF': self.DF if self.DF is not None else np.zeros((3,4)),
                    'jc': self._compute_jacobi_constant(self.ic),
                    'converged': True,
                    'iterations': iteration
                }
                return result, True
            
            if error > 1e6 or np.isnan(error):
                print("Diverging")
                X = best_X.copy()
                damping *= 0.5
                if damping < 0.01:
                    return None, False
                continue
            
            # Build Jacobian (3x4)
            # Constraints: [y, vx, vz] at final time
            # Free vars: [x0, z0, vy0, T]
            
            # Time derivatives
            vy_final = x_final[4]
            
            # Accelerations at final state
            Ux, Uy, Uz, _, _, _, _, _, _ = self.potentials(x_final)
            ax_final = 2*vy_final + Ux
            az_final = Uz
            
            # Build Jacobian     
            self.DF = np.array([
                [Phi[1,0], Phi[1,2], Phi[1,4], vy_final],
                [Phi[3,0], Phi[3,2], Phi[3,4], ax_final],
                [Phi[5,0], Phi[5,2], Phi[5,4], az_final]
            ])
            
            # Augment Jacobian if using PALC
            if using_palc:
                DF_augmented = np.vstack([self.DF, palc_args['tangent']])
            else:
                DF_augmented = self.DF
            
            # Newton-Raphson update
            X_prev_iter = X.copy()
            delta_X = np.linalg.pinv(DF_augmented) @ F_augmented
            X = X - damping * delta_X
            
            # Adaptive damping
            if error < best_error:
                damping = min(1.0, damping * 1.1)
        
        print(f"\nFailed to converge in {max_iter} iterations")
        print(f"Best error: {best_error:.3e}")
        return None, False
    
    def ycrossing_events(self):
        self.track_cross = 1
        results = self.propagate()
        if len(results["tevents"][0][:]) > 1:
            self.tf = results["tevents"][0][1]
        else:
            self.tf = results["tevents"][0][0]
    
    def _compute_tangent(self, DF):
        """
        Compute tangent direction from Jacobian null space
        
        For a 3x4 Jacobian, null space is 1D

        Paraneters:
        -------------
        DF: ndarray
            Jacobian of constraint vector F
        
        Returns:
        -------------
        tangent: ndarray
            Tangent along the solution curve
        """
        
        null_space = scipy.linalg.null_space(DF)
        
        if null_space.shape[1] != 1:
            print(f"Warning: Null space dimension = {null_space.shape[1]}, expected 1")
            tangent = null_space[:, 0]
        else:
            tangent = null_space.flatten()
        
        # Normalize
        tangent = tangent / np.linalg.norm(tangent)
        
        return tangent

    def _compute_jacobi_constant(self, state):
        x, y, z, dx, dy, dz = state[:6]
        r1, r2 = self.distances(state)
        U = 0.5 * (x**2 + y**2) + (1 - self.mu)/r1 + self.mu/r2
        V = dx**2 + dy**2 + dz**2
        return 2*U - V
    
    def npc(self, z_target=0.03, step=0.02, max_members=20,
            max_iter=50, tol=1e-12, damping=1.0, verbose=True):
        if verbose:
            print("=" * 60)
            print("Natural Parameter Continuation (NPC) in z0")
            print(f"  Step size : {step:.4f}")
            print(f"  Target z0 : {z_target:.4f}")
            print("=" * 60)
            print("\nConverging seed orbit...")

        result, converged = self.single_shooting(max_iter=max_iter, tol=tol, damping=damping)
        if not converged:
            raise RuntimeError("NPC: failed to converge seed orbit")
        
        if verbose:
            print(f"Seed Converged: z0={result['ic'][2]:.6e}, "
                  f"Period={2*result['tf']:.4f}")
        
        family = [result]

        for i in range(1, max_members + 1):
            prev_ic = family[-1]['ic'].copy()
            prev_tf = family[-1]['tf']

            new_z0 = prev_ic[2] + step
            if new_z0 >= z_target:
                if verbose:
                    print(f"Reached target z0={z_target:.4f}")
                break
            
            self.ic = prev_ic.copy()
            self.ic[2] = new_z0
            self.tf = prev_tf

            result, converged = self.single_shooting(max_iter=max_iter, tol=tol, damping=damping)

            if not converged:
                if verbose:
                    print(f"Not converged, z0={new_z0:.6f}, halving step")
                    step *= 0.5
                if step <1e-4:
                    print(f"Step reduced to 1e-4. Stopping NPC")
                    break
                continue

            family.append(result)

            z0_c = result['ic'][2]
            z_amp = np.max(np.abs(result['states'][:,2]))
            if verbose:
                print(f"Member {len(family)-1:3d}: z0={z0_c:.4f}, "
                      f"Z={z_amp:.4f}, Period={2*result['tf']:.4f}, "
                      f"iters={result['iterations']}")
            
            if result['iterations'] < 4:
                step = min(0.005, step * 1.2)
        
        if verbose:
            print(f"\nNPC complete: {len(family)} members")
            print(f"Final z0: {family[-1]['ic'][2]:.4f}")
        
        return family
        
    def palc(self, family, z_target=2.5, step=0.005, max_members=300,
        max_iter=50, tol=1e-12, damping=1.0, verbose=True):
        if verbose:
            print("=" * 60)
            print("Pseudo Arc-Length Continuation (PALC)")
            print(f"  Step size     : {step:.4f}")
            print(f"  Target Az     : {z_target:.4f}")
            print(f"  Max members   : {max_members}")
            print("=" * 60 + "\n")

        tangent_prev = None

        for i in range(max_members):
            prev_result = family[-1]
            X_prev = prev_result['free_vars']
            DF_prev = prev_result['DF']

            # Compute tangent to the solution curve
            tangent = self._compute_tangent(DF_prev)

            # Ensure consistent direction
            if tangent_prev is None:
                if tangent[1] < 0:
                    tangent *= -1
            else:
                if np.dot(tangent, tangent_prev) < 0:
                    tangent *= -1

            tangent_prev = tangent.copy()

            # Predictor step
            X_predicted = X_prev + step * tangent

            # Update shooter state from predicted free variables
            # free_vars = [x0, z0, ydot0, tf]
            self.ic[0] = X_predicted[0]
            self.ic[1] = 0.0
            self.ic[2] = X_predicted[1]
            self.ic[3] = 0.0
            self.ic[4] = X_predicted[2]
            self.ic[5] = 0.0
            self.tf    = X_predicted[3]

            palc_args = {
                'delta_s' : step,
                'X_prev'  : X_prev,
                'tangent' : tangent,
            }

            result, converged = self.single_shooting(
                max_iter=max_iter, tol=tol, damping=damping, palc_args=palc_args
            )

            if not converged:
                if verbose:
                    print(f"Step {i} failed - halving step to {step*0.5:.4f}")
                step *= 0.5
                if step < 1e-4:
                    if verbose:
                        print("Step reduced to 1e-4 - stopping PALC.")
                    break
                continue

            family.append(result)

            z0   = result['ic'][2]
            z_amp = np.max(np.abs(result['states'][:, 2]))
            if verbose:
                print(f"Member {len(family)-1:3d}: z0={z0:.4f}, "
                      f"Az={z_amp:.4f}, Period={2*result['tf']:.4f}, "
                      f"iters={result['iterations']}")

            # Adaptive step
            if result['iterations'] < 5:
                step = min(0.01, step * 1.2)

            if z_amp >= z_target:
                if verbose:
                    print(f"\nReached target z-amplitude ({z_target:.4f}).")
                break

        if verbose:
            print(f"\nPALC complete: {len(family)} total family members")

        return family

    def stability(self, stm):
        eigvals, _ = np.linalg.eig(stm)
        return (max(np.abs(eigvals)))
