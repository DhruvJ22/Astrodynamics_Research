"""
Created on Jan 8 21:42:49 2022
Updated 20 Mar 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Objective: This file contains a Class that will serve as the Base Class for any CR3BP computation
            It contains method to propagate a state in CR3BP and compute the Jacobi Constant

    Features:
        1. Integrate CR3BP EOMs
        2. Integrate CR3BP EOMs + Compute the State Transition Matrix for each state
        3. Optional events function is added to numerical integrator to track when states move from -y to +y region and/or visa versa
        4. Access solve_ivp->t_eval to save integrated results at set time stamps
        5. Computes accelration, first-derivate of pseudo-potential, and second-derivative of pseudo-potenital terms terms
        6. Compute Jacobi Constant [nd]
        
References
____________
This work heavily relies on the work done by the various past and current members of the Multi-Body Dynamics Research Group and Prof. Kathleen C. Howell
These are some of the referneces that provide a comprehensive brackground and have been the foundation for the work:
1. E. Zimovan, "Characteristics and Design Strategies for Near Rectilinear Halo Orbits Within the Earth-Moon System," M.S., August 2017
2. E. Zimovan Spreen, "Trajectory Design and Targeting for Applications to the Exploration Program in Cislunar Space," Ph.D., May 2021
3. V. Szebehely, "Theory of Orbits: The Restricted Problem of Three Bodies", 1967
4. W. Koon, M. Lo, J. Marsden, S. Ross, "Dynamical Systems, The Three-Body Problem, and Space Mission Design", 2006
"""
import numpy as np
from scipy.integrate import solve_ivp

class cr3bp_model:    
    """ Base class to investigate dynamics in CR3BP -> EOM propagator, JC calculation
    """
    
    def __init__(self, sys_chars_vals, ic=np.zeros((6)), tf=0, teval=None, int_tol=1e-12, stm_bool=0, xcross_cond=0, int_method='DOP853'):
        """
        Constructor
        
        Parameters
        ----------
        sys_chars_vals : object
            object of Class sys_char
        ic : numpy ndarray (6x1), {Can handle all 42 states for CR3BP+STM integration}
            States are defined about the barycenter of the two primaries, P1 and P2
            Initial condition: 6 states to compute a trajectory;
            [0:x0, 1:y0, 2:z0, 3:vx0, 4:vy0, 5:vz0] [non-dimensional] [nd]
            default = [0,0,0,0,0,0]
        tf : float
            Integration time [nd]
            Can be negative or positive, negative => Integration in backwards time, default = 0
        teval: list of float
            Contains time stamps at which the numerically integrated results will be saved
        int_tol : float
            Absolute = Relative Integration Tolerance
            The default is 1e-12.
        stm_bool : boolean, optional
            0: CR3BP EOM Integration
            1: CR3BP + STM EOM Integration
            The default is 0.
        xcross_cond : int, optional
            0 => No y-crossing check
            1 => Events function to check when crossed y axis in any direction
            2 => Events function to check when corssed y axis from -y to +y
            The default is 0.
        int_method : string, optional
            Specify integration scheme: 'DOP853' or 'LSODA'
            The default is 'DOP853'.
        """
        self.sys_chars_vals = sys_chars_vals
        self.mu = sys_chars_vals.mu
        self.ic = ic
        self.tf = tf
        self.teval = teval
        self.int_tol = int_tol
        self.stm_bool = stm_bool
        self.xcross_cond = xcross_cond
        self.int_method = int_method
        self.results = {}
        
    def propagate(self):
        """Numerically Integrate Circular Restricted Three-Body Problem EOMs
    
        Returns
        -------
          results : Dictionary
              't': time history, [nd]
              'states': state history, [:,i] => 6 states @ t = t_i, [nd]
              'yevents': states at prescribed event, [nd]
              'tevents': time stamp of yevents, [nd]
              'stm': STM element history, [:,:,i] => STM @ t = t_i
        """
        
        # Check to see if any state of I.C. is complex and then accordingly set the datatype
        # This is done to ease the implementation of Complex step derivative formulation to compute numerical partials in the future
        if all(np.imag(imag_check) == 0 for imag_check in self.ic):
            datatype = np.float64
            self.ic = np.real(self.ic)
        else:
            datatype = np.complex
            if self.int_method == "LSODA":
                print(self.int_method + " cannot integrate in complex domain")
                return 0    
        self.datatype = datatype
        
        # Check integration scheme
        if self.int_method != "DOP853" and self.int_method != "LSODA":
            print(
                "Please DROP853 or LSODA as the integration schemes, other schemes may not be compatible with the setup"
            )
            return 0
    
        # Accept 6 states and append 36 inital states if stm_bool = 1
        # Runs even if all 42 sates are given
        if len(self.ic) == 6 and self.stm_bool == 1:
            self.ic = np.concatenate(
                (self.ic, np.identity(6).flatten())
            )  # Appends the IC for STM: IC[6x1] + I_6x6
        elif len(self.ic) == 42:
            self.stm_bool = 1  # To make sure stm_bool is set to 1 if all 42 states are passed
        if len(self.ic) == 6:
            self.stm_bool == 0
        elif len(self.ic) != 6 and len(self.ic) != 42:
            print(
                "Initial conditions are neither of length 6 nor 42, recheck the input, len is "
                + str(len(self.ic))
            )
            return 0
    
        # Events function setup
        t0 = 0  # By default set initial integration time
    
        def xcross(t, y):
            """Track y position state for events function during integration"""
            return y[1]
    
        if self.xcross_cond == 1:  # Track events when crossising -y to +y region or +y to -y region, that is crossing XZ or XY plane
            xcross.terminal = False
            xcross.direction = 0
        elif self.xcross_cond == 2:  # Track events when crossing -y to +y region
            xcross.terminal = True
            xcross.direction = 1
        else:  # Does not track any events
            xcross = None

        # Numerical Integration
        fun = solve_ivp(self.__Nondim_DE_CR3BP_STM, [t0, self.tf], self.ic, method=self.int_method, t_eval = self.teval, events=xcross, rtol=self.int_tol, atol=self.int_tol)    
       
        # Save data to dictionary
        self.results = self.__save_prop_data_cr3p(fun, self.stm_bool, datatype)
    
        return self.results
    
    # Skeleton function for CR3BP + STM Numerical Integration
    def __Nondim_DE_CR3BP_STM(self, t, state_stm):
        """
        Describes the CR3BP EOM + STM and fed into a numerical integrator

        Parameters
        ----------
        t : flaot
            time @ which EOM are evaluated
        state_stm : numpy ndarray (6x1) or (42x1), float64/complex 128
             6 states + 36 elements of 6x6 STM
        Returns
        -------
        dstate_stm : numpy ndarray (6x1) or (42x1), float64/complex 128
            Time derivative of 6 states + 36 elements of 6x6 STM
        """
        
        dstate_stm = np.empty((len(state_stm),), dtype=self.datatype)
        
        dist_p1_p3, dist_p2_p3 = self.rel_dist_cr3bp(state_stm)  # Relative position vectors        
        _,_,_, ax, ay, az = self.ui_partials_acc_cr3bp(state_stm[0:6])
        
        # CR3BP: 6 states
        dstate_stm[0] = state_stm[3]
        dstate_stm[1] = state_stm[4]
        dstate_stm[2] = state_stm[5]
        dstate_stm[3] = ax
        dstate_stm[4] = ay
        dstate_stm[5] = az

        # STM: 36 States
        if len(state_stm) == 42:
            Uxx, Uyy, Uzz, Uxy, Uxz, Uyz = self.uii_partials_cr3bp(state_stm[0:6])# Hetian of  Pesudo-potenital of CR3BP
            At = np.array(
                [
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [Uxx, Uxy, Uxz, 0, 2, 0],
                    [Uxy, Uyy, Uyz, -2, 0, 0],
                    [Uxz, Uyz, Uzz, 0, 0, 0],
                ]
            )

            phi = np.zeros((6, 6), dtype=self.datatype)

            # Assingment of x[] for phi, STM
            count = 6
            for i in range(6):
                for j in range(6):
                    phi[i][j] = state_stm[count]
                    count = count + 1

            phi_d = np.matmul(At, phi)  # STM_dot = A(t)*STM

            dstate_stm[6:42] = phi_d.flatten()

        return dstate_stm
    
    
    def __save_prop_data_cr3p(self, fun, stm_bool, datatype):
        """Saves Numerical Integration results to a dictionary
    
        Parameters
        ----------
        fun : Object
            solve_ivp object
        stm_bool : boolean, optional
            0: CR3BP EOM Integration
            1: CR3BP + STM EOM Integration
        datatype : datatype
            Can be np.float64 or np.complex
    
        Returns
        -------
        results : Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
    
        """
        # Save the 6 states
        t = fun.t
        x = fun.y[0]
        y = fun.y[1]
        z = fun.y[2]
        vx = fun.y[3]
        vy = fun.y[4]
        vz = fun.y[5]
    
        # Save events function results
        tevents = fun.t_events
        yevents = fun.y_events
    
        # Save STM elements
        if stm_bool == 1:
            count = 0
            stm_vals = np.zeros((6, 6, len(t)), dtype=datatype)
            for t1 in range(6):
                for t2 in range(6):
                    for t3 in range(len(t)):
                        stm_vals[t1, t2, t3] = fun.y[count + 6][t3]
                    count = count + 1
    
        results = {}    
        results["t"] = np.asarray(t)
        results["states"] = np.concatenate(
            [
                [np.asarray(x)],
                [np.asarray(y)],
                [np.asarray(z)],
                [np.asarray(vx)],
                [np.asarray(vy)],
                [np.asarray(vz)],
            ]
        ).T
        if stm_bool == 1:
            results["stm"] = stm_vals
        results["tevents"] = tevents
        results["yevents"] = yevents
    
        return results
    
    
    def rel_dist_cr3bp(self, state=None):
        """Compute distance between a satellite(P3) defined in P1-P2 barycenter
            to P1 and P2 in CR3BP
    
        Parameters
        ----------
        states : numpy ndarray (6x1)
            States are defined about the barycenter of the two primaries, P1 and P2
            Default is None
        Returns
        -------
        dist_p1_p3 : float64/complex128
            distance between P3 and P1
        dist_p2_p3 : float64/complex128
            distance between P3 and P2
        """
        if state is None:
            state = self.ic
        
        dist_p1_p3 = ((state[0] + self.mu) ** 2 + state[1] ** 2 + state[2] ** 2) ** 0.5
        dist_p2_p3 = ((state[0] - 1 + self.mu) ** 2 + state[1] ** 2 + state[2] ** 2) ** 0.5
    
        return dist_p1_p3, dist_p2_p3
    
    
    def uii_partials_cr3bp(self, state=None):
        """Compute second-derivate of the pseudo-potenital of the CR3BP EOMs
    
        Parameters
        ----------
        state : numpy ndarray (6x1)
            state are defined about the barycenter of the two primaries, P1 and P2
            Default is None
    
        Returns
        -------
        Uxx, Uyy, Uzz, Uxy, Uxz, Uyz: float64/complex128
            Second-derivaitves of the pseudo-potenial of CR3BP
    
        """
        if state is None:
            state = self.ic
            
        dist_p1_p3, dist_p2_p3 = self.rel_dist_cr3bp(state)
        
        one_minus_mu = 1 - self.mu
        x_plus_mu = state[0] + self.mu
        x_minus_1_plus_mu = state[0] - 1 + self.mu
        
        one_minus_mu_dist_p1_p3_3 = one_minus_mu / dist_p1_p3**3
        mu_dist_p2_p3_3 = self.mu / dist_p2_p3**3
        
        # Second order partials of U for A(t) matrix (6x6)        
        Uxx = (
            1 - one_minus_mu_dist_p1_p3_3
            - mu_dist_p2_p3_3
            + 3 * one_minus_mu * x_plus_mu ** 2 / dist_p1_p3**5
            + 3 * self.mu * x_minus_1_plus_mu ** 2 / dist_p2_p3**5
        )
        Uyy = (
            1 - one_minus_mu_dist_p1_p3_3
            - mu_dist_p2_p3_3
            + 3 * one_minus_mu * state[1] ** 2 / dist_p1_p3**5
            + 3 * self.mu * state[1] ** 2 / dist_p2_p3**5
        )
        Uzz = (
            - one_minus_mu_dist_p1_p3_3
            - mu_dist_p2_p3_3
            + 3 * one_minus_mu * state[2] ** 2 / dist_p1_p3**5
            + 3 * self.mu * state[2] ** 2 / dist_p2_p3**5
        )
        Uxy = (
            3 * one_minus_mu * x_plus_mu * state[1] / dist_p1_p3**5
            + 3 * self.mu * x_minus_1_plus_mu * state[1] / dist_p2_p3**5
        )
        Uxz = (
            3 * one_minus_mu * x_plus_mu * state[2] / dist_p1_p3**5
            + 3 * self.mu * x_minus_1_plus_mu * state[2] / dist_p2_p3**5
        )
        Uyz = (
            3 * one_minus_mu * state[1] * state[2] / dist_p1_p3**5
            + 3 * self.mu * state[1] * state[2] / dist_p2_p3**5
        )
    
        return Uxx, Uyy, Uzz, Uxy, Uxz, Uyz
    
    
    def ui_partials_acc_cr3bp(self, state=None):
        """Compute first-derivateive of pseudo-potenital terms of the CR3BP EOMs and acceleration terms
        
        Parameters
        ----------
        state : numpy ndarray (6x1)
            state are defined about the barycenter of the two primaries, P1 and P2
            Default is None
    
        Returns
        -------
        Ux, Uy, Uz, ax, ay, az: float64/complex128
            First-derivaitves of the pseudo-potenial of CR3BP and the accelration components
        """
        if state is None:
            state = self.ic
        
        dist_p1_p3, dist_p2_p3 = self.rel_dist_cr3bp(state)
        
        one_minus_mu = 1 - self.mu # 1-mu
        x_plus_mu = state[0] + self.mu # x+mu
        x_minus_1_plus_mu = state[0] - 1 + self.mu # x-1+mu
        
        one_minus_mu_dist_p1_p3_3 = one_minus_mu / dist_p1_p3**3 # (1-mu)/d13^3
        mu_dist_p2_p3_3 = self.mu / dist_p2_p3**3 # mu/d23^3

        # Calculate Ux, Uy, Uz, ax, ay, az
        Ux = state[0] - one_minus_mu_dist_p1_p3_3 * x_plus_mu - mu_dist_p2_p3_3 * x_minus_1_plus_mu
        Uy = state[1] - one_minus_mu_dist_p1_p3_3 * state[1] - mu_dist_p2_p3_3 * state[1]
        Uz = -one_minus_mu_dist_p1_p3_3 * state[2] - mu_dist_p2_p3_3 * state[2]
        ax = 2 * state[4] + Ux
        ay = -2 * state[3] + Uy
        az = Uz
    
        return Ux, Uy, Uz, ax, ay, az    
    
    def JC(self, state=None):
        """Computes Jacobi Constant/Jacobi Integral [nd], CR3BP integral constant
        Can handle COMPLEX inputs

        Parameters
        ----------
            state : numpy ndarray (6x1)
            state are defined about the barycenter of the two primaries, P1 and P2
        
        Returns
        -------
            JC: Jacobi Constant, [nd]
        """
        
        if state is None:
            state = self.ic
        
        dist_p1_p3, dist_p2_p3 = self.rel_dist_cr3bp(state)

        JC = (state[0] ** 2
            + state[1] ** 2
            + 2 * (1 - self.mu) / dist_p1_p3
            + 2 * self.mu / dist_p2_p3
            - np.sqrt(state[3] ** 2 + state[4] ** 2 + state[5] ** 2) ** 2
        )  # Not used linalg.norm for |v|^2 as it is not compatiable with complex inputs
        
        return JC