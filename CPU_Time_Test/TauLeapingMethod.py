import random
import math
import numpy as np

def tau_leaping_method(x0, rates, stoich, times, options=None):
    """
    Simulate a stochastic system using the Tau Leaping Method.
    
    Parameters:
        x0 : list
            Initial state vector.
        rates : function
            Function that computes reaction propensities given the state and time.
        stoich : dict
            Dictionary with stoichiometry matrix under key "nu".
        times : list
            List of time points.
        options : dict, optional
            Additional options for simulation. Expects "dt" as the time step.
    
    Returns:
        Z : list of lists
            State of the system at each time step.
        t : list
            Time points corresponding to states.
    """
    X0 = x0[:]
    nu = stoich["nu"]
    t_final = times[-1]
    tau = options.get("dt", 1e-3) if options else 1e-3
    
    Nt = math.floor(t_final / tau)
    Z = [[0] * (Nt + 1) for _ in range(len(X0))]
    t = [0] * (Nt + 1)
    
    for i in range(len(X0)):
        Z[i][0] = X0[i]
    
    for i in range(Nt):
        # Compute propensities
        props = rates([Z[j][i] for j in range(len(X0))], t[i])
        props = [max(0, p) for p in props]  # Ensure no negative propensities
        
        # Generate Poisson variates
        # Y = [poisson(p * tau) if p > 0 else 0 for p in props]
        Y = [np.random.poisson(p * tau) if p > 0 else 0 for p in props]
        
        # Update copy numbers
        for j in range(len(X0)):
            Z[j][i + 1] = Z[j][i] + sum(nu[k][j] * Y[k] for k in range(len(nu)))
        
        # Update time
        t[i + 1] = t[i] + tau
    
    return Z, t


# def poisson(lmbda):
#     """Generate a Poisson random variable."""
#     L = math.exp(-lmbda)
#     k = 0
#     p = 1
#     while p > L:
#         k += 1
#         p *= random.uniform(0, 1)
#     return k - 1
