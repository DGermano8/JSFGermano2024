import numpy as np

def tau_leaping_method(x0, rates, stoich, times, options=None):
    """
    Simulate a stochastic system using the Tau Leaping Method.
    
    Parameters:
        x0 : array-like
            Initial state vector.
        rates : function
            Function that computes reaction propensities given the state and time.
        stoich : dict
            Dictionary with stoichiometry matrix under key "nu".
        times : array-like
            Array of time points.
        options : dict, optional
            Additional options for simulation. Expects "dt" as the time step.
    
    Returns:
        Z : ndarray
            Array of system states at each time step.
        t : ndarray
            Array of time points corresponding to states.
    """
    X0 = np.array(x0)
    nu = np.array(stoich["nu"])
    t_final = times[-1]
    tau = options.get("dt", 1e-3) if options else 1e-3
    
    Nt = int(np.floor(t_final / tau))
    Z = np.zeros((len(X0), Nt + 1))
    t = np.zeros(Nt + 1)
    
    Z[:, 0] = X0
    
    for i in range(Nt):
        # Compute propensities
        props = rates(Z[:, i], t[i])
        props = np.maximum(props, 0)  # Ensure no negative propensities
        
        # Generate Poisson variates
        Y = np.random.poisson(props * tau)
        
        # Update copy numbers
        Z[:, i + 1] = Z[:, i] + nu.T @ Y
        
        # Update time
        t[i + 1] = t[i] + tau
    
    return Z, t
