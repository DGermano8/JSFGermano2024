import numpy as np

def gillespie_direct_method(x0, rates, stoich, times, options=None):
    """
    Simulate a stochastic system using the Gillespie Direct Method.
    
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
            Additional options for simulation.
    
    Returns:
        X : ndarray
            Array of states over time.
        TauArr : ndarray
            Array of time points corresponding to states.
    """
    X0 = np.array(x0)
    nu = np.array(stoich["nu"])
    
    t_final = times[-1]
    
    n_rates, n_compartments = nu.shape
    
    # Initialize solution arrays
    X = np.zeros((n_compartments, 2 * len(times)))
    X[:, 0] = X0
    TauArr = np.zeros(2 * len(times))
    iters = 0
    
    while TauArr[iters] <= t_final:
        # Compute propensities
        props = rates(X[:, iters], TauArr[iters])
        
        # Sample exponential waiting time
        dt = np.random.exponential(1 / np.sum(props))
        
        # Check if the simulation is finished
        iters += 1
        if TauArr[iters - 1] + dt > t_final:
            break
        
        # Sample the next reaction event
        j = np.random.choice(n_rates, p=props / np.sum(props))
        
        # Update copy numbers
        X[:, iters] = X[:, iters - 1] + nu[j, :]
        TauArr[iters] = TauArr[iters - 1] + dt
        
        # Expand arrays if necessary
        if iters >= len(TauArr) - 1:
            X = np.hstack((X, np.zeros((n_compartments, len(TauArr)))))
            TauArr = np.hstack((TauArr, np.zeros(len(TauArr))))
    
    # Trim excess entries
    X = X[:, :iters]
    TauArr = TauArr[:iters]
    
    return X, TauArr
