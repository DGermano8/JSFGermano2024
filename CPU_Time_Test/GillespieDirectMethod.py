import random

def gillespie_direct_method(x0, rates, stoich, times, options=None):
    """
    Simulate a stochastic system using the Gillespie Direct Method.
    
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
            Additional options for simulation.
    
    Returns:
        X : list of lists
            State of the system at each time step.
        TauArr : list
            Time points corresponding to states.
    """
    X0 = x0[:]
    nu = stoich["nu"]
    
    t_final = times[-1]
    
    n_rates = len(nu)
    n_compartments = len(X0)
    
    # Initialize solution arrays
    X = [[0] * (2 * len(times)) for _ in range(n_compartments)]
    TauArr = [0] * (2 * len(times))
    for i in range(n_compartments):
        X[i][0] = X0[i]
    iters = 0
    
    while TauArr[iters] <= t_final:
        # Compute propensities
        props = rates([X[i][iters] for i in range(n_compartments)], TauArr[iters])
        total_propensity = sum(props)
        
        if total_propensity == 0:
            break
        
        # Sample exponential waiting time
        dt = random.expovariate(total_propensity)
        
        # Check if the simulation is finished
        if TauArr[iters] + dt > t_final:
            break
        
        iters += 1
        
        # Sample the next reaction event
        r = random.uniform(0, total_propensity)
        cumulative = 0
        j = -1
        for k in range(n_rates):
            cumulative += props[k]
            if r < cumulative:
                j = k
                break
        
        # Update copy numbers
        for i in range(n_compartments):
            X[i][iters] = X[i][iters - 1] + nu[j][i]
        TauArr[iters] = TauArr[iters - 1] + dt
        
        # Expand arrays if necessary
        if iters >= len(TauArr) - 1:
            for i in range(n_compartments):
                X[i].extend([0] * len(TauArr))
            TauArr.extend([0] * len(TauArr))
    
    # Trim excess entries
    for i in range(n_compartments):
        X[i] = X[i][:iters]
    TauArr = TauArr[:iters]
    
    return X, TauArr
