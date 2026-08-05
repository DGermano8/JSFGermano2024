import gillespy2
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt
import jsf as jsf

from GillespieDirectMethod import gillespie_direct_method



class BirthDeathModel_GP2(gillespy2.Model):
    def __init__(self, x0=10, growth_rate=1, death_rate=0.5, t_max=5, switch_th=100):
        super().__init__(name="BirthDeathModel")
        
        # Parameters
        self.add_parameter(gillespy2.Parameter(name="k_birth", expression=growth_rate ))
        self.add_parameter(gillespy2.Parameter(name="k_death", expression=death_rate))
        
        # Species
        # self.add_species(gillespy2.Species(name="X", initial_value=x0, switch_min=switch_th, switch_tol=0.1))
        self.add_species(gillespy2.Species(name="X", initial_value=x0, constant=False, boundary_condition=False, mode='dynamic', allow_negative_populations=False,  switch_min=switch_th ))
        
        # Reactions
        self.add_reaction(gillespy2.Reaction(name="birth", reactants={'X': 1}, products={'X': 2}, rate=self.listOfParameters['k_birth']))
        self.add_reaction(gillespy2.Reaction(name="death", reactants={'X': 1}, products={}, rate=self.listOfParameters['k_death']))
        
        self.timespan(np.linspace(0, t_max, 1001))

def jsf_process(x0=10, growth_rate=1, death_rate=0.5, t_max=5, switch_th=100):
    # Implement the JSF process
    # Define the Parameters
    pgrowth = growth_rate
    pdeath = death_rate

    # Define the stoichiometry
    reactant_matrix =[[1],
                    [1]]

    product_matrix = [[2],
                    [0]]

    # Define the rates
    rates = lambda x, t: [pgrowth  * x[0],
                        pdeath * x[0]]

    stoich = {
            "nu": [ [a - b for a, b in zip(r1, r2)]
                    for r1, r2 in zip(product_matrix, reactant_matrix) ],
            "DoDisc": [0],
            "nuReactant": reactant_matrix,
            "nuProduct": product_matrix,
            }

    my_opts = {
                "EnforceDo": [0],
                "dt": 0.01,
                "SwitchingThreshold": [switch_th]
            }
    x0_v = [x0]

    sim_jsf = jsf.jsf(x0_v, rates, stoich, t_max, config=my_opts, method="operator-splitting")
    return sim_jsf
    

x0=1000
switch_th=200
growth_rate=0.5
death_rate =1.0
t_max=10
numb_of_traj=2

# initialise summary statistics of size numb_of_traj
summary_time_JSF = np.zeros(numb_of_traj)
summary_cpu_JSF = np.zeros(numb_of_traj)

for i in range(numb_of_traj):
    start_time = time.perf_counter()
    jsf_sim = jsf_process(x0,growth_rate, death_rate, t_max, switch_th)
    end_time = time.perf_counter()

    summary_cpu_JSF[i] = end_time - start_time
    

    time_reached_100_JSF = None
    for t, x in zip(jsf_sim[1], jsf_sim[0][0]):
        if x == 100:
            time_reached_100_JSF = t
            break
    summary_time_JSF[i] = time_reached_100_JSF
    

summary_time_GP2 = np.zeros(numb_of_traj)
summary_cpu_GP2 = 0

model = BirthDeathModel_GP2(x0, growth_rate, death_rate, t_max, switch_th)

start_time = time.perf_counter()
results_tau_hyb = model.run(algorithm="Tau-Hybrid", number_of_trajectories=numb_of_traj )

end_time = time.perf_counter()
summary_cpu_GP2 = (end_time - start_time)/numb_of_traj

for i in range(numb_of_traj):
    time_reached_100_GP2 = None
    x_p = 0
    for t, x in zip(results_tau_hyb[i]['time'], results_tau_hyb[i]['X']):
        if x <= 100 and x_p >= 100:
            time_reached_100_GP2 = t
            break
        x_p = x
    summary_time_GP2[i] = time_reached_100_GP2
  

# results_ode= model.run(algorithm="ODE")


# creat a subplot

fig, axs = plt.subplots(1, 2)
fig.suptitle('Comparison of JSF and GP2')
# plot the point estimate of the mean of summary_cpu_GP2 and summary_cpu_JSF
axs[0].bar(['GP2', 'JSF'], [summary_cpu_GP2, np.mean(summary_cpu_JSF)])
axs[0].set_xlabel('Method')
axs[0].set_ylabel('CPU Time')
axs[0].set_title('CPU Time')

# plot the histogram of summary_time_GP2 and summary_time_JSF
axs[1].hist(summary_time_GP2, bins=20, alpha=0.5, label='GP2')
axs[1].hist(summary_time_JSF, bins=20, alpha=0.5, label='JSF')
axs[1].legend(loc='upper right')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Frequency')
axs[1].set_title('Time to reach 100')

plt.show()
