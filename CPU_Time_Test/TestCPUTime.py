import pandas as pd
import random
import matplotlib.pyplot as plt
import jsf
import time

import pdb

from GillespieDirectMethod import gillespie_direct_method
from TauLeapingMethod import tau_leaping_method

# Define the initial conditions
S0 = 1*pow(10, 2)
I0 = 20
R0 = 0

# Define the Parameters
pbeta = 2/7
pgamma = 1/7
pwane = 1/365
pmu =    1/(85*365)
pkappa = 1/(85*365)

# Define the final time
t_max = 2*365

x0 = [S0, I0, R0]

# Define the stoichiometry
reactant_matrix =[[1, 1, 0],
                  [0, 1, 0],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [0, 0, 1]]

product_matrix = [[0,2,0],
                  [0,0,1],
                  [2,0,0],
                  [1,1,0],
                  [1,0,1],
                  [0,0,0],
                  [0,0,0],
                  [0,0,0],
                  [1,0,0]]

# Define the rates
rates = lambda x, t: [pbeta  * x[0] * x[1] / ( x[0] + x[1] + x[2] ),
                      pgamma * x[1],
                      pkappa * x[0],
                      pkappa * x[1],
                      pkappa * x[2],
                      pmu    * x[0],
                      pmu    * x[1],
                      pmu    * x[2],
                      pwane  * x[2]]

stoich = {
         "nu": [ [a - b for a, b in zip(r1, r2)]
                for r1, r2 in zip(product_matrix, reactant_matrix) ],
         "DoDisc": [0, 0, 0],
         "nuReactant": reactant_matrix,
         "nuProduct": product_matrix,
         }

switch_th = 3
my_opts = {
            "EnforceDo": [0, 0, 0],
            "dt": 0.005,
            "SwitchingThreshold": [pow(10,switch_th), pow(10,switch_th), pow(10,switch_th)]
           }


# track the CPU time for each method
shift_pop = 2
pop_it = 5
runs = 3

jsf_times = [[0 for _ in range(runs)] for _ in range(pop_it)]
gillespie_times = [[0 for _ in range(runs)] for _ in range(pop_it)]
tau_times = [[0 for _ in range(runs)] for _ in range(pop_it)]

jsf_avg = [0 for _ in range(pop_it)]
gillespie_avg = [0 for _ in range(pop_it)]
tau_avg = [0 for _ in range(pop_it)]

for pop_size in range(0, pop_it):
    S0 = 1*pow(10, pop_size+shift_pop)
    I0 = 2
    R0 = 0
    x0 = [S0, I0, R0]

    jsf_times = [0 for _ in range(runs)]
    gillespie_times = [0 for _ in range(runs)]
    tau_times = [0 for _ in range(runs)]

    for i in range(0, runs):  # Run multiple times to track CPU time over all runs
        start_jsf = time.process_time()
        sim_jsf = jsf.jsf(x0, rates, stoich, t_max, config=my_opts, method="operator-splitting")
        finish_jsf = time.process_time()
        jsf_times[i] = (finish_jsf - start_jsf)

        start_gillespie = time.process_time()
        sim_gil = gillespie_direct_method(x0, rates, stoich, [0, t_max], options=None)
        finish_gillespie = time.process_time()
        gillespie_times[i] = (finish_gillespie - start_gillespie)

        start_tau = time.process_time()
        sim_tau = tau_leaping_method(x0, rates, stoich, [0, t_max], options={"dt": 0.001})
        finish_tau = time.process_time()
        tau_times[i] = (finish_tau - start_tau)

    jsf_avg[pop_size] = sum(jsf_times) / len(jsf_times)
    gillespie_avg[pop_size] = sum(gillespie_times) / len(gillespie_times)
    tau_avg[pop_size] = sum(tau_times) / len(tau_times)

    print(f"Population Size: 10^{pop_size+shift_pop}")
    print("|  JSF  |  Gil  |  Tau  |")
    print("|-------|-------|-------|")
    print(f"| {jsf_avg[pop_size]:.3f} | {gillespie_avg[pop_size]:.3f} | {tau_avg[pop_size]:.3f} |", end="\n\n")




# plot the results with scatter plot

plt.scatter([pow(10, i+shift_pop) for i in range(pop_it)], jsf_avg, label="JSF")
plt.scatter([pow(10, i+shift_pop) for i in range(pop_it)], gillespie_avg, label="Gillespie")
plt.scatter([pow(10, i+shift_pop) for i in range(pop_it)], tau_avg, label="Tau")
plt.xlabel("Population Size")
plt.ylabel("Average CPU Time")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()



# print(f"Population Size: 10^{pop_size}")
# print(f"Average JSF CPU Time: {sum(jsf_times) / len(jsf_times)} seconds")
# print(f"Average Gillespie CPU Time: {sum(gillespie_times) / len(gillespie_times)} seconds")
# print(f"Average Tau CPU Time: {sum(tau_times) / len(tau_times)} seconds")
# print("-" * 50)


# # Plot the results using subplots
# fig, axs = plt.subplots(3, 1, figsize=(10, 8))

# axs[0].plot(sim_jsf[1], sim_jsf[0][0], label="S_JSf")
# axs[0].plot(sim_gil[1], sim_gil[0][0], label="S_Gil")
# axs[0].plot(sim_tau[1], sim_tau[0][0], label="S_Tau")
# axs[0].set_ylabel("Susceptible")
# axs[0].legend()

# axs[1].plot(sim_jsf[1], sim_jsf[0][1], label="I_JSf")
# axs[1].plot(sim_gil[1], sim_gil[0][1], label="I_Gil")
# axs[1].plot(sim_tau[1], sim_tau[0][1], label="I_Tau")
# axs[1].set_ylabel("Infected")
# axs[1].legend()

# axs[2].plot(sim_jsf[1], sim_jsf[0][2], label="R_JSf")
# axs[2].plot(sim_gil[1], sim_gil[0][2], label="R_Gil")
# axs[2].plot(sim_tau[1], sim_tau[0][2], label="R_Tau")
# axs[2].set_xlabel("Time")
# axs[2].set_ylabel("Recovered")
# axs[2].legend()

# plt.tight_layout()
# plt.show()
