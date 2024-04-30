
import numpy as np
import scipy.stats
import pypfilt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as matplotlib
# matplotlib.use('QtAgg') 
import pandas as pd
import plotnine as p9
from plotnine import *
import pdb
import random

# this seed will give us a trajectory that goes extinct, making it hard to to ODE inference
random.seed(6)

def run_simulation(instance):
    """
    Simulates a model using pypfilt, and returns a DataFrame of results.

    :param instance: A simulation instance containing settings and parameters.
    :return: A pandas DataFrame containing the simulation results.
    """
    num_reps = instance.settings['num_replicates']
    my_obs_tables = pypfilt.simulate_from_model(instance, particles = num_reps)
    # pdb.set_trace()
    sim_df = pd.DataFrame(my_obs_tables['I'])
    sim_df['particle'] = np.tile(
        np.arange(num_reps),
        sim_df.shape[0] // num_reps
    )

    outDit = 'outputs'

    # write this data frame to a csv file
    sim_df.to_csv(outDit + '/sim_df.csv', index = False)

    # write this data to a ssv file
    sim_df.to_csv(outDit + '/sim_df.ssv', sep = ' ', index = False)

    return sim_df


def simulation_plot(sim_df, name):
    """
    Generates a plot using ggplot to visualize the simulation data.

    :param sim_df: A pandas DataFrame containing the simulation results.
    :param name: The title of the plot.
    :return: A ggplot object representing the plot.
    """
    return (ggplot()
            + geom_line(data = sim_df,
                        mapping = aes(x = "time",
                                      y = "value",
                                      group = "particle"),
                        alpha = 0.3)
            # + scale_y_sqrt(limits = (0, 600))
            + labs(title = name,
                   y = "Population size",
                   x = "Time")
            + theme_bw())


def main():
    """
    Main function that loads instances from a scenario file, runs
    simulations, and saves the plots as a PDF.

    :return: None
    """
    scenario_file = 'sirs-simulation.toml'
    instance_dict = {x.scenario_id: x for x
                     in pypfilt.load_instances(scenario_file)}
    plots = []
    for key, inst in instance_dict.items():
        sim_name = inst.settings['simulation_name']
        plot_path = inst.settings['plot_path']
        sim_df = run_simulation(inst)
        p = simulation_plot(sim_df, sim_name)
        p.save(plot_path)
        plots.append(p)

    p9.save_as_pdf_pages(plots, filename = "outputs/demo-simulations.pdf")


if __name__ == '__main__':
    main()
