import numpy as np
import scipy.stats
import pypfilt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as matplotlib
matplotlib.use('QtAgg')
import pandas as pd
import plotnine as p9
from plotnine import *
import pdb
import tiv
from importlib import reload

reload(tiv)

def run_simulation(instance):
    """
    Runs a simulation from a given instance.

    Parameters
    ----------
    instance : pypfilt.Instance
        The instance to run the simulation from.

    Returns
    -------
    sim_df : pandas.DataFrame
        A dataframe containing the simulation results.
    """
    num_reps = instance.settings['num_replicates']
    my_obs_tables = pypfilt.simulate_from_model(
        instance, particles = num_reps
    )
    num_summaries = round(my_obs_tables['V'].size / num_reps)

    sim_df = pd.DataFrame({
        'time': my_obs_tables['V']['time'],
        'V': my_obs_tables['V']['value'],
        'T': my_obs_tables['T']['value'],
        'replicate': np.tile(np.arange(num_reps), num_summaries)
    })
    return sim_df


def read_data(filename, patient_num):
    """
    Reads the data from the data file.

    Parameters
    ----------
    filename : str
        The name of the data file.
    patient_num : int
        The patient number.

    Returns
    -------
    tcid_df : pandas.DataFrame
        A dataframe containing the data.
    """
    tcid_df = pd.read_csv("data/tcid.csv")
    tcid_df = tcid_df[tcid_df["patient"] == patient_num]
    tcid_df["log10_tcid"] = (tcid_df["log10_tcid"]
                             .astype(str)
                             .apply(lambda x:
                                    float(x[2:])
                                    if x.startswith("<=")
                                    else float(x)))
    tcid_df["is_truncated"] = tcid_df["log10_tcid"] <= 0.5
    return tcid_df


def write_simulation_data(sim_df, patient_num):
    """
    Writes the simulation results to a file.

    Parameters
    ----------
    sim_df : pandas.DataFrame
        A dataframe containing the simulation results.
    patient_num : int
        The patient number.
    """
    for replicate in sim_df["replicate"].unique():
        # zero-pad the replicate number to 3 digits
        output_ssv = f"out/simulated-data-patient-{patient_num}-replicte-{replicate:03d}.ssv"
        tmp = sim_df[sim_df["replicate"] == replicate]
        tmp = tmp.assign(value = tmp['V'].round(1).clip(lower=0.5))
        tmp = tmp[['time', 'value']]
        tmp = tmp[tmp['time'] > 0.0]
        tmp.to_csv(output_ssv,
                   sep = " ",
                   index = False,
                   header = True)


def write_simulation_plot(sim_df, tcid_df, patient_num):
    """
    Plots the simulation results.

    Parameters
    ----------
    sim_df : pandas.DataFrame
        A dataframe containing the simulation results.
    tcid_df : pandas.DataFrame
        A dataframe containing the raw patient data.
    patient_num : int
        The patient number.

    Returns
    -------
    plot : plotnine.ggplot
        The plot.
    """
    plot_df = pd.melt(sim_df,
                      id_vars = ['time', 'replicate'],
                      value_vars = ['V', 'T'])
    plot_title = f"Viral load (patient {tcid_df['patient'].unique()[0]})"
    p9 = (ggplot()
          + geom_line(
              data = plot_df[plot_df["variable"] == "T"],
              mapping = aes(x = "time",
                            y = "value")
          )
          + geom_line(
              data = plot_df[plot_df["variable"] == "V"],
              mapping = aes(x = "time",
                            y = "value",
                            group = "replicate"),
              alpha = 0.2)
          + geom_point(
              data = tcid_df,
              mapping = aes(x = "day",
                            y = "log10_tcid",
                            shape = "is_truncated"),
              color = "red"
          )
          + labs(title = plot_title,
                 y = "",
                 x = "Days post infection")
          + scale_y_continuous(limits = [-4, 10],
                               breaks = np.linspace(-4, 8, 7))
          + theme_bw()
          + theme(legend_position = "none"))

    p9.save(
        f"out/baccam-fit-{patient_num}.png",
        height = 4.1, width = 5.8 # A6
    )


def main():
    """
    """
    insts = list(pypfilt.load_instances('tiv-simulation.toml'))

    for inst in insts:
        patient_num = inst.settings['patient_number']
        sim_df = run_simulation(inst)
        write_simulation_data(sim_df, patient_num)
        data_df = read_data("data/tcid.csv", patient_num)
        write_simulation_plot(sim_df, data_df, patient_num)
    return None


if __name__ == '__main__':
    main()
