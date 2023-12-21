# run-inference-from-simulation.py
#
# program: Inference on simulated data
#
# programmer: Alexander E. Zarebski
# date: 2023-08-29
#
# description: A simple example of parameter estimation using
# simulated data based on parameters from Baccam et al (2006).
#
#
from typing import List, Dict, Any, Optional
import numpy as np
import scipy.stats              # type: ignore
import pypfilt                  # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.image as mpimg # type: ignore
import matplotlib as matplotlib # type: ignore
matplotlib.use('QtAgg')
import pandas as pd
import plotnine as p9
from plotnine import ggplot, geom_rect, aes, geom_ribbon, geom_point, scale_y_log10, scale_x_continuous, labs, theme_bw, geom_vline
import pdb
from src.inf import *

# import the src.inf package again using the reload method from the imports module
from importlib import reload
import src.inf
reload(src.inf)


def const_params_from_prior(prior : Dict) -> Dict[str, float]:
    """
    Extract the constant parameters used in the simulation from the
    prior.
    """
    param_names = ["lnV0", "beta", "p", "c", "gamma"]
    is_const = lambda n: prior[n]["name"] == "constant"
    return {p : prior[p]["args"]["value"]
            for p in param_names if is_const(p) }


def main():
    out_dir = "out"
    in_toml = "tiv-inference-from-simulation.toml"
    #
    # Define parameter and variable names
    #
    param_names = ['lnV0', 'beta', 'p', 'c', 'gamma']
    inst_dict = {
        x.scenario_id: x
        for x in pypfilt.load_instances(in_toml)
       }
    obs_ssv = inst_dict['simulation'].settings['sim_output_file']
    #
    # Simulate some observations
    #
    sim_params = const_params_from_prior(inst_dict['simulation'].settings["prior"])
    sim_result = pypfilt.simulate_from_model(inst_dict['simulation'])
    obs_df = pd.DataFrame(sim_result['V'])
    obs_df.to_csv(obs_ssv, sep = ' ', index = False)
    #
    # Run the particle filter over simulated data
    #
    inst_dict['inference'].settings['observations']['V']['file'] = "out/simulated-data-for-inference.ssv"
    inf_ctx = inst_dict['inference'].build_context()
    inf_results = tiv_run_inference(inf_ctx)
    #
    # Plot the prior-posterior distributions for parameters
    #
    end_time = inf_results['end_time']
    mrgs = inf_results['marginals']
    pst_param_df = inf_results['posterior_param_df']
    for param in param_names:
        dd = pst_param_df[pst_param_df['name'] == param]
        plt_df = plottable_model_cis(dd)
        param_p9 = param_plt_p9(plt_df, sim_params[param], mrgs[param],
                                param)
        param_p9.save(f"{out_dir}/demo-param-{param}-histogram.png",
                      height = 5.8, width = 8.3)
    #
    # Plot the state trajectory
    #
    pst_state_df = inf_results['posterior_state_df']
    plt_df = pst_state_df[pst_state_df['name'] == 'V']
    plt_df_obs = obs_df.copy()
    plt_df_obs['y'] = 10**plt_df_obs['value']
    state_p9 = state_plt_p9(plt_df, plt_df_obs)
    state_p9.save(f"{out_dir}/demo-state-trajectory.png",
            height = 4.1, width = 5.8)


if __name__ == "__main__":
    main()
