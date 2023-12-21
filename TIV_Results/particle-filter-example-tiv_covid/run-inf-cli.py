# run-inf-cli.py
#
# program: Command line interface to inference
#
# programmer: Alexander E. Zarebski
# date: 2023-09-11
#
# description: A simple example of parameter estimation for the TIV
# model using data from a CSV file passed in as a command line
# argument along with a TOML to specify the inference context.
#
# example:
#
# $ python run-inf-cli.py --out out --input_toml cli-tiv-demo.toml --obs_ssv out/simulated-data-for-inference.ssv
#
#
import argparse
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
from src.inf import plottable_model_cis, param_plt_p9, state_plt_p9, tiv_run_inference
import os

def read_cli_args() -> Dict[str, Any]:
    """
    Read the command line arguments and return them as a dictionary.
    """
    parser = argparse.ArgumentParser(description="ClI for using pypfilt with the TIV model.")
    parser.add_argument('--out', required=True, help='Output directory')
    parser.add_argument('--input_toml', required=True, help='Input TOML file')
    parser.add_argument('--obs_ssv', required=True, help='Input SSV file with data')
    parser.add_argument('--param_plots', action='store_true', help='Generate parameter plots')
    parser.add_argument('--state_plot', action='store_true', help='Generate state plot')
    parser.add_argument('--record_summary', action='store_true', help='Summarise the estimates')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()
    return vars(args)


def main():
    cli_args = read_cli_args()
    out_dir = cli_args["out"]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        if cli_args['verbose']:
            print(f"Created output directory: {out_dir}")
            
    param_names = ['lnV0', 'beta', 'p', 'c', 'gamma']
    if cli_args['verbose']:
        print(f"Reading instances from TOML file: {cli_args['input_toml']}")
    inst = list(pypfilt.load_instances(cli_args["input_toml"]))[0]
    if cli_args['verbose']:
        print(f"Reading observations from SSV file: {cli_args['obs_ssv']}")
    inst.settings['observations']['V']['file'] = cli_args['obs_ssv']
    inf_ctx = inst.build_context()
    if cli_args['verbose']:
        print("Running inference")
    inf_results = tiv_run_inference(inf_ctx)

    if cli_args['param_plots']:
        if cli_args['verbose']:
            print("Generating parameter plots")
        end_time = inf_results['end_time']
        mrgs = inf_results['marginals']
        pst_param_df = inf_results['posterior_param_df']
        for param in param_names:
            if cli_args['verbose']:
                print(f"\tGenerating plot for parameter: {param}")
            dd = pst_param_df[pst_param_df['name'] == param]
            plt_df = plottable_model_cis(dd)
            param_p9 = param_plt_p9(plt_df, None, mrgs[param],
                                    param)
            param_p9.save(f"{out_dir}/demo-param-{param}-histogram.png",
                          height = 5.8, width = 8.3)

    if cli_args['state_plot']:
        if cli_args['verbose']:
            print("Generating state plot")
        pst_state_df = inf_results['posterior_state_df']
        plt_df = pst_state_df[pst_state_df['unit'] == 'V']
        plt_df['ymin'] = np.power(10, plt_df['ymin'])
        plt_df['ymax'] = np.power(10, plt_df['ymax'])
        plt_df_obs = pd.read_csv(cli_args['obs_ssv'], sep = ' ')
        plt_df_obs['y'] = 10**plt_df_obs['value']
        state_p9 = state_plt_p9(plt_df, plt_df_obs)
        state_p9.save(f"{out_dir}/demo-state-trajectory.png",
                      height = 4.1, width = 5.8)

    if cli_args['record_summary']:
        if cli_args['verbose']:
            print("Summarising the parameter estimates")
        pt_df = inf_results['posterior_point_ests']
        pt_df.to_csv(f"{out_dir}/demo-point-estimates.csv")

if __name__ == "__main__":
    main()
