
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
import random

# set random seed
random.seed(1)

patient_list = ['432192', '443108', '444332', '444391', '445602', '451152']

for ii in range(0, len(patient_list)):
    try:
        inst = list(pypfilt.load_instances("config/cli-refractory-tiv-jsf.toml"))[0]
        patient_id = '/' + patient_list[ii]

        out_dir = 'outputs/' + patient_id + '/' +  inst.settings['components']['model'] + '_' + str(inst.settings['filter']['particles'])
        cli_args = {'obs_ssv': 'data/' + patient_id + '.ssv'}

        # make output directory if it doesn't exist
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        inst.settings['observations']['V']['file'] = 'data/' + patient_id + '.ssv'

        inf_ctx = inst.build_context()
        # inf_results = tiv_run_inference(inf_ctx)
        end_time = inf_ctx.settings['time']['until']
        prior = inf_ctx.settings['prior']

        param_names = ['lnV0','beta', 'phi', 'rho', 'delta', 'pi']
        has_prior = lambda n: prior[n]["name"] != "constant"
        mrgs = {p : prior[p]
                    for p in param_names if has_prior(p) }

        forecast_time = 14

        fit_result = pypfilt.forecast(inf_ctx, [forecast_time], filename=None)
        pst_df = pd.DataFrame(fit_result.estimation.tables['model_cints'])
        pst_state_df = pd.DataFrame(fit_result.forecasts[forecast_time].tables['backcasts'])
        pst_state_df = pst_state_df[['unit', 'time', 'prob','ymin', 'ymax']]
        pst_param_df = pst_df[pst_df['time'] == forecast_time]
        pst_param_df = pst_param_df[pst_param_df['name'].isin(param_names)]
        pst_param_df = pst_param_df[['prob','ymin', 'ymax', 'name']]

        snpsht_df = pd.DataFrame(fit_result.estimation.tables['snapshot'])
        snpsht_df.to_csv(f"{out_dir}/snpsht_df.csv", index=False)

        betaEst = pst_df[pst_df['name'] == 'beta']
        betaEst.to_csv(f"{out_dir}/betaEst.csv", index=False)
        phiEst = pst_df[pst_df['name'] == 'phi']
        phiEst.to_csv(f"{out_dir}/phiEst.csv", index=False)
        rhoEst = pst_df[pst_df['name'] == 'rho']
        rhoEst.to_csv(f"{out_dir}/rhoEst.csv", index=False)
        deltaEst = pst_df[pst_df['name'] == 'delta']
        deltaEst.to_csv(f"{out_dir}/deltaEst.csv", index=False)
        piEst = pst_df[pst_df['name'] == 'pi']
        piEst.to_csv(f"{out_dir}/piEst.csv", index=False)
        lnV0Est = pst_df[pst_df['name'] == 'lnV0']
        lnV0Est.to_csv(f"{out_dir}/lnV0Est.csv", index=False)

        covmat = pd.DataFrame(fit_result.forecasts[forecast_time].tables['covariencetable'])
        covmat.to_csv(f"{out_dir}/covmat.csv", index=False)

        prd_state_df = pd.DataFrame(fit_result.forecasts[forecast_time].tables['forecasts'])
        prd_state_df = prd_state_df[['unit', 'time', 'prob','ymin', 'ymax']]

        pst_ext_prdf = pd.DataFrame(fit_result.estimation.tables['extinction_prob'])
        prd_ext_prdf = pd.DataFrame(fit_result.forecasts[forecast_time].tables['extinction_prob'])
        ext_prdf = pd.concat([pst_ext_prdf, prd_ext_prdf])
        ext_prdf.to_csv(f"{out_dir}/ext_prdf.csv", index=False)

        pst_point_df = pd.DataFrame(fit_result.estimation.tables['point_ests'])
        pst_point_df = pst_point_df[pst_point_df['time'] == forecast_time]
        pst_point_df = pst_point_df[['prob', 'ymin', 'ymax', 'name']]

        inf_results = {'posterior_state_df': pst_state_df,
                    'posterior_param_df': pst_param_df,
                    'posterior_point_ests': pst_point_df,
                    'end_time': forecast_time,
                    'marginals': mrgs}

            # out_dir = 'ODE_Model'
            # out_dir = 'JSF_Model'

        end_time = inf_results['end_time']
        mrgs = inf_results['marginals']
        pst_param_df = inf_results['posterior_param_df']
        for param in param_names:
                # if cli_args['verbose']:
                #     print(f"\tGenerating plot for parameter: {param}")
            dd = pst_param_df[pst_param_df['name'] == param]
            plt_df = plottable_model_cis(dd)
            param_p9 = param_plt_p9(plt_df, None, mrgs[param],
                                    param)
            param_p9.save(f"{out_dir}/demo-param-{param}-histogram.png",
                        height = 5.8, width = 8.3)

        # cli_args = {'obs_ssv': 'data/patient-1-censored.ssv'}
        pst_state_df = inf_results['posterior_state_df']

        state_df = pd.concat([pst_state_df, prd_state_df])

        plt_df = state_df[state_df['unit'] == 'V']
        plt_df['ymin'] = np.power(10, plt_df['ymin'])
        plt_df['ymax'] = np.power(10, plt_df['ymax'])
        plt_df_obs = pd.read_csv(cli_args['obs_ssv'], sep = ' ')
        plt_df_obs['y'] = 10**plt_df_obs['value']
        # write the pandas dataframe to csv
        plt_df.to_csv(f"{out_dir}/plt_df.csv", index=False)
        state_p9 = state_plt_p9(plt_df, plt_df_obs)
        state_p9.save(f"{out_dir}/demo-state-trajectory.png",
                            height = 4.1, width = 5.8)
    
    except:
        print('Error')