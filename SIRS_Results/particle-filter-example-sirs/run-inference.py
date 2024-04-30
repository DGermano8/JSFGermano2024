
import argparse
from typing import List, Dict, Any, Optional
import numpy as np
import scipy.stats              # type: ignore
import pypfilt                  # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.image as mpimg # type: ignore
import matplotlib as matplotlib # type: ignore
# matplotlib.use('QtAgg')
import pandas as pd
import plotnine as p9
from plotnine import ggplot, geom_rect, aes, geom_ribbon, geom_point, scale_y_log10, scale_x_continuous, labs, theme_bw, geom_vline
import pdb
from inf import plottable_model_cis, param_plt_p9, state_plt_p9, sirs_run_inference
import os
import random

random.seed(2)

# inst = list(pypfilt.load_instances("sirs-inference.toml"))[0]
inst = list(pypfilt.load_instances("sirs-inference.toml"))[0]

outDir = 'outputs/' 

inf_ctx = inst.build_context()
# inf_results = sirs_run_inference(inf_ctx)

state_names = ['S', 'I', 'R']
param_names = ["betaCoef", "gammaCoef", "omegaCoef"]
prior = inf_ctx.settings['prior']

# out_dir = outDir +  inst.settings['components']['model'] + '_' + str(inst.settings['filter']['particles']) + '_04'

out_dir = outDir +  inst.settings['components']['model'] + '_particle=' + str(inst.settings['filter']['particles'])  + '_intMag=' + str(prior['intervention_magnitude']['args']['value'])

# make output directory if it doesn't exist
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

inst.settings['observations']['I']['file'] = outDir + '/sim_df.ssv'

has_prior = lambda n: prior[n]["name"] != "constant"
mrgs = {p : prior[p]
            for p in param_names if has_prior(p) }
end_time = inf_ctx.settings['time']['until']
# fit_result = pypfilt.fit(inf_ctx, filename=None)
#                                       list of forecast times (time to start foecast, until time "until" in .toml file)
forecast_time = 100

extinction_time = end_time

fit_result = pypfilt.forecast(inf_ctx, [forecast_time], filename=None)

# pdb.set_trace()
pst_df = pd.DataFrame(fit_result.estimation.tables['model_cints'])
pst_state_df = pd.DataFrame(fit_result.forecasts[forecast_time].tables['backcasts'])
pst_state_df = pst_state_df[['unit', 'time', 'prob','ymin', 'ymax']]
pst_param_df = pst_df[pst_df['time'] == forecast_time]
pst_param_df = pst_param_df[pst_param_df['name'].isin(param_names)]
pst_param_df = pst_param_df[['prob','ymin', 'ymax', 'name']]

prd_state_df = pd.DataFrame(fit_result.forecasts[forecast_time].tables['forecasts'])
prd_state_df = prd_state_df[['unit', 'time', 'prob','ymin', 'ymax']]

betaEst = pst_df[pst_df['name'] == 'betaCoef']
betaEst.to_csv(f"{out_dir}/betaEst.csv", index=False)
gammaEst = pst_df[pst_df['name'] == 'gammaCoef']
gammaEst.to_csv(f"{out_dir}/gammaEst.csv", index=False)
omegaEst = pst_df[pst_df['name'] == 'omegaCoef']
omegaEst.to_csv(f"{out_dir}/omegaEst.csv", index=False)

ext_prdf = pd.DataFrame(fit_result.forecasts[forecast_time].tables['extinction_prob'])
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
# param_names = ["betaCoef", "gammaCoef", "omegaCoef", "kappaCoef", "muCoef"]
param_names = ["betaCoef", "gammaCoef", "omegaCoef",]

end_time = inf_results['end_time']

mrgs = inf_results['marginals']
pst_param_df = inf_results['posterior_param_df']
pst_param_df.to_csv(f"{out_dir}/pst_param_df.csv", index=False)

for param in param_names:
    # if cli_args['verbose']:
    #     print(f"\tGenerating plot for parameter: {param}")
    dd = pst_param_df[pst_param_df['name'] == param]
    plt_df = plottable_model_cis(dd)
    plt_df.to_csv(f"{out_dir}/plt_df_{param}.csv", index=False)
    param_p9 = param_plt_p9(plt_df, None, mrgs[param],
                            param)
    param_p9.save(f"{out_dir}/demo-param-{param}-histogram.png",
                  height = 5.8, width = 8.3)

cli_args = {'obs_ssv': outDir + '/sim_df.ssv'}
pst_state_df = inf_results['posterior_state_df']
state_df = pd.concat([pst_state_df, prd_state_df])

plt_df = state_df[state_df['unit'] == 'I']
plt_df['ymin'] = plt_df['ymin']
plt_df['ymax'] = plt_df['ymax']
plt_df_obs = pd.read_csv(cli_args['obs_ssv'], sep = ' ')
plt_df_obs['y'] = plt_df_obs['value']
# write the pandas dataframe to csv
plt_df.to_csv(f"{out_dir}/plt_df.csv", index=False)
state_p9 = state_plt_p9(plt_df, plt_df_obs)
state_p9.save(f"{out_dir}/demo-state-trajectory.png",
                      height = 4.1, width = 5.8)


# this is where we get a snapshot of the parameter distribution at the final time
# and then write each of the estimated parameters to a separate csv file
snpsht_df = pd.DataFrame(fit_result.estimation.tables['snapshot'])

# input_csv = 'outputs4/' + patient_list[ii] + '/src.tiv.RefractoryCellModel_JSF_6000/snpsht_df.csv'
# out_dir = 'out/' + patient_list[ii] + '/src.tiv.RefractoryCellModel_JSF_6000'

params_df = snpsht_df
final_time = forecast_time
params_df = params_df[params_df["time"] == final_time]
params_df = params_df[param_names]

for i in range(len(param_names)):
    (params_df[param_names[i]]).to_csv(out_dir + '/' + param_names[i]+'_distribution.csv', index=False)  
