from importlib import reload
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
import src.tiv
reload(src.tiv)


def posterior_dataframes(results):
    """
    Posterior dataframes from the results object.

    :param results: Results object from the inference.
    :return: Dictionary of dataframes.
    """
    post_cints_df = pd.DataFrame(results.estimation.tables['model_cints'])

    def subset_and_order(n):
        tmp = post_cints_df[post_cints_df['name'] == n].copy()
        tmp['prob'] = pd.Categorical(
            tmp['prob'],
            categories=tmp['prob'].unique(),
            ordered=True
        )
        return tmp
    return {'V0_df': subset_and_order('V0'),
            'beta_df': subset_and_order('beta'),
            'p_df': subset_and_order('p'),
            'c_df': subset_and_order('c'),
            'gamma_df': subset_and_order('gamma'),
            'state_T_df': subset_and_order('T'),
            'state_I_df': subset_and_order('I'),
            'state_V_df': subset_and_order('V'),
            'observations': pd.DataFrame(results.obs)}


#
#
#   ***********************
#   *                     *
#   *  Run the inference  *
#   *                     *
#   * ******************* *
#
#
scenario_file = 'tiv-inference-from-simulation.toml'
assert os.path.exists(scenario_file)
instances = list(pypfilt.load_instances(scenario_file))

inf_inst = instances[0]
obs_ssv = inf_inst.settings['observations']['V']['file']
assert os.path.exists(obs_ssv)
time_scale = inf_inst.time_scale()

# # Carry out the inference using the inference instances and
# # generate plots to summarise the results.
# for inf_inst in instances[1:]:
output_id = inf_inst.settings['patient_number']
fcst_time = inf_inst.settings['forecast_time']
ctx = inf_inst.build_context()
results = pypfilt.forecast(ctx, [fcst_time], filename=None)
posterior = posterior_dataframes(results)
#
#
#   ************************
#   *                      *
#   *  Plot the model fit  *
#   *                      *
#   ************************
#
#
v2obs = lambda v: np.log10(v).round(1).clip(lower=0.5)

plt_est_df = posterior['state_V_df']
plt_est_df = plt_est_df[plt_est_df['prob'] == 0]
plt_est_df = plt_est_df.assign(tcid = v2obs(plt_est_df['ymin']))
plt_est_df = plt_est_df[['time', 'tcid']]

plt_ci_df = posterior['state_V_df']
plt_ci_df = plt_ci_df[plt_ci_df['prob'] == 95]
plt_ci_df = plt_ci_df.assign(tcid_lower = v2obs(plt_ci_df['ymin']),
                             tcid_upper = v2obs(plt_ci_df['ymax']))
plt_ci_df = plt_ci_df[['time', 'tcid_lower', 'tcid_upper']]

data_df = posterior['observations']

(ggplot()
    + geom_ribbon(
        data = plt_ci_df,
        mapping = aes(x = "time",
                        ymin = "tcid_lower",
                        ymax = "tcid_upper"),
        alpha = 0.2)
 + geom_line(
        data = plt_est_df,
        mapping = aes(x = "time",
                        y = "tcid"))
 + geom_point(
        data = data_df,
        mapping = aes(x = "time",
                        y = "value"),
     color = "red")
 + theme_bw()
)
