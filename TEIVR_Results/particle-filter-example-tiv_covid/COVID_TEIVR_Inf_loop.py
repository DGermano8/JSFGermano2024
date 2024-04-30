
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
import pickle
import multiprocessing
import multiprocessing


# set random seed
random.seed(1)

patient_list = ['432192', '443108', '444332', '444391', '445602', '451152']

def process_patient(ii):
    try:
        # Your existing code here
        inst = list(pypfilt.load_instances("config/cli-refractory-tiv-jsf.toml"))[0]
        patient_id = '/' + patient_list[ii]

        out_dir = 'outputs4/' + patient_id + '/' +  inst.settings['components']['model'] + '_' + str(inst.settings['filter']['particles'])
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

        forecast_time = 10

        fit_result = pypfilt.forecast(inf_ctx, [forecast_time], filename=None)

        with open(f"{out_dir}/fit_result.pkl", "wb") as file_to_pickle:
            pickle.dump(fit_result, file_to_pickle)
    
    except Exception as e:
        print(e)
        print('Error')

if __name__ == '__main__':
    num_cores = 3  # Specify the number of cores to use

    pool = multiprocessing.Pool(processes=num_cores)
    pool.map(process_patient, range(0, len(patient_list)))
    pool.close()
    pool.join()


    