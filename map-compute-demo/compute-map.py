import numpy as np
from scipy import stats
import pandas as pd
import plotnine as p9
from plotnine import *

# Put all the input/output files here so that they are easy to find
# and change later.

input_csv = "snpsht_df.csv"
output_png_func = lambda a, b: f"out/scatter_{a}_{b}.png"
output_csv = "out/estimate_df.csv"

# There is a lot of data in the snapshot dataframe, but we only need
# the parameter values at the final time point. So, we can filter out
# the rest.

params_df = pd.read_csv(input_csv)
final_time = params_df["time"].max()
params_df = params_df[params_df["time"] == final_time]
param_names = ["lnV0", "beta", "delta", "pi", "phi", "rho"]
params_df = params_df[param_names]

# We can approximate the MAP by fitting a KDE and selecting the
# particular that maximises that density. This is crude but it
# definitely gives us a point within the support and assuming the
# particles have well explored the space it should be reasonable.

params_mat = params_df.to_numpy().T
kernel = stats.gaussian_kde(params_mat)
kernal_pdf_vals = kernel(params_mat)
idx = kernal_pdf_vals.argmax()

map_dict = {var: [params_mat[i, idx]] for i, var in enumerate(param_names)}
map_dict["method"] = "map"
map_df = pd.DataFrame(map_dict)


mean_dict = {var: [params_df[var].mean()] for var in param_names}
mean_dict["method"] = "mean"
mean_df = pd.DataFrame(mean_dict)

# Always a good idea to get a recording of the actual values to refer
# to later. But the plots are also a good way to see how the results
# differ between the element-wise mean the and MAP.

summary_df = pd.concat([map_df, mean_df])
summary_df.to_csv(output_csv, index=False)


for i in range(len(param_names)):
    for j in range(i+1, len(param_names)):
        scatter_p9 = (p9.ggplot(params_df, p9.aes(x=param_names[i], y=param_names[j]))
                      + p9.geom_point()
                      + p9.geom_point(mean_df, size=5, color='red')
                      + p9.geom_point(map_df, size=5, color='blue')
                      + p9.theme_bw())
        scatter_p9.save(output_png_func(param_names[i], param_names[j]),
                        height=4, width=6, dpi=300)
