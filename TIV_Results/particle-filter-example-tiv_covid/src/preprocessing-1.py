import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as matplotlib
matplotlib.use('QtAgg')
import pandas as pd
import plotnine as p9
from plotnine import *

# --------------------------------------------------------------------
# Read the data into a sensible dataframe.
# --------------------------------------------------------------------

tcid_df = pd.read_csv("data/tcid.csv")
tcid_df["log10_tcid"] = (tcid_df["log10_tcid"]
                         .astype(str)
                         .apply(lambda x: float(x[2:]) if x.startswith("<=") else float(x)))
tcid_df["is_truncated"] = tcid_df["log10_tcid"] <= 0.5

# --------------------------------------------------------------------
# Generate a nice plot to check it is in a good format.
# --------------------------------------------------------------------

data_p9 = (ggplot()
           + geom_point(
               data = tcid_df,
               mapping = aes(x = "day",
                             y = "log10_tcid",
                             colour = "is_truncated")
           )
           + geom_hline(yintercept = 0.5, linetype = "dashed")
           + facet_wrap("patient")
           + scale_y_continuous(
               limits = [-4, 8],
               breaks = [2 * i for i in range(-2, 5)]
           )
           + labs(x = "Days post infection",
                  y = "Viral titre (log10 TCID50/ml)")
           + theme_bw()
           + theme(legend_position = "none"))

data_p9.save("out/data-plot.png",
      height = 5.8, width = 8.3)

# --------------------------------------------------------------------
# Create data frames of values for each patient.
# --------------------------------------------------------------------

for patient_number, group_df in tcid_df.groupby("patient"):
    out_ssv = f"out/patient-{patient_number}.ssv"
    data_df = group_df[["day", "log10_tcid"]]
    data_df.rename(columns = {"day": "time", "log10_tcid": "value"},
                   inplace = True)
    data_df["time"] = data_df["time"].astype(float)
    data_df.to_csv(out_ssv, sep = " ", index = False)
