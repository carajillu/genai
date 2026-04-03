from dockflow.analyze.plots import hist_mode_center, get_overlapping_distriubutions

import pandas as pd, matplotlib.pyplot as plt, pickle as pkl, numpy as np, glob, os

files=["round_1/docking/filter_bydock.csv",
       "round_2/docking/filter_bydock.csv",
       "round_3/docking/filter_bydock.csv",
       "round_4/docking/filter_bydock.csv",
       "round_5/docking/filter_bydock.csv"]

file_keys=["Iteration 1","Iteration 2","Iteration 3","Iteration 4","Iteration 5"]

# DOcking score
get_overlapping_distriubutions(files=files,file_keys=file_keys,plot_key="average_score",reference_value=None,x_min=-8.75,x_max=-6.75, title_key="Docking Score (kcal/mol)")
get_overlapping_distriubutions(files=files,file_keys=file_keys,plot_key="TanimotoDistance (raw)",reference_value=None,x_min=0.375,x_max=0.575, title_key="Tanimoto Distance")
get_overlapping_distriubutions(files=files,file_keys=file_keys,plot_key="SA score (raw)",reference_value=None,x_min=1.50,x_max=3.0, title_key="SA score")
get_overlapping_distriubutions(files=files,file_keys=file_keys,plot_key="QED (raw)",reference_value=None,x_min=0.725,x_max=1.00, title_key="QED")
get_overlapping_distriubutions(files=files,file_keys=file_keys,plot_key="average_csa",reference_value=None,x_min=520,x_max=630, title_key="Contact Surface Area ($\AA^{2}$)",scale=100)
