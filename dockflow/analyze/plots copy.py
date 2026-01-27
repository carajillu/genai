from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 15,                 # adjust to 8–10 depending on figure size
    "font.style": "normal",
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "normal",
    "mathtext.default": "it",        # italic math symbols (good for variables)
})


def hist_mode_center(
    df,
    key: str,
    bins: int = 100,
    ax: Optional[plt.Axes] = None,
    reference_value: Optional[float] = None,
    show_mode_line: bool = False,
    x_min: float = None,
    x_max: float = None, 
    label: float=None,
    marker: str="o",
) -> float:
    """
    Plot a smooth density line of df[key], and return the
    center of the bin with the highest count.
    """
    s = df[key].dropna().to_numpy()
    if s.size == 0:
        raise ValueError(f"Column '{key}' has no numeric data after dropping NaNs.")
    
    # Decide plotting range
    xmin_data, xmax_data = s.min(), s.max()

    xmin_plot = x_min if x_min is not None else xmin_data
    xmax_plot = x_max if x_max is not None else xmax_data

    # Compute histogram (counts, bin edges)
    counts, edges = np.histogram(s, bins=bins)
    mids = (edges[:-1] + edges[1:]) / 2.0
    bin_width = edges[1] - edges[0]

    # Mode bin center (from histogram)
    mode_idx = int(np.argmax(counts))
    mode_center = float((edges[mode_idx] + edges[mode_idx + 1]) / 2.0)

    ax = ax or plt.gca()

    # Optional mode line
    if show_mode_line:
        ax.axvline(mode_center, linestyle=":", linewidth=2)

    # Optional reference line
    if reference_value is not None:
        ax.axvline(reference_value, linestyle=":", linewidth=2,label="LSN", color="purple")

    
    from scipy.stats import gaussian_kde  # type: ignore
    xs = np.linspace(xmin_plot, xmax_plot, 512)
    density = gaussian_kde(s)(xs)
    ax.plot(xs, density, linewidth=5,label=label,marker=marker,markevery=0.05,ms=10)

    """
    except Exception:
        # Smoothed histogram fallback
        window = 5 if bins >= 5 else max(3, bins // 2 * 2 + 1)  # odd window
        kernel = np.ones(window) / window
        smooth = np.convolve(counts, kernel, mode="same")
        area = np.trapz(smooth, mids)
        ys = smooth / area if area > 0 else smooth
        
        ax.plot(mids, ys, linewidth=2,label=label,marker=marker,markevery=0.05,)
    """
    return mode_center



def get_overlapping_distriubutions(files: list[str], file_keys: list[str], plot_key: str, title_key: str, reference_value: float=None, x_min: float=None, x_max: float=None, scale: float=1.0):
    """
    Given a list of csv files, it loads each of those as a pandas dataframe (sep=",", header=True)
    It constructs a new df with the column plot_key of each dataframe, now relabled as file_keys[i]
    """
    assert len(files)==len(file_keys),"variables files and file_keys must have the same length"
    markers=["o","s","D","h","v"]
    fig, ax = plt.subplots(figsize=(6,5))
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)


    for i in range(len(files)):
        df_i=pd.read_csv(files[i])
        df_i=df_i[df_i.average_score!=np.inf] # clean failed dockings
        df_i[plot_key]=df_i[plot_key]*scale # just for CSA which we want in angstroms squared
        if i>0:
           reference_value=None
        hist_mode_center(df_i, key=plot_key, bins=100, ax=ax, reference_value=reference_value,x_min=x_min,x_max=x_max,label=file_keys[i],marker=markers[i])
    
    # Force x-limits if requested
    if (x_min is not None) or (x_max is not None):
        ax.set_xlim(x_min,x_max)

    ax.set_xlabel(title_key)
    ax.set_ylabel("Density")

    #ax.legend()
    plt.savefig(f"{plot_key}", dpi=300, bbox_inches="tight")
    plt.show()

    return
        
    

