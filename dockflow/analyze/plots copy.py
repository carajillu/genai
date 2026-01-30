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
    

def kde_1d(
    data: np.ndarray,
    low: float,
    high: float,
    spacing: float,
    buffer: float = 0,
    bandwidth: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute a 1D KDE for `data` and return (x, density) on a fixed grid.

    Parameters
    ----------
    data : np.ndarray
        1D array-like of samples (NaNs/inf are ignored).
    low, high : float
        Inclusive bounds of the evaluation grid (x-axis).
    spacing : float
        Spacing between grid points (must be > 0).
    bandwidth : float, optional
        KDE bandwidth. If None, use library default (scipy) or Silverman's rule (fallback).
    buffer : float, optional
        Buffer to add to the low and high limits. If None, no buffer is added.

    Returns
    -------
    x : np.ndarray
        Grid points from low to high (inclusive where possible).
    y : np.ndarray
        KDE values at x; integrates to ~1 over the real line (the finite grid captures most mass).

    Raises
    ------
    ValueError
        If inputs are invalid or there is no finite data.
    """
    if spacing <= 0:
        raise ValueError("`spacing` must be > 0.")
    if high <= low:
        raise ValueError("`high` must be greater than `low`.")

    s = np.asarray(data).ravel()
    s = s[np.isfinite(s)]
    if s.size == 0:
        raise ValueError("No finite data in `data` after removing NaN/inf.")

    x = np.arange(low - buffer, high + buffer + spacing/2, spacing)  # include high when aligned

    # Try scipy's gaussian_kde first
    try:
        from scipy.stats import gaussian_kde  # type: ignore
        kde = gaussian_kde(s, bw_method=bandwidth)
        y = kde(x)
        return x, y
    except Exception:
        # Fallback: Gaussian KDE with Silverman's rule-of-thumb bandwidth
        print("Using fallback: Gaussian KDE with Silverman's rule-of-thumb bandwidth")
        n = s.size
        if bandwidth is None:
            if n > 1:
                std = s.std(ddof=1)
                iqr = np.subtract(*np.percentile(s, [75, 25]))
                sigma = min(std, iqr / 1.349) if (std > 0 and iqr > 0) else max(std, iqr / 1.349)
                h = 0.9 * sigma * n ** (-1/5) if sigma > 0 else (np.ptp(s) or 1.0) * n ** (-1/5) * 0.9
            else:
                h = 1.0  # arbitrary fallback for a single point
        else:
            h = float(bandwidth)
        h = max(h, np.finfo(float).eps)

        # Vectorized Gaussian kernels
        # y_j = (1 / (n*h*sqrt(2π))) * sum_i exp(-(x_j - s_i)^2 / (2h^2))
        diffs = (x[:, None] - s[None, :]) / h
        y = np.exp(-0.5 * diffs**2).sum(axis=1) / (n * h * np.sqrt(2 * np.pi))
        return x, y


def get_kde(dfs:list[pd.DataFrame],df_names:list[str],key:str,limits:list[float]=[None,None],spacing:float=0.01,buffer:float=0):
    """
    Get the KDE of a list of DataFrames.
    """
    # check that all dfs have a column named key
    for df in dfs:
        if key not in df.columns:
            raise ValueError(f"DataFrame {df} does not have a column named {key}")
    #if limits are not defined, get the maximum and minimum values across column {key} in all dataframes in {dfs}, then set them as the limits
    if limits[0] is None:
        limits[0] = min([df[key].min() for df in dfs])
    if limits[1] is None:
        limits[1] = max([df[key].max() for df in dfs])
    limits[0] = limits[0]-buffer
    limits[1] = limits[1]+buffer
    # create a list of x values
    x = np.arange(limits[0], limits[1], spacing)
    # create a list of y values
    kdes = pd.DataFrame()
    kdes["x"] = x
    for i in range(len(dfs)):
        df = dfs[i]
        kde_i = kde_1d(df[key].to_numpy(), limits[0], limits[1], spacing)[1]
        if len(kde_i)-len(x)==1:
            kde_i = kde_i[:-1]
        elif len(kde_i)-len(x)==-1:
            kde_i = np.concatenate([kde_i, [0]])
        elif len(kde_i) != len(x):
            raise ValueError(f"Length of KDE for {df_names[i]} ({len(kde_i)}) is not equal to length of x ({len(x)}).")
        kdes[df_names[i]] = kde_i
    return kdes

def plot_kde(kde_df: pd.DataFrame,line_width:float=2,colors:list[str]=None,save_name:str=None,xlabel:str="x",ylabel:str="Density",fontsize:float=20,legend: bool=True):
    if kde_df.keys()[0] != "x":
        raise ValueError("First column of kde_df must be 'x'. This probably isn't the right type of dataframe.")
    ax=kde_df.plot.line("x","round_1",linewidth=line_width,color=colors[0],fontsize=fontsize,legend=legend)
    for i in range(2,len(kde_df.keys())):
        kde_df.plot.line("x",kde_df.keys()[i],linewidth=line_width,color=colors[i],ax=ax,fontsize=fontsize,legend=legend)
    ax.set_xlabel(xlabel,fontweight="bold",fontsize=fontsize,labelpad=1)
    ax.set_ylabel(ylabel,fontweight="bold",fontsize=fontsize)
    if legend:
        plt.legend(fontsize=fontsize)
    if save_name is not None:
        plt.savefig(save_name,dpi=300,transparent=True)
    else:
        plt.show()

'''
Example usage:
import pandas as pd, matplotlib.pyplot as plt, glob, os, numpy as np; from dockflow.analyze.plots import hist_mode_center,get_kde,plot_kde
dfs,df_names=[],[]
for i in range(1,6):
    df_i=pd.read_csv(f"round_{i}/docking/filtered_bydock.csv")
    df_i.average_csa=df_i.average_csa*100
    dfs.append(df_i)
    df_names.append(f"round_{i}")
colors = [None,"#0072B2", "#E69F00", "#56B4E9", "#009E73", "#D55E00"]
tanimoto_kde=get_kde(dfs,df_names,key="TanimotoDistance (raw)",spacing=0.001,buffer=0.1)
plot_kde(tanimoto_kde,colors=colors,save_name="Tanimoto.png",xlabel="Tanimoto Distance",fontsize=13,line_width=4,legend=False)
'''

def find_highest(lst: list[float], n: int, reverse: bool = False) -> list[int]:
    """
    Returns the indices of the n highest values in lst.
    If reverse is True, returns the indices of the n lowest values in lst.
    """
    if n < 0:
        raise ValueError("n must be non-negative")

    n = min(n, len(lst))

    # (value, index), sort by value (desc for highest, asc for lowest)
    # tie-break by index for deterministic output
    pairs = sorted(
        enumerate(lst),
        key=lambda x: (x[1], x[0]) if reverse else (-x[1], x[0])
    )
    
    return [idx for idx, _ in pairs[:n]]
    

def pretty_plot(df,key_x,label_x,key_y,label_y,label_series,n_mark=None,mark_label=None,save_name=None):
    was_interactive = plt.isinteractive()
    plt.ioff()  # stop auto-draw while building the plot
    fig, ax = plt.subplots(figsize=(6,5))
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
    ax.plot(df[key_x],df[key_y],linewidth=1,label=label_series)
    if n_mark is not None:
       top_idx = find_highest(df[key_y],n_mark)
       print(top_idx)
       print(df[key_x].iloc[top_idx],df[key_y].iloc[top_idx])
       ax.plot(df[key_x][top_idx],df[key_y][top_idx],"o", color="red", markersize=5,label=mark_label)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    ax.set_xlim(0,500)
    ax.set_ylim(400,1000)
    ax.legend(loc="lower left")
    if save_name is not None:
       plt.savefig(save_name,  dpi=300, bbox_inches="tight") 
    plt.show()
    plt.close(fig)
    return

filename="contacts_3.csv"
rep_number=4
rep=pd.read_csv(filename)
rep["csa_angstrom2"]=rep["Total_Surface_Area_nm2"]*100
rep["time_ns"]=rep["Time_ps"]/1000
label_x="Time(ns)"
label_y="Contact Surface Area ($\AA^{2}$)"
plot_label=f"Replicate {rep_number}"
save_name=f"rep_{rep_number}"
pretty_plot(df=rep,key_x="time_ns",label_x=label_x,key_y="csa_angstrom2",label_y=label_y,label_series=plot_label,n_mark=9,mark_label="Selected for docking",save_name=save_name)

def compare_barplot(df, label_col="label", title="Comparison (holo vs apo)",
                    bar_width=0.5, save_name=None):

    df = df.copy()

    # residues in the same order as df columns
    residue_order = [c for c in df.columns if c != label_col]

    long = df.melt(id_vars=label_col, var_name="residue", value_name="value")
    long["residue"] = pd.Categorical(long["residue"], categories=residue_order, ordered=True)

    # pivot, then reindex rows to enforce the order
    wide = long.pivot(index="residue", columns=label_col, values="value").reindex(residue_order)

    ax = wide.plot(kind="bar", width=bar_width, figsize=(10, 5))

    ax.set_ylabel(r"RMSF ($\AA$)")
    ax.legend()
    plt.xticks(rotation=0, ha="center")
    plt.tight_layout()

    if save_name is not None:
        plt.savefig(save_name, dpi=300, bbox_inches="tight")

    plt.show()
    return wide