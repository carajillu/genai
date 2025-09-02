from typing import Optional
import numpy as np
import matplotlib.pyplot as plt

def hist_mode_center(
    df,
    key: str,
    bins: int = 100,
    ax: Optional[plt.Axes] = None,
    save_name: str=None
) -> float:
    """
    Plot a histogram of df[key] with `bins` bins, overlay a smooth density line,
    and return the center of the bin with the highest count.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.
    key : str
        Column name to analyze.
    bins : int, default=100
        Number of histogram bins.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates a new figure/axes.

    Returns
    -------
    float
        The value at the center of the most populated bin.
    """
    s = df[key].dropna().to_numpy()
    if s.size == 0:
        raise ValueError(f"Column '{key}' has no numeric data after dropping NaNs.")

    # Compute histogram (counts, bin edges)
    counts, edges = np.histogram(s, bins=bins)
    mids = (edges[:-1] + edges[1:]) / 2.0
    bin_width = edges[1] - edges[0]

    # Mode bin center
    mode_idx = int(np.argmax(counts))
    mode_center = float((edges[mode_idx] + edges[mode_idx + 1]) / 2.0)

    # Plot
    ax = ax or plt.gca()
    ax.hist(s, bins=bins, alpha=0.4, edgecolor="black", linewidth=0.5)

    # Try KDE (SciPy); fall back to smoothed histogram if SciPy missing
    try:
        from scipy.stats import gaussian_kde  # type: ignore
        xs = np.linspace(s.min(), s.max(), 512)
        kde = gaussian_kde(s)
        density = kde(xs)  # integrates to 1 over x
        # Scale density to histogram counts for visual comparability
        scaled = density * s.size * bin_width
        ax.plot(xs, scaled, linewidth=2)
    except Exception:
        # Simple moving-average smoothing of counts
        window = 5 if bins >= 5 else max(3, bins // 2 * 2 + 1)  # odd window
        kernel = np.ones(window) / window
        smooth = np.convolve(counts, kernel, mode="same")
        ax.plot(mids, smooth, linewidth=2)

    # Mark the mode bin center
    ax.axvline(mode_center, linestyle="--")
    ax.set_xlabel(key)
    ax.set_ylabel("Count")
    ax.set_title(f"Histogram & mode bin center of '{key}'")
    plt.show()
    if save_name is not None:
       plt.savefig(save_name,dpi=300)

    return mode_center
