import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def qvalue_comparison(
    datasets,
    dataset_labels,
    q_label: str = "q",
    decoy_label: str = "is decoy",
    fdr_thresholds = None,
    log_scale: bool = True,
    title: str = "",
    ax = None,
):
    """
    Plot identification count in function of q-value threshold for multiple
    datasets.
    Parameters
    ----------
    datasets : List[pd.DataFrame]
        list of datasets in the form of `pandas.DataFrame`s
    dataset_labels : List[str]
        list of dataset labels to use in figure legend
    q_label : str, optional
        label of column containing each PSM's q-value, default `q`
    decoy_label : str, optional
        label of column marking decoy PSMs as True and target PSMs as False, default
        `is decoy`
    fdr_thresholds : List[float], optional
        list of FDR thresholds to plot as vertical, dotted lines
    log_scale : bool
        plot x-axis (q-values) in log scale or not.
    ax : matplotlib Axes, optional
        axes object to draw the plot onto, otherwise uses the current axes.
    Returns
    -------
    ax : matplotlib Axes
    """
    fig = plt.figure()
    if not isinstance(datasets, list):
        raise TypeError("`datasets` should be of type `list`.")
    if not datasets:
        raise ValueError("`datasets` cannot be empty.")

    if not fdr_thresholds:
        fdr_thresholds = [0.01, 0.001]

    if ax is None:
        ax = plt.gca()

    max_count = 0

    for i, rerec in enumerate(datasets):
        # Cumulatively count target IDs at each FDR threshold
        df = (
            rerec[~rerec[decoy_label]]
            .reset_index(drop=True)
            .sort_values(q_label, ascending=True)
            .copy()
        )
        df["count"] = (~df[decoy_label]).cumsum()

        # Plot counts
        if dataset_labels:
            label = dataset_labels[i]
        ax.plot(df[q_label], df["count"], label=label, alpha=0.5)

        # Get maximum count, required for vertical line height
        tmp_max = np.max(df["count"])
        if tmp_max > max_count:
            max_count = tmp_max

    # Plot FDR thresholds as dotted lines
    if fdr_thresholds:
        for fdr in fdr_thresholds:
            ax.plot(
                [fdr] * 2,
                np.linspace(0, max_count, 2),
                linestyle="--",
                color="black",
            )

    # Figure labels and legend
    ax.set_xlim(0.00001, 1)
    ax.set_ylabel("Number of identified spectra")
    ax.set_title(title)
    if log_scale:
        ax.set_xlabel("FDR threshold (log scale)")
        ax.set_xscale("log")
    else:
        ax.set_xlabel("FDR threshold")

    ax.legend()
    fig.set_size_inches(12, 10)

    return ax

def read_pout_file(target_pout, decoy_pout):
    """Read target and decoy pout files and combine into single pandas DataFrame."""

    target_pout = pd.read_csv(target_pout, sep="\t")
    decoy_pout = pd.read_csv(decoy_pout, sep="\t")
    target_pout["Label"] = 1
    decoy_pout["Label"] = -1
    pout = pd.concat([target_pout, decoy_pout])

    pout_qvalues = pout[["PSMId", "score", "q-value", "Label", "peptide"]].rename(
        columns={"q-value": "q", "Label": "is decoy"}
    )
    pout_qvalues["is decoy"] = pout["Label"] == -1

    return pout_qvalues