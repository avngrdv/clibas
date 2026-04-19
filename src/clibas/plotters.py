"""
Matplotlib-based plotting utilities for sequencing data analysis.

Provides static plot generation for sequence length distributions, quality
scores, library convergence, UMAP embeddings, and other analysis outputs.
Plots are saved as PNG and SVG files.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial"]


def _save_figs(fig, basename=None):
    """Save figure as PNG and SVG files.

    Args:
        fig: Matplotlib figure object.
        basename (str, optional): Output path without extension.
    """
    if basename is not None:
        # save png and svg, and close the file
        svg = basename + ".svg"
        png = basename + ".png"
        fig.savefig(svg, bbox_inches="tight")
        fig.savefig(png, bbox_inches="tight")
        plt.close()
    return


# Just a container for FastqParser-related plotters
class SequencingData:
    """
    Static plots for FASTQ sequencing data analysis.

    Collection of plotting functions for visualizing sequencing data metrics
    including length distributions, quality scores, conservation, and token
    frequencies.
    """

    def L_distribution(X, Y, where=None, basename=None):
        """
        Plot sequence length distribution histogram.

        Args:
            X (ndarray): Sequence lengths (x-axis values).

            Y (ndarray): Counts for each length (y-axis values).

            where (str, optional): Dataset type ('dna' or 'pep') for labeling.

            basename (str, optional): Output path without extension.

        Example:
            >>> L, counts = np.unique(lengths, return_counts=True)
            >>> SequencingData.L_distribution(L, counts, where='pep',
            ...                               basename='output/length_dist')
        """
        if not where:
            where = ""

        fig = plt.figure(figsize=(9, 3), dpi=300)
        ax = fig.add_subplot(111)
        plt.bar(X, Y, color="#1571da")

        ax.set_ylim(0, 1.02 * np.max(Y))
        ax.set_xlim(np.min(X), np.max(X) + 1)
        ax.set_xticks(np.linspace(np.min(X), np.max(X) + 1, 10))
        ax.set_xticklabels(np.linspace(np.min(X), np.max(X) + 1, 10, dtype=int))

        ax.set_xlabel("Sequence length", fontsize=14)
        ax.set_ylabel("Count", fontsize=14)
        ax.tick_params(axis="both", which="major", labelsize=12)

        title = f"Distribution of sequence lengths in {where} dataset"
        ax.set_title(title, fontsize=14, y=1.04)

        _save_figs(fig, basename=basename)
        return

    def dataset_convergence(C, shannon, where, basename):
        """
        Plot sequence-level convergence with Shannon entropy.

        Generates a log-scale plot showing cumulative sequence count vs
        percentile, annotated with normalized Shannon entropy.

        Args:
            C (ndarray): Sorted sequence counts.

            shannon (float): Normalized Shannon entropy value.

            where (str): Dataset type ('dna' or 'pep').

            basename (str): Output path without extension.
        """
        y = np.sort(C)
        x = 100 * np.divide(np.arange(y.size), y.size)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
        plt.plot(x, y, lw=2.5, c="#1571da", antialiased=True)

        ax.set_xlim(0, 101)
        ax.set_xticks(np.arange(0, 125, 25))

        ax.set_yscale("log")

        ax.set_ylabel(f"{where} sequence count", fontsize=14)
        ax.set_xlabel("Sequence percentile", fontsize=14)
        ax.set_title(f"Sequence-level convergence of {where} dataset", fontsize=16)
        plt.text(
            x=2,
            y=y.max(),
            size=12,
            s=f"normalized Shannon entropy: {shannon:1.4f}",
            horizontalalignment="left",
            verticalalignment="center",
        )
        plt.grid(
            lw=0.5,
            ls="--",
            c="slategrey",
            alpha=0.2,
            dash_capstyle="round",
            dash_joinstyle="round",
            antialiased=True,
        )
        _save_figs(fig, basename=basename)
        return

    def conservation(conservation, where, basename):
        """
        Plot position-wise sequence conservation. Deprecated.

        Args:
            conservation (ndarray): Conservation values per position.

            where (str): Dataset type ('dna' or 'pep').

            basename (str): Output path without extension.
        """
        fig = plt.figure(figsize=(9, 3), dpi=300)
        ax = fig.add_subplot(111)

        ax.plot(np.arange(len(conservation)) + 1, conservation, lw=3.5, c="#1571da")

        ax.set_ylim(0, np.ceil(np.max(conservation)))
        ax.set_xlim(0.5, len(conservation) + 0.5)

        from matplotlib.ticker import MaxNLocator

        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True, prune=None))

        ax.set_xlabel("Position index", fontsize=14)
        ax.set_ylabel("Conservation, norm", fontsize=14)
        ax.tick_params(axis="both", which="major", labelsize=12)

        title = f"Position-wise sequence conservation plot for {where} dataset"
        ax.set_title(title, fontsize=14, y=1.04)

        _save_figs(fig, basename=basename)
        return

    def tokenwise_frequency(freq, yticknames, where=None, loc=None, basename=None):
        """
        Plot position-wise token frequency heatmap.

        Generates a heatmap showing frequency of each token (amino acid, base)
        at each position in sequences.

        Args:
            freq (ndarray): Frequency matrix (n_tokens, sequence_length).

            yticknames (list): Token labels for y-axis.

            where (str, optional): Dataset type for labeling.

            loc (str or list, optional): Region specification for labeling.

            basename (str, optional): Output path without extension.
        """
        if not where:
            where = ""

        if where == "dna":
            ylabel = "Base"
        elif where == "pep":
            ylabel = "Amino acid"
        else:
            ylabel = "Token"

        figsize = (1 + freq.shape[1] / 2, freq.shape[0] / 2)
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300)

        norm = mpl.colors.Normalize(vmin=0, vmax=np.max(freq))
        c = ax.pcolormesh(
            freq, cmap=plt.cm.Blues, norm=norm, edgecolors="w", linewidths=4
        )
        cbar = fig.colorbar(c, ax=ax)
        cbar.ax.set_ylabel("Frequency", rotation=-90, va="bottom", fontsize=22)
        cbar.ax.tick_params(labelsize=20)

        # set ticks
        ax.set_xticks(np.arange(freq.shape[1]) + 0.5)
        ax.set_yticks(np.arange(freq.shape[0]) + 0.5)
        ax.set_xticklabels(np.arange(freq.shape[1]) + 1)
        ax.set_yticklabels(yticknames)

        # set labels
        if loc is not None:
            ax.set_xlabel(f"Position inside library region(s) {loc}", fontsize=21)

        ax.set_ylabel(ylabel, fontsize=21)
        ax.tick_params(axis="both", which="major", labelsize=16)
        ax.set_title(f"Position-wise frequency map for {where} dataset", fontsize=25)

        _save_figs(fig, basename=basename)
        return

    def Q_score_summary(avg, std, loc, basename):
        """
        Plot quality score statistics across positions.

        Plots mean quality score with standard deviation bands.

        Args:
            avg (ndarray): Average quality scores per position.

            std (ndarray): Standard deviation per position.

            loc (str): Region specification for labeling.

            basename (str): Output path without extension.
        """
        fig = plt.figure(figsize=(18, 6), dpi=300)
        ax = fig.add_subplot(111)
        plt.plot(avg, lw=4, c="#3b61b1")
        plt.plot(avg + std, lw=1, c="#0091b5")
        plt.plot(avg - std, lw=1, c="#0091b5")
        ax.fill_between(
            np.arange(len(avg)), avg - std, avg + std, color="#0091b5", alpha=0.15
        )

        ax.set_ylim(0, np.nanmax(avg + std) + 5)
        ax.set_xlim(0, avg.size)
        ax.set_xticks(np.linspace(0, avg.size, 10))
        ax.set_xticklabels(np.linspace(1, avg.size + 1, 10, dtype=int))

        ax.set_xlabel(f"{loc} region(s) index", fontsize=30)
        ax.tick_params(axis="both", which="major", labelsize=25)
        ax.set_ylabel("Q, average log score", fontsize=30)

        title = "Q-score plot"
        ax.set_title(title, fontsize=34, y=1.04)

        _save_figs(fig, basename=basename)
        return


class Analysis:
    """
    Static plots for data analysis outputs.

    Plotting functions for dimensionality reduction, clustering, and
    sequence logos.
    """

    def UMAP_HDBSCAN(
        Y,
        labels,
        C=None,
        sample_name=None,
        basename=None,
        total_reads=None,
        show_annotations=False,
    ):
        """
        Plot UMAP embedding with HDBSCAN cluster assignments.

        Creates a 2D scatter plot of UMAP coordinates colored by cluster,
        with point sizes proportional to sequence counts.

        Args:
            Y (ndarray): UMAP coordinates (n_sequences, 2).

            labels (ndarray): Cluster assignments (0 for noise).

            C (ndarray, optional): Sequence counts for sizing points.

            total_reads(int, optional): Total number of reads in the sample.

            sample_name (str, optional): Sample name for title.

            basename (str, optional): Output path without extension.

            show_annotations (bool): If True, annotate cluster centers with
                cluster numbers. Default is False.
        """
        from matplotlib.colors import ListedColormap

        cmap = ListedColormap(sns.color_palette("husl", labels.max() + 5))

        fig = plt.figure(figsize=(8, 8), dpi=300)
        ax = fig.add_subplot(111)

        # Replicate _compute_marker_size logic exactly
        min_px = 2
        max_px = 75
        scaling = 0.33

        if C is not None:
            if total_reads is None:
                total_reads = C.sum()  # fallback
            f = C / total_reads
            sizes = min_px + (max_px - min_px) * np.power(f, scaling)
            sizes = sizes**2
        else:
            sizes = 5500 * np.power(np.divide(1, Y.shape[0]), 0.53)

        if not sample_name:
            sample_name = "unnamed sample"

        noise = labels == 0
        plt.scatter(
            Y[:, 0][noise][::-1],
            Y[:, 1][noise][::-1],
            alpha=0.7,
            c="#070D0D",
            marker="o",
            edgecolors="none",
            s=sizes[noise][::-1],
        )

        aint_noise = labels != 0
        plt.scatter(
            Y[:, 0][aint_noise][::-1],
            Y[:, 1][aint_noise][::-1],
            alpha=0.7,
            c=labels[aint_noise][::-1],
            cmap=cmap,
            marker="o",
            edgecolors="none",
            s=sizes[aint_noise][::-1],
        )

        if show_annotations:
            for cluster in labels:
                if not cluster == -1:
                    x_coord = np.average(Y[:, 0][labels == cluster])
                    y_coord = np.average(Y[:, 1][labels == cluster])
                    plt.text(
                        x=x_coord,
                        y=y_coord,
                        s=f"{cluster}",
                        size=15,
                        weight="bold",
                        alpha=0.3,
                        color="#323232",
                        horizontalalignment="center",
                        verticalalignment="center",
                    )
        plt.axis("off")

        title = f"umap/hdbscan: {sample_name}"
        ax.set_title(title, fontsize=20, y=1.04)

        _save_figs(fig, basename=basename)
        return

    def ClusteringHyperParams(
        min_clusters, min_samples, scores, sample_name=None, basename=None
    ):
        """
        Plot HDBSCAN hyperparameter optimization results.

        Generates heatmap showing clustering scores across different
        min_cluster_size and min_samples combinations.

        Args:
            min_clusters (ndarray): min_cluster_size values tested.

            min_samples (ndarray): min_samples values tested.

            scores (ndarray): Clustering scores (2D grid).

            sample_name (str, optional): Sample name for title.

            basename (str, optional): Output path without extension.
        """
        fig = plt.figure(figsize=(7, 6), dpi=300)
        ax = fig.add_subplot(111)
        c = plt.pcolor(scores, cmap=sns.color_palette("mako", as_cmap=True))
        fig.colorbar(c, ax=ax)

        ax.set_xlabel("min_sample", fontsize=16)
        ax.set_ylabel("min_cluster", fontsize=16)

        ax.set_xticks(np.arange(len(min_samples)) + 0.5)
        ax.set_yticks(np.arange(len(min_clusters)) + 0.5)

        ax.set_xticklabels(min_samples)
        ax.set_yticklabels(min_clusters)

        ax.tick_params(axis="both", which="major", labelsize=14)

        title = f"hdbscan clustering scores: {sample_name}"
        ax.set_title(title, fontsize=18, y=1.02)

        _save_figs(fig, basename=basename)
        return

    def seq_logo(X=None, C=None, freq=None, alphabet=None, basename=None):
        """
        Generate sequence logo from sequences or frequency matrix.

        Creates information content-based sequence logo using logomaker.
        Can generate from raw sequences with counts or pre-computed frequencies.

        Args:
            X (ndarray, optional): Sequences (2D array). Required if freq=None.

            C (ndarray, optional): Sequence counts for weighting. Required if freq=None.

            freq (DataFrame, optional): Pre-computed frequency matrix. If provided,
                X and C are ignored.

            alphabet (list): Token alphabet for columns.

            basename (str, optional): Output path without extension.

        Note:
            Either provide (X, C, alphabet) or (freq, alphabet), not both.

        Example:
            >>> #from sequences
            >>> seq_logo(X=sequences, C=counts, alphabet=['A','C','D','E'],
            ...          basename='output/logo')
            >>>
            >>> #from frequency matrix
            >>> seq_logo(freq=freq_df, alphabet=['A','C','D','E'],
            ...          basename='output/logo')
        """
        import logomaker

        from clibas.misc import freq_to_information_content

        if freq is None:
            # compute weighted frequency matrix
            seq_len = X.shape[-1]
            freq = pd.DataFrame(0.0, index=range(seq_len), columns=alphabet)
            for seq, w in zip(X, C, strict=False):
                for i, aa in enumerate(seq):
                    if aa:
                        freq.loc[i, aa] += w

            freq = freq.div(freq.sum(axis=1), axis=0)

        # this list covers both aas and bases; I think should be robust enough?
        token_colors = {
            "A": "#64b964",
            "C": "#ffed6f",
            "D": "#e31a1c",
            "E": "#ff7f00",
            "F": "#6a3d9a",
            "G": "#b15928",
            "H": "#1f78b4",
            "I": "#33a02c",
            "K": "#fb9a99",
            "L": "#cab2d6",
            "M": "#ffff99",
            "N": "#a6cee3",
            "P": "#fdbf6f",
            "Q": "#b2df8a",
            "R": "#1f78b4",
            "S": "#e31a1c",
            "T": "#ff7f00",
            "V": "#33a02c",
            "W": "#6a3d9a",
            "Y": "#6a3d9a",
        }
        # fallback for custom aas - just color them black:
        for col in alphabet:
            if col not in token_colors:
                token_colors[col] = "#323232"

        # make the logo
        IC = freq_to_information_content(freq).fillna(0)
        logo = logomaker.Logo(IC, font_name="Consolas", color_scheme=token_colors)

        logo.ax.set_xlabel("Position", fontsize=14)
        logo.ax.set_ylabel("Information (bits)", fontsize=14)

        logo.ax.patch.set_edgecolor("#aaaaaa")
        logo.ax.patch.set_linewidth(0.25)

        for axis in ["top", "bottom", "left", "right"]:
            logo.ax.spines[axis].set_linewidth(0.25)
            logo.ax.spines[axis].set_color("#aaaaaa")

        logo.ax.tick_params(
            axis="both", which="major", width=0.25, color="#aaaaaa", labelsize=12
        )

        _save_figs(logo.fig, basename=basename)
        return


class Miscellaneous:
    def single_linkage_dendrogram(link, labels=None, basename=None):
        """
        Plot hierarchical clustering dendrogram using single linkage.

        Generates a dendrogram visualization from a linkage matrix, typically
        derived from hierarchical clustering (e.g., Ward linkage). Displays
        hierarchical relationships among tokens or sequence features.

        Args:
            link (ndarray): Linkage matrix from hierarchical clustering.

            labels (list, optional): Labels for the dendrogram leaves.

            basename (str, optional): Output path without extension for saving plots.

        Example:
            >>> from scipy.cluster.hierarchy import linkage
            >>> Z = linkage(X, method='ward')
            >>> Miscellaneous.single_linkage_dendrogram(Z, labels=feature_names,
            ...                                         basename='output/dendrogram')
        """

        from scipy.cluster.hierarchy import dendrogram

        dims = ((link.shape[0] + 1) * 0.3, 3)
        fig = plt.figure(figsize=dims, dpi=300)
        ax = fig.add_subplot(111)
        dendrogram(link, labels=labels, p=50, ax=ax)

        ax.set_xlabel("Tokens", fontsize=16)
        ax.set_ylabel("Height", fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=14)

        title = "Ward linkage dendrogram for the feature matrix"
        ax.set_title(title, fontsize=18, y=1.02)

        if basename is not None:
            # save png and svg, and close the file
            svg = basename + ".svg"
            png = basename + ".png"
            fig.savefig(svg, bbox_inches="tight")
            fig.savefig(png, bbox_inches="tight")
            plt.close()

        return
