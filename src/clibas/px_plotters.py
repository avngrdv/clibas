"""
Interactive Plotly-based plotting utilities for UMAP/HDBSCAN analysis.

Provides interactive HTML dashboards with linked views for exploring
dimensionality reduction and clustering results. Includes UMAP scatter plots,
cluster statistics, frequency heatmaps, and sequence logos.
"""

import os

import numpy as np
import pandas as pd

import clibas.plotters as P

try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

except ImportError:
    msg = "Failed to import plotly packages. Please install clibas with `pip install clibas[ml]`. . ."
    raise ImportError(msg)


def _save_hmtl(fig, fname=None):
    """
    Save Plotly figure as HTML file.

    Args:
        fig: Plotly figure object.

        fname (str, optional): Output filename. '.html' added if not present.
    """
    if fname:
        if not fname.endswith(".html"):
            fname += ".html"

        fig.write_html(fname, include_plotlyjs="cdn", full_html=True)
    return


def _compute_marker_size(Q, cluster_mask=None, min_px=2, max_px=70, scaling=0.33):
    """
    Calculate marker sizes based on sequence counts.

    Args:
        Q (dict): Analysis results with 'C' (counts) and 'total_reads'.

        cluster_mask (ndarray, optional): Boolean mask for filtering.

        min_px (int): Minimum marker size in pixels. Default is 2.

        max_px (int): Maximum marker size in pixels. Default is 70.

        scaling (float): Power scaling factor. Default is 0.33.

    Returns:
        ndarray: Marker sizes in pixels.
    """
    C = Q["C"]
    if cluster_mask is not None:
        C = C[cluster_mask]

    f = C / Q["total_reads"]
    return min_px + (max_px - min_px) * np.power(f, scaling)


def _get_palette(clusters, bw=False):
    """
    Generate color palette for clusters.

    Args:
        clusters (ndarray): Cluster labels.

        bw (bool): If True, returns black/white palette. Default is False.

    Returns:
        dict: Mapping of cluster labels (as strings) to hex colors.
    """
    import seaborn as sns

    clusters = np.unique(clusters)
    palette = sns.color_palette("husl", clusters.size).as_hex()

    if 0 in clusters:
        palette[0] = "#323232"

    if bw:
        palette = ["#323232" for c in palette]

    return {
        str(cluster): color for cluster, color in zip(clusters, palette, strict=False)
    }


def _get_heatmap(X, labels, cluster=None, alphabet=None):
    """
    Compute frequency heatmap for sequences.

    Args:
        X (ndarray): Sequences (2D array).

        labels (ndarray): Cluster assignments.

        cluster (int, optional): If specified, compute for this cluster only.

        alphabet (array-like): Token alphabet.

    Returns:
        ndarray: Frequency matrix (n_tokens, sequence_length).
    """
    from clibas.misc import get_freqs

    if cluster is not None:
        freq = get_freqs(X[labels == cluster], alphabet=alphabet)
    else:
        freq = get_freqs(X, alphabet=alphabet)

    return freq


def umap_dashboard(Q, palette=None, cluster_mask=None):
    """
    Create interactive UMAP scatter plot.

    Generates Plotly scatter plot of UMAP embedding with hover information
    showing sequence rank, sequence, cluster, and count.

    Args:
        Q (dict): Analysis results with keys 'X', 'Y', 'labels', 'C', 'total_reads'.

        palette (dict, optional): Color mapping for clusters.

        cluster_mask (ndarray, optional): Boolean mask for filtering sequences.

    Returns:
        plotly.graph_objects.Figure: Interactive UMAP plot.
    """
    if cluster_mask is None:
        cluster_mask = np.ones(len(Q["X"]), dtype=bool)

    L = Q["X"].shape[1]
    seq = Q["X"][cluster_mask].view(f"S{L}").ravel().astype(f"U{L}")
    size = _compute_marker_size(Q, cluster_mask)

    if palette is None:
        palette = _get_palette(Q["labels"], bw=False)

    df = pd.DataFrame(
        {
            "umap1": Q["Y"][cluster_mask, 0],
            "umap2": Q["Y"][cluster_mask, 1],
            "seq": seq,
            "size": size,
            "cluster": Q["labels"][cluster_mask],
            "count": Q["C"][cluster_mask],
        }
    )
    df["cluster"] = df["cluster"].astype(str)
    df["rank"] = df["count"].rank(ascending=False, method="min").astype(int)

    fig = px.scatter(
        df,
        x="umap1",
        y="umap2",
        opacity=0.7,
        color="cluster",
        hover_data={"seq": False, "count": False, "size": False, "cluster": False},
        color_discrete_map=palette,
    )
    for trace in fig.data:
        cluster_name = trace.name
        cluster_df = df[df["cluster"] == cluster_name]
        trace.update(
            marker=dict(size=cluster_df["size"], line=dict(width=0)),
            customdata=cluster_df[["rank", "seq", "cluster", "count"]].values,
            hovertemplate='<b><span style="color: #905c54">%{customdata[0]}:</span> '
            '<span style="font-family: Consolas; color: #323232">%{customdata[1]}</span></b><br>'
            '<b style="color: #905c54">cluster:</b> '
            '<span style="color: #323232">%{customdata[2]}</span><br>'
            '<b style="color: #905c54">count:</b> '
            '<span style="color: #323232">%{customdata[3]}</span><br>'
            "<extra></extra>",
        )
    return fig


def barplot_dashboard(Q, palette=None, cluster=None, all_cluster_names=None):
    """
    Create interactive cluster statistics bar plot.

    Args:
        Q (dict): Analysis results with 'cluster_summary' DataFrame.

        palette (dict, optional): Color mapping for clusters.

        cluster (int, optional): If specified, highlights only this cluster.

        all_cluster_names (list, optional): All cluster names for x-axis.

    Returns:
        plotly.graph_objects.Figure: Interactive bar plot.
    """
    if cluster is None:
        df = pd.DataFrame(
            {
                "cluster": Q["cluster_summary"]["Cluster number"].astype(str),
                "score": Q["cluster_summary"]["Cluster score"],
                "purity": Q["cluster_summary"]["Cluster purity"],
                "size": Q["cluster_summary"]["Cluster size"],
            }
        )
    else:
        cluster_info = Q["cluster_summary"][
            Q["cluster_summary"]["Cluster number"] == cluster
        ].iloc[0]

        df = pd.DataFrame(
            {
                "cluster": all_cluster_names,
                "score": [
                    cluster_info["Cluster score"] if str(c) == str(cluster) else 0
                    for c in all_cluster_names
                ],
                "purity": [
                    cluster_info["Cluster purity"] if str(c) == str(cluster) else 0
                    for c in all_cluster_names
                ],
                "size": [
                    cluster_info["Cluster size"] if str(c) == str(cluster) else 0
                    for c in all_cluster_names
                ],
            }
        )

    fig = px.bar(
        df,
        x="cluster",
        y="score",
        color="cluster",
        color_discrete_map=palette,
        hover_data={"cluster": False, "score": False, "purity": False, "size": False},
    )
    for i, trace in enumerate(fig.data):
        cluster_val = df.iloc[i]["cluster"]
        cluster_data = df[df["cluster"] == cluster_val].iloc[0]

        trace.update(
            marker=dict(line=dict(width=0)),
            hovertemplate='<b style="color: #905c54">Cluster:</b> <span style="font-family: Arial; color: #323232">'
            + str(cluster_data["cluster"])
            + "</span><br>"
            + '<b style="color: #905c54">Score:</b> <span style="font-family: Arial; color: #323232">'
            + f"{cluster_data['score']:.2f}"
            + "</span><br>"
            + '<b style="color: #905c54">Purity:</b> <span style="font-family: Arial; color: #323232">'
            + f"{cluster_data['purity']:.2f}"
            + "</span><br>"
            + '<b style="color: #905c54">Size:</b> <span style="font-family: Arial; color: #323232">'
            + str(cluster_data["size"])
            + "</span><br>"
            + "<extra></extra>",
        )
    return fig


def freq_dashboard(freq, alphabet=None):
    """
    Create interactive frequency heatmap.

    Args:
        freq (ndarray): Frequency matrix (n_tokens, sequence_length).

        alphabet (list): Token labels.

    Returns:
        plotly.graph_objects.Figure: Interactive heatmap.
    """
    x_dim_size = freq.shape[1]

    if alphabet is not None and len(alphabet) > 0:
        if isinstance(alphabet[0], bytes):
            alphabet = [
                a.decode("utf-8") if isinstance(a, bytes) else a for a in alphabet
            ]

    fig = go.Figure(
        data=go.Heatmap(
            z=freq,
            x=list(range(1, x_dim_size + 1)),
            y=list(alphabet),
            colorscale="Blues",
            hovertemplate='<b style="color: #905c54">freq:</b> <span style="font-family: Consolas; color: #323232">%{z:.2f}</span><br>'
            + "<extra></extra>",
            colorbar=dict(
                len=200 / 600,
                thickness=15,
                x=1.02,
                y=1.025,
                yanchor="top",
                outlinecolor="#aaaaaa",
                outlinewidth=0.25,
                ticks="outside",
                tickcolor="#aaaaaa",
                ticklen=4,
                tickwidth=0.25,
            ),
        )
    )
    return fig


def four_panel_dashboard(Q, alphabet=None, fname=None):
    """
    Create comprehensive 4-panel interactive analysis dashboard.

    Generates an HTML dashboard with four linked panels:
    1. UMAP scatter plot (top-left)
    2. Frequency heatmap (top-right)
    3. Cluster score bar plot (bottom-left)
    4. Sequence logo (bottom-right)

    Includes dropdown menu to switch between overall view and individual
    cluster views. Sequence logos are generated as SVG files and embedded.

    Args:
        Q (dict): Analysis results containing:
            - 'X': Sequences (2D array)
            - 'Y': UMAP coordinates (n_sequences, 2)
            - 'labels': Cluster assignments
            - 'C': Sequence counts
            - 'cluster_summary': DataFrame with cluster statistics
            - 'total_reads': Total number of reads
            - 'name': Sample name

        alphabet (list): Token alphabet for plots.

        fname (str): Output HTML filename without extension.

    Returns:
        plotly.graph_objects.Figure: Four-panel dashboard figure.

    Note:
        Creates a 'logos' subdirectory in the output directory containing
        SVG sequence logos for overall and per-cluster views.

    Example:
        >>> four_panel_dashboard(
        ...     Q=analysis_results,
        ...     alphabet=['A','C','D','E','F','G'],
        ...     fname='output/sample_dashboard'
        ... )
    """
    clusters = np.unique(Q["labels"])
    palette = _get_palette(np.unique(Q["labels"]), bw=False)

    logo_dir = os.path.join(os.path.dirname(fname), "logos")
    if not os.path.isdir(logo_dir):
        os.makedirs(logo_dir)

    # overall data
    freq_data = {}
    freq_data["overall"] = _get_heatmap(
        Q["X"].astype("U1"), Q["labels"], cluster=None, alphabet=alphabet
    )
    P.Analysis.seq_logo(
        X=Q["X"].astype("U1"),
        C=Q["C"],
        freq=None,
        alphabet=alphabet,
        basename=os.path.join(logo_dir, "overall"),
    )
    # individual clusters
    for c in clusters:
        mask = Q["labels"] == c
        freq_data[f"cluster_{c}"] = _get_heatmap(
            Q["X"].astype("U1"), Q["labels"], cluster=c, alphabet=alphabet
        )
        P.Analysis.seq_logo(
            X=Q["X"][mask].astype("U1"),
            C=Q["C"][mask],
            freq=None,
            alphabet=alphabet,
            basename=os.path.join(logo_dir, f"cluster_{c}"),
        )
    fig = make_subplots(
        rows=2,
        cols=2,
        row_heights=[600, 200],
        column_widths=[0.5, 0.5],
        vertical_spacing=0.1,
        horizontal_spacing=0.1,
        specs=[
            [{"type": "scatter"}, {"type": "heatmap"}],
            [{"type": "bar"}, {"type": "xy"}],
        ],
    )
    # store umap x,y lims
    umap = umap_dashboard(Q, palette=palette)
    umap_xlim = [Q["Y"][:, 0].min(), Q["Y"][:, 0].max()]
    umap_ylim = [Q["Y"][:, 1].min(), Q["Y"][:, 1].max()]

    # add umap traces (overall)
    for trace in umap.data:
        fig.add_trace(trace, row=1, col=1)

    # add umap traces (by cluster)
    for cluster in clusters:
        cluster_mask = Q["labels"] == cluster
        cluster_umap = umap_dashboard(Q, palette=palette, cluster_mask=cluster_mask)

        for trace in cluster_umap.data:
            trace.visible = False
            fig.add_trace(trace, row=1, col=1)

    bars = barplot_dashboard(Q, palette=palette)
    all_cluster_names = [str(c) for c in clusters]

    for trace in bars.data:
        fig.add_trace(trace, row=2, col=1)

    for cluster in clusters:
        cluster_bar = barplot_dashboard(
            Q, palette=palette, cluster=cluster, all_cluster_names=all_cluster_names
        )
        for trace in cluster_bar.data:
            trace.visible = False
            fig.add_trace(trace, row=2, col=1)

    x_dim_size = freq_data["overall"].shape[1]
    for key in ["overall"] + [f"cluster_{c}" for c in clusters]:
        freq = freq_data[key]
        visible = key == "overall"

        heatmap_trace = go.Heatmap(
            z=freq,
            x=list(range(1, x_dim_size + 1)),
            y=list(alphabet),
            colorscale="Blues",
            visible=visible,
            hovertemplate='<b style="color: #905c54">freq:</b> <span style="font-family: Consolas; color: #323232">%{z:.2f}</span><br>'
            + "<extra></extra>",
            colorbar=dict(
                len=200 / 600,
                thickness=15,
                x=1.005,
                y=1.02,
                yanchor="top",
                outlinecolor="#aaaaaa",
                outlinewidth=0.25,
                ticks="outside",
                tickcolor="#aaaaaa",
                ticklen=4,
                tickwidth=0.25,
            ),
        )
        fig.add_trace(heatmap_trace, row=1, col=2)

    for i, key in enumerate(["overall"] + [f"cluster_{c}" for c in clusters]):
        svg_path = f"logos/{key}.svg"
        visible = key == "overall"

        fig.add_layout_image(
            dict(
                source=svg_path,
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.2,
                sizex=0.5,
                sizey=0.25,
                xanchor="left",
                yanchor="top",
                visible=visible,
                name=f"logo_{key}",
            )
        )
    n_umap_overall = len(umap.data)
    n_umap_per_cluster = len(cluster_umap.data)
    n_bar_overall = len(bars.data)
    n_bar_per_cluster = len(cluster_bar.data)
    n_clusters = len(clusters)

    def make_vis_list(selected_idx):
        visible = []

        # UMAP traces
        visible.extend(
            [True] * n_umap_overall if selected_idx == -1 else [False] * n_umap_overall
        )
        for i in range(n_clusters):
            visible.extend(
                [True] * n_umap_per_cluster
                if i == selected_idx
                else [False] * n_umap_per_cluster
            )

        # bar traces
        visible.extend(
            [True] * n_bar_overall if selected_idx == -1 else [False] * n_bar_overall
        )
        for i in range(n_clusters):
            visible.extend(
                [True] * n_bar_per_cluster
                if i == selected_idx
                else [False] * n_bar_per_cluster
            )

        # freq heatmap traces
        visible.extend([i == (selected_idx + 1) for i in range(n_clusters + 1)])
        return visible

    def make_img_vis(selected_idx):
        return [i == (selected_idx + 1) for i in range(n_clusters + 1)]

    def upd_img_vis(visible_list):
        return [
            dict(img.to_plotly_json(), visible=vis)
            for img, vis in zip(fig.layout.images, visible_list, strict=False)
        ]

    buttons = []
    buttons.append(
        dict(
            label="Overall",
            method="update",
            args=[
                {"visible": make_vis_list(-1)},
                {"images": upd_img_vis(make_img_vis(-1))},
            ],
        )
    )
    for i, cluster in enumerate(clusters):
        buttons.append(
            dict(
                label=f"Cluster {cluster}",
                method="update",
                args=[
                    {"visible": make_vis_list(i)},
                    {"images": upd_img_vis(make_img_vis(i))},
                ],
            )
        )
    # umap axes
    fig.update_xaxes(
        showgrid=False,
        showticklabels=False,
        title="umap1",
        showline=True,
        linewidth=1,
        linecolor="#aaaaaa",
        mirror=True,
        range=umap_xlim,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        showgrid=False,
        showticklabels=False,
        title="umap2",
        showline=True,
        linewidth=1,
        linecolor="#aaaaaa",
        mirror=True,
        range=umap_ylim,
        scaleanchor="x1",
        scaleratio=1,
        row=1,
        col=1,
    )
    # barplot axes
    fig.update_xaxes(
        title="Cluster (hover for details)",
        showgrid=False,
        showline=True,
        linewidth=0.25,
        linecolor="#aaaaaa",
        mirror=True,
        title_font=dict(size=14),
        categoryorder="array",
        categoryarray=all_cluster_names,
        row=2,
        col=1,
    )
    fig.update_yaxes(
        title="Cluster score",
        showgrid=False,
        showline=True,
        linewidth=0.25,
        linecolor="#aaaaaa",
        mirror=True,
        range=[0, Q["cluster_summary"]["Cluster score"].max() + 0.5],
        title_font=dict(size=14),
        row=2,
        col=1,
    )
    # frequency heatmap axes
    fig.update_xaxes(
        title="Position",
        side="bottom",
        showgrid=False,
        showline=True,
        linewidth=0.25,
        linecolor="#aaaaaa",
        mirror=True,
        title_font=dict(size=14),
        ticks="outside",
        tickmode="linear",
        tick0=1,
        dtick=1,
        tickcolor="#aaaaaa",
        ticklen=4,
        tickwidth=0.25,
        row=1,
        col=2,
    )
    fig.update_yaxes(
        title="Token",
        showgrid=False,
        showline=True,
        linewidth=0.25,
        linecolor="#aaaaaa",
        mirror=True,
        title_font=dict(size=14),
        ticks="outside",
        tickmode="linear",
        tick0=1,
        dtick=1,
        tickcolor="#aaaaaa",
        ticklen=4,
        tickwidth=0.25,
        row=1,
        col=2,
    )
    # logo image axes
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        row=2,
        col=2,
        visible=False,
    )
    fig.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        scaleanchor="x4",
        visible=False,
        scaleratio=1,
        row=2,
        col=2,
    )
    fig.update_layout(
        height=800,
        width=1100,
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="Arial"),
        showlegend=False,
        hoverlabel=dict(bgcolor="white", font_family="Arial"),
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=0.44,
                xanchor="left",
                y=1.15,
                yanchor="top",
                bgcolor="white",
                bordercolor="#aaaaaa",
                borderwidth=0.25,
                font=dict(size=12),
            )
        ],
    )
    annotations = [
        dict(
            text="UMAP projection",
            xref="x domain",
            yref="y domain",
            x=0,
            y=1.05,
            showarrow=False,
            font=dict(size=14, family="Arial", color="#323232", weight="bold"),
            row=1,
            col=1,
        ),
        dict(
            text="Cluster scores",
            xref="x domain",
            yref="y domain",
            x=0,
            y=1.14,
            showarrow=False,
            font=dict(size=14, family="Arial", color="#323232", weight="bold"),
            row=2,
            col=1,
        ),
        dict(
            text="Sequence conservation (unweighed)",
            xref="x domain",
            yref="y domain",
            x=0,
            y=1.05,
            showarrow=False,
            font=dict(size=14, family="Arial", color="#323232", weight="bold"),
            row=1,
            col=2,
        ),
    ]
    for annot in annotations:
        row = annot.pop("row")
        col = annot.pop("col")
        fig.add_annotation(annot, row=row, col=col)

    x0, x1 = fig.layout["xaxis4"].domain
    y0, y1 = fig.layout["yaxis4"].domain

    fig.add_annotation(
        text="Sequence logo (weighed)",
        xref="paper",
        yref="paper",
        x=(x0 + x1) / 2 - 0.065,
        y=y1 - 0.042,
        showarrow=False,
        font=dict(size=14, family="Arial", color="#323232", weight="bold"),
    )
    _save_hmtl(fig, fname=fname)
    return fig


def single_manifold_embedding_dashboard(tup, fname=None, cols=2):
    """
    Create multi-sample UMAP comparison dashboard.

    Generates side-by-side UMAP plots for multiple samples embedded on a
    single shared manifold, enabling direct visual comparison of different
    datasets.

    Args:
        tup (list): List of analysis result dictionaries, each containing:
            - 'X': Sequences
            - 'Y': UMAP coordinates (on shared manifold)
            - 'labels': Cluster assignments
            - 'C': Sequence counts
            - 'name': Sample name

        fname (str, optional): Output HTML filename without extension.

        cols (int): Number of columns in grid layout. Default is 2.

    Returns:
        plotly.graph_objects.Figure: Multi-panel comparison dashboard.

    Note:
        All samples should be embedded using single_manifold=True in the
        umap_hdbscan_analysis operation to ensure coordinates are comparable.

    Example:
        >>> #sfter running analysis with single_manifold=True
        >>> single_manifold_embedding_dashboard(
        ...     tup=sample_results_list,
        ...     fname='output/manifold_comparison',
        ...     cols=2
        ... )
    """
    n = len(tup)
    rows = int(np.ceil(n / cols))

    x_min = [d["Y"][:, 0].min() for d in tup]
    x_max = [d["Y"][:, 0].max() for d in tup]
    y_min = [d["Y"][:, 1].min() for d in tup]
    y_max = [d["Y"][:, 1].max() for d in tup]

    x_range = (min(x_min) - 0.25, max(x_max) + 0.25)
    y_range = (min(y_min) - 0.25, max(y_max) + 0.25)

    fig = make_subplots(
        rows=rows,
        cols=cols,
        horizontal_spacing=0.02,
        vertical_spacing=0.05,
        subplot_titles=[d["name"] for i, d in enumerate(tup)],
    )

    for i, d in enumerate(tup):
        row = i // cols + 1
        col = i % cols + 1

        palette = _get_palette(d["labels"], bw=True)
        umap_fig = umap_dashboard(d, palette=palette)

        for trace in umap_fig.data:
            fig.add_trace(trace, row=row, col=col)

        fig.update_xaxes(
            showgrid=False,
            showticklabels=False,
            range=x_range,
            linecolor="#aaaaaa",
            linewidth=0.5,
            scaleanchor=f"y{row}{col}",
            scaleratio=1,
            row=row,
            mirror=True,
            col=col,
        )
        fig.update_yaxes(
            showgrid=False,
            showticklabels=False,
            range=y_range,
            linecolor="#aaaaaa",
            linewidth=0.5,
            row=row,
            mirror=True,
            col=col,
        )
    panel_size = 600
    total_width = panel_size * cols * 0.93
    total_height = panel_size * rows

    fig.update_layout(
        height=total_height,
        width=total_width,
        margin=dict(l=20, r=10, t=60, b=20),
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="Arial", color="#323232", weight="bold", size=14),
        showlegend=False,
        hoverlabel=dict(bgcolor="white", font_family="Arial"),
    )
    _save_hmtl(fig, fname=fname)
    return fig
