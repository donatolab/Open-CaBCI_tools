import argparse
import csv
import json
import os
import sys
import tempfile
import textwrap
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib")
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import networkx as nx
import numpy as np
import seaborn as sns

from bmitools import Calcium
from bmitools.utils import compute_dff0
from bmitools.utils.correlation.correlation import get_corr2
from bmitools.utils.network.network import (
    generate_graph_from_connected_nodes,
    get_degree_distribution,
    load_correlation_matrix,
)


DEFAULT_DATA = (
    REPOSITORY_ROOT
    / "test_data"
    / "DON-014371"
    / "20230411Rec"
    / "calibration"
    / "suite2p"
    / "plane0"
)
DEFAULT_OUTPUT = REPOSITORY_ROOT / "examples" / "outputs"

TOKENS = {
    "surface": "#FCFCFD",
    "panel": "#FFFFFF",
    "ink": "#1F2430",
    "muted": "#6F768A",
    "grid": "#E6E8F0",
    "axis": "#D7DBE7",
}
BLUE = {
    "xlight": "#EAF1FE",
    "light": "#CEDFFE",
    "base": "#A3BEFA",
    "mid": "#5477C4",
    "dark": "#2E4780",
}
ORANGE = {
    "xlight": "#FFEDDE",
    "light": "#FFBDA1",
    "base": "#F0986E",
    "mid": "#CC6F47",
    "dark": "#804126",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run representative Open-CaBCI_tools analyses on bundled real data."
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--threshold", type=float, default=2.0)
    parser.add_argument("--network-threshold", type=float, default=0.15)
    return parser.parse_args()


def use_chart_theme():
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": TOKENS["surface"],
            "figure.edgecolor": "none",
            "savefig.facecolor": TOKENS["surface"],
            "savefig.edgecolor": "none",
            "axes.facecolor": TOKENS["panel"],
            "axes.edgecolor": TOKENS["axis"],
            "axes.labelcolor": TOKENS["ink"],
            "axes.grid": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": TOKENS["grid"],
            "grid.linewidth": 0.8,
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "sans-serif"],
            "font.monospace": ["DejaVu Sans Mono", "monospace"],
            "patch.linewidth": 1.0,
        },
    )


def add_figure_header(fig, title, subtitle):
    fig.text(
        0.07,
        0.985,
        textwrap.fill(title, 88),
        ha="left",
        va="top",
        fontsize=15,
        fontweight="semibold",
        color=TOKENS["ink"],
    )
    fig.text(
        0.07,
        0.945,
        textwrap.fill(subtitle, 130),
        ha="left",
        va="top",
        fontsize=9,
        color=TOKENS["muted"],
    )


def validate_input(data_dir):
    required = ["F.npy", "Fneu.npy", "spks.npy", "iscell.npy", "ops.npy", "stat.npy"]
    missing = [name for name in required if not (data_dir / name).is_file()]
    if missing:
        raise FileNotFoundError("Missing Suite2p files: " + ", ".join(missing))


def compute_pairwise_correlations(binary_activity):
    cell_count = binary_activity.shape[0]
    correlation = np.eye(cell_count, dtype=np.float32)
    pvalues = np.zeros((cell_count, cell_count), dtype=np.float32)
    rows = []
    for first in range(cell_count):
        for second in range(first + 1, cell_count):
            result, _ = get_corr2(
                binary_activity[first],
                binary_activity[second],
                zscore=False,
            )
            value = float(result[0]) if np.isfinite(result[0]) else 0.0
            pvalue = float(result[1]) if np.isfinite(result[1]) else 1.0
            correlation[first, second] = value
            correlation[second, first] = value
            pvalues[first, second] = pvalue
            pvalues[second, first] = pvalue
            rows.append((first, second, value, pvalue))
            rows.append((second, first, value, pvalue))
    return correlation, pvalues, np.asarray(rows, dtype=np.float64)


def prepare_analysis(data_dir, output_dir, threshold, network_threshold):
    validate_input(data_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    calcium = Calcium()
    calcium.data_dir = str(data_dir)
    calcium.verbose = True
    calcium.load_suite2p()
    calcium.load_footprints()

    ops = np.load(data_dir / "ops.npy", allow_pickle=True).item()
    source_cell_ids = np.load(data_dir / "source_cell_ids.npy")
    sample_rate = float(ops.get("fs", calcium.sample_rate))
    corrected = calcium.F - 0.7 * calcium.Fneu
    _, dff_time_first = compute_dff0(corrected.T)
    dff = np.asarray(dff_time_first.T, dtype=np.float32)

    calcium.data_dir = str(output_dir)
    binary_path = output_dir / "binarized.npy"
    binary_path.unlink(missing_ok=True)
    binary_activity, _ = calcium.binarize(dff, thresh=threshold)
    binary_activity = np.asarray(binary_activity, dtype=np.uint8)
    np.save(binary_path, binary_activity)

    bin_frames = 10
    usable_frames = dff.shape[1] - dff.shape[1] % bin_frames
    scale = np.std(dff[:, :usable_frames], axis=1, keepdims=True)
    standardized = (dff[:, :usable_frames] - np.mean(dff[:, :usable_frames], axis=1, keepdims=True)) / np.maximum(scale, 1e-6)
    population_bins = standardized.reshape(
        standardized.shape[0], -1, bin_frames
    ).mean(axis=2)
    calcium.recompute_PCA = True
    pca, embedding = calcium.compute_PCA(
        population_bins.T,
        suffix1="demo_",
        recompute=True,
        save=True,
    )

    correlation, pvalues, pairs = compute_pairwise_correlations(binary_activity)
    pair_path = output_dir / "correlation_pairs.npy"
    np.save(pair_path, pairs)
    adjacency = load_correlation_matrix(str(pair_path), corr_thresh=network_threshold)
    np.fill_diagonal(adjacency, 0)
    graph = generate_graph_from_connected_nodes(adjacency)
    graph = nx.Graph(graph)
    graph.remove_edges_from(nx.selfloop_edges(graph))
    graph.remove_nodes_from(list(nx.isolates(graph)))
    degree_bins, degree_counts = get_degree_distribution(graph)

    results_path = output_dir / "demo_results.npz"
    np.savez_compressed(
        results_path,
        dff=dff,
        binary_activity=binary_activity,
        pca_embedding=embedding,
        pca_explained_variance_ratio=pca.explained_variance_ratio_,
        correlation_matrix=correlation,
        correlation_pvalues=pvalues,
        adjacency_matrix=adjacency,
        source_cell_ids=source_cell_ids,
        sample_rate_hz=sample_rate,
        pca_bin_frames=bin_frames,
        degree_bins=degree_bins,
        degree_counts=degree_counts,
    )

    with (output_dir / "network_edges.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["cell_a", "cell_b", "source_cell_a", "source_cell_b", "correlation", "pvalue"])
        for first, second in sorted(graph.edges()):
            writer.writerow(
                [
                    first,
                    second,
                    int(source_cell_ids[first]),
                    int(source_cell_ids[second]),
                    float(correlation[first, second]),
                    float(pvalues[first, second]),
                ]
            )

    analysis = {
        "calcium": calcium,
        "ops": ops,
        "source_cell_ids": source_cell_ids,
        "sample_rate": sample_rate,
        "dff": dff,
        "binary": binary_activity,
        "embedding": embedding,
        "explained": pca.explained_variance_ratio_,
        "correlation": correlation,
        "adjacency": adjacency,
        "graph": graph,
    }
    return analysis, results_path


def save_figure(fig, output_dir, stem):
    for suffix in ("png", "svg"):
        output_path = output_dir / f"{stem}.{suffix}"
        fig.savefig(
            output_path,
            dpi=180,
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
        )
        if suffix == "svg":
            lines = output_path.read_text().splitlines()
            output_path.write_text("\n".join(line.rstrip() for line in lines) + "\n")
    plt.close(fig)


def validate_analysis(analysis):
    binary = analysis["binary"]
    if binary.ndim != 2 or binary.shape[0] < 2 or binary.shape[1] < 2:
        raise RuntimeError(f"Unexpected binary activity shape: {binary.shape}")
    if not np.isin(binary, [0, 1]).all():
        raise RuntimeError("Binarized activity contains values outside {0, 1}")
    for name in ("dff", "embedding", "correlation"):
        if not np.isfinite(analysis[name]).all():
            raise RuntimeError(f"{name} contains non-finite values")
    if analysis["embedding"].shape[0] < 20:
        raise RuntimeError("PCA embedding is too short for population-geometry analysis")


def plot_binarization(analysis, output_dir):
    calcium = analysis["calcium"]
    dff = analysis["dff"]
    binary = analysis["binary"]
    sample_rate = analysis["sample_rate"]
    ops = analysis["ops"]
    source_ids = analysis["source_cell_ids"]
    duration = dff.shape[1] / sample_rate

    fig = plt.figure(figsize=(13, 8.2), facecolor=TOKENS["surface"])
    grid = fig.add_gridspec(2, 2, width_ratios=(0.9, 1.45), hspace=0.34, wspace=0.28)
    ax_fov = fig.add_subplot(grid[:, 0])
    ax_trace = fig.add_subplot(grid[0, 1])
    ax_raster = fig.add_subplot(grid[1, 1])

    mean_image = np.asarray(ops["meanImg"])
    low, high = np.percentile(mean_image, [1, 99.7])
    ax_fov.imshow(mean_image, cmap="gray", vmin=low, vmax=high)
    for contour in calcium.contours:
        ax_fov.plot(contour[:, 0], contour[:, 1], color=ORANGE["base"], linewidth=0.7, alpha=0.85)
    ax_fov.set_title("Real field of view and selected ROIs", loc="left", fontsize=10, color=TOKENS["ink"])
    ax_fov.set_xlabel("x (pixels)")
    ax_fov.set_ylabel("y (pixels)")
    ax_fov.grid(False)

    event_counts = binary.sum(axis=1)
    trace_ids = np.argsort(event_counts)[-6:][::-1]
    time = np.arange(dff.shape[1]) / sample_rate
    scale = max(float(np.percentile(np.abs(dff[trace_ids]), 98)), 1e-3)
    for offset, cell_id in enumerate(trace_ids):
        trace = np.clip(dff[cell_id] / scale, -1.2, 2.0) + offset * 2.4
        ax_trace.plot(time, trace, color=BLUE["dark"], linewidth=0.75)
        event_times = time[binary[cell_id].astype(bool)]
        ax_trace.scatter(
            event_times,
            np.full(event_times.shape, offset * 2.4 - 0.55),
            s=5,
            color=ORANGE["base"],
            edgecolor=ORANGE["dark"],
            linewidth=0.25,
        )
    ax_trace.set_yticks(
        np.arange(len(trace_ids)) * 2.4,
        [f"ROI {int(source_ids[cell_id])}" for cell_id in trace_ids],
    )
    ax_trace.set_xlim(0, duration)
    ax_trace.set_xlabel("Time (s)")
    ax_trace.set_title("Fluorescence traces and detected events", loc="left", fontsize=10, color=TOKENS["ink"])
    ax_trace.grid(axis="x", linestyle=":")
    ax_trace.grid(axis="y", visible=False)

    raster_cmap = sns.blend_palette(
        [TOKENS["panel"], BLUE["light"], BLUE["dark"]], as_cmap=True
    )
    sns.heatmap(
        binary,
        ax=ax_raster,
        cmap=raster_cmap,
        cbar=False,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
    )
    tick_seconds = np.linspace(0, duration, 6)
    ax_raster.set_xticks(tick_seconds / duration * binary.shape[1])
    ax_raster.set_xticklabels([f"{value:.0f}" for value in tick_seconds])
    ax_raster.set_xlabel("Time (s)")
    ax_raster.set_ylabel("ROI")
    ax_raster.set_title("Binarized population activity", loc="left", fontsize=10, color=TOKENS["ink"])

    add_figure_header(
        fig,
        "Real calcium traces are converted into sparse population events",
        f"DON-014371 · 20230411Rec calibration · {dff.shape[0]} classified cells · {dff.shape[1]} frames · {duration:.1f} s",
    )
    fig.subplots_adjust(top=0.88, left=0.07, right=0.97, bottom=0.08)
    save_figure(fig, output_dir, "calcium_binarization")


def plot_population_geometry(analysis, output_dir, network_threshold):
    embedding = analysis["embedding"]
    explained = analysis["explained"]
    correlation = analysis["correlation"]
    graph = analysis["graph"]
    calcium = analysis["calcium"]
    sample_rate = analysis["sample_rate"]
    source_ids = analysis["source_cell_ids"]

    fig = plt.figure(figsize=(13, 5.2), facecolor=TOKENS["surface"])
    grid = fig.add_gridspec(1, 3, wspace=0.34)
    ax_pca = fig.add_subplot(grid[0, 0])
    ax_corr = fig.add_subplot(grid[0, 1])
    ax_graph = fig.add_subplot(grid[0, 2])

    time_seconds = np.arange(embedding.shape[0]) * 10 / sample_rate
    time_cmap = sns.blend_palette(
        [BLUE["xlight"], BLUE["light"], BLUE["base"], BLUE["dark"]], as_cmap=True
    )
    ax_pca.plot(embedding[:, 0], embedding[:, 1], color=BLUE["light"], linewidth=0.7, zorder=1)
    points = ax_pca.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=time_seconds,
        cmap=time_cmap,
        s=14,
        edgecolor=BLUE["dark"],
        linewidth=0.2,
        zorder=2,
    )
    colorbar = fig.colorbar(points, ax=ax_pca, fraction=0.046, pad=0.04)
    colorbar.set_label("Time (s)")
    ax_pca.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax_pca.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax_pca.set_title("Population trajectory", loc="left", fontsize=10, color=TOKENS["ink"])
    ax_pca.grid(linestyle=":")

    diverging = sns.blend_palette(
        [ORANGE["dark"], ORANGE["xlight"], TOKENS["panel"], BLUE["xlight"], BLUE["dark"]],
        as_cmap=True,
    )
    sns.heatmap(
        correlation,
        ax=ax_corr,
        cmap=diverging,
        vmin=-0.5,
        vmax=0.5,
        center=0,
        square=True,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={"label": "Pearson r", "shrink": 0.72},
    )
    ax_corr.set_xlabel("ROI")
    ax_corr.set_ylabel("ROI")
    ax_corr.set_title("Pairwise event correlation", loc="left", fontsize=10, color=TOKENS["ink"])

    positions = {}
    for cell_id in graph.nodes:
        stat = calcium.stat[cell_id]
        positions[cell_id] = (float(np.mean(stat["xpix"])), -float(np.mean(stat["ypix"])))
    degrees = dict(graph.degree())
    if graph.number_of_nodes():
        node_colors = [degrees[node] for node in graph.nodes]
        nx.draw_networkx_edges(
            graph,
            positions,
            ax=ax_graph,
            edge_color=BLUE["light"],
            width=0.8,
            alpha=0.65,
        )
        nodes = nx.draw_networkx_nodes(
            graph,
            positions,
            ax=ax_graph,
            node_color=node_colors,
            cmap=time_cmap,
            node_size=55,
            edgecolors=BLUE["dark"],
            linewidths=0.6,
        )
        labels = {
            node: str(int(source_ids[node]))
            for node in sorted(degrees, key=degrees.get, reverse=True)[:5]
        }
        nx.draw_networkx_labels(graph, positions, labels=labels, ax=ax_graph, font_size=6)
        colorbar = fig.colorbar(nodes, ax=ax_graph, fraction=0.046, pad=0.04)
        colorbar.set_label("Degree")
        colorbar.locator = mticker.MaxNLocator(integer=True)
        colorbar.update_ticks()
    else:
        ax_graph.text(0.5, 0.5, "No significant edges", ha="center", va="center", transform=ax_graph.transAxes)
    ax_graph.set_title("Spatial functional network", loc="left", fontsize=10, color=TOKENS["ink"])
    ax_graph.set_xlabel("x (pixels)")
    ax_graph.set_ylabel("y (pixels)")
    ax_graph.grid(False)

    add_figure_header(
        fig,
        "Low-dimensional population geometry and significant co-activity structure",
        f"PCA uses 10-frame bins; network edges require Pearson r ≥ {network_threshold:.2f} and p < 0.10; node positions are ROI centroids",
    )
    fig.subplots_adjust(top=0.83, left=0.07, right=0.98, bottom=0.13)
    save_figure(fig, output_dir, "population_geometry")


def write_summary(analysis, output_dir, threshold, network_threshold):
    binary = analysis["binary"]
    explained = analysis["explained"]
    correlation = analysis["correlation"]
    graph = analysis["graph"]
    upper = correlation[np.triu_indices_from(correlation, k=1)]
    summary = {
        "cells": int(binary.shape[0]),
        "frames": int(binary.shape[1]),
        "sample_rate_hz": float(analysis["sample_rate"]),
        "duration_seconds": float(binary.shape[1] / analysis["sample_rate"]),
        "binarization_threshold_std": float(threshold),
        "active_samples": int(binary.sum()),
        "active_fraction": float(binary.mean()),
        "median_events_per_cell": float(np.median(binary.sum(axis=1))),
        "pca_pc1_variance_fraction": float(explained[0]),
        "pca_pc2_variance_fraction": float(explained[1]),
        "median_pairwise_correlation": float(np.median(upper)),
        "network_threshold": float(network_threshold),
        "network_nodes": int(graph.number_of_nodes()),
        "network_edges": int(graph.number_of_edges()),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    args = parse_args()
    data_dir = args.data.resolve()
    output_dir = args.output.resolve()
    use_chart_theme()
    analysis, results_path = prepare_analysis(
        data_dir,
        output_dir,
        args.threshold,
        args.network_threshold,
    )
    validate_analysis(analysis)
    plot_binarization(analysis, output_dir)
    plot_population_geometry(analysis, output_dir, args.network_threshold)
    summary = write_summary(
        analysis,
        output_dir,
        args.threshold,
        args.network_threshold,
    )
    print(json.dumps(summary, indent=2))
    print(f"Saved numerical results: {results_path}")
    print(f"Saved figures: {output_dir / 'calcium_binarization.png'}")
    print(f"Saved figures: {output_dir / 'population_geometry.png'}")


if __name__ == "__main__":
    main()
