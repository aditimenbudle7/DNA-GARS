#!/usr/bin/env python3
"""
Plot DNA-GARS benchmark results.

Usage:
    python3 scripts/plot.py [results_dir]     (default: results)

Reads search.csv, align.csv and reconstruct.csv written by the benchmark menu
and saves PNG graphs next to them. Requires matplotlib (pip install matplotlib).
"""
import csv
import os
import sys
from collections import defaultdict

try:
    import matplotlib

    matplotlib.use("Agg")  # no window needed, just write files
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("matplotlib is required: pip install matplotlib")

MARKERS = ["o", "s", "^", "D", "v"]
STYLES = ["-", "--", "-.", ":", "-"]


def read_csv(path):
    if not os.path.exists(path):
        print(f"skip: {path} not found")
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def save(fig, out_dir, name):
    path = os.path.join(out_dir, name)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"saved {path}")


def line_chart(series, title, xlabel, ylabel, log_x=True, log_y=True):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, (name, points) in enumerate(series.items()):
        # log axes cannot show 0, so drop non-positive points
        points = sorted(p for p in points if p[1] > 0 or not log_y)
        if not points:
            continue
        ax.plot([p[0] for p in points], [p[1] for p in points],
                marker=MARKERS[i % len(MARKERS)], linestyle=STYLES[i % len(STYLES)],
                linewidth=1.5, markersize=5, label=name)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    return fig


def plot_search(rows, out_dir):
    scenarios = sorted({r["scenario"] for r in rows})
    for scenario in scenarios:
        for metric, label in (("time_ms", "Time (ms)"), ("comparisons", "Comparisons")):
            series = defaultdict(list)
            for r in rows:
                if r["scenario"] != scenario:
                    continue
                value = float(r[metric])
                name = r["algorithm"]
                if metric == "comparisons" and name == "Rabin-Karp":
                    # Rabin-Karp's work is mostly hash comparisons, so count both
                    value += float(r.get("hash_checks", 0) or 0)
                    name = "Rabin-Karp (hash + char)"
                series[name].append((int(r["n"]), value))
            fig = line_chart(series, f"Pattern search ({scenario} genome): {label.lower()} vs n",
                             "Genome length n", label)
            save(fig, out_dir, f"search_{scenario}_{metric}.png")


def plot_align(rows, out_dir):
    series = defaultdict(list)
    for r in rows:
        series[r["algorithm"]].append((int(r["cells"]), float(r["time_ms"])))
    fig = line_chart(series, "Alignment: time vs DP cells (n x m)", "DP cells (n x m)", "Time (ms)")
    save(fig, out_dir, "align_time.png")


def plot_reconstruct(rows, out_dir):
    series = defaultdict(list)
    lengths = defaultdict(dict)
    for r in rows:
        k = int(r["k"])
        series[r["algorithm"]].append((k, float(r["time_ms"])))
        lengths[k][r["algorithm"]] = int(r["superstring_len"])
    fig = line_chart(series, "Reconstruction: time vs number of fragments k",
                     "Fragments k", "Time (ms)", log_x=False)
    save(fig, out_dir, "reconstruct_time.png")

    ratio = [(k, v["Greedy"] / v["DP-Bitmask"]) for k, v in sorted(lengths.items())
             if "Greedy" in v and "DP-Bitmask" in v]
    if ratio:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.bar([k for k, _ in ratio], [x for _, x in ratio])
        ax.axhline(1.0, color="black", linewidth=1)
        ax.set_ylim(bottom=0.9)
        ax.set_title("Greedy approximation: superstring length / optimal length")
        ax.set_xlabel("Fragments k")
        ax.set_ylabel("Greedy / optimal")
        ax.grid(True, axis="y", alpha=0.3)
        save(fig, out_dir, "reconstruct_greedy_ratio.png")


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "results"
    if not os.path.isdir(out_dir):
        sys.exit(f"folder not found: {out_dir}")

    search = read_csv(os.path.join(out_dir, "search.csv"))
    align = read_csv(os.path.join(out_dir, "align.csv"))
    recon = read_csv(os.path.join(out_dir, "reconstruct.csv"))

    if search:
        plot_search(search, out_dir)
    if align:
        plot_align(align, out_dir)
    if recon:
        plot_reconstruct(recon, out_dir)
    if not (search or align or recon):
        print("No CSV files found. Run the benchmarks from the DNA-GARS menu first.")


if __name__ == "__main__":
    main()
