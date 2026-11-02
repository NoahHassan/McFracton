"""Combine the *_bins.txt files of independent chains (one per seed) into observables with error bars.

    python scripts/analyse_bins.py results/spiderweb_corrected_job=123_seed=*_bins.txt
    python scripts/analyse_bins.py <files> --out table.tsv --plot figure.png

A bins file is written by mcf_run when bin_size > 0: per temperature, one line per bin of consecutive
measurement sweeps, holding the bin means of every observable X and of X^2 (measured after every sweep).

Per chain and temperature the bins are merged into --blocks consecutive blocks of equal length, and
every quantity gets a jackknife error over those blocks. That error is only right if a block is much
longer than the autocorrelation time, which is why the script also prints tau_int of every observable
from the same blocks (they differ: the helicity moduli can be much slower than the energy). The
quantities:

    E/N         mean energy per site
    C/N         specific heat per site, (<E^2> - <E>^2) / (T^2 N)
    HelicityB, HelicityExx, HelicityExy, InstantonDensity
                means of the per-configuration estimators the code measures

Across chains the result is the mean of the per-chain estimates, and its error (err) is the standard
error of that mean over chains: independent chains need no assumption about autocorrelation, so this is
the error bar to quote. err_jk is the within-chain jackknife error combined over chains,
sqrt(sum err_i^2) / n; the two agreeing is the check that the blocks are long enough.

N is linear_size^2 * temporal_size, read from the file header.
"""
import argparse
import math
import sys
from collections import defaultdict

import numpy as np


def read_bins(path):
    """Header values and {temperature: rows} of one bins file; rows as an array per column name."""
    header = {}
    with open(path) as f:
        lines = f.read().splitlines()
    data_lines = []
    names = None
    for line in lines:
        if line.startswith("#"):
            parts = line[1:].split("\t")
            parts = [p.strip() for p in parts]
            # "# key\tvalue\tkey\tvalue ..." lines; prose comment lines have no tab and are skipped
            if len(parts) >= 2 and len(parts) % 2 == 0:
                for k, v in zip(parts[0::2], parts[1::2]):
                    header[k] = v
            continue
        if names is None:
            names = line.split("\t")
            continue
        if line.strip():
            data_lines.append(line.split("\t"))
    if names is None:
        raise ValueError(f"{path}: no column header")
    values = np.array([[float(x) for x in row] for row in data_lines]) if data_lines else np.zeros((0, len(names)))
    columns = {name: values[:, i] for i, name in enumerate(names)}
    by_t = defaultdict(list)
    for i, t in enumerate(columns["T"]):
        by_t[t].append(i)
    per_t = {t: {name: col[idx] for name, col in columns.items()} for t, idx in by_t.items()}
    return header, names, per_t


def blocks_of(rows, n_blocks, discard):
    """Merge bins into n_blocks blocks of (almost) equal sweep count; returns (weights, {name: block means})."""
    n = rows["n_sweeps"]
    first = int(round(discard * len(n)))
    keep = slice(first, None)
    n = n[keep]
    if len(n) < 2:
        return None, None
    n_blocks = min(n_blocks, len(n))
    # Split the bin sequence into contiguous groups with about the same number of bins each.
    edges = np.linspace(0, len(n), n_blocks + 1).round().astype(int)
    weights = np.array([n[a:b].sum() for a, b in zip(edges[:-1], edges[1:])])
    means = {}
    for name, col in rows.items():
        if name in ("T", "bin", "n_sweeps"):
            continue
        col = col[keep]
        means[name] = np.array([(col[a:b] * n[a:b]).sum() for a, b in zip(edges[:-1], edges[1:])]) / weights
    return weights, means


def estimates(weights, means, temperature, n_sites, observables):
    """Every derived quantity from (weighted) block means."""
    w = weights / weights.sum()
    avg = {name: float((m * w).sum()) for name, m in means.items()}
    out = {}
    e, e2 = avg["Energy"], avg["Energy_sq"]
    out["E/N"] = e / n_sites
    out["C/N"] = (e2 - e * e) / (temperature ** 2 * n_sites)
    for name in observables:
        if name != "Energy":
            out[name] = avg[name]
    return out


def jackknife(weights, means, temperature, n_sites, observables):
    """Estimates and jackknife errors from the blocks, leaving one block out at a time."""
    full = estimates(weights, means, temperature, n_sites, observables)
    k = len(weights)
    leave_out = []
    for i in range(k):
        mask = np.arange(k) != i
        leave_out.append(estimates(weights[mask], {n: m[mask] for n, m in means.items()},
                                   temperature, n_sites, observables))
    errors = {}
    for q in full:
        vals = np.array([lo[q] for lo in leave_out])
        errors[q] = math.sqrt((k - 1) / k * ((vals - vals.mean()) ** 2).sum())
    return full, errors


def tau_int(weights, means, name):
    """tau_int of one observable, in sweeps, from the blocks: block length * var(block means) / (2 var(one sweep)).

    Only an estimate, and one that saturates at about a block length: when it comes out comparable to
    the block length, the blocks are too short for that observable and its jackknife error is too small.
    """
    w = weights / weights.sum()
    x = float((means[name] * w).sum())
    var_single = float((means[name + "_sq"] * w).sum()) - x * x
    if var_single <= 0 or len(weights) < 2:
        return float("nan")
    var_blocks = float(((means[name] - x) ** 2 * w).sum()) * len(weights) / (len(weights) - 1)
    return weights.mean() * var_blocks / (2.0 * var_single)


def half_difference(rows, discard, name):
    """Mean of the second half of the kept bins minus the mean of the first half (sweep-weighted)."""
    n = rows["n_sweeps"]
    first = int(round(discard * len(n)))
    n, x = n[first:], rows[name][first:]
    if len(n) < 2:
        return float("nan")
    half = len(n) // 2
    return float((x[half:] * n[half:]).sum() / n[half:].sum() - (x[:half] * n[:half]).sum() / n[:half].sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", help="bins files, one per chain")
    parser.add_argument("--blocks", type=int, default=20, help="blocks per chain and temperature (default 20)")
    parser.add_argument("--discard", type=float, default=0.0,
                        help="fraction of the bins at each temperature dropped as extra thermalization (default 0)")
    parser.add_argument("--out", help="write the table here instead of to stdout")
    parser.add_argument("--plot", help="write a figure (png/pdf) of every quantity against T")
    args = parser.parse_args()

    chains = []
    n_sites = None
    observables = None
    for path in args.files:
        header, names, per_t = read_bins(path)
        sites = int(float(header["linear_size"])) ** 2 * int(float(header["temporal_size"]))
        obs = [n for n in names[3:] if not n.endswith("_sq")]
        if n_sites is None:
            n_sites, observables = sites, obs
        elif sites != n_sites or obs != observables:
            sys.exit(f"{path}: lattice or observables differ from {args.files[0]}")
        chains.append((path, per_t))

    quantities = ["E/N", "C/N"] + [o for o in observables if o != "Energy"]
    temperatures = sorted({t for _, per_t in chains for t in per_t})

    rows = []
    for t in temperatures:
        per_chain, per_chain_err, taus, sweeps, block_lengths, drifts = [], [], [], [], [], []
        for path, per_t in chains:
            if t not in per_t:
                continue
            weights, means = blocks_of(per_t[t], args.blocks, args.discard)
            if weights is None:
                continue
            drifts.append(half_difference(per_t[t], args.discard, "Energy") / n_sites)
            est, err = jackknife(weights, means, t, n_sites, observables)
            per_chain.append(est)
            per_chain_err.append(err)
            taus.append({o: tau_int(weights, means, o) for o in observables})
            sweeps.append(weights.sum())
            block_lengths.append(weights.mean())
        n = len(per_chain)
        if n == 0:
            continue
        row = {"T": t, "n_chains": n, "sweeps_per_chain": float(np.mean(sweeps)),
               "block_length": float(np.mean(block_lengths))}
        for o in observables:
            row["tau " + o] = float(np.nanmean([tau[o] for tau in taus]))
        # Equilibration check: in equilibrium the two halves of a chain differ only by noise, which
        # averages out over chains. A drift shared by the chains (z of several sigma) means the
        # temperature was not thermalized and its numbers, the specific heat above all, are not
        # equilibrium values.
        d = np.array(drifts)
        row["drift E/N"] = float(d.mean())
        row["drift_z"] = float(d.mean() / (d.std(ddof=1) / math.sqrt(n))) if n > 1 and d.std() > 0 else float("nan")
        for q in quantities:
            vals = np.array([c[q] for c in per_chain])
            row[q] = float(vals.mean())
            row["err " + q] = float(vals.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
            row["err_jk " + q] = float(math.sqrt(sum(e[q] ** 2 for e in per_chain_err)) / n)
        rows.append(row)

    columns = ["T", "n_chains", "sweeps_per_chain", "block_length", "drift E/N", "drift_z"]
    for q in quantities:
        columns += [q, "err " + q, "err_jk " + q]
    columns += ["tau " + o for o in observables]
    for row in rows:
        if abs(row["drift_z"]) > 3:
            print(f"warning: T = {row['T']:.6g} drifts (E/N second half - first half = {row['drift E/N']:.3g}, "
                  f"z = {row['drift_z']:.1f}): not thermalized, raise the thermalization or use --discard",
                  file=sys.stderr)
    out = open(args.out, "w") if args.out else sys.stdout
    out.write(f"# {len(chains)} chains, N = {n_sites} sites, {args.blocks} blocks per chain and temperature, "
              f"discard = {args.discard}\n")
    out.write("# err: standard error over chains (quote this); err_jk: within-chain jackknife, combined\n")
    out.write("# tau <X>: integrated autocorrelation time of X in sweeps, mean over chains; it cannot come out\n"
              "# much above block_length, and when it gets close to it the blocks are too short for X\n"
              "# drift E/N: second half minus first half of each chain, mean over chains; drift_z: that in units\n"
              "# of its standard error. |drift_z| > 3 means the temperature was not thermalized\n")
    out.write("\t".join(c.replace(" ", "_") for c in columns) + "\n")
    for row in rows:
        out.write("\t".join(f"{row[c]:.10g}" for c in columns) + "\n")
    if args.out:
        out.close()

    if args.plot:
        plot(rows, quantities, args.plot, len(chains), n_sites)


def plot(rows, quantities, path, n_chains, n_sites):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
    plt.rcParams.update({
        "font.size": 9, "figure.dpi": 150, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.edgecolor": AXIS,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 9, "axes.titleweight": "bold",
    })
    t = np.array([r["T"] for r in rows])
    n_cols = 3
    n_rows = math.ceil(len(quantities) / n_cols)
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(11, 3.1 * n_rows), sharex=True, squeeze=False)
    for i, q in enumerate(quantities):
        a = ax[i // n_cols, i % n_cols]
        y = np.array([r[q] for r in rows])
        e = np.array([r["err " + q] for r in rows])
        a.errorbar(t, y, yerr=e, color="#2a78d6", marker="o", markersize=3.5, linewidth=1.4,
                   elinewidth=0.9, capsize=2, markeredgecolor=SURFACE, markeredgewidth=0.6)
        a.set_title(q)
    for j in range(len(quantities), n_rows * n_cols):
        ax[j // n_cols, j % n_cols].axis("off")
    for a in ax[-1]:
        a.set_xlabel(r"temperature $T$")
    fig.suptitle(f"{n_chains} chains, N = {n_sites} sites; error bars: standard error over chains",
                 fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path)


if __name__ == "__main__":
    main()
