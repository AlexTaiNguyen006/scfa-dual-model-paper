#!/usr/bin/env python3
"""Generate all manuscript figures (Fig 1–10) from results CSVs."""

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .utils import build_paths

# plot colors
COLORS = {
    "acetate":    "#1b9e77",
    "propionate": "#d95f02",
    "butyrate":   "#7570b3",
    "glucose":    "#e7298a",
    "oxygen":     "#66a61e",
    "co2":        "#a65628",
    "objective":  "#e6ab02",
    "baseline":   "#999999",
    "delta":      "#d62728",
}

MODEL_COLORS = {
    "Recon3D":   "#3a86ff",
    "Human-GEM": "#ef476f",
}

COND_ORDER = ["StachysDose_High", "StachysDose_Mid", "StachysDose_Low"]
XLABELS = ["High", "Mid", "Low"]

FIG_W, FIG_H = 7, 4.5
DPI = 300

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
})


def _save(fig, paths, name):
    """Save figure as both PNG and PDF."""
    fig.savefig(paths.figs_dir / f"{name}.png")
    fig.savefig(paths.figs_dir / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}.png / .pdf")


def main():
    paths = build_paths()
    df = pd.read_csv(paths.results / "merged_dose_scfa_host.csv")

    # sort
    order_map = {c: i for i, c in enumerate(COND_ORDER)}
    df["_ord"] = df["condition"].map(order_map)
    df = df.sort_values("_ord").reset_index(drop=True)
    x = np.arange(len(df))

    # Fig 1: SCFA availability by dose
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for scfa, col_name in [("acetate", "acetate_mmol_gDW_hr"),
                            ("propionate", "propionate_mmol_gDW_hr"),
                            ("butyrate", "butyrate_mmol_gDW_hr")]:
        ax.plot(x, df[col_name], "o-", color=COLORS[scfa],
                label=scfa.capitalize(), lw=2, markersize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_xlabel("Stachyose Dose")
    ax.set_ylabel("SCFA Availability (mmol/gDW/hr)")
    ax.set_title("SCFA Availability by Dose Condition")
    ax.legend()
    _save(fig, paths, "Fig1")

    # Fig 2: SCFA molar ratios
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    total = (df["acetate_mmol_gDW_hr"] + df["propionate_mmol_gDW_hr"]
             + df["butyrate_mmol_gDW_hr"])
    frac_ac  = df["acetate_mmol_gDW_hr"] / total
    frac_ppa = df["propionate_mmol_gDW_hr"] / total
    frac_but = df["butyrate_mmol_gDW_hr"] / total

    bar_w = 0.55
    ax.bar(x, frac_ac, bar_w, color=COLORS["acetate"], label="Acetate")
    ax.bar(x, frac_ppa, bar_w, bottom=frac_ac,
           color=COLORS["propionate"], label="Propionate")
    ax.bar(x, frac_but, bar_w, bottom=frac_ac + frac_ppa,
           color=COLORS["butyrate"], label="Butyrate")
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_ylabel("Molar Fraction")
    ax.set_title("SCFA Molar Ratio by Condition")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper right")
    _save(fig, paths, "Fig2")

    # Fig 3: Dual-model ATPM comparison
    mc = pd.read_csv(paths.results / "model_comparison.csv")
    cond_order_mc = ["StachysDose_High", "StachysDose_Mid", "StachysDose_Low"]
    labels_mc = ["High", "Mid", "Low"]
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    w = 0.35
    for i, mdl in enumerate(["Recon3D", "Human-GEM"]):
        sub = mc[mc["model"] == mdl].copy()
        sub["_o"] = sub["condition"].map({c: j for j, c in enumerate(cond_order_mc)})
        sub = sub.sort_values("_o")
        xpos = np.arange(len(sub)) + (i - 0.5) * w
        ax.bar(xpos, sub["objective_value"], w, label=mdl,
               color=MODEL_COLORS[mdl], edgecolor="black", linewidth=0.5)
    ax.set_xticks(np.arange(len(labels_mc)))
    ax.set_xticklabels(labels_mc)
    ax.set_xlabel("Stachyose Dose")
    ax.set_ylabel("ATPM (mmol/gDW/hr)")
    ax.set_title("Dual-Model ATPM Comparison")
    ax.legend()
    _save(fig, paths, "Fig3")

    # Fig 4: Host exchange fluxes
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    flux_series = [
        ("Glucose",    "glucose_flux",    COLORS["glucose"]),
        ("Acetate",    "acetate_flux",    COLORS["acetate"]),
        ("Propionate", "propionate_flux", COLORS["propionate"]),
        ("Butyrate",   "butyrate_flux",   COLORS["butyrate"]),
        ("O\u2082",    "oxygen_flux",     COLORS["oxygen"]),
        ("CO\u2082",   "co2_flux",        COLORS["co2"]),
    ]
    for label, col, color in flux_series:
        if col in df.columns:
            ax.plot(x, df[col], "o-", color=color, label=label,
                    lw=2, markersize=6)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_ylabel("Flux (mmol/gDW/hr)")
    ax.set_title("Exchange Fluxes by Condition")
    ax.legend(ncol=2, loc="lower left", fontsize=9)
    _save(fig, paths, "Fig4")

    # Fig 5: Pathway flux heatmap
    pw_cols = [c for c in df.columns if c.startswith("pathway_")]
    if pw_cols:
        pw = df[pw_cols].copy()
        pw.index = XLABELS
        pw.columns = [c.replace("pathway_", "") for c in pw_cols]
        pw = pw.loc[:, (pw != 0).any() & pw.notna().any()]

        if not pw.empty:
            ncols = pw.shape[1]
            fig_width = max(7, ncols * 0.9)
            fig, ax = plt.subplots(figsize=(fig_width, 4.5))
            im = ax.imshow(pw.values.astype(float), aspect="auto",
                           cmap="RdYlBu_r")
            ax.set_xticks(np.arange(ncols))
            ax.set_xticklabels(pw.columns, rotation=45, ha="right")
            ax.set_yticks(np.arange(len(pw)))
            ax.set_yticklabels(pw.index)
            ax.set_title("Pathway Fluxes by Condition")
            vmax = pw.values.max()
            for i in range(len(pw)):
                for j in range(ncols):
                    v = pw.iloc[i, j]
                    if pd.notna(v) and v != 0:
                        txt_color = "white" if v > vmax * 0.6 else "black"
                        ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                                fontsize=9, color=txt_color)
            plt.colorbar(im, ax=ax, shrink=0.8)
            fig.tight_layout()
            _save(fig, paths, "Fig5")

    # Fig 6: Sensitivity analysis
    sa_path = paths.results / "sensitivity_analysis.csv"
    if sa_path.exists():
        sa = pd.read_csv(sa_path)
        sweeps = sa["sweep"].unique()
        n_sweeps = len(sweeps)
        fig, axes = plt.subplots(1, n_sweeps, figsize=(4 * n_sweeps, 4.5),
                                  squeeze=False)
        for idx, sweep_name in enumerate(sweeps):
            ax = axes[0, idx]
            for mdl, color in MODEL_COLORS.items():
                sub = sa[(sa["sweep"] == sweep_name) & (sa["model"] == mdl)]
                if not sub.empty:
                    ax.plot(sub["sweep_value"], sub["ATPM"], "o-",
                            color=color, label=mdl, lw=2, markersize=5)
            ax.set_xlabel(f"{sweep_name.capitalize()} (mmol/gDW/hr)")
            ax.set_ylabel("ATPM (mmol/gDW/hr)")
            ax.set_title(f"{sweep_name.capitalize()} sweep")
            ax.legend(fontsize=8)
        fig.suptitle("One-at-a-Time Sensitivity Analysis", fontsize=13,
                     fontweight="bold")
        fig.tight_layout()
        _save(fig, paths, "Fig6")

    # Fig 7: Propionate rescue analysis
    rescue_path = paths.results / "table_rescue_constrained.csv"
    if rescue_path.exists():
        rc = pd.read_csv(rescue_path)
        conds = rc["Condition"].tolist()
        x_rc = np.arange(len(conds))
        w = 0.2
        COL_ORIG = "#3a86ff"
        COL_RESC = "#06d6a0"
        COL_HGEM = "#ef476f"

        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 5))

        # Panel A: Propionate flux
        dose_mask = rc["Condition"] != "Baseline"
        rc_dose = rc[dose_mask]
        x_d = np.arange(len(rc_dose))
        ax_a.bar(x_d - w, rc_dose["Recon3D_orig_PPA"].abs(), w,
                 color=COL_ORIG, label="Recon3D (original)",
                 edgecolor="black", linewidth=0.5)
        constrained_ppa_col = "Rescued_constrained_PPA"
        if constrained_ppa_col in rc_dose.columns:
            ax_a.bar(x_d, rc_dose[constrained_ppa_col].abs(), w,
                     color=COL_RESC, label="Recon3D + PPCOACm",
                     edgecolor="black", linewidth=0.5)
        ax_a.bar(x_d + w, rc_dose["Human_GEM_PPA"].abs(), w,
                 color=COL_HGEM, label="Human-GEM",
                 edgecolor="black", linewidth=0.5)
        ax_a.set_xticks(x_d)
        ax_a.set_xticklabels(rc_dose["Condition"].tolist())
        ax_a.set_xlabel("Dose Condition")
        ax_a.set_ylabel("Propionate uptake (mmol gDW\u207b\u00b9 hr\u207b\u00b9)")
        ax_a.set_title("A   Propionate exchange flux", loc="left",
                        fontweight="bold")
        ax_a.legend(fontsize=8)

        # Panel B: ATPM
        ax_b.bar(x_rc - w, rc["Recon3D_orig_ATPM"], w,
                 color=COL_ORIG, label="Recon3D (original)",
                 edgecolor="black", linewidth=0.5)
        ax_b.bar(x_rc, rc["Rescued_constrained_ATPM"], w,
                 color=COL_RESC, label="Recon3D + PPCOACm (constrained)",
                 edgecolor="black", linewidth=0.5)
        ax_b.bar(x_rc + w, rc["Human_GEM_ATPM"], w,
                 color=COL_HGEM, label="Human-GEM",
                 edgecolor="black", linewidth=0.5)
        ax_b.set_xticks(x_rc)
        ax_b.set_xticklabels(conds)
        ax_b.set_xlabel("Dose Condition")
        ax_b.set_ylabel("Predicted ATPM (mmol gDW\u207b\u00b9 hr\u207b\u00b9)")
        ax_b.set_title("B   Predicted ATPM flux", loc="left",
                        fontweight="bold")
        ax_b.legend(fontsize=8)
        fig.tight_layout()
        _save(fig, paths, "Fig7")

    # Fig 8: SCFA ratio sensitivity
    rs_path = paths.results / "ratio_sensitivity.csv"
    if rs_path.exists():
        rs = pd.read_csv(rs_path)
        fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
        labels = []
        for mdl, color in MODEL_COLORS.items():
            sub = rs[rs["model"] == mdl]
            if not sub.empty:
                ax.barh(np.arange(len(sub)) + (0.2 if mdl == "Human-GEM" else -0.2),
                        sub["ATPM"], 0.35, color=color, label=mdl,
                        edgecolor="black", linewidth=0.5)
                if mdl == "Recon3D":
                    labels = sub["ratio_profile"].tolist()
        if labels:
            ax.set_yticks(np.arange(len(labels)))
            ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("ATPM (mmol/gDW/hr)")
        ax.set_title("SCFA Ratio Sensitivity")
        ax.legend()
        fig.tight_layout()
        _save(fig, paths, "Fig8")

    # Fig 9: pFBA comparison
    pfba_path = paths.results / "pfba_comparison.csv"
    if pfba_path.exists():
        pf = pd.read_csv(pfba_path)
        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 5))
        cond_labels = {"StachysDose_Low": "Low", "StachysDose_Mid": "Mid",
                       "StachysDose_High": "High"}

        # Panel A: ATPM — FBA vs pFBA
        for mdl, color in MODEL_COLORS.items():
            sub = pf[pf["model"] == mdl].copy()
            sub["label"] = sub["condition"].map(cond_labels)
            sub = sub.sort_values("label")
            xp = np.arange(len(sub))
            ax_a.plot(xp, sub["ATPM_FBA"], "o-", color=color, label=f"{mdl} FBA", lw=2)
            ax_a.plot(xp, sub["ATPM_pFBA"], "s--", color=color,
                      label=f"{mdl} pFBA", lw=2, alpha=0.7)
            if mdl == "Recon3D":
                ax_a.set_xticks(xp)
                ax_a.set_xticklabels(sub["label"].tolist())
        ax_a.set_ylabel("ATPM (mmol/gDW/hr)")
        ax_a.set_title("A   FBA vs pFBA objective", loc="left",
                        fontweight="bold")
        ax_a.legend(fontsize=8)

        # Panel B: Total flux reduction
        for mdl, color in MODEL_COLORS.items():
            sub = pf[pf["model"] == mdl].copy()
            sub["label"] = sub["condition"].map(cond_labels)
            sub = sub.sort_values("label")
            xp = np.arange(len(sub))
            ax_b.bar(xp + (0.2 if mdl == "Human-GEM" else -0.2),
                     sub["flux_reduction_pct"], 0.35, color=color, label=mdl,
                     edgecolor="black", linewidth=0.5)
        ax_b.set_xticks(np.arange(len(cond_labels)))
        ax_b.set_xticklabels(["High", "Low", "Mid"])
        ax_b.set_ylabel("Total Flux Reduction (%)")
        ax_b.set_title("B   Flux reduction by pFBA", loc="left",
                        fontweight="bold")
        ax_b.legend(fontsize=8)
        fig.tight_layout()
        _save(fig, paths, "Fig9")

    # Fig 10: Multi-threshold FVA
    fva_path = paths.results / "fva_multi_threshold.csv"
    if fva_path.exists():
        fva = pd.read_csv(fva_path)
        atpm_fva = fva[fva["reaction"] == "ATPM"].copy()
        fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
        cond_labels = {"StachysDose_Low": "Low", "StachysDose_Mid": "Mid",
                       "StachysDose_High": "High"}
        for mdl, color in MODEL_COLORS.items():
            for cond, marker in zip(["StachysDose_Low", "StachysDose_Mid",
                                      "StachysDose_High"],
                                     ["o", "s", "D"]):
                sub = atpm_fva[(atpm_fva["model"] == mdl) &
                               (atpm_fva["condition"] == cond)]
                if sub.empty:
                    continue
                sub = sub.sort_values("fraction_optimum")
                label = f"{mdl} {cond_labels[cond]}"
                ax.fill_between(sub["fraction_optimum"],
                                sub["fva_minimum"], sub["fva_maximum"],
                                alpha=0.15, color=color)
                ax.plot(sub["fraction_optimum"], sub["fva_minimum"],
                        marker + "-", color=color, lw=1.5, markersize=5,
                        label=label)
                ax.plot(sub["fraction_optimum"], sub["fva_maximum"],
                        marker + "--", color=color, lw=1.5, markersize=5,
                        alpha=0.6)
        ax.set_xlabel("Fraction of Optimum")
        ax.set_ylabel("ATPM Range (mmol/gDW/hr)")
        ax.set_title("FVA Solution Space Across Thresholds")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        _save(fig, paths, "Fig10")

    print("done — all figures generated")


if __name__ == "__main__":
    main()

