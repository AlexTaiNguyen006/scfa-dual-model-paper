#!/usr/bin/env python3


import pandas as pd
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .utils import build_paths, load_config

#colorblind-friendly palette
PAL = {
    "pipeline": "#0072b2",
    "agora2":   "#e69f00",
    "micom":    "#009e73",
    "lit":      "#cc79a7",
}

SOURCE_COLORS = {
    "AGORA2_western":    PAL["agora2"],
    "AGORA2_highfiber":  PAL["agora2"],
    "AGORA2_lowfiber":   PAL["agora2"],
    "Baldini2019_avg":   PAL["agora2"],
    "MICOM_starchtype":  PAL["micom"],
    "MICOM_fos":         PAL["micom"],
    "literature_invitro": PAL["lit"],
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
})


def load_community_scfa(path):
    
    df = pd.read_csv(path, comment="#")
    #normalize column names
    df.columns = [c.strip().lower() for c in df.columns]
    for col in ["acetate", "propionate", "butyrate"]:
        if col not in df.columns:
            raise ValueError(f"Community SCFA file missing column: {col}")
    return df


def load_pipeline_conditions(results_dir):
    
    path = results_dir / "scfa_inputs_canonical.csv"
    df = pd.read_csv(path)
    return df


def cross_validate(community_df, pipeline_df):
    
    pipe_vals = pipeline_df[["acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr",
                              "butyrate_mmol_gDW_hr"]].values
    pipe_conds = pipeline_df["condition"].values

    rows = []
    for _, crow in community_df.iterrows():
        cvec = np.array([crow["acetate"], crow["propionate"], crow["butyrate"]])
        dists = np.linalg.norm(pipe_vals - cvec, axis=1)
        nearest_idx = int(np.argmin(dists))
        nearest_cond = pipe_conds[nearest_idx]
        nearest_dist = dists[nearest_idx]

        #check if within range
        total_community = crow["acetate"] + crow["propionate"] + crow["butyrate"]
        pipe_total = pipe_vals[nearest_idx].sum()

        rows.append({
            "source": crow.get("source", "unknown"),
            "diet_context": crow.get("diet_context", ""),
            "community_acetate": crow["acetate"],
            "community_propionate": crow["propionate"],
            "community_butyrate": crow["butyrate"],
            "community_total": total_community,
            "nearest_condition": nearest_cond.replace("StachysDose_", ""),
            "pipeline_total": pipe_total,
            "euclidean_distance": round(nearest_dist, 3),
            "within_pipeline_range": "Yes" if nearest_dist < pipe_total * 0.5 else "Approximate",
        })

    return pd.DataFrame(rows)


def fig_cross_validation(community_df, pipeline_df, figs_dir):
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    #panel a: total scfa comparison
    #pipeline conditions as horizontal bands
    cond_labels = {"StachysDose_Low": "Low", "StachysDose_Mid": "Mid",
                   "StachysDose_High": "High"}
    for _, row in pipeline_df.iterrows():
        total = (row["acetate_mmol_gDW_hr"] + row["propionate_mmol_gDW_hr"]
                 + row["butyrate_mmol_gDW_hr"])
        label = cond_labels[row["condition"]]
        ax1.axhline(total, color=PAL["pipeline"], ls="--", lw=1.2, alpha=0.6)
        ax1.text(len(community_df) - 0.3, total + 0.3, f"Pipeline: {label}",
                 fontsize=8, color=PAL["pipeline"], ha="right")

    #community points
    for i, (_, crow) in enumerate(community_df.iterrows()):
        total = crow["acetate"] + crow["propionate"] + crow["butyrate"]
        color = SOURCE_COLORS.get(crow.get("source", ""), "#666")
        ax1.bar(i, total, color=color, alpha=0.8, width=0.7)

    ax1.set_xticks(range(len(community_df)))
    labels = community_df["source"].apply(lambda s: s.replace("_", "\n")).tolist()
    ax1.set_xticklabels(labels, fontsize=7, rotation=0)
    ax1.set_ylabel("Total SCFA (mmol gDW$^{-1}$ hr$^{-1}$)")
    ax1.set_title("A) Community SCFA Output vs Pipeline Conditions")

    #legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PAL["agora2"], label="AGORA2"),
        Patch(facecolor=PAL["micom"], label="MICOM"),
        Patch(facecolor=PAL["lit"], label="Literature"),
        plt.Line2D([0], [0], color=PAL["pipeline"], ls="--", label="Pipeline doses"),
    ]
    ax1.legend(handles=legend_elements, fontsize=8)

    #panel b: scfa ratio comparison
    #show acetate:propionate:butyrate ratios as stacked bars
    for i, (_, crow) in enumerate(community_df.iterrows()):
        total = crow["acetate"] + crow["propionate"] + crow["butyrate"]
        if total == 0:
            continue
        ac_frac = crow["acetate"] / total
        pp_frac = crow["propionate"] / total
        bt_frac = crow["butyrate"] / total
        ax2.bar(i, ac_frac, color="#0072b2", width=0.7)
        ax2.bar(i, pp_frac, bottom=ac_frac, color="#e69f00", width=0.7)
        ax2.bar(i, bt_frac, bottom=ac_frac + pp_frac, color="#009e73", width=0.7)

    #add pipeline ratios
    n = len(community_df)
    for j, (_, row) in enumerate(pipeline_df.iterrows()):
        total = (row["acetate_mmol_gDW_hr"] + row["propionate_mmol_gDW_hr"]
                 + row["butyrate_mmol_gDW_hr"])
        ac_frac = row["acetate_mmol_gDW_hr"] / total
        pp_frac = row["propionate_mmol_gDW_hr"] / total
        bt_frac = row["butyrate_mmol_gDW_hr"] / total
        idx = n + j
        ax2.bar(idx, ac_frac, color="#0072b2", width=0.7, alpha=0.5, hatch="//")
        ax2.bar(idx, pp_frac, bottom=ac_frac, color="#e69f00", width=0.7,
                alpha=0.5, hatch="//")
        ax2.bar(idx, bt_frac, bottom=ac_frac + pp_frac, color="#009e73",
                width=0.7, alpha=0.5, hatch="//")

    all_labels = (community_df["source"]
                  .apply(lambda s: s.replace("_", "\n")).tolist()
                  + [cond_labels[c].replace("Dose_", "") + "\n(pipeline)"
                     for c in pipeline_df["condition"]])
    ax2.set_xticks(range(len(all_labels)))
    ax2.set_xticklabels(all_labels, fontsize=7)
    ax2.set_ylabel("Molar Fraction")
    ax2.set_title("B) SCFA Ratio: Community Models vs Pipeline")
    ax2.set_ylim(0, 1.05)

    #ratio legend
    ratio_legend = [
        Patch(facecolor="#0072b2", label="acetate"),
        Patch(facecolor="#e69f00", label="propionate"),
        Patch(facecolor="#009e73", label="butyrate"),
    ]
    ax2.legend(handles=ratio_legend, fontsize=8)

    fig.tight_layout()
    fig.savefig(figs_dir / "fig6_microbiome_crossvalidation.png")
    fig.savefig(figs_dir / "fig6_microbiome_crossvalidation.pdf")
    plt.close(fig)
    print("  fig6_microbiome_crossvalidation.png / .pdf")


def main():




    paths = build_paths()

    #load community scfa data
    community_path = paths.root / "data" / "inputs" / "agora2_community_scfa.csv"
    if not community_path.exists():
        print("  no community SCFA file found — skipping microbiome integration")
        return

    print("step 02a: microbiome integration & cross-validation")
    community_df = load_community_scfa(community_path)
    print(f"  loaded {len(community_df)} community SCFA predictions")

    pipeline_df = load_pipeline_conditions(paths.results)

    #cross-validate
    xval = cross_validate(community_df, pipeline_df)
    xval_path = paths.results / "microbiome_crossvalidation.csv"
    xval.to_csv(xval_path, index=False)
    print(f"  wrote cross-validation to {xval_path}")
    print(xval[["source", "community_total", "nearest_condition",
                "euclidean_distance"]].to_string(index=False))

    #generate figure
    fig_cross_validation(community_df, pipeline_df, paths.figs_dir)

    #summary statistics
    comm_totals = community_df["acetate"] + community_df["propionate"] + community_df["butyrate"]
    pipe_totals = (pipeline_df["acetate_mmol_gDW_hr"]
                   + pipeline_df["propionate_mmol_gDW_hr"]
                   + pipeline_df["butyrate_mmol_gDW_hr"])
    print(f"\n  community SCFA range: {comm_totals.min():.1f} – {comm_totals.max():.1f} mmol/gDW/hr")
    print(f"  pipeline condition range: {pipe_totals.min():.1f} – {pipe_totals.max():.1f} mmol/gDW/hr")
    overlap = (comm_totals.min() < pipe_totals.max()) and (comm_totals.max() > pipe_totals.min())
    print(f"  ranges overlap: {'YES' if overlap else 'NO'}")


if __name__ == "__main__":
    main()
