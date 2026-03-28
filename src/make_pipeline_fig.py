"""Generate pipeline schematic figure (Figure S1).

Three-panel layout:
  A — Pipeline: Dietary Input → SCFA Translation → Dual-Model FBA → Predicted ATPM
  B — Robustness: Sensitivity Analysis, Solution Robustness, Pathway Rescue
  C — Validation: External Validation against published ATP yields
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

# ── Colors ────────────────────────────────────────────────────────
COL_INPUT  = "#C5CAE9"   # light indigo  — input data
COL_COMP   = "#B2DFDB"   # light teal    — computation
COL_OUTPUT = "#FFF9C4"   # light yellow  — output
COL_VALID  = "#F3E5F5"   # light purple  — validation / robustness
EDGE       = "#455A64"   # dark blue-grey
ARROW_COL  = "#37474F"
SECTION_COL = "#90A4AE"  # section labels

# ── Helper ────────────────────────────────────────────────────────
def _box(ax, cx, cy, w, h, title, subtitle, facecolor, fontsize_t=11,
         fontsize_s=7.5):
    """Draw a rounded box with bold title + italic subtitle."""
    rect = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                          boxstyle="round,pad=0.12", facecolor=facecolor,
                          edgecolor=EDGE, linewidth=1.3, zorder=2)
    ax.add_patch(rect)
    ax.text(cx, cy + 0.18, title, ha='center', va='center',
            fontsize=fontsize_t, fontweight='bold', family='serif', zorder=3)
    if subtitle:
        ax.text(cx, cy - 0.22, subtitle, ha='center', va='center',
                fontsize=fontsize_s, fontstyle='italic', color='#555555',
                family='serif', zorder=3)


def _arrow(ax, x1, y1, x2, y2):
    """Draw a simple arrow between two points."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->,head_width=0.25,head_length=0.12",
                                color=ARROW_COL, lw=1.8),
                zorder=1)


# ── Figure setup ──────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 9.5))
ax.set_xlim(-0.5, 10.5)
ax.set_ylim(-0.5, 10)
ax.axis('off')

BW, BH = 2.0, 1.1   # box width / height
BWs = 2.4            # wider boxes for section B

# ══════════════════════════════════════════════════════════════════
# Section A — Pipeline
# ══════════════════════════════════════════════════════════════════
ax.text(0.0, 9.5, "A", fontsize=16, fontweight='bold', family='serif')
ax.text(0.35, 9.5, "Pipeline", fontsize=11, fontstyle='italic',
        color=SECTION_COL, family='serif', va='center')

ya = 8.5
xs_a = [1.2, 3.7, 6.2, 8.8]

_box(ax, xs_a[0], ya, BW, BH, "Dietary Input",
     "S. affinis stachyose\n(25 / 50 / 100 g tubers)", COL_INPUT)
_box(ax, xs_a[1], ya, BW, BH, "SCFA Translation",
     "Ac : Pr : Bu\n\u2248 65 : 23 : 13", COL_INPUT)
_box(ax, xs_a[2], ya, BW, BH, "Dual-Model FBA",
     "Recon3D + Human-GEM\n(ATPM objective)", COL_COMP)
_box(ax, xs_a[3], ya, BW, BH, "Predicted ATPM",
     "Dose-dependent ATP\nmaintenance flux", COL_OUTPUT)

for i in range(3):
    _arrow(ax, xs_a[i] + BW/2 + 0.05, ya, xs_a[i+1] - BW/2 - 0.05, ya)

# ══════════════════════════════════════════════════════════════════
# Section B — Robustness
# ══════════════════════════════════════════════════════════════════
ax.text(0.0, 6.8, "B", fontsize=16, fontweight='bold', family='serif')
ax.text(0.35, 6.8, "Robustness", fontsize=11, fontstyle='italic',
        color=SECTION_COL, family='serif', va='center')

yb = 5.5
xs_b = [1.8, 5.0, 8.2]

_box(ax, xs_b[0], yb, BWs, BH, "Sensitivity Analysis",
     "OAT sweeps (6 params)\nRatio sensitivity profiles", COL_VALID)
_box(ax, xs_b[1], yb, BWs, BH, "Solution Robustness",
     "FVA (90\u2013100% thresholds)\npFBA flux minimisation", COL_VALID)
_box(ax, xs_b[2], yb, BWs, BH, "Pathway Rescue",
     "PPCOACm diagnosis\nModel convergence test", COL_VALID)

# Arrows from pipeline row A down to robustness row B
for xb in xs_b:
    _arrow(ax, 6.2, ya - BH/2 - 0.05, xb, yb + BH/2 + 0.05)

# ══════════════════════════════════════════════════════════════════
# Section C — Validation
# ══════════════════════════════════════════════════════════════════
ax.text(0.0, 3.8, "C", fontsize=16, fontweight='bold', family='serif')
ax.text(0.35, 3.8, "Validation", fontsize=11, fontstyle='italic',
        color=SECTION_COL, family='serif', va='center')

yc = 2.5
_box(ax, 5.0, yc, 3.6, 1.3, "External Validation",
     "Predicted vs. published ATP yields\n(Bergman 1990; R\u00e9m\u00e9sy et al. 1986)",
     COL_VALID, fontsize_t=12, fontsize_s=8.5)

# ══════════════════════════════════════════════════════════════════
# Legend
# ══════════════════════════════════════════════════════════════════
legend_y = 0.3
legend_items = [
    (COL_INPUT, "Input data"),
    (COL_COMP,  "Computation"),
    (COL_OUTPUT, "Output"),
    (COL_VALID, "Validation / Robustness"),
]
lx = 1.0
for col, label in legend_items:
    rect = FancyBboxPatch((lx - 0.55, legend_y - 0.25), 1.8, 0.5,
                          boxstyle="round,pad=0.08", facecolor=col,
                          edgecolor=EDGE, linewidth=1.0, zorder=2)
    ax.add_patch(rect)
    ax.text(lx + 0.35, legend_y, label, ha='center', va='center',
            fontsize=8, family='serif', zorder=3)
    lx += 2.4

plt.tight_layout()

# ── Save ──────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
out_png = os.path.join(ROOT, "outputs", "figs", "FigS1.png")
out_pdf = os.path.join(ROOT, "outputs", "figs", "FigS1.pdf")
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved {out_png}")
print(f"Saved {out_pdf}")
