"""Generate manuscript Figure S1 pipeline schematic."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

# Colors tuned to match the requested reference style.
COL_INPUT = "#AFC0D6"
COL_COMP = "#E7DDB3"
COL_OUTPUT = "#BBB5DA"
COL_VALID = "#96DAB9"
EDGE = "#3A4C66"
ARROW = "#3A4C66"
SECTION = "#8EA1B1"

# ── Helper ────────────────────────────────────────────────────────
def _box(ax, cx, cy, w, h, title, subtitle, facecolor, fontsize_t=13, fontsize_s=8.5):
    """Draw a rounded box with title and subtitle."""
    rect = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                          boxstyle="round,pad=0.12", facecolor=facecolor,
                 edgecolor=EDGE, linewidth=1.5, zorder=2)
    ax.add_patch(rect)
    ax.text(cx, cy + 0.14, title, ha='center', va='center',
         fontsize=fontsize_t, fontweight='bold', family='DejaVu Sans',
         color='#1F2E42', zorder=3)
    if subtitle:
     ax.text(cx, cy - 0.22, subtitle, ha='center', va='center',
          fontsize=fontsize_s, fontstyle='italic', color='#4C5F75',
          family='DejaVu Sans', zorder=3)


def _arrow(ax, x1, y1, x2, y2):
    """Draw a directional connector."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->,head_width=0.22,head_length=0.1",
                                color=ARROW, lw=1.7),
                zorder=1)


# Figure setup
fig, ax = plt.subplots(figsize=(9.7, 7.7))
ax.set_xlim(0, 14)
ax.set_ylim(0, 10)
ax.axis('off')

BW, BH = 2.6, 1.0

# Section A
ax.text(0.2, 9.55, "A", fontsize=20, fontweight='bold', family='DejaVu Sans', color='#384A5E')
ax.text(0.2, 9.15, "Pipeline", fontsize=12, fontstyle='italic', family='DejaVu Sans', color=SECTION)

ya = 8.25
xs_a = [2.1, 5.3, 8.5, 11.7]

_box(ax, xs_a[0], ya, BW, BH, "Dietary Input",
     "S. affinis stachyose\n(25 / 50 / 100 g tubers)", COL_INPUT)
_box(ax, xs_a[1], ya, BW, BH, "SCFA Translation",
     "Ac : Pr : Bu\n\u2248 65 : 23 : 13", COL_COMP)
_box(ax, xs_a[2], ya, BW, BH, "Dual-Model FBA",
     "Recon3D  +  Human-GEM\n(ATPM objective)", COL_COMP)
_box(ax, xs_a[3], ya, BW, BH, "Predicted ATPM",
     "Dose-dependent ATP\nmaintenance flux", COL_OUTPUT)

for i in range(3):
    _arrow(ax, xs_a[i] + BW/2 + 0.05, ya, xs_a[i+1] - BW/2 - 0.05, ya)

# Section B
ax.text(0.2, 6.8, "B", fontsize=20, fontweight='bold', family='DejaVu Sans', color='#384A5E')
ax.text(0.2, 6.4, "Robustness", fontsize=12, fontstyle='italic', family='DejaVu Sans', color=SECTION)

yb = 5.7
xs_b = [3.1, 6.9, 10.7]

_box(ax, xs_b[0], yb, 3.0, BH, "Sensitivity Analysis",
     "OAT sweeps (6 params)\nRatio sensitivity profiles", COL_VALID)
_box(ax, xs_b[1], yb, 3.0, BH, "Solution Robustness",
     "FVA (90\u2013100% thresholds)\npFBA flux minimisation", COL_VALID)
_box(ax, xs_b[2], yb, 3.0, BH, "Pathway Rescue",
     "PPCOACm diagnosis\nModel convergence test", COL_VALID)

# Connections from Dual-Model FBA to robustness modules.
for xb in xs_b:
    _arrow(ax, xs_a[2], ya - BH/2 - 0.02, xb, yb + BH/2 + 0.02)

# Section C
ax.text(0.2, 3.9, "C", fontsize=20, fontweight='bold', family='DejaVu Sans', color='#384A5E')
ax.text(0.2, 3.5, "Validation", fontsize=12, fontstyle='italic', family='DejaVu Sans', color=SECTION)

yc = 2.9
_box(
    ax, 7.0, yc, 7.2, 1.1, "External Validation",
    "Predicted vs. published ATP yields\n(Bergman 1990; Remesy et al. 1986)",
    COL_VALID, fontsize_t=18, fontsize_s=12
)
_arrow(ax, xs_b[1], yb - BH/2 - 0.02, 7.0, yc + 0.58)

# Legend
legend_y = 1.0
legend_items = [
    (COL_INPUT, "Input data"),
    (COL_COMP,  "Computation"),
    (COL_OUTPUT, "Output"),
    (COL_VALID, "Validation / Robustness"),
]
lx = 1.7
for col, label in legend_items:
    rect = FancyBboxPatch((lx - 1.0, legend_y - 0.25), 2.2, 0.5,
                          boxstyle="round,pad=0.08", facecolor=col,
                          edgecolor=EDGE, linewidth=1.0, zorder=2)
    ax.add_patch(rect)
    ax.text(lx + 0.1, legend_y, label, ha='center', va='center',
            fontsize=12, family='DejaVu Sans', zorder=3)
    lx += 3.1

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
