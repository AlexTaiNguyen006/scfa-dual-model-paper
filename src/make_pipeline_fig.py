"""Generate pipeline schematic figure (Figure S2)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(figsize=(10, 2.8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 3)
ax.axis('off')

astyle = dict(arrowstyle="->,head_width=0.3,head_length=0.15",
              color="#2C5F8A", lw=2)

boxes = [
    (1.0, "S. affinis\nstachyose\n(dietary input)"),
    (3.2, "Colonic\nfermentation\n(SCFA production)"),
    (5.4, "Portal\ndelivery\n(SCFA uptake vectors)"),
    (7.6, "Hepatic FBA\n(Recon3D / Human-GEM\nATPM objective)"),
    (9.4, "Predicted\nATPM"),
]

for x, txt in boxes:
    ax.text(x, 1.5, txt, ha='center', va='center', fontsize=8.5,
            fontweight='bold', family='serif',
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#E8F0FE",
                      edgecolor="#2C5F8A", linewidth=1.5))

arrow_pairs = [(1.85, 2.35), (4.05, 4.55), (6.25, 6.75), (8.45, 8.85)]
for x1, x2 in arrow_pairs:
    ax.annotate('', xy=(x2, 1.5), xytext=(x1, 1.5), arrowprops=astyle)

ax.text(5.0, 0.25,
        "Three dose scenarios (Low / Mid / High)  \u00b7  "
        "Dual-model comparison  \u00b7  "
        "Sensitivity + FVA + Rescue analysis",
        ha='center', va='center', fontsize=7.5, fontstyle='italic',
        color='#555555', family='serif')

plt.tight_layout()
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
out = os.path.join(ROOT, "outputs", "figs",
                   "figS2_pipeline.png")
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved {out}")
