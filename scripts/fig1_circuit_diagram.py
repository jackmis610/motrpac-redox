"""
Figure 1: Vandiver–Neufer thermodynamic circuit schematic
Creates figures/fig1_neufer_circuit.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

fig, ax = plt.subplots(figsize=(13, 8.5))
ax.set_xlim(0, 13); ax.set_ylim(0, 8.5)
ax.axis("off")
fig.patch.set_facecolor("white")

# ── Palette ────────────────────────────────────────────────
C_ETS    = "#2166AC"   # blue – ETS complexes
C_ATP    = "#762A83"   # purple – ATP synthase
C_NNT    = "#8C510A"   # brown – NNT
C_GSH    = "#E08214"   # orange – glutathione
C_TRX    = "#01665E"   # teal – thioredoxin
C_QPOOL  = "#4DAC26"   # green – Q-pool feeders
C_MEM    = "#CCCCCC"   # membrane grey
C_ARROW  = "#444444"
C_H2O2   = "#D6604D"   # red – ROS
C_TEXT   = "#222222"

def box(ax, x, y, w, h, color, text, fontsize=8.5, text_color="white", lw=1.2, style="round,pad=0.1"):
    rect = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle=style, linewidth=lw,
                          edgecolor="white", facecolor=color, zorder=3)
    ax.add_patch(rect)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            color=text_color, fontweight="bold", zorder=4, linespacing=1.4)

def arrow(ax, x1, y1, x2, y2, color=C_ARROW, lw=1.6, style="-|>", label=None, label_x=None, label_y=None, fs=7.5):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color,
                                lw=lw, mutation_scale=14), zorder=2)
    if label:
        lx = label_x if label_x else (x1+x2)/2
        ly = label_y if label_y else (y1+y2)/2
        ax.text(lx, ly, label, ha="center", va="center", fontsize=fs,
                color=color, fontstyle="italic", zorder=5,
                bbox=dict(fc="white", ec="none", pad=1.5))

# ══════════════════════════════════════════════════════════════
# IMM band
# ══════════════════════════════════════════════════════════════
imm = mpatches.FancyBboxPatch((0.4, 3.5), 9.2, 0.9,
                               boxstyle="round,pad=0.1", linewidth=0,
                               facecolor=C_MEM, alpha=0.35, zorder=1)
ax.add_patch(imm)
ax.text(0.65, 3.95, "IMM", ha="center", va="center", fontsize=7,
        color="grey", fontstyle="italic")

# ══════════════════════════════════════════════════════════════
# MATRIX side labels
# ══════════════════════════════════════════════════════════════
ax.text(0.22, 5.5, "Matrix", ha="center", va="center", fontsize=8,
        color="grey", fontstyle="italic", rotation=90)
ax.text(0.22, 3.0, "IMS", ha="center", va="center", fontsize=8,
        color="grey", fontstyle="italic", rotation=90)
ax.axhline(3.45, xmin=0.03, xmax=0.74, color=C_MEM, lw=0.8, alpha=0.6)
ax.axhline(4.45, xmin=0.03, xmax=0.74, color=C_MEM, lw=0.8, alpha=0.6)

# ══════════════════════════════════════════════════════════════
# Q-pool feeders (left column, matrix)
# ══════════════════════════════════════════════════════════════
ax.text(1.2, 7.7, "Q-pool feeders", ha="center", fontsize=7.5,
        color=C_QPOOL, fontweight="bold")
for i, g in enumerate(["FADH₂\n(SDH, ETF,\nGPD2, DHODH)",
                        "NADH\n(CI)"]):
    box(ax, 1.2, 7.1 - i*1.05, 1.55, 0.75, C_QPOOL, g, fontsize=7, text_color="white")

# Arrow: feeders → ubiquinone pool
ax.text(1.2, 5.55, "Ubiquinone\npool (CoQ)", ha="center", va="center",
        fontsize=7.5, color=C_QPOOL, fontweight="bold",
        bbox=dict(fc="#E8F4E8", ec=C_QPOOL, lw=0.8, boxstyle="round,pad=0.2"))
arrow(ax, 1.2, 5.95, 1.2, 5.78, color=C_QPOOL, lw=1.2)
arrow(ax, 1.2, 6.45, 1.2, 5.78, color=C_QPOOL, lw=1.2)

# ══════════════════════════════════════════════════════════════
# ETS complexes on membrane
# ══════════════════════════════════════════════════════════════
for cx, label in [(2.1, "CI\n(NDUF*)"),
                  (3.2, "CIII\n(UQCR*)"),
                  (4.3, "CIV\n(COX*)"),
                  (5.3, "Cyt c")]:
    box(ax, cx, 3.95, 0.85, 0.72, C_ETS, label, fontsize=7.5)

# Proton arrows (upward = pumping into IMS)
for cx in [2.1, 3.2, 4.3]:
    ax.annotate("", xy=(cx, 3.45), xytext=(cx, 3.95 - 0.36),
                arrowprops=dict(arrowstyle="-|>", color=C_ETS, lw=1.4,
                                mutation_scale=10), zorder=2)
ax.text(3.2, 2.9, "H⁺  →  Δψ + ΔpH  =  ΔGpmf", ha="center", va="center",
        fontsize=8, color=C_ETS, fontstyle="italic")

# CoQ → CIII electron flow
arrow(ax, 1.65, 5.55, 2.65, 4.32, color=C_QPOOL, lw=1.2,
      style="-|>", label="e⁻", label_x=2.1, label_y=5.05)
# Cyt c connecting CIII → CIV
arrow(ax, 3.65, 3.95, 3.75, 3.95, color=C_ETS, lw=1.0)
# O₂ → H₂O at CIV
ax.text(4.3, 5.1, "O₂ → H₂O", ha="center", fontsize=7.5,
        color=C_ETS, fontstyle="italic")
arrow(ax, 4.3, 4.85, 4.3, 4.32, color=C_ETS, lw=1.0, style="-|>")

# H₂O₂ leak arrow (from CIII)
ax.annotate("", xy=(3.2, 5.7), xytext=(3.2, 4.32),
            arrowprops=dict(arrowstyle="-|>", color=C_H2O2,
                            lw=1.3, linestyle="dashed",
                            mutation_scale=12), zorder=2)
ax.text(3.55, 5.1, "H₂O₂\n(leak)", ha="left", fontsize=7.5,
        color=C_H2O2, fontweight="bold")

# ══════════════════════════════════════════════════════════════
# ATP synthase
# ══════════════════════════════════════════════════════════════
box(ax, 6.25, 3.95, 0.95, 0.72, C_ATP, "ATP\nsynthase\n(ATP5*)", fontsize=7)
ax.annotate("", xy=(6.25, 4.32), xytext=(6.25, 3.45),
            arrowprops=dict(arrowstyle="-|>", color=C_ATP,
                            lw=1.5, mutation_scale=12), zorder=2)
arrow(ax, 6.25, 4.68, 6.25, 5.35, color=C_ATP, lw=1.4,
      label="ΔGATP", label_x=6.65, label_y=5.05)
ax.text(6.25, 5.65, "ATP\nproduction", ha="center", fontsize=7.5,
        color=C_ATP, fontweight="bold")

# ══════════════════════════════════════════════════════════════
# NNT
# ══════════════════════════════════════════════════════════════
box(ax, 7.55, 3.95, 0.95, 0.72, C_NNT, "NNT", fontsize=10)
ax.annotate("", xy=(7.55, 4.32), xytext=(7.55, 3.45),
            arrowprops=dict(arrowstyle="-|>", color=C_NNT,
                            lw=1.5, mutation_scale=12), zorder=2)
# NNT arrow to NADPH pool
arrow(ax, 7.55, 4.68, 7.55, 5.35, color=C_NNT, lw=1.4,
      label="ΔGNADPH", label_x=8.1, label_y=5.0)

# NADH → NADPH label below NNT
ax.text(7.55, 3.2, "NADH → NADPH", ha="center", fontsize=7,
        color=C_NNT, fontstyle="italic")
ax.text(7.55, 2.85, "(uses ΔGpmf)", ha="center", fontsize=7,
        color=C_NNT, fontstyle="italic")

# NADPH pool
ax.text(7.55, 5.7, "Mitochondrial\nNADPH pool", ha="center", va="center",
        fontsize=8, color=C_NNT, fontweight="bold",
        bbox=dict(fc="#F5E8CC", ec=C_NNT, lw=0.9, boxstyle="round,pad=0.25"))

# ══════════════════════════════════════════════════════════════
# GSH circuit
# ══════════════════════════════════════════════════════════════
box(ax, 9.4, 6.5, 1.7, 0.85, C_GSH,
    "GSH system\nGPx1/4 · GSR\nGCLC · GCLM · Glrx2", fontsize=7)
arrow(ax, 7.95, 5.7, 8.7, 6.5, color=C_NNT, lw=1.3,
      label="NADPH", label_x=8.15, label_y=6.25)
# H₂O₂ → GSH
arrow(ax, 3.7, 5.7, 8.55, 6.5, color=C_H2O2, lw=1.1,
      label="H₂O₂", label_x=6.2, label_y=6.3)

# ══════════════════════════════════════════════════════════════
# TRX circuit
# ══════════════════════════════════════════════════════════════
box(ax, 9.4, 5.3, 1.7, 0.85, C_TRX,
    "TRX system\nTxnrd2 · Txn2\nPrdx3 · Prdx5", fontsize=7)
arrow(ax, 7.95, 5.55, 8.55, 5.3, color=C_NNT, lw=1.3)
arrow(ax, 3.7, 5.7, 8.55, 5.3, color=C_H2O2, lw=1.1)

# Outcome
ax.text(11.25, 5.9, "Proteome\nredox tone\n(~17,000 Cys)", ha="center", va="center",
        fontsize=8, color=C_TEXT, fontweight="bold",
        bbox=dict(fc="#F0F0F0", ec="#999999", lw=0.8, boxstyle="round,pad=0.3"))
arrow(ax, 10.27, 6.5, 10.82, 6.1, color=C_GSH, lw=1.1)
arrow(ax, 10.27, 5.3, 10.82, 5.7, color=C_TRX, lw=1.1)

# ══════════════════════════════════════════════════════════════
# SOD
# ══════════════════════════════════════════════════════════════
box(ax, 5.5, 5.9, 1.1, 0.65, "#B35806", "SOD1/2\n(O₂•⁻→H₂O₂)", fontsize=7)
arrow(ax, 4.3, 5.7, 5.0, 5.9, color=C_H2O2, lw=1.0)
arrow(ax, 5.9, 5.9, 6.4, 6.5, color=C_H2O2, lw=1.0)

# ══════════════════════════════════════════════════════════════
# Bottom thermodynamic flow panel
# ══════════════════════════════════════════════════════════════
y0 = 0.9
ax.add_patch(mpatches.FancyBboxPatch((0.5, 0.3), 11.8, 1.05,
             boxstyle="round,pad=0.1", facecolor="#F7F7F7",
             edgecolor="#CCCCCC", lw=1.0))
flow_items = [
    (1.2,  y0, "ΔG_redox\n(substrate\noxidation)",   C_QPOOL),
    (3.0,  y0, "ΔGpmf\n(proton-motive\nforce)",       C_ETS),
    (5.1,  y0, "ΔG_ATP\n(ATP\nsynthesis)",              C_ATP),
    (7.2,  y0, "ΔG_NADPH\n(NNT-driven\nregeneration)", C_NNT),
    (9.4,  y0, "GSH/TRX\nbuffering\ncapacity",          C_GSH),
    (11.35,y0, "Proteome\nredox\ntone",                  C_TEXT),
]
for x, y, lbl, col in flow_items:
    box(ax, x, y, 1.55, 0.72, col, lbl, fontsize=7.5)
for i in range(len(flow_items)-1):
    x1 = flow_items[i][0]   + 0.78
    x2 = flow_items[i+1][0] - 0.78
    ax.annotate("", xy=(x2, y0), xytext=(x1, y0),
                arrowprops=dict(arrowstyle="-|>", color="#888888",
                                lw=1.4, mutation_scale=13))

# ══════════════════════════════════════════════════════════════
# Title and legend
# ══════════════════════════════════════════════════════════════
ax.text(6.5, 8.25,
        "Vandiver–Neufer Thermodynamic Framework: Mitochondrial Redox Circuit",
        ha="center", va="center", fontsize=11, fontweight="bold", color=C_TEXT)

legend_items = [
    mpatches.Patch(fc=C_ETS,    label="ETS complexes (CI/CIII/CIV, Cyt c)"),
    mpatches.Patch(fc=C_ATP,    label="ATP synthase (ATP5*)"),
    mpatches.Patch(fc=C_NNT,    label="NNT (ΔGpmf → NADPH)"),
    mpatches.Patch(fc=C_GSH,    label="Glutathione system"),
    mpatches.Patch(fc=C_TRX,    label="Thioredoxin system"),
    mpatches.Patch(fc=C_QPOOL,  label="Q-pool feeder enzymes"),
    mpatches.Patch(fc=C_H2O2,   label="ROS (H₂O₂, O₂•⁻)"),
]
ax.legend(handles=legend_items, loc="upper right",
          bbox_to_anchor=(12.9, 8.4), fontsize=7,
          framealpha=0.9, edgecolor="#CCCCCC",
          ncol=1, handlelength=1.2, handleheight=0.9)

plt.tight_layout(pad=0.2)
plt.savefig("figures/fig1_neufer_circuit.png", dpi=220,
            bbox_inches="tight", facecolor="white")
plt.savefig("figures/fig1_neufer_circuit.pdf",
            bbox_inches="tight", facecolor="white")
print("Saved figures/fig1_neufer_circuit.png")
print("Saved figures/fig1_neufer_circuit.pdf")
