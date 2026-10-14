"""Подсчёт числа образцов на каждый тип ткани и построение барплота."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "meta_tissues_hidden.csv"
OUT = ROOT / "paper" / "tissue_sample_counts.png"

meta = pd.read_csv(META)
counts = meta["Tissue_Number"].value_counts().sort_index()

fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(counts.index.astype(str), counts.values, color="#4C72B0",
              edgecolor="black", linewidth=0.8, width=0.6)
ax.set_xlabel("Tissue number", fontsize=13)
ax.set_ylabel("Number of samples", fontsize=13)
ax.set_title("Samples per tissue type", fontsize=15, fontweight="bold")
ax.tick_params(axis="both", labelsize=12)
ax.grid(axis="y", alpha=0.3, linestyle="--")
ax.set_axisbelow(True)
ax.set_ylim(0, counts.max() * 1.15)
for bar, v in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width() / 2, v + 40, str(v),
            ha="center", va="bottom", fontsize=12, fontweight="bold")
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close(fig)

print("Counts per tissue:")
print(counts.to_string())
print(f"Total samples: {counts.sum()}")
print(f"Saved plot: {OUT}")