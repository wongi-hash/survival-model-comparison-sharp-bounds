"""Figure 2: direct identification and separate-risk subtraction across sensitivity values for one
configuration per application (identified intervals: thick bars; 95% simultaneous envelopes: thin lines).
The configurations were chosen to illustrate differences in directional support; all configurations are in
tables/applications_fixed_configurations.csv.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
g = pd.read_csv(RES / "rotterdam_gbsg/direct_vs_separate.csv"); g = g[g.seed == 62]
t = pd.read_csv(RES / "rotterdam_temporal/seed52.csv")
n = pd.read_csv(RES / "nwts/seed62.csv")
rg = pd.read_csv(RES / "rotterdam_gbsg/retrospective_target_contrast.csv").set_index("seed").loc[62, "contrast"]
rt = pd.read_csv(RES / "rotterdam_temporal/retrospective_target_contrast.csv").set_index("seed").loc[52, "contrast"]
rn = pd.read_csv(RES / "nwts/retrospective_target_contrast.csv").set_index("seed").loc[62, "contrast"]
panels = [("(A) Rotterdam → GBSG\n(other hospital), configuration 62", g, rg),
          ("(B) Rotterdam 1978–88 → 1989–93\n(later patients), configuration 52", t, rt),
          ("(C) NWTS-3 → NWTS-4\n(next study), configuration 62", n, rn)]
fig, axes = plt.subplots(1, 3, figsize=(7.5, 5.4))
for ax, (title, d, retro) in zip(axes, panels):
    d = d.sort_values("gamma"); y = np.arange(len(d))
    for off, pre, col, lab in ((-0.17, "dir", "#1f3b73", "Direct"), (0.17, "sep", "#c0504d", "Separate")):
        L, U = d[f"L_{pre}"], d[f"U_{pre}"]
        lo, hi = d[f"{pre}_sim_lower"], d[f"{pre}_sim_upper"]
        ax.hlines(y + off, lo, hi, color=col, lw=1.3, label=f"{lab}: 95% simultaneous envelope")
        ax.hlines(y + off, L, U, color=col, lw=6, alpha=0.6, label=f"{lab}: identified set")
        act = d[f"action_{pre}"].to_numpy()
        for yi, a, h in zip(y + off, act, hi):
            if a != "DEFER": ax.text(h, yi, " C", color=col, va="center", ha="left", fontsize=8.5, fontweight="bold")
    ax.axvline(0, color="0.3", lw=0.8)
    ax.axvline(retro, color="0.5", lw=0.9, ls=":", label="Later-observed IPCW estimate")
    ax.set_yticks(y); ax.set_yticklabels([f"{v:g}" for v in d.gamma]); ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_title(title, fontsize=9); ax.tick_params(labelsize=8.5); ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.spines[["top", "right"]].set_visible(False)
axes[1].set_xlabel("Brier-risk difference (candidate − reference)", fontsize=10)
axes[0].set_ylabel("Sensitivity bound γ", fontsize=10)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8, columnspacing=1.2, handlelength=1.8)
fig.tight_layout(rect=(0, 0.09, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(Path(__file__).with_name(f"Figure2_applications.{ext}"), dpi=300, bbox_inches="tight")
print("Figure 2 written")
