"""Run both gravity models and plot their differences; requires matplotlib/numpy."""
import argparse
import csv
import io
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binary", default="build/orbit-access")
parser.add_argument("--output", default="docs/orbit-comparison.png")
args = parser.parse_args()
models = {}
for model in ("two-body", "j2"):
    result = subprocess.run([args.binary, "propagate", "86400", "60", model],
                            check=True, capture_output=True, text=True)
    rows = list(csv.DictReader(io.StringIO(result.stdout)))
    models[model] = {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}

fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, layout="constrained")
for model, color in (("two-body", "#2166ac"), ("j2", "#b35806")):
    data = models[model]
    position = np.array([data[f"{axis}_km"] for axis in "xyz"]).T
    velocity = np.array([data[f"v{axis}_km_s"] for axis in "xyz"]).T
    momentum = np.cross(position, velocity)
    node = np.degrees(np.unwrap(np.arctan2(momentum[:, 0], -momentum[:, 1])))
    hours = data["time_s"] / 3600
    axes[0].plot(hours, node, label=model, color=color)
    axes[1].plot(hours, np.linalg.norm(position, axis=1), label=model, color=color)
    energy = data["energy_km2_s2"]
    axes[2].plot(hours, (energy - energy[0]) / abs(energy[0]), color=color)
axes[0].set(ylabel="Ascending node (deg)", title="One day from the same 7000-km, 51.6° initial orbit")
axes[0].legend()
axes[1].set(ylabel="Geocentric radius (km)")
axes[2].set(ylabel="Relative energy change", xlabel="Elapsed time (hours)")
for axis in axes:
    axis.grid(alpha=0.2)
    axis.spines[["top", "right"]].set_visible(False)
output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=160)
summary = {}
for model, data in models.items():
    summary[model] = {"max_relative_energy_change": float(np.max(np.abs(
        (data["energy_km2_s2"] - data["energy_km2_s2"][0]) / data["energy_km2_s2"][0])))}
output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary))
