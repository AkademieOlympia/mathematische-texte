import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent

files = [
    "axes_5e8.csv",
    "axes_5e9.csv",
    "axes_5e10.csv",
]

rows = []

for f in files:
    df = pd.read_csv(ROOT / f)
    N = int(df["N"].iloc[0])

    for dist, g in df.groupby("dist"):
        mean_d16 = g["d16"].mean()
        mean_d25 = g["d25"].mean()
        mean_d34 = g["d34"].mean()
        mean_purity = g["purity"].mean()
        mean_deff = g["d_eff"].mean()

        abs_means = {
            "d16": abs(mean_d16),
            "d25": abs(mean_d25),
            "d34": abs(mean_d34),
        }
        dominant = max(abs_means, key=abs_means.get)

        rows.append(
            {
                "N": N,
                "dist": dist,
                "pairs": len(g),
                "mean_d16": mean_d16,
                "mean_d25": mean_d25,
                "mean_d34": mean_d34,
                "mean_purity": mean_purity,
                "mean_d_eff": mean_deff,
                "dominant_from_mean": dominant,
            }
        )

out = pd.DataFrame(rows).sort_values(["N", "dist"])
out.to_csv(ROOT / "axes_aggregate.csv", index=False)

print(out.to_string(index=False))
