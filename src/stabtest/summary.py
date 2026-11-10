"""Collect the headline numbers quoted in the README into ``summary.json``."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def _ms(x: pd.Series) -> dict[str, float]:
    return {"mean": float(x.mean()), "std": float(x.std(ddof=1)) if len(x) > 1 else 0.0}


def summarise(res: Path) -> dict:
    """Read the CSVs in ``res`` and write ``res / 'summary.json'``."""
    out: dict = {}
    c = pd.read_csv(res / "curves.csv")
    out["curves"] = {
        m: {int(k): _ms(g["error"]) for k, g in d.groupby("copies")} for m, d in c.groupby("method")
    }
    h = pd.read_csv(res / "heatmap.csv")
    out["heatmap"] = {
        f"dephase={a},readout={b}": {
            "copies_needed": _ms(g["copies_needed"]),
            "agnostic_asymptotic_error": _ms(g["agnostic_asymptotic_error"]),
            "null_accept": _ms(g["null_accept"]),
            "alt_accept": _ms(g["alt_accept"]),
        }
        for (a, b), g in h.groupby(["dephase", "readout"])
    }
    s = pd.read_csv(res / "scaling.csv").replace(np.inf, np.nan)
    out["scaling"] = {
        f"n={n}|{z}|{m}": _ms(g["copies_needed"]) for (n, z, m), g in s.groupby(["n", "noise", "method"])
    }
    ratios = {}
    for (n, z), g in s.groupby(["n", "noise"]):
        a = g[g["method"].str.startswith("Pauli")].set_index("seed")["copies_needed"]
        b = g[g["method"].str.startswith("Bell")].set_index("seed")["copies_needed"]
        if len(a) and a.notna().all():
            ratios[f"n={n}|{z}"] = _ms(a / b.loc[a.index])
    out["shadow_over_bell_copy_ratio"] = ratios
    p = pd.read_csv(res / "critical_noise.csv")
    out["critical_noise"] = {int(n): _ms(g["p_crit"]) for n, g in p.groupby("n")}
    w = pd.read_csv(res / "waterfall.csv")
    piv = w.pivot(index="seed", columns="step", values="error")[list(dict.fromkeys(w["step"]))]
    out["waterfall_error"] = {k: _ms(piv[k]) for k in piv.columns}
    out["waterfall_delta"] = {k: _ms(piv.diff(axis=1)[k]) for k in piv.columns[1:]}
    (res / "summary.json").write_text(json.dumps(out, indent=2))
    return out
