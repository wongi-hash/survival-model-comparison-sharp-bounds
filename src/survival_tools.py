"""Shared helpers for the applications and the repeated-split analysis.

- event-stratified development/calibration split;
- ridge Cox predictors (Efron ties, interpolated Breslow baseline) evaluated at a fixed horizon;
- IPCW estimate of the target Brier-risk difference with a marginal Kaplan-Meier censoring distribution.
"""
import hashlib
from pathlib import Path

import numpy as np

from cox_ridge import fit_cox, h0_tau_interp, predict as _cox_predict

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORK = ROOT / "work"        # patient-level intermediate files (created locally, not distributed)
RESULTS = ROOT / "results"  # aggregate outputs


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_sha256(path, expected):
    if sha256(path) != expected:
        raise RuntimeError(f"SHA-256 mismatch for {path}; see data/README.md")


def stratified_split(event, seed):
    """Event-stratified random split: one third of each event stratum goes to calibration."""
    rng = np.random.default_rng(seed)
    cal = []
    for k in (0, 1):
        idx = rng.permutation(np.flatnonzero(event == k))
        cal.extend(idx[: int(round(len(idx) / 3))])
    cal = np.sort(np.asarray(cal))
    dev = np.setdiff1d(np.arange(len(event)), cal)
    return dev, cal


def fit_ridge_cox(x, time, event, tau, pen=0.01):
    """Ridge Cox model (objective -loglik/n + 0.5*pen*||b||^2 on standardized covariates, Efron ties)
    with the Breslow cumulative baseline hazard linearly interpolated at tau."""
    lock = fit_cox(x, time, event, pen=pen, ties="efron", scale="mean")
    lock["H0_tau"] = h0_tau_interp(lock, x, time, event, tau)
    return lock


def predict(lock, x):
    """Predicted survival probability at the horizon used in fit_ridge_cox."""
    return _cox_predict(lock, lock["H0_tau"], x)


def ipcw_brier_contrast(t, e, q0, q1, tau):
    """Candidate-minus-reference Brier risk at tau (see ipcw_brier_scores)."""
    b0, b1 = ipcw_brier_scores(t, e, q0, q1, tau)
    return b1 - b0


def ipcw_brier_scores(t, e, q0, q1, tau):
    """IPCW Brier risks of the reference and candidate at tau with a marginal Kaplan-Meier censoring
    distribution: G(t-) for events at or before tau and G(tau) for patients followed beyond tau."""
    t = np.asarray(t, float); e = np.asarray(e, int)
    o = np.argsort(t, kind="mergesort"); ts = t[o]; cs = 1 - e[o]
    ut, inv = np.unique(ts, return_inverse=True)
    n_at = len(ts) - np.searchsorted(ts, ut, side="left")
    dc = np.bincount(inv, weights=cs, minlength=len(ut))
    surv = np.cumprod(1 - dc / n_at)

    def g_left(q):
        i = np.searchsorted(ut, q, side="left") - 1; g = np.ones(len(q)); ok = i >= 0; g[ok] = surv[i[ok]]; return g

    def g_at(q):
        i = np.searchsorted(ut, q, side="right") - 1; g = np.ones(len(q)); ok = i >= 0; g[ok] = surv[i[ok]]; return g

    w = np.zeros(len(t)); died = (t <= tau) & (e == 1); alive = t > tau
    w[died] = 1 / np.maximum(g_left(t[died]), 1e-8)
    w[alive] = 1 / np.maximum(g_at(np.full(alive.sum(), tau)), 1e-8)
    y = alive.astype(float)
    return float(np.mean(w * (y - q0) ** 2)), float(np.mean(w * (y - q1) ** 2))


def retrospective_interval(t, e, q0, q1, tau, seed, boot=999):
    """Point estimate and percentile 95% interval from `boot` resamples of target patients."""
    est = ipcw_brier_contrast(t, e, q0, q1, tau)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(boot):
        i = rng.integers(0, len(t), len(t))
        bs.append(ipcw_brier_contrast(t[i], e[i], q0[i], q1[i], tau))
    lo, hi = np.quantile(bs, [0.025, 0.975])
    return est, float(lo), float(hi)
