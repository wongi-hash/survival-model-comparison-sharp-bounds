"""Ridge-penalized Cox proportional hazards model in numpy.

Objective: -loglik/n + 0.5*pen*||beta||^2 (scale="mean") or -loglik + 0.5*pen*||beta||^2 (scale="sum"),
on covariates standardized with their sample mean and standard deviation. Ties: Efron (default) or Breslow.
"""
import numpy as np


def fit_cox(x, time, event, pen=0.01, ties="efron", scale="mean"):
    mu, sd = x.mean(0), x.std(0); sd[sd < 1e-12] = 1.0; z = (x - mu) / sd
    n, p = z.shape; ut = np.unique(time[event == 1]); beta = np.zeros(p)
    groups = [(np.flatnonzero(time >= t), np.flatnonzero((time == t) & (event == 1))) for t in ut]
    c = (1.0 / n) if scale == "mean" else 1.0
    for _ in range(100):
        r = np.exp(z @ beta); g = np.zeros(p); H = np.zeros((p, p))
        for R, D in groups:
            d = len(D); zR, rR = z[R], r[R]; S0 = rR.sum(); S1 = rR @ zR; S2 = (zR * rR[:, None]).T @ zR
            zD, rD = z[D], r[D]; T0 = rD.sum(); T1 = rD @ zD; T2 = (zD * rD[:, None]).T @ zD
            g += zD.sum(0)
            for l in range(d):
                f = (l / d) if ties == "efron" else 0.0
                den = S0 - f * T0; m = (S1 - f * T1) / den
                g -= m; H -= (S2 - f * T2) / den - np.outer(m, m)
        g = c * g - pen * beta; H = c * H - pen * np.eye(p)
        step = np.linalg.solve(H, g); beta = beta - step
        if np.max(np.abs(step)) < 1e-11: break
    return {"mu": mu, "sd": sd, "beta": beta}


def h0_tau_interp(lock, x, time, event, tau):
    """Breslow cumulative baseline hazard on all distinct observed times, linearly interpolated at tau."""
    risk = np.exp(((x - lock["mu"]) / lock["sd"]) @ lock["beta"])
    ut = np.unique(time); H = []; h = 0.0
    for t in ut:
        d = ((time == t) & (event == 1)).sum()
        if d: h += d / risk[time >= t].sum()
        H.append(h)
    return float(np.interp(tau, ut, np.asarray(H)))


def predict(lock, h0, x):
    """Survival probability exp{-H0(tau) exp(x'beta)}, restricted to [1e-4, 1 - 1e-4]."""
    return np.clip(np.exp(-h0 * np.exp(((x - lock["mu"]) / lock["sd"]) @ lock["beta"])), 1e-4, 1 - 1e-4)
