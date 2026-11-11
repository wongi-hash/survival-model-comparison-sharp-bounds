"""Data-generating mechanism and estimators for the simulation study.

Source covariate X ~ N(0, 1); target covariate X ~ N(TXSHIFT, 1).
Fixed predictions q0(x) = expit(-0.15 + 0.45x), q1(x) = expit(0.10 + 0.80x).
Source survival beyond tau = 1 given Z = (q0, q1): m_S(z) = expit(-0.35 + 1.25 q0 - 0.70 q1).
Event times: Weibull with shape WK and S(tau | z) = m_S(z).
Censoring times: exponential with rate lambda_C * exp(ALPHA_C' z_std), z_std = (q0, q1) standardized.
Estimation: conditional Cox-Breslow censoring model given Z, IPCW-weighted logistic outcome model.
"""
import hashlib

import numpy as np

GAMMAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5]); TAU = 1.0; WK = 1.5
TXSHIFT = 0.5; ALPHA_C = np.array([1.0, -0.75]); G_MIN = 0.02
BOOT = 199


def expit(x): return 1.0/(1.0+np.exp(-np.clip(np.asarray(x, float), -35, 35)))
def logit(p): p = np.clip(np.asarray(p, float), 1e-6, 1-1e-6); return np.log(p/(1-p))
def qv(x): x = np.asarray(x, float); return expit(-0.15+0.45*x), expit(0.10+0.80*x)
def psrc(q0, q1): return expit(-0.35+1.25*q0-0.70*q1)
def stable_seed(*parts): return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def endpoints(q0, q1, ps):
    """Direct identified-interval endpoints at each gamma."""
    a = q1*q1-q0*q0; b = q1-q0; lp = logit(ps); out = np.empty((len(GAMMAS), 2))
    for j, g in enumerate(GAMMAS):
        pm, pp = expit(lp-g), expit(lp+g)
        out[j, 0] = np.mean(a-2*b*np.where(b >= 0, pp, pm))
        out[j, 1] = np.mean(a-2*b*np.where(b >= 0, pm, pp))
    return out


def calib_lambda(target_cens):
    """Baseline censoring rate giving the target overall censoring proportion in the source population."""
    rng = np.random.default_rng(20260907); x = rng.normal(0, 1, 500_000); q0, q1 = qv(x); ps = psrc(q0, q1)
    lam = -np.log(np.clip(ps, 1e-8, 1))/TAU**WK; t = (-np.log(rng.random(len(x)))/np.clip(lam, 1e-10, None))**(1/WK)
    z = np.column_stack([q0, q1]); z = (z-z.mean(0))/z.std(0); u = -np.log(rng.random(len(x))); mult = np.exp(z@ALPHA_C)
    lo, hi = 1e-4, 20.0
    for _ in range(60):
        mid = (lo+hi)/2; c = u/(mid*mult)
        if np.mean(c < t) < target_cens: lo = mid
        else: hi = mid
    return (lo+hi)/2


def population_truth():
    """Population direct endpoints from 2,000,000 target profiles."""
    rng = np.random.default_rng(990011); x = rng.normal(TXSHIFT, 1, 2_000_000); q0, q1 = qv(x); ps = psrc(q0, q1)
    return endpoints(q0, q1, ps)


TRUE_EP = population_truth()


def _cox(time, ce, z):
    time = np.asarray(time, float); ce = np.asarray(ce, int); z = np.asarray(z, float)
    mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0; x = (z-mu)/sd; et = np.unique(time[ce == 1])
    if len(et) == 0: return mu, sd, np.zeros(x.shape[1]), np.array([]), np.array([])
    asc = np.argsort(time); ta = time[asc]; xa = x[asc]; i0 = np.searchsorted(ta, et, 'left')
    d = np.array([int(((time == t) & (ce == 1)).sum()) for t in et], float)
    sumx = np.array([x[(time == t) & (ce == 1)].sum(0) for t in et]); xaxa = xa[:, :, None]*xa[:, None, :]
    beta = np.zeros(x.shape[1])
    for _ in range(40):
        ra = np.exp(np.clip(xa@beta, -25, 25)); r0 = np.cumsum(ra[::-1])[::-1]
        r1 = np.cumsum((ra[:, None]*xa)[::-1], 0)[::-1]; r2 = np.cumsum((ra[:, None, None]*xaxa)[::-1], 0)[::-1]
        S0 = r0[i0]; m = r1[i0]/S0[:, None]
        info = np.sum(d[:, None, None]*(r2[i0]/S0[:, None, None]-m[:, :, None]*m[:, None, :]), 0)
        sc = sumx.sum(0)-np.sum(d[:, None]*m, 0); step = np.linalg.pinv(info, rcond=1e-10)@sc
        mx = np.max(np.abs(step))
        if mx > 2: step *= 2/mx
        beta += step
        if np.max(np.abs(step)) < 1e-8: break
    ra = np.exp(np.clip(xa@beta, -25, 25)); r0 = np.cumsum(ra[::-1])[::-1]; inc = d/r0[i0]
    return mu, sd, beta, et, inc


def cipcw(time, event, z, tau):
    """Fixed-horizon IPCW weights from a Cox-Breslow censoring model given z."""
    time = np.asarray(time, float); event = np.asarray(event, int); z = np.asarray(z, float)
    mu, sd, beta, et, dh = _cox(time, 1-event, z); ph = np.exp(np.clip(((z-mu)/sd)@beta, -25, 25)); cum = np.cumsum(dh)
    def h0(q, left):
        idx = np.searchsorted(et, np.asarray(q, float), side=('left' if left else 'right'))-1
        o = np.zeros(len(np.atleast_1d(q))); ok = idx >= 0; o[ok] = cum[np.minimum(idx[ok], len(cum)-1)]; return o
    gtl = np.clip(np.exp(-h0(time, True)*ph), G_MIN, 1.0); gtau = np.clip(np.exp(-h0(np.repeat(tau, len(time)), False)*ph), G_MIN, 1.0)
    y = (time > tau).astype(float); w = np.zeros(len(time)); w[time > tau] = 1/gtau[time > tau]
    f = (time <= tau) & (event == 1); w[f] = 1/gtl[f]
    return y, w


def wlog(x, y, w):
    """Weighted logistic regression by Newton-Raphson."""
    xd = np.column_stack([np.ones(len(x)), x]); b = np.zeros(xd.shape[1])
    for _ in range(60):
        p = expit(xd@b); ww = np.clip(w*p*(1-p), 1e-10, None)
        step = np.linalg.pinv((xd*ww[:, None]).T@xd, rcond=1e-10)@(xd.T@(w*(y-p)))
        mx = np.max(np.abs(step))
        if mx > 3: step *= 3/mx
        b += step
        if np.max(np.abs(step)) < 1e-8: break
    return b


def fitps(ct, ce, cq0, cq1, tq0, tq1):
    """Fitted source survival probabilities m_S at the target profiles."""
    z = np.column_stack([cq0, cq1]); y, w = cipcw(ct, ce, z, TAU); mu = z.mean(0); sd = z.std(0); sd[sd < 1e-8] = 1.0
    beta = wlog((z-mu)/sd, y, w); zt = np.column_stack([tq0, tq1])
    return np.clip(expit(beta[0]+((zt-mu)/sd)@beta[1:]), 1e-4, 1-1e-4)
