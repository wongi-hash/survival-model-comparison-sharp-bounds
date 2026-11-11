"""Helpers for the Rotterdam -> GBSG application (predictors fitted with lifelines)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

TAU = 365.25 * 5
CAP = 365.25 * 7
SEEDS = [42, 52, 62, 72]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def source_outcome(raw):
    time = np.where(raw["recur"].to_numpy() == 1, raw["rtime"].to_numpy(), raw["dtime"].to_numpy()).astype(float)
    event = ((raw["recur"].to_numpy() == 1) | (raw["death"].to_numpy() == 1)).astype(int)
    event = np.where(time <= CAP, event, 0)
    time = np.minimum(time, CAP)
    return time, event


def target_outcome(raw):
    time = np.minimum(raw["rfstime"].to_numpy(dtype=float), CAP)
    event = np.where(raw["rfstime"].to_numpy(dtype=float) <= CAP, raw["status"].to_numpy(dtype=int), 0)
    return time, event


def feature_frame(raw, cohort):
    age10 = raw["age"].to_numpy(dtype=float) / 10.0
    if cohort == "rotterdam":
        size = raw["size"].astype(str)
        size_mid = size.eq("20-50").astype(float).to_numpy()
        size_large = size.eq(">50").astype(float).to_numpy()
    elif cohort == "gbsg":
        size = raw["size"].to_numpy(dtype=float)
        size_mid = ((size > 20) & (size <= 50)).astype(float)
        size_large = (size > 50).astype(float)
    else:
        raise ValueError(cohort)
    nodes = raw["nodes"].to_numpy(dtype=float)
    if np.any(nodes <= 0):
        raise ValueError("feature_frame requires node-positive rows")
    return pd.DataFrame({
        "age3": age10 ** 3,
        "age3log": (age10 ** 3) * np.log(age10),
        "meno": raw["meno"].to_numpy(dtype=float),
        "size_mid": size_mid,
        "size_large": size_large,
        "nodes_inv_sqrt": 1.0 / np.sqrt(nodes),
        "hormon": raw["hormon"].to_numpy(dtype=float),
        "er": raw["er"].to_numpy(dtype=float),
    })


def fit_scaler(frame, cols):
    mu = frame[cols].mean(axis=0)
    sd = frame[cols].std(axis=0, ddof=0).replace(0, 1.0)
    return mu, sd


def apply_scaler(frame, cols, mu, sd):
    return (frame[cols] - mu) / sd


def fit_predictor(train_features, train_time, train_event, cols, penalizer=0.01):
    mu, sd = fit_scaler(train_features, cols)
    x = apply_scaler(train_features, cols, mu, sd).copy()
    x["time"] = train_time
    x["event"] = train_event
    model = CoxPHFitter(penalizer=penalizer)
    model.fit(x, duration_col="time", event_col="event", show_progress=False)
    return {"model": model, "cols": list(cols), "mu": mu, "sd": sd}


def predict_survival(lock, features, tau=TAU):
    x = apply_scaler(features, lock["cols"], lock["mu"], lock["sd"])
    pred = lock["model"].predict_survival_function(x, times=[tau]).iloc[0].to_numpy(dtype=float)
    return np.clip(pred, 1e-4, 1 - 1e-4)


def serialize_lock(lock):
    model = lock["model"]
    return {
        "columns": lock["cols"],
        "mean": {k: float(v) for k, v in lock["mu"].items()},
        "sd": {k: float(v) for k, v in lock["sd"].items()},
        "coef": {k: float(v) for k, v in model.params_.items()},
        "penalizer": float(model.penalizer),
        "baseline_survival_tau": float(model.baseline_survival_.reindex(model.baseline_survival_.index.union([TAU])).sort_index().ffill().loc[TAU].iloc[0]),
    }
