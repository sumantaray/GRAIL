"""Population metrics. Inputs to WD and SWD are size-matched by the caller."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist, jensenshannon, pdist
from scipy.stats import false_discovery_control, mannwhitneyu, spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score


def correlation(x, y):
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(spearmanr(x, y).statistic)


def wasserstein(x, y):
    if len(x) != len(y) or len(x) == 0:
        raise ValueError("Exact uniform W1 requires equal nonempty samples in this implementation.")
    cost = cdist(x, y, metric="euclidean")
    row, col = linear_sum_assignment(cost)
    return float(cost[row, col].mean())


def projected_w1(x, y, directions):
    if len(x) != len(y):
        raise ValueError("SWD diagnostic expects matched population sizes.")
    return float(np.abs(np.sort(x @ directions, axis=0) - np.sort(y @ directions, axis=0)).mean())


def fit_markers(a, b, per_direction):
    delta = b.mean(0) - a.mean(0)
    p = mannwhitneyu(a, b, axis=0, method="asymptotic").pvalue
    q = false_discovery_control(np.nan_to_num(p, nan=1.0))
    result = []
    for sign in (1, -1):
        keep = np.flatnonzero(sign * delta > 0)
        order = np.lexsort((keep, -sign * delta[keep], q[keep]))
        result.extend(keep[order[:per_direction]])
    return np.asarray(result, dtype=int)


def marker_trend(observed, generated, markers):
    obs = np.array([x.mean(axis=0)[markers] for x in observed])
    pred = np.array([x.mean(axis=0)[markers] for x in generated])
    scores = np.array([correlation(x, y) for x, y in zip(obs.T, pred.T)])
    valid = np.isfinite(scores)
    return (float(scores[valid].mean()) if valid.any() else np.nan), int(valid.sum())


def reconstruction(ae, anchors, groups, seed):
    from .model import encode
    if min(anchors.shape) < 64:
        raise ValueError("The reconstruction audit needs at least 64 training anchors and 64 genes.")
    pca = PCA(n_components=64, svd_solver="full", random_state=seed).fit(anchors)
    rows = []
    for stage, x in groups.items():
        predictions = {"mean": np.broadcast_to(anchors.mean(0), x.shape),
                       "pca64": pca.inverse_transform(pca.transform(x)),
                       "autoencoder": encode(ae.decoder, encode(ae.encoder, x))}
        for method, pred in predictions.items():
            error = (pred - x) ** 2
            rows.append(dict(stage=stage, method=method, n=len(x),
                             gene_rmse=np.median(np.sqrt(error.mean(0))),
                             cell_rmse=np.median(np.sqrt(error.mean(1)))))
    return rows


def geometry(projected, grid):
    centers = np.array([projected[t].mean(0) for t in grid], dtype=np.float64)
    steps = np.diff(centers, axis=0)
    norms = np.linalg.norm(steps, axis=1)
    denom = norms[:-1] * norms[1:]
    valid = denom > 1e-12
    curv = np.linalg.norm(np.diff(steps, axis=0), axis=1)[valid] / denom[valid]
    return dict(TL=norms.sum(), chord=np.linalg.norm(centers[-1] - centers[0]),
                AC=curv.mean() if len(curv) else np.nan, valid_curvature_terms=int(valid.sum()))


def population(real, generated, proj_real, proj_gen, seed, types=None, annotator=None):
    apd = pdist(proj_real).mean()
    tri = np.triu_indices(proj_real.shape[1], k=1)
    row = dict(hvg_spearman=correlation(real.mean(0), generated.mean(0)),
               diversity_ratio=pdist(proj_gen).mean() / apd if apd > 0 else np.nan,
               covariance_spearman=correlation(np.cov(proj_real.T)[tri], np.cov(proj_gen.T)[tri]))
    if len(real) >= 5:
        x, y = np.vstack([real, generated]), np.repeat([0, 1], len(real))
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        models = {"lr": LogisticRegression(solver="liblinear", max_iter=1000, random_state=seed),
                  "rf": RandomForestClassifier(n_estimators=100, random_state=seed, n_jobs=1)}
        for name, model in models.items():
            row[f"auroc_{name}"] = cross_val_score(model, x, y, cv=cv, scoring="roc_auc").mean()
    if annotator is not None and types is not None:
        predicted = annotator.predict(generated)
        labels = np.union1d(types, np.union1d(annotator.classes_, predicted))
        real_p = [(types == name).mean() for name in labels]
        gen_p = [(predicted == name).mean() for name in labels]
        row["cell_type_jsd_base2"] = jensenshannon(real_p, gen_p, base=2) ** 2
    return row
