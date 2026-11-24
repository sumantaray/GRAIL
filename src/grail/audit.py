"""Evaluate a completed seed with its training-only feature map and frozen networks."""
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

from . import metrics
from .model import encode, generate, pair_codes


def evaluate(data, X, t, train, test, qc, ia, ib, ae, gen, cfg, seed, variant, external=None):
    rng = np.random.default_rng(seed)
    endpoints = {s: test[t[test] == s] for s in (0.0, 1.0)}
    a, b = X[ia], X[ib]
    anchors = np.vstack([a, b])
    za, zb = encode(ae.encoder, a), encode(ae.encoder, b)
    markers = metrics.fit_markers(X[train[t[train] == 0]], X[train[t[train] == 1]], cfg["markers_per_direction"])
    stages = sorted(set(t[qc]) - {0., 1.})
    grid = sorted(set(cfg["grid"]) | set(stages) | {0., 1.})
    rows = {"reconstruction": metrics.reconstruction(ae, anchors, {s: X[i] for s,i in endpoints.items()}, seed)}
    all_preds, all_latent, paired_ids = {}, {}, {}
    for strategy in ("nearest", "random", "centroid"):
        dest, partners = pair_codes(za, zb, strategy, np.random.default_rng(seed + 1))
        paths = {s: (1 - s) * za + s * dest for s in grid}
        all_latent[strategy] = paths
        all_preds[strategy] = {s: generate(gen, z) for s,z in paths.items()}
        if partners is not None:
            paired_ids[strategy] = data.obs_names[ib[partners]].tolist()
        if strategy == "nearest":
            all_preds["linear_expression"] = {s: (1 - s) * a + s * b[partners] for s in grid}
            all_preds["ae_decoder"] = {s: encode(ae.decoder, z) for s,z in paths.items()}
    if external:
        for name, pops in external.items():
            if name in all_preds or not set(grid).issubset(pops):
                raise ValueError(f"External method {name} needs a unique name and all grid coordinates {grid}.")
            all_preds[name] = pops
    pca = PCA(n_components=min(10, X.shape[1], len(train)), random_state=seed).fit(X[train])
    zmean, zscale = np.vstack([za, zb]).mean(0), np.vstack([za, zb]).std(0)
    zscale[zscale == 0] = 1
    directions = rng.normal(size=(za.shape[1], 256))
    directions /= np.linalg.norm(directions, axis=0)
    clf = LogisticRegression(max_iter=1000, solver="liblinear", random_state=seed).fit(X[train], t[train])
    annotation = None
    if "cell_type" in data.obs:
        ytrain = data.obs["cell_type"].iloc[train]
        valid = ytrain.notna().to_numpy()
        if ytrain[valid].nunique() > 1:
            annotation = LogisticRegression(max_iter=1000, random_state=seed).fit(X[train[valid]], ytrain[valid].astype(str))
    rows.update({key: [] for key in ("wd", "population", "sensitivity", "latent", "endpoint_wd", "mgtc", "geometry", "ccm", "branch")})
    cap = cfg["evaluation_cells"]
    # Stage controls use one set of pairs and one real sample at each stage.
    limits = [len(ia), cap] + [len(v) for v in endpoints.values()]
    limits += [int(np.sum(qc & (t == s))) for s in stages]
    for method, pops in (external or {}).items():
        limits.extend(len(pops[s]) for s in grid)
    n = min(limits)
    if n < 2:
        raise ValueError("Evaluation requires at least two eligible cells per compared population.")
    pairs = rng.choice(len(ia), n, replace=False)
    reference = {s: rng.choice(ids, n, replace=False) for s,ids in endpoints.items()}
    reference.update({s: rng.choice(np.flatnonzero(qc & (t == s)), n, replace=False) for s in stages})
    real = {s: X[ids] for s,ids in reference.items()}
    projected = {s: pca.transform(v) for s,v in real.items()}
    gen_proj = {}
    for method, pops in all_preds.items():
        matched = {s: (pop[rng.choice(len(pop), n, replace=False)] if external and method in external else pop[pairs])
                   for s,pop in pops.items()}
        gen_proj[method] = {s: pca.transform(pop) for s,pop in matched.items()}
        if len(stages) >= 2:
            score, valid = metrics.marker_trend([real[s] for s in stages], [matched[s] for s in stages], markers)
            rows["mgtc"].append(dict(method=method, mgtc=score, valid_markers=valid))
        for s in grid:
            rows["ccm"].append(dict(method=method, t=s, mean_p_stage_b=clf.predict_proba(pops[s])[:,1].mean()))
        for s in stages:
            rows["wd"].append(dict(method=method, t=s, n=n, wd=metrics.wasserstein(gen_proj[method][s], projected[s])))
            if method == "nearest" or method in (external or {}):
                labels = None
                if annotation is not None:
                    y = data.obs["cell_type"].iloc[reference[s]]
                    if y.notna().all():
                        labels = y.astype(str).to_numpy()
                scores = metrics.population(real[s], matched[s], projected[s], gen_proj[method][s], seed, labels, annotation)
                record = dict(t=s, n=n, **scores)
                if external:
                    record["method"] = method
                rows["population"].append(record)
    zpath = all_latent["nearest"]
    ref_projected = {s: (zpath[s][pairs] - zmean) / zscale for s in grid}
    gp = gen_proj["nearest"]
    def best_other(values, s, distance, target):
        dist = [distance(values[v], target) for v in stages if v != s]
        return min(dist) if dist else np.nan
    for s in stages:
        target = projected[s]
        fakes = all_preds["nearest"]
        rows["sensitivity"].append(dict(t=s, n=n,
            paired_rms=np.sqrt(np.square(fakes[s][pairs] - fakes[0][pairs]).mean(1)).mean(),
            wd_matched=metrics.wasserstein(gp[s], target), wd_fixed_a=metrics.wasserstein(gp[0], target),
            wd_other=best_other(gp, s, metrics.wasserstein, target), wd_real_b=metrics.wasserstein(projected[1], target)))
        zs = (encode(ae.encoder, real[s]) - zmean) / zscale
        def swd(x, y):
            return metrics.projected_w1(x, y, directions)
        rows["latent"].append(dict(t=s, n=n, swd_matched=swd(ref_projected[s], zs), swd_a=swd(ref_projected[0], zs),
                                   swd_b=swd(ref_projected[1], zs), swd_other=best_other(ref_projected, s, swd, zs)))
    for s in endpoints:
        rows["endpoint_wd"].append(dict(reference_t=s, n=n, wd_generated_t0=metrics.wasserstein(gp[0], projected[s])))
    geom_grid = sorted(cfg["geometry_grid"])
    common = [X[train]]
    for method in all_preds.values():
        for s in geom_grid:
            vals = method[s]
            common.append(vals[rng.choice(len(vals), min(16, len(vals)), replace=False)])
    geom_map = PCA(n_components=pca.n_components_, random_state=seed).fit(np.vstack(common))
    for method, pops in all_preds.items():
        coords = {s: geom_map.transform(pops[s]) for s in geom_grid}
        rows["geometry"].append(dict(method=method, **metrics.geometry(coords, geom_grid)))
    if "branch" in data.obs:
        terminal = train[t[train] == 1]
        y = data.obs["branch"].iloc[terminal]
        valid = y.notna().to_numpy()
        if y[valid].nunique() >= 2:
            bc = LogisticRegression(max_iter=1000, random_state=seed).fit(X[terminal[valid]], y[valid].astype(str))
            all_b = data.obs.loc[qc & (t == 1), "branch"].dropna().astype(str).to_numpy()
            pred = bc.predict(all_preds["nearest"][1])
            for k in bc.classes_:
                rows["branch"].append(dict(branch=k, reference_count=(all_b == k).sum(), reference_n=len(all_b),
                                           reference_percent=100 * (all_b == k).mean(), generated_count=(pred == k).sum(),
                                           generated_n=len(pred), generated_percent=100 * (pred == k).mean()))
    audit_ids = dict(reference_ids={str(s): data.obs_names[ids].tolist() for s,ids in reference.items()},
                     paired_sample_ids=data.obs_names[ia[pairs]].tolist(), paired_b_ids=paired_ids,
                     marker_indices=markers.tolist(), variant=variant)
    return rows, all_preds["nearest"], audit_ids
