"""The fixed-pair, varying-weight worked example from Table 3."""
import numpy as np
from scipy.stats import spearmanr


def calculate():
    t = np.array([0.25, 0.50, 0.75, 1.0])
    a = np.array([[1., 8., 3.], [3., 4., 5.]])
    b = np.array([[5., 4., 6.], [7., 0., 4.]])
    w = np.array([0.25, 0.75])
    pairs = np.array([(1 - s) * a + s * b for s in t])
    fixed = w @ pairs
    delta = w @ (b - a)
    output = {}
    for g, gene in enumerate(["increasing", "decreasing", "flat"]):
        score = spearmanr(fixed[:, g], np.sign(delta[g])*t).statistic if delta[g] else np.nan
        output[gene] = (fixed[:,g], score)
    changing = np.array([[0.10,0.90],[0.90,0.10]] * 2)
    values = (changing * pairs[:,:,0]).sum(axis=1)
    output["varying weights, increasing"] = values, spearmanr(values, t).statistic
    return output


if __name__ == "__main__":
    for label, (values, dc) in calculate().items():
        print(label, values.tolist(), "excluded" if np.isnan(dc) else f"DC={dc:.3f}")
