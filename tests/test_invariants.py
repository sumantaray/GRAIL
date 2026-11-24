"""Regression checks for the reviewer's fixed-weight, leakage and inference questions."""
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import torch

from grail.data import feature_space, select_anchors
from grail.metrics import wasserstein, projected_w1
from grail.model import Generator, generate, pair_codes
from scripts.directional_check import calculate


class RevisionChecks(unittest.TestCase):
    def test_affine_example(self):
        result = calculate()
        self.assertAlmostEqual(result["increasing"][1], 1.)
        self.assertAlmostEqual(result["decreasing"][1], 1.)
        self.assertTrue(np.isnan(result["flat"][1]))
        self.assertAlmostEqual(result["varying weights, increasing"][1], .6)

    def test_intermediates_cannot_select_features(self):
        rng = np.random.default_rng(42)
        t = np.repeat([0., .5, 1.], 20)
        raw = rng.poisson(3, (60,100)).astype("float32")
        data = ad.AnnData(sparse.csr_matrix(raw), obs=pd.DataFrame({"t": t}, index=[f"c{i}" for i in range(60)]))
        cfg = dict(min_genes=4,min_cells=3,hvg=66,target_sum=10000,test_fraction=.2)
        x, train, test, _, spec = feature_space(data,cfg,42)
        # Perturb both kinds of held-out cells, preserving the endpoint labels.
        raw[t==.5] = rng.poisson(64,raw[t==.5].shape)
        raw[test] = rng.poisson(99,raw[test].shape)
        data.X = sparse.csr_matrix(raw)
        xx, other_train, other_test, _, second = feature_space(data,cfg,42)
        self.assertEqual(spec,second)
        np.testing.assert_array_equal(x[train],xx[other_train])
        np.testing.assert_array_equal(test,other_test)
        self.assertFalse(set(train)&set(test))
        self.assertTrue(np.all(t[train]!=.5))

    def test_pairing_and_deterministic_readout(self):
        torch.manual_seed(42)
        rng = np.random.default_rng(42)
        x = rng.normal(size=(80,64)).astype("float32")
        ix = select_anchors(x,40,rng)
        self.assertEqual(len(set(ix)),40)
        z, ids = pair_codes(x[:4], x[4:],"nearest",rng)
        np.testing.assert_array_equal(z,x[4:][ids])
        model = Generator(20).eval()
        actual = generate(model,x)
        with torch.no_grad():
            expected = model(torch.from_numpy(x),torch.zeros_like(torch.from_numpy(x))).numpy()
        np.testing.assert_array_equal(actual,expected)
        np.testing.assert_array_equal(actual,generate(model,x))

    def test_known_distances(self):
        x = np.array([[0.],[3.]])
        self.assertAlmostEqual(wasserstein(x, x[::-1]+2),2.)
        self.assertAlmostEqual(projected_w1(x, x+2, np.ones((1,3))),2.)


if __name__ == "__main__":
    unittest.main()
