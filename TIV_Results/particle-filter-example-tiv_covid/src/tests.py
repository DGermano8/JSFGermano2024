import unittest
import scipy.stats as stats
import numpy as np

import statsmodels.stats.weightstats as smws

class TestStatsModelsDescrStatsWQuantiles(unittest.TestCase):
    """
    The descriptive stats test with weights and quantiles. Check that this works with a toy example.
    """
    def setUp(self):
        self.x = np.array([1.0, 2.0])
        self.w = np.array([2.0, 3.0])
        self.q = np.array([0.35, 0.39, 0.41, 0.45])
        self.expected_mean = 1.6
        self.expected_quantiles = np.array([1.0, 1.0, 2.0, 2.0])
        self.StatsObj = smws.DescrStatsW(self.x, weights = self.w)

    def test_mean(self):
        self.assertAlmostEqual(self.StatsObj.mean, self.expected_mean)

    def test_quantiles(self):
        self.assertListEqual(list(self.StatsObj.quantile(self.q)), list(self.expected_quantiles))


class TestScipyLognormal(unittest.TestCase):
    """
    The scipy.stats.lognorm API requires specifying mu and sigma in a
    strange way.
    """

    def setUp(self):
        self.pdf = lambda x, mu, sig: np.exp(- (np.log(x) - mu)**2 / (2 * sig**2)) / (x * sig * np.sqrt(2 * np.pi))
        self.x = np.array([1E-6, 1.0, 2.0, 5.0, 6.789])
        self.th1 = {'mu': 1.0, 'sigma': 1.0}
        self.th2 = {'mu': 2.0, 'sigma': 2.3}
        self.parameter_sets = [self.th1, self.th2]

    def test_pdf(self):
        for th in self.parameter_sets:
            scipy_vals = stats.lognorm.pdf(self.x, th['sigma'], scale = np.exp(th['mu']))
            wrong_vals = stats.lognorm.pdf(self.x, s = th['sigma'], loc = th['mu'])
            my_vals = self.pdf(self.x, th['mu'], th['sigma'])

            for ix in range(len(self.x)):
                self.assertAlmostEqual(scipy_vals[ix], my_vals[ix])
                self.assertNotEqual(wrong_vals[ix], my_vals[ix])

    def test_frozen(self):
        for th in self.parameter_sets:
            frozen_dist = stats.lognorm(th['sigma'], scale = np.exp(th['mu']))
            scipy_vals = frozen_dist.pdf(self.x)
            my_vals = self.pdf(self.x, th['mu'], th['sigma'])

            for ix in range(len(self.x)):
                self.assertAlmostEqual(scipy_vals[ix], my_vals[ix])


if __name__ == '__main__':
    unittest.main()
