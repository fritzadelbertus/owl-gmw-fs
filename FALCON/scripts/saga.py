"""
This is a light version of SAGA: https://github.com/PQShield/SAGA
"""

# Estimators for moments
# Data management
from copy import deepcopy

# Math functions
from math import ceil, exp, log, sqrt
import re

# For debugging purposes
import time

import matplotlib.pyplot as plt

# Numpy stuff
# For HZ multivariate test, used in scipy.spatial.distance.mahalanobis
# doornik hansen
from numpy import (
    array,
    corrcoef,
    cov,
    diag,
    floor,
    mean,
    power,
    set_printoptions,
    sort,
    tile,
    transpose,
)
from numpy import log as nplog
from numpy import sqrt as npsqrt
from numpy.linalg import eigh, inv, matrix_rank
import pandas

# qqplot
import scipy.stats as stats

# Statistical (normality) tests
# Distributions
from scipy.stats import chi2, chisquare, kurtosis, moment, skew

# Tailcut rate
tau = 14
# Minimal size of a bucket for the chi-squared test (must be >= 5)
chi2_bucket = 10
# Minimal p-value
pmin = 0.001
# Print options
set_printoptions(precision=4)


def gaussian(x, mu, sigma):
    """
    Gaussian function of center mu and "standard deviation" sigma.
    """
    return exp(-((x - mu) ** 2) / (2 * (sigma**2)))


def make_gaussian_pdt(mu, sigma):
    """
    Make the probability distribution table (PDT) of a discrete Gaussian.
    The output is a dictionary.
    """
    # The distribution is restricted to [-zmax, zmax).
    zmax = int(ceil(tau * sigma))
    pdt = dict()
    for z in range(int(floor(mu)) - zmax, int(ceil(mu)) + zmax):
        pdt[z] = gaussian(z, mu, sigma)
    gauss_sum = sum(pdt.values())
    for z in pdt:
        pdt[z] /= gauss_sum
    return pdt


class UnivariateSamples:
    """
    Class for computing statistics on univariate Gaussian samples.
    """

    def __init__(self, mu, sigma, list_samples):
        """
        Input:
        - the expected center mu of a discrete Gaussian over Z
        - the expected standard deviation sigma of a discrete Gaussian over Z
        - a list of samples defining an empiric distribution

        Output:
        - the means of the expected and empiric distributions
        - the standard deviations of the expected and empiric distributions
        - the skewness of the expected and empiric distributions
        - the kurtosis of the expected and empiric distributions
        - a chi-square test between the two distributions
        """
        zmax = int(ceil(tau * sigma))
        # Expected center standard variation.
        self.exp_mu = mu
        self.exp_sigma = sigma
        # Number of samples
        self.nsamples = len(list_samples)
        self.histogram = dict()
        self.outlier = 0
        # Initialize histogram
        start = int(floor(mu)) - zmax
        end = int(ceil(mu)) + zmax
        for z in range(start, end):
            self.histogram[z] = 0
        for z in list_samples:
            # Detect and count outliers (samples not in [-zmax, zmax))
            if z not in self.histogram:
                self.outlier += 1
            # Fill histogram according to the samples
            else:
                self.histogram[z] += 1
        # Empiric mean, variance, skewness, kurtosis and standard deviation
        self.mean = sum(list_samples) / self.nsamples
        self.variance = moment(list_samples, 2)
        self.skewness = skew(list_samples)
        self.kurtosis = kurtosis(list_samples)
        self.stdev = sqrt(self.variance)
        # Chi-square statistic and p-value
        self.chi2_stat, self.chi2_pvalue = self.chisquare()
        # Final assessment: the dataset is valid if:
        # - the chi-square p-value is higher than pmin
        # - there is no outlier
        self.is_valid = True
        self.is_valid &= self.chi2_pvalue > pmin
        self.is_valid &= self.outlier == 0

    def __repr__(self):
        """
        Print the sample statistics in a readable form.
        """
        rep = "\n"
        rep += f"Testing a Gaussian sampler with center = {self.exp_mu} and sigma = {self.exp_sigma}\n"
        rep += f"Number of samples: {self.nsamples}\n\n"
        rep += "Moments  |   Expected     Empiric\n"
        rep += "---------+----------------------\n"
        rep += f"Mean:    |   {self.exp_mu:.5f}      {self.mean:.5f}\n"
        rep += f"St. dev. |   {self.exp_sigma:.5f}      {self.stdev:.5f}\n"
        rep += f"Skewness |   {0:.5f}      {self.skewness:.5f}\n"
        rep += f"Kurtosis |   {0:.5f}      {self.kurtosis:.5f}\n"
        rep += "\n"
        rep += f"Chi-2 statistic:   {self.chi2_stat}\n"
        rep += f"Chi-2 p-value:     {self.chi2_pvalue}   (should be > {pmin})\n"
        rep += "\n"
        rep += f"How many outliers? {self.outlier}"
        rep += "\n\n"
        rep += f"Is the sample valid? {self.is_valid}"
        return rep

    def chisquare(self):
        """
        Run a chi-square test to compare the expected and empiric distributions
        """
        # We construct two histograms:
        # - the expected one (exp_histogram)
        # - the empirical one (histogram)
        histogram = deepcopy(self.histogram)
        # The chi-square test require buckets to have enough elements,
        # so we aggregate samples in the left and right tails in two buckets
        exp_histogram = make_gaussian_pdt(self.exp_mu, self.exp_sigma)
        obs = list(histogram.values())
        exp = list(exp_histogram.values())
        z = 0
        while 1:
            if z >= len(exp) - 1:
                break
            while (z < len(exp) - 1) and (exp[z] < chi2_bucket / self.nsamples):
                obs[z + 1] += obs[z]
                exp[z + 1] += exp[z]
                obs.pop(z)
                exp.pop(z)
            z += 1
        obs[-2] += obs[-1]
        exp[-2] += exp[-1]
        obs.pop(-1)
        exp.pop(-1)
        exp = [round(prob * self.nsamples) for prob in exp]
        diff = self.nsamples - sum(exp_histogram.values())
        exp_histogram[int(round(self.exp_mu))] += diff
        res = chisquare(obs, f_exp=exp)
        return res


class MultivariateSamples:
    """
    Class for computing statistics on multivariate Gaussian samples
    """

    def __init__(self, sigma, list_samples):
        """
        Input:
        - sigma: an expected standard deviation
        - list_samples: a list of (expected) multivariate samples

        Output:
        - univariates[]: a list of UnivariateSamples objects (one / coordinate)
        - covariance: an empiric covariance matrix
        - DH, AS, PO, PA: statistics and p-values for the Doornik-Hansen test
        - dc_pvalue: a p-value for our custom covariance-based test
        """
        # Parse the signatures and store them
        self.nsamples = len(list_samples)
        self.dim = len(list_samples[0])
        self.data = pandas.DataFrame(list_samples)
        # Expected center and standard deviation
        self.exp_mu = 0
        self.exp_si = sigma
        # Testing sphericity
        # For each coordinate, perform an univariate analysis
        self.univariates = [None] * self.dim
        for i in range(self.dim):
            self.univariates[i] = UnivariateSamples(0, sigma, self.data[i])
        self.nb_gaussian_coord = sum(
            (self.univariates[i].chi2_pvalue > pmin) for i in range(self.dim)
        )
        # Estimate the (normalized) covariance matrix
        self.covariance = cov(self.data.transpose()) / (self.exp_si**2)
        self.DH, self.AS, self.PO, self.PA = doornik_hansen(self.data)
        self.dc_pvalue = diagcov(self.covariance, self.nsamples)

    def __repr__(self):
        """
        Print the sample statistics in a readable form.
        """
        rep = "\n"
        rep += f"Testing a centered multivariate Gaussian of dimension = {self.dim} and sigma = {self.exp_si:.3f}\n"
        rep += f"Number of samples: {self.nsamples}\n"
        rep += "\n"
        rep += "The test checks that the data corresponds to a multivariate Gaussian, by doing the following:\n"
        rep += "1 - Print the covariance matrix (visual check). One can also plot\n"
        rep += "    the covariance matrix by using self.show_covariance()).\n"
        rep += "2 - Perform the Doornik-Hansen test of multivariate normality.\n"
        rep += f"    The p-value obtained should be > {pmin}\n"
        rep += "3 - Perform a custom test called covariance diagonals test.\n"
        rep += "4 - Run a test of univariate normality on each coordinate\n"
        rep += "\n"
        rep += f"1 - Covariance matrix ({self.dim} x {self.dim}):\n{self.covariance}\n"
        rep += "\n"
        if self.nsamples < 4 * self.dim:
            rep += "Warning: it is advised to have at least 8 times more samples than the dimension n.\n"
        rep += f"2 - P-value of Doornik-Hansen test:                {self.PO:.4f}\n"
        rep += "\n"
        rep += f"3 - P-value of covariance diagonals test:          {self.dc_pvalue:.4f}\n"
        rep += "\n"
        rep += f"4 - Gaussian coordinates (w/ st. dev. = sigma)?    {self.nb_gaussian_coord} out of {self.dim}\n"
        return rep

    def show_covariance(self):
        """
        Visual representation of the covariance matrix
        """
        plt.imshow(self.covariance, interpolation="nearest")
        plt.show()


def doornik_hansen(data):
    """
    Perform the Doornik-Hansen test
    (https://doi.org/10.1111/j.1468-0084.2008.00537.x)

    This computes and transforms multivariate variants of the skewness
    and kurtosis, then computes a chi-square statistic on the results.
    """
    data = pandas.DataFrame(data)
    data = deepcopy(data)

    n = len(data)
    p = len(data.columns)
    # R is the correlation matrix, a scaling of the covariance matrix
    # R has dimensions dim * dim
    R = corrcoef(data.transpose())
    L, V = eigh(R)
    for i in range(p):
        if L[i] <= 1e-12:
            L[i] = 0
        if L[i] > 1e-12:
            L[i] = 1 / sqrt(L[i])
    L = diag(L)

    if matrix_rank(R) < p:
        V = pandas.DataFrame(V)
        G = V.loc[:, (L != 0).any(axis=0)]
        data = data.dot(G)
        ppre = p
        p = data.size / len(data)
        raise ValueError(
            "NOTE:Due that some eigenvalue resulted zero, \
                          a new data matrix was created. Initial number \
                          of variables = ",
            ppre,
            ", were reduced to = ",
            p,
        )
        R = corrcoef(data.transpose())
        L, V = eigh(R)
        L = diag(L)

    means = [list(data.mean())] * n
    stddev = [list(data.std(ddof=0))] * n

    Z = (data - pandas.DataFrame(means)) / pandas.DataFrame(stddev)
    Zp = Z.dot(V)
    Zpp = Zp.dot(L)
    st = Zpp.dot(transpose(V))

    # skew is the multivariate skewness (dimension dim)
    # kurt is the multivariate kurtosis (dimension dim)
    skew = mean(power(st, 3), axis=0)
    kurt = mean(power(st, 4), axis=0)

    # Transform the skewness into a standard normal z1
    n2 = n * n
    b = 3 * (n2 + 27 * n - 70) * (n + 1) * (n + 3)
    b /= (n - 2) * (n + 5) * (n + 7) * (n + 9)
    w2 = -1 + sqrt(2 * (b - 1))
    d = 1 / sqrt(log(sqrt(w2)))
    y = skew * sqrt((w2 - 1) * (n + 1) * (n + 3) / (12 * (n - 2)))
    # Use numpy log/sqrt as math versions dont have array input
    z1 = d * nplog(y + npsqrt(y * y + 1))

    # Transform the kurtosis into a standard normal z2
    d = (n - 3) * (n + 1) * (n2 + 15 * n - 4)
    a = (n - 2) * (n + 5) * (n + 7) * (n2 + 27 * n - 70) / (6 * d)
    c = (n - 7) * (n + 5) * (n + 7) * (n2 + 2 * n - 5) / (6 * d)
    k = (n + 5) * (n + 7) * (n * n2 + 37 * n2 + 11 * n - 313) / (12 * d)
    al = a + (skew**2) * c
    chi = (kurt - 1 - (skew**2)) * k * 2
    z2 = (((chi / (2 * al)) ** (1 / 3)) - 1 + 1 / (9 * al)) * npsqrt(9 * al)
    kurt -= 3

    # omnibus normality statistic
    DH = z1.dot(z1.transpose()) + z2.dot(z2.transpose())
    AS = n / 6 * skew.dot(skew.transpose()) + n / 24 * kurt.dot(kurt.transpose())
    # degrees of freedom
    v = 2 * p
    # p-values
    PO = 1 - chi2.cdf(DH, v)
    PA = 1 - chi2.cdf(AS, v)

    return DH, AS, PO, PA


def diagcov(cov_mat, nsamples):
    """
    This test studies the population covariance matrix.
    Suppose it is of this form:
     ____________
    |     |     |
    |  1  |  3  |
    |_____|_____|
    |     |     |
    |     |  2  |
    |_____|_____|

    The test will first compute sums of elements on diagonals of 1, 2 or 3,
    and store them in the table diagsum of size 2 * dim:
    - First (dim / 2) lines = means of each diag. of 1 above leading diag.
    - Following (dim / 2) lines = means of each diag. of 2 above leading diag.
    - Following (dim / 2) lines = means of each diag. of 3 above leading diag.
    - Last (dim / 2) lines = means of each diag. of 3 below leading diag.

    We are making the assumption that each cell of the covariance matrix
    follows a normal distribution of variance 1 / n. Assuming independence
    of each cell in a diagonal, each diagonal sum of k elements should
    follow a normal distribution of variance k / n (hence of variance
    1 after normalization by n / k).

    We then compute the sum of the squares of all elements in diagnorm.
    If is supposed to look like a chi-square distribution
    """
    dim = len(cov_mat)
    n0 = dim // 2
    diagsum = [0] * (2 * dim)
    for i in range(1, n0):
        diagsum[i] = sum(cov_mat[j][i + j] for j in range(n0 - i))
        diagsum[i + n0] = sum(cov_mat[n0 + j][n0 + i + j] for j in range(n0 - i))
        diagsum[i + 2 * n0] = sum(cov_mat[j][n0 + i + j] for j in range(n0 - i))
        diagsum[i + 3 * n0] = sum(cov_mat[j][n0 - i + j] for j in range(n0 - i))
    # Diagnorm contains the normalized sums, which should be normal
    diagnorm = diagsum[:]
    for i in range(1, n0):
        nfactor = sqrt(nsamples / (n0 - i))
        diagnorm[i] *= nfactor
        diagnorm[i + n0] *= nfactor
        diagnorm[i + 2 * n0] *= nfactor
        diagnorm[i + 3 * n0] *= nfactor

    # Each diagnorm[i + _ * n0] should be a random normal variable
    chistat = sum(elt**2 for elt in diagnorm)
    pvalue = 1 - chi2.cdf(chistat, df=4 * (n0 - 1))
    return pvalue


def parse_multivariate_file(filename):
    """
    Parse a file containing several multivariate samples.

    Input:
    - the file name of a file containing k lines
      - each line corresponds to a multivariate sample
      - the samples are all assumed to be from the same distribution

    Output:
    - sigma: the expected standard deviation of the samples
    - data: a Python list of length k, containing all the samples
    """
    with open(filename) as f:
        sigma = 0
        data = []
        while True:
            # Parse each line
            line = f.readline()
            if not line:
                break  # EOF
            sample = re.split(", |,\n", line)
            sample = [int(elt) for elt in sample[:-1]]
            data += [sample]
            sigma += sum(elt**2 for elt in sample)
        # sigma is the expected sigma based on the samples
        sigma = sqrt(sigma / (len(data) * len(data[0])))
    return (sigma, data)


def test_sig(n=128, nb_sig=1000, perturb=False, level=0):
    """
    Test signatures output by a Python implementation of Falcon.
    This test allow to perturb the FFT by setting the rightmost node
    of the FFT tree (of the private key) to 0. One can check that, at
    least for moderate levels (0 to 4), the test will end up detecting
    (via diagcov) that the signatures output do not follow the correct
    distribution.

    Input:
    - n: the degree of the ring
    - nb_sig: number of signatures
    - perturb: if set to 1, one node in the FFT tree is set to 0
    - level: determines which node (the rightmost one at a given level)
      is set to 0
    """
    start = time.time()
    # Generate a private key
    sk = falcon.SecretKey(n)
    # Perturb the FFT tree
    if perturb is True:
        # Check that the level is less than the FFT tree depth
        assert (1 << level) < n
        u, k = sk.T_fft, n
        # Find the node
        for _ in range(level):
            u = u[2]
            k >>= 1
        # Zero-ize the node
        u[0] = [0] * k
    end = time.time()
    print(f"Took {end - start:.2f} seconds to generate the private key.")

    # Compute signatures
    message = "0"
    start = time.time()
    list_signatures = [sk.sign(message, reject=False) for _ in range(nb_sig)]
    # Strip away the nonces and concatenate the s_1's and s_2's
    list_signatures = [sig[1][0] + sig[1][1] for sig in list_signatures]
    end = time.time()
    print(f"Took {end - start:.2f} seconds to generate the samples.")
    # Perform the statistical test
    start = time.time()
    samples_data = MultivariateSamples(sk.sigma, list_signatures)
    end = time.time()
    print(f"Took {end - start:.2f} seconds to run a statistical test.")
    return sk, samples_data


#######################
# Supplementary Stuff #
#######################


def qqplot(data):
    """
    https://www.itl.nist.gov/div898/handbook/eda/section3/qqplot.htm
    """
    data = pandas.DataFrame(data)
    data = deepcopy(data)

    S = cov(data.transpose(), bias=1)
    n = len(data)
    p = len(data.columns)

    means = [list(data.mean())] * n
    difT = data - pandas.DataFrame(means)
    Dj = diag(difT.dot(inv(S)).dot(difT.transpose()))
    Y = data.dot(inv(S)).dot(data.transpose())
    Ytdiag = array(pandas.DataFrame(diag(Y.transpose())))
    Djk = -2 * Y.transpose()
    Djk += tile(Ytdiag, (1, n)).transpose()
    Djk += tile(Ytdiag, (1, n))
    Djk_quick = []
    for i in range(n):
        Djk_quick += list(Djk.values[i])

    chi2_random = chi2.rvs(p - 1, size=len(Dj))
    chi2_random = sort(chi2_random)
    r2 = stats.linregress(sort(Dj), sort(chi2_random))[2] ** 2
    plt.title("R-Squared = %0.20f" % r2, fontsize=9)
    plt.suptitle("QQ plot for Multivariate Normality", fontweight="bold", fontsize=12)

    plt.savefig("qqplot.eps", format="eps", bbox_inches="tight", pad_inches=0)
    plt.show()
