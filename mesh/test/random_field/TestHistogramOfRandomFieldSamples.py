import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


# Change this file name to the name of your random field (relative to $CHASTE_TEST_OUTPUT)
# This should match the name you used in TestSamplesFromCachedRandomField::TestOuputSamplesToFile()
file_name = 'CachedRandomFields/xy_0.000_1.000_2.000_5.000_10_12_0_0_0.800_0.500.rfg'

grf = np.genfromtxt(os.path.join(os.environ['CHASTE_TEST_OUTPUT'], file_name + '.check'), delimiter=',')
print(grf)
plt.figure(1, figsize=(16, 10), dpi=72)
plt.subplot(1, 1, 1)

mu, sigma = 0, 1

plt.axis([-3*sigma, 3*sigma, 0, 0.5])
plt.grid(True)
n, bins, patches = plt.hist(grf, bins=50, density=True, facecolor='green', alpha=0.75)

# add a 'best fit' line
y1 = norm.pdf(bins, mu, sigma)
l1 = plt.plot(bins, y1, 'b--', linewidth=5)

plt.show()