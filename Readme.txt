
#This package demonstrates the methods developed in Laha et al. (2020) for shape-constrained tests of stochastic
#dominance and shape-cnstrained estimation of the Hellinger-distance.
# The R function for the shape-constrained tests: SDNN
# Estimator of Hellinger distance: hd.uni (unimodal density),
# hd.lc (log-concave density), hd.lc.sm (smoothed log-concave density).
# Use the R function hell.ci to compute the corresponding Weld-Type confidence intervals.
# See the manual SDNNtests_0.0.0.9000.pdf for more details.

#To install the package SDNNtests type in your R console the following:
library(devtools)
install_github('nilanjanalaha/SDNNtests')

#To search for the R function SDNN, try
??SDNN

#To search for the R function hd.uni, try
??hd.uni

