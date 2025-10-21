This is an arbitrary masked version of the keccakf-1600 permutation. It contains the 5 methods, theta, rho, pi, chi and iota. Each of these methods has been designed for
side-channel security and has been evaluated at first order masking with TVLA order 1-4 for the linnear operations theta, rho, pi and iota and at masking orders 1-3 with 
TVLA orders 1-4 for the only non-linnear method chi. This is not a claim of security, simply that it is expected not to leak information with a T value greater than 4.5.

1st order evaluation 50k vs 50k fixed vs random
2nd order evaluation 200k vs 200k fixed vs random
3rd order evaluation 500k vs 500k fixed vs random
