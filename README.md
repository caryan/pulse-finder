# pulse-finder

GRAPE style quantum optimal control. This code is an archived version of an implementation the GRAPE
algorithm [10.1016/j.jmr.2004.11.004](https://doi.org/10.1016/j.jmr.2004.11.004). For current implementations consider using the optimal control functionality builtin
into [QuTiP](http://qutip.org/docs/4.0.2/guide/guide-control.html) or more modern approaches that
leverage automatic differentiation frameworks such as
[GRAPE-Tensorflow](https://github.com/SchusterLab/quantum-optimal-control).

This code base was used extensively in quantum information processing demonstrations in both liquid-
and solid-state NMR and also ESR. Some features of this implementation are incoherent averaging over
uncertainty in Hamiltonian parameters for experimental robustness, conjugate gradients to speed up
convergence, and timestep derivatives to find time optimal control. Some further details can be found in
[10.1103/PhysRevA.78.012328](https://doi.org/10.1103/PhysRevA.78.012328).

This code was used in the following publications:

1. Moussa, O., da Silva, M.P., Ryan, C.A. and Laflamme, R., 2012. Practical experimental certification of computational quantum gates using a twirling procedure. Physical review letters, 109(7), 070504. [10.1103/PhysRevLett.109.070504](https://doi.org/10.1103/PhysRevLett.109.070504)
2. Zhang, Y., Ryan, C.A., Laflamme, R. and Baugh, J., 2011. Coherent control of two nuclear spins using the anisotropic hyperfine interaction. Physical review letters, 107(17), 170503.
[doi.org/10.1103/PhysRevLett.107.170503](https://doi.org/10.1103/PhysRevLett.107.170503)
3. Moussa, O., Baugh, J., Ryan, C.A. and Laflamme, R., 2011. Demonstration of sufficient control for two rounds of quantum error correction in a solid state ensemble quantum information processor. Physical review letters, 107(16), 160501. [10.1103/PhysRevLett.107.160501](https://doi.org/10.1103/PhysRevLett.107.160501)
4. Souza, A.M., Zhang, J., Ryan, C.A. and Laflamme, R., 2011. Experimental magic state distillation for fault-tolerant quantum computing. Nature communications, 2, 169. [10.1038/ncomms1166](https://doi.org/10.1038/ncomms1166)
5. Moussa, O., Ryan, C.A., Cory, D.G. and Laflamme, R., 2010. Testing contextuality on quantum ensembles with one clean qubit. Physical review letters, 104(16), 160501. [10.1103/PhysRevLett.104.160501](https://doi.org/10.1103/PhysRevLett.104.160501)
6. Ryan, C.A., Laforest, M. and Laflamme, R., 2009. Randomized benchmarking of single-and multi-qubit control in liquid-state NMR quantum information processing. New Journal of Physics, 11(1), 013034. [0.1088/1367-2630/11/1/013034](https://doi.org/10.1088/1367-2630/11/1/013034)
7. Ryan, C.A., Negrevergne, C., Laforest, M., Knill, E. and Laflamme, R., 2008. Liquid-state nuclear magnetic resonance as a testbed for developing quantum control methods. Physical Review A, 78(1), 012328. [10.1103/PhysRevA.78.012328](https://doi.org/10.1103/PhysRevA.78.012328)
8. Ryan, C.A., Moussa, O., Baugh, J. and Laflamme, R., 2008. Spin based heat engine: demonstration of multiple rounds of algorithmic cooling. Physical review letters, 100(14), 140501. [10.1103/PhysRevLett.100.140501](https://doi.org/10.1103/PhysRevLett.100.140501)
9. Emerson, J., Silva, M., Moussa, O., Ryan, C., Laforest, M., Baugh, J., Cory, D.G. and Laflamme, R., 2007. Symmetrized characterization of noisy quantum processes. Science, 317(5846), pp.1893-1896. [10.1126/science.1145699](https://doi.org/10.1126/science.1145699)
