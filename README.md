DQMC-Hubbard Model
==================

[![forthebadge](https://forthebadge.com/images/badges/works-on-my-machine.svg)](https://forthebadge.com)

C++ implementation of determinant Quantum Monte Carlo `(DQMC)` simulation applied to half-filled Hubbard model.

Prerequisite
------------
1. g++ and cmake installed.
2. Eigen library, linked to Intel MKL.
3. BLAS and Lapack for SVD stabilization.

Accomplished
------------
1. Cyclic update (sweep back and forth).
2. Stable sweeping process with high computational efficiency.
3. Support equal-time and time-displaced measurements of physical quantities.
4. Support attractive and repulsive hubbard interaction.
5. Program works well in condition of moderate model parameters. 

TODO
----
1. Extract fermion spectrum function by numerical analytic continuation (e.g. Stochastic Analytic Continuation, `SAC`).
2. Determine the critical temperature of superconducting transition in negative-U cases.

References
----------
1. James Gubernatis, Naoki Kawashima, Philipp Werner  
   Quantum Monte Carlo Methods: Algorithms for Lattice Models  
   Cambridge University Press, 2016. [DOI](https://doi.org/10.1017/CBO9780511902581)
2. H. FehskeR. SchneiderA. Weiße  
   Computational Many-Particle Physics  
   Springer, 2008. [DOI](https://doi.org/10.1007/978-3-540-74686-7)
3. Douglas J. Scalapino, Steven R. White, and Shoucheng Zhang  
   Insulator, metal, or superconductor: The criteria  
   Phys. Rev. B 47, 7995. [DOI](https://doi.org/10.1103/PhysRevB.47.7995)
4. Xiao Yan Xu, Kai Sun, Yoni Schattner, Erez Berg, and Zi Yang Meng  
   Non-Fermi Liquid at (2 + 1) D Ferromagnetic Quantum Critical Point  
   Phys. Rev. X 7, 031058. [DOI](https://doi.org/10.1103/PhysRevX.7.031058)
   