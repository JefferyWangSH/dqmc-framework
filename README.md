# DQMC-Hubbard Model [![forthebadge](https://forthebadge.com/images/badges/works-on-my-machine.svg)](https://forthebadge.com)

C++ implementation of determinant Quantum Monte Carlo `(DQMC)` simulation applied to half-filled fermion Hubbard model is presented in this repository. 

Some preliminary supports for doped cases are also realized by simple reweighing of configuration wights.

---

## Installation ##

### Prerequisite ###

* gcc/g++ and cmake ( version > 3.21.1 ) installed.
* Boost C++ libraries ( version > 1.71 ) installed.
* Eigen library ( version > 3.3.7 ), linked to Intel Math Kernel Library ( MKL ).
* `mkl_lapacke.h` for high-accuracy numerical stabilization using SVD decomposition.

### Usage ###

1. Download source codes from github: 
    ``` shell
    $ git clone https://github.com/JefferyWangSH/DQMC-Hubbard-Model.git dqmc-hubbard
    ```
2. Creat empty directory `build` for program compilation:
    ``` shell
    $ mkdir ./dqmc-hubbard/build && cd ./dqmc-hubbard/build
    ```
3. Cmake and make:
    ``` shell
    $ cmake .. && make
    ```
4. Run the code and using command line option `--help` to see helping messages:
    ``` shell
    $ ./dqmc_hubbard --help
    ```


## Features ##

1. Program works well in condition of moderate model parameters.
2. Cyclic update ( sweep back and forth ).
3. Stable sweeping procedure with high computational efficiency.
4. Support checkerboard break-ups with efficient linear algebra.
5. Support attractive and repulsive hubbard interaction.
6. Support slightly doping away from half-filled case. 
7. Support equal-time and time-displaced measurements of physical observables.


## References ##

1. H. Fehske, R. Schneider, A. Weiße, Computational Many-Particle Physics, *Springer*, 2008. [DOI](https://doi.org/10.1007/978-3-540-74686-7)
2. James Gubernatis, Naoki Kawashima, Philipp Werner, Quantum Monte Carlo Methods: Algorithms for Lattice Models, *Cambridge University Press*, 2016. [DOI](https://doi.org/10.1017/CBO9780511902581)
3. Xiao Yan Xu, *Tutorial*: Solving square lattice Hubbard Model with DQMC, 2016. [DOI](http://ziyangmeng.iphy.ac.cn/files/teaching/SummerSchoolSimpleDQMCnoteXYX201608.pdf)
4. C. N. Varney, C. R. Lee, Z. J. Bai et al., Quantum Monte Carlo study of the two-dimensional fermion Hubbard model, *Phys. Rev. B 80, 075116*, 2009. [DOI](https://doi.org/10.1103/PhysRevB.80.075116)
5. Thereza Paiva, Raimundo R. dos Santos, R. T. Scalettar et al., Critical temperature for the two-dimensional attractive Hubbard model, *Phys. Rev. B 69, 184501*, 2004. [DOI](https://doi.org/10.1103/PhysRevB.69.184501)
6. Xiao Yan Xu, Kai Sun, Yoni Schattner et al., Non-Fermi Liquid at (2 + 1) D Ferromagnetic Quantum Critical Point, *Phys. Rev. X 7, 031058*. [DOI](https://doi.org/10.1103/PhysRevX.7.031058)
7. Douglas J. Scalapino, Steven R. White, and Shoucheng Zhang, Insulator, metal, or superconductor: The criteria, *Phys. Rev. B 47, 7995.* [DOI](https://doi.org/10.1103/PhysRevB.47.7995)


## License & Support ##

The code is open source under GPL-3.0 License. 

If any questions, feel free to contact me via email 17307110117@fudan.edu.cn.