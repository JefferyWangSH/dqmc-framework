# DQMC-Hubbard Model [![forthebadge](https://forthebadge.com/images/badges/works-on-my-machine.svg)](https://forthebadge.com)

In this repository, we present a C++ implementation of determinant Quantum Monte Carlo `(DQMC)` simulation applied to the half-filled fermion Hubbard model.

Some preliminary supports for doped cases are also realized by simple reweighting of configuration wights.

---

## Installation ##

### Prerequisite ###

* `gcc/g++` `( version > 4.8.5 )` and `cmake` `( version > 2.8.12 )` installed.
* `Boost C++ libraries` `( version > 1.71 )` installed.
* `Eigen library` `( version > 3.4.0 )` providing a user-friendly interface of matrices.
* `Intel Math Kernel Library (MKL)` for high-accuracy linear algebra and numerical stabilization.
* `Intel implementation of Message Passing Interface (Intel MPI)` for large range of distributed parallel.

### Usage ###

1. Download source codes from github.
    ``` shell
    $ # download the source code
    $ git clone https://github.com/JefferyWangSH/DQMC-Hubbard-Model.git {PROGRAM_ROOT}
    ```
2. Enter the [`build/`](build/) directory and run [`runcmake.sh`](build/runcmake.sh) which will analyse the program using cmake.
    ``` shell
    $ # initialize cmake
    $ cd {PROGRAM_ROOT}/build && ./runcmake.sh
    ```
3. Enter the [`run/`](run/) directory and compile the codes. 
    ``` shell
    $ # build the program
    $ cd {PROGRAM_ROOT}/run && ./make.sh
    ```
4. Run the script [`run/batch.sh`](run/batch.sh) to start the simulation if the program is successfully build. ( We use `Slurm` for the management of program tasks, hence the simulation parameters should be edited in the script [`run/run.sh`](run/run.sh). )
    ```shell
    $ # start simulation
    $ ./batch.sh

    $ # edit the simulation params in `run.sh`
    $ vim run.sh
    ```
5. One can always use command line option `--help` to see helping messages:
    ``` shell
    $ # see helping messages
    $ mpirun {PROGRAM_ROOT}/build/dqmc_hubbard --help
    ```


## Features ##

1. Program works well in condition of moderate model parameters.
2. Stable sweeping procedure with high computational efficiency.
3. Support distributed parallel using MPI.  
4. Support both attractive and repulsive hubbard interaction.
5. Support simulation of doped systems by reweighting.
6. Support equal-time and dynamical measurements of physical observables.
7. Support checkerboard break-ups with efficient linear algebra.


## References ##

1. H. Fehske, R. Schneider, A. Weiße, Computational Many-Particle Physics, *Springer*, 2008. [doi](https://doi.org/10.1007/978-3-540-74686-7)
2. James Gubernatis, Naoki Kawashima, Philipp Werner, Quantum Monte Carlo Methods: Algorithms for Lattice Models, *Cambridge University Press*, 2016. [doi](https://doi.org/10.1017/CBO9780511902581)
3. Xiao Yan Xu, *Tutorial*: Solving square lattice Hubbard Model with DQMC, 2016. [doi](http://ziyangmeng.iphy.ac.cn/files/teaching/SummerSchoolSimpleDQMCnoteXYX201608.pdf)
4. C. N. Varney, C. R. Lee, Z. J. Bai et al., Quantum Monte Carlo study of the two-dimensional fermion Hubbard model, *Phys. Rev. B 80, 075116*, 2009. [doi](https://doi.org/10.1103/PhysRevB.80.075116)
5. Thereza Paiva, Raimundo R. dos Santos, R. T. Scalettar et al., Critical temperature for the two-dimensional attractive Hubbard model, *Phys. Rev. B 69, 184501*, 2004. [doi](https://doi.org/10.1103/PhysRevB.69.184501)
6. Xiao Yan Xu, Kai Sun, Yoni Schattner et al., Non-Fermi Liquid at (2 + 1) D Ferromagnetic Quantum Critical Point, *Phys. Rev. X 7, 031058*. [doi](https://doi.org/10.1103/PhysRevX.7.031058)
7. Douglas J. Scalapino, Steven R. White, and Shoucheng Zhang, Insulator, metal, or superconductor: The criteria, *Phys. Rev. B 47, 7995.* [doi](https://doi.org/10.1103/PhysRevB.47.7995)


## License & Support ##

The code is open source under GPL-3.0 License. 

If any questions, feel free to contact me via email 17307110117@fudan.edu.cn.