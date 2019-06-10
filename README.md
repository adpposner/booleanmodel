# booleanmodel
Source for Boolean model (direct solution &amp; simulation) for associated manuscript

The model used consists of two parts - the solver subroutines, written in Fortran 2008 (only because of the use of `iso_c_binding`), and the interface, written in C99 with a quad-precision extension. They have been optimized for use with the Intel compiler collection (icc/ifort), all 2018-2019 versions, and Intel MKL 2018-2019 versions, but have been built successfully with gcc/gfortran as well.

These have been configured to work on both mainline Linux flavors and Mac OS. Provided that one has Intel MKL set up and has suitable compilers, compilation and execution are fairly straightforward.

One is strongly encouraged to use the [Intel Link Line Advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor) to configure the appropriate include and link flags.

## Solver subroutines
Located in the *fort* directory. The makefile should be configured with appropriate entries for `FC` and `AR`. The values for `FINCLUDEFLAGS` and `MKLFLAGS` may be acceptable for default MKL installations. All object files are built in the `fbuild` subdirectory, while all module files are built in the `fmod` subdirectories by default. Additionally, they are composited by default into `fortlibs.a`; this is entirely optional.

The values of `FFLAGS` are currently ifort-specific, but can all be removed to run w/gcc. Note the use of the 128-bit quad precision kind parameter qp - this is both ifort- and gfortran-compatible. These items aside, building these items should be as simple as calling `make` in the `fort` subdirectory, once appropriate changes have been made. 

Building within the fort subdirectory generates a small demo binary called **fgensln**, which can be used to check the basic routines. Parameters can be entered manually into *sys_solve_prog.f90*, and all intermediate operators will be written to the `fort/testoutput` subdirectory.

Intel MKL is not specifically necessary for the compilation of these units. All of the external dependencies are concentrated in `blas_interface.f90` - its name is self-explanatory. The only external functions used are `dgemm,dtrmm` and `geev` (95), and the interface actually down- and up-converts between quad and double precision. This conversion does not appreciably decrease the accuracy of the results, as most of the floating-point error stems from the binomial and hypergeometric probabilities; the LA routines are used relatively late in the entire process. To maintain quad precision, see [here](https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=2739) for info on recompiling LAPACK accordingly.

## Interface layer
The interface files are written in C, almost entirely C99. The one caveat is the use of the non-standard `__float128` type, which is both Intel & GNU C compatible - one should still ensure that it is a valid 128-bit floating point type, as it is used in a few places solely for assignment. 

Most of the interface is specified in `f_interface.c` and its associated header file. In the manuscript, file `generatevars.c` uses a Sobol quasirandom number generator to iterate through a range of variables, whose limits can be set in `generatevars.c`. This can add to runtime, however, as the interface attempts to reuse any workspace variables which are not altered from run to run - e.g. the *M* operator is only regenerated if the connectivity coefficients *m_p_nz, m_deg_low* or *m_deg_high* are changed. Also in `generatevars.c` one can change the number of different miRNA levels tested through the macros `N_MICRO_*`. 

The program's main entry point is in `call_gensln.c`, which has a simple `demo` function that can be used to see how the workspace variables are changed between iterations. Parameters in the demo function can be altered, as can the output file for results.

Alternatively, the `runSolver` method can be used, with the number of iterations and output file specified directly.

## Results files
The output of the solver is a flat tab-delimited text file with a header specifying the particular entries. Generally, it looks like:

|nMess | nMicro | tf_deg_low | tf_deg_high| m_deg_low |  m_deg_high | tf_p_nz | defect | pmr | m_p_nz | rho | p_zero | p_one | entropy_rate | expression_mean | expression_std | condent_0 | ... | condent_150 |
|----|----|---|----|---|----|----|---|----|---|----|----|---|----|---| --- | --- | --- | --- |
| 130 | 0 | 3 | 20 | 14 | 45 | 0.2 | 0.3 | 1.5 | 0.4 | 0.85 | 0.01 | 0.01 | H_0 | &mu; | &sigma; | h_0 | ... | h_150|
| 130 | 10 | 3 | 20 | 14 | 45 | 0.2 | 0.3 | 1.5 | 0.4 | 0.85 | 0.01 | 0.01 | H_0 | &mu; | &sigma; | h_0 | ... | h_150|

where the value nMicro changes first. Regardless of what `nMess` is set to, 151 entries will be listed, with entries from `nMess+1` onwards being equal to 0. This output is fairly self-explanatory and easy to work with.

### Simulation
I am currently cleaning up/improving the readability of the CUDA-based associated simulation. Its source will be posted shortly. As one might expect, it does depend on the CUDA toolkit, v9.0 and above.

### Caveats
Beware of certain hard-coded values such as *MESSMAX* & *MICROMAX* (both on Fortran and C sides). These can be altered as needed, but obviously must be >= *nMess*, *nMicro* respectively. Also note that a failure to find a suitable stationary distribution will report as a mean expression of 0.0; origins can be probed by substituting in the `sys_solve_cplte` method in `f_interface.c`.
### Questions
The model is given here primarily for review purposes. It should be easy to build and test for someone with some experience with the methods and compilation process. The flipside of this is that it *is not* intended to be universally-accessible software. It is quite compute-intensive and does require some more involved mathematical routines, and does not necessarily report user-friendly errors. If you have any questions or concerns, please feel free to reach out and I will be glad to help you with whatever I can. 
