femex v.1.1.0
=============

Research use for coding with PDE/FEM in 2D. Tested on Matlab 2013a with gcc-4.8 and gfortran-4.4.

System Requirements:
* g++/gcc-4.7+ (for c++11)
* gfortran-4.7+ and gfortran-4.4
* Matlab 2013a+(gcc-4.4+ compatible)

if you would like to use AGMG-3.0 as solver, ``gfortran-4.4`` is also needed. 
Thanks to ``ABI``(http://www.fortran90.org/src/faq.html#abi), gfortran 4.7+ can link libraries on compiled files from gfortran 4.4, since the library is backward compatible.

ILUPACK is already included for a choice of solvers.



How to use:
-----------
* Modify ``Makefile.in`` to set ``MATLAB_ROOT`` to the directory which includes ``bin``.

* Simply add femex directory to current workspace path.

* Follow the example to set up boundary conditions.

 
Note:
------
* Since ``femex`` is built with Matlab's ``mex`` interface, it will need license for use which is not included.

