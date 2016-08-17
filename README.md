femex v.2.1.0 
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
* For those who hadn't change the default gcc/gfortran .so files within Matlab, you should PRELOAD the system's gcc/gfortran libs first before starting Matlab.
* ``femex`` only works in shell mode, window mode will report TLS problem.

In your shell, type in something like:

```
export LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/4.8/libstdc++.so:/usr/lib/gcc/x86_64-linux-gnu/4.8/libgcc_s.so:/usr/lib/gcc/x86_64-linux-gnu/4.8/libgfortran.so
```
then start Matlab using

```
matlab -nodesktop -nosplash
```

* Back up Matlab's out-of-date ``libstdc++`` dynamic library files, try to use system's library.

* Modify ``Makefile.in`` to set ``MATLAB_ROOT`` to the directory which includes ``bin``.

* Simply add femex directory to current workspace path.

* Follow the example to set up boundary conditions.

 
Note:
------
* Since ``femex`` is built with Matlab's ``mex`` interface, it will need license for use which is not included.

