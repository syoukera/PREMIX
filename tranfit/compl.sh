gfortran -c tranfit.f
gfortran -c cklib.f
gfortran -c dmath.f -fPIC -fno-range-check

gfortran -o tranfite tranfit.o cklib.o dmath.o