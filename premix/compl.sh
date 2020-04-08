#!/bin/bash

gfortran -c dmath.f -fPIC -fno-range-check
gfortran -c twopnt.f
gfortran -c tranlib.f
gfortran -c cklib.f
gfortran -c driv.f
gfortran -c premix.f

gfortran -o premixe premix.o driv.o cklib.o tranlib.o twopnt.o dmath.o