#!/bin/bash
ifort -c thermalstart.f90
ifort -o Q-test-yamauchi  Q-test-yamauchi.f90 thermalstart.o
