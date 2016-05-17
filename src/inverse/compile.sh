#!/bin/bash
#ifort -o kernel kernel.f90 -lm
ifort -o Qmu-lsqrinv  Qmu-lsqrinv.f90 -lm
