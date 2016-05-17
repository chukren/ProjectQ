ifort -132 -save -o  love love.for ltable-love.for
ifort -132 -save -o  rayleigh rayleigh.for ltable-rayleigh.for
ifort -132 -save -o  separate separate.for
ifort -132 -save -o  read-eigenfunction   read-eigenfunction.for 
ifort -o loveq loveQ.f90
ifort -o rayleighq rayleighQ.f90


