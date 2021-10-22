#!/bin/bash

ll=4
lt=80
b=4.0
u=-4.0
mu=0.0
nbin=20
nsweep=100
cb="false"
equal_measure="true"
dynamic_measure="true"

../build/dqmc_hubbard --ll=${ll} --lt=${lt} --beta=${b} --u=${u} --mu=${mu} --checkerboard=${cb} --eqtime=${equal_measure} --dynamic=${dynamic_measure} --nbin=${nbin} --nsweep=${nsweep}

exit 0