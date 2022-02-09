touch *.*
rm ./plate_advancing.tst -rf
CC='mpicc -D_MPI=4' make plate_advancing.tst
