touch *.*
rm ./60d_plate_advancing_simulation.tst -rf
CC='mpicc -D_MPI=4' make 60d_plate_advancing_simulation.tst
