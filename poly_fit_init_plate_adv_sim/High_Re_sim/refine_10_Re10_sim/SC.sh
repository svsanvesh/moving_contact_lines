touch *.*
rm ./70d_plate_adv_air_water_r10_Re10_init_long_run.tst
CC='mpicc -D_MPI=4' make 70d_plate_adv_air_water_r10_Re10_init_long_run.tst
