touch *.*
rm ./30d_plate_adv_air_water_r10_Re10_polyfit_init_long_run_high_ST.tst
CC='mpicc -D_MPI=4' make 30d_plate_adv_air_water_r10_Re10_polyfit_init_long_run_high_ST.tst
