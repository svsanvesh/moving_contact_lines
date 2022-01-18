touch *.*
rm ./65d_plate_adv_air_water_r12_Re0-0677_polyfit_init_long_run.tst
CC='mpicc -D_MPI=4' make 65d_plate_adv_air_water_r12_Re0-0677_polyfit_init_long_run.tst
