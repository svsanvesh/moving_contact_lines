touch *.*
rm ./30d_plate_adv_air_water_r10_Re0-0677_polyfit_init_long_run.tst -rf
CC='mpicc -D_MPI=12' make 30d_plate_adv_air_water_r10_Re0-0677_polyfit_init_long_run.tst
