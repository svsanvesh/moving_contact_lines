touch *.*
rm ./30d_plate_adv_air_water_r10_Re0-3045_poly_fit_init.tst -rf
CC='mpicc -D_MPI=4' make 30d_plate_adv_air_water_r10_Re0-3045_poly_fit_init.tst
