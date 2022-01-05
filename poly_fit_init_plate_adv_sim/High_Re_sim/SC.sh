touch *.*
rm ./30d_plate_adv_air_water_r9_Re10_polyfit_init.tst
CC='mpicc -D_MPI=4' make 30d_plate_adv_air_water_r9_Re10_polyfit_init.tst
