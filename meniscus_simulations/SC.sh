touch *.*
rm ./30d_plate_adv_air_water_r9_Re0-3045_poly_fit_init.tst
rm ./30d_plate_adv_air_water_r9_Re0-3045_poly_fit_init/measures
rm ./30d_plate_adv_air_water_r9_Re0-3045_poly_fit_init/xyzu*
CC='mpicc -D_MPI=4' make 30d_plate_adv_air_water_r9_Re0-3045_poly_fit_init.tst
