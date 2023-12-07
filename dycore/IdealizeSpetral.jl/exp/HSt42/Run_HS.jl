using JGCM

include("HS.jl")

#############################################################
end_day = 2
spinup_day = 0

physics_params = Dict{String,Float64}("σ_b"=>0.7, "k_f" => 1.0, "k_a" => 1.0/40.0, "k_s" => 1.0/4.0, "ΔT_y" => 65.0, "Δθ_z" => 10.0) ### 60.0
op_man = Atmos_Spectral_Dynamics_Main(physics_params, end_day, spinup_day)
# Finalize_Output!(op_man, "RH50_500_100_final.dat", "RH50_500_100.dat")
# Finalize_Output!(op_man, "HS_front_RH80_PR0_PRRELAX86400_with_w_all_final.dat", "HS_front_RH80_PR0_PRRELAX86400_with_w_all.dat")

# Finalize_Output!(op_man, "RH80_test_f0_final.dat", "RH80_test_f0_all.dat")
# Finalize_Output!(op_man, "100day_RH50_PR10_test_qp_equal_q_q_ref_t_eq_all.dat", "100day_RH50_PR10_test_qp_equal_q_q_ref_t_eq_final.dat")
# Finalize_Output!(op_man, "500day_test_all.dat", "500day_test_final.dat")
# Finalize_Output!(op_man, "PR10_400day_test_ice_all.dat", "400day_test_ice_final.dat")
# Finalize_Output!(op_man, "PR10_200day_test_ice_all.dat", "800day_test_ice_final.dat")
# Finalize_Output!(op_man, "PR20_200day_test_ice_all.dat", "PR20_400day_test_ice_final.dat")
# Finalize_Output!(op_man, "PR30_200day_test_ice_all.dat", "PR30_400day_test_ice_final.dat")

# Finalize_Output!(op_man, "400day_RH50_PR10_test_new_tq_all.dat", "400day_RH50_PR10_test_new_tq_final.dat")

# Finalize_Output!(op_man, "100day_test_h_advection_final.dat", "100day_test_h_advection_all.dat")
# Finalize_Output!(op_man, "100day_test_RH80_original_h_advection_final.dat", "100day_test_RH80_original_h_advection_all.dat")
# Finalize_Output!(op_man, "day200_hope_test_final.dat", "day200_hope_test_all.dat")

# Finalize_Output!(op_man, "day50_qv_v_test_final.dat", "day50_qv_v_test_all.dat")
# Finalize_Output!(op_man, "day200_qv_v_test_final.dat", "day200_qv_v_test_all.dat")

Finalize_Output!(op_man, "test_final.dat", "test_all.dat")
# Finalize_Output!(op_man, "1205_25day_factor123_tracers_c_final.dat", "1205_25day_factor123_tracers_c_all.dat")

# Finalize_Output!(op_man, "1205_25day_factor1_with_tracers_c_and_factor3_final.dat", "1205_25day_factor1_with_tracers_c_and_factor3_all.dat")

# Finalize_Output!(op_man, "1204_25day_factor123_final.dat", "1204_25day_factor123_all.dat")
# Finalize_Output!(op_man, "1203_25day_factor123_final.dat", "1203_25day_factor123_all.dat")

# Finalize_Output!(op_man, "revise_geopot_final.dat", "revise_geopot_all.dat")
# Finalize_Output!(op_man, "200day_revise_geopot_final.dat", "200day_revise_geopot_all.dat")

# Finalize_Output!(op_man, "25day_factor2_new_final.dat", "25day_factor2_new_all.dat")
# Finalize_Output!(op_man, "25day_from_dry_final.dat", "25day_from_dry_all.dat")
# Finalize_Output!(op_man, "50day_from_dry_final.dat", "50day_from_dry_all.dat")
# Finalize_Output!(op_man, "100day_from_dry_final.dat", "100day_from_dry_all.dat")

# Finalize_Output!(op_man, "correction2_100day_from_dry_final.dat", "correction2_100day_from_dry_all.dat")
# Finalize_Output!(op_man, "correction2_300day_from_dry_final.dat", "correction2_300day_from_dry_all.dat")

# Finalize_Output!(op_man, "300day_from_dry_final.dat", "300day_from_dry_all.dat")

# Finalize_Output!(op_man, "25day_L10_final.dat", "25day_L10_all.dat")
# Finalize_Output!(op_man, "100day_L10_final.dat", "100day_L10_all.dat")

# Finalize_Output!(op_man, "200day_factor2_new_final.dat", "200day_factor2_new_all.dat")
# Finalize_Output!(op_man, "400day_factor2_new_final.dat", "400day_factor2_new_all.dat")
# Finalize_Output!(op_man, "800day_factor2_new_final.dat", "800day_factor2_new_all.dat")

# Finalize_Output!(op_man, "800day_finalcorection_final.dat", "800day_finalcorection_all.dat")

# Finalize_Output!(op_man, "25day_pqpz_final.dat", "25day_pqpz_all.dat")
# Finalize_Output!(op_man, "200day_pqpz_final.dat", "200day_pqpz_all.dat")
# Finalize_Output!(op_man, "800day_pqpz_final.dat", "800day_pqpz_all.dat")




# Finalize_Output!(op_man, "50day_ice_factor3_final.dat", "50day_ice_factor3_all.dat")
# Finalize_Output!(op_man, "100day_ice_factor3_final.dat", "100day_ice_factor3_all.dat")
# Finalize_Output!(op_man, "200day_ice_factor3_final.dat", "200day_ice_factor3_all.dat")

# Finalize_Output!(op_man, "25day_test_final.dat", "25day_test_all.dat")
# Finalize_Output!(op_man, "25day_nonrefactor_test_final.dat", "25day_nonrefactor3_test_all.dat")

# Finalize_Output!(op_man, "300day_test_final.dat", "300day_test_all.dat")
# Finalize_Output!(op_man, "500day_test_final.dat", "500day_test_all.dat")
# Finalize_Output!(op_man, "800day_test_final.dat", "800day_test_all.dat")
# Finalize_Output!(op_man, "100day_test_final.dat", "100day_test_all.dat")


# Finalize_Output!(op_man, "day200_hope_test_final.dat", "day200_hope_test_all.dat")
# Finalize_Output!(op_man, "day800_hope_test_final.dat", "day800_hope_test_all.dat")




# Sigma_Zonal_Mean_Contourf(op_man, "Contourf")

"""
try 6 hour output end_day = 300, spinup = 100
WARNING:
remember to edit initial field before run the model!!!


day0_500: end_day = 500, spinup = 100

day501_1000: end_day = 500, spinup = 0

day1001_1500: end_day = 500, spinup = 0

day1501_2000: end_day = 500, spinup = 0

"""






