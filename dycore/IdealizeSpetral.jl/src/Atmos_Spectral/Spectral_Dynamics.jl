using Statistics
using Interpolations
export Compute_Corrections_Init, Compute_Corrections!, Four_In_One!, Spectral_Dynamics!, Get_Topography!, Spectral_Initialize_Fields!, Spectral_Dynamics_Physics!, Atmosphere_Update!



function Compute_Corrections_Init(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3}, grid_ps_p::Array{Float64, 3}, grid_t_p::Array{Float64, 3}, 
    grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, grid_δt::Array{Float64, 3},  
    Δt::Int64, grid_energy_temp::Array{Float64, 3}, grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_δtracers::Array{Float64,3}, grid_tracers_all::Array{Float64,3})
    
    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    
    sum_tracers_p = 0.

    if (do_mass_correction) 
        mean_ps_p = Area_Weighted_Global_Mean(mesh, grid_ps_p)
    end
    
    if (do_energy_correction) 
        # due to dissipation introduced by the forcing
        cp_air, grav = atmo_data.cp_air, atmo_data.grav
        grid_energy_temp  .=  0.5*((grid_u_p + Δt*grid_δu).^2 + (grid_v_p + Δt*grid_δv).^2) + cp_air*(grid_t_p + Δt*grid_δt)
        mean_energy_p = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_p)
        ###
    end

    if (do_water_correction)
        # error("water correction has not implemented")
        mean_moisture_p     =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_p, grid_ps_p)

    end
    
    return mean_ps_p, mean_energy_p, mean_moisture_p 
end 

function Compute_Corrections!(semi_implicit::Semi_Implicit_Solver, vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    mean_ps_p::Float64, mean_energy_p::Float64, mean_moisture_p::Float64,
    grid_u_n::Array{Float64, 3}, grid_v_n::Array{Float64, 3},
    grid_energy_temp::Array{Float64, 3}, grid_ps_p::Array{Float64, 3},grid_ps_c::Array{Float64, 3},
    grid_ps_n::Array{Float64, 3}, spe_lnps_n::Array{ComplexF64, 3}, 
    grid_t_n::Array{Float64, 3}, spe_t_n::Array{ComplexF64, 3},
    grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_tracers_n::Array{Float64, 3}, 
    grid_t::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_p_half::Array{Float64, 3}, grid_z_full::Array{Float64, 3}, grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3},
    grid_geopots::Array{Float64, 3}, grid_w_full::Array{Float64,3}, grid_t_p::Array{Float64, 3}, dyn_data::Dyn_Data, grid_δt::Array{Float64,3}, factor1::Array{Float64,3}, factor2::Array{Float64,3})#,  grid_tracers_diff::Array{Float64, 3})


    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    
    
    if (do_mass_correction) 
        mean_ps_n = Area_Weighted_Global_Mean(mesh, grid_ps_n)
        mass_correction_factor = mean_ps_p/mean_ps_n
        grid_ps_n .*= mass_correction_factor
        #P00 = 1 
        spe_lnps_n[1,1,1] += log(mass_correction_factor)
    end
    
    if (do_energy_correction) 
        cp_air, grav = atmo_data.cp_air, atmo_data.grav
        
        grid_energy_temp .=  0.5*(grid_u_n.^2 + grid_v_n.^2) + cp_air*grid_t_n
        mean_energy_n         = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_n)
        
        temperature_correction = grav*(mean_energy_p - mean_energy_n)/(cp_air*mean_ps_p)
        #@info grav, mean_energy_p , mean_energy_n, cp_air, mean_ps_p
        grid_t_n       .+= temperature_correction
        spe_t_n[1,1,:] .+= temperature_correction
    end

    # @info mean_ps_p, mean_energy_p, mass_correction_factor, temperature_correction
    ### By CJY 0517
    nλ = mesh.nλ
    nθ = mesh.nθ
    nd = mesh.nd
    grav = atmo_data.grav
    integrator = semi_implicit.integrator
    Δt = Get_Δt(integrator)

    if (do_water_correction) 
        # ### original ###
        # factor11             = zeros(((128,64,20)))
        # factor11            .= factor1 .+ factor2
        # @info maximum(factor11), minimum(factor11)

        # if maximum(factor11) != 0. #&& minimum(factor11) != 0.
        #     grid_tracers_n[grid_tracers_n .< 0.] .= 0.
        #     mean_factor11_n      =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, factor11, grid_ps_n)
        #     mean_moisture_n      =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        #     factor11           .*= (mean_moisture_p-mean_moisture_n)/mean_factor11_n
        #     grid_tracers_n     .+=  factor11
        #     grid_tracers_n[grid_tracers_n .< 0.] .= 0.
        #     @info "add"
        #     @info maximum(factor1), minimum(factor1)
        #     @info maximum(factor2), minimum(factor2)

        # else
        #     grid_tracers_n[grid_tracers_n .< 0.] .= 0.
        #     mean_moisture_n       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        #     grid_tracers_n     .*= mean_moisture_p ./ mean_moisture_n 
        #     grid_tracers_n[grid_tracers_n .< 0.] .= 0.
        #     @info "plus"
        # end
        ###

        # mean_factor_all_n      =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, factor_all, grid_ps_n)
        # mean_moisture_n       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        # factor_all           .*= (mean_moisture_p-mean_moisture_n)/mean_factor_all_n
        # grid_tracers_n     .+=  factor_all
        # grid_tracers_n[grid_tracers_n .< 0.] .= 0.
        # mean_factor_all_n      =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, factor_all, grid_ps_n)
        # factor_all           .*= (mean_moisture_p-mean_moisture_n)/mean_factor_all_n

        grid_tracers_n[grid_tracers_n .< 0.] .= 0.
        mean_moisture_n       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        grid_tracers_n     .*= mean_moisture_p ./ mean_moisture_n 
        grid_tracers_n[grid_tracers_n .< 0.] .= 0.

        mean_moisture_n  =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        
        ### 10/30 
        @info "#### mass correction:", (mean_moisture_n - mean_moisture_p)
        return mean_moisture_n
    end
    
end 


"""
compute vertical mass flux and velocity 
grid_M_half[:,:,k+1] = downward mass flux/per unit area across the K+1/2
grid_w_full[:,:,k]   = dp/dt vertical velocity 

update residuals
grid_δps[:,:,k]  += -∑_{r=1}^nd Dr = -∑_{r=1}^nd ∇(vrΔp_r)
grid_δt[:,:,k]  += κTw/p 
(grid_δu[:,:,k], grid_δv[:,:,k]) -= RT ∇p/p 

!  cell boundary. This is the "vertical velocity" in the hybrid coordinate system.
!  When vertical coordinate is pure sigma: grid_M_half = grid_ps*d(sigma)/dt
"""

function Four_In_One!(vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
    grid_div::Array{Float64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
    grid_ps::Array{Float64,3},  grid_Δp::Array{Float64,3}, grid_lnp_half::Array{Float64,3}, grid_lnp_full::Array{Float64,3}, grid_p_full::Array{Float64,3},
    grid_dλ_ps::Array{Float64,3}, grid_dθ_ps::Array{Float64,3}, 
    grid_t::Array{Float64,3}, 
    grid_M_half::Array{Float64,3}, grid_w_full::Array{Float64,3}, 
    grid_δu::Array{Float64,3}, grid_δv::Array{Float64,3}, grid_δps::Array{Float64,3}, grid_δt::Array{Float64,3}, grid_δtracers::Array{Float64,3})
    
    rdgas, cp_air = atmo_data.rdgas, atmo_data.cp_air
    nd, bk = vert_coord.nd, vert_coord.bk
    Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
    vert_difference_option = vert_coord.vert_difference_option
    
    kappa = rdgas / cp_air
    
    # dmean_tot = ∇ ∑_{k=1}^{nd} vk Δp_k = ∑_{k=1}^{nd} Dk
    nλ, nθ, _ = size(grid_ps)
    dmean_tot = zeros(Float64, nλ, nθ)
    Δlnp_p = zeros(Float64, nλ, nθ)
    Δlnp_m = zeros(Float64, nλ, nθ)
    Δlnp = zeros(Float64, nλ, nθ)
    x1 = zeros(Float64, nλ, nθ)
    dlnp_dλ = zeros(Float64, nλ, nθ)
    dlnp_dθ = zeros(Float64, nλ, nθ)
    dmean = zeros(Float64, nλ, nθ)
    x5 = zeros(Float64, nλ, nθ)
        
    if (vert_difference_option == "simmons_and_burridge") 
        for k = 1:nd
        Δp = grid_Δp[:,:,k]
        
        Δlnp_p .= grid_lnp_half[:,:,k + 1] - grid_lnp_full[:,:,k]
        Δlnp_m .= grid_lnp_full[:,:,k]   - grid_lnp_half[:,:,k]
        Δlnp   .= grid_lnp_half[:,:,k + 1] - grid_lnp_half[:,:,k]
        
        # angular momentum conservation 
        #    ∇p_k/p =  [(lnp_k - lnp_{k-1/2})∇p_{k-1/2} + (lnp_{k+1/2} - lnp_k)∇p_{k+1/2}]/Δpk
        #         =  [(lnp_k - lnp_{k-1/2})B_{k-1/2} + (lnp_{k+1/2} - lnp_k)B_{k+1/2}]/Δpk * ∇ps
        #         =  x1 * ∇ps
        x1 .= (bk[k] * Δlnp_m + bk[k + 1] * Δlnp_p ) ./ Δp
        
        dlnp_dλ .= x1 .* grid_dλ_ps[:,:,1]
        dlnp_dθ .= x1 .* grid_dθ_ps[:,:,1]
        
        
        # (grid_δu, grid_δv) -= RT ∇p/p 
        grid_δu[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dλ
        grid_δv[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dθ
        
        # dmean = ∇ (vk Δp_k) =  divk Δp_k + vk  Δbk[k] ∇ p_s
        dmean .= grid_div[:,:,k] .* Δp + Δbk[k] * (grid_u[:,:,k] .* grid_dλ_ps[:,:,1] + grid_v[:,:,k] .* grid_dθ_ps[:,:,1])
        
    
        # energy conservation for temperature
        # w/p = dlnp/dt = ∂lnp/∂t + dσ ∂lnp/∂σ + v∇lnp
        # dσ ∂ξ_k/∂σ = [M_{k+1/2}(ξ_k+1/2 - ξ_k) + M_{k-1/2}(ξ_k - ξ_k-1/2)]/Δp_k
        # weight the same way (TODO)
        # vertical advection operator (M is the downward speed)
        # dσ ∂lnp_k/∂σ = [M_{k+1/2}(lnp_k+1/2 - lnp_k) + M_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
        # ∂lnp/∂t = 1/p ∂p/∂t = [∂p/∂t_{k+1/2}(lnp_k+1/2 - lnp_k) + ∂p/∂t_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
        # As we know
        # ∂p/∂t_{k+1/2} = -∑_{r=1}^k Dr - M_{k+1/2}
        
        # ∂lnp/∂t + dσ ∂lnp/∂σ =  [(-∑_{r=1}^k Dr)(lnp_k+1/2 - lnp_k) + (-∑_{r=1}^{k-1} Dr)(lnp_k - lnp_k-1/2)]/Δp_k
        #                      = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k
        
        x5 .= -(dmean_tot .* Δlnp + dmean .* Δlnp_p) ./ Δp .+ grid_u[:,:,k] .* dlnp_dλ + grid_v[:,:,k] .* dlnp_dθ
        # grid_δt += κT w/p
        grid_δt[:,:,k] .+=  kappa * grid_t[:,:,k] .* x5
        # grid_w_full = w
        grid_w_full[:,:,k] .= x5 .* grid_p_full[:,:,k]
        # update dmean_tot to ∑_{r=1}^k ∇(vrΔp_r)
        dmean_tot .+= dmean
        # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
        grid_M_half[:,:,k + 1] .= -dmean_tot
        end
        
    else
        error("vert_difference_option ", vert_difference_option, " is not a valid value for option")
        
    end
    # ∂ps/∂t = -∑_{r=1}^nd ∇(vrΔp_r) = -dmean_tot
    grid_δps[:,:,1] .-= dmean_tot
    
    for k = 1:nd-1
        # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
        grid_M_half[:,:,k+1] .+= dmean_tot * bk[k+1]
    end
    
    grid_M_half[:,:,1] .= 0.0
    grid_M_half[:,:,nd + 1] .= 0.0
end 



"""
The governing equations are
∂div/∂t = ∇ × (A, B) - ∇^2E := f^d                    
∂lnps/∂t= (-∑_k div_k Δp_k + v_k ∇ Δp_k)/ps := f^p    
∂T/∂t = -(u,v)∇T - dσ∂T∂σ + κTw/p + J:= f^t           
Φ = f^Φ                                               

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 

We have 
δdiv = f^d - I^d + I^d
δlnps = f^p - I^p + I^p
δT = f^t - I^t + I^t

"""

function Spectral_Dynamics!(mesh::Spectral_Spherical_Mesh,  vert_coord::Vert_Coordinate, 
    atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
    semi_implicit::Semi_Implicit_Solver)
    
    # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor     = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv     = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt             = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    
    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n


    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps = dyn_data.grid_lnps
    
    grid_div, grid_absvor, grid_vor = dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half = dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half
    
    grid_energy_full, spe_energy = dyn_data.grid_energy_full, dyn_data.spe_energy
    
    # By CJY2
    spe_tracers_n = dyn_data.spe_tracers_n
    spe_tracers_c = dyn_data.spe_tracers_c
    spe_tracers_p = dyn_data.spe_tracers_p 
    
    grid_tracers_n = dyn_data.grid_tracers_n
    grid_tracers_c = dyn_data.grid_tracers_c
    grid_tracers_p = dyn_data.grid_tracers_p 
    
    grid_tracers_diff = dyn_data.grid_tracers_diff
    
    spe_δtracers   = dyn_data.spe_δtracers  
    grid_δtracers  = dyn_data.grid_δtracers 
    grid_tracers_full = dyn_data.grid_tracers_full
    
    ### 11/07
    grid_z_full = dyn_data.grid_z_full
    grid_z_half = dyn_data.grid_z_half
    ###
    grid_w_full = dyn_data.grid_w_full
    ### 
    add_water   = dyn_data.add_water
    ###
    # todo !!!!!!!!
    #  grid_q = grid_t
    ###
    grav = atmo_data.grav
    integrator = semi_implicit.integrator
    Δt = Get_Δt(integrator)
    factor1 = dyn_data.factor1 
    factor2 = dyn_data.factor2 
    factor3 = dyn_data.factor3  
    # factor4 = dyn_data.factor4  

    grid_z_full = dyn_data.grid_z_full
    grid_z_half = dyn_data.grid_z_half
    grid_δtracers = dyn_data.grid_δtracers 

    K_E = dyn_data.K_E
    pqpz = dyn_data.pqpz

    # grid_ps = dyn_data.grid_ps_c

    ###############################################################################
    ###
    # original 
    # pressure difference
    grid_Δp = dyn_data.grid_Δp
    # temporary variables
    grid_δQ = dyn_data.grid_d_full1
    ### 11/12
    qv_global_intergral = dyn_data.qv_global_intergral
        
    # incremental quantities
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt

    integrator = semi_implicit.integrator
    Δt = Get_Δt(integrator)

    # Calculate latent heat and modify qv_current
    # HS_forcing_water_vapor!(grid_tracers_c,  grid_δtracers, grid_t, grid_δt, grid_p_full)

    mean_ps_p, mean_energy_p, mean_moisture_p = Compute_Corrections_Init(vert_coord, mesh, atmo_data,
    grid_u_p, grid_v_p, grid_ps_p, grid_t_p, 
    grid_δu, grid_δv, grid_δt,  
    Δt, grid_energy_full, grid_tracers_p, grid_tracers_c, grid_δtracers, grid_tracers_full)
    
    # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    ###
    ###

    ##
    C_E = 0.0044
    Lv = 2.5*10^6.
    Rv = 461.
    Rd = 287.
    cp = 1004.
    # # ### factor3
    ### use n

    # grid_δtracers .-= factor3 ./(2. .* Δt)
    ### try
    # @info maximum(grid_u)
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt
    grid_tracers_diff_new = HS_forcing_water_vapor!(semi_implicit, grid_tracers_n,  grid_t_n, grid_δt, grid_p_full, grid_u, grid_v, factor3, grid_δtracers, grid_tracers_c, grid_t)
    grid_tracers_diff    .= grid_tracers_diff_new
    grid_tracers_c[grid_tracers_c .< 0] .= 0 

    # grid_tracers_c[grid_tracers_c .< 0] .= 0 
    ### 11/08
    
    V_c  = zeros(((128,64,20)))
    V_c .= (grid_u[:,:,:].^2 .+ grid_v[:,:,:].^2).^0.5
    ### add moisture at surface following paper
    ### ∂q_a/∂t = C_E * V_a * (q_sat,a - q_a) ./ z_a 
    ### factor1
    # cal rho
    rho = zeros(((128,64,20)))
    for i in 1:20
        rho[:,:,i] .=  grid_p_half[:,:,i] ./ Rd ./ (grid_t[:,:,i])
    end
    # rho_s = zeros(((128,64,1)))
    # rho_s[:,:,1] .=  grid_ps_n[:,:,1] ./ Rd ./ (grid_t[:,:,20])

    # cal za
    tv = zeros(((128,64,1)))
    za = zeros(((128,64,1)))
    tv[:,:,1] .= grid_t[:,:,20] .* (1. .+ 0.608 .* grid_tracers_c[:,:,20])

    za[:,:,1] .= (Rd .* tv./9.81 .* (log.(grid_ps[:,:,1]) ./ ((grid_p_full[:,:,20] .+ grid_p_half[:,:,21]) ./ 2.) ) ./2)
    ### 

    grid_tracers_c_max  = deepcopy(grid_tracers_c)
    grid_tracers_c_max .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) )) ./ (grid_p_full .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) )) 

    ### 1205 1205_25day_factor1_with_tracers_c_and_factor3_final.dat test factor1 use grid_tracers_c and grid_tracers_c_max and add on grid_δtracers, without factor2
    factor1[:,:,20] .=  C_E .* V_c[:,:,20] .* (grid_tracers_c_max[:,:,20] .- min.(grid_tracers_c[:,:,20], grid_tracers_c_max[:,:,20])) ./ za[:,:,1] 
    # grid_tracers_c[:,:,20] .+= factor1[:,:,20]
    grid_δtracers[:,:,20] .+= factor1[:,:,20] ./(2. .* Δt)
    # ###
        ### precipitation, add water then surface flux
    # # eddy diffusivity coefficient, K_E
    # spe_δtracers   .= 0.
    # grid_δtracers  .= 0.
    """
    V_a = V_c[:,:,20]
    for i in 17:21
        K_E[:,:,i] .= C_E .* V_a .* za[:,:,1]
    end
    K_E[:,:, 1:16] .= C_E .* V_a .* za[:,:,1] .* exp.(-((grid_p_half[:,:,17] .- grid_p_half[:,:,1:16]) ./ 10000.).^2)
    ### cal PBL Scheme
    rpdel  = zeros(((128,64,20))) ### = 1 / (p^n_{+} - p^n_{-}) , which p^_{-} mean upper layer
    for i in 1:20
        rpdel[:,:,i] .= 1. ./ (grid_p_half[:,:,i+1] .- grid_p_half[:,:,i])
    end

    CA     = zeros(((128,64,20)))
    CC     = zeros(((128,64,20)))
    CE     = zeros(((128,64,20+1)))
    CF     = zeros(((128,64,20+1)))

    for k in 1:19
        CA[:,:,k]   .= (rpdel[:,:,k]   .* 2. .* Δt .* grav.^2 .* K_E[:,:,k+1]   .* rho[:,:,k+1].^2 
                       ./ (grid_p_full[:,:,k+1] .- grid_p_full[:,:,k]))
        CC[:,:,k+1] .= (rpdel[:,:,k+1] .* 2. .* Δt .* grav.^2 .* K_E[:,:,k+1]   .* rho[:,:,k+1].^2
                       ./ (grid_p_full[:,:,k+1] .- grid_p_full[:,:,k]))
    end
    # @info maximum(rpdel) ## OK
    # @info maximum(rho) ## OK
    # @info maximum(grid_p_full[:,:,18+1] .- grid_p_full[:,:,18]) # # OK
    
    CA[:,:,20]   .= 0.
    CC[:,:, 1]   .= 0.
    CE[:,:,21]   .= 0.
    CF[:,:,21]   .= 0.
    # @info minimum(CA)
    # @info minimum(CC)


    for k in 20:-1:1
        CE[:,:,k]    .= CC[:,:,k] ./ (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1])
        CF[:,:,k]    .= ((grid_tracers_c[:,:,k] .+ CA[:,:,k] .* CF[:,:,k+1])
                        ./ (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1]))
    end
    # @info maximum(CE)
    # @info maximum(CF)

    # first calculate the updates at the top model level
    # grid_δtracers[:,:,1] .+= (CF[:,:,1] .- grid_tracers_c[:,:,1]) ./ (2. .* Δt)
    ### WARNING factor1 just factor, so it did  ./ ./ (2. .* Δt). 
    ### So did factor2
    factor2[:,:,1]        .= (CF[:,:,1] .- grid_tracers_c[:,:,1]) #./ (2. .* Δt)  # because CE at top = 0
    grid_tracers_c[:,:,1] .= CF[:,:,1] 
    # Loop over the remaining level
    for k in 2:19
        grid_δtracers[:,:,k]  .+= (CE[:,:,k] .* grid_tracers_c[:,:,k-1] .+ CF[:,:,k] .- grid_tracers_c[:,:,k]) ./ (2. .* Δt)
        factor2[:,:,k]         .= (CE[:,:,k] .* grid_tracers_c[:,:,k-1] .+ CF[:,:,k] .- grid_tracers_c[:,:,k])  #./ (2. .* Δt)
        grid_tracers_c[:,:,k]  .=  CE[:,:,k] .* grid_tracers_c[:,:,k-1] .+ CF[:,:,k]
    end
    @info maximum(CE), minimum(CE)
    @info maximum(CF), minimum(CF)
    """
    ###
    # compute ∇ps = ∇lnps * ps
    Compute_Gradients!(mesh, spe_lnps_c,  grid_dλ_ps, grid_dθ_ps)
    grid_dλ_ps .*= grid_ps
    grid_dθ_ps .*= grid_ps


    
    # compute grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, 
    # except the contributions from geopotential or vertical advection
    Four_In_One!(vert_coord, atmo_data, grid_div, grid_u, grid_v, grid_ps, 
    grid_Δp, grid_lnp_half, grid_lnp_full, grid_p_full,
    grid_dλ_ps, grid_dθ_ps, 
    grid_t, 
    grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, grid_δtracers)

    Compute_Geopotential!(vert_coord, atmo_data, 
    grid_lnp_half, grid_lnp_full,  
    grid_t, 
    grid_geopots, grid_geopot_full, grid_geopot_half, grid_tracers_n)

    grid_δlnps .= grid_δps ./ grid_ps
    Trans_Grid_To_Spherical!(mesh, grid_δlnps, spe_δlnps)

    # compute vertical advection, todo  finite volume method 
    Vert_Advection!(vert_coord, grid_u, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δu  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_v, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δv  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_t, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δt  .+= grid_δQ
    ### By CJY2
    Vert_Advection!(vert_coord, grid_tracers_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme,  grid_δQ)
    grid_δtracers .+= grid_δQ 
    ### spectral tracers need to be done first 
    Add_Horizontal_Advection!(mesh, spe_tracers_c, grid_u, grid_v, grid_δtracers) 
    Trans_Grid_To_Spherical!(mesh, grid_δtracers, spe_δtracers)
    Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
    Filtered_Leapfrog!(integrator, spe_δtracers, spe_tracers_p, spe_tracers_c, spe_tracers_n)
    Trans_Spherical_To_Grid!(mesh, spe_tracers_n, grid_tracers_n)
    Add_Horizontal_Advection!(mesh, spe_t_c, grid_u, grid_v, grid_δt)

    ### original grid_tracers_diff
    # grid_tracers_diff_new = HS_forcing_water_vapor!(semi_implicit, grid_tracers_n,  grid_t_n, grid_δt, grid_p_full, grid_u, grid_v, factor3)
    # grid_tracers_diff    .= grid_tracers_diff_new
    ### 10/30
    # @info maximum(grid_tracers_diff)
    ###
    Trans_Grid_To_Spherical!(mesh, grid_δt, spe_δt)
   
    
    grid_absvor = dyn_data.grid_absvor
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    
    
    grid_δu .+=  grid_absvor .* grid_v
    grid_δv .-=  grid_absvor .* grid_u
    
    
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)

    grid_energy_full .= grid_geopot_full .+ 0.5 * (grid_u.^2 + grid_v.^2)
    Trans_Grid_To_Spherical!(mesh, grid_energy_full, spe_energy)
    Apply_Laplacian!(mesh, spe_energy)
    spe_δdiv .-= spe_energy
    
    
    
    Implicit_Correction!(semi_implicit, vert_coord, atmo_data,
    spe_div_c, spe_div_p, spe_lnps_c, spe_lnps_p, spe_t_c, spe_t_p, 
    spe_δdiv, spe_δlnps, spe_δt)
    
    Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
    Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
    Compute_Spectral_Damping!(integrator, spe_t_c, spe_t_p, spe_δt)
    # ### By CJY2
    # Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
    # ###
        
    Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
    Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
    Filtered_Leapfrog!(integrator, spe_δlnps, spe_lnps_p, spe_lnps_c, spe_lnps_n)
    Filtered_Leapfrog!(integrator, spe_δt, spe_t_p, spe_t_c, spe_t_n)
    
    
    Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_n, grid_lnps)
    grid_ps_n .= exp.(grid_lnps)
    Trans_Spherical_To_Grid!(mesh, spe_t_n, grid_t_n) 
    # spe_δtracers   .= 0.
    # grid_δtracers  .= 0.
    # factor1_loc, factor2_loc, factor3_loc, K_E_loc, rho_loc,
    # @info maximum(factor1), minimum(factor1)
    # @info maximum(factor2), minimum(factor2)
     mean_moisture_n_loc = Compute_Corrections!(semi_implicit, vert_coord, mesh, atmo_data, mean_ps_p, mean_energy_p,mean_moisture_p, 
        grid_u_n, grid_v_n,
        grid_energy_full, grid_ps_p,grid_ps,
        grid_ps_n, spe_lnps_n, 
        grid_t_n, spe_t_n, 
        grid_tracers_p, grid_tracers_c, grid_tracers_n,
        grid_t, grid_p_full, grid_p_half, grid_z_full, grid_u_p, grid_v_p, grid_geopots, grid_w_full, grid_t_p, dyn_data, grid_δt, factor1, factor2)

    # factor1 .= factor1_loc
    # factor2 .= factor2_loc
    # factor3 .= factor3_loc

    add_water .= factor1 .+ factor2

    # K_E .= K_E_loc
    # pqpz .= pqpz_loc

    # rho .= rho_loc

    qv_global_intergral .= mean_moisture_n_loc
    ###

    day_to_sec = 86400
    if (integrator.time%day_to_sec == 0)
        # dyn_data.grid_tracers_c[dyn_data.grid_tracers_c .< 0] .= 0
        @info "Day: ", div(integrator.time,day_to_sec), " Max |U|,|V|,|P|,|T|,|qv|: ", maximum(abs.(dyn_data.grid_u_c)), maximum(abs.(dyn_data.grid_v_c)), maximum(dyn_data.grid_p_full), maximum(dyn_data.grid_t_c), maximum(dyn_data.grid_tracers_c), maximum(dyn_data.grid_tracers_diff)
        @info "Day: ", div(integrator.time,day_to_sec), " Min |U|,|V|,|P|,|T|,|qv|: ", minimum(abs.(dyn_data.grid_u_c)), minimum(abs.(dyn_data.grid_v_c)), minimum(dyn_data.grid_p_full), minimum(dyn_data.grid_t_c), minimum(dyn_data.grid_tracers_c)
    end

    Time_Advance!(dyn_data)

    # @info "min dyn.grid_tracers_c" minimum(dyn_data.grid_tracers_n)
    
    #@info "sec: ", integrator.time+1200, sum(abs.(grid_u_n)), sum(abs.(grid_v_n)), sum(abs.(grid_t_n)) , sum(abs.(grid_ps_n))
    #@info "max: ", maximum(abs.(grid_u_n)), maximum(abs.(grid_v_n)), maximum(abs.(grid_t_n)) , maximum(abs.(grid_ps_n))
    #@info "loc", grid_u_n[100,30,10],  grid_t_n[100,30,10], grid_u_n[1,32,1],  grid_t_n[1,32,1]
    
    #@assert(maximum(grid_u) <= 100.0 && maximum(grid_v) <= 100.0)


    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    
    return 
end 

function Get_Topography!(grid_geopots::Array{Float64, 3})
    # grid_geopots .= 0.0
    ### 2023/10/25
    read_file = load("300day_dry_run_all.dat")
    grid_geopots .= read_file["grid_geopots_xyzt"][:,:,1,300]
    return
end 

function Spectral_Initialize_Fields!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, sea_level_ps_ref::Float64, init_t::Float64, grid_geopots::Array{Float64,3}, dyn_data::Dyn_Data, Δt::Int64, physics_params::Dict{String, Float64})
    
    spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
    spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
    grid_u, grid_v, grid_ps, grid_t = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
    grid_u_p, grid_v_p, grid_ps_p, grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
    
    grid_lnps,  grid_vor, grid_div =  dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
    
    grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
    
    ### By CJY2
    spe_tracers_n = dyn_data.spe_tracers_n
    spe_tracers_c = dyn_data.spe_tracers_c
    spe_tracers_p = dyn_data.spe_tracers_p 
        
    grid_tracers_n = dyn_data.grid_tracers_n
    grid_tracers_c = dyn_data.grid_tracers_c
    grid_tracers_p = dyn_data.grid_tracers_p 
    
    grid_tracers_diff = dyn_data.grid_tracers_diff
    ###
    ### 10/31
    grid_δu, grid_δv, grid_δps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
    grid_u_p, grid_v_p,  grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_t_p
    grid_p_half, grid_p_full = dyn_data.grid_p_half, dyn_data.grid_p_full
    
    ### 11/08
    grid_z_full = dyn_data.grid_z_full
    grid_w_full = dyn_data.grid_w_full

    
    # ### 11/02
    read_file  = load("300day_dry_run_all.dat")
    initial_day = 300
    # # ### 11/08
    grid_z_full .= read_file["grid_z_full_xyzt"][:,:,:,initial_day]
    grid_w_full .= read_file["grid_w_full_xyzt"][:,:,:,initial_day]
    # ###
    grid_δu       .= read_file["grid_δu_xyzt"][:,:,:,initial_day]
    grid_δv       .= read_file["grid_δv_xyzt"][:,:,:,initial_day]
    grid_δt       .= read_file["grid_δt_xyzt"][:,:,:,initial_day]
    grid_δps       .= read_file["grid_δps_xyzt"][:,:,1,initial_day]

    grid_p_half .= read_file["grid_p_half_xyzt"][:,:,:,initial_day]
    grid_p_full .= read_file["grid_p_full_xyzt"][:,:,:,initial_day]

    
    
    rdgas = atmo_data.rdgas
    # grid_t    .= init_t
    ### 2023/10/27
    grid_t    .= read_file["grid_t_c_xyzt"][:,:,:,initial_day] 
    ###
    # dΦ/dlnp = -RT    Δp = -ΔΦ/RT
    grid_geopots .= read_file["grid_geopots_xyzt"][:,:,:,initial_day]
    # grid_lnps .= log(sea_level_ps_ref) .- grid_geopots / (rdgas * init_t)
    grid_lnps .= read_file["grid_lnps_xyzt"][:,:,1,initial_day]
    # grid_ps   .= exp.(grid_lnps)
    grid_ps    .= read_file["grid_ps_xyzt"][:,:,1,initial_day]
    

    # By CJY
    # spe_div_c .= 0.0
    # spe_vor_c .= 0.0

    # # initial perturbation
    num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
    
    # initial_perturbation = 1.0e-7/sqrt(2.0)
    # initial vorticity perturbation used in benchmark code
    # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
    
    # for k = nd-2:nd
    #     spe_vor_c[2,5,k] = initial_perturbation
    #     spe_vor_c[6,9,k] = initial_perturbation
    #     spe_vor_c[2,4,k] = initial_perturbation  
    #     spe_vor_c[6,8,k] = initial_perturbation
    # end
    

    ##
    spe_vor_c[:,:,:] .= read_file["spe_vor_c_xyzt"][:,:,:,initial_day]
    spe_div_c[:,:,:] .= read_file["spe_div_c_xyzt"][:,:,:,initial_day]
    grid_u[:,:,:]    .= read_file["grid_u_c_xyzt"][:,:,:,initial_day]
    grid_v[:,:,:]    .= read_file["grid_v_c_xyzt"][:,:,:,initial_day]  
    
    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
    
    # initial spectral fields (and spectrally-filtered) grid fields
    Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
    Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)

    Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps_c)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_c,  grid_lnps)

    # By CJY ### grid_ps .= exp.(grid_lnps)
    #grid_ps .= read_file["grid_ps_xyzt"][:,:,1,90]
    grid_ps .= exp.(grid_lnps)
    
    Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)

    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
    
    Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)

    #update pressure variables for hs forcing
    grid_p_half   .= read_file["grid_p_half_xyzt"][:,:,:,initial_day]
    grid_Δp       .= read_file["grid_Δp_xyzt"][:,:,:,initial_day]
    grid_lnp_half .= read_file["grid_lnp_half_xyzt"][:,:,:,initial_day]
    grid_p_full   .= read_file["grid_p_full_xyzt"][:,:,:,initial_day]
    grid_lnp_full .= read_file["grid_lnp_full_xyzt"][:,:,:,initial_day]
        
    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
    grid_lnp_half, grid_p_full, grid_lnp_full)
    
    
    spe_vor_p  .= spe_vor_c
    spe_div_p  .= spe_div_c
    spe_lnps_p .= spe_lnps_c
    spe_t_p    .= spe_t_c


    grid_u_p   .= grid_u
    grid_v_p   .= grid_v
    grid_ps_p  .= grid_ps
    grid_t_p   .= grid_t

    # Tracer initialization
    initial_RH      = 0.8
    Lv              = 2.5*10^6.
    Rv              = 461.
    ### 2023/10/25
    grid_tracers_c .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) .* initial_RH)) ./ (grid_p_full .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) .* initial_RH))     
    qv0 = 0.018
    # λc = mesh.λc # lon
    θc = mesh.θc # lat
    phi_hw = 2 * pi / 9 * deg2rad(40)
    p_hw = 30000.
    phi = LinRange(-90,90,64)
    p0 = 100000.
    # for k in 1:20
    #     for j in 1:64
    #        for i in 1:128
    #            grid_tracers_c[i,j,k] = qv0 * exp(-((grid_p_full[i,j,k]/grid_ps[i,j,1] - 1.)*(p0/p_hw))^2) * exp(-((deg2rad(phi[j]))/phi_hw)^4) 
    #        end            
    #     end
    # end
    # grid_tracers_c[:,:,1] .= 0.
    ###
    # grid_tracers_c .= read_file["grid_tracers_c_xyz1t"][:,:,:,initial_day]
    ###
    
    ###
    Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
    Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
    
    
    # By CJY2
    spe_tracers_p  .= spe_tracers_c
    grid_tracers_p .= grid_tracers_c
    

     
end 


function Spectral_Dynamics_Physics!(atmo_data::Atmo_Data, mesh::Spectral_Spherical_Mesh, dyn_data::Dyn_Data, Δt::Int64, physics_params::Dict{String, Float64})
    grid_δu, grid_δv, grid_δps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
    grid_u_p, grid_v_p,  grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_t_p
    grid_p_half, grid_p_full = dyn_data.grid_p_half, dyn_data.grid_p_full
    grid_t_eq = dyn_data.grid_t_eq

    grid_δtracers = dyn_data.grid_δtracers
    spe_δtracers  = dyn_data.spe_δtracers

    
    grid_δps .= 0.0

    spe_δtracers   .= 0.
    grid_δtracers  .= 0.
    
    HS_Forcing!(atmo_data, Δt, mesh.sinθ, grid_u_p, grid_v_p, grid_p_half, grid_p_full, grid_t_p, grid_δu, grid_δv,
    grid_t_eq, grid_δt, physics_params)

    
end


function Atmosphere_Update!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, semi_implicit::Semi_Implicit_Solver, 
                            dyn_data::Dyn_Data, physcis_params::Dict{String, Float64})

    Δt = Get_Δt(semi_implicit.integrator)
    Spectral_Dynamics_Physics!(atmo_data, mesh,  dyn_data, Δt, physcis_params)
    Spectral_Dynamics!(mesh,  vert_coord , atmo_data, dyn_data, semi_implicit)
    ### original
    grid_ps , grid_Δp, grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_ps_c,  dyn_data.grid_Δp, dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full 
    
    grid_t = dyn_data.grid_t_c
    grid_geopots, grid_z_full, grid_z_half = dyn_data.grid_geopots, dyn_data.grid_z_full, dyn_data.grid_z_half

    ### 1201
    grid_tracers_n = dyn_data.grid_tracers_n
        
    Compute_Pressures_And_Heights!(atmo_data, vert_coord,     
    grid_ps, grid_geopots, grid_t, 
    grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full, grid_z_full, grid_z_half, grid_tracers_n)

    return
end 


function HS_forcing_water_vapor!(semi_implicit::Semi_Implicit_Solver, grid_tracers_n::Array{Float64, 3},  grid_t_n::Array{Float64, 3}, grid_δt::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_u::Array{Float64, 3},  grid_v::Array{Float64, 3}, factor3::Array{Float64, 3}, grid_δtracers::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3},grid_t::Array{Float64, 3})

    integrator = semi_implicit.integrator
    Δt = Get_Δt(integrator)
    cp  = 1004.
    Lv = 2.5*10^6.
    Rd = 287.04
    Rv = 461.
    # C_E = 0.0044
    grid_tracers_diff      = zeros(size(grid_tracers_c)...)
    grid_tracers_c_max     = zeros(size(grid_tracers_c)...)
    
    grid_tracers_c_max  = deepcopy(grid_tracers_c)
    grid_tracers_c_max .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) )) ./ (grid_p_full .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) )) 
    grid_tracers_diff  .= (max.(grid_tracers_c, grid_tracers_c_max) .- grid_tracers_c_max) 
    grid_tracers_c .-= (max.(grid_tracers_c, grid_tracers_c_max) .- grid_tracers_c_max) 
    #grid_δtracers .-= (max.(grid_tracers_n, grid_tracers_n_max) .- grid_tracers_n_max) ./ (2 .* Δt)
    


    # FIXME
    day_to_sec = 86400.
    L = 0.1
    ### 2023/10/28
    # grid_tracers_diff  .= max.(0, grid_tracers_c .- grid_tracers_c_max) #./ (1 .+ Lv ./ cp .* Lv .* grid_tracers_c_max ./ Rv ./ grid_t .^2) #./ (2*Δt)
    # # ### ice
    # Ls = 2.83 * 10^6.
    # grid_t_minus   = zeros(size(grid_tracers_c)...)
    # grid_t_minus_final    = zeros(size(grid_tracers_c)...)
    # grid_t_minus_final2   = zeros(size(grid_tracers_c)...)

    
    # grid_t_plus   = zeros(size(grid_tracers_c)...)
    # grid_t_plus_final    = zeros(size(grid_tracers_c)...)
    # grid_t_plus_final2   = zeros(size(grid_tracers_c)...)
    
    # grid_t_minus .= min.(0., grid_t.-273.15)
    # grid_t_minus_final  .= ifelse.(grid_t_minus .< 0, grid_t_minus, grid_t_minus .= 0)
    # grid_t_minus_final2 .= ifelse.(grid_t_minus_final .> 0, grid_t_minus_final, grid_t_minus .= 1)
    # grid_δt  .+= (grid_tracers_diff .* grid_t_minus_final2 .* Ls ./ cp)./day_to_sec .* L
    
    # grid_t_plus .= max.(0., grid_t.-273.15)
    # grid_t_plus_final  .= ifelse.(grid_t_plus .> 0, grid_t_plus, grid_t_plus .= 0)
    # grid_t_plus_final2 .= ifelse.(grid_t_plus_final .< 0, grid_t_plus, grid_t_plus .= 1)
    # grid_δt  .+= (grid_tracers_diff .* grid_t_plus_final2 .* Lv ./ cp)./day_to_sec .* L  
    ###
    ### only water condense
    # grid_δt  .+= (grid_tracers_diff  .* Lv ./ cp)./day_to_sec .* L
    # @info maximum(grid_δt)

    ### let grid_tracers_diff be the C (condensation rate)
    ### original
    #grid_tracers_diff  .= factor3 #./ (1 .+ Lv ./ cp .* Lv .* grid_tracers_c_max ./ Rv ./ grid_t .^2) #./ (2*Δt)
    factor3 .= grid_tracers_diff
    grid_δt  .+= (grid_tracers_diff .* Lv ./ cp) ./day_to_sec .* L 
    ###
    return grid_tracers_diff
end


function e_to_qv_2D(e::Array{Float64, 2}, grid_P::Array{Float64, 2})
    # if minimum(grid_P - 0.378 * e) .< 0
    #     @info "warning: e >> P"
    #     return (e .* 0.622) ./ (grid_P)
    # end
    return (e .* 0.622) ./ (grid_P - 0.378 * e)
end

function qv_to_e_2D(qv::Array{Float64, 2}, grid_P::Array{Float64, 2})
    return (qv .* grid_P) ./ (0.378 .* qv .+ 0.622)
end


# function Spectral_Dynamics_Main()
#   # the decay of a sinusoidal disturbance to a zonally symmetric flow 
#   # that resembles that found in the upper troposphere in Northern winter.
#   name = "Spectral_Dynamics"
#   #num_fourier, nθ, nd = 63, 96, 20
#   num_fourier, nθ, nd = 42, 64, 20
#   #num_fourier, nθ, nd = 21, 32, 20
#   num_spherical = num_fourier + 1
#   nλ = 2nθ
  
#   radius = 6371000.0
#   omega = 7.292e-5
#   sea_level_ps_ref = 1.0e5
#   init_t = 264.0
  
#   # Initialize mesh
#   mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
#   θc, λc = mesh.θc,  mesh.λc
#   cosθ, sinθ = mesh.cosθ, mesh.sinθ
  
#   vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
#   # Initialize atmo_data
#   do_mass_correction = true
#   do_energy_correction = true
#   do_water_correction = false
  
#   use_virtual_temperature = false
#   atmo_data = Atmo_Data(name, nλ, nθ, nd, do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature, sinθ, radius,  omega)
  
#   # Initialize integrator
#   damping_order = 4
#   damping_coef = 1.15741e-4
#   robert_coef  = 0.04 
  
#   implicit_coef = 0.5
#   day_to_sec = 86400
#   start_time = 0
#   end_time = 2*day_to_sec  #
#   Δt = 1200
#   init_step = true
  
#   integrator = Filtered_Leapfrog(robert_coef, 
#   damping_order, damping_coef, mesh.laplacian_eig,
#   implicit_coef, Δt, init_step, start_time, end_time)
  
#   ps_ref = sea_level_ps_ref
#   t_ref = fill(300.0, nd)
#   wave_numbers = mesh.wave_numbers
#   semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data,
#   integrator, ps_ref, t_ref, wave_numbers)
  
#   # Initialize data
#   dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)
  
  
#   NT = Int64(end_time / Δt)
  
#   Get_Topography!(dyn_data.grid_geopots)
  
#   Spectral_Initialize_Fields!(mesh, atmo_data, vert_coord, sea_level_ps_ref, init_t,
#   dyn_data.grid_geopots, dyn_data)
  

#   Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#   Update_Init_Step!(semi_implicit)
#   integrator.time += Δt 
#   for i = 2:NT

#     Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#     integrator.time += Δt
#     @info integrator.time

#   end
  
# end


# #Spectral_Dynamics_Main()

