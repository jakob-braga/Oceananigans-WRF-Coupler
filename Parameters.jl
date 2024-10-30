using Random

module OceanParameters
    #=
    --------------------------
    - domain, grid, and time -
    --------------------------
    =# 
    ## grid
    const Nx = 512
    const Ny = 64
    const Nz = 64
    ## domain [m]
    const Lx = 307200
    const dx = Lx/Nx
    const Ly = 38400
    const dy = Ly/Ny
    const Lz = 1024
    const dz = Lz/Nz
    ## time stepping [s]
    const Nt = 5
    const dt = 60

    #=
    -------------
    - constants -
    -------------
    =#
    ## density
    const ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
    const N2 = (3.7e-3)^2
    const g = 9.8
    ## temperature
    const cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater
    const dTdz = 0.01 # K m⁻¹
    const d_theta = 2.0
    ## wind stress
    const cᴰ = 1.6e-3 # dimensionless drag coefficient
    const ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

    #=
    ---------------
    - BCs and ICs -
    ---------------
    =#
    const cj = 6
    const alpha = 2e-4
    const ellx = 38400*3
    const ellz = 1024
    ## temperature
    @inline function bottom_temperature_flux(x, y, t)
        return 0.0
    end
    @inline function initial_temperature(x, y, z)
        return 300 + d_theta*sech(cj/ellx*(x-Lx/2))^2*sech(cj/ellz*z)^2 + N2/(g*alpha)*z
    end
    ## velocities [m/s]
    @inline function initial_u(x, y, z)
        return 1e-4*randn()
    end
    @inline function initial_v(x, y, z)
        coeff = -2*g*alpha*d_theta*ellz/(1e-4*ellx)
        x_ = cj/ellx*(x-Lx/2)
        z_ = cj/ellz*z
        return coeff*sech(x_)^2*tanh(x_)*(tanh(z_)+1) + 1e-3*randn()
    end
end

module AtmospshereParameters
    ## TODO add wrf specific parameters here like background winds, temperature feild, etc.
end
