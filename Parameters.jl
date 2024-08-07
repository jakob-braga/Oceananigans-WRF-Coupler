using Random

module OceanParameters
    #=
    --------------------------
    - domain, grid, and time -
    --------------------------
    =# 
    ## grid
    const Nx = 64
    const Ny = 64
    const Nz = 32
    ## domain [m]
    const Lx = 4096
    const dx = Lx/Nx
    const Ly = 4096
    const dy = Ly/Ny
    const Lz = 512
    const dz = Lz/Nz
    ## time stepping [s]
    const Nt = 20
    const dt = 10

    #=
    -------------
    - constants -
    -------------
    =#
    ## density
    const ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
    ## temperature
    const Qʰ = 200.0  # W m⁻², surface _heat_ flux
    const cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater
    const dTdz = 0.01 # K m⁻¹
    ## wind stress
    const cᴰ = 2.5e-3 # dimensionless drag coefficient
    const ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

    #=
    ---------------
    - BCs and ICs -
    ---------------
    =#
    ## temperature
    @inline function top_temperature_flux(x, y, t)
        return Qʰ / (ρₒ * cᴾ)
    end
    @inline function bottom_temperature_flux(x, y, t)
        return 0
    end
    @inline function initial_temperature(x, y, z)
        return 300 + 0.01*z + 5*sech(50/Lx*(x-Lx/2))*(tanh(10/Lz*(z))+1)
    end
    ## velocities [m/s]
    @inline function initial_u(x, y, z)
        return 0.0 #1e-6*randn()
    end
end

module AtmospshereParameters
    ## TODO add wrf specific parameters here like background winds, temperature feild, etc.
end
