module OceananigansManager
    
    using NCDatasets
    using Printf
    
    using Oceananigans
    using Oceananigans: write_output!
    using SeawaterPolynomials
    
    import ..OceanParameters as OP
    
    
    global model
    global writer

    ###############################
    ## Managing Shared Variables ##
    ###############################
    function SST()
        return model.tracers.T[1:OP.Nx,1:OP.Ny,OP.Nz]
    end
    
    #####################
    ## Model Specifics ##
    #####################
    
    function construct_model(_u_stress::Function, _v_stress::Function)
        ## grid
        grid = RectilinearGrid(
        size = (OP.Nx, OP.Ny, OP.Nz),
        x = (0, OP.Lx),
        y = (0, OP.Ly),
        z = (-OP.Lz, 0)
        )

        ## boundary conditions
        T_bcs = FieldBoundaryConditions(
            top = FluxBoundaryCondition(OP.top_temperature_flux),
            bottom = GradientBoundaryCondition(OP.bottom_temperature_flux)
        )
        u_bcs = FieldBoundaryConditions(
            top = FluxBoundaryCondition(_u_stress)
        )
        v_bcs = FieldBoundaryConditions(
            top = FluxBoundaryCondition(_v_stress)
        )

        ## buoyancy
        eos = SeawaterPolynomials.TEOS10EquationOfState()
        buoyancy = SeawaterBuoyancy(equation_of_state=eos)

        ## actual model
        global model = NonhydrostaticModel(;
            grid,
            #buoyancy,
            advection = UpwindBiasedFifthOrder(),
            timestepper = :RungeKutta3,
            tracers = (:T, :S),
            coriolis = FPlane(f=1e-4),
            closure = AnisotropicMinimumDissipation(),
            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs)
        )

        ## initial conditions
        set!(model, u=OP.initial_u, T=OP.initial_temperature)
    end


    ############
    ## Writer ##
    ############
    function construct_writer(_filename::String)
        global writer = NetCDFOutputWriter(
                model, 
                merge(model.velocities, model.tracers),
                filename = _filename,
                schedule = IterationInterval(1),
                overwrite_existing = true
        )
        ## also write the initial state
        write_output!(writer, model)
    end


    ###################
    ## Time Stepping ##
    ###################

    @inline function time_step()
        time_step!(model, OP.dt)
        write_output!(writer, model)
    end

end