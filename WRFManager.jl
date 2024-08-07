module WRFManager
    
    using Glob
    using NCDatasets
    
    import ..OceanParameters as OP
    
    ## useful constants
    const wrfin_filename = "wrfinput_d01"
    const wrfrst_filename = "wrfrst_d01_0001-01-01_00:00:00"
    const wrfout_filename = "wrfout_d01_0001-01-01_00:00:00"
    
    ## dynamic variables, fancy (hacky) way to define a global variable
    ## so i can use this namespace as a class
    global const U10 = Ref{Array{Float32, 3}}(zeros((1, 1, 1)))
    global const V10 = Ref{Array{Float32, 3}}(zeros((1, 1, 1)))


    ###############################
    ## Managing Shared Variables ##
    ###############################
    
    @inline function update_data()
        wrf_data = NCDataset(string("./wrf/", wrfrst_filename), "r")
        global U10[] = wrf_data["U10"][:]
        global V10[] = wrf_data["V10"][:]
        close(wrf_data)
    end

    @inline function u_stress(x, y, t)
        i = max(1, min(OP.Nx, Int(trunc(x/OP.dx))))
        j = max(1, min(OP.Ny, Int(trunc(y/OP.dy))))
        tmp_u10 = U10[][i, j, 1]
        return -OP.ρₐ/OP.ρₒ*OP.cᴰ*tmp_u10*abs(tmp_u10)
    end

    @inline function v_stress(x, y, t)
        i = max(1, min(OP.Nx, Int(trunc(x/OP.dx))))
        j = max(1, min(OP.Ny, Int(trunc(y/OP.dy))))
        tmp_v10 = V10[][i, j, 1]
        return -OP.ρₐ/OP.ρₒ*OP.cᴰ*tmp_v10*abs(tmp_v10)
    end


    ################################
    ## editing wrf specific files ##
    ################################

    @inline function edit_namelist(_variable::String, _value::String)
        (tmppath, tmpio) = mktemp()
        open("./wrf/namelist.input") do io
            for line in eachline(io, keep=true) # keep so the new line isn't chomped
                if length(line) >= length(_variable)+2 # +2 for spaces
                    if occursin(" "*_variable*" ", line)
                        eq_idx = findfirst("=", line)
                        line = replace(line, line[eq_idx[1]:end] => "= "*_value*",\n")
                    end
                end
                write(tmpio, line)
            end
        end
        close(tmpio)
        mv(tmppath, "./wrf/namelist.input", force=true)
    end
    
    @inline function edit_wrfin(_variable::String, _initial_condition::Function)
        wrfin_data = NCDataset(string("./wrf/", wrfin_filename), "a")
        wrfin_size = size(wrfin_data[_variable])
        dx = wrfin_data.attrib["DX"]
        dy = wrfin_data.attrib["DY"]
        initial_val = zeros(Float32, wrfin_size)
        for i = 1:wrfin_size[1]
            for j = 1:wrfin_size[2]
                x = (i-1)*dx
                y = (j-1)*dy
                z = 0
                initial_val[i, j, 1] = _initial_condition(x, y, z)
            end
        end
        wrfin_data[_variable][:] = initial_val
        close(wrfin_data)
    end
    
    @inline function edit_wrfrst(_variable::String, _data)
        wrfrst_data = NCDataset(string("./wrf/", wrfrst_filename), "a")
        wrfrst_data[_variable][:, :, 1] = _data
        close(wrfrst_data)
    end
    

    ####################################
    ## Cleaning up and renaming files ##
    ####################################

    @inline function rename_outputs(; _iteration::Integer = 0)
        ## depending on the iteration you will get the first or second iteration
        if _iteration == 0
            rst_idx = 1
        else
            rst_idx = 2
        end
        mv(glob("./wrf/wrfout_*")[1], string("./wrf/wrfoutput", lpad(_iteration, 4, "0")), force=true)
        mv(glob("./wrf/wrfrst*")[rst_idx], string("./wrf/", wrfrst_filename), force=true)
    end

    @inline function clean_directory()
        wrfrst_files = glob("./wrf/wrfrst*")
        wrfout_files = glob("./wrf/wrfout_*")
        foreach(f -> rm(f), wrfrst_files)
        foreach(f -> rm(f), wrfout_files)
    end
    
    
    #################
    ## Running WRF ##
    #################

    @inline function initialize_wrf()
        ## edit the namelist in wrf
        edit_namelist("restart", ".false.")
        edit_namelist("run_seconds", string(OP.dt))
        edit_namelist("history_interval_s", string(OP.dt))
        edit_namelist("restart_interval_s", string(OP.dt))
        edit_namelist("e_we", string(OP.Nx + 1))
        edit_namelist("e_sn", string(OP.Ny + 1))
        edit_namelist("dx", string(trunc(OP.dx)))
        edit_namelist("dy", string(trunc(OP.dy)))
        ## actually run the ideal.exe
        cd("./wrf/")
        run(`./ideal.exe`)
        cd("..")
        ## set the SST (TSK) for WRF based on the ocean temperature given
        edit_wrfin("TSK", OP.initial_temperature)
        ## Run WRF for the duration of a single Ocean time step
        print("Running WRF inital time step\n")
        run_wrf()
        print("WRF initial time step complete")
        ## lastly make sure wrf is set to restart now
        edit_namelist("restart", ".true.")
        ## clean up the files
        rename_outputs(_iteration=0)
        ## based off wrfoutput we define the windstress
        update_data()
    end

    @inline function run_wrf()
        cd("./wrf/")
        run(`mpirun ./wrf.exe`)
        cd("..")
    end




end