module WRFManager
    
    using Glob
    using NCDatasets
    
    import ..OceanParameters as OP
    
    ## useful constants
    const wrfin_filename = "wrfinput_d01"
    const wrfrst_filename = "wrfrst_d01_0001-01-01_00:00:00"
    const wrfout_filename = "wrfout_d01_0001-01-01_00:00:00"

    ###############################
    ## Managing Shared Variables ##
    ###############################

    @inline function get_U10()
        wrf_data = NCDataset(string("./wrf/", wrfrst_filename), "r")
        tmp_U10 = wrf_data["U10"][:, :, 1]
        close(wrf_data)
        return tmp_U10
    end

    @inline function get_V10()
        wrf_data = NCDataset(string("./wrf/", wrfrst_filename), "r")
        tmp_V10 = wrf_data["V10"][:, :, 1]
        close(wrf_data)
        return tmp_V10
    end

    @inline function get_UST()
        wrf_data = NCDataset(string("./wrf/", wrfrst_filename), "r")
        tmp_UST = wrf_data["UST"][:, :, 1]
        close(wrf_data)
        return tmp_UST
    end

    @inline function get_HFX()
        wrf_data = NCDataset(string("./wrf/", wrfrst_filename), "r")
        tmp_HFX = wrf_data["HFX"][:, :, 1]
        close(wrf_data)
        return tmp_HFX
    end
    
    ## TODO clean this up so its not just returning an array, this should be a struct
    @inline function get_coupled_vars()
        wrf_data = NCDataset(string("./wrf/", wrfrst_filename), "r")
        tmp_U10 = wrf_data["U10"][:, :, 1]
        tmp_V10 = wrf_data["V10"][:, :, 1]
        tmp_UST = wrf_data["UST"][:, :, 1]
        tmp_HFX = wrf_data["HFX"][:, :, 1]
        close(wrf_data)
        return (tmp_U10, tmp_V10, tmp_UST, tmp_HFX)
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
        ## construct surface velocities on the ocean based off the initializations in the parameters config
    end

    @inline function run_wrf()
        cd("./wrf/")
        run(`mpirun ./wrf.exe`)
        cd("..")
    end




end
