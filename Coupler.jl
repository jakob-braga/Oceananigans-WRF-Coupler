using Glob
using NCDatasets
using Printf

using Oceananigans
using Oceananigans: write_output!
using Oceananigans.Units: minute, minutes, hour
using SeawaterPolynomials

include("./Parameters.jl")
include("./OceananigansManager.jl")
include("./WRFManager.jl")

import .OceanParameters as OP
import .WRFManager as WFM
import .OceananigansManager as OCM


########################
## WRF Initialization ##
########################

print("--Initializing Wrf--\n")

## remove all the possible wrfrst files from previous runs
WFM.clean_directory()
## initialize wrf based on Parameters.jl
WFM.initialize_wrf()

print("--WRF Initialization Complete--\n\n")

################################
## Ocean Model Initialization ##
################################

print("--Initializing Oceananigans Model--")

## construct the model based on WRF wind stress
OCM.construct_model(WFM.u_stress, WFM.v_stress)
## construct writer
OCM.construct_writer("./ocean_output.nc")

print("--Initialization of Oceananigans Model Complete--\n\n")

###################
## Time Stepping ##
###################

for n = 1:OP.Nt
    ## run the ocean model
    OCM.time_step()
    ## update WRF variables
    WFM.edit_wrfrst("TSK", OCM.SST())
    ## run wrf
    WFM.run_wrf()
    WFM.rename_outputs(_iteration=n)
    ## update wrf manager and wind stress
    WFM.update_data()
end
