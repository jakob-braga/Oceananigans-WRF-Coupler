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


###########################################
## Managing Significant Shared Variables ##
###########################################
mutable struct CoupledVars
    U10
    V10
    UST
    HFX
    CoupledVars() = new()
end

function initialize_coupled_vars!(_vars)
    _vars.U10 = zeros((OP.Nx, OP.Ny))
    _vars.V10 = zeros((OP.Nx, OP.Ny))
    _vars.UST = zeros((OP.Nx, OP.Ny))
    _vars.HFX = zeros((OP.Nx, OP.Ny))
end

function update_coupled_vars!(_vars)
    _vars.U10 = WFM.get_U10()
    _vars.V10 = WFM.get_V10()
    _vars.UST = WFM.get_UST()
    _vars.HFX = WFM.get_HFX()
end

coupled_vars = CoupledVars()
initialize_coupled_vars!(coupled_vars)

####################################################
## Defining Boundary Conditions from Coupled Vars ##
####################################################

@inline function u_stress(x, y, z)
    i = max(1, min(OP.Nx, Int(trunc(x/OP.dx))))
    j = max(1, min(OP.Ny, Int(trunc(y/OP.dy))))
    u10 = coupled_vars.U10[i, j]
    v10 = coupled_vars.V10[i, j]
    ust = coupled_vars.UST[i, j]
    return -OP.ρₐ/OP.ρₒ*ust^2*u10/sqrt(u10^2 + v10^2)
end

@inline function v_stress(x, y, t)
    i = max(1, min(OP.Nx, Int(trunc(x/OP.dx))))
    j = max(1, min(OP.Ny, Int(trunc(y/OP.dy))))
    u10 = coupled_vars.U10[i, j]
    v10 = coupled_vars.V10[i, j]
    ust = coupled_vars.UST[i, j]
    return -OP.ρₐ/OP.ρₒ*ust^2*v10/sqrt(u10^2 + v10^2)
end

@inline function T_flux(x, y, t)
    i = max(1, min(OP.Nx, Int(trunc(x/OP.dx))))
    j = max(1, min(OP.Ny, Int(trunc(y/OP.dy))))
    hfx = coupled_vars.HFX[i, j]
    return hfx/(OP.ρₒ*OP.cᴾ)
end

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

print("--Initializing Oceananigans Model--\n")

## construct the model based on coupled boundary conditions
OCM.construct_model(
    u_stress,
    v_stress,
    T_flux
)
## construct writer
OCM.construct_writer("./ocean_output.nc")

print("--Initialization of Oceananigans Model Complete--\n\n")

###################
## Time Stepping ##
###################

for n = 1:OP.Nt
    ## give a progress message every 10 frames
    if mod(n, 10) == 0
        @printf("Iteration: %04d", n)
    end
    ## run the ocean model
    OCM.time_step()
    ## update WRF variables
    WFM.edit_wrfrst("TSK", OCM.get_SST())
    ## run wrf
    WFM.run_wrf()
    WFM.rename_outputs(_iteration=n)
    ## update coupled vars
    update_coupled_vars!(coupled_vars)
end
