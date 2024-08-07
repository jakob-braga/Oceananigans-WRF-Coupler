------------------
-- Dependancies --
------------------

Glob
NCDatasets
Oceananigans
Printf
Random
SeawaterPolynomials


--------------
-- Overview --
--------------

I wrote this like an object oriented program because thats all I know, given that julia does not
support this it means I defined some global variables within namespaces and treated the modules
as a class, I am sorry.

This is the general idea of what each module/"class" does

Parameters.jl - user input simulation details, initial conditions, and non shared boundary conditions
WRFManager.jl - handles all WRF related things from file management to storing the wind stress
OceananigansManager.jl - handles all Oceananigans related things, and also stores the model and writer
Coupler.jl - the script that is run to run both models


-----------------
-- Quick Start --
-----------------

1) edit Parameters.jl to set user defined fields*
2) run the Coupler.jl script
3) oceananigans outputs are in ./ocean_output.nc and wrf outputs are in the ./wrf/ directory
4) the wrf outputs are organized as wrfoutputxxxx where xxxx is the iteration number
   you can something like nco wrfoutput* -O wrfoutput to concatinate everything to one file
   
*The WRF executable I shipped this with is picky about grid size, it likes having at least a 10x10
grid in the horizontal for each cpu core being used, if you run this on a machine with multiple cores
you may get an error with wrf or mpi, try running the script with a single core or up the grid resolution


