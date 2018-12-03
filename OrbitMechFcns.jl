"""
THINGS TO LOOK UP

ode solvers -> define time array
Make structure non-immutable -> mutable struct
Rotation Matricies ? made my own for now
plotting some of this stuff (EarthPlot, ground tracks ...)


FUNCTIONS INCLUDED
    * Orbital Element Structure
    * OE2ECI
    * GravitationalForce
    * Rotation Matricies
    * ECI2OE
    * anom2E, E2anom
    * E2M, M2E
STILL TO ADD
    * Timing (GMST,UT1,...)
    * ECI2ECEF
    * ECEF2GEOCEN/GEODED
    * SatPropagator
    * Lambert solver ?..

    * transit 290 function ...

Look into
    * plotting (earthPlot, ground tracks,M_Map package)
    * Other ODE solvers

"""

# Packages that are useful
using LinearAlgebra
using DifferentialEquations
using Plots
using CSV

# Read in all functions from their respective files
include("anomalytransformations.jl")
include("OEstruct.jl")
include("RotMat.jl")
include("GravitationalForcePropagator.jl")
include("OE_ECItransformations.jl")
include("TimeFunctions.jl")
include("MeanMotionProp.jl")
include("EarthGroundPlot.jl")
include("GeocentricConversion.jl")
include("GroundRange.jl")
