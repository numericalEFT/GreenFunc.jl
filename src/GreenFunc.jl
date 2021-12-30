module GreenFunc
using StaticArrays, Lehmann, CompositeGrids
# Write your package code here.

abstract type TimeDomain end
abstract type ImTime <: TimeDomain end
abstract type ReTime <: TimeDomain end
abstract type ImFreq <: TimeDomain end
abstract type ReFreq <: TimeDomain end
abstract type DLRFreq <: TimeDomain end

export TimeDomain, ImTime, ReTime, ImFreq, ReFreq, DLRFreq

include("green/Green.jl")
export Green2DLR, toTau, toMatFreq, toDLR, get



end
