"""
    TemperatureRegime

Thermodynamic context carried by a Green's function independently of its
physical mesh. Real-time and real-frequency meshes do not determine a
temperature, so their Green's functions must choose a regime explicitly.
"""
abstract type TemperatureRegime end

"""
    ZeroTemperature()

Exact ground-state (`T = 0`) Green-function context. This is not represented
as `FiniteTemperature(Inf)` because Matsubara grids and thermal kernels do not
have a regular `β = Inf` representation.
"""
struct ZeroTemperature <: TemperatureRegime end

"""
    FiniteTemperature(β)

Finite-temperature equilibrium context at inverse temperature `β > 0`.
"""
struct FiniteTemperature{T<:Real} <: TemperatureRegime
    β::T
    function FiniteTemperature(β::T) where {T<:Real}
        β isa Bool && throw(ArgumentError("inverse temperature β must be a real number, not Bool"))
        isfinite(β) || throw(ArgumentError("inverse temperature β must be finite"))
        β > zero(β) || throw(ArgumentError("inverse temperature β must be positive"))
        return new{T}(β)
    end
end

inverse_temperature(::ZeroTemperature) = nothing
inverse_temperature(temperature::FiniteTemperature) = temperature.β

Base.show(io::IO, ::ZeroTemperature) = print(io, "ZeroTemperature()")
Base.show(io::IO, temperature::FiniteTemperature) =
    print(io, "FiniteTemperature(", temperature.β, ")")

Base.:(==)(::ZeroTemperature, ::ZeroTemperature) = true
Base.:(==)(left::FiniteTemperature, right::FiniteTemperature) = left.β == right.β
Base.isequal(::ZeroTemperature, ::ZeroTemperature) = true
Base.isequal(left::FiniteTemperature, right::FiniteTemperature) =
    isequal(left.β, right.β)
Base.hash(::ZeroTemperature, seed::UInt) = hash(:zero_temperature, seed)
Base.hash(temperature::FiniteTemperature, seed::UInt) =
    hash((:finite_temperature, temperature.β), seed)

function _temperature_from_meshes(meshes)
    imaginary_meshes = filter(mesh -> mesh isa Union{ImTime,ImFreq,DLRFreq}, meshes)
    isempty(imaginary_meshes) && return nothing

    β = first(imaginary_meshes).β
    for mesh in Iterators.drop(imaginary_meshes, 1)
        isapprox(mesh.β, β; atol=0, rtol=sqrt(eps(Float64))) ||
            throw(ArgumentError("all imaginary-time/frequency meshes must have the same inverse temperature"))
    end
    return FiniteTemperature(β)
end

function _resolve_temperature(meshes, temperature::Union{Nothing,TemperatureRegime})
    mesh_temperature = _temperature_from_meshes(meshes)
    if isnothing(mesh_temperature)
        isnothing(temperature) && throw(ArgumentError(
            "temperature is required when no imaginary-time/frequency mesh fixes β; " *
            "pass ZeroTemperature() or FiniteTemperature(β)",
        ))
        return temperature
    end

    isnothing(temperature) && return mesh_temperature
    temperature isa ZeroTemperature && throw(ArgumentError(
        "imaginary-time/frequency meshes are finite-temperature domains and cannot use ZeroTemperature()",
    ))
    temperature isa FiniteTemperature || throw(ArgumentError(
        "imaginary-time/frequency meshes require a FiniteTemperature context",
    ))
    isapprox(temperature.β, mesh_temperature.β;
        atol=0, rtol=sqrt(eps(Float64))) || throw(ArgumentError(
        "temperature β=$(temperature.β) does not match imaginary mesh β=$(mesh_temperature.β)",
    ))
    return temperature
end
