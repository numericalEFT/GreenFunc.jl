"""
    thermal_distribution(energy, statistics, temperature)

Equilibrium Fermi or Bose occupation at an energy measured relative to the
chemical potential. `ZeroTemperature()` is evaluated as the exact `β -> Inf`
limit without constructing a Matsubara mesh or storing `β = Inf`.

At a fermionic step exactly at zero the symmetric value `1/2` is returned. The
bosonic distribution is singular there and throws a `DomainError`.
"""
function thermal_distribution(energy::Real, statistics::Bool, ::ZeroTemperature)
    isfinite(energy) || throw(ArgumentError("energy must be finite"))
    if statistics === FERMION
        iszero(energy) && return 0.5
        return energy < zero(energy) ? 1.0 : 0.0
    end
    iszero(energy) && throw(DomainError(energy,
        "the zero-temperature Bose distribution is singular at zero energy"))
    return energy < zero(energy) ? -1.0 : 0.0
end

function thermal_distribution(energy::Real, statistics::Bool,
                              temperature::FiniteTemperature)
    isfinite(energy) || throw(ArgumentError("energy must be finite"))
    exponent = temperature.β * energy
    if statistics === FERMION
        return exponent >= zero(exponent) ?
            exp(-exponent) / (1 + exp(-exponent)) :
            1 / (1 + exp(exponent))
    end
    iszero(exponent) && throw(DomainError(exponent,
        "the Bose distribution is singular at zero energy"))
    return inv(expm1(exponent))
end

function _one_plus_exchange_distribution(energy::Real, statistics::Bool,
                                         ::ZeroTemperature)
    if statistics === FERMION
        iszero(energy) && return 0.5
        return energy < zero(energy) ? 0.0 : 1.0
    end
    iszero(energy) && throw(DomainError(energy,
        "the zero-temperature Bose distribution is singular at zero energy"))
    return energy < zero(energy) ? 0.0 : 1.0
end

function _one_plus_exchange_distribution(energy::Real, statistics::Bool,
                                         temperature::FiniteTemperature)
    exponent = temperature.β * energy
    if statistics === FERMION
        return exponent >= zero(exponent) ?
            1 / (1 + exp(-exponent)) :
            exp(exponent) / (1 + exp(exponent))
    end
    iszero(exponent) && throw(DomainError(exponent,
        "the Bose distribution is singular at zero energy"))
    return -inv(expm1(-exponent))
end

function _single_refreq_axis(gf::Gf)
    axes = findall(mesh -> mesh isa ReFreq, gf.mesh)
    length(axes) == 1 || throw(ArgumentError(
        "equilibrium spectral relations require exactly one ReFreq physical mesh",
    ))
    return only(axes)
end

function _scale_physical_axis(data, factors, axis::Int)
    shape = ntuple(dim -> dim == axis ? length(factors) : 1, ndims(data))
    return data .* reshape(factors, shape)
end

function _real_axis_gf(gf::Gf, data, component::Symbol)
    ma = MeshArray(; mesh=gf.fullmesh, data=data, dtype=eltype(data))
    return Gf(ma;
        target_ndim=gf.target_ndim,
        statistics=gf.statistics,
        component=component,
        temperature=gf.temperature,
        target_labels=gf.target_labels,
        metadata=gf.metadata,
    )
end

function _target_adjoint_data(gf::Gf)
    gf.target_ndim in (0, 2) || throw(ArgumentError(
        "retarded/advanced conversion requires a scalar or matrix target",
    ))
    gf.target_ndim == 0 && return conj.(gf.data)
    gf.target_shape[1] == gf.target_shape[2] || throw(DimensionMismatch(
        "retarded/advanced conversion requires a square target matrix",
    ))

    output = similar(gf.data)
    physical_dims = size(gf.data)[3:end]
    for index in CartesianIndices(physical_dims)
        physical_index = Tuple(index)
        output[:, :, physical_index...] .=
            adjoint(@view gf.data[:, :, physical_index...])
    end
    return output
end

"""
    advanced_from_retarded(retarded) -> Gf

Construct the advanced component pointwise as the target-space adjoint of a
retarded real-frequency Green's function. Scalar targets are complex
conjugated; matrix targets are conjugate-transposed at every physical point.
"""
function advanced_from_retarded(retarded::Gf)
    retarded.component === :retarded || throw(ArgumentError(
        "advanced_from_retarded requires component=:retarded",
    ))
    any(mesh -> mesh isa ReFreq, retarded.mesh) || throw(ArgumentError(
        "advanced_from_retarded requires a ReFreq physical mesh",
    ))
    return _real_axis_gf(retarded, _target_adjoint_data(retarded), :advanced)
end

"""
    spectral_from_retarded(retarded; constraint=:hermitian, atol=0,
                           rtol=sqrt(eps())) -> SpectralDensity

Construct the normalized spectral density
`ρ(ω) = i [G^R(ω) - G^A(ω)] / (2π)`. For matrix targets the adjoint is taken in
target space; this deliberately differs from elementwise `-imag.(G^R)/π` for
complex off-diagonal entries.
"""
function spectral_from_retarded(retarded::Gf;
                                constraint::Symbol=:hermitian,
                                atol::Real=0,
                                rtol::Real=sqrt(eps(Float64)))
    advanced = advanced_from_retarded(retarded)
    data = (im / (2π)) .* (retarded.data .- advanced.data)
    spectral = _real_axis_gf(retarded, data, :spectral)
    return SpectralDensity(spectral; constraint=constraint, atol=atol, rtol=rtol)
end

function _component_from_spectral(spectral::SpectralDensity, component::Symbol)
    gf = parent(spectral)
    physical_axis = _single_refreq_axis(gf)
    frequencies = collect(gf.mesh[physical_axis])
    exchange_sign = gf.statistics === FERMION ? -1 : 1

    factors = if component === :lesser
        occupations = thermal_distribution.(frequencies, Ref(gf.statistics),
            Ref(gf.temperature))
        -2π * im .* exchange_sign .* occupations
    elseif component === :greater
        weights = _one_plus_exchange_distribution.(frequencies,
            Ref(gf.statistics), Ref(gf.temperature))
        -2π * im .* weights
    elseif component === :keldysh
        occupations = thermal_distribution.(frequencies, Ref(gf.statistics),
            Ref(gf.temperature))
        -2π * im .* (1 .+ 2 .* exchange_sign .* occupations)
    else
        throw(ArgumentError("unsupported equilibrium component $component"))
    end

    data_axis = gf.target_ndim + physical_axis
    data = _scale_physical_axis(gf.data, factors, data_axis)
    return _real_axis_gf(gf, data, component)
end

"""
    lesser_from_spectral(spectral) -> Gf

Construct the equilibrium lesser component from a spectral density `ρ(ω)`.
The convention is `G<(ω) = -2πi ξ n_ξ(ω) ρ(ω)`, where `ξ=-1` for
fermions and `ξ=+1` for bosons.
"""
lesser_from_spectral(spectral::SpectralDensity) =
    _component_from_spectral(spectral, :lesser)

"""
    greater_from_spectral(spectral) -> Gf

Construct the equilibrium greater component using
`G>(ω) = -2πi [1 + ξ n_ξ(ω)] ρ(ω)`.
"""
greater_from_spectral(spectral::SpectralDensity) =
    _component_from_spectral(spectral, :greater)

"""
    keldysh_from_spectral(spectral) -> Gf

Construct `G^K = G^> + G^<` from an equilibrium spectral density.
"""
keldysh_from_spectral(spectral::SpectralDensity) =
    _component_from_spectral(spectral, :keldysh)

"""
    equilibrium_components(spectral)

Return the equilibrium `lesser`, `greater`, and `keldysh` components generated
from `spectral` under its exact zero- or finite-temperature context.
"""
function equilibrium_components(spectral::SpectralDensity)
    return (
        lesser=lesser_from_spectral(spectral),
        greater=greater_from_spectral(spectral),
        keldysh=keldysh_from_spectral(spectral),
    )
end
