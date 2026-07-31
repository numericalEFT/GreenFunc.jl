"""
    Gf(ma::MeshArray; target_ndim=0, statistics, component=:generic,
       temperature=nothing, target_labels=nothing, metadata=NamedTuple())

Semantic wrapper around a [`MeshArray`](@ref). The first `target_ndim` axes are
target-space axes; all remaining axes are physical meshes. No target axes are
inferred from an existing `MeshArray`.
"""
struct Gf{A<:MeshArray,T<:TemperatureRegime,L,M<:NamedTuple}
    array::A
    target_ndim::Int
    statistics::Bool
    component::Symbol
    temperature::T
    target_labels::L
    metadata::M
end

_target_shape(gf::Gf) = ntuple(i -> size(gf.array, i), gf.target_ndim)
_physical_meshes(gf::Gf) = gf.array.mesh[(gf.target_ndim + 1):end]

function _validate_target_labels(labels, target_shape)
    isnothing(labels) && return nothing
    labels isa Tuple || throw(ArgumentError(
        "target_labels must be a tuple containing one label collection per target axis",
    ))
    length(labels) == length(target_shape) || throw(DimensionMismatch(
        "target_labels has $(length(labels)) axes, expected $(length(target_shape))",
    ))
    for (axis, (axis_labels, n)) in enumerate(zip(labels, target_shape))
        applicable(length, axis_labels) ||
            throw(ArgumentError("target labels for axis $axis must have a length"))
        length(axis_labels) == n || throw(DimensionMismatch(
            "target labels for axis $axis have length $(length(axis_labels)), expected $n",
        ))
    end
    return labels
end

function _validate_mesh_statistics(meshes, statistics)
    for mesh in meshes
        if mesh isa Union{ImTime,ImFreq,DLRFreq}
            mesh.isFermi == statistics || throw(ArgumentError(
                "statistics does not match temporal mesh $(typeof(mesh))",
            ))
        end
    end
    return nothing
end

const _REAL_AXIS_COMPONENTS =
    (:retarded, :advanced, :lesser, :greater, :keldysh, :spectral)

function _validate_component_meshes(meshes, component::Symbol)
    has_real_axis = any(mesh -> mesh isa Union{ReFreq,ReTime}, meshes)
    has_imaginary_axis = any(mesh -> mesh isa Union{ImTime,ImFreq,DLRFreq}, meshes)

    component in _REAL_AXIS_COMPONENTS && !has_real_axis && throw(ArgumentError(
        "component=$component requires a ReFreq or ReTime physical mesh",
    ))
    component === :matsubara && !has_imaginary_axis && throw(ArgumentError(
        "component=:matsubara requires an ImTime, ImFreq, or DLRFreq physical mesh",
    ))
    has_real_axis && component === :generic && throw(ArgumentError(
        "real-axis Green's functions require an explicit component, for example :retarded or :spectral",
    ))
    return nothing
end

function Gf(ma::MeshArray;
    target_ndim::Integer=0,
    statistics::Bool,
    component::Symbol=:generic,
    temperature::Union{Nothing,TemperatureRegime}=nothing,
    target_labels=nothing,
    metadata::NamedTuple=NamedTuple())
    0 <= target_ndim < ndims(ma) || throw(ArgumentError(
        "target_ndim must be between 0 and $(ndims(ma) - 1); at least one physical mesh is required",
    ))
    eltype(ma.data) <: Number ||
        throw(ArgumentError("Gf data must have a numeric element type"))

    target_ndim = Int(target_ndim)
    labels = _validate_target_labels(target_labels,
        ntuple(i -> size(ma, i), target_ndim))
    meshes = ma.mesh[(target_ndim + 1):end]
    _validate_mesh_statistics(meshes, statistics)
    _validate_component_meshes(meshes, component)
    resolved_temperature = _resolve_temperature(meshes, temperature)
    return Gf{typeof(ma),typeof(resolved_temperature),typeof(labels),typeof(metadata)}(
        ma, target_ndim, statistics, component, resolved_temperature, labels, metadata,
    )
end

"""
    Gf(meshes...; target_shape=(), data=nothing, dtype=ComplexF64, statistics,
       component=:generic, temperature=nothing, target_labels=nothing,
       metadata=NamedTuple())

Create a Green's function with target axes followed by physical mesh axes.
For example, `target_shape=(2, 2)` creates a matrix-valued function and retains
all complex off-diagonal entries in `data`.
"""
function Gf(meshes...;
    target_shape::Tuple=(),
    data::Union{Nothing,AbstractArray}=nothing,
    dtype=ComplexF64,
    statistics::Bool,
    component::Symbol=:generic,
    temperature::Union{Nothing,TemperatureRegime}=nothing,
    target_labels=nothing,
    metadata::NamedTuple=NamedTuple())
    isempty(meshes) && throw(ArgumentError("at least one physical mesh is required"))
    all(n -> n isa Integer && n > 0, target_shape) ||
        throw(ArgumentError("target_shape entries must be positive integers"))
    target_axes = map(Base.OneTo, target_shape)
    ma = MeshArray(target_axes..., meshes...; dtype=dtype, data=data)
    return Gf(ma;
        target_ndim=length(target_shape), statistics=statistics,
        component=component, temperature=temperature,
        target_labels=target_labels, metadata=metadata,
    )
end

Base.parent(gf::Gf) = gf.array
Base.size(gf::Gf) = size(gf.array)
Base.axes(gf::Gf) = axes(gf.array)
Base.length(gf::Gf) = length(gf.array)
Base.eltype(::Type{<:Gf{A}}) where {A} = eltype(A)
Base.getindex(gf::Gf, inds...) = getindex(gf.array, inds...)
Base.setindex!(gf::Gf, value, inds...) = setindex!(gf.array, value, inds...)
Base.eachindex(gf::Gf) = eachindex(gf.array)

function Base.getproperty(gf::Gf, name::Symbol)
    name === :data && return getfield(gf, :array).data
    name === :mesh && return _physical_meshes(gf)
    name === :fullmesh && return getfield(gf, :array).mesh
    name === :dims && return getfield(gf, :array).dims
    name === :target_shape && return _target_shape(gf)
    return getfield(gf, name)
end

Base.propertynames(::Gf, private::Bool=false) = private ?
    (:array, :target_ndim, :statistics, :component, :temperature,
        :target_labels, :metadata,
        :data, :mesh, :fullmesh, :dims, :target_shape) :
    (:data, :mesh, :target_shape, :statistics, :component, :temperature,
        :target_labels, :metadata)

function Base.similar(gf::Gf, ::Type{T}=eltype(gf.array)) where {T}
    return Gf(similar(gf.array, T);
        target_ndim=gf.target_ndim, statistics=gf.statistics,
        component=gf.component, temperature=gf.temperature,
        target_labels=gf.target_labels, metadata=gf.metadata,
    )
end

function Base.copy(gf::Gf)
    return Gf(MeshArray(; mesh=gf.fullmesh, data=copy(gf.data), dtype=eltype(gf.data));
        target_ndim=gf.target_ndim, statistics=gf.statistics,
        component=gf.component, temperature=gf.temperature,
        target_labels=gf.target_labels, metadata=gf.metadata,
    )
end

Base.show(io::IO, gf::Gf) = print(io,
    "Gf with target shape $(gf.target_shape), physical mesh $(typeof(gf.mesh)), ",
    "statistics = $(gf.statistics ? :fermion : :boson), component = $(gf.component), ",
    "temperature = $(gf.temperature)",
)

function _rewrap(gf::Gf, ma::MeshArray)
    return Gf(ma;
        target_ndim=gf.target_ndim, statistics=gf.statistics,
        component=gf.component, temperature=gf.temperature,
        target_labels=gf.target_labels, metadata=gf.metadata,
    )
end

function _parent_mesh_dim(gf::Gf, dim::Integer)
    1 <= dim <= length(gf.mesh) || throw(DimensionMismatch(
        "physical mesh dimension must be between 1 and $(length(gf.mesh))",
    ))
    return gf.target_ndim + Int(dim)
end

function _apply_gf_transform(transform, gf::Gf, args...;
    dim::Union{Nothing,Integer}=nothing, kwargs...)
    ma = if isnothing(dim)
        # Preserve the MeshArray transform's own mesh-discovery default.
        transform(parent(gf), args...; kwargs...)
    else
        transform(parent(gf), args...; dim=_parent_mesh_dim(gf, dim), kwargs...)
    end
    return _rewrap(gf, ma)
end


dlr_to_imfreq(gf::Gf, args...; kwargs...) =
    _apply_gf_transform(dlr_to_imfreq, gf, args...; kwargs...)
dlr_to_imtime(gf::Gf, args...; kwargs...) =
    _apply_gf_transform(dlr_to_imtime, gf, args...; kwargs...)
imfreq_to_dlr(gf::Gf; kwargs...) =
    _apply_gf_transform(imfreq_to_dlr, gf; kwargs...)
imtime_to_dlr(gf::Gf; kwargs...) =
    _apply_gf_transform(imtime_to_dlr, gf; kwargs...)
to_dlr(gf::Gf; kwargs...) = _apply_gf_transform(to_dlr, gf; kwargs...)
to_imfreq(gf::Gf; kwargs...) = _apply_gf_transform(to_imfreq, gf; kwargs...)
to_imtime(gf::Gf; kwargs...) = _apply_gf_transform(to_imtime, gf; kwargs...)

"""
    BlockGf(pairs::Pair{Symbol,<:Gf}...)

Ordered collection of named Green's-function blocks. Blocks may have different
target shapes, but must share statistics, component, and physical meshes.
"""
struct BlockGf{N,B<:Tuple}
    names::NTuple{N,Symbol}
    blocks::B
end


function _same_temporal_family(a, b)
    a isa ImTime && return b isa ImTime
    a isa ImFreq && return b isa ImFreq
    a isa DLRFreq && return b isa DLRFreq
    a isa ReFreq && return b isa ReFreq
    a isa ReTime && return b isa ReTime
    return !(b isa TemporalGrid)
end

_same_mesh_value(a::Number, b::Number) = isapprox(a, b)
_same_mesh_value(a, b) = isequal(a, b)

function _same_mesh(a, b)
    if a isa TemporalGrid || b isa TemporalGrid
        a isa TemporalGrid && b isa TemporalGrid || return false
        _same_temporal_family(a, b) || return false
    else
        typeof(a) == typeof(b) || return false
    end
    length(a) == length(b) || return false
    for property in (:β, :Euv, :rtol, :symmetry, :isFermi)
        if hasproperty(a, property) || hasproperty(b, property)
            hasproperty(a, property) && hasproperty(b, property) || return false
            isequal(getproperty(a, property), getproperty(b, property)) || return false
        end
    end
    return all(_same_mesh_value(x, y) for (x, y) in zip(a, b))
end

function _same_physical_meshes(left::Tuple, right::Tuple)
    length(left) == length(right) || return false
    return all(_same_mesh(a, b) for (a, b) in zip(left, right))
end

function BlockGf(pairs::Pair{Symbol,<:Gf}...)
    isempty(pairs) && throw(ArgumentError("BlockGf requires at least one block"))
    names = map(first, pairs)
    allunique(names) || throw(ArgumentError("BlockGf names must be unique"))
    blocks = map(last, pairs)
    reference = first(blocks)
    for (name, block) in zip(names[2:end], blocks[2:end])
        block.statistics == reference.statistics ||
            throw(ArgumentError("block $name has inconsistent statistics"))
        block.component == reference.component ||
            throw(ArgumentError("block $name has inconsistent component"))
        block.temperature == reference.temperature ||
            throw(ArgumentError("block $name has inconsistent temperature"))
        _same_physical_meshes(block.mesh, reference.mesh) ||
            throw(ArgumentError("block $name has inconsistent physical meshes"))
    end
    return BlockGf{length(names),typeof(blocks)}(names, blocks)
end

Base.length(blocks::BlockGf) = length(blocks.blocks)
Base.keys(blocks::BlockGf) = blocks.names
Base.values(blocks::BlockGf) = blocks.blocks
Base.getindex(blocks::BlockGf, i::Integer) = blocks.blocks[i]
function Base.getindex(blocks::BlockGf, name::Symbol)
    i = findfirst(isequal(name), blocks.names)
    isnothing(i) && throw(KeyError(name))
    return blocks.blocks[i]
end
Base.iterate(blocks::BlockGf, state::Int=1) = state > length(blocks) ? nothing :
    (blocks.names[state] => blocks.blocks[state], state + 1)
Base.pairs(blocks::BlockGf) = blocks

function Base.getproperty(blocks::BlockGf, name::Symbol)
    name === :statistics && return first(getfield(blocks, :blocks)).statistics
    name === :component && return first(getfield(blocks, :blocks)).component
    name === :temperature && return first(getfield(blocks, :blocks)).temperature
    name === :mesh && return first(getfield(blocks, :blocks)).mesh
    return getfield(blocks, name)
end

Base.show(io::IO, blocks::BlockGf) = print(io,
    "BlockGf with blocks $(collect(blocks.names)), component = $(blocks.component)",
)

"""
    SpectralDensity(gf::Gf; constraint=:hermitian, atol=0, rtol=sqrt(eps()))

Validated scalar- or matrix-valued spectral density on a mesh containing
[`ReFreq`](@ref). `constraint=:hermitian` checks every scalar or complete target matrix;
`constraint=:psd` additionally checks positive semidefiniteness. Complex
off-diagonal entries are retained and checked against their conjugate partners.
"""
struct SpectralDensity{G<:Gf}
    gf::G
    constraint::Symbol
    atol::Float64
    rtol::Float64
end

function _validate_spectral_matrices(gf::Gf, constraint, atol, rtol)
    all(isfinite, gf.data) || throw(ArgumentError("spectral-density data must be finite"))
    constraint === :none && return nothing

    if gf.target_ndim == 0
        for index in eachindex(gf.data)
            value = gf.data[index]
            isapprox(value, conj(value); atol=atol, rtol=rtol) ||
                throw(ArgumentError("scalar spectral-density value at linear index $index is not real"))
            if constraint === :psd
                tolerance = atol + rtol * abs(value)
                real(value) >= -tolerance || throw(ArgumentError(
                    "scalar spectral-density value at linear index $index is negative",
                ))
            end
        end
        return nothing
    end

    physical_dims = size(gf.data)[3:end]
    for index in CartesianIndices(physical_dims)
        physical_index = Tuple(index)
        matrix = @view gf.data[:, :, physical_index...]
        isapprox(matrix, adjoint(matrix); atol=atol, rtol=rtol) ||
            throw(ArgumentError("spectral-density matrix at physical index $physical_index is not Hermitian"))
        if constraint === :psd
            scale = maximum(abs, matrix)
            tolerance = atol + rtol * scale
            minimum(eigvals(Hermitian(Matrix(matrix)))) >= -tolerance ||
                throw(ArgumentError("spectral-density matrix at physical index $physical_index is not positive semidefinite"))
        end
    end
    return nothing
end

function SpectralDensity(gf::Gf;
    constraint::Symbol=:hermitian,
    atol::Real=0,
    rtol::Real=sqrt(eps(Float64)))
    constraint in (:none, :hermitian, :psd) ||
        throw(ArgumentError("constraint must be :none, :hermitian, or :psd"))
    gf.component === :spectral ||
        throw(ArgumentError("a SpectralDensity requires component=:spectral"))
    atol >= 0 || throw(ArgumentError("atol must be nonnegative"))
    rtol >= 0 || throw(ArgumentError("rtol must be nonnegative"))
    any(mesh -> mesh isa ReFreq, gf.mesh) ||
        throw(ArgumentError("a SpectralDensity requires a ReFreq physical mesh"))
    gf.target_ndim in (0, 2) || throw(ArgumentError(
        "a SpectralDensity requires a scalar or matrix target",
    ))
    gf.target_ndim == 0 || gf.target_shape[1] == gf.target_shape[2] ||
        throw(DimensionMismatch("spectral-density target matrix must be square"))

    _validate_spectral_matrices(gf, constraint, atol, rtol)
    return SpectralDensity{typeof(gf)}(gf, constraint, Float64(atol), Float64(rtol))
end

Base.parent(spectral::SpectralDensity) = spectral.gf
Base.size(spectral::SpectralDensity) = size(spectral.gf)
Base.length(spectral::SpectralDensity) = length(spectral.gf)
Base.getindex(spectral::SpectralDensity, inds...) = getindex(spectral.gf, inds...)

function Base.getproperty(spectral::SpectralDensity, name::Symbol)
    name in (:data, :mesh, :target_shape, :statistics, :component,
        :temperature, :target_labels, :metadata) &&
        return getproperty(getfield(spectral, :gf), name)
    return getfield(spectral, name)
end

Base.show(io::IO, spectral::SpectralDensity) = print(io,
    "SpectralDensity with target shape $(spectral.target_shape), constraint = $(spectral.constraint)",
)
