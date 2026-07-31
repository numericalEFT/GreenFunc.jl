struct DescendingTestGrid <: CompositeGrids.AbstractGrid{Float64}
    grid::Vector{Float64}
    size::Int
end
DescendingTestGrid(points::Vector{Float64}) = DescendingTestGrid(points, length(points))
Base.size(grid::DescendingTestGrid) = size(grid.grid)
Base.getindex(grid::DescendingTestGrid, i::Int) = grid.grid[i]

@testset "Real-domain meshes" begin
    ω = ReFreq(-2.0, 2.0, 5)
    @test collect(ω) == [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test length(ω) == 5
    @test locate(ω, 0.0) == 3
    @test volume(ω) ≈ sum(volume(ω, i) for i in eachindex(ω))

    ωrev = ReFreq([2.0, 1.0, 0.0, -1.0])
    @test collect(ωrev) == [2.0, 1.0, 0.0, -1.0]
    @test collect(ωrev.grid) == [-1.0, 0.0, 1.0, 2.0]

    reversed_grid = DescendingTestGrid([2.0, 1.0, 0.0])
    ωgridrev = ReFreq(reversed_grid)
    @test collect(ωgridrev) == [2.0, 1.0, 0.0]
    @test collect(ωgridrev.grid) == [0.0, 1.0, 2.0]

    t = ReTime(; window=(0.0, 1.0), n_t=3)
    @test collect(t) == [0.0, 0.5, 1.0]
    @test ReFreq(; window=(-1.0, 1.0), n_w=3)[2] == 0.0
    @test occursin("-5.0", sprint(show, ReFreq(-5.0, -1.0, 5)))

    @test_throws ArgumentError ReFreq([0.0, 0.0, 1.0])
    @test_throws ArgumentError ReTime([0.0, 0.5, 0.5])
    @test_throws ArgumentError ReFreq([0.0, 2.0, 1.0])
end

@testset "Green-function semantic wrappers" begin
    mesh = ReFreq(-1.0, 1.0, 3)
    scalar_ma = MeshArray(mesh; dtype=ComplexF64, data=ComplexF64[1, 2, 3])
    scalar = Gf(scalar_ma; statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    @test scalar.target_ndim == 0
    @test scalar.target_shape == ()
    @test scalar.mesh == (mesh,)
    @test parent(scalar) === scalar_ma

    data = zeros(ComplexF64, 2, 2, length(mesh))
    data[1, 2, :] .= 1 .+ 2im
    data[2, 1, :] .= 3 .- 4im
    matrix_gf = Gf(mesh;
        target_shape=(2, 2), data=data, statistics=FERMION,
        component=:retarded, temperature=ZeroTemperature(),
        target_labels=((:a, :b), (:a, :b)),
        metadata=(source=:test,),
    )
    @test matrix_gf.target_ndim == 2
    @test matrix_gf.target_shape == (2, 2)
    @test matrix_gf.fullmesh[1:2] == (Base.OneTo(2), Base.OneTo(2))
    @test matrix_gf.data[1, 2, 2] == 1 + 2im
    @test matrix_gf.data[2, 1, 2] == 3 - 4im
    @test matrix_gf.metadata.source == :test
    for preserved in (copy(matrix_gf), similar(matrix_gf))
        @test preserved.statistics == matrix_gf.statistics
        @test preserved.component == matrix_gf.component
        @test preserved.temperature isa ZeroTemperature
        @test preserved.target_labels == matrix_gf.target_labels
        @test preserved.metadata == matrix_gf.metadata
    end
    @test copy(matrix_gf).data == matrix_gf.data

    @test_throws DimensionMismatch Gf(mesh;
        target_shape=(2, 2), statistics=FERMION,
        component=:retarded, temperature=ZeroTemperature(),
        target_labels=((:a,), (:a, :b)),
    )
    imtime = ImTime(2.0, FERMION; grid=[0.0, 1.0, 2.0])
    @test_throws ArgumentError Gf(imtime; statistics=BOSON)

    imfreq = ImFreq(10.0, FERMION)
    imfreq_data = reshape(
        ComplexF64[inv(im * ω - 1) for ω in imfreq],
        1, 1, length(imfreq),
    )
    dlr_gf = to_dlr(Gf(imfreq;
        target_shape=(1, 1), data=imfreq_data, statistics=FERMION,
        component=:matsubara,
        target_labels=((:orbital,), (:orbital,)),
        metadata=(source=:transform_test,),
    ))
    @test dlr_gf isa Gf
    @test dlr_gf.target_shape == (1, 1)
    @test only(dlr_gf.mesh) isa DLRFreq
    @test dlr_gf.component === :matsubara
    @test dlr_gf.statistics === FERMION
    @test dlr_gf.temperature isa FiniteTemperature
    @test inverse_temperature(dlr_gf.temperature) == 10.0
    @test dlr_gf.target_labels == ((:orbital,), (:orbital,))
    @test dlr_gf.metadata == (source=:transform_test,)

    multi_mesh_data = repeat(
        reshape(imfreq_data, 1, 1, 1, length(imfreq)), 1, 1, 2, 1,
    )
    multi_mesh_gf = Gf(1:2, imfreq;
        target_shape=(1, 1), data=multi_mesh_data, statistics=FERMION,
        component=:matsubara,
    )
    multi_mesh_dlr = to_dlr(multi_mesh_gf; dim=2)
    @test multi_mesh_dlr.mesh[1] == 1:2
    @test multi_mesh_dlr.mesh[2] isa DLRFreq
    @test_throws DimensionMismatch to_dlr(multi_mesh_gf; dim=3)

    up = Gf(mesh; target_shape=(2, 2), data=data,
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    down = Gf(mesh; target_shape=(1, 1),
        data=zeros(ComplexF64, 1, 1, length(mesh)),
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    blocks = BlockGf(:up => up, :down => down)
    @test collect(keys(blocks)) == [:up, :down]
    @test blocks[:up] === up
    @test first(blocks) == (:up => up)
    arbitrary_mesh = ReFreq(collect(mesh))
    arbitrary_up = Gf(arbitrary_mesh; target_shape=(2, 2), data=copy(data),
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    @test BlockGf(:uniform => up, :arbitrary => arbitrary_up) isa BlockGf
    retime_up = Gf(ReTime(collect(mesh)); target_shape=(2, 2), data=copy(data),
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    @test_throws ArgumentError BlockGf(:frequency => up, :time => retime_up)
    @test_throws KeyError blocks[:missing]
    @test_throws ArgumentError BlockGf(:up => up, :duplicate => Gf(
        ReFreq(-2.0, 2.0, 3); target_shape=(2, 2),
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature(),
    ))
    @test_throws ArgumentError BlockGf(:up => up, :down => Gf(
        mesh; target_shape=(2, 2), statistics=FERMION, component=:advanced,
        temperature=ZeroTemperature(),
    ))
    finite_up = Gf(mesh; target_shape=(2, 2), data=copy(data),
        statistics=FERMION, component=:retarded,
        temperature=FiniteTemperature(4.0))
    @test_throws ArgumentError BlockGf(:zero => up, :finite => finite_up)
end

@testset "Temperature regimes and mesh inference" begin
    @test_throws ArgumentError FiniteTemperature(0.0)
    @test_throws ArgumentError FiniteTemperature(-1.0)
    @test_throws ArgumentError FiniteTemperature(Inf)
    @test_throws ArgumentError FiniteTemperature(NaN)
    @test inverse_temperature(ZeroTemperature()) === nothing
    @test inverse_temperature(FiniteTemperature(3.0)) == 3.0

    iw = ImFreq(7.0, FERMION)
    inferred = Gf(iw; statistics=FERMION, component=:matsubara)
    @test inferred.temperature isa FiniteTemperature
    @test inferred.temperature.β == 7.0
    explicit = Gf(iw; statistics=FERMION, component=:matsubara,
        temperature=FiniteTemperature(7.0))
    @test explicit.temperature.β == 7.0
    @test_throws ArgumentError Gf(iw; statistics=FERMION,
        component=:matsubara, temperature=FiniteTemperature(8.0))
    @test_throws ArgumentError Gf(iw; statistics=FERMION,
        component=:matsubara, temperature=ZeroTemperature())
    @test_throws ArgumentError Gf(ImTime(7.0, FERMION), ImFreq(8.0, FERMION);
        statistics=FERMION, component=:matsubara)

    rw = ReFreq(-1.0, 1.0, 3)
    @test_throws ArgumentError Gf(rw; statistics=FERMION,
        temperature=ZeroTemperature())
    @test_throws ArgumentError Gf(rw; statistics=FERMION,
        component=:retarded)
    zero_real = Gf(rw; statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    finite_real = Gf(rw; statistics=FERMION, component=:retarded,
        temperature=FiniteTemperature(5.0))
    @test zero_real.temperature isa ZeroTemperature
    @test finite_real.temperature.β == 5.0
end

@testset "SpectralDensity validation" begin
    mesh = ReFreq(-1.0, 1.0, 3)
    data = zeros(ComplexF64, 2, 2, length(mesh))
    for i in eachindex(mesh)
        data[:, :, i] .= [2.0 1.0im; -1.0im 2.0]
    end
    gf = Gf(mesh; target_shape=(2, 2), data=data,
        statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature())
    spectral = SpectralDensity(gf; constraint=:psd)
    @test spectral.data[1, 2, 1] == 1im
    @test parent(spectral) === gf

    nonhermitian = copy(data)
    nonhermitian[2, 1, 2] = 1im
    @test_throws ArgumentError SpectralDensity(Gf(mesh;
        target_shape=(2, 2), data=nonhermitian,
        statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature(),
    ))

    indefinite = copy(data)
    indefinite[:, :, 1] .= [1.0 2.0; 2.0 1.0]
    indefinite_gf = Gf(mesh; target_shape=(2, 2), data=indefinite,
        statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature())
    @test SpectralDensity(indefinite_gf; constraint=:hermitian) isa SpectralDensity
    @test_throws ArgumentError SpectralDensity(indefinite_gf; constraint=:psd)

    @test_throws ArgumentError SpectralDensity(Gf(ReTime(0.0, 1.0, 3);
        target_shape=(2, 2), data=data,
        statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature(),
    ))
    @test_throws DimensionMismatch SpectralDensity(Gf(mesh;
        target_shape=(2, 3), statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature(),
    ))

    retarded = Gf(mesh; target_shape=(2, 2), data=data,
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
    @test_throws ArgumentError SpectralDensity(retarded)

    scalar_gf = Gf(mesh; data=ComplexF64[0.5, 1.0, 2.0],
        statistics=FERMION, component=:spectral,
        temperature=FiniteTemperature(2.0))
    scalar_spectral = SpectralDensity(scalar_gf; constraint=:psd)
    @test scalar_spectral.target_shape == ()
    @test scalar_spectral.data == scalar_gf.data
    @test_throws ArgumentError SpectralDensity(Gf(mesh;
        data=ComplexF64[0.5, 1.0im, 2.0], statistics=FERMION,
        component=:spectral, temperature=ZeroTemperature()))
    negative_scalar = Gf(mesh; data=ComplexF64[0.5, -0.1, 2.0],
        statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature())
    @test SpectralDensity(negative_scalar; constraint=:hermitian) isa SpectralDensity
    @test_throws ArgumentError SpectralDensity(negative_scalar; constraint=:psd)
end

@testset "Retarded, advanced, and spectral conversion" begin
    mesh = ReFreq(-2.0, 2.0, 5)
    epsilon, eta = 0.3, 0.2
    retarded_data = ComplexF64[inv(ω - epsilon + im * eta) for ω in mesh]
    retarded = Gf(mesh; data=retarded_data, statistics=FERMION,
        component=:retarded, temperature=ZeroTemperature(),
        metadata=(model=:single_pole,))
    advanced = advanced_from_retarded(retarded)
    spectral = spectral_from_retarded(retarded; constraint=:psd)
    @test advanced.component === :advanced
    @test advanced.temperature isa ZeroTemperature
    @test advanced.metadata == retarded.metadata
    @test advanced.data ≈ conj.(retarded.data)
    @test spectral.component === :spectral
    @test spectral.temperature isa ZeroTemperature
    @test spectral.data ≈ -imag.(retarded.data) ./ π

    hermitian_part = ComplexF64[0.2 0.3+0.4im; 0.3-0.4im -0.1]
    density = ComplexF64[1.0 0.2im; -0.2im 0.7]
    matrix_data = Array{ComplexF64}(undef, 2, 2, length(mesh))
    for i in eachindex(mesh)
        matrix_data[:, :, i] .= hermitian_part .- im * π .* density
    end
    matrix_retarded = Gf(mesh; target_shape=(2, 2), data=matrix_data,
        statistics=FERMION, component=:retarded,
        temperature=FiniteTemperature(4.0),
        target_labels=((:a, :b), (:a, :b)))
    matrix_advanced = advanced_from_retarded(matrix_retarded)
    matrix_spectral = spectral_from_retarded(matrix_retarded; constraint=:psd)
    for i in eachindex(mesh)
        @test matrix_advanced.data[:, :, i] ≈ adjoint(matrix_retarded.data[:, :, i])
        @test matrix_spectral.data[:, :, i] ≈ density
    end
    @test matrix_spectral.temperature.β == 4.0
    @test matrix_spectral.target_labels == matrix_retarded.target_labels
end

@testset "Equilibrium spectral relations" begin
    @test thermal_distribution(-1.0, FERMION, ZeroTemperature()) == 1.0
    @test thermal_distribution(0.0, FERMION, ZeroTemperature()) == 0.5
    @test thermal_distribution(1.0, FERMION, ZeroTemperature()) == 0.0
    @test_throws DomainError thermal_distribution(0.0, BOSON, ZeroTemperature())
    @test_throws DomainError thermal_distribution(0.0, BOSON, FiniteTemperature(2.0))
    @test thermal_distribution(-1.0, BOSON, ZeroTemperature()) == -1.0
    @test thermal_distribution(1.0, BOSON, ZeroTemperature()) == 0.0

    stable_temperature = FiniteTemperature(1.0)
    for statistics in (FERMION, BOSON)
        occupations = thermal_distribution.((-1000.0, 1000.0), Ref(statistics),
            Ref(stable_temperature))
        @test all(isfinite, occupations)
    end
    @test thermal_distribution(-1000.0, FERMION, stable_temperature) == 1.0
    @test thermal_distribution(1000.0, FERMION, stable_temperature) == 0.0
    @test thermal_distribution(-1000.0, BOSON, stable_temperature) == -1.0
    @test thermal_distribution(1000.0, BOSON, stable_temperature) == 0.0

    zero_mesh = ReFreq([-2.0, 0.0, 2.0])
    zero_spectral = SpectralDensity(Gf(zero_mesh; data=ones(ComplexF64, 3),
        statistics=FERMION, component=:spectral,
        temperature=ZeroTemperature()); constraint=:psd)
    zero_components = equilibrium_components(zero_spectral)
    @test zero_components.lesser.data ≈ ComplexF64[2π * im, π * im, 0]
    @test zero_components.greater.data ≈ ComplexF64[0, -π * im, -2π * im]
    @test zero_components.keldysh.data ≈
        zero_components.greater.data .+ zero_components.lesser.data

    beta = 1.7
    frequencies = [-2.0, -0.5, 0.5, 2.0]
    density = ComplexF64[0.4, 0.8, 1.2, 0.6]
    for statistics in (FERMION, BOSON)
        spectral = SpectralDensity(Gf(ReFreq(frequencies); data=density,
            statistics=statistics, component=:spectral,
            temperature=FiniteTemperature(beta)); constraint=:psd)
        components = equilibrium_components(spectral)
        xi = statistics === FERMION ? -1 : 1
        fdt = statistics === FERMION ?
            tanh.(beta .* frequencies ./ 2) :
            coth.(beta .* frequencies ./ 2)
        @test components.keldysh.data ≈
            components.greater.data .+ components.lesser.data
        @test components.greater.data ≈
            xi .* exp.(beta .* frequencies) .* components.lesser.data
        @test components.keldysh.data ≈
            fdt .* (components.greater.data .- components.lesser.data)
        @test im .* (components.greater.data .- components.lesser.data) ./ (2π) ≈
            spectral.data
        for component in values(components)
            @test component.temperature.β == beta
            @test component.statistics == statistics
        end
    end

    boson_zero = SpectralDensity(Gf(ReFreq([-2.0, 2.0]);
        data=ones(ComplexF64, 2), statistics=BOSON, component=:spectral,
        temperature=ZeroTemperature()); constraint=:psd)
    boson_components = equilibrium_components(boson_zero)
    @test boson_components.lesser.data ≈ ComplexF64[2π * im, 0]
    @test boson_components.greater.data ≈ ComplexF64[0, -2π * im]

    descending_frequencies = [2.0, 0.5, -0.5, -2.0]
    descending_density = ComplexF64[1.0, 2.0, 3.0, 4.0]
    descending = SpectralDensity(Gf(ReFreq(descending_frequencies);
        data=descending_density, statistics=FERMION, component=:spectral,
        temperature=FiniteTemperature(beta)); constraint=:psd)
    descending_lesser = lesser_from_spectral(descending)
    expected_occupations = thermal_distribution.(descending_frequencies,
        Ref(FERMION), Ref(FiniteTemperature(beta)))
    @test collect(only(descending.mesh)) == descending_frequencies
    @test descending_lesser.data ≈ 2π * im .* expected_occupations .* descending_density
end
