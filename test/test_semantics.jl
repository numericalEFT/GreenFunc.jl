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
    scalar = Gf(scalar_ma; statistics=FERMION)
    @test scalar.target_ndim == 0
    @test scalar.target_shape == ()
    @test scalar.mesh == (mesh,)
    @test parent(scalar) === scalar_ma

    data = zeros(ComplexF64, 2, 2, length(mesh))
    data[1, 2, :] .= 1 .+ 2im
    data[2, 1, :] .= 3 .- 4im
    matrix_gf = Gf(mesh;
        target_shape=(2, 2), data=data, statistics=FERMION,
        component=:retarded,
        target_labels=((:a, :b), (:a, :b)),
        metadata=(source=:test,),
    )
    @test matrix_gf.target_ndim == 2
    @test matrix_gf.target_shape == (2, 2)
    @test matrix_gf.fullmesh[1:2] == (Base.OneTo(2), Base.OneTo(2))
    @test matrix_gf.data[1, 2, 2] == 1 + 2im
    @test matrix_gf.data[2, 1, 2] == 3 - 4im
    @test matrix_gf.metadata.source == :test
    @test copy(matrix_gf).data == matrix_gf.data

    @test_throws DimensionMismatch Gf(mesh;
        target_shape=(2, 2), statistics=FERMION,
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
    ))
    @test dlr_gf isa Gf
    @test dlr_gf.target_shape == (1, 1)
    @test only(dlr_gf.mesh) isa DLRFreq

    multi_mesh_data = repeat(
        reshape(imfreq_data, 1, 1, 1, length(imfreq)), 1, 1, 2, 1,
    )
    multi_mesh_gf = Gf(1:2, imfreq;
        target_shape=(1, 1), data=multi_mesh_data, statistics=FERMION,
    )
    multi_mesh_dlr = to_dlr(multi_mesh_gf; dim=2)
    @test multi_mesh_dlr.mesh[1] == 1:2
    @test multi_mesh_dlr.mesh[2] isa DLRFreq
    @test_throws DimensionMismatch to_dlr(multi_mesh_gf; dim=3)

    up = Gf(mesh; target_shape=(2, 2), data=data,
        statistics=FERMION, component=:retarded)
    down = Gf(mesh; target_shape=(1, 1),
        data=zeros(ComplexF64, 1, 1, length(mesh)),
        statistics=FERMION, component=:retarded)
    blocks = BlockGf(:up => up, :down => down)
    @test collect(keys(blocks)) == [:up, :down]
    @test blocks[:up] === up
    @test first(blocks) == (:up => up)
    arbitrary_mesh = ReFreq(collect(mesh))
    arbitrary_up = Gf(arbitrary_mesh; target_shape=(2, 2), data=copy(data),
        statistics=FERMION, component=:retarded)
    @test BlockGf(:uniform => up, :arbitrary => arbitrary_up) isa BlockGf
    retime_up = Gf(ReTime(collect(mesh)); target_shape=(2, 2), data=copy(data),
        statistics=FERMION, component=:retarded)
    @test_throws ArgumentError BlockGf(:frequency => up, :time => retime_up)
    @test_throws KeyError blocks[:missing]
    @test_throws ArgumentError BlockGf(:up => up, :duplicate => Gf(
        ReFreq(-2.0, 2.0, 3); target_shape=(2, 2),
        statistics=FERMION, component=:retarded,
    ))
    @test_throws ArgumentError BlockGf(:up => up, :down => Gf(
        mesh; target_shape=(2, 2), statistics=FERMION, component=:advanced,
    ))
end

@testset "SpectralDensity validation" begin
    mesh = ReFreq(-1.0, 1.0, 3)
    data = zeros(ComplexF64, 2, 2, length(mesh))
    for i in eachindex(mesh)
        data[:, :, i] .= [2.0 1.0im; -1.0im 2.0]
    end
    gf = Gf(mesh; target_shape=(2, 2), data=data,
        statistics=FERMION, component=:spectral)
    spectral = SpectralDensity(gf; constraint=:psd)
    @test spectral.data[1, 2, 1] == 1im
    @test parent(spectral) === gf

    nonhermitian = copy(data)
    nonhermitian[2, 1, 2] = 1im
    @test_throws ArgumentError SpectralDensity(Gf(mesh;
        target_shape=(2, 2), data=nonhermitian,
        statistics=FERMION, component=:spectral,
    ))

    indefinite = copy(data)
    indefinite[:, :, 1] .= [1.0 2.0; 2.0 1.0]
    indefinite_gf = Gf(mesh; target_shape=(2, 2), data=indefinite,
        statistics=FERMION, component=:spectral)
    @test SpectralDensity(indefinite_gf; constraint=:hermitian) isa SpectralDensity
    @test_throws ArgumentError SpectralDensity(indefinite_gf; constraint=:psd)

    @test_throws ArgumentError SpectralDensity(Gf(ReTime(0.0, 1.0, 3);
        target_shape=(2, 2), data=data,
        statistics=FERMION, component=:spectral,
    ))
    @test_throws DimensionMismatch SpectralDensity(Gf(mesh;
        target_shape=(2, 3), statistics=FERMION, component=:spectral,
    ))
end
