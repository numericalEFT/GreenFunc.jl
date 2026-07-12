# Green's Functions

`Gf` adds target-space, statistics, component, and temperature semantics without
changing the underlying `MeshArray` layout. Target axes are always first and are
only created when `target_shape` is explicit. Real-axis functions must state both
their component and temperature regime explicitly:

```julia
mesh = ReFreq(-4.0, 4.0, 401)
data = zeros(ComplexF64, 2, 2, length(mesh))
data[1, 2, :] .= 0.1im
data[2, 1, :] .= -0.1im
g = Gf(mesh; target_shape=(2, 2), data=data,
       statistics=FERMION, component=:retarded,
       temperature=ZeroTemperature(),
       target_labels=((:up, :down), (:up, :down)))
```

Use `ZeroTemperature()` for an exact ground-state context and
`FiniteTemperature(β)` for finite inverse temperature `β > 0`. `β = Inf` is
deliberately not used to represent zero temperature. `ImTime`, `ImFreq`, and
`DLRFreq` meshes infer `FiniteTemperature(mesh.β)`; an explicitly supplied
finite temperature must agree with every imaginary mesh. Real-frequency and
real-time meshes cannot infer a temperature. They also require an explicit
component such as `:retarded`, `:advanced`, `:lesser`, `:greater`, `:keldysh`, or
`:spectral`.

For `Gf` transforms, an explicit `dim` indexes `g.mesh` (the physical meshes),
not the leading target axes stored in the parent `MeshArray`. `copy`, `similar`,
and mesh transforms preserve statistics, component, temperature, target labels,
and metadata.

Named blocks preserve insertion order and require compatible physical meshes,
statistics, components, and identical temperature regimes:

```julia
blocks = BlockGf(:up => g, :down => copy(g))
```

`SpectralDensity` accepts scalar targets or square matrix targets on `ReFreq`.
Its wrapped `Gf` must have `component=:spectral`. Scalar values must be real;
matrix values must be Hermitian at every frequency. Use
`constraint=:hermitian` (the default), `:psd`, or `:none`.

## Retarded and spectral components

`advanced_from_retarded` constructs the advanced component by complex
conjugation for a scalar target and by the adjoint at each physical point for a
matrix target. `spectral_from_retarded` uses
`ρ(ω) = i [Gᴿ(ω) - Gᴬ(ω)] / (2π)` and validates the result:

```julia
ω = ReFreq(-4.0, 4.0, 401)
gR = Gf(ω; data=ComplexF64[inv(w - 0.3 + 0.1im) for w in ω],
        statistics=FERMION, component=:retarded,
        temperature=ZeroTemperature())
gA = advanced_from_retarded(gR)
ρ = spectral_from_retarded(gR; constraint=:psd)
```

## Equilibrium components

`thermal_distribution(energy, statistics, temperature)` evaluates stable Fermi
or Bose occupations. At zero temperature the fermionic value at the step is
`1/2`; the Bose distribution at zero energy is singular and throws a
`DomainError`.

Given a `SpectralDensity`, `lesser_from_spectral`, `greater_from_spectral`, and
`keldysh_from_spectral` construct equilibrium real-frequency components.
`equilibrium_components(ρ)` returns all three as a named tuple. With
`ξ = -1` for fermions and `ξ = +1` for bosons, the conventions are

```math
G^<(\omega) = -2\pi i\,\xi n_\xi(\omega)\rho(\omega),\qquad
G^>(\omega) = -2\pi i\,[1+\xi n_\xi(\omega)]\rho(\omega),
```

```math
G^K = G^> + G^<,\qquad
\rho = \frac{i}{2\pi}(G^>-G^<),\qquad
G^> = \xi e^{\beta\omega}G^<.
```

For a descending `ReFreq`, data and generated thermal factors retain that same
user-visible frequency order.

```@autodocs
Modules = [GreenFunc]
```
