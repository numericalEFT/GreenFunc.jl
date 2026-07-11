# Green's Functions

`Gf` adds target-space and statistics semantics without changing the underlying
`MeshArray` layout. Target axes are always first and are only created when
`target_shape` is explicit:

```julia
mesh = ReFreq(-4.0, 4.0, 401)
data = zeros(ComplexF64, 2, 2, length(mesh))
data[1, 2, :] .= 0.1im
data[2, 1, :] .= -0.1im
g = Gf(mesh; target_shape=(2, 2), data=data,
       statistics=FERMION, component=:retarded,
       target_labels=((:up, :down), (:up, :down)))
```

For `Gf` transforms, an explicit `dim` indexes `g.mesh` (the physical meshes),
not the leading target axes stored in the parent `MeshArray`.

Named blocks preserve insertion order and require compatible physical meshes,
statistics, and components:

```julia
blocks = BlockGf(:up => g, :down => copy(g))
```

`SpectralDensity` validates full complex target matrices at every physical-mesh
point. Use `constraint=:hermitian` (the default), `:psd`, or `:none`.

```@autodocs
Modules = [GreenFunc]
```
