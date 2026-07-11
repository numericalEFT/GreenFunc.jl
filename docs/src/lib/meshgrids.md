# MeshGrids

In addition to `ImTime`, `ImFreq`, and `DLRFreq`, real-axis meshes are
available as `ReTime` and `ReFreq`. They accept either explicit sorted points
or a uniform window:

```julia
ω = ReFreq(-5.0, 5.0, 1001)
t = ReTime(; window=(0.0, 20.0), n_points=401)
```

```@autodocs
Modules = [GreenFunc.MeshGrids]
```
