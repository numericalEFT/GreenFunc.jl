var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GreenFunc","category":"page"},{"location":"#GreenFunc","page":"Home","title":"GreenFunc","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GreenFunc.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GreenFunc]","category":"page"},{"location":"#GreenFunc.dlr_to_imfreq-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Any}} where {T, N, MT}","page":"Home","title":"GreenFunc.dlr_to_imfreq","text":"function dlr_to_imfreq(obj::MeshArray, ngrid=nothing; dim=nothing)\n\nTransform a Green's function in DLR space to Matsubara frequency space.  #Arguements\n\n'obj': Function in DLR space\nngrid: The Matsubara-frequency grid which the function transforms into. Default value is the Matsubara-frequency grid from the ImFreq constructor.\ndim: The dimension of the temporal mesh. Default value is the first ImFreq mesh.\n\n\n\n\n\n","category":"method"},{"location":"#GreenFunc.dlr_to_imtime-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Any}} where {T, N, MT}","page":"Home","title":"GreenFunc.dlr_to_imtime","text":"function dlr_to_imtime(obj::MeshArray, tgrid=nothing; dim=nothing)\n\nTransform a Green's function in DLR space to the imaginary-time space.  #Arguements\n\n'obj': Function in DLR space\ntgrid: The imaginary-time grid which the function transforms into. Default value is the imaginary-time grid from the ImFreq constructor.\ndim: The dimension of the temporal mesh. Default value is the first ImTime mesh.\n\n\n\n\n\n","category":"method"},{"location":"#GreenFunc.imfreq_to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"Home","title":"GreenFunc.imfreq_to_dlr","text":"function imfreq_to_dlr(obj::MeshArray; dim=nothing, rtol=1e-12, sym=:none)\n\nCalculate the DLR sepctral density of a Matsubara-frequency Green's function. #Arguements\n\n'obj': Function in the Matsubara-frequency space\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImFreq.\nrtol: The relative tolerance of the DLR transform. Default value is 1e-12.\nsym: The symmetry of the Green's function, :none, :ph or :pha. Default value is :none.\n\n\n\n\n\n","category":"method"},{"location":"#GreenFunc.imtime_to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"Home","title":"GreenFunc.imtime_to_dlr","text":"function imtime_to_dlr(obj::MeshArray; dim=nothing, rtol=1e-12, sym=:none)\n\nCalculate the DLR sepctral density of an imaginary-time Green's function. #Arguements\n\n'obj': Function in the imaginary-time space\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime.\nrtol: The relative tolerance of the DLR transform. Default value is 1e-12.\nsym: The symmetry of the Green's function, :none, :ph or :pha. Default value is :none.\n\n\n\n\n\n","category":"method"},{"location":"#GreenFunc.to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"Home","title":"GreenFunc.to_dlr","text":"function to_dlr(obj::MeshArray; dim=nothing, rtol=1e-12, sym=:none)\n\nCalculate the DLR sepctral density of an imaginary-time or Matsubara-frequency Green's function. #Arguements\n\n'obj': Function in the imaginary-time space or in the Matsubara-frequency\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime or ImFreq.\nrtol: The relative tolerance of the DLR transform. Default value is 1e-12.\nsym: The symmetry of the Green's function, :none, :ph or :pha. Default value is :none.\n\n\n\n\n\n","category":"method"}]
}
