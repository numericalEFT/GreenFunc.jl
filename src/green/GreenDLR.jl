"""
General Green's function. 
"""

"""
Green's function on a multi-dimensional mesh plus one in-built Discrete Lehmann Representation.
#Parameters:
- 'T': type of data
- 'TType': type of time domain, TType<:TimeDomain
- 'TGT': type of time grid
- 'SGT': type of space grid

#Members:
"""
mutable struct GreenDLR{T<:Any,Type<:TimeDomain,TGT,SGT}<:AbstractArray
    name::Symbol
    dlrGrid::DLRGrid

    #########   Mesh   ##############
    
    timeGrid::TGT
    mesh::SGT
    inner_state::Array{::Int, 1}
    ###########     data   ###########
    data::Array{T}


    function GreenDLR{T}(name::Symbol, timeType::TT, β, isFermi::Bool, Euv, spaceGrid; inner_state=[1,],
        timeSymmetry::Symbol = :none, rtol = 1e-8, kwargs...
    ) where {T<:Any,TT<:TimeDomain}
       
        gnew = new{T,TT,typeof(timeGrid),typeof(mesh)}(
            name,  dlrGrid,
            timeGrid,
             mesh, inner_state,
            data)
        return gnew
    end
end

"""
def __getitem__(self, key):

        # First case : g[:] = RHS ... will be g << RHS
        if key == self._full_slice:
            return self

        # Only one argument. Must be a mesh point, idx or slicing rank1 target space
        if not isinstance(key, tuple):
            if isinstance(key, (MeshPoint, Idx)):
                return self.data[key.linear_index if isinstance(key, MeshPoint) else self._mesh.index_to_linear(key.idx)]
            else: key = (key,)

        # If all arguments are MeshPoint, we are slicing the mesh or evaluating
        if all(isinstance(x, (MeshPoint, Idx)) for x in key):
            assert len(key) == self.rank, "wrong number of arguments in [ ]. Expected %s, got %s"%(self.rank, len(key))
            return self.data[tuple(x.linear_index if isinstance(x, MeshPoint) else m.index_to_linear(x.idx) for x,m in zip(key,self._mesh._mlist))]

        # If any argument is a MeshPoint, we are slicing the mesh or evaluating
        elif any(isinstance(x, (MeshPoint, Idx)) for x in key):
            assert len(key) == self.rank, "wrong number of arguments in [[ ]]. Expected %s, got %s"%(self.rank, len(key))
            assert all(isinstance(x, (MeshPoint, Idx, slice)) for x in key), "Invalid accessor of Greens function, please combine only MeshPoints, Idx and slice"
            assert self.rank > 1, "Internal error : impossible case" # here all == any for one argument
            mlist = self._mesh._mlist 
            for x in key:
                if isinstance(x, slice) and x != self._full_slice: raise NotImplementedError("Partial slice of the mesh not implemented")
            # slice the data
            k = tuple(x.linear_index if isinstance(x, MeshPoint) else m.index_to_linear(x.idx) if isinstance(x, Idx) else x for x,m in zip(key,mlist)) + self._target_rank * (slice(0, None),)
            dat = self._data[k]
            # list of the remaining lists
            mlist = [m for i,m in filter(lambda tup_im : not isinstance(tup_im[0], (MeshPoint, Idx)), zip(key, mlist))]
            assert len(mlist) > 0, "Internal error" 
            mesh = MeshProduct(*mlist) if len(mlist)>1 else mlist[0]
            sing = None 
            r = Gf(mesh = mesh, data = dat)
            r.__check_invariants()
            return r

        # In all other cases, we are slicing the target space
        else : 
            assert self.target_rank == len(key), "wrong number of arguments. Expected %s, got %s"%(self.target_rank, len(key))

            # Assume empty indices (scalar_valued)
            ind = GfIndices([])

            # String access: transform the key into a list integers
            if all(isinstance(x, str) for x in key):
                warnings.warn("The use of string indices is deprecated", DeprecationWarning)
                assert self._indices, "Got string indices, but I have no indices to convert them !"
                key_tpl = tuple(self._indices.convert_index(s,i) for i,s in enumerate(key)) # convert returns a slice of len 1

            # Slicing with ranges -> Adjust indices
            elif all(isinstance(x, slice) for x in key): 
                key_tpl = tuple(key)
                ind = GfIndices([ v[k]  for k,v in zip(key_tpl, self._indices.data)])

            # Integer access
            elif all(isinstance(x, int) for x in key):
                key_tpl = tuple(key)

            # Invalid Access
            else:
                raise NotImplementedError("Partial slice of the target space not implemented")

            dat = self._data[ self._rank*(slice(0,None),) + key_tpl ]
            r = Gf(mesh = self._mesh, data = dat, indices = ind)

            r.__check_invariants()
            return r
"""
function Base.getindex(green::GreenDLR,key...)


end

"""
def __setitem__(self, key, val):

    # Only one argument and not a slice. Must be a mesh point, Idx
    if isinstance(key, (MeshPoint, Idx)):
    self.data[key.linear_index if isinstance(key, MeshPoint) else self._mesh.index_to_linear(key.idx)] = val

        # If all arguments are MeshPoint, we are slicing the mesh or evaluating
        elif isinstance(key, tuple) and all(isinstance(x, (MeshPoint, Idx)) for x in key):
            assert len(key) == self.rank, "wrong number of arguments in [ ]. Expected %s, got %s"%(self.rank, len(key))
        self.data[tuple(x.linear_index if isinstance(x, MeshPoint) else m.index_to_linear(x.idx) for x,m in zip(key,self._mesh._mlist))] = val

    else:
        self[key] << val
"""
function Base.setindex!(green::GreenDLR; kwargs...)

end

function Base.view(green::GreenDLR,inds...)
    

end

function Base.getproperty(obj::GreenDLR{T,TT,TGT,SGT}, sym::Symbol) where {T,TT,TGT,SGT}
    if sym === :isFermi
        return obj.dlrGrid.isFermi
    elseif sym === :β
        return obj.dlrGrid.β
    elseif sym === :timeType
        return TT
    elseif sym === :timeSymmetry
        return obj.dlrGrid.symmetry
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

function Base.size(green::GreenDLR)
    return (green.color, green.color, size(green.spaceGrid), size(green.timeGrid))
end

"""
   def __lshift__(self, A):
        A can be two things:
          * G << any_init will init the GFBloc with the initializer
          * G << g2 where g2 is a GFBloc will copy g2 into self
        if isinstance(A, Gf):
            if self is not A: # otherwise it is useless AND does not work !!
                assert self.mesh == A.mesh, "Green function meshes are not compatible:\n  %s\nand\n  %s" % (self.mesh, A.mesh)
                self.copy_from(A)
        elif isinstance(A, lazy_expressions.LazyExpr): # A is a lazy_expression made of GF, scalars, descriptors
            A2 = descriptors.convert_scalar_to_const(A)
            def e_t (x):
                if not isinstance(x, descriptors.Base): return x
                tmp = self.copy()
                x(tmp)
                return tmp
            self.copy_from (lazy_expressions.eval_expr_with_context(e_t, A2) )
        elif isinstance(A, lazy_expressions.LazyExprTerminal): #e.g. g<< SemiCircular (...)
            self << lazy_expressions.LazyExpr(A)
        elif descriptors.is_scalar(A): #in the case it is a scalar ....
            self << lazy_expressions.LazyExpr(A)
        else:
            raise NotImplemented
                    return self
"""
function Base.:<<(green::GreenDLR, greenSrc::GreenDLR)
end


"""
 def copy(self) : 
        Deep copy of the Greens function.
            Returns
            -------
                G : Gf
            Copy of self.
        return Gf (mesh = self._mesh.copy(), 
                   data = self._data.copy(), 
                   indices = self._indices.copy(), 
                   name = self.name)

"""
function Base.deepcopy(green::GreenDLR)

end


function Base.:-(green::GreenDLR)

end

function Base.:+(green::GreenDLR)

end

function Base.:+(greenL::GreenDLR, greenR::GreenDLR)

end

function Base.:-(greenL::GreenDLR, greenR::GreenDLR)

end

function Base.:*(greenL::GreenDLR, greenR::GreenDLR)

end

"""
 def copy_from(self, another):
        Copy the data of another Greens function into self.
        self._mesh.copy_from(another.mesh)
        assert self._data.shape == another._data.shape, "Shapes are incompatible: " + str(self._data.shape) + " vs " + str(another._data.shape)
        self._data[:] = another._data[:]
        self._indices = another._indices.copy()
        self.__check_invariants()
"""

function copy_from(green::GreenDLR, greenSrc::GreenDLR)

end


"""
def density(self, *args, **kwargs):
    rCompute the density matrix of the Greens function
        Parameters
        ----------
        beta : float, optional
            Used for finite temperature density calculation with ``MeshReFreq``.
        Returns
        -------
        density_matrix : ndarray
            Single particle density matrix with shape ``target_shape``.
        Notes
        -----
        Only works for single mesh Greens functions with a, Matsubara,
        real-frequency, or Legendre mesh.
return gf_fnt.density(self, *args, **kwargs)
"""
function density(green::GreenDLR; kwargs...)

end
"""
def total_density(self, *args, **kwargs):
        Compute total density.
            Returns
-------
    density : float
Total density of the Greens function.
    Notes
    -----
        Only implemented for single mesh Greens function with a,
            Matsubara, real-frequency, or Legendre mesh.
               
        return np.trace(gf_fnt.density(self, *args, **kwargs))
"""
function total_density(green::GreenDLR; kwargs...)

end

"""
   def transpose(self): 
        Take the transpose of a matrix valued Greens function.
            Returns
            -------
                G : Gf (copy)
            The transpose of the Greens function.
                Notes
                -----
                    Only implemented for single mesh matrix valued Greens functions.
                   

        # FIXME Why this assert ?
        #assert any( (isinstance(self.mesh, x) for x in [meshes.MeshImFreq, meshes.MeshReFreq])), "Method invalid for this Gf"

        assert self.rank == 1, "Transpose only implemented for single mesh Greens functions"
        assert self.target_rank == 2, "Transpose only implemented for matrix valued Greens functions"

        d = np.transpose(self.data.copy(), (0, 2, 1))
        return Gf(mesh = self.mesh, data= d, indices = self.indices.transpose())
"""

function transpose(green:GreenDLR)


end
"""
    def from_L_G_R(self, L, G, R):
        rMatrix transform of the target space of a matrix valued Greens function.
        Sets the current Greens function :math:`g_{ab}` to the matrix transform of :math:`G_{cd}`
        using the left and right transform matrices :math:`L_{ac}` and :math:`R_{db}`.
        .. math::
            g_{ab} = \sum_{cd} L_{ac} G_{cd} R_{db}
        Parameters
        ----------
        L : (a, c) ndarray
            Left side transform matrix.
        G : Gf matrix valued target_shape == (c, d)
            Greens function to transform.
        R : (d, b) ndarray
            Right side transform matrix.
        Notes
        -----
        Only implemented for Greens functions with a single mesh.
       

        assert self.rank == 1, "Only implemented for Greens functions with one mesh"
        assert self.target_rank == 2, "Matrix transform only valid for matrix valued Greens functions"

        assert len(L.shape) == 2, "L needs to be two dimensional"
        assert len(R.shape) == 2, "R needs to be two dimensional"

        assert L.shape[1] == G.target_shape[0], "Dimension mismatch between L and G"
        assert R.shape[0] == G.target_shape[1], "Dimension mismatch between G and R"

        assert L.shape[0] == self.target_shape[0], "Dimension mismatch between L and self"
        assert R.shape[1] == self.target_shape[1], "Dimension mismatch between R and self"

        if not L.strides == sorted(L.strides):
            L = L.copy(order='C')

        if not R.strides == sorted(R.strides):
            R = R.copy(order='C')

        wrapped_aux.set_from_gf_data_mul_LR(self.data, L, G.data, R)
"""
function from_L_G_R(self,L,G::GreenDLR,R)

end





"""
    function toTau(green::Green2DLR, targetGrid =  green.dlrGrid.τ)
Convert Green's function to τ space by Fourier transform.
If green is already in τ space then it will be interpolated to the new grid.

#Arguements
- 'green': Original Green's function
- 'targetGrid': Grid of outcome Green's function. Default: DLR τ grid
"""
function toTau(green::Green2DLR, targetGrid = green.dlrGrid.τ)

    if targetGrid isa AbstractVector
        targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
    end

    # do nothing if the domain and the grid remain the same
    if green.timeType == ImTime && length(green.timeGrid.grid) ≈ length(targetGrid.grid) && green.timeGrid.grid ≈ targetGrid.grid
        return green
    end
    if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
        return green
    end


    if (green.timeType == ImTime)
        dynamic = tau2tau(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == ImFreq)
        dynamic = matfreq2tau(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == DLRFreq)
        dynamic = dlr2tau(green.dlrGrid, green.dynamic, targetGrid.grid; axis = 4)
    end

    return Green2DLR{eltype(dynamic)}(
        green.name, IMTIME, green.β, green.isFermi, green.dlrGrid.Euv, green.spaceGrid, green.color;
        timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, rtol = green.dlrGrid.rtol,
        dynamic = dynamic, instant = green.instant)
end

"""
    function toMatFreq(green::Green2DLR, targetGrid =  green.dlrGrid.n)
Convert Green's function to matfreq space by Fourier transform.
If green is already in matfreq space then it will be interpolated to the new grid.

#Arguements
- 'green': Original Green's function
- 'targetGrid': Grid of outcome Green's function. Default: DLR n grid
"""
function toMatFreq(green::Green2DLR, targetGrid = green.dlrGrid.n)

    if targetGrid isa AbstractVector
        targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
    end

    # do nothing if the domain and the grid remain the same
    if green.timeType == ImFreq && length(green.timeGrid.grid) ≈ length(targetGrid.grid) && green.timeGrid.grid ≈ targetGrid.grid
        return green
    end
    if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
        return green
    end


    if (green.timeType == ImFreq)
        dynamic = matfreq2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == ImTime)
        dynamic = tau2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == DLRFreq)
        dynamic = dlr2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid; axis = 4)
    end

    return Green2DLR{eltype(dynamic)}(
        green.name, IMFREQ, green.β, green.isFermi, green.dlrGrid.Euv, green.spaceGrid, green.color;
        timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, rtol = green.dlrGrid.rtol,
        dynamic = dynamic, instant = green.instant)

end

"""
    function toDLR(green::Green2DLR)
Convert Green's function to dlr space.

#Arguements
- 'green': Original Green's function
"""
function toDLR(green::Green2DLR)

    # do nothing if the domain and the grid remain the same
    if green.timeType == DLRFreq
        return green
    end
    if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
        return green
    end


    if (green.timeType == ImTime)
        dynamic = tau2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == ImFreq)
        dynamic = matfreq2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid; axis = 4)
    end

    return Green2DLR{eltype(dynamic)}(
        green.name, DLRFREQ, green.β, green.isFermi, green.dlrGrid.Euv, green.spaceGrid, green.color;
        timeSymmetry = green.timeSymmetry, timeGrid = green.dlrGrid.ω, rtol = green.dlrGrid.rtol,
        dynamic = dynamic, instant = green.instant)

end



"""
    function dynamic(green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, color1::Int, color2::Int, timeMethod::TM , spaceMethod::SM) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid,TM,SM}
Find value of Green's function's dynamic part at given color and k/x by interpolation.
Interpolation method is by default depending on the grid, but could also be chosen to be linear.
#Argument
- 'green': Green's function
- 'time': Target τ/ω_n point
- 'space': Target k/x point
- 'color1': Target color1
- 'color2': Target color2
- 'timeMethod': Method of interpolation for time
- 'spaceMethod': Method of interpolation for space 
"""
function _interpolation(TIM, SIM; green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, color1::Int, color2::Int) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid}
    # for double composite
    if isempty(green.data)
        error("Dynamic Green's function can not be empty!")
    else
        spaceNeighbor = CompositeGrids.Interp.findneighbor(SIM, green.spaceGrid, space)
        println(TIM)
        if green.timeType == ImFreq && TIM != DLRInterp
            timeGrid = (green.timeGrid.grid * 2 .+ 1) * π / green.β
            comTimeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(timeGrid)}(timeGrid)
            comTime = (2*time+1)*π/green.β
        else
            comTimeGrid = green.timeGrid
            comTime = time
        end

        timeNeighbor = CompositeGrids.Interp.findneighbor(TIM, comTimeGrid, comTime)
        data_slice = view(green.data, color1, color2, spaceNeighbor.index, timeNeighbor.index)
        data_slice_xint = CompositeGrids.Interp.interpsliced(spaceNeighbor,data_slice, axis=1)
        result = CompositeGrids.Interp.interpsliced(timeNeighbor,data_slice_xint, axis=1)
    end
    return result
end

function _interpolation( ::LinearInterp , ::LinearInterp
    ;green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, color1::Int, color2::Int,
    ) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid}
    # for double composite and double linear
    if isempty(green.data)
        error("Dynamic Green's function can not be empty!")
    else
        if green.timeType == ImFreq
            timeGrid = (green.timeGrid.grid * 2 .+ 1) * π / green.β
            comTimeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(timeGrid)}(timeGrid)            
            comTime = (2*time+1)*π/green.β
        else
            comTimeGrid = green.timeGrid
            comTime = time
        end
        data_slice = view(green.data, color1, color2, :,:)
        result = CompositeGrids.Interp.linear2D(data_slice, green.spaceGrid, comTimeGrid,space,comTime)
    end
    return result
end


function _interpolation(::DLRInterp, SIM; green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, color1::Int, color2::Int) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid}
    # for composite space and dlr time
    if isempty(green.data)
        error("Dynamic Green's function can not be empty!")
    else
        spaceNeighbor = CompositeGrids.Interp.findneighbor(SIM, green.spaceGrid, space)
        data_slice = view(green.data, color1, color2, spaceNeighbor.index,:)
        data_slice_xint = CompositeGrids.Interp.interpsliced(spaceNeighbor,data_slice, axis=1)
        if green.timeType == ImFreq
            result = (matfreq2matfreq(green.dlrGrid, data_slice_xint, [time,], green.timeGrid.grid))[1]
        elseif green.timeType == ImTime
            result = (tau2tau(green.dlrGrid, data_slice_xint, [time,], green.timeGrid.grid))[1]
        end
    end
    return result
end

interpolation(;timeMethod::TM, spaceMethod::SM ,green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space , color1::Int, color2::Int) where {TM,SM,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid,DT,TT} = _interpolation(InterpMethod(TGT,TM), InterpMethod(SGT, SM); green, time, space, color1, color2)

function interpolation(green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, color1::Int=1, color2::Int=color1, timeMethod::TM=DEFAULTINTERP , spaceMethod::SM=DEFAULTINTERP) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid,TM,SM}
    return  interpolation(; timeMethod = timeMethod, spaceMethod = spaceMethod, green = green, time=time, space=space, color1=color1, color2 =color2)
end

function interpolation(green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, timeMethod::TM, spaceMethod::SM) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid,TM,SM}
    return  interpolation(; timeMethod = timeMethod, spaceMethod = spaceMethod, green = green, time=time, space=space, color1=1, color2 =1)
end
