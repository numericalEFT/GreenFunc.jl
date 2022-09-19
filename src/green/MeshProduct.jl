"""
Cartisian product of 1 dimensional meshes 
"""

"""

#Members:
"""
struct MeshProduct{MT}
    mlist::MT
    function MeshProduct(vargs...
        )
        #@assert all(v -> (v isa Mesh), vargs) "all arguments should variables"
        mprod = Tuple(v for v in vargs)
        mnew = new{typeof(mprod)}(mprod)
        return mnew
    end
end


#TODO:return the sum of length of all meshes
function Base.length(obj::MeshProduct)
    len = 1
    for l in length.(obj.mlist)
        len = len * l
    end
    return len
end

#TODO:return the size of Ith mesh in mlist 
function Base.size(obj::MeshProduct, I::Int)
    return length.(obj.mlist)[I]
end

#TODO:return the size of all meshes in a tuple 
function Base.size(obj::MeshProduct)
    return length.(obj.mlist)
end

function rank(obj::MeshProduct)
    return length(obj.mlist)
end

#TODO: find the linearindex I corresponding to the given index
function index_to_linear(obj::MeshProduct, index...)
    bn = reverse(Tuple(prod(size(obj)[1:n-1]) for (n,sz) in enumerate(size(obj))))
    li = 1
    for (i,id) in enumerate(index)
      li = li+(id-1)*bn[i]
    end
    return li
end

#TODO: find the index corresponding to linearindex I
function linear_to_index(obj::MeshProduct, I::Int)
    d =rank(obj)
    bn = reverse(Tuple(prod(size(obj)[1:n-1]) for (n,sz) in enumerate(size(obj))))
    index = zeros(Int32,d)
    index[1] = (I-1) ÷ bn[1] +1
    for k in 2:d
        index[k]=((I-1)%bn[k-1])÷bn[k]+1 
    end
    return Tuple(index)
end



#TODO:for all n meshes in mlist, return [..., (mlist[i])[index[i]], ...] 
function Base.getindex(obj::MeshProduct, index...)
    return Tuple(obj.mlist[i][id] for (i,id) in enumerate(index))
end
#TODO:find the index corresponds to linearindex I, then for all n meshes in mlist, return [..., (mlist[i])[index[i]], ...] 
function Base.getindex(obj::MeshProduct, I::Int)
    return getindex(meshprod,linear_to_index(meshprod,I)[1])
end

#TODO:return the sliced pieces of 
function Base.view(obj::MeshProduct,inds...)
    return 1
end


Base.firstindex(obj::MeshProduct) = 1
Base.lastindex(obj::MeshProduct) = length(obj)
# iterator
Base.iterate(obj::MeshProduct) = (obj[1],1)
Base.iterate(obj::MeshProduct, state) = (state>=length(obj)) ? nothing : (obj[state+1],state+1)


#TODO:nice print
function Base.show(obj::MeshProduct)
    return 1
end


"""
 All meshes in mlist should have locate and volume functions. Here in meshproduct we just delegate these functions to the meshes, and return the proper array of returned values.
"""
function locate(obj::MeshProduct, index...)
    return Tuple(locate(obj, index[mi]) for (mi, m) in enumerate(obj))
end

function volume(obj::MeshProduct, index...)
    return reduce(*, volume(obj, index[mi]) for (mi, m) in enumerate(obj))
end

function locate(obj::MeshProduct, I::Int)
    index = linear_to_index(obj,I)
    return (locate(obj, index[mi]) for (mi, m) in enumerate(obj))
end

function volume(obj::MeshProduct, I::Int)
    index = linear_to_index(obj,I)
    return reduce(*, volume(obj, index[mi]) for (mi, m) in enumerate(obj))
end

