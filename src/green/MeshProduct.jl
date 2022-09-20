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
    return 1
end

#TODO:return the size of Ith mesh in mlist 
function Base.size(obj::MeshProduct, I::Int)
    return 1
end

#TODO:return the size of all meshes in a tuple 
function Base.size(obj::MeshProduct)
    return 1
end

function rank(obj::MeshProduct)
    return length(mlist)
end

#TODO: find the linearindex I corresponding to given index
function index_to_linear(obj::MeshProduct, index...)
    return 1
end

#TODO: find the index corresponding to linearindex I
function linear_to_index(obj::MeshProduct)
    return 1
end



#TODO:for all n meshes in mlist, return [..., (mlist[i])[index[i]], ...] 
function Base.getindex(obj::MeshProduct, index...)
    return 1
end
#TODO:find the index corresponds to linearindex I, then for all n meshes in mlist, return [..., (mlist[i])[index[i]], ...] 
function Base.getindex(obj::MeshProduct, I::Int)
    return 1
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

