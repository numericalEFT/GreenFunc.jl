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
Base.length(obj::MeshProduct) = reduce(*, length(m) for m in obj.mlist)

#TODO:return the size of Ith mesh in mlist 
Base.size(obj::MeshProduct, I::Int) = length(obj.mlist[I])

#TODO:return the size of all meshes in a tuple 
Base.size(obj::MeshProduct)=Tuple(length.(obj.mlist))

rank(obj::MeshProduct) = length(obj.mlist)

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
Base.getindex(obj::MeshProduct, index...) = Tuple(obj.mlist[i][id] for (i,id) in enumerate(index))

#TODO:find the index corresponds to linearindex I, then for all n meshes in mlist, return [..., (mlist[i])[index[i]], ...] 
function Base.getindex(obj::MeshProduct, I::Int)
    index =linear_to_index(obj,I)
    return Tuple(obj.mlist[i][id] for (i,id) in enumerate(index))
end

#TODO:return the sliced pieces of 
# function Base.view(obj::MeshProduct,inds...)
#     return Tuple(view(obj, i) for i in inds)
#     #return 1
# end
# t[1] --> view of the first mesh


Base.firstindex(obj::MeshProduct) = 1
Base.lastindex(obj::MeshProduct) = length(obj)
# iterator
Base.iterate(obj::MeshProduct) = (obj[1],1)
Base.iterate(obj::MeshProduct, state) = (state>=length(obj)) ? nothing : (obj[state+1],state+1)


#TODO:nice print
Base.show(obj::MeshProduct)=string("MeshProduct of:\n",join((repr(x) for x in obj.mlist),",\n"))


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

