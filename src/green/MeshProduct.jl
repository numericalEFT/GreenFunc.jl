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
        mprod = Tuple(vargs...)
        mnew = new{typeof(mprod)}(
        mprod)
        return mnew
    end
end

function Base.length(obj::MeshProduct)

end


function Base.size(obj::MeshProduct, dim::Int)


end

function rank(obj::MeshProduct)


end



function Base.getindex(obj::MeshProduct, indices...)


end

function Base.getindex(obj::MeshProduct, index::Int)


end

function Base.view(obj::MeshProduct,inds...)
    

end

function Base.setindex!(obj::MeshProduct, index::Int)


end


Base.firstindex(obj::MeshProduct) = 1
Base.lastindex(obj::MeshProduct) = length(obj)
# iterator
Base.iterate(obj::MeshProduct) = (obj[1],1)
Base.iterate(obj::MeshProduct, state) = (state>=length(obj)) ? nothing : (obj[state+1],state+1)


function Base.iterate(obj::MeshProduct)
    return Iterators.product(obj.mlist)
end

"""
    def copy(self):
        
Deep copy

        return self.__class__(*[x.copy() for x in self._mlist])
"""

function Base.deepcopy(obj::MeshProduct)

end


"""
def copy_from(self, another):
Deep copy

        assert self.rank == another.rank, "copy_from requires the same rank for meshes"
        return self.__class__(*[x.copy_from(y) for x,y in zip(self._mlist, another._mlist)])

"""
function copy_from(obj::MeshProduct, src::MeshProduct)

end


function index_to_linear(obj::MeshProduct)

end

function Base.show(obj::MeshProduct)

end
